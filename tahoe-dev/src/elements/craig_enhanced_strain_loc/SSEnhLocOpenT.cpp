#include "SSEnhLocOpenT.h"
#include "ShapeFunctionT.h"
#include "SSSolidMatT.h"

/* materials lists */
#include "SSSolidMatList1DT.h"
#include "SSSolidMatList2DT.h"
#include "SSSolidMatList3DT.h"
#include "SSMatSupportT.h"
#include "ParameterContainerT.h"
#include "ModelManagerT.h"
#include "OutputSetT.h"

#include "DetCheckT.h"
#include "dTensor4DT.h"
#include "math.h"

using namespace Tahoe;

/*initialize static variables */
bool SSEnhLocOpenT::fLocalizationHasBegun = false;
bool SSEnhLocOpenT::fSeedElementsSet = false;
double SSEnhLocOpenT::fDetAMin = 1.0e99;
int SSEnhLocOpenT::fLeastDetEle = -1;

/* constructor */
SSEnhLocOpenT::SSEnhLocOpenT(const ElementSupportT& support):
  SmallStrainT(support),
  fOpenBand(NULL)
{
  SmallStrainT::SetName("small_strain_enh_loc_open");
  	//jump_out.open("jump.info");
	fEdgeOfBandElements.Free();
	fEdgeOfBandCoords.Free();
}

/* implementation of the ParameterInterfaceT interface */
void SSEnhLocOpenT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SmallStrainT::DefineParameters(list);

	/*PARAMETERS FOR Open model*/
	list.AddParameter(fFrictionCoeff, "friction_coefficient");
	list.AddParameter(fAlpha_zeta, "shear_tensile_displacement_ratio_alpha_delta");
	list.AddParameter(fAlpha_sigma, "shear_tensile_strength_ratio_alpha_sigma");
	list.AddParameter(fZeta_star, "critical_slip_distance_zeta_star");

	ParameterT multiBand(fMultiBand, "Allow_Multiple_Bands");
	multiBand.SetDefault(false);
	list.AddParameter(multiBand);
	
	
	list.AddParameter(fBVPType, "BVP_type");
	list.AddParameter(fFirstElementToLocalize, "First_element_to_localize");
	
	double dummy;
	
	list.AddParameter(dummy, "First_component_of_normal_if_prefailed");
	list.AddParameter(dummy, "Second_component_of_normal_if_prefailed");
	list.AddParameter(dummy, "First_component_of_slip_dir_if_prefailed");
	list.AddParameter(dummy, "Second_component_of_slip_dir_if_prefailed");

	//ParameterT initialSlipRate(fInitialSlipRate, "Initial_Slip_Rate_neg_one_for_V_star");
	//initialSlipRate.SetDefault(-1.0);
	//list.AddParameter(initialSlipRate);
}

/* information about subordinate parameter lists */
void SSEnhLocOpenT::DefineSubs(SubListT& sub_list) const
{	
  /* inherited */
  SmallStrainT::DefineSubs(sub_list);
}

void SSEnhLocOpenT::TakeParameterList(const ParameterListT& list)
{
  const char caller[] = "SmallStrainT::TakeParameterList";
  /* strain displacement option before calling
     SolidElementT::TakeParameterList */
  int b = list.GetParameter("strain_displacement");
  fStrainDispOpt = (b == kStandardB) ? kStandardB : kMeanDilBbar;
  
  /* inherited */
  SolidElementT::TakeParameterList(list);

  /* dimension workspace */
  fGradU.Dimension(NumSD());        
  if (fStrainDispOpt == kMeanDilBbar) {
    fLocDispTranspose.Dimension(fLocDisp.Length());
    fMeanGradient.Dimension(NumSD(), NumElementNodes());
  }        
  
  /* offset to class needs flags */
  fNeedsOffset = fMaterialNeeds[0].Length();
  
  /* set material needs */
  for (int i = 0; i < fMaterialNeeds.Length(); i++)
    {      
      /* needs array */
      ArrayT<bool>& needs = fMaterialNeeds[i];
      
      /* resize array */
      needs.Resize(needs.Length() + 2, true);
      
      /* casts are safe since class contructs materials list */
      ContinuumMaterialT* pcont_mat = (*fMaterialList)[i];
      SSSolidMatT* mat = (SSSolidMatT*) pcont_mat;
      
      /* collect needs */
      needs[fNeedsOffset + kstrain     ] = true; //mat->Need_Strain();
      needs[fNeedsOffset + kstrain_last] = true; //mat->Need_Strain_last();

      /* consistency */
      needs[kNeedDisp] = needs[kNeedDisp] || needs[fNeedsOffset +
						   kstrain];
      needs[KNeedLastDisp] = needs[KNeedLastDisp] || needs[fNeedsOffset
							   + kstrain_last];
    }

  /* what's needed */
  bool need_strain = true;
  bool need_strain_last = true;
  for (int i = 0; i < fMaterialNeeds.Length(); i++) {
    const ArrayT<bool>& needs = fMaterialNeeds[i];
    need_strain = need_strain || needs[fNeedsOffset + kstrain];
    need_strain_last = need_strain_last || needs[fNeedsOffset +
						 kstrain_last];
  }
  
  /* allocate strain list */
  if (need_strain) {
    fStrain_List.Dimension(NumIP());
    for (int i = 0; i < NumIP(); i++)
      fStrain_List[i].Dimension(NumSD());
  }
  
  /* allocate "last" strain list */
  if (need_strain_last) {
    fStrain_last_List.Dimension(NumIP());
    for (int i = 0; i < NumIP(); i++)
      fStrain_last_List[i].Dimension(NumSD());
  }


  /*PARAMETERS FOR Open Model*/	
  fFrictionCoeff = list.GetParameter("friction_coefficient");
  fAlpha_zeta = list.GetParameter("shear_tensile_displacement_ratio_alpha_delta");
  fAlpha_sigma = list.GetParameter("shear_tensile_strength_ratio_alpha_sigma");
  fZeta_star = list.GetParameter("critical_slip_distance_zeta_star");
  
  fMultiBand = list.GetParameter("Allow_Multiple_Bands");
  fBVPType = list.GetParameter("BVP_type");
  fFirstElementToLocalize = list.GetParameter("First_element_to_localize");
  //fFirstElementToLocalize -= 1;	

  if (fBVPType == kPreFailed)
    {
      fPreFailedNormal.Dimension(2);
      fPreFailedSlipDir.Dimension(2);

      fPreFailedNormal[0] = list.GetParameter("First_component_of_normal_if_prefailed");
      fPreFailedNormal[1] = list.GetParameter("Second_component_of_normal_if_prefailed");
      fPreFailedSlipDir[0] = list.GetParameter("First_component_of_slip_dir_if_prefailed");
      fPreFailedSlipDir[1] = list.GetParameter("Second_component_of_slip_dir_if_prefailed");
	  
	  fPreFailedNormal.UnitVector();
	  fPreFailedSlipDir.UnitVector();
	  
	  //PreFailElements();
    }
}

/* calculate the internal force contribution ("-k*d") */
void SSEnhLocOpenT::FormKd(double constK)
{

  if(!IsElementTraced())
  {
    SmallStrainT::FormKd(constK);
	cout << "not traced\n";
	}  
  else 
      {
	cout << "SSEnhLocOpenT::FormKd\n";
	  
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/* collect incremental heat */
	bool need_heat = fElementHeat.Length() == fShapes->NumIP();
	
	SetGlobalShape(); //do I need to write new version - think no
	dArrayT jumpIncr = CalculateJump();//fOpenBand -> Jump();
	/* record jump in fOpenBand */
	fOpenBand -> StoreJump(jumpIncr);
	jumpIncr -= fOpenBand -> LastJump();
	
	fShapes->TopIP();
	while (fShapes->NextIP())
	  {
	    /* strain displacement matrix */
	    if (fStrainDispOpt == kMeanDilBbar)
	      Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
	    else
	      Set_B(fShapes->Derivatives_U(), fB);

	    /* B^T * Cauchy stress */
		dSymMatrixT strainIncr = fStrain_List [CurrIP()];
		strainIncr -= fStrain_last_List [CurrIP()]; 
		  
		dSymMatrixT gradActiveTensorJumpIncr =
		FormGradActiveTensorJumpIncr(jumpIncr, NumSD(), CurrIP());
		gradActiveTensorJumpIncr.ScaleOffDiagonal(0.5);      	

		strainIncr.AddScaled(-1.0, gradActiveTensorJumpIncr);
			
		dSymMatrixT stressIncr(NumSD());
		stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);
		//cout << "CurrIP() = " << CurrIP() << endl;
		stressIncr += fOpenBand->Stress_List(CurrIP());
	    fB.MultTx(stressIncr, fNEEvec);
	    
	    /* accumulate */
	    fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);
	    
	    /* incremental heat generation */
	    if (need_heat) 
	      fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	  }     
      }
}


void SSEnhLocOpenT::FormStiffness(double constK)
{

  if (fStrainDispOpt == kMeanDilBbar)
    cout << "Warning SSEnhLocCraigT::FormStiffness, b-bar integration not implemented for enhanced strain element, postlocalization results may be garbage.";

  if (!IsElementTraced() || fOpenBand->BandState() == kInactive )
  {
    /* form stiffness in standard way */
    SmallStrainT::FormStiffness(constK);

  }
  else 
  {

	/* matrix format */
	dMatrixT::SymmetryFlagT format = dMatrixT::kWhole;

	/* integrate element stiffness */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	int ndof = NumDOF();
	int nen = NumElementNodes();
	int nedof = nen*ndof;//# of element dof
	double area = 0.0;
	dArrayT jump = fOpenBand->Jump();
	
	if (jump[0] == 0.0 )
	{ 	
		double k_zeta_zeta = 0.0;
		dArrayT k_d_zeta(nedof), k_zeta_d(nedof);
		dArrayT k_d_zeta_work(nedof), k_zeta_d_work(nedof);
		k_d_zeta = 0.0;
		k_zeta_d = 0.0;

		dSymMatrixT gradActiveTensorShearDir(ndof), dGdSigma(ndof), nTensorn(ndof);
		dMatrixT fLHSWork(nedof),fDfB((fCurrMaterial->ce_ijkl().Rows()),nedof);

		dGdSigma = FormdGdSigma(ndof);
	
		fShapes->TopIP();
		while ( fShapes->NextIP() )
		{
			area += (*Det)*(*Weight);
			double scale = constK*(*Det++)*(*Weight++);
	
			//form ke_dd - assume elastic
			/* strain displacement matrix */
			if (fStrainDispOpt == kMeanDilBbar)
				Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
			else
				Set_B(fShapes->Derivatives_U(), fB);

			/* get D matrix */
			fD.SetToScaled(scale, fCurrMaterial->ce_ijkl());

			/* multiply b(transpose) * db, taking account of symmetry, */
			/* and accumulate in elstif */
			fDfB.MultAB(fD, fB);
			fLHSWork.MultATB(fB, fDfB);
			fLHS +=fLHSWork;
			/* update nMatrixT to get rid of fLHSwork */

			//form k_d_zeta
			gradActiveTensorShearDir =
			FormGradActiveTensorShearDir(ndof, CurrIP());

			fDfB.MultTx(gradActiveTensorShearDir, k_d_zeta_work);
			k_d_zeta += k_d_zeta_work;

			//form k_zeta_d
			fDfB.MultTx(dGdSigma, k_zeta_d_work);
			k_zeta_d += k_zeta_d_work;
		}
	
		k_zeta_d *= 1.0/area;

		if (fOpenBand->BandState() == kElastic)
			k_zeta_zeta = ElasticBalanceGradient(jump)(1,1);
		else 
			k_zeta_zeta = PlasticBalanceGradient(jump)(1,1);

		k_zeta_zeta *= -1.0;

	fLHS.Outer(k_d_zeta, k_zeta_d, -1.0/k_zeta_zeta, dMatrixT::kAccumulate);
    }
	else
	{
		dMatrixT k_zeta_zeta(ndof); 
		dMatrixT k_d_zeta(nedof, ndof), k_zeta_d(ndof, nedof);
		dMatrixT k_d_zeta_work(nedof, ndof), k_zeta_d_work(ndof, nedof);
		k_zeta_zeta = 0.0;
		k_d_zeta = 0.0;
		k_zeta_d = 0.0;

		dMatrixT gradActiveTensorBandDirs(dSymMatrixT::NumValues(ndof), ndof);
		dMatrixT dGsdSigma(dSymMatrixT::NumValues(ndof), ndof);
		dMatrixT fLHSWork(nedof),fDfB((fCurrMaterial->ce_ijkl().Rows()),nedof);

		dGsdSigma = FormdGsdSigma(ndof);
	
		fShapes->TopIP();
		while ( fShapes->NextIP() )
		{
			area += (*Det)*(*Weight);
			double scale = constK*(*Det++)*(*Weight++);
	
			//form ke_dd - assume elastic
			/* strain displacement matrix */
			if (fStrainDispOpt == kMeanDilBbar)
				Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
			else
				Set_B(fShapes->Derivatives_U(), fB);

			/* get D matrix */
			fD.SetToScaled(scale, fCurrMaterial->ce_ijkl());

			/* multiply b(transpose) * db, taking account of symmetry, */
			/* and accumulate in elstif */
			fDfB.MultAB(fD, fB);
			fLHSWork.MultATB(fB, fDfB);
			fLHS +=fLHSWork;
			/* update nMatrixT to get rid of fLHSwork */

			//form k_d_zeta
			gradActiveTensorBandDirs =
			FormGradActiveTensorBandDirs(ndof, CurrIP());

			//fDfB.MultTx(gradActiveTensorBandDirs, k_d_zeta_work); //change to MultAB?????
			k_d_zeta_work.MultATB(fDfB, gradActiveTensorBandDirs);
			k_d_zeta += k_d_zeta_work;

			//form k_zeta_d
			fDfB.MultAB(dGsdSigma, k_zeta_d_work);
			k_zeta_d += k_zeta_d_work;
		}
	
		k_zeta_d *= 1.0/area;

		if (fOpenBand->BandState() == kElastic)
			k_zeta_zeta = ElasticBalanceGradient(jump);
		else 
			k_zeta_zeta = PlasticBalanceGradient(jump);

		//k_zeta_zeta *= -1.0;
		
		double det = k_zeta_zeta(0,0) * k_zeta_zeta(1,1)
				-  k_zeta_zeta(0,1) * k_zeta_zeta(1,0);
					
		dMatrixT k_zeta_zetaInverse(2);
		k_zeta_zetaInverse(0,0) = k_zeta_zeta(1,1);
		k_zeta_zetaInverse(0,1) = -1.0* k_zeta_zeta(0,1);
		k_zeta_zetaInverse(1,0) = -1.0* k_zeta_zeta(1,0);
		k_zeta_zetaInverse(1,1) = k_zeta_zeta(0,0);
		k_zeta_zetaInverse /= det;

	fLHS.MultABC(k_d_zeta, k_zeta_zetaInverse, k_zeta_d, dMatrixT::kWhole, dMatrixT::kAccumulate); 
	}
  }
}
  

dArrayT SSEnhLocOpenT::CalculateJump()
{

	dArrayT jump = fOpenBand -> LastJump();
	double normalStress = NormalStress(jump); 
	if (normalStress < 0.0 &&
	    ShearStress(jump) < -1.0 * fFrictionCoeff * normalStress)
		{
			fOpenBand -> SetBandState(kInactive);
			return jump;
		}

   if (fOpenBand -> Cohesion() == 0.0)
   {
	  fOpenBand -> SetBandState(kDamage);
	  return CalculatePlasticJump();
	}
	
	jump = CalculateElasticJump();

    //cout << "jump = " << jump << endl;

	if (EffectiveTraction(jump) > fOpenBand->Cohesion())
	{

		jump = CalculatePlasticJump();
		//cout << "jump = " << jump << endl;
		fOpenBand -> SetBandState(kDamage);
	}	
	else
	{
		cout << "jump is elastic\n";
		fOpenBand -> SetBandState(kElastic);
	}

	//cout << "fOpenBand = " << fOpenBand << endl;

		
	return jump;
}

double SSEnhLocOpenT::EffectiveTraction(dArrayT jump)
{
	double normalStress = NormalStress(jump);
	double shearStress = ShearStress(jump);
	
	return sqrt(normalStress*normalStress*fAlpha_sigma*fAlpha_sigma + shearStress*shearStress);
}

dArrayT SSEnhLocOpenT::CalculateElasticJump()
{

  double newtonTol = 1.0e-12;
  double slipRateTol = 1.0e-12;

  
  int newtonCounter = 0;
  int maxIter = 30;

  dArrayT jump = fOpenBand->LastJump();
  
  dArrayT tractionBalance = ElasticTractionBalance(jump);
  double initialNormalBalance = tractionBalance[0];
  double initialShearBalance = tractionBalance[1];

  while ((fabs(tractionBalance[0]/initialNormalBalance) > newtonTol && fabs(tractionBalance[0]) > newtonTol)
			|| (fabs(tractionBalance[1]/initialShearBalance) > newtonTol && fabs(tractionBalance[1]) > newtonTol))
    {
      /* if too many iterations, stop and cut load step */
      if(newtonCounter++ >maxIter)
		{
 			cout << "SSEnhLocOpenT::CalculateJumpIncrement, Newton iteration did not converge, cutting load step\n";
			throw ExceptionT::kGeneralFail;
		}
		
      /* update jump increment via Newton iteration */
      dMatrixT balanceGradient = ElasticBalanceGradient(jump);
	  
	  /* generate inverse */
	  double det = balanceGradient(0,0) * balanceGradient(1,1)
					-  balanceGradient(0,1) * balanceGradient(1,0);
					
	  /* if band is closed and in compression, opening jump is 0, solve 1x1 for shear */
	  if (det == 0.0)
	  {
		jump[0] = 0.0;
		jump[1] -= tractionBalance[1]/balanceGradient(1,1);
	  }
	  /* otherwise solve 2x2 for both degrees of freedom */
	  else
	  {
		dMatrixT balanceGradientInverse(2);
		balanceGradientInverse(0,0) = balanceGradient(1,1);
		balanceGradientInverse(0,1) = -1* balanceGradient(0,1);
		balanceGradientInverse(1,0) = -1* balanceGradient(1,0);
		balanceGradientInverse(1,1) = balanceGradient(0,0);
		balanceGradientInverse /= det;
	 
	    dArrayT jumpIncrement(NumSD());
		balanceGradientInverse.Multx(tractionBalance, jumpIncrement);
	 
		jump -=  jumpIncrement;
		if (jump[0] < 0.0)
			jump[0] = 0.0;
	   }
	 tractionBalance = ElasticTractionBalance(jump);

      //slipRate -= yieldFn/DPhidSlipRate(slipRate, jumpIncrement, thetaNew);

      /*update increment of ISV */
      //jumpIncrement = JumpIncrement(slipRate);
      //thetaNew = ThetaNew(slipRate);

      /*reform DeltaG*/
      //yieldFn = Phi(slipRate, jumpIncrement, thetaNew);
    }
  
/* store? */
//store current jump, cohesion?
return jump;
 
}      

dArrayT SSEnhLocOpenT::ElasticTractionBalance(dArrayT jump)
{
	dArrayT tractionBalance(2);
	double normalStress = NormalStress(jump);
	double shearStress = ShearStress(jump);
	double effectiveJump = EffectiveJump(jump);
	
	if (fOpenBand->Cohesion() == fOpenBand->InitialCohesion() && jump == 0.0)
	{
		tractionBalance = 0.0;
		return tractionBalance;
	}
	
	/* normal traction balance */
	if (normalStress <= 0.0)
	{
	    if (jump[0] == 0.0)
			tractionBalance[0] = 0.0;
		else 
		tractionBalance[0] =  normalStress - 
				jump[0] * fAlpha_zeta * fOpenBand->Cohesion() * fOpenBand->InitialCohesion() * fZeta_star
								/ (fAlpha_sigma * (fOpenBand->InitialCohesion() - fOpenBand->Cohesion()));
		
		tractionBalance[1] = shearStress - 
							jump[1] * fOpenBand->Cohesion() * fOpenBand->InitialCohesion() * fZeta_star
							/  (fOpenBand->InitialCohesion() - fOpenBand->Cohesion()); 
		
	    double friction =  normalStress * fFrictionCoeff;
		if (tractionBalance[1] + friction < 0.0)
			tractionBalance[1] = 0.0;
		else
			tractionBalance += friction;		
	}
	else
	{
		tractionBalance[0] =  normalStress - 
				jump[0] * fAlpha_zeta * fOpenBand->Cohesion() * fOpenBand->InitialCohesion() * fZeta_star
								/ (fAlpha_sigma * (fOpenBand->InitialCohesion() - fOpenBand->Cohesion()));
								
		tractionBalance[1] = shearStress - 
							jump[1] * fOpenBand->Cohesion() * fOpenBand->InitialCohesion() * fZeta_star
							/  (fOpenBand->InitialCohesion() - fOpenBand->Cohesion()); 
	}			
	
	return tractionBalance;
}

dMatrixT SSEnhLocOpenT::ElasticBalanceGradient(dArrayT jump)
{
	dMatrixT balanceGradient(2);
	
	double normalStress = NormalStress(jump);
	double shearStress = ShearStress(jump);
	double effectiveJump = EffectiveJump(jump);
	
	
	if (normalStress <= 0.0 && jump[0] == 0.0)
	{
     	balanceGradient(0,0) = 0.0;
		balanceGradient(1,0) = 0.0;

	  
	  balanceGradient(1,0) = DShearStressDNormalSlip(); //normal slip does not occur - no effect on friction?
	  balanceGradient(1,1) = DShearStressDShearSlip() -
				fOpenBand->Cohesion() * fOpenBand->InitialCohesion() * fZeta_star
				/ (fOpenBand->InitialCohesion() - fOpenBand->Cohesion())
				+ fFrictionCoeff * DNormalStressDShearSlip(); //if no slip, won't come back here?
		
	}
	else
	{
		balanceGradient(0,0) = DNormalStressDNormalSlip() - 
				fAlpha_zeta * fOpenBand->Cohesion() * fOpenBand->InitialCohesion() * fZeta_star
				/ (fAlpha_sigma * (fOpenBand->InitialCohesion() - fOpenBand->Cohesion()));
		balanceGradient(0,1) = DNormalStressDShearSlip();
		balanceGradient(1,0) = DShearStressDNormalSlip();
		balanceGradient(1,1) = DShearStressDShearSlip() -
				fOpenBand->Cohesion() * fOpenBand->InitialCohesion() * fZeta_star
				/ (fOpenBand->InitialCohesion() - fOpenBand->Cohesion());
	}
	return balanceGradient;
}


double SSEnhLocOpenT::DNormalStressDNormalSlip()
{
    int ndof = NumDOF();
	dArrayT gradActive = AvgGradActive(NumDOF());
	dArrayT normal = fOpenBand -> Normal();
	
	/* 4-way contraction */
	return -1.0 * QuadContraction(fCurrMaterial->ce_ijkl(), normal, normal, normal, gradActive);
}

double SSEnhLocOpenT::DNormalStressDShearSlip()
{
    int ndof = NumDOF();
	dArrayT gradActive = AvgGradActive(NumDOF());
	dArrayT normal = fOpenBand -> Normal();
	dArrayT shearDir = fOpenBand -> ShearDir();
	
	/* 4-way contraction */
	return -1.0 * QuadContraction(fCurrMaterial->ce_ijkl(), normal, normal, shearDir, gradActive);
}

double SSEnhLocOpenT::DShearStressDNormalSlip()
{
    int ndof = NumDOF();
	dArrayT gradActive = AvgGradActive(NumDOF());
	dArrayT normal = fOpenBand -> Normal();
	dArrayT shearDir = fOpenBand -> ShearDir();
	
	/* 4-way contraction */
	return -1.0 * QuadContraction(fCurrMaterial->ce_ijkl(), normal, shearDir, normal, gradActive);
}

double SSEnhLocOpenT::DShearStressDShearSlip()
{
    int ndof = NumDOF();
	dArrayT gradActive = AvgGradActive(NumDOF());
	dArrayT normal = fOpenBand -> Normal();
	dArrayT shearDir = fOpenBand -> ShearDir();
	
	/* 4-way contraction */
	return -1.0 *  QuadContraction(fCurrMaterial->ce_ijkl(), normal, shearDir, shearDir, gradActive);
}

double SSEnhLocOpenT::QuadContraction(dMatrixT c_ijkl, 
		dArrayT vector1, dArrayT vector2, dArrayT vector3, dArrayT vector4)
{
    int ndof = NumDOF(); //change to vector1 length
	dTensor4DT C(ndof, ndof, ndof, ndof);
	
	C.ConvertTangentFrom2DTo4D(C, fCurrMaterial->ce_ijkl());
	
	double d = 0.0;
	
	for (int i = 0; i < ndof; i++)
		for (int j = 0; j < ndof; j++)
			for (int k = 0; k < ndof; k++)
				for (int l = 0; l < ndof; l++)
					d += C(i,j,k,l) * vector1[i] * vector2[j] * vector3[k] * vector4[l];
	return d;
}


dArrayT SSEnhLocOpenT::AvgGradActive(int ndof)
{
  iAutoArrayT activeNodes = fOpenBand->ActiveNodes();
  int A;
  dArrayT grad_f(ndof);

  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();      
  double area = 0.0;
  dArrayT avgGradActive(ndof);
  avgGradActive = 0.0;

  fShapes->TopIP();
  while (fShapes->NextIP())
    {
      double scale = (*Det++)*(*Weight++);
      area += scale;

	  //fShapes->SetIP(ip);
	  fShapes->SetIP(CurrIP());
      Set_B(fShapes->Derivatives_U(), fB);

      for (int i=0; i<ndof; i++)
      {
	  grad_f = 0.0;
      activeNodes.Top();

      while(activeNodes.Next())      
		{
			A = activeNodes.Current();  
			grad_f[i] += fB(i, (A)*ndof +i);
		}
	   }
	   grad_f *= scale;
	   avgGradActive += grad_f; 
    }
  avgGradActive /= area;
  return avgGradActive;		
}

dArrayT SSEnhLocOpenT::CalculatePlasticJump()
{
	 double newtonTol = 1.0e-12;
  double slipRateTol = 1.0e-12;

  
  int newtonCounter = 0;
  int maxIter = 30; 

  dArrayT jump = fOpenBand->LastJump();
  
  dArrayT tractionBalance = PlasticTractionBalance(jump);
  double initialNormalBalance = tractionBalance[0];
  double initialShearBalance = tractionBalance[1];

  cout << "tractionBalance = " << tractionBalance << endl;

  while ((fabs(tractionBalance[0]/initialNormalBalance) > newtonTol && fabs(tractionBalance[0]) > newtonTol)
			|| (fabs(tractionBalance[1]/initialShearBalance) > newtonTol && fabs(tractionBalance[1]) > newtonTol))
    {
	  cout << "hi\n";
	
      /* if too many iterations, stop and cut load step */
      if(newtonCounter++ >maxIter)
		{
			cout << "SSEnhLocOpenT::CalculateJumpIncrement, Newton iteration did not converge, cutting load step\n";
			throw ExceptionT::kGeneralFail;
		}
		
      /* update jump increment via Newton iteration */
      dMatrixT balanceGradient = PlasticBalanceGradient(jump);
	  
	  cout << "balanceGradient = " << balanceGradient << endl; 
	  
	  /* generate inverse */
	  double det = balanceGradient(0,0) * balanceGradient(1,1)
					-  balanceGradient(0,1) * balanceGradient(1,0);
					
	  /* if band is closed and in compression, opening jump is 0, solve 1x1 for shear */
	  if (jump[0] <= 0.0 && NormalStress(jump) <= 0.0)
	  {
		jump[0] = 0.0;
		jump[1] -= tractionBalance[1]/balanceGradient(1,1);
	  }
	  /* otherwise solve 2x2 for both degrees of freedom */
	  else
	  {
		dMatrixT balanceGradientInverse(2);
		balanceGradientInverse(0,0) = balanceGradient(1,1);
		balanceGradientInverse(0,1) = -1.0 * balanceGradient(0,1);
		balanceGradientInverse(1,0) = -1.0 * balanceGradient(1,0);
		balanceGradientInverse(1,1) = balanceGradient(0,0);
		balanceGradientInverse /= det;

	    dArrayT jumpIncrement(NumSD());
		balanceGradientInverse.Multx(tractionBalance, jumpIncrement);
	 
		jump -=  jumpIncrement;
		if (jump[0] < 0.0)
			jump[0] = 0.0;
			
		cout << "jump = " << jump << endl;
	   }
	 tractionBalance = PlasticTractionBalance(jump);
    }
	return jump;
}

dArrayT SSEnhLocOpenT::PlasticTractionBalance(dArrayT jump)
{
		double cohesion = CurrentCohesion(jump);
		
		if (cohesion <= 0.0)
			return CohesionlessTractionBalance(jump);
			
			dArrayT tractionBalance(2);
	double normalStress = NormalStress(jump);
	double shearStress = ShearStress(jump);
	double effectiveJump = EffectiveJump(jump);
	
	
	if (fOpenBand->Cohesion() == fOpenBand->InitialCohesion() && jump == 0.0)
	{
	
		if (normalStress <= 0.0)
		{
			tractionBalance[0] = 0.0;
			tractionBalance[1] = shearStress - normalStress * fFrictionCoeff - fOpenBand -> Cohesion();
		}
		else
		{
		
		    double beta = fAlpha_zeta * shearStress / (fAlpha_sigma * normalStress);
			if (beta == 0.0)
			{
				tractionBalance[0] = normalStress 
					- fOpenBand -> Cohesion() * fAlpha_zeta / fAlpha_sigma;
				tractionBalance[1] = 0.0;
			}
			else
			{
				double work = 1.0 / (sqrt((fAlpha_zeta *fAlpha_zeta/ (beta * beta)) + 1));
				tractionBalance[0] = normalStress 
					- fOpenBand -> Cohesion() * fAlpha_zeta * work/ (fAlpha_sigma * beta);
				tractionBalance[1] = shearStress - fOpenBand -> Cohesion() * work;
			}
		}
		
		return tractionBalance;
	}
	
	/* normal traction balance */
	if (normalStress <= 0.0)
	{
	    if (jump[0] == 0.0)
			tractionBalance[0] = 0.0;
		else 
		tractionBalance[0] =  normalStress - 
				jump[0] * fAlpha_zeta * cohesion * fOpenBand->InitialCohesion() * fZeta_star
								/ (fAlpha_sigma * (fOpenBand->InitialCohesion() - cohesion));
		
		tractionBalance[1] = shearStress - 
							jump[1] * cohesion * fOpenBand->InitialCohesion() * fZeta_star
							/  (fOpenBand->InitialCohesion() - cohesion); 
		
	    double friction =  normalStress * fFrictionCoeff;
		if (tractionBalance[1] + friction < 0.0)
			tractionBalance[1] = 0.0;
		else
			tractionBalance += friction;		
	}
	else
	{
		tractionBalance[0] =  normalStress - 
				jump[0] * fAlpha_zeta * cohesion * fOpenBand->InitialCohesion() * fZeta_star
								/ (fAlpha_sigma * (fOpenBand->InitialCohesion() - cohesion));
								
		tractionBalance[1] = shearStress - 
							jump[1] * cohesion * fOpenBand->InitialCohesion() * fZeta_star
							/  (fOpenBand->InitialCohesion() - cohesion); 
	}			
	
	return tractionBalance;
		
		  
}

/* returns current cohesion given jump */
double SSEnhLocOpenT::CurrentCohesion(dArrayT jump)
{
    double effectiveJump = EffectiveJump(jump);
	
	if (effectiveJump >= fZeta_star)
		return 0.0;
		
	return fOpenBand->InitialCohesion() * (1.0 - effectiveJump)/effectiveJump;
}

double SSEnhLocOpenT::EffectiveJump(dArrayT jump)
{
	return sqrt(fAlpha_zeta * fAlpha_zeta * jump[0] * jump[0] + jump[1] * jump[1]);
}

dMatrixT SSEnhLocOpenT::PlasticBalanceGradient(dArrayT jump)
{
    dMatrixT balanceGradient(2);
	double cohesion = CurrentCohesion(jump);
		
	if (cohesion <= 0.0)
		return CohesionlessBalanceGradient(jump);
		
	double normalStress = NormalStress(jump);
	double shearStress = ShearStress(jump);
	double effectiveJump = EffectiveJump(jump);
	
	//add in plastic increments 
	if (normalStress <= 0.0 && jump[0] == 0.0)
	{
	  balanceGradient(0,0) = 0.0;
	  balanceGradient(1,0) = 0.0;

	  
	  balanceGradient(1,0) = DShearStressDNormalSlip(); //normal slip does not occur - no effect on friction?
	  balanceGradient(1,1) = DShearStressDShearSlip() -
				cohesion * fOpenBand->InitialCohesion() * fZeta_star
				/ (fOpenBand->InitialCohesion() - cohesion)
				+ fFrictionCoeff * DNormalStressDShearSlip(); 
				 
	  balanceGradient(1,1) += 2.0 * fOpenBand->InitialCohesion() * jump[1] * jump[1] 
							/ (effectiveJump * effectiveJump * effectiveJump);
	}
	else
	{
		balanceGradient(0,0) = DNormalStressDNormalSlip() - 
				fAlpha_zeta * cohesion * fOpenBand->InitialCohesion() * fZeta_star
				/ (fAlpha_sigma * (fOpenBand->InitialCohesion() - cohesion));
		balanceGradient(0,1) = DNormalStressDShearSlip();
		balanceGradient(1,0) = DShearStressDNormalSlip();
		balanceGradient(1,1) = DShearStressDShearSlip() -
				cohesion * fOpenBand->InitialCohesion() * fZeta_star
				/ (fOpenBand->InitialCohesion() - cohesion);
				
		balanceGradient(0,0) += 2.0 * fAlpha_zeta * fOpenBand->InitialCohesion() * jump[0] * jump[0] 
							/ (effectiveJump * effectiveJump * effectiveJump);
		balanceGradient(0,1) += 2.0 * fAlpha_zeta * fOpenBand->InitialCohesion() * jump[0] * jump[1] 
							/ (fAlpha_sigma * effectiveJump * effectiveJump * effectiveJump);
		balanceGradient(1,0) += 2.0 * fAlpha_zeta * fOpenBand->InitialCohesion() * jump[1] * jump[0] 
							/ ( effectiveJump * effectiveJump * effectiveJump);
		balanceGradient(1,1) += 2.0 * fOpenBand->InitialCohesion() * jump[1] * jump[1] 
							/ (effectiveJump * effectiveJump * effectiveJump);
	}
	return balanceGradient;
}


dArrayT SSEnhLocOpenT::CohesionlessTractionBalance(dArrayT jump)
{
	dArrayT tractionBalance(2);
	double normalStress = NormalStress(jump);
	
	
	/* normal traction balance */
	if (normalStress <= 0.0 && jump[0])
	{
		tractionBalance[0] = 0.0;
	
		tractionBalance[1] = ShearStress(jump) + normalStress * fFrictionCoeff;
	}
	else
	{
		tractionBalance[0] = normalStress;
		tractionBalance[1] = ShearStress(jump);
	}			
	
	return tractionBalance;
}

dMatrixT SSEnhLocOpenT::CohesionlessBalanceGradient(dArrayT jump)
{
	dMatrixT balanceGradient(2);
	double normalStress = NormalStress(jump);
	
	
	/* normal traction balance */
	if (normalStress <= 0.0 && jump[0])
	{
     	balanceGradient(0,0) = 0.0;
		balanceGradient(1,0) = 0.0;
	  
		balanceGradient(1,0) = DShearStressDNormalSlip(); //normal slip does not occur - no effect on friction?
		balanceGradient(1,1) = DShearStressDShearSlip() +
				 fFrictionCoeff * DNormalStressDShearSlip(); //if no slip, won't come back here?
	}
	else
	{
		balanceGradient(0,0) = DNormalStressDNormalSlip(); 
		balanceGradient(0,1) = DNormalStressDShearSlip();
		balanceGradient(1,0) = DShearStressDNormalSlip();
		balanceGradient(1,1) = DShearStressDShearSlip();
	}			
	
	return balanceGradient;
}

dSymMatrixT SSEnhLocOpenT::StressIncrOnBand(dArrayT jumpIncrement)
{
  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();      
  double area = 0.0;
  int ndof = NumDOF();
  dSymMatrixT gradActiveTensorJumpIncr(ndof);
  dSymMatrixT avgStressIncr(ndof);
  avgStressIncr = 0.0;

  fShapes->TopIP();
  while (fShapes->NextIP())
    {
      double scale = (*Det++)*(*Weight++);
      area += scale;
       
	  /* for no jump increment */ 
      dSymMatrixT strainIncr = fStrain_List [CurrIP()];
      strainIncr -= fStrain_last_List [CurrIP()];

      /* add effect of jump */
	  //strainIncr.AddScaled( /* Tensor(GradActive, jumpIncr */);  
	  
      gradActiveTensorJumpIncr = FormGradActiveTensorJumpIncr(jumpIncrement, ndof, CurrIP());
      //gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);
      strainIncr.AddScaled(-1.0, gradActiveTensorJumpIncr);

      dSymMatrixT stressIncr(NumSD());
      stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);

      avgStressIncr.AddScaled(scale, stressIncr);
    }
  avgStressIncr/=area;

  return avgStressIncr;
}

dSymMatrixT SSEnhLocOpenT::LastStressOnBand()
{
  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();      
  double area = 0.0;
  dSymMatrixT avgStress(NumDOF());
  avgStress = 0.0;

  fShapes->TopIP();
  while (fShapes->NextIP())
    {
      double scale = (*Det++)*(*Weight++);
      area += scale;

      avgStress.AddScaled(scale, fOpenBand->Stress_List(CurrIP())); //fOpenBand???
    }
  avgStress/=area;

  return avgStress;
}


double SSEnhLocOpenT::NormalStress(dArrayT jump)
{
  dSymMatrixT currStress = LastStressOnBand();
  dArrayT jumpIncr = jump;
  jumpIncr -= fOpenBand->LastJump();
  
  currStress += StressIncrOnBand(jumpIncr);

  return currStress.MultmBn(fOpenBand->Normal(), fOpenBand->Normal());
}

double SSEnhLocOpenT::ShearStress(dArrayT jump)
{
  dSymMatrixT currStress = LastStressOnBand();
  
  dArrayT jumpIncr = jump;
  jumpIncr -= fOpenBand->LastJump();
  currStress += StressIncrOnBand(jumpIncr);

  return currStress.MultmBn(fOpenBand->ShearDir(), fOpenBand->Normal());
}

/*------end math functions----------------------------------------*/

/* can get rid of slipDir, perpSlipDir later */
OpenBandT* SSEnhLocOpenT::FormNewBand(dArrayT normal, dArrayT slipDir,
				   dArrayT perpSlipDir, dArrayT coords, double area)
{
  ArrayT<dSymMatrixT> stressList;
  stressList.Dimension(NumIP());

  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();
  //double area = 0.0;

  /*for residual cohesion*/
  double normalStress = 0.0;
  double shearStress = 0.0;

  Det    = fShapes->IPDets();
  Weight = fShapes->IPWeights();

  fShapes->TopIP();
  for (int i = 0; i < NumIP(); i++)
    {
      fShapes -> NextIP();
      stressList[i].Dimension(NumSD());
      stressList[i] = fCurrMaterial -> s_ij();

      double scale = (*Det++)*(*Weight++);

      normalStress += scale * stressList[i].MultmBn(normal, normal);
      shearStress += scale * stressList[i].MultmBn(perpSlipDir, normal); 
    }

  normalStress/= area;
  shearStress = fabs(shearStress)/area;

  double residCohesion;
  
  if (normalStress < 0.0)
  {
	residCohesion = shearStress;
  	residCohesion += normalStress * fFrictionCoeff;
  }
  else
  {
    double weightedNormalStress = fAlpha_sigma * normalStress;
	residCohesion = sqrt(shearStress * shearStress + weightedNormalStress * weightedNormalStress);
  }

  if (fBVPType == kPreFailed && ElementSupport().Time() == 0.0)
    residCohesion = 0.0;

return new OpenBandT(normal, perpSlipDir, coords, residCohesion, stressList, this);
}

void SSEnhLocOpenT::CloseStep()
{

cout << "\n CloseStep \n \n";

/*update traced elements */ 
if (fLocalizationHasBegun)
	{
	  Top();
      while (NextElement())
	{
	  if (IsElementTraced())
	    {
		SetGlobalShape();	

	    	/* loop over integration points */
	    	fShapes->TopIP();
	      while (fShapes->NextIP())
			{    
	      	dSymMatrixT strainIncr = fStrain_List [CurrIP()];
	      	strainIncr -= fStrain_last_List [CurrIP()]; 
	      
		    dArrayT jumpIncr = fOpenBand -> Jump();
			jumpIncr -= fOpenBand -> LastJump();
		  
	     	dSymMatrixT gradActiveTensorJumpIncr =
			FormGradActiveTensorJumpIncr(jumpIncr, NumSD(), CurrIP());
	      	gradActiveTensorJumpIncr.ScaleOffDiagonal(0.5);      	

	      	strainIncr.AddScaled(-1.0, gradActiveTensorJumpIncr);
			
	     	dSymMatrixT stressIncr(NumSD());
	      	stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);
	      	fOpenBand -> IncrementStress(stressIncr, CurrIP());
			
	      	}
		  /* Update Cohesion */
		  if (fOpenBand -> BandState() == kDamage)
			fOpenBand -> UpdateCohesion(1.0 - EffectiveJump(fOpenBand -> Jump() )/fZeta_star);
			
	      fOpenBand -> CloseStep();
	      	
			//cout << "residual cohesion = " << fBand->ResidualCohesion();
	      	//cout << "fBand->JumpIncrement = " << fBand->JumpIncrement();
	      	//cout << ", fBand->Jump() = " << fBand->Jump() << endl;	
	    }			
	}
	}

 if (fBVPType == kPreFailed && ElementSupport().Time() == 0.0)
   {
     //cout << "time = 0.0\n" << flush; 
     PreFailElements();
   }

bool finishedTracing = false;
	
/* check for newly localized elements on existing bands */
//fJustChecking = false;
fEdgeOfBandElements.Top();
fEdgeOfBandCoords.Top();

while(!finishedTracing)
{
while(fEdgeOfBandElements.Next())
	{
	  fEdgeOfBandCoords.Next();
	  cout << "fEdgeOfBandElements.Current() = " << flush << fEdgeOfBandElements.Current() << endl << flush;
	  //cout << "fEdgeOfBandCoords.Current() = " << fEdgeOfBandCoords.Current() << endl;	  
	  GetElement(fEdgeOfBandElements.Current());
	  if (IsElementLocalized())
		TraceElement();
	}

//reset to check for elements not on an existing band
fDetAMin = 1.0e99;
finishedTracing = true;
//fEdgeOfBandElements.Previous();
fEdgeOfBandElements.Current(fEdgeOfBandElements.Position() -1);

/* check for newly localized elements not on existing bands */
//fJustChecking = true;
	  //cout << "hi hi\n";
if (!fSeedElementsSet)
if (!fLocalizationHasBegun || fMultiBand) 
    {
	  cout << "hi\n";
      //choose first element then let band progress
 
      Top();
      while (NextElement())
	{
		  //cout << "1\n";
		if (!IsElementTraced())
		{
		    //cout << "CurrElementNumber = " << CurrElementNumber() << "\n";
			GetElement(CurrElementNumber());

			if (IsElementLocalized())
				finishedTracing = false;
		}
	}
			if (!finishedTracing)
			{
				fLocalizationHasBegun = true;
				//fEdgeOfBandElements.Current(0); //should be at
	  
	  			cout << "fLeastDetEle = " << fLeastDetEle << endl;
				GetElement(fLeastDetEle);
				//fEdgeOfBandCoords.Free();
				fEdgeOfBandElements.Append(fLeastDetEle);
				fEdgeOfBandCoords.Append(Centroid());	  
				//fEdgeOfBandElements.Top();
				//fEdgeOfBandCoords.Top();
			}
    }
	
	
}


 /* uncomment these two lines to restrict band propagation after initial
    localization */
//if (fLocalizationHasBegun)
//	fSeedElementsSet = true;

  SmallStrainT::CloseStep();

}

bool SSEnhLocOpenT::IsElementTraced()
{
 	 //int elementNumber = CurrElementNumber();
 	 return IsElementTraced(CurrElementNumber());
}
  
bool SSEnhLocOpenT::IsElementTraced(int elementNumber)
{
  bool isTraced = fTracedElements.HasKey(elementNumber);

  if (isTraced)
    LoadBand(elementNumber);

  return isTraced;
}
	
bool SSEnhLocOpenT::IsElementLocalized()
{
  /* newly localized elements only */
  if (IsElementTraced())
	return false;

  bool locCheck = false;
  double detA, detAMin  = 1.0e99;
  AutoArrayT <dArrayT> normals;
  AutoArrayT <dArrayT> slipDirs;
  //AutoArrayT <dArrayT> bestNormals;
  //AutoArrayT <dArrayT> bestSlipDirs;

  /* loop over integration points */
  fShapes->TopIP();
  while ( fShapes->NextIP() )
    {
      //if (CurrIP() == 0 && CurrElementNumber() == 0)
      //cout << "stress = " << fCurrMaterial->s_ij() << endl;

      //checker.SetfStructuralMatSupport(*fSSMatSupport);
      if (fCurrMaterial->IsLocalized(normals, slipDirs, detA))
	{
	  locCheck = true;
	  normals.Top();
	  while (normals.Next())

	  if (detA < detAMin)
	    detAMin = detA;

	  //cout << "detAMin = " << flush << detAMin << flush << endl;
	  //cout << "fDetAMin = " << fDetAMin << endl << endl;
	
      if ((fBVPType == kNonhomogeneous && detAMin < fDetAMin) || (fBVPType !=
		kNonhomogeneous && CurrElementNumber() == fFirstElementToLocalize))
      {
	  
		fDetAMin = detAMin;
		fLeastDetEle = CurrElementNumber();
		//cout << "detAMin = " << detAMin << endl;
		//cout << "fDetAMin = " << fDetAMin << endl;
		//cout << "fLeastDetEle = " << fLeastDetEle << endl;
		//save best normals and slipdirs here?
      }
	}
	}

  return locCheck;
}

bool SSEnhLocOpenT::TraceElement()
{
  cout << "Tracing element " << CurrElementNumber() << endl;

  bool locCheck = false;
  double detA, detAMin  = 1.0e99;
  AutoArrayT <dArrayT> normals;
  AutoArrayT <dArrayT> slipDirs;
  AutoArrayT <dArrayT> bestNormals;
  AutoArrayT <dArrayT> bestSlipDirs;

  /* loop over integration points */
  fShapes->TopIP();
  while ( fShapes->NextIP() )
    {
      //if (CurrIP() == 0 && CurrElementNumber() == 0)
      //cout << "stress = " << fCurrMaterial->s_ij() << endl;

      //checker.SetfStructuralMatSupport(*fSSMatSupport);
      if (fCurrMaterial->IsLocalized(normals, slipDirs, detA))
	{
	  locCheck = true;
	  normals.Top();
	  while (normals.Next())

	  if (detA < detAMin)
	    {
	      detAMin = detA;
	      bestNormals.Dimension(normals);
	      normals.CopyInto(bestNormals);
	      bestSlipDirs.Dimension(slipDirs);
	      slipDirs.CopyInto(bestSlipDirs);
	    }
	}
    }

  if (locCheck)
    {
      ChooseNormals(bestNormals, bestSlipDirs);
      /* remove element from list of Edge elements*/
      fEdgeOfBandElements.DeleteAt(fEdgeOfBandElements.Position());
      fEdgeOfBandElements.Current(fEdgeOfBandElements.Position() -1);
      fEdgeOfBandCoords.DeleteAt(fEdgeOfBandCoords.Position());  
      fEdgeOfBandCoords.Current(fEdgeOfBandCoords.Position() - 1);
     }

  return locCheck;
}

void SSEnhLocOpenT::ChooseNormals(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipDirs)
{
  cout << "choosing normals\n";
  normals.Top();
  slipDirs.Top();

  int ndof = NumDOF();

  dArrayT normal(ndof);
  dArrayT slipDir(ndof);
  dArrayT perpSlipDir(ndof);

  //determine strain on band
  //should we use elastic strain?
  dMatrixT avgGradU(ndof);
  avgGradU = 0.0;

  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();
  double area = 0.0;

  fShapes->TopIP();
  while (fShapes-> NextIP())
    {
      double scale = (*Det++)*(*Weight++);
      area += scale;
      fShapes->GradU(fLocDisp, fGradU, CurrIP());
      avgGradU.AddScaled(scale, fGradU);
      //cout << "fGradU = \n" << fGradU << endl;
    }
  avgGradU /= area;

  double prod, maxProd = -1.0;

  while(normals.Next())
    {
      //cout << "current normal = \n" << normals.Current() << endl; 
      slipDirs.Next();
      //prod = fabs( avgGradU.MultmBn(normals.Current(), slipDirs.Current()));
      prod = fabs( avgGradU.MultmBn(slipDirs.Current(), normals.Current()));
	  //cout << "prod = " << prod << endl;
      //prod = (normals.Current() [1]) * (normals.Current() [2]);
      if (prod > maxProd)
	{
	  //cout << "best normal = \n" << normals.Current() << endl; 
	  normal = normals.Current(); 
	  slipDir = slipDirs.Current();
          maxProd = prod;
	}
    }
	
	
 //1 For propagating band normal in same direction
 #if 0
   OpenBandT* fBandTemp = fOpenBand;
 
   //get element neighbors array
   ModelManagerT& model = ElementSupport().ModelManager();
  iArray2DT neighbors;
  ArrayT<StringT> ids;
  
  ElementBlockIDs(ids);  
  //cout << "hi " << flush;

  model.ElementNeighbors(ids, neighbors);

  //cout << "neighbors =\n" << neighbors << endl;

  //2D
  int numSides;
  int numSidesFound = 0;

  switch (GeometryCode())
{
 case GeometryT::kQuadrilateral:
   {
     numSides = 4;
     break;
   }
 case GeometryT::kTriangle: 
   {
     numSides = 3;
     break;
   }
 default:
   {
     cout << "SSEnhLocCraigT::AddNewEdgeElements, geometry not implemented. \n" << flush;
     throw ExceptionT::kGeneralFail;
   }
}

   int elementNumber = CurrElementNumber();
   //search for traced neighbor
    for(int i = 0; i < numSides; i++)
   		if (neighbors(elementNumber,i) != -1 && IsElementTraced(neighbors(elementNumber ,i)))
         {   
			//load that element - done by IsElementTraced
   
			//get that normal
			normal = fOpenBand -> Normal();
         }
 
 
  fOpenBand = fBandTemp;
 #endif

  //make sure slip direction is dilatant
  if (normal.Dot(normal, slipDir)<0.0)
    slipDir *= -1.0;

  perpSlipDir = slipDir;
  perpSlipDir.AddScaled(-1.0*slipDir.Dot(slipDir, normal), normal);
  perpSlipDir.UnitVector();

/*
  if (fNoBandDilation)
    {
      slipDir = perpSlipDir;
    }
*/

  fOpenBand = FormNewBand(normal, slipDir, perpSlipDir, fEdgeOfBandCoords.Current(), area);

  fTracedElements.Insert(CurrElementNumber(), fOpenBand);
  AddNewEdgeElements(CurrElementNumber());
}

void SSEnhLocOpenT::AddNewEdgeElements(int elementNumber)
{
  ModelManagerT& model = ElementSupport().ModelManager();
  iArray2DT neighbors;
  ArrayT<StringT> ids;
  
  ElementBlockIDs(ids);  
  //cout << "hi " << flush;

  model.ElementNeighbors(ids, neighbors);

  //cout << "neighbors =\n" << neighbors << endl;

  //2D
  int numSides;
  int numSidesFound = 0;
  //int numElementSides = 0;
  
 
  
  iAutoArrayT activeNodes = fOpenBand->ActiveNodes();
  //activeNodes.Free();
  //activeNodes.Dimension(fBand->ActiveNodes());
  //fBand->ActiveNodes().CopyInto(activeNodes);



  switch (GeometryCode())
{
 case GeometryT::kQuadrilateral:
   {
     numSides = 4;
     break;
   }
 case GeometryT::kTriangle: 
   {
     numSides = 3;
     break;
   }
 default:
   {
     cout << "SSEnhLocCraigT::AddNewEdgeElements, geometry not implemented. \n" << flush;
     throw ExceptionT::kGeneralFail;
   }
}


  LocalArrayT nodalCoords = InitialCoordinates();
  dArrayT nodalCoord1(NumSD()), nodalCoord2(NumSD()); //coords a particular node

  //cout << "fBand = " << fBand << endl;
  
  //temporary workaround- IsElementTraced loads the element, changing the coordinates
  OpenBandT* fBandTemp = fOpenBand;

  for(int i = 0; i < numSides; i++)
    //if (((fBand->ActiveNodes()).HasValue((i+1) % numSides) && !((fBand->ActiveNodes()).HasValue(i)))
	//|| (!((fBand->ActiveNodes()).HasValue((i+1) % numSides)) && (fBand->ActiveNodes()).HasValue(i)))
	if ((activeNodes.HasValue((i+1) % numSides) && !(activeNodes.HasValue(i)))
	|| (!(activeNodes.HasValue((i+1) % numSides)) && activeNodes.HasValue(i)))
      {

	  //cout << "fBand = " << fBand << endl;	
      fOpenBand = fBandTemp;
	  
	
		if (!(neighbors(elementNumber,i) == -1 || IsElementTraced(neighbors(elementNumber ,i))))
	  	{
	    	//get coords
			fOpenBand = fBandTemp;
			  //cout << "fBand = " << fBand << endl;
	    	dArrayT localizedEleCoord = fOpenBand -> Coords();

		    for (int j = 0; j < nodalCoords.MinorDim(); j++)
		      {		
				nodalCoord1 [j] = nodalCoords(i,j);
				nodalCoord2 [j] = nodalCoords((i+1) % numSides, j); //replace numSides w/ num Ele nodes!
	      		}
		    /*
	      	cout << CurrElementNumber() << endl;
	      	cout << "nodalCoords =\n" << nodalCoords;
	      	cout << "neighbors =\n" << neighbors << endl;	
	      	cout << "nodalCoord1 = \n" << nodalCoord1 << endl;	
	      	cout << "nodalCoord2 = \n" << nodalCoord2 << endl;
		    */


	    	dArrayT interceptCoords = InterceptCoords(localizedEleCoord,
						      nodalCoord1, nodalCoord2);
	    	//stick element in fEdgeOfBandElements
	    	//EdgeOfBandElement* newElement = new EdgeOfBandElement;
	    	//newElement->elementNumber = CurrElementNumber();
	    	//newElement -> coords = interceptCoords;
	    	if(fEdgeOfBandElements.AppendUnique(neighbors(elementNumber,i)))
		  {
		    cout << "Successfully appended element " <<
		      neighbors(elementNumber,i) << endl << "coords =\n" << interceptCoords << endl;
		    fEdgeOfBandCoords.Append(interceptCoords);
		  }
	  	} 
		if (++numSidesFound > 1) 
	  		break;
      }
  cout << "numSidesFound = " << numSidesFound <<endl;

}

dArrayT SSEnhLocOpenT::InterceptCoords(dArrayT& localizedEleCoord,
dArrayT& nodalCoord1, dArrayT& nodalCoord2)
{
  //assumes straight sides
  dArrayT sideVector = nodalCoord2;
  sideVector -= nodalCoord1;

  dArrayT perpSlipDir = fOpenBand -> ShearDir();

  //cout << "Interceptcoords: nodalCoord1 =\n" << nodalCoord1 << "\nNodalCoord2 =\n" << nodalCoord2
  //     << "\nlocalizedEleCoord =\n" << localizedEleCoord << "\nperpSlipDir =\n" << perpSlipDir; 
	   
  double alpha = sideVector[0] * (localizedEleCoord[1] - nodalCoord1[1]) -
  sideVector[1] * (localizedEleCoord[0] - nodalCoord1[0]);
  alpha /= sideVector[1] * perpSlipDir[0] - sideVector[0] *
  perpSlipDir[1];

  dArrayT interceptCoord = localizedEleCoord;
  interceptCoord.AddScaled(alpha, perpSlipDir);

  return interceptCoord;
}


dArrayT SSEnhLocOpenT::Centroid()
{
  dArrayT centroid(NumDOF());

  double area = 0.0;
  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();
  dArrayT coords(NumDOF());

  centroid = 0.0;

  fShapes->TopIP();
  while ( fShapes->NextIP() )
    {
      double scale = (*Det++)*(*Weight++);
      area += scale;
      fShapes->IPCoords(coords);
      centroid.AddScaled(scale,coords);
    }
  centroid *= 1.0/area;
  
  return centroid;
 
}

void SSEnhLocOpenT::LoadBand(int elementNumber)
{
   fOpenBand = fTracedElements[elementNumber];

  //SSEnhLocCraigT::LoadBand(elementNumber);
  //change to static cast?
  //fOpenBand = dynamic_cast<OpenBandT*> (fOpenBand);
}

/* current element operations */
void SSEnhLocOpenT::GetElement(int elementNumber)
{
  /* inherited */
  //bool result = ContinuumElementT::NextElement();
  fElementCards.Current(elementNumber);  
  //cout << "elementNumber =" << elementNumber << endl;

  /* get material pointer */
  ContinuumMaterialT* pcont_mat = (*fMaterialList)[CurrentElement().MaterialNumber()];
      
  /* cast is safe since class contructs materials list */
  fCurrMaterial = (SolidMaterialT*) pcont_mat;
  
  //SetLocalX(fLocInitCoords);
  SetGlobalShape();
 }

void SSEnhLocOpenT::PreFailElements()
{
	fLocalizationHasBegun = true;
	
	AutoArrayT<dArrayT> normals;
	AutoArrayT<dArrayT> slipDirs;
	
	normals.Append(fPreFailedNormal);
	slipDirs.Append(fPreFailedSlipDir);

	GetElement(fFirstElementToLocalize);
	fEdgeOfBandCoords.Append(Centroid());
	fEdgeOfBandCoords.Current(0);		
	ChooseNormals(normals, slipDirs);
	fEdgeOfBandCoords.DeleteAt(0);
  
	fEdgeOfBandElements.Top();
	fEdgeOfBandCoords.Top();
	
	while(fEdgeOfBandElements.Next())
	    {
	      cout << "checking for localization of element " <<
	      	fEdgeOfBandElements.Current() << endl; 

		  fEdgeOfBandCoords.Next();
	      GetElement(fEdgeOfBandElements.Current());    
		  ChooseNormals(normals, slipDirs);	
	
      /* remove element from list of Edge elements*/
      fEdgeOfBandElements.DeleteAt(fEdgeOfBandElements.Position());
      fEdgeOfBandElements.Current(fEdgeOfBandElements.Position() -1);
      fEdgeOfBandCoords.DeleteAt(fEdgeOfBandCoords.Position());  
      fEdgeOfBandCoords.Current(fEdgeOfBandCoords.Position() - 1);	
	}
}


dSymMatrixT SSEnhLocOpenT::FormdGdSigma(int ndof)
{
  dSymMatrixT dGdSigma(ndof), work(ndof);
  dMatrixT dGNonSym(ndof);

  dGNonSym.Outer(fOpenBand->Normal(), fOpenBand-> ShearDir());
  dGdSigma.Symmetrize(dGNonSym);
 
  work.Outer(fOpenBand->Normal());
  dGdSigma.AddScaled(fFrictionCoeff, work);
 
  /* use vector values for mult */
  dGdSigma.ScaleOffDiagonal(2.0);

  return dGdSigma;
}

dMatrixT SSEnhLocOpenT::FormdGsdSigma(int ndof)
{
  dSymMatrixT dGdSigma(ndof), work(ndof);
  dMatrixT dGNonSym(ndof), dGs(dSymMatrixT::NumValues(ndof), ndof);

  dGNonSym.Outer(fOpenBand->Normal(), fOpenBand->ShearDir());
  dGdSigma.Symmetrize(dGNonSym);
 
  work.Outer(fOpenBand->Normal());
 
  /* use vector values for mult */
  work.ScaleOffDiagonal(2.0);
  dGdSigma.ScaleOffDiagonal(2.0);

  for (int i=0; i < dSymMatrixT::NumValues(ndof); i++)
  {
	dGs(i, 0) = work[i];
	dGs(i, 1) = dGdSigma[i];
  }

  return dGs;
}

dSymMatrixT SSEnhLocOpenT::FormGradActiveTensorJumpIncr(dArrayT jumpIncr, int ndof, int ip)
{
  dSymMatrixT G(ndof);
  dMatrixT G_NonSym(ndof);
  iAutoArrayT activeNodes = fOpenBand->ActiveNodes();
  int A;
  dArrayT grad_f(ndof);
  grad_f = 0.0;

  fShapes->SetIP(ip);
  Set_B(fShapes->Derivatives_U(), fB);

  for (int i=0; i<ndof; i++)
    {
      activeNodes.Top();

      while(activeNodes.Next())      
	{
	  A = activeNodes.Current();  
	  grad_f[i] += fB(i, (A)*ndof +i);
	}
    }

  G_NonSym.Outer(grad_f, jumpIncr);
  G.Symmetrize(G_NonSym);

  /* use vector values for mult */
  //G.ScaleOffDiagonal(2.0);

  return G;
}

dSymMatrixT SSEnhLocOpenT::FormGradActiveTensorShearDir(int ndof, int ip)
{
  dSymMatrixT G(ndof);
  dMatrixT G_NonSym(ndof);
  iAutoArrayT activeNodes = fOpenBand->ActiveNodes();
  int A;
  dArrayT grad_f(ndof);
  grad_f = 0.0;

  fShapes->SetIP(ip);
  Set_B(fShapes->Derivatives_U(), fB);

  for (int i=0; i<ndof; i++)
    {
      activeNodes.Top();

      while(activeNodes.Next())      
	{
	  A = activeNodes.Current();  
	  grad_f[i] += fB(i, (A)*ndof +i);
	}
    }

  G_NonSym.Outer(grad_f, fOpenBand->ShearDir());
  G.Symmetrize(G_NonSym);

  /* use vector values for mult */
  G.ScaleOffDiagonal(2.0);

  return G;
}

dMatrixT SSEnhLocOpenT::FormGradActiveTensorBandDirs(int ndof, int ip)
{
  dSymMatrixT G(ndof);
  dMatrixT G_NonSym(ndof);
  iAutoArrayT activeNodes = fOpenBand->ActiveNodes();
  int A, vert_length =  dSymMatrixT::NumValues(ndof);
  G.Symmetrize(G_NonSym);

  /* use vector values for mult */
  G.ScaleOffDiagonal(2.0);;
  dArrayT grad_f(ndof);
  grad_f = 0.0;
  

  fShapes->SetIP(ip);
  Set_B(fShapes->Derivatives_U(), fB);

  for (int i=0; i<ndof; i++)
    {
      activeNodes.Top();

      while(activeNodes.Next())      
	{
	  A = activeNodes.Current();  
	  grad_f[i] += fB(i, (A)*ndof +i);
	}
    }

  G_NonSym.Outer(grad_f, fOpenBand->Normal());
  G.Symmetrize(G_NonSym);

  /* use vector values for mult */
  G.ScaleOffDiagonal(2.0);
  
  dMatrixT gatbd(vert_length, ndof);

  for (int i=0; i < vert_length; i++)
    gatbd(0, i) = G[i];

  G_NonSym.Outer(grad_f, fOpenBand->ShearDir());
  G.Symmetrize(G_NonSym);

  /* use vector values for mult */
  G.ScaleOffDiagonal(2.0);  
  for (int i=0; i < vert_length; i++)
    gatbd(1, i) = G[i];

  return gatbd;
}