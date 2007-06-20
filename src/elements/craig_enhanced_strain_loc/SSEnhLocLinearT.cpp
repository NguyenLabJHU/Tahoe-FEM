 #include "SSEnhLocLinearT.h"
#include "ShapeFunctionT.h"
#include "SSSolidMatT.h"
#include "dTensor4DT.h"

/* materials lists */
#include "SSSolidMatList1DT.h"
#include "SSSolidMatList2DT.h"
#include "SSSolidMatList3DT.h"
#include "SSMatSupportT.h"
#include "ParameterContainerT.h"
#include "ModelManagerT.h"
#include "OutputSetT.h"

#include "DetCheckT.h"
#include <math.h>

using namespace Tahoe;

/*initialize static variables */
bool SSEnhLocLinearT::fLocalizationHasBegun = false;
bool SSEnhLocLinearT::fSeedElementsSet = false;
double SSEnhLocLinearT::fDetAMin = 1.0e99;
int SSEnhLocLinearT::fLeastDetEle = -1;

/* constructor */
SSEnhLocLinearT::SSEnhLocLinearT(const ElementSupportT& support):
	SmallStrainT(support),
	fBand(NULL)
	//fLeastDetEle(-1)
{
	SmallStrainT::SetName("small_strain_enh_loc_linear");
	jump_out.open("jump.info");
	fEdgeOfBandElements.Free();
	fEdgeOfBandCoords.Free();
}

/* destructor */
/*
SSEnhLocLinearT::~SSEnhLocLinearT(void)
{}
*/

/* implementation of the ParameterInterfaceT interface */
void SSEnhLocLinearT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SmallStrainT::DefineParameters(list);

	/*PARAMETERS FOR ENHANCED STRAIN*/
	list.AddParameter(fH_delta_0, "Post-Localization_softening_parameter_H_Delta"); 
	list.AddParameter(fNoBandDilation, "Disallow_Dilation_on_Band");
	list.AddParameter(fLocalizedFrictionCoeff, "Localized_Friction_Coefficient");
	//list.AddParameter(fMultiBand, "Allow_Multiple_Bands");
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

    //list.AddParameter(fPreFailedNormal[0], "First_component_of_normal_if_prefailed");
	//list.AddParameter(fPreFailedNormal[1], "Second_component_of_normal_if_prefailed");
	//list.AddParameter(fPreFailedSlipDir[0], "First_component_of_slip_dir_if_prefailed");
	//list.AddParameter(fPreFailedSlipDir[1], "Second_component_of_slip_dir_if_prefailed");
	
}

/* information about subordinate parameter lists */
void SSEnhLocLinearT::DefineSubs(SubListT& sub_list) const
{	
  /* inherited */
  SmallStrainT::DefineSubs(sub_list);
}

void SSEnhLocLinearT::TakeParameterList(const ParameterListT& list)
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

  /*PARAMETERS FOR ENHANCED STRAIN*/
  fH_delta_0 = list.GetParameter("Post-Localization_softening_parameter_H_Delta"); 
  fNoBandDilation = list.GetParameter("Disallow_Dilation_on_Band");
  fLocalizedFrictionCoeff = list.GetParameter("Localized_Friction_Coefficient");
  fMultiBand = list.GetParameter("Allow_Multiple_Bands");
  fBVPType = list.GetParameter("BVP_type");
  fFirstElementToLocalize = list.GetParameter("First_element_to_localize");
  fFirstElementToLocalize -= 1;	//adjust to internal Tahoe numbering

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

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct a new material support and return a pointer */
MaterialSupportT* SSEnhLocLinearT::NewMaterialSupport(MaterialSupportT* p) const
{
  return SmallStrainT::NewMaterialSupport(p);
}


/* return a pointer to a new material list */
MaterialListT* SSEnhLocLinearT::NewMaterialList(const StringT& name, int size)
{
  return SmallStrainT::NewMaterialList(name, size);
}


/* calculate the internal force contribution ("-k*d") */
void SSEnhLocLinearT::FormKd(double constK)
{

  if(!IsElementTraced() ||  (fBand -> JumpIncrement(0) == 0.0 && fBand -> JumpIncrement(1) == 0.0))
    SmallStrainT::FormKd(constK);
  else 
      {
	//cout << "SSEnhLocLinearT::FormKd\n";
	  
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/* collect incremental heat */
	bool need_heat = fElementHeat.Length() == fShapes->NumIP();
	
	fShapes->TopIP();
	while (fShapes->NextIP())
	  {
	    /* strain displacement matrix */
	    if (fStrainDispOpt == kMeanDilBbar)
	      Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
	    else
	      Set_B(fShapes->Derivatives_U(), fB);

	    /* Conforming Strain Increment */
	    dSymMatrixT strainIncr = fStrain_List [CurrIP()];
	    strainIncr -= fStrain_last_List [CurrIP()];

	    /* modify strain increment for avg jump Increment */
	    dSymMatrixT gradActiveTensorFlowDir =
	    FormGradActiveTensorFlowDir(NumSD(), CurrIP());
	    gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);

		//double avgJumpIncr = 0.0;
		//for (int i = 0; i < fBand -> NumSurfaceIPs(); i ++)
		//	avgJumpIncr += fBand -> JumpIncrement(i);
		//avgJumpIncr /= fBand -> NumSurfaceIPs();	

		//cout << "strainIncr = \n" << strainIncr << endl;
		//cout << "gradActiveTensorFlowDir = \n" << gradActiveTensorFlowDir << endl;
		//cout << "avgJumpIncr = " << avgJumpIncr << endl;
		//cout << "JumpIncrAtBulkIP(CurrIP()) = " << JumpIncrAtBulkIP(CurrIP()) << endl;
		
	    strainIncr.AddScaled(-1.0*JumpIncrAtBulkIP(CurrIP()), gradActiveTensorFlowDir);

		//cout << "strainIncr = \n" << strainIncr << endl;
		
		/* modify strain increment for jump increment gradient */
		dSymMatrixT perpSlipDirOuter(NumSD());
		dArrayT perpSlipDir = fBand -> PerpSlipDir(0);
		perpSlipDirOuter.Outer(perpSlipDir);
		
		dMatrixT ce_inverse = fCurrMaterial->ce_ijkl();
		ce_inverse.Inverse(); 
		double denom = QuadContraction(ce_inverse, perpSlipDir, perpSlipDir, perpSlipDir, perpSlipDir);

		double linearModifier;
		if (fBand -> IsBulkIPActive(CurrIP()))
			linearModifier = 1.0; 
		else
			linearModifier = 0.0;
				
		linearModifier -= F_hAtBulkIP(CurrIP());		
		linearModifier = (fBand -> JumpIncrement(1) - fBand -> JumpIncrement(0))/
								(fBand -> DistanceBetweenIPs() * denom);
								
		//cout << "linearModifier = " << linearModifier << endl;
		
		/*
		double linearModifier;
		if (fBand -> IsBulkIPActive(CurrIP()))
			linearModifier = 1.0; 
		else
			linearModifier = 0.0;
				
		linearModifier -= F_hAtBulkIP(CurrIP());
		linearModifier /= fBand -> DistanceBetweenIPs();
		linearModifier *= (fBand -> JumpIncrement(1) - fBand -> JumpIncrement(0));
		
	    strainIncr.AddScaled(1.0*linearModifier, perpSlipDirOuter);
		*/
		
		//cout << "strainIncr = \n" << strainIncr << endl;
		
		/* stress increment */
	    dSymMatrixT stressIncr(NumSD());
	    stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);
	    stressIncr.AddScaled(linearModifier, perpSlipDirOuter);
		stressIncr += fBand->Stress_List(CurrIP());
	    fB.MultTx(stressIncr, fNEEvec);
	    
	    /* accumulate */
	    fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);
		
		//cout << "fRHS = \n" << fRHS << endl;
	    
	    /* incremental heat generation */
	    if (need_heat) 
	      fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	  }     
      }
}

/* form the element stiffness matrix */
void SSEnhLocLinearT::FormStiffness(double constK)
{
  if (fStrainDispOpt == kMeanDilBbar)
    cout << "Warning SSEnhLocLinearT::FormStiffness, b-bar integration not implemented for enhanced strain element, postlocalization results may be garbage.";

  if (!IsElementTraced() || (fBand->JumpIncrement(0) == 0.0 && fBand -> JumpIncrement(1) == 0.0))//todo - more elegant?
    {
      /* form stiffness in standard way */
      SmallStrainT::FormStiffness(constK);
    }
  else //if already localized, use localized stiffness routine
    {
	/* matrix format */
	dMatrixT::SymmetryFlagT format = dMatrixT::kWhole;

	/* integrate element stiffness */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	int ndof = NumDOF();
	int nen = NumElementNodes();
	int nedof = nen*ndof;//# of element dof
	double k_zeta_zeta, area = 0.0;
	dArrayT k_d_zeta(nedof), k_zeta_d(nedof);
	dArrayT k_d_zeta_work(nedof), k_zeta_d_work(nedof);


	dSymMatrixT gradActiveTensorFlowDir(ndof), dGdSigma(ndof);
	dMatrixT fLHSWork(nedof),fDfB((fCurrMaterial->ce_ijkl().Rows()),nedof);

	dGdSigma = FormdGdSigma(ndof);

    /* form ke_dd */
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
                //fLHS.MultATB(fB, fD, format, dMatrixT::kAccumulate);	
	}

    /* loop over and form other matrices for each band ip */
	dSymMatrixT perpSlipDirOuter(NumSD());
    perpSlipDirOuter.Outer(fBand -> PerpSlipDir(0));
	/*scale for vector use */
    perpSlipDirOuter.ScaleOffDiagonal(2.0);

	for (int bandIP = 0; bandIP < fBand -> NumSurfaceIPs(); bandIP++)
	if (fBand -> JumpIncrement(bandIP) != 0.0)
	{
		k_d_zeta = 0.0;
		k_zeta_d = 0.0;
		k_zeta_zeta = 0.0;
		double modifier;

		/*reset IP Dets and Weights */
		Det    = fShapes->IPDets();
		Weight = fShapes->IPWeights();
	
		fShapes->TopIP();
		while ( fShapes->NextIP() )
		{
			//form k_d_zeta
			gradActiveTensorFlowDir =
			FormGradActiveTensorFlowDir(ndof, CurrIP());

			/*  modify to account for linear band */
			double bandIPcoord = fBand -> IPBandCoord(CurrIP());
			
			double interpolant;
			if (bandIP == 0)
				interpolant = 0.5 - bandIPcoord/(fBand -> DistanceBetweenIPs());
			else if (bandIP == 1 )
				interpolant = 0.5 + bandIPcoord/(fBand -> DistanceBetweenIPs());
			else 
			{
				cout << "SSEnhLocLinearT::FormStiffness: invalid band IP" << flush;
				throw ExceptionT::kGeneralFail;
			}
			//cout << "interpolant = " << interpolant << endl;
			
			gradActiveTensorFlowDir *= interpolant;
			
			
			if (fBand -> IsBulkIPActive(CurrIP()))
				modifier = 1.0; 
			else
			    modifier = 0.0;
				
			modifier -= F_hAtBulkIP(CurrIP());
			modifier /= fBand -> DistanceBetweenIPs();
			
			if (bandIP == 1)
				modifier *= -1.0;
			
			/*
			gradActiveTensorFlowDir.AddScaled(modifier, perpSlipDirOuter);
		    */
			
			dMatrixT ce_inverse = fCurrMaterial->ce_ijkl();
			ce_inverse.Inverse(); 
			dArrayT perpSlipDir = fBand -> PerpSlipDir(0);
			double denom = QuadContraction(ce_inverse, perpSlipDir, perpSlipDir, perpSlipDir, perpSlipDir);
		
			modifier /= (fBand -> DistanceBetweenIPs() * denom);
		
			double scale = constK*(*Det++)*(*Weight++);
		
			/* strain displacement matrix */
			if (fStrainDispOpt == kMeanDilBbar)
				Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
			else
				Set_B(fShapes->Derivatives_U(), fB);

			/* get D matrix */
			fD.SetToScaled(scale, fCurrMaterial->ce_ijkl());

			/* multiply b(transpose) * db, taking account of symmetry, */
			fDfB.MultAB(fD, fB);
			
			//cout << "fB = \n" << fB << endl;

			fDfB.MultTx(gradActiveTensorFlowDir, k_d_zeta_work);
			k_d_zeta += k_d_zeta_work;
			
			/* reset */
			k_d_zeta_work = 0.0;
			fB.Multx(perpSlipDirOuter, k_d_zeta_work);
			k_d_zeta_work *= scale * modifier;
			
			k_d_zeta += k_d_zeta_work; 
			
        }

		//form k_zeta_d
		fD.SetToScaled(constK, fCurrMaterial->ce_ijkl());

		/* multiply b(transpose) * db, taking account of symmetry, */
		dMatrixT bB = B_at_surfaceIP(bandIP);	
		fDfB.MultAB(fD, bB);
		
		//cout << "bB = \n" << bB << endl;
		
		//fDfB.MultTx(dGdSigma, k_zeta_d_work, dMatrixT::kOverwrite);
		fDfB.MultTx(dGdSigma, k_zeta_d_work);
		k_zeta_d += k_zeta_d_work;
		
		gradActiveTensorFlowDir = FormGradActiveTensorFlowDirAtBandIP(ndof, bandIP);

        //cout << "gradActiveTensorFlowDir =\n" << gradActiveTensorFlowDir << endl;
        //cout << "dGdSigma =\n" << dGdSigma << endl;		
		//form k_zeta_zeta
		k_zeta_zeta = fD.MultmBn(dGdSigma,gradActiveTensorFlowDir);

		//k_zeta_d *= 1.0/area;
		//k_zeta_zeta *= 1.0/area;
		
		/*
		cout << "k_d_zeta = " << k_d_zeta << endl;
		cout << "k_zeta_d  = " << k_zeta_d << endl;
		cout << "k_zeta_zeta = " << k_zeta_zeta << endl;
		*/
		
		k_zeta_zeta += fBand->EffectiveSoftening(bandIP);
		//cout << "fLHS =\n" << fLHS << endl;
		fLHS.Outer(k_d_zeta, k_zeta_d, -1.0/k_zeta_zeta, dMatrixT::kAccumulate);
    }
   //cout << "Element Number " << CurrElementNumber() << ". InitialCoodinates() =\n" << InitialCoordinates() << endl;

	}
			//cout << "fLHS =\n" << fLHS << endl;
}

/* compute the measures of strain/deformation over the element */
void SSEnhLocLinearT::SetGlobalShape(void)
{
  SmallStrainT::SetGlobalShape();
  
  /* subtract band deformation */
  if (IsElementTraced())
    {
      int ndof = NumDOF();
      dSymMatrixT gradActiveTensorFlowDir(ndof);
      
      if (fStrainDispOpt == kMeanDilBbar)
	cout << "Warning - B-bar not implemented for localized element\n";
      else
	{
	  int material_number = CurrentElement().MaterialNumber();
	  const ArrayT<bool>& needs = fMaterialNeeds[material_number];
	  
	  for (int i = 0; i < fBand -> NumSurfaceIPs(); i++)
	  {
	    //cout << "i = " << i << endl;
		double jumpIncrement = CalculateJumpIncrement(i);
	  }

#if 0
	  /* loop over integration points again */
	  for (int i = 0; i < NumIP(); i++)
	    {
	      gradActiveTensorFlowDir =
		FormGradActiveTensorFlowDir(ndof, i);
	      
	      /*change shear strains back to matrix values */
	      /* vector values are used when created */
	      gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);
				
	      /* deformation gradient */
	      if (needs[fNeedsOffset + kstrain])
		{			  
		  fStrain_List[i].AddScaled(-(fBand->Jump()+
		     jumpIncrement),gradActiveTensorFlowDir);
		}
		
	      /* "last" deformation gradient */ //is this right?
	      if (needs[fNeedsOffset + kstrain_last])
		{
		  fStrain_last_List[i].AddScaled(-(fBand->Jump()),
						 gradActiveTensorFlowDir);
		}
	    }
#endif
	}
    }
}

/***********************************************************************
 * Protected
 ***********************************************************************/



double SSEnhLocLinearT::CalculateJumpIncrement(int bandIP)
{
	int ndof = NumDOF();
  
	/* calculate conforming Strain increment*/
	dSymMatrixT conformingStrainIncrement = ConformingStrainIncrAtCoord(fBand -> SurfaceIPCoords(bandIP));
  
    //cout << "conformingStrainIncrement = \n" << conformingStrainIncrement << endl;
  
	/*determine trial traction*/
	dArrayT workingTraction = fBand -> SurfaceIPTraction(bandIP); 
  
	dSymMatrixT trialStressIncrement(ndof);
	trialStressIncrement.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), conformingStrainIncrement);
  
	dArrayT trialTractionIncrement(ndof);
	trialStressIncrement.Multx(fBand -> Normal(), trialTractionIncrement);
  
	workingTraction.AddScaled(1.0, trialTractionIncrement);	
	
	if (!IsBandActive(workingTraction, bandIP))
	{
		//cout << "Band IP inactive, bandIP = " << bandIP << "jumpIncrement = 0.0" << endl;
		fBand -> StoreJumpIncrement(0.0, bandIP);
		return 0.0;
	}
  
	/* claculate denominator */
	dSymMatrixT dGdSigma = FormdGdSigma(ndof);
	dSymMatrixT gradActiveTensorFlowDir(ndof);
    
	/* calculate flow direction at IP */
	gradActiveTensorFlowDir = FormGradActiveTensorFlowDirAtBandIP(ndof, bandIP);
	
	//cout << "gradActiveTensorFlowDir =\n" << gradActiveTensorFlowDir << endl;
	
	fD.SetToScaled(1.0, fCurrMaterial->ce_ijkl());
	dArrayT dGfD(fD.Rows());
	fD.MultTx(dGdSigma, dGfD);
  
	double jumpIncrement;
	double denom = gradActiveTensorFlowDir.Dot(dGfD,gradActiveTensorFlowDir);
	
	//cout << "workingTraction = " << workingTraction << endl;
	//cout << "denom = " << denom << endl; 
	double normalTraction = workingTraction.Dot(workingTraction, fBand -> Normal());
	double shearTraction = workingTraction.Dot(workingTraction, fBand -> PerpSlipDir(bandIP));
	
	//cout << "shearTraction = " << shearTraction << ", normalTraction = " << normalTraction << endl;
	//cout << "Resid cohesion = " << fBand -> ResidualCohesion(bandIP) << endl;
	
	/* assuming there is still cohesion softening softening */
	if (fBand -> ResidualCohesion(bandIP) > 0.0)
	{
  
		/* calculate jumpIncrement */
		jumpIncrement = shearTraction + fLocalizedFrictionCoeff * normalTraction
					- fBand -> ResidualCohesion(bandIP);
		//cout << "bandIP = " << bandIP << ", jumpIncrement = " << jumpIncrement << endl;
		//cout << "Resid cohesion = " << fBand -> ResidualCohesion(bandIP) << endl;
	
		jumpIncrement /= (denom + fH_delta_0);
	
		/* check to see if cohesion softening done */
		if (fBand -> ResidualCohesion(bandIP) > -1.0 * fH_delta_0 * jumpIncrement)
		{
			fBand->SetEffectiveSoftening(fH_delta_0, bandIP); 
			fBand -> StoreJumpIncrement(jumpIncrement, bandIP);
			//cout << "bandIP = " << bandIP << ", jumpIncrement = " << jumpIncrement << endl;
			return jumpIncrement;
		}
	}
   
	/* if softening done, recalculate jump increment*/
	/* todo - verify that this works for cohesion transition -seems to, but I'm skeptical */ 
   
	/* calculate jumpIncrement for the case of completed softening */
	jumpIncrement = shearTraction + fLocalizedFrictionCoeff * normalTraction;
	jumpIncrement /= denom;
   
	fBand->SetEffectiveSoftening(0.0, bandIP); 
	fBand -> StoreJumpIncrement(jumpIncrement, bandIP);
	cout << "Softening done, bandIP = " << bandIP << "jumpIncrement = " << jumpIncrement << endl;
	return jumpIncrement;  
}

bool SSEnhLocLinearT::IsBandActive(dArrayT workingTraction, int bandIP)
{
  //return true;

  double normalTraction = workingTraction.Dot(workingTraction, fBand -> Normal());
  double shearTraction = workingTraction.Dot(workingTraction, fBand -> PerpSlipDir(bandIP));  

  
  /* check for yielding */
  if (shearTraction < 0.0)
  {
	fBand -> FlipSlipDir(bandIP);
    shearTraction *= -1.0;
  }
  
  if (shearTraction > fLocalizedFrictionCoeff * -1.0 * normalTraction 
		+ fBand -> ResidualCohesion(bandIP))
  {
    fBand -> SetActive(bandIP, true);
	return true;
  }
  
  fBand -> SetActive(bandIP, false);
  return false;

}


void SSEnhLocLinearT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
               const iArrayT& e_codes, dArray2DT& e_values)
{

//cout << "ComputeOutput \n";
//SmallStrainT::ComputeOutput(n_codes, n_values, e_codes, e_values);
//return;

  cout << "ComputeOutput \n";


/* number of output values */
    int n_out = n_codes.Sum();
    int e_out = e_codes.Sum();
        
    int n_simo, n_extrap;
    if (qUseSimo) {
        n_simo = n_out - n_codes[iNodalDisp] - n_codes[iNodalCoord];
        n_extrap = n_codes[iNodalDisp] + n_codes[iNodalCoord];
    } else {
        n_simo = 0;
        n_extrap = n_out;
    }
            
    /* nothing to output */
    if (n_out == 0 && e_out == 0) return;

    /* dimensions */
    int nsd = NumSD();
    int ndof = NumDOF();
    int nen = NumElementNodes();
    int nnd = ElementSupport().NumNodes();
    int nstrs = fB.Rows();

    /* reset averaging workspace */
    ElementSupport().ResetAverage(n_extrap);

    /* allocate element results space */
    e_values.Dimension(NumElements(), e_out);

    /* nodal work arrays */
    dArray2DT nodal_space(nen, n_extrap);
    dArray2DT nodal_all(nen, n_extrap);
    dArray2DT coords, disp;
    dArray2DT nodalstress, princstress, matdat;
    dArray2DT energy, speed;
    dArray2DT Poynting;

    /* ip values */
    dSymMatrixT cauchy((nstrs != 4) ? nsd : dSymMatrixT::k3D_plane), nstr_tmp;
    dArrayT ipmat(n_codes[iMaterialData]), ipenergy(1);
    dArrayT ipspeed(nsd), ipprincipal(nsd), ipPoynting(nsd);

    /* set shallow copies */
    double* pall = nodal_space.Pointer();
    coords.Alias(nen, n_codes[iNodalCoord], pall); pall += coords.Length();
    disp.Alias(nen, n_codes[iNodalDisp], pall)   ; pall += disp.Length();

    /* work space for Poynting vector */
    dMatrixT stress_mat, F_inv, PbyJ;
    if (n_codes[iPoyntingVector]) {
        stress_mat.Dimension(nsd);
        F_inv.Dimension(nsd);
        PbyJ.Dimension(nsd);
    }

    /* workspaces for Simo */
    int simo_offset = coords.MinorDim() + disp.MinorDim();
    dArray2DT simo_space(nen,qUseSimo ? n_simo : 0);
    dArray2DT simo_all(nen,qUseSimo ? n_simo : 0);
    dArray2DT simoNa_bar(nen,qUseSimo ? 1 : 0);
    dArray2DT simo_force;
    dArray2DT simo_mass;
    iArrayT simo_counts;
        
    if (!qUseSimo) 
    {
        nodalstress.Alias(nen, n_codes[iNodalStress], pall); pall += nodalstress.Length();
        princstress.Alias(nen, n_codes[iPrincipal], pall)  ; pall += princstress.Length();
        energy.Alias(nen, n_codes[iEnergyDensity], pall)   ; pall += energy.Length();
        speed.Alias(nen, n_codes[iWaveSpeeds], pall)       ; pall += speed.Length();
        matdat.Alias(nen, n_codes[iMaterialData], pall)    ; pall += matdat.Length();
        Poynting.Alias(nen, n_codes[iPoyntingVector], pall);
    }
    else
    {
        simo_force.Dimension(ElementSupport().NumNodes(),qUseSimo ? n_simo : 0);
        simo_mass.Dimension(ElementSupport().NumNodes(),qUseSimo ? 1 : 0);
        simo_counts.Dimension(ElementSupport().NumNodes());

        pall = simo_space.Pointer();
        nodalstress.Alias(nen, n_codes[iNodalStress], pall); pall += nodalstress.Length();
        princstress.Alias(nen, n_codes[iPrincipal], pall)  ; pall += princstress.Length();
        energy.Alias(nen, n_codes[iEnergyDensity], pall)   ; pall += energy.Length();
        speed.Alias(nen, n_codes[iWaveSpeeds], pall)       ; pall += speed.Length();
        matdat.Alias(nen, n_codes[iMaterialData], pall)    ; pall += matdat.Length();
        Poynting.Alias(nen, n_codes[iPoyntingVector], pall);

        simo_mass = 0.;
        simo_force = 0.;
        simo_counts = 0;
    }
        
    /* element work arrays */
    dArrayT element_values(e_values.MinorDim());
    pall = element_values.Pointer();
    dArrayT centroid, ip_centroid, ip_mass;
    dArrayT ip_coords(nsd);
    if (e_codes[iCentroid])
    {
        centroid.Alias(nsd, pall); pall += nsd;
        ip_centroid.Dimension(nsd);
    }
    if (e_codes[iMass]) {
        ip_mass.Alias(NumIP(), pall); 
        pall += NumIP();
    }
    double w_tmp, ke_tmp;
    double mass;
    double& strain_energy = (e_codes[iStrainEnergy]) ? *pall++ : w_tmp;
    double& kinetic_energy = (e_codes[iKineticEnergy]) ? *pall++ : ke_tmp;
    dArrayT linear_momentum, ip_velocity;
    if (e_codes[iLinearMomentum])
    {
        linear_momentum.Alias(ndof, pall); pall += ndof;
        ip_velocity.Dimension(ndof);
    }
    else if (e_codes[iKineticEnergy]) ip_velocity.Dimension(ndof);

    dArray2DT ip_stress;
    if (e_codes[iIPStress])
    {
        ip_stress.Alias(NumIP(), e_codes[iIPStress]/NumIP(), pall);
        pall += ip_stress.Length();
    }
    dArray2DT ip_material_data;
    if (e_codes[iIPMaterialData])
    {
        ip_material_data.Alias(NumIP(), e_codes[iIPMaterialData]/NumIP(), pall);
        pall += ip_material_data.Length();
        ipmat.Dimension(ip_material_data.MinorDim());
    }

    /* check that degrees are displacements */
    int interpolant_DOF = InterpolantDOFs();

    bool is_axi = Axisymmetric();
    double Pi2 = 2.0*acos(-1.0);

    Top();
    while (NextElement())
        if (CurrentElement().Flag() != ElementCardT::kOFF)
        {
            /* initialize */
            nodal_space = 0.0;
            simo_space = 0.;
            simo_all = 0.;
            simoNa_bar = 0.;

            /* global shape function values */
            SetGlobalShape();
			
			/**********************************************************************/
	     
            /* output jump for element */
			jump_out.open_append("jump.info");
			if(IsElementTraced())
			    /* todo - currently only outputs first jump - update */
				jump_out << setw(16) << fBand->JumpIncrement(0) + fBand->Jump(0) << " ";
			else
			    jump_out << setw(16) << "0.0 ";	
	    


            /* collect nodal values */
            if (e_codes[iKineticEnergy] || e_codes[iLinearMomentum] || n_codes[iPoyntingVector]) {
                if (fLocVel.IsRegistered())
                    SetLocalU(fLocVel);
                else
                    fLocVel = 0.0;
            }

            /* coordinates and displacements all at once */
            if (n_codes[iNodalCoord]) fLocInitCoords.ReturnTranspose(coords);
            if (n_codes[ iNodalDisp]) {
                if (interpolant_DOF)
                    fLocDisp.ReturnTranspose(disp);
                else
                    NodalDOFs(CurrentElement().NodesX(), disp);
            }
 
            /* initialize element values */
            mass = strain_energy = kinetic_energy = 0;
            if (e_codes[iCentroid]) centroid = 0.0;
            if (e_codes[iLinearMomentum]) linear_momentum = 0.0;
            const double* j = fShapes->IPDets();
            const double* w = fShapes->IPWeights();

            /* integrate */
            dArray2DT Na_X_ip_w;
            fShapes->TopIP();
            while (fShapes->NextIP())
            {
                /* density may change with integration point */
                double density = fCurrMaterial->Density();

                /* element integration weight */
                double ip_w = (*j++)*(*w++);
                if (is_axi) {
                    fShapes->IPCoords(ip_coords);
                    ip_w *= Pi2*ip_coords[0]; /* radius is the x1 coordinate */
                }

                if (qUseSimo || qNoExtrap)
                {
                    Na_X_ip_w.Dimension(nen,1);
                    if (qUseSimo)
                    {
                        const double* Na_X = fShapes->IPShapeX();
                        Na_X_ip_w = ip_w;
                        for (int k = 0; k < nen; k++)
                            Na_X_ip_w(k,0) *= *Na_X++;
                        simoNa_bar += Na_X_ip_w;
                    }
                    else
                        for (int k = 0; k < nen; k++)
                            Na_X_ip_w(k,0) = 1.;
                }

                /* get Cauchy stress */
                //const dSymMatrixT stress(NumSD());                
                //const dSymMatrixT& stress = fCurrMaterial->s_ij();
                //dSymMatrixT stress(NumSD());
                
        // #if 1
                if (IsElementTraced())
          {
            dSymMatrixT strainIncr = fStrain_List [CurrIP()];
            strainIncr -= fStrain_last_List [CurrIP()];
            
            dSymMatrixT gradActiveTensorFlowDir =
              FormGradActiveTensorFlowDir(NumSD(), CurrIP());
            gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);
            
			/* todo - update this, currently uses only first jump */
            strainIncr.AddScaled(-1.0*fBand->JumpIncrement(0), gradActiveTensorFlowDir);
            dSymMatrixT stressIncr(NumSD());
            stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);
            stressIncr += fBand->Stress_List(CurrIP());
            //stress = stressIncr;
            //cout << "stressIncr = " << stressIncr << endl;
            cauchy.Translate(stressIncr);
          }
                else
          {
            const dSymMatrixT& stress = fCurrMaterial->s_ij();
            //cout << "stress = " << stress << endl;
            cauchy.Translate(stress);
            //cout << "cauchy = \n" << cauchy << endl;
          }    
                 
        //#endif
        //cout << "hi\n";    
        //cauchy.Translate(stress);
        //cout << "cauchy = \n" << cauchy << endl;


                /* stresses */
                if (n_codes[iNodalStress]) {        
                    if (qNoExtrap)
                        for (int k = 0; k < nen; k++)
                            nodalstress.AddToRowScaled(k,Na_X_ip_w(k,0),cauchy);
                    else
                        fShapes->Extrapolate(cauchy, nodalstress);
                }

                if (e_codes[iIPStress]) {
                    double* row = ip_stress(fShapes->CurrIP());
                    nstr_tmp.Set(nsd, row);
                    nstr_tmp = cauchy;
                    row += cauchy.Length();
                    nstr_tmp.Set(nsd, row);
                    fCurrMaterial->Strain(nstr_tmp);
                }

                /* wave speeds */
                if (n_codes[iWaveSpeeds])
                {
                    /* acoustic wave speeds */
                    fCurrMaterial->WaveSpeeds(fNormal, ipspeed);
                    if (qNoExtrap)
                        for (int k = 0; k < nen; k++)
                            speed.AddToRowScaled(k,Na_X_ip_w(k,0),ipspeed);
                        else
                            fShapes->Extrapolate(ipspeed, speed);
                }

                /* principal values - compute principal before smoothing */
                if (n_codes[iPrincipal])
                {
                    /* compute eigenvalues */
                    cauchy.PrincipalValues(ipprincipal);
                    if (qNoExtrap)
                        for (int k = 0; k < nen; k++)
                            princstress.AddToRowScaled(k,Na_X_ip_w(k,0),ipprincipal);
                    else
                        fShapes->Extrapolate(ipprincipal, princstress);        
                }

                /* strain energy density */
                if (n_codes[iEnergyDensity] || e_codes[iStrainEnergy] || n_codes[iPoyntingVector])
                {
                    double ip_strain_energy = fCurrMaterial->StrainEnergyDensity();

                    /* nodal average */
                    if (n_codes[iEnergyDensity])
                    {
                        ipenergy[0] = ip_strain_energy;
                        if (qNoExtrap)
                            for (int k = 0; k < nen; k++)
                                energy.AddToRowScaled(k,Na_X_ip_w(k,0),ipenergy);
                        else
                            fShapes->Extrapolate(ipenergy,energy);
                    }

                    /* integrate over element */
                    if (e_codes[iStrainEnergy])
                        strain_energy += ip_w*ip_strain_energy;
                        
                    /* Poynting vector */
                    if (n_codes[iPoyntingVector]) {
                        ipPoynting = 0.0;
                        ipPoynting[0] += fv_ss*ip_strain_energy;
                    }
                }

                /* material stuff */
                if (n_codes[iMaterialData] || e_codes[iIPMaterialData])
                {
                    /* compute material output */
                    fCurrMaterial->ComputeOutput(ipmat);

                    /* store nodal data */
                    if (n_codes[iMaterialData])
                    {
                        if (qNoExtrap)
                            for (int k = 0; k < nen; k++)
                                matdat.AddToRowScaled(k,Na_X_ip_w(k,0),ipmat);
                        else 
                            fShapes->Extrapolate(ipmat, matdat);
                    }

                    /* store element data */
                    if (e_codes[iIPMaterialData]) ip_material_data.SetRow(fShapes->CurrIP(), ipmat);
                }

                /* mass averaged centroid */
                if (e_codes[iCentroid] || e_codes[iMass])
                {
                    /* mass */
                    mass += ip_w*density;

                    /* integration point mass */
                    if (e_codes[iMass]) ip_mass[fShapes->CurrIP()] = ip_w*density;

                    /* moment */
                    if (e_codes[iCentroid]) {
                        fShapes->IPCoords(ip_centroid);
                        centroid.AddScaled(ip_w*density, ip_centroid);
                    }
                }

                /* kinetic energy/linear momentum */
                if (e_codes[iKineticEnergy] || e_codes[iLinearMomentum] || n_codes[iPoyntingVector])
                {
                    /* velocity at integration point */
                    fShapes->InterpolateU(fLocVel, ip_velocity);                                        
                    double ke_density = 0.5*density*dArrayT::Dot(ip_velocity, ip_velocity);

                    /* kinetic energy */
                    if (e_codes[iKineticEnergy])
                        kinetic_energy += ip_w*ke_density;
        
                    /* linear momentum */
                    if (e_codes[iLinearMomentum])
                        linear_momentum.AddScaled(ip_w*density, ip_velocity);
                
                    /* Poynting vector */
                    if (n_codes[iPoyntingVector]) {
                        ipPoynting[0] += fv_ss*ke_density;

                        /* is finite strain */
                        
                        #if 0
                        FSSolidMatT* fs_mat = TB_DYNAMIC_CAST(FSSolidMatT*, fCurrMaterial);
                        if (fs_mat)
                        {
                            /* compute PK1/J */
                            cauchy.ToMatrix(stress_mat);
                            F_inv.Inverse(fs_mat->F());
                            PbyJ.MultABT(stress_mat, F_inv);
                        
                            PbyJ.Multx(ip_velocity, ipPoynting, 1.0/F_inv.Det(), dMatrixT::kAccumulate);
                        } 
                        else { /* small strain */
                        
                        #endif
                        
                            cauchy.Multx(ip_velocity, ipPoynting, 1.0, dMatrixT::kAccumulate);
                        
                        //}
                        
                        /* extrapolate */
                        ipPoynting *= -1; /* positive toward crack tip */
                        fShapes->Extrapolate(ipPoynting, Poynting);
                    }
                }
            }

            /* copy in the cols */
            int colcount = 0;
            nodal_all.BlockColumnCopyAt(disp       , colcount); colcount += disp.MinorDim();
            nodal_all.BlockColumnCopyAt(coords     , colcount); colcount += coords.MinorDim();

            if (!qUseSimo)
            {
                if (qNoExtrap)
                {
                    double nip(fShapes->NumIP());
                    nodalstress /= nip;
                    princstress /= nip;
                    energy /= nip;
                    speed /= nip;
                    matdat /= nip;
                }
                nodal_all.BlockColumnCopyAt(nodalstress, colcount); colcount += nodalstress.MinorDim();
                nodal_all.BlockColumnCopyAt(princstress, colcount); colcount += princstress.MinorDim();
                nodal_all.BlockColumnCopyAt(energy     , colcount); colcount += energy.MinorDim();
                nodal_all.BlockColumnCopyAt(speed      , colcount); colcount += speed.MinorDim();
                nodal_all.BlockColumnCopyAt(matdat     , colcount); colcount += matdat.MinorDim();
                nodal_all.BlockColumnCopyAt(Poynting   , colcount); colcount += Poynting.MinorDim();
            }
            else
            {        
                colcount = 0;
                simo_all.BlockColumnCopyAt(nodalstress, colcount); colcount += nodalstress.MinorDim();
                simo_all.BlockColumnCopyAt(princstress, colcount); colcount += princstress.MinorDim();
                simo_all.BlockColumnCopyAt(energy     , colcount); colcount += energy.MinorDim();
                simo_all.BlockColumnCopyAt(speed      , colcount); colcount += speed.MinorDim();
                simo_all.BlockColumnCopyAt(matdat     , colcount); colcount += matdat.MinorDim();
                simo_all.BlockColumnCopyAt(Poynting   , colcount); colcount += Poynting.MinorDim();

                iArrayT currIndices = CurrentElement().NodesX();
                simo_force.Accumulate(currIndices,simo_all);
                simo_mass.Accumulate(currIndices,simoNa_bar);
                for (int i = 0; i < currIndices.Length(); i++)
                    simo_counts[currIndices[i]]++;
            }

            /* accumulate - extrapolation done from ip's to corners => X nodes */
            ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_all);

            /* element values */
            if (e_codes[iCentroid]) centroid /= mass;

            /* store results */
            e_values.SetRow(CurrElementNumber(), element_values);
        }

    /* get nodally averaged values */
    const OutputSetT& output_set = ElementSupport().OutputSet(fOutputID);
    const iArrayT& nodes_used = output_set.NodesUsed();
    dArray2DT extrap_values(nodes_used.Length(), n_extrap);
    extrap_values.RowCollect(nodes_used, ElementSupport().OutputAverage());

    int tmpDim = extrap_values.MajorDim();
    n_values.Dimension(tmpDim,n_out);
    n_values.BlockColumnCopyAt(extrap_values,0);
    if (qUseSimo)
    {        
        int rowNum = 0;
//      iArrayT nodes_used(tmpDim);
        dArray2DT tmp_simo(tmpDim, n_simo);
        for (int i = 0; i < simo_force.MajorDim(); i++)
            if (simo_counts[i] > 0) {
//              nodes_used[rowNum] = i;
                simo_force.ScaleRow(i, 1./simo_mass(i,0));
                tmp_simo.SetRow(rowNum, simo_force(i));
                rowNum++;
            }

        /* collect final values */
        n_values.BlockColumnCopyAt(tmp_simo, simo_offset);

        /* write final values back into the averaging workspace */
        if (extrap_values.MinorDim() != n_values.MinorDim()) {
            ElementSupport().ResetAverage(n_values.MinorDim());
            ElementSupport().AssembleAverage(nodes_used, n_values);
        }
    }

  /* inherited */
  //SmallStrainT::ComputeOutput(n_codes, n_values, e_codes, e_values);
  
  jump_out << endl;  

}


void SSEnhLocLinearT::CloseStep(void)
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
				/* Conforming Strain Increment */
				dSymMatrixT strainIncr = fStrain_List [CurrIP()];
				strainIncr -= fStrain_last_List [CurrIP()];
				
				//cout << "strainIncr =\n" << strainIncr; 

				/* modify strain increment for avg jump Increment */
				dSymMatrixT gradActiveTensorFlowDir =
				FormGradActiveTensorFlowDir(NumSD(), CurrIP());
				gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);

				//double avgJumpIncr = 0.0;
				//for (int i = 0; i < fBand -> NumSurfaceIPs(); i ++)
				//	avgJumpIncr += fBand -> JumpIncrement(i);
				//avgJumpIncr /= fBand -> NumSurfaceIPs();	
		
				strainIncr.AddScaled(-1.0*JumpIncrAtBulkIP(CurrIP()), gradActiveTensorFlowDir);
		
				/* modify strain increment for jump increment gradient */
				dSymMatrixT perpSlipDirOuter(NumSD());
				perpSlipDirOuter.Outer(fBand -> PerpSlipDir(0));
		
		
				/*
				double linearModifier;
				if (fBand -> IsBulkIPActive(CurrIP()))
					linearModifier = 1.0; 
				else
					linearModifier = 0.0;
				
				linearModifier -= F_hAtBulkIP(CurrIP());
				linearModifier /= fBand -> DistanceBetweenIPs();
				linearModifier *= (fBand -> JumpIncrement(1) - fBand -> JumpIncrement(0));
				*/
		
				dMatrixT ce_inverse = fCurrMaterial->ce_ijkl();
				ce_inverse.Inverse();
				dArrayT perpSlipDir = fBand -> PerpSlipDir(0); 
				double denom = QuadContraction(ce_inverse, perpSlipDir, perpSlipDir, perpSlipDir, perpSlipDir);

				double linearModifier;
				if (fBand -> IsBulkIPActive(CurrIP()))
					linearModifier = 1.0; 
				else
					linearModifier = 0.0;
				
				linearModifier -= F_hAtBulkIP(CurrIP());		
				linearModifier = (fBand -> JumpIncrement(1) - fBand -> JumpIncrement(0))/
								(fBand -> DistanceBetweenIPs() * denom);
		
				//strainIncr.AddScaled(-1.0*linearModifier, perpSlipDirOuter);		
	      
				//cout << "strainIncr =\n" << strainIncr;
		  
				dSymMatrixT stressIncr(NumSD());
				stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);
				stressIncr.AddScaled(linearModifier, perpSlipDirOuter);
				fBand -> IncrementStress(stressIncr, CurrIP());
				//cout << "stressIncr = \n" << stressIncr << endl;
	      	}
			
			/* update band traction */
			for (int bandIP = 0; bandIP < fBand -> NumSurfaceIPs(); bandIP ++)
			{
				dSymMatrixT strainIncr = ConformingStrainIncrAtCoord(fBand -> SurfaceIPCoords(bandIP));
				
				//cout << "strainIncr =\n" << strainIncr;  
				
				dSymMatrixT gradActiveTensorFlowDir =  
					FormGradActiveTensorFlowDirAtBandIP(NumSD(), bandIP);
				gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);
				strainIncr.AddScaled(-1.0 * fBand -> JumpIncrement(bandIP), gradActiveTensorFlowDir);
				 
				//cout << "strainIncr =\n" << strainIncr; 
				 
				dSymMatrixT stressIncr(NumSD());
				stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);
				
				dArrayT tractionIncr(NumSD());
				stressIncr.Multx(fBand -> Normal(), tractionIncr);
				//cout << "tractionIncr = \n" << tractionIncr << endl;
				
				fBand -> IncrementTractionAtBandIP(tractionIncr, bandIP);				
			}
			
			fBand -> CloseStep();
	      	
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
	  //cout << "hi\n";
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
				//if (fBVPType == kHomogeneous && !fLocalizationHasBegun)
				//	fLeastDetEle = fFirstElementToLocalize;
			
				fLocalizationHasBegun = true;
	  
	  			//cout << "fLeastDetEle = " << fLeastDetEle << endl;
				GetElement(fLeastDetEle);
				//fEdgeOfBandCoords.Free();
				fEdgeOfBandElements.Append(fLeastDetEle);
				fEdgeOfBandCoords.Append(Centroid());	  
				//fEdgeOfBandElements.Top();
				//fEdgeOfBandCoords.Top();
			}
    }
}

 //cout << "hi\n" << flush; 



 /* uncomment these two lines to restrict band propagation after initial
    localization */
if (fLocalizationHasBegun)
	fSeedElementsSet = true;

  SmallStrainT::CloseStep();
 }


/***********************************************************************
 * Protected
 ***********************************************************************/

/* current element operations */
void SSEnhLocLinearT::GetElement(int elementNumber)
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




//move to surface mat model?
dSymMatrixT SSEnhLocLinearT::FormdGdSigma(int ndof)
{
  dSymMatrixT dGdSigma(ndof), work(ndof);
  dMatrixT dGNonSym(ndof);

//todo - figure out how to choose which perpSLipDir
  dGNonSym.Outer(fBand->Normal(), fBand->PerpSlipDir(0));
  dGdSigma.Symmetrize(dGNonSym);
  work.Outer(fBand->Normal());
  dGdSigma.AddScaled(fLocalizedFrictionCoeff, work);

  /* use vector values for mult */
  dGdSigma.ScaleOffDiagonal(2.0);

  return dGdSigma;
}

dSymMatrixT SSEnhLocLinearT::FormGradActiveTensorFlowDir(int ndof, int ip)
{
  dSymMatrixT G(ndof);
  dMatrixT G_NonSym(ndof);
  iAutoArrayT activeNodes = fBand->ActiveNodes();
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
	  //cout << "active node " << A << endl;
	  grad_f[i] += fB(i, (A)*ndof +i);
	}
    }

  /* todo - update this to take whichever slipdir */
  G_NonSym.Outer(grad_f, fBand->SlipDir(0));
  G.Symmetrize(G_NonSym);

  /* use vector values for mult */
  G.ScaleOffDiagonal(2.0);

  return G;
}

dSymMatrixT SSEnhLocLinearT::FormGradActiveTensorFlowDirAtBandIP(int ndof, int bandIP)
{
  dSymMatrixT G(ndof);
  dMatrixT G_NonSym(ndof);
  iAutoArrayT activeNodes = fBand->ActiveNodes();
  int A;
  dArrayT grad_f(ndof);
  grad_f = 0.0;


  
	/* transform coordinates to Parent Domain (Local Coordinates)*/
	dArrayT localCoords(NumSD());
	/*
	if (!fShapes -> ParentDomain().MapToParentDomain(InitialCoordinates(), fBand -> SurfaceIPCoords(bandIP) , localCoords))
	{
		cout << "SSEnhLocLinearT::FormGradACtiveTensorFlowDirAtBandIP, failed to map to local coordinates\n " << flush;
		throw ExceptionT::kGeneralFail;
	}*/
	
	fShapes -> ParentDomain().MapToParentDomain(InitialCoordinates(), fBand -> SurfaceIPCoords(bandIP) , localCoords);
	

  /* Get Local-to-Global transformation */
  dArrayT shapeFunctions;
  dArray2DT shapeFunctionGradients, globalShapeFunctionGradients;
  dMatrixT jacobian(NumSD());
  
  fShapes -> GradU(InitialCoordinates(), jacobian, 
	localCoords, shapeFunctions, shapeFunctionGradients); 

  //cout << "jacobian =\n" << jacobian << endl;
	
  jacobian.Inverse();

  /* Calculate Global Derivatives */
  
  //cout << "jacobian =\n" << jacobian << endl;
  //cout << "shapeFunctionGradients = \n" << shapeFunctionGradients << endl;
	  
  fShapes -> TransformDerivatives(jacobian, shapeFunctionGradients, globalShapeFunctionGradients);
  //calls for dArray2DTs for 2nd 2 args */

  //cout << "globalShapeFunctionGradients = \n" << globalShapeFunctionGradients << endl;

  for (int i=0; i<ndof; i++)
    {
      activeNodes.Top();

      while(activeNodes.Next())      
	{
	  A = activeNodes.Current();  
	  grad_f[i] += globalShapeFunctionGradients(i,A);//think this is correct order, see ShapeFunctionT:GradU
	}
    }

  G_NonSym.Outer(grad_f, fBand->SlipDir(bandIP));
  G.Symmetrize(G_NonSym);

  /* use vector values for mult */
  G.ScaleOffDiagonal(2.0);

  return G;
}



bool SSEnhLocLinearT::IsElementTraced()
{
 	 //int elementNumber = CurrElementNumber();
 	 return IsElementTraced(CurrElementNumber());
}
  
bool SSEnhLocLinearT::IsElementTraced(int elementNumber)
{
  bool isTraced = fTracedElements.HasKey(elementNumber);

  if (isTraced)
    LoadBand(elementNumber);

  return isTraced;
}

void SSEnhLocLinearT::LoadBand(int elementNumber)
{
  fBand = fTracedElements[elementNumber];
}

bool SSEnhLocLinearT::IsElementLocalized()
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
      //cout << "hi \n";

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

bool SSEnhLocLinearT::TraceElement()
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

void SSEnhLocLinearT::PreFailElements()
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

void SSEnhLocLinearT::ChooseNormals(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipDirs)
{
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
	  
	  //temp -overwrite - todo - delete
	  //prod = normals.Current() [0] * 7.646947e-01 - normals.Current() [0] * 6.443928e-01;
	  
      if (prod > maxProd)
	{
	  //cout << "best normal = \n" << normals.Current() << endl; 
	  normal = normals.Current(); 
	  slipDir = slipDirs.Current();
          maxProd = prod;
	}
    }

 //1 For propagating band normal in same direction
 #if 1
   LinearBandT* fBandTemp = fBand;
 
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
     cout << "SSEnhLocLinearT::AddNewEdgeElements, geometry not implemented. \n" << flush;
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
			normal = fBand -> Normal();
         }
 
 
  fBand = fBandTemp;
 #endif



  //make sure slip direction is dilatant
  if (normal.Dot(normal, slipDir)<0.0)
    slipDir *= -1.0;

  //normal *= -1.0;
  //slipDir *= -1.0;

  perpSlipDir = slipDir;
  perpSlipDir.AddScaled(-1.0*slipDir.Dot(slipDir, normal), normal);
  perpSlipDir.UnitVector();

  if (fNoBandDilation)
    {
      slipDir = perpSlipDir;
    }

  //very temp
  /*
  normal[0] = .8;
  normal[1] = .6;
  slipDir[0] = .6;
  slipDir[1] = -.8;
  perpSlipDir = slipDir;
	*/

  fBand = FormNewBand(normal, slipDir, perpSlipDir, fEdgeOfBandCoords.Current(), area);
  
  //cout << "fBand->PerpSlipDir() = " << fBand->PerpSlipDir() << endl;
  //cout << "fBand->Coords() = " << fBand->Coords() << endl;
  //cout << "fBand = " << fBand << endl;
  //cout << "1 " << flush;

  fTracedElements.Insert(CurrElementNumber(), fBand);

  //cout << "fBand = " << fBand << endl;
  //cout << "2 " << flush;

  AddNewEdgeElements(CurrElementNumber());
  
  //cout << "3 " << flush;
}

LinearBandT* SSEnhLocLinearT::FormNewBand(dArrayT normal, dArrayT slipDir,
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

  double residCohesion = shearStress + normalStress * fLocalizedFrictionCoeff;

  if (fBVPType == kPreFailed && ElementSupport().Time() == 0.0)
    residCohesion = 0.0;

  // todo remove stress average and calculation of resid cohesion
  cout << "residCohesion = " << residCohesion << endl;

return new LinearBandT(normal, slipDir, perpSlipDir, coords, fH_delta_0, stressList, this);
}

void SSEnhLocLinearT::AddNewEdgeElements(int elementNumber)
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
  
 
  
  iAutoArrayT activeNodes = fBand->ActiveNodes();
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
     cout << "SSEnhLocLinearT::AddNewEdgeElements, geometry not implemented. \n" << flush;
     throw ExceptionT::kGeneralFail;
   }
}


  LocalArrayT nodalCoords = InitialCoordinates();
  dArrayT nodalCoord1(NumSD()), nodalCoord2(NumSD()); //coords a particular node

  //cout << "fBand = " << fBand << endl;
  
  //temporary workaround- IsElementTraced loads the element, changing the coordinates
  LinearBandT* fBandTemp = fBand;

  for(int i = 0; i < numSides; i++)
    //if (((fBand->ActiveNodes()).HasValue((i+1) % numSides) && !((fBand->ActiveNodes()).HasValue(i)))
	//|| (!((fBand->ActiveNodes()).HasValue((i+1) % numSides)) && (fBand->ActiveNodes()).HasValue(i)))
	if ((activeNodes.HasValue((i+1) % numSides) && !(activeNodes.HasValue(i)))
	|| (!(activeNodes.HasValue((i+1) % numSides)) && activeNodes.HasValue(i)))
      {

	  //cout << "fBand = " << fBand << endl;	
      fBand = fBandTemp;
	  
	
		if (!(neighbors(elementNumber,i) == -1 || IsElementTraced(neighbors(elementNumber ,i))))
	  	{
	    	//get coords
			fBand = fBandTemp;
			  //cout << "fBand = " << fBand << endl;
	    	dArrayT localizedEleCoord = fBand -> Coords();

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

dArrayT SSEnhLocLinearT::InterceptCoords(dArrayT& localizedEleCoord,
dArrayT& nodalCoord1, dArrayT& nodalCoord2)
{
  //assumes straight sides
  dArrayT sideVector = nodalCoord2;
  sideVector -= nodalCoord1;

  dArrayT perpSlipDir = fBand -> PerpSlipDir(0);//at this stage doesn't matter which

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


dArrayT SSEnhLocLinearT::Centroid()
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

/*----------------------------------------------------------------------------------*/
/*------------------ Strain Functions for Surface Gauss Points ---------------------*/
/*----------------------------------------------------------------------------------*/

dSymMatrixT SSEnhLocLinearT::ConformingStrainAtCoord(dArrayT coord)
{
	/* transform coordinates to Parent Domain (Local Coordinates)*/
	dArrayT localCoords(NumSD());
	/*
	if (!(fShapes -> ParentDomain().MapToParentDomain(InitialCoordinates(),coord, localCoords)))
	{
		cout << "SSEnhLocLinearT::StrainAtCoord, failed to map to local coordinates\n " << flush;
		throw ExceptionT::kGeneralFail;
	}*/
	
	fShapes -> ParentDomain().MapToParentDomain(InitialCoordinates(),coord, localCoords);
	
	/* Calculate Strain at Coordinates */
	dMatrixT grad_U_local(NumSD()), grad_U_global(NumSD()); 
	dArray2DT DNa; //DNa - shape function gradients wrt local coords
    dArrayT Na;	// Na - shape functions
	
	/* Displacements() in ContinuumElementT */
	fShapes -> GradU(Displacements(), grad_U_local, localCoords, Na, DNa); 
	
	/* Get Local-to-Global transformation */
	dMatrixT jacobian(NumSD());
	fShapes -> GradU(InitialCoordinates(), jacobian, localCoords, Na, DNa); 
	
	jacobian.Inverse();
	
	//fShapes -> TransformDerivatives(jacobian, grad_U_local, grad_U_global);
	grad_U_global.MultAB(grad_U_local, jacobian);
	
	/* Symmetrize */
	dSymMatrixT strain(NumSD());
	strain.Symmetrize(grad_U_global);
	
	return strain;
}

dSymMatrixT SSEnhLocLinearT::ConformingStrainIncrAtCoord(dArrayT coord)
{
    dSymMatrixT incr = ConformingStrainAtCoord(coord);

	/* transform coordinates to Parent Domain (Local Coordinates)*/
	dArrayT localCoords(NumSD());
	/*
	if (!fShapes -> ParentDomain().MapToParentDomain(InitialCoordinates(),coord, localCoords))
	{
		cout << "SSEnhLocLinearT::StrainAtCoord, failed to map to local coordinates\n " << flush;
		throw ExceptionT::kGeneralFail;
	}*/
	
	fShapes -> ParentDomain().MapToParentDomain(InitialCoordinates(),coord, localCoords);
	
	/* Calculate Strain at Coordinates */
	dMatrixT grad_U_local(NumSD()), grad_U_global(NumSD()); 
	dArray2DT DNa; //DNa - shape function gradients wrt local coords
    dArrayT Na;	// Na - shape functions
	
	
	/* LastDisplacements() in SolidElementT */
	fShapes -> GradU(LastDisplacements(), grad_U_local, localCoords, Na, DNa); 
	
	/* Get Local-to-Global transformation */
	dMatrixT jacobian(NumSD());
	fShapes -> GradU(InitialCoordinates(), jacobian, localCoords, Na, DNa); 
	
	jacobian.Inverse();
	
	grad_U_global.MultAB(grad_U_local, jacobian);
	//fShapes -> TransformDerivatives(jacobian, grad_U_local, grad_U_global);
	
	/* Symmetrize */
	dSymMatrixT strainLast(NumSD());
	strainLast.Symmetrize(grad_U_global);
	
	incr -= strainLast;
	
	return incr;
}

/* usually evaluated at bulk ip coordinates */
double SSEnhLocLinearT::F_hAtBulkIP(int ip)
{
	double f_h = 0.0;
	const double *shapeFunctions = fShapes -> IPShapeX(ip);
	iAutoArrayT activeNodes = fBand -> ActiveNodes();
	
	/* loop over active nodes*/
	activeNodes.Top();
    //fBand -> ActiveNodes().Top();

	while(activeNodes.Next())      
		f_h += shapeFunctions [activeNodes.Current()]; 
	//while(fBand -> ActiveNodes().Next())      
	//	f_h += shapeFunctions [fBand -> ActiveNodes().Current()];
	
	
	return f_h;
}

dMatrixT SSEnhLocLinearT::B_at_surfaceIP(int bandIP)
{
	/* transform coordinates to Parent Domain (Local Coordinates)*/
	dArrayT localCoords(NumSD());
	dArrayT coord = fBand -> SurfaceIPCoords(bandIP);
	
	/*
	if (!fShapes -> ParentDomain().MapToParentDomain(InitialCoordinates(),coord, localCoords))
	{
		cout << "SSEnhLocLinearT::StrainAtCoord, failed to map to local coordinates\n " << flush;
		throw ExceptionT::kGeneralFail;
	}*/
	
	fShapes -> ParentDomain().MapToParentDomain(InitialCoordinates(),coord, localCoords);
	
	/* Get Local-to-Global transformation */
  dArrayT shapeFunctions;
  dArray2DT shapeFunctionGradients, globalShapeFunctionGradients;
  dMatrixT jacobian(NumSD());
  
  fShapes -> GradU(InitialCoordinates(), jacobian, 
	localCoords, shapeFunctions, shapeFunctionGradients); 
	
  jacobian.Inverse();

  /* Calculate Global Derivatives */
  fShapes -> TransformDerivatives(jacobian, shapeFunctionGradients, globalShapeFunctionGradients);
  


  dMatrixT bB;
  int nnd = shapeFunctionGradients.MinorDim();
  //double* pB = bB.Pointer();
		
  /* arrange into B matrix - follows SolidElementT*/
  if (shapeFunctionGradients.MajorDim() == 2)
  {
		 bB.Dimension(3, nnd*2);
		 double* pB = bB.Pointer();
         const double* pNax = globalShapeFunctionGradients(0);
         const double* pNay = globalShapeFunctionGradients(1);
         for (int i = 0; i < nnd; i++)
         {             /* see Hughes (2.8.20) */
             *pB++ = *pNax;
             *pB++ = 0.0;
             *pB++ = *pNay;
 
             *pB++ = 0.0;
             *pB++ = *pNay++;
             *pB++ = *pNax++;
         }
   }	 
   else
   {
     cout << "SSEnhLocLinearT::B_at_surfaceIP, only 2 dimensions "<< endl << flush;
	 throw ExceptionT::kGeneralFail;
   }
   
  return bB;
}

double SSEnhLocLinearT::JumpIncrAtBulkIP(int ip)
{
	double coord = fBand -> IPBandCoord(ip);
	double distance = fBand -> DistanceBetweenIPs();
	
	return (0.5 - coord/distance) * fBand -> JumpIncrement(0)
		+ (0.5 + coord/distance) * fBand -> JumpIncrement(1);
}

double SSEnhLocLinearT::QuadContraction(dMatrixT c_ijkl, 
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