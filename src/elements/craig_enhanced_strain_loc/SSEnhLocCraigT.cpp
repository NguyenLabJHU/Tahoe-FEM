/* $Id: SSEnhLocCraigT.cpp,v 1.10 2005-03-17 21:35:54 cfoster Exp $ */
#include "SSEnhLocCraigT.h"
#include "ShapeFunctionT.h"
#include "SSSolidMatT.h"

/* materials lists */
#include "SSSolidMatList1DT.h"
#include "SSSolidMatList2DT.h"
#include "SSSolidMatList3DT.h"
#include "SSMatSupportT.h"
#include "ParameterContainerT.h"
#include "ModelManagerT.h"

#include "DetCheckT.h"


using namespace Tahoe;

/* constructor */
SSEnhLocCraigT::SSEnhLocCraigT(const ElementSupportT& support):
	SmallStrainT(support),
	fBand(NULL)
{
	SmallStrainT::SetName("small_strain_enh_loc_craig");
}

/* destructor */
/*
SSEnhLocCraigT::~SSEnhLocCraigT(void)
{
	delete fSSMatSupport;
}
*/

/* implementation of the ParameterInterfaceT interface */
void SSEnhLocCraigT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SmallStrainT::DefineParameters(list);


	/*PARAMETERS FOR ENHANCED STRAIN*/
	list.AddParameter(fH_delta_0, "Post-Localization_softening_parameter_H_Delta"); 
	list.AddParameter(fNoBandDilation, "Disallow_Dilation_on_Band");
	list.AddParameter(fLocalizedFrictionCoeff, "Localized_Friction_Coefficient");
}

/* information about subordinate parameter lists */
void SSEnhLocCraigT::DefineSubs(SubListT& sub_list) const
{	
	/* inherited */
  //SolidElementT::DefineSubs(sub_list);
  SmallStrainT::DefineSubs(sub_list);

	/* element block/material specification */
	//sub_list.AddSub("small_strain_enh_loc_craig_element_block", ParameterListT::OnePlus);
}

void SSEnhLocCraigT::TakeParameterList(const ParameterListT& list)
{
  //SmallStrainT::TakeParameterList(list);

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
}

/* extract the list of material parameters */
void SSEnhLocCraigT::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{

  SmallStrainT::CollectMaterialInfo(all_params, mat_params);

}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct a new material support and return a pointer */
MaterialSupportT* SSEnhLocCraigT::NewMaterialSupport(MaterialSupportT* p) const
{
  return SmallStrainT::NewMaterialSupport(p);
}


/* return a pointer to a new material list */
MaterialListT* SSEnhLocCraigT::NewMaterialList(const StringT& name, int size)
{
  return SmallStrainT::NewMaterialList(name, size);
}


/* calculate the internal force contribution ("-k*d") */
void SSEnhLocCraigT::FormKd(double constK)
{


  /*
  if (fInitialModulus == 0.0)
    {
      fInitialModulus.Dimension(NumDOF());
      fInitialModulus = fCurrMaterial->c_ijkl();
    }
  */

  if(!IsElementTraced())
    SmallStrainT::FormKd(constK);
  else 
      {
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

	    /* B^T * Cauchy stress */
	    dSymMatrixT strainIncr = fStrain_List [CurrIP()];
	    strainIncr -= fStrain_last_List [CurrIP()];
	    //	cout << "strainIncr =\n" << strainIncr << endl;
	    dSymMatrixT stressIncr(NumSD());
	    stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);
	    stressIncr += fBand->Stress_List(CurrIP());
	    fB.MultTx(stressIncr, fNEEvec);
	    
	    /* accumulate */
	    fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);
	    
	    /* incremental heat generation */
	    if (need_heat) 
	      fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	  }     
      }
}

/* form the element stiffness matrix */
void SSEnhLocCraigT::FormStiffness(double constK)
{


  /*
  double ck_area = 0.0;
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

  fShapes->TopIP();
  while (fShapes->NextIP())
    ck_area += (*Det++)*(*Weight++);



  cout << "Element # " << CurrElementNumber() << "has area " << ck_area <<
  endl;
  */

  if (fStrainDispOpt == kMeanDilBbar)
    cout << "Warning SSEnhLocCraigT::FormStiffness, b-bar integration not implemented for enhanced strain element, postlocalization results may be garbage.";

  if (!IsElementTraced() || !fBand->IsActive() )
    {
    /* form stiffness in standard way */
    SmallStrainT::FormStiffness(constK);

    }
  else //if already localized, use localized stiffness routine
    {

      //cout << "constK =\n" << constK << endl;

	/* matrix format */
	dMatrixT::SymmetryFlagT format = dMatrixT::kWhole;

	/* integrate element stiffness */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	int ndof = NumDOF();
	int nen = NumElementNodes();
	int nedof = nen*ndof;//# of element dof
	double k_zeta_zeta = 0.0, area = 0.0;
	dArrayT k_d_zeta(nedof), k_zeta_d(nedof);
	dArrayT k_d_zeta_work(nedof), k_zeta_d_work(nedof);
	k_d_zeta = 0.0;
	k_zeta_d = 0.0;

	dSymMatrixT gradActiveTensorFlowDir(ndof), dGdSigma(ndof);
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
                //fLHS.MultATB(fB, fD, format, dMatrixT::kAccumulate);	

		//form k_d_zeta
		gradActiveTensorFlowDir = FormGradActiveTensorFlowDir(ndof);

		//k_d_zeta_work = 0.0;
		//fDfB.MultTx(gradActiveTensorFlowDir, k_d_zeta_work,::dMatrixT::kOverwrite);
		fDfB.MultTx(gradActiveTensorFlowDir, k_d_zeta_work);
		k_d_zeta += k_d_zeta_work;
		//cout << "k_d_zeta =\n" << k_d_zeta << endl;

		//form k_zeta_d
		//fDfB.MultTx(dGdSigma, k_zeta_d_work, dMatrixT::kOverwrite);
		fDfB.MultTx(dGdSigma, k_zeta_d_work);
		k_zeta_d += k_zeta_d_work;
		//cout << "k_zeta_d =\n" << k_zeta_d << endl;

		//form k_zeta_zeta
		k_zeta_zeta +=
		fD.MultmBn(dGdSigma,gradActiveTensorFlowDir);
	}

	k_d_zeta *= 1.0/area;

	k_zeta_zeta *= 1.0/area;
	k_zeta_zeta += fBand->EffectiveSoftening();
	/*k_zeta_zeta *= (k_zeta_zeta + fBand->H_delta())/(k_zeta_zeta + (1
									-
									fBand-> EffectiveSoftening()) * fBand->H_delta());*/

	//k_zeta_zeta *= -1.0;

	fLHS.Outer(k_d_zeta, k_zeta_d, -1.0/k_zeta_zeta, dMatrixT::kAccumulate);
	  }
}

/* compute the measures of strain/deformation over the element */
void SSEnhLocCraigT::SetGlobalShape(void)
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

	    double jumpIncrement = CalculateJumpIncrement();

	    /* loop over integration points again */
	    for (int i = 0; i < NumIP(); i++)
	      {
		gradActiveTensorFlowDir = FormGradActiveTensorFlowDir(ndof);
		
		/*change shear strains back to matrix values */
		/* vector values are used when created */
		gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);
		
		
		/* deformation gradient */
		if (needs[fNeedsOffset + kstrain])
		  {			  
		    // if (i==0)
		    // cout << "fStrain_List[" << i << "] = " <<
		    //	fStrain_List[i] << endl;
		    
		    fStrain_List[i].AddScaled(-(fBand->Jump()
						+ jumpIncrement),
					      gradActiveTensorFlowDir);

		    //if (i==0)
		    //  cout << "fStrain_List[" << i << "] = " <<
		    //	fStrain_List[i] << endl;
		  }
		
		/* "last" deformation gradient */ //is this right?
		if (needs[fNeedsOffset + kstrain_last])
		  {
		    
		    //if (i==0)
		    //    cout << "fStrain_last_List[" << i << "] = " <<
		    //	fStrain_last_List[i] << endl;
		    
		    fStrain_last_List[i].AddScaled(-(fBand->Jump()),
						   gradActiveTensorFlowDir);
		    
		    //if (i==0)
		    //cout << "fStrain_last_List[" << i << "] = " <<
		    //	fStrain_last_List[i] << endl;
		    
		  }
	      }
	  }
    }
}

/***********************************************************************
 * Protected
 ***********************************************************************/



double SSEnhLocCraigT::CalculateJumpIncrement()
{
  if (!IsBandActive())
    return 0.0;
  
  int ndof = NumDOF();
  dSymMatrixT dGdSigma = FormdGdSigma(ndof);
  dSymMatrixT gradActiveTensorFlowDir(ndof);
    
  //fD.SetToScaled(1.0, HookeanMatT::Modulus());
  fD.SetToScaled(1.0, fCurrMaterial->ce_ijkl());
  dArrayT dGfD(fD.Rows());
  fD.MultTx(dGdSigma, dGfD);
	    
  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();	    
  double area = 0.0;
  double jumpIncrement = 0.0;
  double jumpWork = 0.0;

  /* loop over integration points */
  for (int i = 0; i < NumIP(); i++)
    {
      //area += (*Det)*(*Weight);
      double scale = (*Det++)*(*Weight++);
      area += scale;
      dSymMatrixT strainIncr = fStrain_List [i];
      strainIncr -= fStrain_last_List [i];
      strainIncr.ScaleOffDiagonal(2.0);

      gradActiveTensorFlowDir = FormGradActiveTensorFlowDir(ndof);
      
      //cout << "scale = " << scale << endl;
      //cout << "dGfD = " << dGfD << endl << endl;
      //cout << "strainIncr = " << strainIncr << endl;
      
      jumpIncrement += scale * strainIncr.Dot(dGfD, strainIncr);
      jumpWork += scale * gradActiveTensorFlowDir.Dot(dGfD,gradActiveTensorFlowDir);		  
    }
  
  //cout << "area = " << area << endl ;
  //cout << "jumpWork = " << jumpWork << endl;		
  
  
  //fJumpIncrement /= area;
  //fJumpIncrement += fBand->Jump()*jumpWork; //already incorporated into strain??
  //jumpWork = -1.0*fabs(jumpWork);

  jumpIncrement /= (jumpWork + area * fBand->H_delta());

  cout << "jumpIncrement = " << jumpIncrement <<endl;

  double trialDeltaResidCohesion = -1.0*fabs(jumpIncrement)*fBand->H_delta();
  /* check to see that residual cohsion does not drop below 0, adjust if nec */
  if (fBand->ResidualCohesion() < trialDeltaResidCohesion)
    {
      cout << fBand->ResidualCohesion() << " " << trialDeltaResidCohesion
	   << endl;


      //full step with no softening
      jumpIncrement *= (jumpWork + area * fBand -> H_delta())/jumpWork;
      // *fraction of step over which no softening occurs
      jumpIncrement *= (trialDeltaResidCohesion -
			fBand->ResidualCohesion())/trialDeltaResidCohesion;
      double jumpIncrementSign = jumpIncrement/fabs(jumpIncrement); 
      //plus amount to get cohesion to zero
      jumpIncrement -= jumpIncrementSign * (fBand->ResidualCohesion())/(fBand-> H_delta());
      
      //fBand->SetEffectiveSoftening(-1.0*fBand->ResidualCohesion()/fabs(jumpIncrement));
      fBand->SetEffectiveSoftening(0.0);
      //fBand -> SetEffectiveSoftening(fBand->ResidualCohesion()/trialDeltaResidCohesion);
    }
  else
    {
      fBand->SetEffectiveSoftening(fBand->H_delta()); 
    } 

  fBand -> StoreJumpIncrement(jumpIncrement);

  cout << "jumpIncrement = " << jumpIncrement <<endl;
  //cout << "fBand->Jump() = " << fBand->Jump() << endl <<endl;

  return jumpIncrement;
}

bool SSEnhLocCraigT::IsBandActive()
{
  /*calculate average stress assuming no jump*/
  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();
	
  double area = 0.0;
  double normalStress = 0.0;
  double shearStress = 0.0;

  /* collect incremental heat */
  //bool need_heat = fElementHeat.Length() == fShapes->NumIP();
  
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
      //	cout << "strainIncr =\n" << strainIncr << endl;
      dSymMatrixT stressIncr(NumSD());
      stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);
      stressIncr += fBand->Stress_List(CurrIP());
      
      area += (*Det)*(*Weight);
      double scale = (*Det++)*(*Weight++);

      normalStress += scale * stressIncr.MultmBn(fBand->Normal(),fBand-> Normal());
      shearStress += scale * stressIncr.MultmBn(fBand->PerpSlipDir(),fBand-> Normal()); 
    }

  if (shearStress < 0.0)
    {
      /* align slip direction with shear stress direction to get 
	 correct yield surface */
      fBand-> FlipSlipDir();
      shearStress *= -1.0;
    }

  normalStress/= area;
  shearStress = shearStress/area;

  double neededCohesion = shearStress + normalStress * fLocalizedFrictionCoeff;
 
  if (fBand-> ResidualCohesion() < neededCohesion)
    {
      fBand-> SetActive(true);
      return true;
    }
  else
    {
      fBand -> SetActive(false);
      return false;
    }
}


void SSEnhLocCraigT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
			   const iArrayT& e_codes, dArray2DT& e_values)
{

  SolidElementT::ComputeOutput(n_codes, n_values, e_codes, e_values);

}


void SSEnhLocCraigT::CloseStep(void)
{

  Top();
  while (NextElement())
    {
      
      /* If element has not localized yet, check */
      if (!IsElementTraced())
	{
	  if (IsElementLocalized())
	    {
	      cout << "hi\n";
	      // isLocalized = true;
	      /* grab stresses */
	      
	      /* temp - this works only for single element */
	      
	      //replae to } with AllcoateBand or move to IsLocalized
	      /*
	      fStress_List.Dimension(NumIP());
	      fShapes->TopIP();
	      for (int i = 0; i < NumIP(); i++)
		{
		  fShapes -> NextIP();
		  cout << "s_ij = " << fCurrMaterial -> s_ij() << endl << flush;
		  fStress_List[i].Dimension(NumSD());
		  fStress_List[i] = fCurrMaterial -> s_ij();
		}
	      */
	    }
	}
      else
	{
	  //fBand -> CloseStep();

	  
	  fBand->IncrementJump();
	  cout << "JumpIncrement = " << fBand -> JumpIncrement() << endl;
	  cout << "Jump = " << fBand -> Jump() << endl; 


	  fShapes->TopIP();
	  while (fShapes->NextIP())
	    {
	      dSymMatrixT strainIncr = fStrain_List [CurrIP()];
	      strainIncr -= fStrain_last_List [CurrIP()]; 
	      dSymMatrixT stressIncr(NumSD());
	      stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);
	      fBand -> IncrementStress(stressIncr, CurrIP());
		//fStress_List[CurrIP()] += stressIncr; 
	    }
	  
	  fBand->UpdateCohesion();


	}
    }
  SmallStrainT::CloseStep();		
}

/***********************************************************************
 * Protected
 ***********************************************************************/

//move to surface mat model?
dSymMatrixT SSEnhLocCraigT::FormdGdSigma(int ndof)
{
  dSymMatrixT dGdSigma(ndof), work(ndof);
  dMatrixT dGNonSym(ndof);

  //cout << "normal =\n" << normal << endl;
  //cout << "slipdir =\n" <<  

  dGNonSym.Outer(fBand->Normal(), fBand->PerpSlipDir());
  dGdSigma.Symmetrize(dGNonSym);
  work.Outer(fBand->Normal());
  dGdSigma.AddScaled(fLocalizedFrictionCoeff, work);

  /* use vector values for mult */
  dGdSigma.ScaleOffDiagonal(2.0);

  return dGdSigma;
}

dSymMatrixT SSEnhLocCraigT::FormGradActiveTensorFlowDir(int ndof)
{
  dSymMatrixT G(ndof);
  dMatrixT G_NonSym(ndof);
  iAutoArrayT activeNodes = fBand->ActiveNodes();
  int A;
  dArrayT grad_f(ndof);
  grad_f = 0.0;

  for (int i=0; i<ndof; i++)
    {
      activeNodes.Top();

      while(activeNodes.Next())      
	{
	  A = activeNodes.Current();  
	  //grad_f[i] += fB(i, (A-1)*ndof +i);
	  grad_f[i] += fB(i, (A)*ndof +i);

	  //cout << "A =\n" << A << endl;
	  //cout << "grad_f =\n" << grad_f << endl;
	}

    }

  //cout << "SlipDir =\n" << fBand->SlipDir() << endl;

  G_NonSym.Outer(grad_f, fBand->SlipDir());
  G.Symmetrize(G_NonSym);

  /* use vector values for mult */
  G.ScaleOffDiagonal(2.0);
  //G.Symmetrize(G_NonSym.Outer(grad_f,m));

  return G;
}

bool SSEnhLocCraigT::IsElementTraced()
{
  int elementNumber = CurrElementNumber();
  bool isTraced = fTracedElements.HasKey(elementNumber);

  if (isTraced)
    LoadBand(elementNumber);

  return isTraced;
}

void SSEnhLocCraigT::LoadBand(int elementNumber)
{
  fBand = fTracedElements[elementNumber];
}

bool SSEnhLocCraigT::IsElementLocalized()
{
  bool locCheck = false;
  double detA, detAMin  = 1.0e99;
  AutoArrayT <dArrayT> normals;
  AutoArrayT <dArrayT> slipDirs;
  AutoArrayT <dArrayT> bestNormals;
  AutoArrayT <dArrayT> bestSlipDirs;
  
  //cannot reach element card during ComputeOutput
  
  /* loop over integration points */
  fShapes->TopIP();
  while ( fShapes->NextIP() )
    {
      
      // DetCheckT checker(fCurrMaterial->s_ij(), fCurrMaterial->c_ijkl(), fInitialModulus);
      
      /*is this necessary? */
      //checker.SetfStructuralMatSupport(*fSSMatSupport);
      if (fCurrMaterial->IsLocalized(normals, slipDirs, detA))
	{
	  locCheck = true;
	  cout << "detA = " << detA << endl;
	  normals.Top();
	  while (normals.Next())
	    cout << "normal (IsElementLocalized) = \n " << normals.Current() << endl;

	  if (detA < detAMin)
	    {
	      detAMin = detA;
	      bestNormals.Dimension(normals);
	      normals.CopyInto(bestNormals);
	      bestSlipDirs.Dimension(slipDirs);
	      slipDirs.CopyInto(bestSlipDirs);
	    }
	}
      //normals.Free();
      //slipDirs.Free();
    }



  if (locCheck)
    ChooseNormals(bestNormals, bestSlipDirs);



  return locCheck;
}

void SSEnhLocCraigT::ChooseNormals(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipDirs)
{

  /*
  normals.Top();
  while(normals.Next())
    cout << "normal = \n" << normals.Current() << endl;
  */

  normals.Top();
  slipDirs.Top();

  normals.Next();
  slipDirs.Next();

  // cout << "normals.Current = \n " << flush << normals.Current() << endl;

  // implement how to choose
  dArrayT normal = normals.Current();
  dArrayT slipDir = slipDirs.Current();
  dArrayT perpSlipDir;

  // cout << "normal = \n" << normal;
  //cout << "slipDir = \n" << slipDir; 

  //make sure slip direction is dilatant
  if (normal.Dot(normal, slipDir)<0.0)
    slipDir *= -1.0;

  perpSlipDir = slipDir;
  perpSlipDir.AddScaled(-1.0*slipDir.Dot(slipDir, normal), normal);
  perpSlipDir.UnitVector();

  if (fNoBandDilation)
    {
      slipDir = perpSlipDir;
    }



  //get centroid - 
  //later this can be point on edge if neighboring element is localized
  dArrayT centroid = Centroid();

  //temp
  //slipDir *= -1.0;
  //perpSlipDir *= -1.0;


  ArrayT<dSymMatrixT> stressList;
  stressList.Dimension(NumIP());

  /*for residual cohesion*/
  double normalStress = 0.0;
  double shearStress = 0.0;
  double area = 0.0;

  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();

  fShapes->TopIP();
  for (int i = 0; i < NumIP(); i++)
    {
      fShapes -> NextIP();
      //cout << "s_ij = " << fCurrMaterial -> s_ij() << endl << flush;
      stressList[i].Dimension(NumSD());
      stressList[i] = fCurrMaterial -> s_ij();

      area += (*Det)*(*Weight);
      double scale = (*Det++)*(*Weight++);

      normalStress += scale * stressList[i].MultmBn(normal, normal);
      shearStress += scale * stressList[i].MultmBn(perpSlipDir, normal); 
    }

  normalStress/= area;
  shearStress = fabs(shearStress)/area;


  double residCohesion = shearStress + normalStress * fLocalizedFrictionCoeff;

  fBand = FormNewBand(normal, slipDir, perpSlipDir, centroid, residCohesion, stressList);
  //fBand = new BandT(slipDir, normal, normal, centroid, this);

  fTracedElements.Insert(CurrElementNumber(), fBand);

  /*
  cout << fBand->Normal() << endl << endl;
  cout << fBand->SlipDir() << endl << endl;
  cout << fBand->PerpSlipDir() << endl << endl;
  cout << centroid << endl << endl;
  cout << fBand->Jump() << endl << endl;
  */

}

BandT* SSEnhLocCraigT::FormNewBand(dArrayT normal, dArrayT slipDir,
dArrayT perpSlipDir, dArrayT coords, double residCohesion, ArrayT<dSymMatrixT>
stressList)
{
return new BandT(normal, slipDir, perpSlipDir, coords, fH_delta_0, residCohesion, stressList, this);
}

dArrayT SSEnhLocCraigT::Centroid()
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
