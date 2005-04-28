/* $Id: SSEnhLocCraigT.cpp,v 1.14 2005-04-28 00:45:55 cfoster Exp $ */
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
#include <math.h>

using namespace Tahoe;

/*initialize static variables */
bool SSEnhLocCraigT::fLocalizationHasBegun = false;
double SSEnhLocCraigT::fDetAMin = 1.0;

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
{}
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
  SmallStrainT::DefineSubs(sub_list);
}

void SSEnhLocCraigT::TakeParameterList(const ParameterListT& list)
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

	    
	    dSymMatrixT gradActiveTensorFlowDir =
	    FormGradActiveTensorFlowDir(NumSD(), CurrIP());
	    gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);

	    strainIncr.AddScaled(-1.0*fBand->JumpIncrement(), gradActiveTensorFlowDir);
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
  if (fStrainDispOpt == kMeanDilBbar)
    cout << "Warning SSEnhLocCraigT::FormStiffness, b-bar integration not implemented for enhanced strain element, postlocalization results may be garbage.";

  if (!IsElementTraced() || !fBand->IsActive() )
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
		gradActiveTensorFlowDir =
		FormGradActiveTensorFlowDir(ndof, CurrIP());

		fDfB.MultTx(gradActiveTensorFlowDir, k_d_zeta_work);
		k_d_zeta += k_d_zeta_work;

		//form k_zeta_d
		//fDfB.MultTx(dGdSigma, k_zeta_d_work, dMatrixT::kOverwrite);
		fDfB.MultTx(dGdSigma, k_zeta_d_work);
		k_zeta_d += k_zeta_d_work;

		//form k_zeta_zeta
		k_zeta_zeta +=
		fD.MultmBn(dGdSigma,gradActiveTensorFlowDir);
	}
	k_zeta_d *= 1.0/area;

	k_zeta_zeta *= 1.0/area;
	k_zeta_zeta += fBand->EffectiveSoftening();
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
      double scale = (*Det++)*(*Weight++);
      area += scale;
      dSymMatrixT strainIncr = fStrain_List [i];
      strainIncr -= fStrain_last_List [i];
      strainIncr.ScaleOffDiagonal(2.0);

      gradActiveTensorFlowDir = FormGradActiveTensorFlowDir(ndof, i);
      
      jumpIncrement += scale * strainIncr.Dot(dGfD, strainIncr);
      jumpWork += scale * gradActiveTensorFlowDir.Dot(dGfD,gradActiveTensorFlowDir);		  
    }
  jumpIncrement /= (jumpWork + area * fBand->H_delta());

  double trialDeltaResidCohesion = -1.0*fabs(jumpIncrement)*fBand->H_delta();
  /* check to see that residual cohsion does not drop below 0, adjust if nec */
  if (fBand->ResidualCohesion() < trialDeltaResidCohesion)
    {
      //full step with no softening
      jumpIncrement *= (jumpWork + area * fBand -> H_delta())/jumpWork;

      // *fraction of step over which no softening occurs
      jumpIncrement *= (trialDeltaResidCohesion -
			fBand->ResidualCohesion())/trialDeltaResidCohesion;
      double jumpIncrementSign = jumpIncrement/fabs(jumpIncrement); 

      //plus amount to get cohesion to zero
      jumpIncrement -= jumpIncrementSign * (fBand->ResidualCohesion())/(fBand-> H_delta());
      
      fBand->SetEffectiveSoftening(0.0);
    }
  else
    {
      fBand->SetEffectiveSoftening(fBand->H_delta()); 
    } 
  fBand -> StoreJumpIncrement(jumpIncrement);

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
  if (fLocalizationHasBegun)
    {
      /*update traced elements */ 
      Top();
      while (NextElement())
	{
	  if (IsElementTraced())
	    {
	      dSymMatrixT strainIncr = fStrain_List [CurrIP()];
	      strainIncr -= fStrain_last_List [CurrIP()]; 
	      
	      dSymMatrixT gradActiveTensorFlowDir =
		FormGradActiveTensorFlowDir(NumSD(), CurrIP());
	      gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);
	      
	      strainIncr.AddScaled(-1.0*fBand->JumpIncrement(), gradActiveTensorFlowDir);
	      
	      dSymMatrixT stressIncr(NumSD());
	      stressIncr.A_ijkl_B_kl(fCurrMaterial->ce_ijkl(), strainIncr);
	      fBand -> IncrementStress(stressIncr, CurrIP());
	      
	      cout << "jumpIncrement = " << fBand->JumpIncrement() << endl;
	      cout << "fBand->Jump() = " << fBand->Jump() << endl;	
	      
	      fBand -> CloseStep();
	    }
	}
      /* check for newly localized elements */
      fEdgeOfBandElements.Top();
      while(fEdgeOfBandElements.Next())
	{
	  GetElement(fEdgeOfBandElements.Current());
	  IsElementLocalized();
	}
    }
  else
    {
      //choose first element then let band progress
      bool localizationHasBegun = false;
      Top();
      while (NextElement())
	{
	  if (IsElementLocalized())
	    localizationHasBegun = true;
	}
      if (localizationHasBegun)
	{
	  fLocalizationHasBegun = true;
	  fEdgeOfBandElements.Top();
	  //localize 1st element?
	  while(fEdgeOfBandElements.Next())
	    {
	      GetElement(fEdgeOfBandElements.Current());    
	      IsElementLocalized();
	    }
	}
    }

  SmallStrainT::CloseStep();		
}


/***********************************************************************
 * Protected
 ***********************************************************************/

/* current element operations */
void SSEnhLocCraigT::GetElement(int elementNumber)
{
  /* inherited */
  //bool result = ContinuumElementT::NextElement();
  fElementCards.Current(elementNumber);  

  /* get material pointer */
  ContinuumMaterialT* pcont_mat = (*fMaterialList)[CurrentElement().MaterialNumber()];
      
  /* cast is safe since class contructs materials list */
  fCurrMaterial = (SolidMaterialT*) pcont_mat;
 }




//move to surface mat model?
dSymMatrixT SSEnhLocCraigT::FormdGdSigma(int ndof)
{
  dSymMatrixT dGdSigma(ndof), work(ndof);
  dMatrixT dGNonSym(ndof);

  dGNonSym.Outer(fBand->Normal(), fBand->PerpSlipDir());
  dGdSigma.Symmetrize(dGNonSym);
  work.Outer(fBand->Normal());
  dGdSigma.AddScaled(fLocalizedFrictionCoeff, work);

  /* use vector values for mult */
  dGdSigma.ScaleOffDiagonal(2.0);

  return dGdSigma;
}

dSymMatrixT SSEnhLocCraigT::FormGradActiveTensorFlowDir(int ndof, int ip)
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
	  grad_f[i] += fB(i, (A)*ndof +i);
	}
    }

  G_NonSym.Outer(grad_f, fBand->SlipDir());
  G.Symmetrize(G_NonSym);

  /* use vector values for mult */
  G.ScaleOffDiagonal(2.0);

  return G;
}


bool SSEnhLocCraigT::IsElementTraced()
{
 	 //int elementNumber = CurrElementNumber();
 	 return IsElementTraced(CurrElementNumber());
}
  
bool SSEnhLocCraigT::IsElementTraced(int elementNumber)
{
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
  
  cout << "hi \n ";

  /* loop over integration points */
  fShapes->TopIP();
  while ( fShapes->NextIP() )
    {
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
    if (fLocalizationHasBegun)
      ChooseNormals(bestNormals, bestSlipDirs);
    else
      if (detAMin < fDetAMin)
      {
	fDetAMin = detAMin;
	fEdgeOfBandElements.Free();
	//EdgeOfBandElement* newElement = new EdgeOfBandElement;
	//newElement->elementNumber = CurrElementNumber();
	dArrayT coords = Centroid();
	fEdgeOfBandElements.Append(CurrElementNumber());
	fEdgeOfBandCoords.Append(coords);
      }

  return locCheck;
}

void SSEnhLocCraigT::ChooseNormals(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipDirs)
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
      cout << "fGradU = \n" << fGradU << endl;
    }
  avgGradU /= area;

  double prod, maxProd = -1.0;

  while(normals.Next())
    {
      slipDirs.Next();
      prod = fabs( avgGradU.MultmBn(normals.Current(), slipDirs.Current()));
      if (prod > maxProd)
	{
	  normal = normals.Current(); 
	  slipDir = slipDirs.Current();
	}
    }

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
  //dArrayT centroid = Centroid();

  ArrayT<dSymMatrixT> stressList;
  stressList.Dimension(NumIP());

  /*for residual cohesion*/
  double normalStress = 0.0;
  double shearStress = 0.0;
  //double area = 0.0;

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

  fBand = FormNewBand(normal, slipDir, perpSlipDir, fEdgeOfBandCoords.Current(), residCohesion, stressList);

  fTracedElements.Insert(CurrElementNumber(), fBand);

  AddNewEdgeElements(CurrElementNumber());
}

BandT* SSEnhLocCraigT::FormNewBand(dArrayT normal, dArrayT slipDir,
dArrayT perpSlipDir, dArrayT coords, double residCohesion, ArrayT<dSymMatrixT>
stressList)
{
return new BandT(normal, slipDir, perpSlipDir, coords, fH_delta_0, residCohesion, stressList, this);
}

void SSEnhLocCraigT::AddNewEdgeElements(int elementNumber)
{
  ModelManagerT& model = ElementSupport().ModelManager();
  iArray2DT neighbors;

  model.ElementNeighbors(model.ElementGroupIDs(), neighbors);

  //2D
  int numSides;
  int numSidesFound = 0;
  int numElementSides = 0;
  iAutoArrayT activeNodes = fBand->ActiveNodes();

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
  dArrayT nodalCoord1, nodalCoord2; //coords a particular node

  for(int i = 0; i < numElementSides; i++)
    if ((activeNodes.HasValue(i+1 % numSides) && !activeNodes.HasValue(i))
	|| (!activeNodes.HasValue(i+1 % numSides) && activeNodes.HasValue(i)))
      {
		if (!(neighbors(elementNumber,i) = -1 || IsElementTraced(neighbors(elementNumber ,i))))
	  	{
	    	//get coords
	    	dArrayT localizedEleCoord = fBand -> Coords();

		    for (int j = 0; j < nodalCoords.MinorDim(); j++)
		      {		
				nodalCoord1 [j] = nodalCoords(i,j);
				nodalCoord2 [j] = nodalCoords(i+1 % numSides, j);
	      		}

	    	dArrayT interceptCoords = InterceptCoords(localizedEleCoord,
						      nodalCoord1, nodalCoord2);
	    	//stick element in fEdgeOfBandElements
	    	//EdgeOfBandElement* newElement = new EdgeOfBandElement;
	    	//newElement->elementNumber = CurrElementNumber();
	    	//newElement -> coords = interceptCoords;
	    	if(fEdgeOfBandElements.AppendUnique(CurrElementNumber()))
		  fEdgeOfBandCoords.Append(interceptCoords);
	  	} 
		if (++numSidesFound > 1) 
	  	break;
      }
  cout << "numSidesFound = " << numSidesFound <<endl;

}

dArrayT SSEnhLocCraigT::InterceptCoords(dArrayT& localizedEleCoord,
dArrayT& nodalCoord1, dArrayT& nodalCoord2)
{
  //assumes staright sides
  dArrayT sideVector = nodalCoord2;
  sideVector -= nodalCoord1;

  dArrayT perpSlipDir = fBand -> PerpSlipDir();

  double alpha = sideVector[1] * (localizedEleCoord[2] - nodalCoord1[2]) -
  sideVector[2] * (localizedEleCoord[1] - nodalCoord1[1]);
  alpha /= sideVector[2] * perpSlipDir[1] - sideVector[1] *
  perpSlipDir[2];

  dArrayT interceptCoord = localizedEleCoord;
  interceptCoord.AddScaled(alpha, perpSlipDir);

  return interceptCoord;
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
