/* $Id: SSEnhLocCraigT.cpp,v 1.5 2005-02-25 03:22:40 cfoster Exp $ */
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
	//	HookeanMatT(NumDOF()),
	//ParameterInterfaceT("element_base"), //only needed for MI w/HookeanMatT
	isLocalized(false),
	isLocalizedTemp(false),
	fBand(NULL),
	fInitialModulus(1)
  //fNeedsOffset(-1),
  //fSSMatSupport(NULL)
{
	SmallStrainT::SetName("small_strain_enh_loc_craig");
	//	HookeanMatT::Dimension(NumDOF());
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
	list.AddParameter(fH_Delta, "Post-Localization_softening_parameter_H_Delta"); 
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
	
	SmallStrainT::TakeParameterList(list);

	/*PARAMETERS FOR ENHANCED STRAIN*/
	fH_Delta = list.GetParameter("Post-Localization_softening_parameter_H_Delta"); 
	fNoBandDilation = list.GetParameter("Disallow_Dilation_on_Band");
	fLocalizedFrictionCoeff = list.GetParameter("Localized_Friction_Coefficient");

	/* "INITIALIZE" PARAMETERS */
	//fInitialModulus = fCurrMaterial->c_ijkl();
	//can't initialize this here fCurrMaterial not set
	fInitialModulus = 0.0;

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

  if (fInitialModulus == 0.0)
    {
      fInitialModulus.Dimension(NumDOF());
      fInitialModulus = fCurrMaterial->c_ijkl();
    }

  if(!isLocalized)
    SmallStrainT::FormKd(constK);

    //cout << "fRHS =\n" << fRHS <<endl;

#if 1
    if(isLocalized)
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
	    stressIncr.A_ijkl_B_kl(fInitialModulus, strainIncr);
	    stressIncr += fStress_List[CurrIP()];
	    fB.MultTx(stressIncr, fNEEvec);
	    
	    /* accumulate */
	    fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);
	    
	    /* incremental heat generation */
	    if (need_heat) 
	      fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	  }     
      }

#endif

}

/* form the element stiffness matrix */
void SSEnhLocCraigT::FormStiffness(double constK)
{


  if (fStrainDispOpt == kMeanDilBbar)
    cout << "Warning SSEnhLocCraigT::FormStiffness, b-bar integration not implemented for enhanced strain element, postlocalization results may be garbage.";

  if (!isLocalized)
    {
    /* form stiffness in standard way */
    SmallStrainT::FormStiffness(constK);

      /* check for localization */
    
    if (IsElementLocalized())
      isLocalizedTemp = true;
    else
      isLocalizedTemp = false;
    

    }
  else //if already localized, use localized stiffness routine
    {

      cout << "constK =\n" << constK << endl;

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
	dMatrixT fLHSWork(nedof),// fDfB((HookeanMatT::Modulus().Rows()),nedof);
	  fDfB((fInitialModulus.Rows()),nedof);

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
		//fD.SetToScaled(scale, HookeanMatT::Modulus());
		fD.SetToScaled(scale, fInitialModulus);

		//cout << "fLHS =\n" << fLHS << endl;		
		//cout << "fD =\n" << fD << endl;
		cout << "fB =\n" << fB << endl;

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
		cout << "k_d_zeta =\n" << k_d_zeta << endl;

		//form k_zeta_d
		//fDfB.MultTx(dGdSigma, k_zeta_d_work, dMatrixT::kOverwrite);
		fDfB.MultTx(dGdSigma, k_zeta_d_work);
		k_zeta_d += k_zeta_d_work;
		cout << "k_zeta_d =\n" << k_zeta_d << endl;

		//form k_zeta_zeta
		k_zeta_zeta +=
		fD.MultmBn(dGdSigma,gradActiveTensorFlowDir);
		cout << "k_zeta_zeta =\n" << k_zeta_zeta << endl;
	}

	cout << "fLHS =\n" << fLHS << endl;
	//cout << "area =" << area << endl;

	k_d_zeta *= 1.0/area;

	k_zeta_zeta *= 1.0/area;
	k_zeta_zeta += fH_Delta;
	//k_zeta_zeta *= -1.0;



	fLHS.Outer(k_d_zeta, k_zeta_d, -1.0/k_zeta_zeta, dMatrixT::kAccumulate);
	//cout << "fLHS =\n" << fLHS << endl;
    }
}

/* compute the measures of strain/deformation over the element */
void SSEnhLocCraigT::SetGlobalShape(void)
{

  SmallStrainT::SetGlobalShape();

  /* subtract band deformation */
  if (isLocalized)
    {
	if (fStrainDispOpt == kMeanDilBbar)
	  cout << "Warning - B-bar not implemented for localized element\n";
	else
	  {

	    int material_number = CurrentElement().MaterialNumber();
	    const ArrayT<bool>& needs = fMaterialNeeds[material_number];

	    int ndof = NumDOF();
	    dSymMatrixT dGdSigma = FormdGdSigma(ndof);
	    dSymMatrixT gradActiveTensorFlowDir(ndof);

	    //cout << "fInitialModulus =\n" << fInitialModulus << endl;
	    //cout << "fD =\n" << fD << endl;


	    //fD.SetToScaled(1.0, HookeanMatT::Modulus());
	    fD.SetToScaled(1.0, fInitialModulus);
	    dArrayT dGfD(fD.Rows());
	    fD.MultTx(dGdSigma, dGfD);

	    const double* Det    = fShapes->IPDets();
	    const double* Weight = fShapes->IPWeights();	    
	    double area = 0.0;
	    fJumpIncrement = 0.0;
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
 
		  fJumpIncrement += scale * strainIncr.Dot(dGfD, strainIncr);
		  jumpWork += scale * gradActiveTensorFlowDir.Dot(dGfD,gradActiveTensorFlowDir);		  
		}
		
		//cout << "area = " << area << endl ;
		//cout << "jumpWork = " << jumpWork << endl;		
		 
		
		//fJumpIncrement /= area;
		//fJumpIncrement += fBand->Jump()*jumpWork; //already incorporated into strain??
		fJumpIncrement /= (jumpWork + area * fH_Delta);
		//fJumpIncrement *= 1.1182;   

		cout << "fJumpIncrement = " << fJumpIncrement <<endl;
		//cout << "fBand->Jump() = " << fBand->Jump() << endl <<endl;

		/* loop over integration points again */
		for (int i = 0; i < NumIP(); i++)
		{
		  gradActiveTensorFlowDir = FormGradActiveTensorFlowDir(ndof);

		  /*change shear strains back to matrix values */
		  gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);


			/* deformation gradient */
			if (needs[fNeedsOffset + kstrain])
			  {			    
				 fStrain_List[i].AddScaled(-(fBand->Jump() + fJumpIncrement), gradActiveTensorFlowDir);
			  }

			/* "last" deformation gradient */ //is this right?
			if (needs[fNeedsOffset + kstrain_last])
			{
				 fStrain_last_List[i].AddScaled(-(fBand->Jump()), gradActiveTensorFlowDir);
			}
		}
	  }
    }
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/** driver for calculating output values */
/* Used to check localization - is there a more appropriate fn? */


void SSEnhLocCraigT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
			   const iArrayT& e_codes, dArray2DT& e_values)
{


  SolidElementT::ComputeOutput(n_codes, n_values, e_codes, e_values);

#if 0

  /* If element has not localized yet, check */
  if (!isLocalized)
    {
      if (isLocalizedTemp)
	{
	isLocalized = true;
	
#if 0

	/* allocate for current strain if material did not */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain])
	  {
	    fStrain_List.Dimension(NumIP());
	    for (int i = 0; i < NumIP(); i++)
	      fStrain_List[i].Dimension(NumSD());
	    
	    needs[kNeedDisp] = true;
	  }
	
	/* allocate for last strain if material did not */
	if (!needs[fNeedsOffset + kstrain_last])
	  {
	    fStrain_last_List.Dimension(NumIP());
	    for (int i = 0; i < NumIP(); i++)
	      fStrain_last_List[i].Dimension(NumSD());
	    
	    needs[KNeedLastDisp] = true;
	  }
			
#endif

	/* grab stresses */
	cout << "NumIP() =" << NumIP() << endl << flush;
	cout << "Density =" << fCurrMaterial->Density() << endl << flush;
	
	/* temp - this works only for single element */
	Top();
	while (NextElement())
	  {

	fStress_List.Dimension(NumIP());
	fShapes->TopIP();
	for (int i = 0; i < NumIP(); i++)
	  {
	    fShapes -> NextIP();
	    cout << "s_ij = " << fCurrMaterial->s_ij() << endl << flush;
	    fStress_List[i].Dimension(NumSD());
	    fStress_List[i] = fCurrMaterial -> s_ij();
	    //fShapes -> NextIP();
	  }    
	  } // end temp - while(NextElement())
	} 
    }

  else
    {
      fBand->IncrementJump(fJumpIncrement);
      cout << "Jump = " << fBand -> Jump() << endl; 
    }

  cout << "isLocalized = " << isLocalized << endl;
  cout << "isLocalizedTemp = " << isLocalizedTemp << endl; 

#endif

}


void SSEnhLocCraigT::CloseStep(void)
{

 /* If element has not localized yet, check */
  if (!isLocalized)
    {
      if (isLocalizedTemp)
	{
	  isLocalized = true;
	
#if 0
	  
	  /* allocate for current strain if material did not */
	  int mat_num = CurrentElement().MaterialNumber();
	  const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	  if (!needs[fNeedsOffset + kstrain])
	    {
	      fStrain_List.Dimension(NumIP());
	      for (int i = 0; i < NumIP(); i++)
		fStrain_List[i].Dimension(NumSD());
	      
	      needs[kNeedDisp] = true;
	    }
	  
	  /* allocate for last strain if material did not */
	  if (!needs[fNeedsOffset + kstrain_last])
	    {
	      fStrain_last_List.Dimension(NumIP());
	      for (int i = 0; i < NumIP(); i++)
		fStrain_last_List[i].Dimension(NumSD());
	      
	      needs[KNeedLastDisp] = true;
	    }
			
#endif

	  /* grab stresses */
	  //cout << "NumIP() =" << NumIP() << endl << flush;

	/* temp - this works only for single element */
	Top();
	while (NextElement())
	  {

	    fStress_List.Dimension(NumIP());
	    fShapes->TopIP();
	    for (int i = 0; i < NumIP(); i++)
	      {
		fShapes -> NextIP();
		cout << "s_ij = " << fCurrMaterial -> s_ij() << endl << flush;
		fStress_List[i].Dimension(NumSD());
		fStress_List[i] = fCurrMaterial -> s_ij();
	      }
	  }
		}
    }
  else
    {
      fBand->IncrementJump(fJumpIncrement);
      cout << "Jump = " << fBand -> Jump() << endl; 
	  
	  Top();
	  while(NextElement())
	  {
	   
	   fShapes -> TopIP();
	   while(fShapes->NextIP())
	   {
		dSymMatrixT strainIncr = fStrain_List [CurrIP()];
		strainIncr -= fStrain_last_List [CurrIP()]; 
		dSymMatrixT stressIncr(NumSD());
		stressIncr.A_ijkl_B_kl(fInitialModulus, strainIncr);
		fStress_List[CurrIP()] += stressIncr;
	   }
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

bool SSEnhLocCraigT::IsElementLocalized()
{
      bool locCheck = false;
      double detA, detAMin  = 1.0;
      AutoArrayT <dArrayT> normals;
      AutoArrayT <dArrayT> slipDirs;
      AutoArrayT <dArrayT> bestNormals;
      AutoArrayT <dArrayT> bestSlipDirs;

      //cannot reach element card during ComputeOutput


  /* loop over integration points */
      fShapes->TopIP();
      while ( fShapes->NextIP() )
    {
      //DetCheckT checker(fCurrMaterial->s_ij(), fCurrMaterial->c_ijkl(), HookeanMatT::Modulus());
      DetCheckT checker(fCurrMaterial->s_ij(), fCurrMaterial->c_ijkl(), fInitialModulus);
      
      /*is this necessary? */
      checker.SetfStructuralMatSupport(*fSSMatSupport);
      if (checker.IsLocalized_SS(normals, slipDirs, detA))
	{
	  locCheck = true;
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
    ChooseNormals(normals, slipDirs);



  return locCheck;
}

void SSEnhLocCraigT::ChooseNormals(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipDirs)
{

  normals.Top();
  slipDirs.Top();

  normals.Next();
  slipDirs.Next();

  // implement how to choose
  dArrayT normal = normals.Current();
  dArrayT slipDir = slipDirs.Current();
  dArrayT perpSlipDir;

  //cout << "normal = \n" << normal;
  //cout << "slipDir = \n" << slipDir; 

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

  fBand = new BandT(normal, slipDir, perpSlipDir, centroid, this);
  //fBand = new BandT(slipDir, normal, normal, centroid, this);

  /*
  cout << fBand->Normal() << endl << endl;
  cout << fBand->SlipDir() << endl << endl;
  cout << fBand->PerpSlipDir() << endl << endl;
  cout << centroid << endl << endl;
  cout << fBand->Jump() << endl << endl;
  */

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
