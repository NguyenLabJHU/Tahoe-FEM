/* $Id: SSEnhLocCraigT.cpp,v 1.1 2004-08-31 16:41:19 cfoster Exp $ */
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
	HookeanMatT(NumDOF()),
	ParameterInterfaceT("element_base"),
	isLocalized(false),
	fBand(NULL)
  //fNeedsOffset(-1),
  //fSSMatSupport(NULL)
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

	/* strain-displacement relation */

	/*
	ParameterT strain_displacement(ParameterT::Enumeration, "strain_displacement");
	strain_displacement.AddEnumeration("standard", kStandardB);
    strain_displacement.AddEnumeration("B-bar", kMeanDilBbar);
    strain_displacement.SetDefault(kStandardB);
	list.AddParameter(strain_displacement);
	*/

	/*ADD PARAMETERS FOR ENHANCED STRAIN HERE */
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


#if 0

/* return the description of the given inline subordinate parameter list */
ParameterInterfaceT* SSEnhLocCraigT::NewSub(const StringT& name) const
{



	if (name == "small_strain_enh_loc_craig_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT(name);
		
		/* list of element block ID's (defined by ElementBaseT) */
		block->AddSub("block_ID_list", ParameterListT::Once);
	
		/* choice of materials lists (inline) */
		block->AddSub("small_strain_enh_loc_craig_material_choice", ParameterListT::Once, true);
	
		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	}
	else /* inherited */
		return SolidElementT::NewSub(name);
}

/* return the description of the given inline subordinate parameter list. */
void SSEnhLocCraigT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "small_strain_enh_loc_craig_material_choice")
	{
		order = ParameterListT::Choice;
		
		/* list of choices */
		sub_lists.AddSub("small_strain_material_1D");
		sub_lists.AddSub("small_strain_material_2D");
		sub_lists.AddSub("small_strain_material_3D");
	}
	else /* inherited */
		SolidElementT::DefineInlineSub(name, order, sub_lists);
}

#endif

void SSEnhLocCraigT::TakeParameterList(const ParameterListT& list)
{
	
	SmallStrainT::TakeParameterList(list);

	/*ADD PARAMETERS FOR ENHACED STRAIN HERE */
	fH_Delta = list.GetParameter("Post-Localization_softening_parameter_H_Delta"); 
	fNoBandDilation = list.GetParameter("Disallow_Dilation_on_Band");
	fLocalizedFrictionCoeff = list.GetParameter("Localized_Friction_Coefficient");

	/* "INITIALIZE" PARAMETERS */

#if 0

  const char caller[] = "SSEnhLocCraigT::TakeParameterList";

	/* strain displacement option before calling SolidElementT::TakeParameterList */
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
		needs[fNeedsOffset + kstrain     ] = mat->Need_Strain();
		needs[fNeedsOffset + kstrain_last] = mat->Need_Strain_last();
		
		/* consistency */
		needs[kNeedDisp] = needs[kNeedDisp] || needs[fNeedsOffset + kstrain];
		needs[KNeedLastDisp] = needs[KNeedLastDisp] || needs[fNeedsOffset + kstrain_last];
	}

	/* what's needed */
	bool need_strain = false;
	bool need_strain_last = false;
	for (int i = 0; i < fMaterialNeeds.Length(); i++) {
		const ArrayT<bool>& needs = fMaterialNeeds[i];
		need_strain = need_strain || needs[fNeedsOffset + kstrain];
		need_strain_last = need_strain_last || needs[fNeedsOffset + kstrain_last];
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

#endif

}

/* extract the list of material parameters */
void SSEnhLocCraigT::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{

  SmallStrainT::CollectMaterialInfo(all_params, mat_params);

#if 0
	const char caller[] = "SSEnhLocCraigT::CollectMaterialInfo";
	
	/* initialize */
	mat_params.Clear();
	
	/* collected material parameters */
	int num_blocks = all_params.NumLists("small_strain_enh_loc_craig_element_block");
	for (int i = 0; i < num_blocks; i++) {

		/* block information */	
		const ParameterListT& block = all_params.GetList("small_strain_enh_loc_craig_element_block", i);
		
		/* resolve material list name */
		if (i == 0) {
			const ParameterListT& mat_list_params = block.GetListChoice(*this, "small_strain_enh_loc_craig_material_choice");
			mat_params.SetName(mat_list_params.Name());
		}
		
		/* collect material parameters */
		const ParameterListT& mat_list = block.GetList(mat_params.Name());
		const ArrayT<ParameterListT>& mat = mat_list.Lists();
		mat_params.AddList(mat[0]);
	}

#endif

}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct a new material support and return a pointer */
MaterialSupportT* SSEnhLocCraigT::NewMaterialSupport(MaterialSupportT* p) const
{

  return SmallStrainT::NewMaterialSupport(p);

#if 0
	/* allocate */
	if (!p) p = new SSMatSupportT(NumDOF(), NumIP());

	/* inherited initializations */
	SolidElementT::NewMaterialSupport(p);
	
	/* set SolidMatSupportT fields */
	SSMatSupportT* ps = TB_DYNAMIC_CAST(SSMatSupportT*, p);
	if (ps) {
		ps->SetLinearStrain(&fStrain_List);
		ps->SetLinearStrain_last(&fStrain_last_List);
	}

	return p;

#endif

}


/* return a pointer to a new material list */
MaterialListT* SSEnhLocCraigT::NewMaterialList(const StringT& name, int size)
{
  //cout << "name = " << name << endl;
  return SmallStrainT::NewMaterialList(name, size);

#if 0

	/* resolve dimension */
	int nsd = -1;
	if (name == "small_strain_material_1D") nsd = 1;
	else if (name == "small_strain_material_2D") nsd = 2;
	else if (name == "small_strain_material_3D") nsd = 3;
	
	/* no match */
	if (nsd == -1) return NULL;

	/* full list */
	if (size > 0)
	{
		/* material support */
		if (!fSSMatSupport) {
			fSSMatSupport = TB_DYNAMIC_CAST(SSMatSupportT*, NewMaterialSupport());
			if (!fSSMatSupport)
				ExceptionT::GeneralFail("SSEnhLocCraigT::NewMaterialList");
		}

		if (nsd == 1)
			return new SSSolidMatList1DT(size, *fSSMatSupport);
		else if (nsd == 2)
			return new SSSolidMatList2DT(size, *fSSMatSupport);
		else if (nsd == 3)
			return new SSSolidMatList3DT(size, *fSSMatSupport);
	}
	else
	{
		if (nsd == 1)
			return new SSSolidMatList1DT;
		else if (nsd == 2)
			return new SSSolidMatList2DT;
		else if (nsd == 3)
			return new SSSolidMatList3DT;
	}
	
	/* no match */
	return NULL;

#endif

}


/* calculate the internal force contribution ("-k*d") */
void SSEnhLocCraigT::FormKd(double constK)
{

  if (isLocalized)
    {

    }
  else
    SmallStrainT::FormKd(constK);

#if 0
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
		fB.MultTx(fCurrMaterial->s_ij(), fNEEvec);

		/* accumulate */
		fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);
		
		/* incremental heat generation */
		if (need_heat) 
			fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
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

      /* check for localization - flag for next time step -moved */
    /*
    if (IsElementLocalized())
      isLocalized = true;
    else
      isLocalized = false;
    */

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
	dSymMatrixT gradActiveTensorFlowDir(ndof), dGdSigma(ndof);
	dMatrixT fLHSWork(nedof), fDfB((HookeanMatT::Modulus().Rows()),nedof);

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
		fD.SetToScaled(scale, HookeanMatT::Modulus());
							
		/* multiply b(transpose) * db, taking account of symmetry, */
		/* and accumulate in elstif */
		fDfB.MultAB(fD, fB);
		fLHSWork.MultATB(fB, fD);
		fLHS +=fLHSWork;
		/* update nMatrix to get rid of fLHSwork */
                //fLHS.MultATB(fB, fD, format, dMatrixT::kAccumulate);	

		//form k_d_zeta
		gradActiveTensorFlowDir = FormGradActiveTensorFlowDir(ndof);

		fDfB.MultTx(gradActiveTensorFlowDir, k_d_zeta, dMatrixT::kAccumulate);

		//form k_zeta_d
		fDfB.MultTx(dGdSigma, k_zeta_d, dMatrixT::kAccumulate);

		//form k_zeta_zeta
		k_zeta_zeta = fD.MultmBn(dGdSigma, gradActiveTensorFlowDir);
		
	}

	k_d_zeta *= 1.0/area;

	k_zeta_zeta *= 1.0/area;
	k_zeta_zeta += fH_Delta;

	fLHS.Outer(k_d_zeta, k_zeta_d, -1.0*k_zeta_zeta, dMatrixT::kAccumulate);

    }


#if 0

	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	/* integrate element stiffness */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/********DEBUG*******/
	bool print = false; 
	int pos = fElementCards.Position(); 
	if (pos == 1&&0)  
	  print = true; 
	/*******************/
	
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{

		double scale = constK*(*Det++)*(*Weight++);
	
		/* strain displacement matrix */
		if (fStrainDispOpt == kMeanDilBbar)
			Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
		else
			Set_B(fShapes->Derivatives_U(), fB);

		/* get D matrix */
		fD.SetToScaled(scale, fCurrMaterial->c_ijkl());
		if (print) cout << "\nmodulus: "<<fCurrMaterial->c_ijkl();
							
		/* multiply b(transpose) * db, taking account of symmetry, */
		/* and accumulate in elstif */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}

#endif

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

	    double ndof = NumDOF();
	    dSymMatrixT dGdSigma = FormdGdSigma(ndof);
	    dSymMatrixT gradActiveTensorFlowDir(ndof);
	    fD.SetToScaled(1.0, HookeanMatT::Modulus());
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
		  area += (*Det)*(*Weight);
		  double scale = (*Det++)*(*Weight++);
		  dSymMatrixT tempStrain = fStrain_List [i];
		  tempStrain.ScaleOffDiagonal(2.0);

		  gradActiveTensorFlowDir = FormGradActiveTensorFlowDir(ndof);

		  fJumpIncrement += scale * tempStrain.Dot(dGfD, tempStrain);
		  jumpWork += gradActiveTensorFlowDir.Dot(dGfD,gradActiveTensorFlowDir);
		  
		}
		fJumpIncrement /= area;
		fJumpIncrement += fBand->Jump()*jumpWork;
		fJumpIncrement /= (jumpWork + fH_Delta);



		/* loop over integration points again */
		for (int i = 0; i < NumIP(); i++)
		{
		  gradActiveTensorFlowDir = FormGradActiveTensorFlowDir(ndof);

		  /*change shear strains back to matrix values */
		  gradActiveTensorFlowDir.ScaleOffDiagonal(0.5);


			/* deformation gradient */
			if (needs[fNeedsOffset + kstrain])
			  {			    

				 fStrain_List[i].AddScaled(fBand->Jump() + fJumpIncrement, gradActiveTensorFlowDir);
			  }

			/* "last" deformation gradient */
			if (needs[fNeedsOffset + kstrain_last])
			{
				 fStrain_List[i].AddScaled(fBand->Jump(), gradActiveTensorFlowDir);
			}
		}

	  }
    }

#if 0
	/* inherited */
	SolidElementT::SetGlobalShape();

	/* material information */
	int material_number = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[material_number];
	
	/* using B-bar */
	if (fStrainDispOpt == kMeanDilBbar)
	{
		/* compute mean of shape function gradients */
		SetMeanGradient(fMeanGradient);

		/* loop over integration points */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/* set B-bar */
			int ip = fShapes->CurrIP();
			Set_B_bar(fShapes->Derivatives_U(ip), fMeanGradient, fB);
	
			/* deformation gradient */
			if (needs[fNeedsOffset + kstrain])
			{
				/* transpose displacement array */
				fLocDisp.ReturnTranspose(fLocDispTranspose);

				/* compute strain using B-bar */
				dSymMatrixT& strain = fStrain_List[ip];
				fB.Multx(fLocDispTranspose, strain);
				strain.ScaleOffDiagonal(0.5);
			}

			/* "last" deformation gradient */
			if (needs[fNeedsOffset + kstrain_last])
			{
				/* transpose displacement array */
				fLocLastDisp.ReturnTranspose(fLocDispTranspose);

				/* compute strain using B-bar */
				dSymMatrixT& strain = fStrain_last_List[ip];
				fB.Multx(fLocDispTranspose, strain);
				strain.ScaleOffDiagonal(0.5);
			}
		}		
	}
	else
	{
		/* loop over integration points */
		for (int i = 0; i < NumIP(); i++)
		{
			/* deformation gradient */
			if (needs[fNeedsOffset + kstrain])
			{
				/* displacement gradient */
				fShapes->GradU(fLocDisp, fGradU, i);

				/* symmetric part */
				 fStrain_List[i].Symmetrize(fGradU);
			}

			/* "last" deformation gradient */
			if (needs[fNeedsOffset + kstrain_last])
			{
				/* displacement gradient */
				fShapes->GradU(fLocLastDisp, fGradU, i);

				/* symmetric part */
				 fStrain_last_List[i].Symmetrize(fGradU);
			}
		}
	} 

#endif
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

  /* If element has not localized yet, check */
  if (!isLocalized)
    {
      if(IsElementLocalized())
	isLocalized = true;
    }
  else
    fBand->IncrementJump(fJumpIncrement);


}

/***********************************************************************
 * Private
 ***********************************************************************/


#if 0
/*private function - not needed - yet */

/* compute mean shape function gradient, Hughes (4.5.23) */
/* for b-bar only ? */
void SSEnhLocCraigT::SetMeanGradient(dArray2DT& mean_gradient) const
{


	int nip = NumIP();
	const double* det = fShapes->IPDets();
	const double*   w = fShapes->IPWeights();

	/* volume */
	double vol = 0.0;
	for (int i = 0; i < nip; i++)
		vol += w[i]*det[i];

	/* initialize */
	mean_gradient = 0.0;			

	/* integrate */
	for (int i = 0; i < nip; i++)
		mean_gradient.AddScaled(w[i]*det[i]/vol, fShapes->Derivatives_U(i)); 


}

#endif


/***********************************************************************
 * Protected
 ***********************************************************************/

//move to surface mat model?
dSymMatrixT SSEnhLocCraigT::FormdGdSigma(double ndof)
{
  dSymMatrixT dGdSigma(ndof), work(ndof);
  dMatrixT dGNonSym(ndof);

  dGNonSym.Outer(fNormal, fBand->PerpSlipDir());
  dGdSigma.Symmetrize(dGNonSym);
  work.Outer(fNormal);
  dGdSigma.AddScaled(fLocalizedFrictionCoeff, work);

  /* use vector values for mult */
  dGdSigma.ScaleOffDiagonal(2.0);

  return dGdSigma;
}



dSymMatrixT SSEnhLocCraigT::FormGradActiveTensorFlowDir(double ndof)
{
  dSymMatrixT G(ndof);
  dMatrixT G_NonSym(ndof);
  iAutoArrayT activeNodes = fBand->ActiveNodes();
  int A;
  dArrayT grad_f(ndof);

  for (int i=0; i<ndof; i++)
    {
      activeNodes.Top();
      do
	{
	  A = activeNodes.Current();  
	  grad_f[i] += fB(i, (A-1)*ndof +i);
	}
      while(activeNodes.Next());
    }

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

  /* loop over integration points */
  for (int i = 0; i < NumIP(); i++)
    {

      DetCheckT checker(fCurrMaterial->s_ij(), fCurrMaterial->c_ijkl(), HookeanMatT::Modulus());
      
      /*is this necessary? */
      checker.SetfStructuralMatSupport(*fSSMatSupport);
      if (checker.IsLocalized_SS(normals, slipDirs, detA))
	{
	  locCheck = true;
	  if (detA < detAMin)
	    {
	      detAMin = detA;
	      normals.CopyInto(bestNormals);
	      slipDirs.CopyInto(bestSlipDirs);
	    }
	}
    }

  if (locCheck)
    ChooseNormals(bestNormals, bestSlipDirs);

  return locCheck;
}

void SSEnhLocCraigT::ChooseNormals(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipDirs)
{

  normals.Top();
  slipDirs.Top();

  // implement how to choose
  dArrayT normal = normals.Current();
  dArrayT slipDir = slipDirs.Current();
  dArrayT perpSlipDir;

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




  fBand = new BandT(normal, slipDir, perpSlipDir, centroid, this);

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
      //area += (*Det)*(*Weight);
      double scale = (*Det++)*(*Weight++);
      area += scale;
      fShapes->IPCoords(coords);
      centroid.AddScaled(scale,coords);

    }
  centroid *= 1.0/area;
  
  return centroid;
 
}
