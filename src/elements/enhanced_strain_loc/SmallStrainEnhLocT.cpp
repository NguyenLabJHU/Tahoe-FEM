/* $Id: SmallStrainEnhLocT.cpp,v 1.4 2005-02-02 21:12:57 raregue Exp $ */
#include "SmallStrainEnhLocT.h"
#include "ShapeFunctionT.h"
#include "SSSolidMatT.h"

/* materials lists */
#include "SSSolidMatList1DT.h"
#include "SSSolidMatList2DT.h"
#include "SSSolidMatList3DT.h"
#include "SSMatSupportT.h"
#include "ParameterContainerT.h"
#include "ModelManagerT.h"

using namespace Tahoe;

/* initialize static variables */
bool SmallStrainEnhLocT::fFirstPass = true;

/* constructor */
SmallStrainEnhLocT::SmallStrainEnhLocT(const ElementSupportT& support):
	SolidElementT(support),
	fNeedsOffset(-1),
	fSSMatSupport(NULL)
{
	SetName("small_strain_enh_loc");
}

/* destructor */
SmallStrainEnhLocT::~SmallStrainEnhLocT(void)
{
	delete fSSMatSupport;
}

/* finalize current step - step is solved */
void SmallStrainEnhLocT::CloseStep(void)
{
	/* inherited */
	SolidElementT::CloseStep();
	
	/* store converged solution */
	fElementLocScalars_last = fElementLocScalars;
	fElementLocSlipDir_last = fElementLocSlipDir;
	fElementLocMuDir_last = fElementLocMuDir;
	fElementLocInternalVars_last = fElementLocInternalVars;
	fElementVolume_last = fElementVolume;
	
	ss_enh_out	<< setw(outputFileWidth) << "\n\n time_step    element " << setw(outputFileWidth) << "?" 
			<< setw(outputFileWidth) << "?" << endl;
	ss_enh_out	<< endl << "----------------------------------------------------------------------------------------------" << endl;

}
	
/* restore last converged state */
GlobalT::RelaxCodeT SmallStrainEnhLocT::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = SolidElementT::ResetStep();
	
	/* restore converged solution */
	fElementLocScalars = fElementLocScalars_last;
	fElementLocSlipDir = fElementLocSlipDir_last;
	fElementLocMuDir = fElementLocMuDir_last;
	fElementLocInternalVars = fElementLocInternalVars_last;
	fElementVolume = fElementVolume_last;
	
	return relax;
}

/* read restart information from stream */
void SmallStrainEnhLocT::ReadRestart(istream& in)
{
	/* inherited */
	SolidElementT::ReadRestart(in);
	
	/* read restart data */
	in >> fElementLocScalars;
	in >> fElementLocSlipDir;
	in >> fElementLocMuDir;
	in >> fElementLocInternalVars;
	in >> fElementVolume;
	
	/* reset last state */
	fElementLocScalars_last = fElementLocScalars;
	fElementLocSlipDir_last = fElementLocSlipDir;
	fElementLocMuDir_last = fElementLocMuDir;
	fElementLocInternalVars_last = fElementLocInternalVars;
	fElementVolume_last = fElementVolume;
}

/* write restart information from stream */
void SmallStrainEnhLocT::WriteRestart(ostream& out) const
{
	/* inherited */
	SolidElementT::WriteRestart(out);
	
	/* write restart data */
	out << fElementLocScalars << '\n';
	out << fElementLocSlipDir << '\n';
	out << fElementLocMuDir << '\n';
	out << fElementLocInternalVars << '\n';
	out << fElementVolume << '\n';
}

/* implementation of the ParameterInterfaceT interface */
void SmallStrainEnhLocT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SolidElementT::DefineParameters(list);

	/* strain-displacement relation */
	ParameterT strain_displacement(ParameterT::Enumeration, "strain_displacement");
	strain_displacement.AddEnumeration("standard", kStandardB);
    strain_displacement.AddEnumeration("B-bar", kMeanDilBbar);
    strain_displacement.SetDefault(kStandardB);
	list.AddParameter(strain_displacement);
	
	/* cohesive surface model constants */
	double cohesion_r, cohesion_p, alpha_c, phi_r, phi_p, alpha_phi, psi_p, alpha_psi;
	list.AddParameter(cohesion_r, "residual_cohesion");
	list.AddParameter(cohesion_p, "peak_cohesion");
	list.AddParameter(alpha_c, "cohesion_softening_coefficient");
	list.AddParameter(phi_r, "residual_friction_angle_rad");
	list.AddParameter(phi_p, "peak_friction_angle_rad");
	list.AddParameter(alpha_phi, "friction_softening_coefficient");
	list.AddParameter(psi_p, "peak_dilation_angle_rad");
	list.AddParameter(alpha_psi, "dilation_softening_coefficient");
	
}

/* information about subordinate parameter lists */
void SmallStrainEnhLocT::DefineSubs(SubListT& sub_list) const
{	
	/* inherited */
	SolidElementT::DefineSubs(sub_list);

	/* element block/material specification */
	sub_list.AddSub("small_strain_enh_loc_element_block", ParameterListT::OnePlus);
}

/* return the description of the given inline subordinate parameter list */
ParameterInterfaceT* SmallStrainEnhLocT::NewSub(const StringT& name) const
{
	if (name == "small_strain_enh_loc_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT(name);
		
		/* list of element block ID's (defined by ElementBaseT) */
		block->AddSub("block_ID_list", ParameterListT::Once);
	
		/* choice of materials lists (inline) */
		block->AddSub("small_strain_enh_loc_material_choice", ParameterListT::Once, true);
	
		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	}
	else /* inherited */
		return SolidElementT::NewSub(name);
}

/* return the description of the given inline subordinate parameter list. */
void SmallStrainEnhLocT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "small_strain_enh_loc_material_choice")
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

void SmallStrainEnhLocT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SmallStrainEnhLocT::TakeParameterList";

	/* strain displacement option before calling SolidElementT::TakeParameterList */
	int b = list.GetParameter("strain_displacement");
	fStrainDispOpt = (b == kStandardB) ? kStandardB : kMeanDilBbar;

	/* inherited */
	SolidElementT::TakeParameterList(list);
	
	/* dimension workspace */
	fGradU.Dimension(NumSD());	
	if (fStrainDispOpt == kMeanDilBbar) {
		fLocDispTranspose.Dimension(fLocDisp.Length());
		//fMeanGradient.Dimension(NumSD(), NumElementNodes());
	}	
	fMeanGradient.Dimension(NumSD(), NumElementNodes());

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
	
	
	/** \name localization element info storage */
	/*@{*/
	/** fixed values */
	fElementLocNormal.Dimension(NumElements(),NumSD());
	fElementLocNormal = 0.0;
	fElementLocTangent.Dimension(NumElements(),NumSD());
	fElementLocTangent = 0.0;
	fElementLocNormal1.Dimension(NumElements(),NumSD());
	fElementLocNormal1 = 0.0;
	fElementLocNormal2.Dimension(NumElements(),NumSD());
	fElementLocNormal2 = 0.0;
	fElementLocNormal3.Dimension(NumElements(),NumSD());
	fElementLocNormal3 = 0.0;
	fElementLocSlipDir1.Dimension(NumElements(),NumSD());
	fElementLocSlipDir1 = 0.0;
	fElementLocSlipDir2.Dimension(NumElements(),NumSD());
	fElementLocSlipDir2 = 0.0;
	fElementLocSlipDir3.Dimension(NumElements(),NumSD());
	fElementLocSlipDir3 = 0.0;
	fElementLocTangent1.Dimension(NumElements(),NumSD());
	fElementLocTangent1 = 0.0;
	fElementLocTangent2.Dimension(NumElements(),NumSD());
	fElementLocTangent2 = 0.0;
	fElementLocTangent3.Dimension(NumElements(),NumSD());
	fElementLocTangent3 = 0.0;
	fElementLocGradEnh.Dimension(NumElements(),NumSD());
	fElementLocGradEnh = 0.0;
	//hardcode element edge number for hex element for now
	int numedges = 12;
	fElementLocEdgeIntersect.Dimension(NumElements(),numedges);
	fElementLocEdgeIntersect = 0.0;
	
	/** variable from time step to time step */
	fElementLocScalars.Dimension(NumElements(),kNUM_SCALAR_TERMS);	
	fElementLocScalars = 0.0;
	fElementLocScalars_last.Dimension(NumElements(),kNUM_SCALAR_TERMS);
	fElementLocScalars_last = 0.0;
	fElementLocSlipDir.Dimension(NumElements(),NumSD());
	fElementLocSlipDir = 0.0;
	fElementLocSlipDir_last.Dimension(NumElements(),NumSD());
	fElementLocSlipDir_last = 0.0;
	fElementLocMuDir.Dimension(NumElements(),NumSD());
	fElementLocMuDir = 0.0;
	fElementLocMuDir_last.Dimension(NumElements(),NumSD());
	fElementLocMuDir_last = 0.0;
	fElementLocInternalVars.Dimension(NumElements(),kNUM_ISV_TERMS);	
	fElementLocInternalVars = 0.0;
	fElementLocInternalVars_last.Dimension(NumElements(),kNUM_ISV_TERMS);
	fElementLocInternalVars_last = 0.0;
	fElementVolume.Dimension(NumElements());
	fElementVolume = 0.0;
	fElementVolume_last.Dimension(NumElements());
	fElementVolume_last = 0.0;
	/*@}*/
	
	/* localization element arrays */
	normals.Dimension(NumSD());
	slipdirs.Dimension(NumSD());
	
	normal1.Dimension(NumSD());
	normal2.Dimension(NumSD());
	normal3.Dimension(NumSD());
	normal_chosen.Dimension(NumSD());
	slipdir1.Dimension(NumSD());
	slipdir2.Dimension(NumSD());
	slipdir3.Dimension(NumSD());
	slipdir_chosen.Dimension(NumSD());
	tangent1.Dimension(NumSD());
	tangent2.Dimension(NumSD());
	tangent3.Dimension(NumSD());
	tangent_chosen.Dimension(NumSD());
	
	/* need to initialize previous volume */
	Top();
	while (NextElement())
	{
		/* inherited - computes gradients and standard 
		 * deformation gradients */
		SolidElementT::SetGlobalShape();

		/* compute mean of shape function gradients */
		double& vol = fElementVolume_last[CurrElementNumber()];
		SetMeanGradient(fMeanGradient, vol);
	}
	
	/* cohesive surface material parameters */
	fCohesiveSurface_Params.Dimension ( kNUM_CS_TERMS );
	fCohesiveSurface_Params[kc_r] = list.GetParameter("residual_cohesion");
	fCohesiveSurface_Params[kc_p] = list.GetParameter("peak_cohesion");
	fCohesiveSurface_Params[kalpha_c] = list.GetParameter("cohesion_softening_coefficient");
	fCohesiveSurface_Params[kphi_r] = list.GetParameter("residual_friction_angle_rad");
	fCohesiveSurface_Params[kphi_p] = list.GetParameter("peak_friction_angle_rad");
	fCohesiveSurface_Params[kalpha_phi] = list.GetParameter("friction_softening_coefficient");
	fCohesiveSurface_Params[kpsi_p] = list.GetParameter("peak_dilation_angle_rad");
	fCohesiveSurface_Params[kalpha_psi] = list.GetParameter("dilation_softening_coefficient");
	
	
	outputPrecision = 10;
	outputFileWidth = outputPrecision + 8;
	
	if (fFirstPass) 
	{
		ss_enh_out.open("ss_enh.info");
		fFirstPass = false;
	}
	else ss_enh_out.open_append("ss_enh.info");
	
}

/* extract the list of material parameters */
void SmallStrainEnhLocT::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{
	const char caller[] = "SmallStrainEnhLocT::CollectMaterialInfo";
	
	/* initialize */
	mat_params.Clear();
	
	/* collected material parameters */
	int num_blocks = all_params.NumLists("small_strain_enh_loc_element_block");
	for (int i = 0; i < num_blocks; i++) {

		/* block information */	
		const ParameterListT& block = all_params.GetList("small_strain_enh_loc_element_block", i);
		
		/* resolve material list name */
		if (i == 0) {
			const ParameterListT& mat_list_params = block.GetListChoice(*this, "small_strain_enh_loc_material_choice");
			mat_params.SetName(mat_list_params.Name());
		}
		
		/* collect material parameters */
		const ParameterListT& mat_list = block.GetList(mat_params.Name());
		const ArrayT<ParameterListT>& mat = mat_list.Lists();
		mat_params.AddList(mat[0]);
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct a new material support and return a pointer */
MaterialSupportT* SmallStrainEnhLocT::NewMaterialSupport(MaterialSupportT* p) const
{
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
}


/* return a pointer to a new material list */
MaterialListT* SmallStrainEnhLocT::NewMaterialList(const StringT& name, int size)
{
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
				ExceptionT::GeneralFail("SmallStrainEnhLocT::NewMaterialList");
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
}


/* choose the normal and slipdir given normals and slipdirs from bifurcation condition */
void SmallStrainEnhLocT::ChooseNormalAndSlipDir(void)
{
	/* current element number */
	int elem = CurrElementNumber();
	loc_flag = fElementLocScalars[elem,kLocFlag];
	
	/* fetch normals and slipdirs for element */
	fElementLocNormal1.RowAlias(elem, normal1);
	fElementLocNormal2.RowAlias(elem, normal2);
	fElementLocNormal3.RowAlias(elem, normal3);
	fElementLocTangent1.RowAlias(elem, tangent1);
	fElementLocTangent2.RowAlias(elem, tangent2);
	fElementLocTangent3.RowAlias(elem, tangent3);
	fElementLocSlipDir1.RowAlias(elem, slipdir1);
	fElementLocSlipDir2.RowAlias(elem, slipdir2);
	fElementLocSlipDir3.RowAlias(elem, slipdir3);
	
	/* loop through nodes as inner product of nodal displacements and each normal */
	/*
	fShapes->TopNode??();
	while (fShapes->NextNode??())
	{
		
	}	
	*/
	
	/*
	normal_chosen = ?;
	slipdir_chosen = ?;
	tangent_chosen = ?;
	*/
	
	/* store chosen normal and slip direction vectors */
	fElementLocNormal.SetRow(elem, normal_chosen);
	fElementLocSlipDir.SetRow(elem, slipdir_chosen);
	fElementLocTangent.SetRow(elem, tangent_chosen);
}

/* given the normal and one point, determine active nodes */
void SmallStrainEnhLocT::DetermineActiveNodes(void)
{
	/* current element number */
	int elem = CurrElementNumber();
	loc_flag = fElementLocScalars[elem,kLocFlag];
	
	/* fetch chosen normal */
	fElementLocNormal.RowAlias(elem, normal_chosen);
	
	/* loop through nodes and determine active nodes */
	
	/* loop through nodes and determine edge intersectons */
	/*
	for (int i = 0; i < numedges; i++)
	{
		//grab associated nodes with edges??
		fElementLocEdgeIntersect[elem,i] = ?;
	}
	*/
	
	/*
	see function SetGlobalShape to see how to grab Grad of shape function and assemble GradEnh
	knowing active nodes
	*/
}

/* calculate the internal force contribution ("-k*d") */
void SmallStrainEnhLocT::FormKd(double constK)
{
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	/* collect incremental heat */
	bool need_heat = fElementHeat.Length() == fShapes->NumIP();

	/* current element number */
	int elem = CurrElementNumber();
	loc_flag = fElementLocScalars[elem,kLocFlag];
	
	/* fetch normal and slipdir for element */
	fElementLocNormal.RowAlias(elem, normal_chosen);
	fElementLocSlipDir.RowAlias(elem, slipdir_chosen);
	fElementLocTangent.RowAlias(elem, tangent_chosen);
	
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
		
		if (loc_flag == 2) 
		{
			//modify stiffness matrix
			
		}
		
		/* incremental heat generation */
		if (need_heat) 
			fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	}	
}

/* form the element stiffness matrix */
void SmallStrainEnhLocT::FormStiffness(double constK)
{
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
	
	/* current element info */
	int elem = CurrElementNumber();
	loc_flag = fElementLocScalars[elem,kLocFlag];
	double vol = fElementVolume[elem];
	
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
		
		/*
		check for localization if element has not localized, or has localized
		but not yet traced (the stress state will change when elements are traced
		and the cohesive surface model is activated)
		*/
		if (loc_flag == 2) 
		{
			//modify stiffness matrix
			
		}
		else
		{
			//check for localization
			bool checkloc = fCurrMaterial->IsLocalized(normals,slipdirs);
			if (checkloc) 
			{
				loc_flag = 1;
				fElementLocScalars[elem,kLocFlag] = loc_flag;
				normals.Top();
				slipdirs.Top();
				int num_normals = normals.Length();
				if (num_normals == 2)
				{
					normals.Next();
					normal1 = normals.Current();
					fElementLocNormal1.SetRow(elem, normal1);
					normals.Next();
					normal2 = normals.Current();
					fElementLocNormal2.SetRow(elem, normal2);
					normal3 = 0.0;
					fElementLocNormal3.SetRow(elem, normal3);
					
					slipdirs.Next();
					slipdir1 = slipdirs.Current();
					fElementLocSlipDir1.SetRow(elem, slipdir1);
					slipdirs.Next();
					slipdir2 = slipdirs.Current();
					fElementLocSlipDir2.SetRow(elem, slipdir2);
					slipdir3 = 0.0;
					fElementLocSlipDir3.SetRow(elem, slipdir3);
					
					//calculate tangent
				}
				else if (num_normals == 3)
				{
					normals.Next();
					normal1 = normals.Current();
					fElementLocNormal1.SetRow(elem, normal1);
					normals.Next();
					normal2 = normals.Current();
					fElementLocNormal2.SetRow(elem, normal2);
					normals.Next();
					normal3 = normals.Current();
					fElementLocNormal3.SetRow(elem, normal3);
					
					slipdirs.Next();
					slipdir1 = slipdirs.Current();
					fElementLocSlipDir1.SetRow(elem, slipdir1);
					slipdirs.Next();
					slipdir2 = slipdirs.Current();
					fElementLocSlipDir2.SetRow(elem, slipdir2);
					slipdirs.Next();
					slipdir3 = slipdirs.Current();
					fElementLocSlipDir3.SetRow(elem, slipdir3);
				
					//calculate tangent
				}
				else
				{
					ExceptionT::GeneralFail("SmallStrainEnhLocT::FormStiffness - incorrect number of normals");
				}
				
				
				ss_enh_out	<< endl << "  " << endl << "normal1: " << setw(outputFileWidth) << normal1[0] 
							<< setw(outputFileWidth) << normal1[1] <<  setw(outputFileWidth) << normal1[2]
							<< setw(outputFileWidth) << "slipdir1: " << setw(outputFileWidth) << slipdir1[0] 
							<< setw(outputFileWidth) << slipdir1[1] <<  setw(outputFileWidth) << slipdir1[2]; 
				ss_enh_out	<< endl << "normal2: " << setw(outputFileWidth) << normal2[0] 
							<< setw(outputFileWidth) << normal2[1] <<  setw(outputFileWidth) << normal2[2]
							<< setw(outputFileWidth) << "slipdir2: " << setw(outputFileWidth) << slipdir2[0] 
							<< setw(outputFileWidth) << slipdir2[1] <<  setw(outputFileWidth) << slipdir2[2]; 
				ss_enh_out	<< endl << "normal3: " << setw(outputFileWidth) << normal3[0] 
							<< setw(outputFileWidth) << normal3[1] <<  setw(outputFileWidth) << normal3[2]
							<< setw(outputFileWidth) << "slipdir3: " << setw(outputFileWidth) << slipdir3[0] 
							<< setw(outputFileWidth) << slipdir3[1] <<  setw(outputFileWidth) << slipdir3[2]; 						


				//do band tracing somewhere else in the element, and only when the stress state is converged

			}
			
		}
	}
}

/* compute the measures of strain/deformation over the element */
void SmallStrainEnhLocT::SetGlobalShape(void)
{
	/* current element number */
	int elem = CurrElementNumber();
	
	/* inherited */
	SolidElementT::SetGlobalShape();

	/* compute mean of shape function gradients */
	double& vol = fElementVolume[elem];
	SetMeanGradient(fMeanGradient, vol);
	
	/* last deformed volume */
	double& vol_last = fElementVolume_last[elem];
	
	/* material information */
	int material_number = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[material_number];
	
	/* using B-bar */
	if (fStrainDispOpt == kMeanDilBbar)
	{
		/* compute mean of shape function gradients */
		//SetMeanGradient(fMeanGradient, vol);

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
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* compute mean shape function gradient, Hughes (4.5.23) */
void SmallStrainEnhLocT::SetMeanGradient(dArray2DT& mean_gradient, double& vol) const
{
	int nip = NumIP();
	const double* det = fShapes->IPDets();
	const double*   w = fShapes->IPWeights();

	/* volume */
	vol = 0.0;
	for (int i = 0; i < nip; i++)
		vol += w[i]*det[i];

	/* initialize */
	mean_gradient = 0.0;			

	/* integrate */
	for (int i = 0; i < nip; i++)
		mean_gradient.AddScaled(w[i]*det[i]/vol, fShapes->Derivatives_U(i));
}

