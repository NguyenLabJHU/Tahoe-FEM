/* $Id: SmallStrainEnhLocT.cpp,v 1.13 2005-03-09 19:27:51 raregue Exp $ */
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

//#include "BasicSupportT.h"

using namespace Tahoe;

/* parameters */
const double Pi = acos(-1.0);
const double smallnum = 1.0e-4;

/* initialize static variables */
bool SmallStrainEnhLocT::fFirstPass = true;
bool SmallStrainEnhLocT::fDeBug = true;
bool SmallStrainEnhLocT::fFirstTrace = false;

/* constructor */
SmallStrainEnhLocT::SmallStrainEnhLocT(const ElementSupportT& support):
	SolidElementT(support),
	fNeedsOffset(-1),
	fSSMatSupport(NULL),
	fK_dd(ElementMatrixT::kNonSymmetric),
	fK_dzeta(ElementMatrixT::kNonSymmetric),
	fK_zetad(ElementMatrixT::kNonSymmetric),
	displ_u(LocalArrayT::kDisp)
{
	SetName("small_strain_enh_loc");
}

/* destructor */
SmallStrainEnhLocT::~SmallStrainEnhLocT(void)
{
	delete fSSMatSupport;
}


/* finalize current step - step is solved */
void SmallStrainEnhLocT::InitStep(void)
{
	/* inherited */
	SolidElementT::InitStep();
	
	/* initialize to converged previous solution */
	fElementLocScalars = fElementLocScalars_last;
	fElementLocSlipDir = fElementLocSlipDir_last;
	fElementLocMuDir = fElementLocMuDir_last;
	fElementLocInternalVars = fElementLocInternalVars_last;
	fElementVolume = fElementVolume_last;
	
	if (fDeBug)
	{
		ss_enh_out	<< endl << "**********************************************************************************************";
		ss_enh_out	<< endl << "**********************************************************************************************" << endl;
		ss_enh_out	<< endl 
					<< setw(outputFileWidth) << "time_step"
					<< endl;
		int step = fSSMatSupport->StepNumber();
		ss_enh_out	<< setw(outputFileWidth) << step
					<< endl;	
	}
			
}


/* finalize current step - step is solved */
void SmallStrainEnhLocT::CloseStep(void)
{
	/* inherited */
	SolidElementT::CloseStep();
	
	int nen, elem_num;
	
	Top();
	while (NextElement())
	{
		/* shape function derivatives, jacobians, local coords */
		SetGlobalShape();
		
		/* get displacements */
        SetLocalU(fLocDisp);
        
        /* get nodal coordinates */
        SetLocalX(fLocInitCoords);
		
		elem_num = CurrElementNumber();
		nen = NumElementNodes();
		loc_flag = fElementLocScalars[kNUM_SCALAR_TERMS*elem_num + kLocFlag];
			
		if ( fabs(loc_flag - 1.0) < smallnum )
		{
			/* choose normal and slip direction based on current element deformation */
			ChooseNormalAndSlipDir(fLocDisp, elem_num, nen);
		}
	
		if ( fabs(loc_flag - 1.0) < smallnum )
		{
			/* determine active nodes and band trace */
			DetermineActiveNodesTrace(fLocInitCoords, elem_num, nen);
		}
		
	} // while next element
	
	/* store converged solution */
	fElementLocScalars_last = fElementLocScalars;
	fElementLocSlipDir_last = fElementLocSlipDir;
	fElementLocMuDir_last = fElementLocMuDir;
	fElementLocInternalVars_last = fElementLocInternalVars;
	fElementVolume_last = fElementVolume;
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

/* see if pair status has changed */
GlobalT::RelaxCodeT SmallStrainEnhLocT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT code = SolidElementT::RelaxSystem();
	
	int nen, elem_num;
	
	Top();
	while (NextElement())
	{
		/* shape function derivatives, jacobians, local coords */
		SetGlobalShape();
		
		elem_num = CurrElementNumber();
		nen = NumElementNodes();
		loc_flag = fElementLocScalars[kNUM_SCALAR_TERMS*elem_num + kLocFlag];
		const ElementCardT& element = CurrentElement();

		/* until element is traced, i.e. loc_flag = 2, check for localization */
		if ( loc_flag < 2.0 && element.IsAllocated() )
		{
			/* initialize element localization data */
			fElementLocScalars[kNUM_SCALAR_TERMS*elem_num + kLocFlag] = 0.0;
			
			fElementLocScalars[kNUM_SCALAR_TERMS*elem_num + kdetAmin] = 1.0e99;
			
			fElementLocScalars[kNUM_SCALAR_TERMS*elem_num + kdissip_max] = 0.0;
			
			/*
			normal1 = 0.0;
			fElementLocNormal1.SetRow(elem_num, normal1);
			normal2 = 0.0;
			fElementLocNormal2.SetRow(elem_num, normal2);
			normal3 = 0.0;
			fElementLocNormal3.SetRow(elem_num, normal3);
			
			slipdir1 = 0.0;
			fElementLocSlipDir1.SetRow(elem_num, slipdir1);
			slipdir2 = 0.0;
			fElementLocSlipDir2.SetRow(elem_num, slipdir2);
			slipdir3 = 0.0;
			fElementLocSlipDir3.SetRow(elem_num, slipdir3);
	
			tangent1 = 0.0;
			fElementLocTangent1.SetRow(elem_num, tangent1);
			tangent2 = 0.0;
			fElementLocTangent2.SetRow(elem_num, tangent2);
			tangent3 = 0.0;
			fElementLocTangent3.SetRow(elem_num, tangent3);
			*/
			
			normal_chosen = 0.0;
			fElementLocNormal.SetRow(elem_num, normal_chosen);
			slipdir_chosen = 0.0;
			fElementLocSlipDir.SetRow(elem_num, slipdir_chosen);
			tangent_chosen = 0.0;
			fElementLocTangent.SetRow(elem_num, tangent_chosen);

			/* check for localization */
			CheckLocalization(elem_num);
		}
	
	} // while next element

	return code;
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
	// there are 3 possible normals, so 3 possible psi's
	fElementLocPsiSet.Dimension(NumElements(),3);
	fElementLocPsiSet = 0.0;
	fElementLocGradEnh.Dimension(NumElements(),NumSD()*NumIP());
	fElementLocGradEnh = 0.0;
	fElementLocGradEnhIP.Dimension(NumIP(),NumSD());
	fElementLocGradEnhIP = 0.0;
	//hardcode element edge number for hex element for now
	int numedges = 12;
	fElementLocEdgeIntersect.Dimension(NumElements(),numedges);
	fElementLocEdgeIntersect = 0.0;
	
	fElementLocStartSurface.Dimension(NumElements(),NumSD());
	fElementLocStartSurface = 0.0;
	
	fElementLocNodesActive.Dimension(NumElements(),NumElementNodes());
	fElementLocNodesActive = 0;
	
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
	
	mu_dir.Dimension(NumSD());
	
	start_surface_vect.Dimension(NumSD());
	
	start_surface_vect_read.Dimension(NumSD());
	
	grad_enh.Dimension(NumSD()*NumIP());
	
	q_isv.Dimension(kNUM_ISV_TERMS);
	h_q.Dimension(kNUM_ISV_TERMS);
	DqDzeta.Dimension(kNUM_ISV_TERMS);
	DhqDzeta.Dimension(kNUM_ISV_TERMS);
	
	DslipdirDzeta.Dimension(NumSD());
	
	F_mun.Dimension(NumSD());
	G_enh.Dimension(NumSD());
	F_nn.Dimension(NumSD());
	
	fDe.Dimension(dSymMatrixT::NumValues(NumSD()));
	
	fK_dd.Dimension(NumElementNodes(),NumElementNodes());
	int dum = 1;
	fK_dzeta.Dimension(NumElementNodes(),dum);
	fK_zetad.Dimension(dum,NumElementNodes());
	
	displ_u.Dimension (NumElementNodes(), NumSD());
	//fDispl->RegisterLocal(displ_u);
	
	node_displ.Dimension(NumSD());
	node_displ = 0.0;
	node_coords.Dimension(NumSD());
	node_coords = 0.0;
	node_coords.Dimension(NumSD());
	node_shape_deriv = 0.0;
	
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
	
	
	/* create output file for debugging */
	if (fDeBug)
	{
		outputPrecision = 10;
		outputFileWidth = outputPrecision + 8;
		
		if (fFirstPass) 
		{
			ss_enh_out.open("ss_enh.info");
			fFirstPass = false;
		}
		else ss_enh_out.open_append("ss_enh.info");
	}
	
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


/* check for localization */
void SmallStrainEnhLocT::CheckLocalization(int& elem)
{
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		/* determine localization */
		//bool checkloc = fCurrMaterial->IsLocalized(normals,slipdirs);
		//bool checkloc = fCurrMaterial->IsLocalized(normals,slipdirs,detAmin);
		bool checkloc = fCurrMaterial->IsLocalized(normals,slipdirs,detAs,dissipations_fact);
		
		/* check if localized and if has minimum determinant for IPs */
		//if (checkloc && detAmin < fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdetAmin])
		/* check if localized */
		if (checkloc)
		{
			/* determine minimum determinant among normals at different IPs */
			detAmin = 1.0e99;
			detAs.Top();
			while (detAs.Next())
			{
				detA_tmp = detAs.Current();
				if (detA_tmp < detAmin) detAmin = detA_tmp;
			}
			
			/* determine maximum dissipation among normals at different IPs */
			dissip_max = 0.0;
			dissipations_fact.Top();
			while (dissipations_fact.Next())
			{
				dissip_tmp = dissipations_fact.Current();
				if (dissip_tmp > dissip_max) dissip_max = dissip_tmp;
			}
			
			/* if maximum dissipation, do further localization calcs */
			//if (dissip_max> fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdissip_max])
			/* if minimum determinant, do further localization calcs */
			if (detAmin < fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdetAmin])
			{
				if (fDeBug)
				{
					const int ip = fShapes->CurrIP()+1;
					ss_enh_out	<< endl 
								<< setw(outputFileWidth) << "element " << elem 
								<< setw(outputFileWidth) <<  "IP" << ip;
				}
				
				loc_flag = 1.0; // localized, not traced
				fElementLocScalars[kNUM_SCALAR_TERMS*elem + kLocFlag] = loc_flag;
				
				fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdetAmin] = detAmin;
				
				fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdissip_max] = dissip_max;
				
				/* calculate tangent unit vectors and write output for debugging */
				normals.Top();
				slipdirs.Top();
				tangents.Free();
				detAs.Top();
				dissipations_fact.Top();
				int num_normals = normals.Length();
				dArrayT dummyt;
				while (normals.Next())
				{
					normal1 = normals.Current();
					slipdirs.Next();
					slipdir1 = slipdirs.Current();
					double product = dArrayT::Dot(normal1,slipdir1);
					psi1 = asin(product);
					fElementLocPsiSet[elem,0] = psi1;
					double sec = 1.0/cos(psi1);
					if (fabs(psi1-Pi/2.0) < smallnum)
					{
						tangent1 = 0.0;
					}
					else 
					{
						tangent1 = slipdir1;
						tangent1 *= sec;
						dummyt = normal1;
						dummyt *= tan(psi1);
						tangent1 -= dummyt;
					}
					tangents.Append(tangent1);
					if (fDeBug)
					{
						detAs.Next();
						detA_tmp = detAs.Current();
						dissipations_fact.Next();
						dissip_tmp = dissipations_fact.Current();
						ss_enh_out	<< endl << "detA: " << setw(outputFileWidth) << detA_tmp << setw(outputFileWidth)
									<< "dissipation: " << setw(outputFileWidth) << dissip_tmp; 			
						ss_enh_out	<< endl << "normal: " << setw(outputFileWidth) << normal1[0] 
									<< setw(outputFileWidth) << normal1[1] <<  setw(outputFileWidth) << normal1[2]
									<< setw(outputFileWidth) << "slipdir: " << setw(outputFileWidth) << slipdir1[0] 
									<< setw(outputFileWidth) << slipdir1[1] <<  setw(outputFileWidth) << slipdir1[2]; 
						ss_enh_out	<< endl << "tangent:" << setw(outputFileWidth) << tangent1[0] 
									<< setw(outputFileWidth) << tangent1[1] <<  setw(outputFileWidth) << tangent1[2];
					}
				} // while normals.Next
				
				/* assign current normals, slipdirs, tangents, and detAs 
				 * as values associated with minimum determinant for element */
				normals_min.Free();
				slipdirs_min.Free();
				tangents_min.Free();
				detAs_min.Free();
				dissipations_fact_min.Free();
				normals_min = normals;
				slipdirs_min = slipdirs;
				tangents_min = tangents;
				detAs_min = detAs;
				dissipations_fact_min = dissipations_fact;
				
			
			
			/*
			if (num_normals == 1)
			{
				normals.Next();
				normal1 = normals.Current();
				fElementLocNormal1.SetRow(elem, normal1);
				normal2 = 0.0;
				fElementLocNormal2.SetRow(elem, normal2);
				normal3 = 0.0;
				fElementLocNormal3.SetRow(elem, normal3);
				
				slipdirs.Next();
				slipdir1 = slipdirs.Current();
				fElementLocSlipDir1.SetRow(elem, slipdir1);
				slipdir2 = 0.0;
				fElementLocSlipDir2.SetRow(elem, slipdir2);
				slipdir3 = 0.0;
				fElementLocSlipDir3.SetRow(elem, slipdir3);
				
				//calculate tangents
				double product = dArrayT::Dot(normal1,slipdir1);
				psi1 = asin(product);
				fElementLocPsiSet[elem,0] = psi1;
				double sec = 1.0/cos(psi1);
				if (fabs(psi1-Pi/2.0) < smallnum)
				{
					tangent1 = 0.0;
				}
				else 
				{
					tangent1 = slipdir1;
					tangent1 *= sec;
					dummyt = normal1;
					dummyt *= tan(psi1);
					tangent1 -= dummyt;
				}
				fElementLocTangent1.SetRow(elem, tangent1);
				
				tangent2 = 0.0;
				fElementLocTangent2.SetRow(elem, tangent2);
				
				tangent3 = 0.0;
				fElementLocTangent3.SetRow(elem, tangent3);
			}
			else if (num_normals == 2)
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
				
				//calculate tangents
				double product = dArrayT::Dot(normal1,slipdir1);
				psi1 = asin(product);
				fElementLocPsiSet[elem,0] = psi1;
				double sec = 1.0/cos(psi1);
				if (fabs(psi1-Pi/2.0) < smallnum)
				{
					tangent1 = 0.0;
				}
				else 
				{
					tangent1 = slipdir1;
					tangent1 *= sec;
					dummyt = normal1;
					dummyt *= tan(psi1);
					tangent1 -= dummyt;
				}
				fElementLocTangent1.SetRow(elem, tangent1);
				
				product = dArrayT::Dot(normal2,slipdir2);
				psi2 = asin(product);
				fElementLocPsiSet[elem,1] = psi2;
				sec = 1.0/cos(psi2);
				if (fabs(psi2-Pi/2.0) < smallnum) 
				{
					tangent2 = 0.0;
				}
				else 
				{
					tangent2 = slipdir2;
					tangent2 *= sec;
					dummyt = normal2;
					dummyt *= tan(psi2);
					tangent2 -= dummyt;
				}
				fElementLocTangent2.SetRow(elem, tangent2);
				
				tangent3 = 0.0;
				fElementLocTangent3.SetRow(elem, tangent3);
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
				double product = dArrayT::Dot(normal1,slipdir1);
				psi1 = asin(product);
				fElementLocPsiSet[elem,0] = psi1;
				double sec = 1.0/cos(psi1);
				if (fabs(psi1-Pi/2.0) < smallnum)
				{
					tangent1 = 0.0;
				}
				else 
				{
					tangent1 = slipdir1;
					tangent1 *= sec;
					dummyt = normal1;
					dummyt *= tan(psi1);
					tangent1 -= dummyt;
				}
				fElementLocTangent1.SetRow(elem, tangent1);
				
				product = dArrayT::Dot(normal2,slipdir2);
				psi2 = asin(product);
				fElementLocPsiSet[elem,1] = psi2;
				sec = 1.0/cos(psi2);
				if (fabs(psi2-Pi/2.0) < smallnum) 
				{
					tangent2 = 0.0;
				}
				else 
				{
					tangent2 = slipdir2;
					tangent2 *= sec;
					dummyt = normal2;
					dummyt *= tan(psi2);
					tangent2 -= dummyt;
				}
				fElementLocTangent2.SetRow(elem, tangent2);
				
				product = dArrayT::Dot(normal3,slipdir3);
				psi3 = asin(product);
				fElementLocPsiSet[elem,2] = psi3;
				// doesn't work
				//if (fabs(psi3-Pi/2.0) < 1.0e-3) tangent3 = 0.0;
				//else tangent3 = slipdir3/cos(psi3) - tan(psi3)*normal3;
				sec = 1.0/cos(psi3);
				if (fabs(psi3-Pi/2.0) < smallnum) 
				{
					tangent3 = 0.0;
				}
				else 
				{
					tangent3 = slipdir3;
					tangent3 *= sec;
					dummyt = normal3;
					dummyt *= tan(psi3);
					tangent3 -= dummyt;
				}
				fElementLocTangent3.SetRow(elem, tangent3);
			}
			else
			{
				ExceptionT::GeneralFail("SmallStrainEnhLocT::FormStiffness - incorrect number of normals");
			}
		
			
			if (fDeBug)
			{
				const int ip = fShapes->CurrIP()+1;
				ss_enh_out	<< endl 
							<< setw(outputFileWidth) << "element " << elem 
							<< setw(outputFileWidth) <<  "IP" << ip;
				ss_enh_out	<< endl << "detAmin: " << setw(outputFileWidth) << detAmin; 			
				ss_enh_out	<< endl << "normal1: " << setw(outputFileWidth) << normal1[0] 
							<< setw(outputFileWidth) << normal1[1] <<  setw(outputFileWidth) << normal1[2]
							<< setw(outputFileWidth) << "slipdir1: " << setw(outputFileWidth) << slipdir1[0] 
							<< setw(outputFileWidth) << slipdir1[1] <<  setw(outputFileWidth) << slipdir1[2]; 
				ss_enh_out	<< endl << "tangent1:" << setw(outputFileWidth) << tangent1[0] 
							<< setw(outputFileWidth) << tangent1[1] <<  setw(outputFileWidth) << tangent1[2];
				ss_enh_out	<< endl << "normal2: " << setw(outputFileWidth) << normal2[0] 
							<< setw(outputFileWidth) << normal2[1] <<  setw(outputFileWidth) << normal2[2]
							<< setw(outputFileWidth) << "slipdir2: " << setw(outputFileWidth) << slipdir2[0] 
							<< setw(outputFileWidth) << slipdir2[1] <<  setw(outputFileWidth) << slipdir2[2]; 
				ss_enh_out	<< endl << "tangent2:" << setw(outputFileWidth) << tangent2[0] 
							<< setw(outputFileWidth) << tangent2[1] <<  setw(outputFileWidth) << tangent2[2];
				ss_enh_out	<< endl << "normal3: " << setw(outputFileWidth) << normal3[0] 
							<< setw(outputFileWidth) << normal3[1] <<  setw(outputFileWidth) << normal3[2]
							<< setw(outputFileWidth) << "slipdir3: " << setw(outputFileWidth) << slipdir3[0] 
							<< setw(outputFileWidth) << slipdir3[1] <<  setw(outputFileWidth) << slipdir3[2];
				ss_enh_out	<< endl << "tangent3:" << setw(outputFileWidth) << tangent3[0] 
							<< setw(outputFileWidth) << tangent3[1] <<  setw(outputFileWidth) << tangent3[2] << endl;
			}
			*/
			
			} // if (detAmin < fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdetAmin])

		} // if (checkloc)
		
	} // while next IP
	
}


/* choose the normal and slipdir given normals and slipdirs from bifurcation condition */
void SmallStrainEnhLocT::ChooseNormalAndSlipDir(LocalArrayT& displ_elem, int& elem, int& nen)
{
	int i, nodeindex;
	double product, sum1, sum2, sum3;

	/* fetch normals and slipdirs for element */
	/*
	fElementLocNormal1.RowCopy(elem, normal1);
	fElementLocNormal2.RowCopy(elem, normal2);
	fElementLocNormal3.RowCopy(elem, normal3);
	fElementLocTangent1.RowCopy(elem, tangent1);
	fElementLocTangent2.RowCopy(elem, tangent2);
	fElementLocTangent3.RowCopy(elem, tangent3);
	fElementLocSlipDir1.RowCopy(elem, slipdir1);
	fElementLocSlipDir2.RowCopy(elem, slipdir2);
	fElementLocSlipDir3.RowCopy(elem, slipdir3);
	*/
	
	/* determine normal by using minimum determinant associated with normal */
	normals_min.Top();
	slipdirs_min.Top();
	tangents_min.Top();
	detAs_min.Top();
	detAmin = 1.0e99;
	while (normals_min.Next())
	{
		detAs_min.Next();
		detA_tmp = detAs_min.Current();
		if (detA_tmp < detAmin)
		{
			detAmin = detA_tmp;
			normal_chosen = normals_min.Current();
			slipdirs_min.Next();
			slipdir_chosen = slipdirs_min.Current();
			tangents_min.Next();
			tangent_chosen = tangents_min.Current();
		}
	}

	/* determine normal based on minimum dissipation calculation */
	/*
	normals_min.Top();
	slipdirs_min.Top();
	tangents_min.Top();
	dissipations_fact_min.Top();
	dissip_max = 0.0;
	while (normals_min.Next())
	{
		dissipations_fact_min.Next();
		dissip_tmp = dissipations_fact_min.Current();
		if (dissip_tmp > dissip_max)
		{
			dissip_max = dissip_tmp;
			normal_chosen = normals_min.Current();
			slipdirs_min.Next();
			slipdir_chosen = slipdirs_min.Current();
			tangents_min.Next();
			tangent_chosen = tangents_min.Current();
		}
	}
	*/

	/* determines normal for shear loading,
	 * which is not applicable for cracking under pure tension, 
	 * or compaction banding under compression */
	//comment out for now
	/*
	sum1 = 0.0;
	sum2 = 0.0;
	sum3 = 0.0;
	for (i=0; i < nen; i++)
	{
		if (NumSD() == 2)
		{
			node_displ[0] = displ_elem[i];
			node_displ[1] = displ_elem[i+nen];
		}
		else if (NumSD() == 3)
		{
			node_displ[0] = displ_elem[i];
			node_displ[1] = displ_elem[i+nen];
			node_displ[2] = displ_elem[i+2*nen];
		}
		
		//while normals.Next
		product = dArrayT::Dot(normal1,node_displ);
		sum1 += product;
		product = dArrayT::Dot(normal2,node_displ);
		sum2 += product;
		product = dArrayT::Dot(normal3,node_displ);
		sum3 += product;
	}
	
	if ( normal3.Magnitude() > smallnum && normal2.Magnitude() > smallnum )
	{
		if ( fabs(sum1) < fabs(sum2) && fabs(sum1) < fabs(sum3))
		{
			normal_chosen = normal1;
			slipdir_chosen = slipdir1;
			tangent_chosen = tangent1;
		}
		else if ( fabs(sum2) < fabs(sum1) && fabs(sum2) < fabs(sum3))
		{
			normal_chosen = normal2;
			slipdir_chosen = slipdir2;
			tangent_chosen = tangent2;
		}
		else if ( fabs(sum3) < fabs(sum1) && fabs(sum3) < fabs(sum2))
		{
			normal_chosen = normal3;
			slipdir_chosen = slipdir3;
			tangent_chosen = tangent3;
		}
		else
		{
			normal_chosen = normal1;
			slipdir_chosen = slipdir1;
			tangent_chosen = tangent1;
		}
	}
	else if ( normal2.Magnitude() > smallnum )
	{
		if ( fabs(sum1) < fabs(sum2) )
		{
			normal_chosen = normal1;
			slipdir_chosen = slipdir1;
			tangent_chosen = tangent1;
		}
		else if ( fabs(sum2) < fabs(sum1) )
		{
			normal_chosen = normal2;
			slipdir_chosen = slipdir2;
			tangent_chosen = tangent2;
		}
		else
		{
			normal_chosen = normal1;
			slipdir_chosen = slipdir1;
			tangent_chosen = tangent1;
		}
	}
	else
	{
		normal_chosen = normal1;
		slipdir_chosen = slipdir1;
		tangent_chosen = tangent1;
	}
	*/
	
	/* store chosen normal and slip direction vectors */
	fElementLocNormal.SetRow(elem, normal_chosen);
	fElementLocSlipDir.SetRow(elem, slipdir_chosen);
	fElementLocTangent.SetRow(elem, tangent_chosen);
	
}

/* given the normal and one point, determine active nodes */
void SmallStrainEnhLocT::DetermineActiveNodesTrace(LocalArrayT& coords_elem, int& elem, int& nen)
{
	/* fetch chosen normal */
	fElementLocNormal.RowCopy(elem, normal_chosen);
	
	/* fetch slip surface starting intersection point */
	if (!fFirstTrace)
	{
		fElementLocStartSurface.SetRow(elem, start_surface_vect_read);
		fFirstTrace = true;
	}
	else
	{
		
	}
	
	fElementLocStartSurface.RowCopy(elem, start_surface_vect);
	/*
	read first start point from input file;
	calculate others based on slip surface intersection with adjacent elements
	*/
	
	dArrayT diff_vector;
	double product;
	
	/* loop through nodes and determine active nodes */
	for (int i=0; i < nen; i++)
	{
		if (NumSD() == 2)
		{
			node_coords[0] = coords_elem[i];
			node_coords[1] = coords_elem[i+nen];
			
			diff_vector.DiffOf(node_coords,start_surface_vect);
			product = dArrayT::Dot(normal_chosen,diff_vector);
			if (product > 0.0)
			{
				fElementLocNodesActive[NumElementNodes()*elem + i] = 1;
				//calculate enhancement in FormKd
			}
		}
		else if (NumSD() == 3)
		{
			node_coords[0] = coords_elem[i];
			node_coords[1] = coords_elem[i+nen];
			node_coords[2] = coords_elem[i+2*nen];
		}
		
	}
	
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
	
	/*
	if (trace)
	{
		loc_flag = 2.0;
		fElementLocScalars[kNUM_SCALAR_TERMS*elem + kLocFlag] = loc_flag;
	}
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
	loc_flag = fElementLocScalars[kNUM_SCALAR_TERMS*elem + kLocFlag];
	int nen = NumElementNodes();
	
	/* fetch normal and slipdir for element */
	fElementLocNormal.RowCopy(elem, normal_chosen);
	fElementLocSlipDir.RowCopy(elem, slipdir_chosen);
	fElementLocTangent.RowCopy(elem, tangent_chosen);
	
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
		
		if ( fabs(loc_flag - 2.0) < smallnum ) 
		{
			const int ip = fShapes->CurrIP(); 
			const dArray2DT& DNa = fShapes->Derivatives_U();

			/* loop through nodes and calculate enhancement function */
			for (int i=0; i < nen; i++)
			{
				DNa.RowCopy(i,node_shape_deriv);
				if ( fElementLocNodesActive[NumElementNodes()*elem + i] == 1 )
				{
					grad_enh += node_shape_deriv;
				}		
			}
			fElementLocGradEnhIP.SetRow(ip, grad_enh);
			
			//modify stiffness matrix
			
		}
		
		/* incremental heat generation */
		if (need_heat) 
			fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	}	
	
	fElementLocGradEnh.SetRow(elem, fElementLocGradEnhIP);

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
	loc_flag = fElementLocScalars[kNUM_SCALAR_TERMS*elem + kLocFlag];
	double vol = fElementVolume[elem];

	/* element has localized and has been traced, thus fetch data to modify the stiffness matrix */
	//if ( fabs(loc_flag - 2.0) < smallnum ) 
	if ( fabs(loc_flag - 1.0) < smallnum ) 
	{
		/* fetch normal and slipdir for element */
		fElementLocNormal.RowCopy(elem, normal_chosen);
		fElementLocSlipDir.RowCopy(elem, slipdir_chosen);
		fElementLocTangent.RowCopy(elem, tangent_chosen);
		fElementLocMuDir.RowCopy(elem, mu_dir);
		
		fElementLocGradEnh.RowCopy(elem, grad_enh);
	
		if (fDeBug)
		{
			ss_enh_out	<< endl << endl << "----------------------------------------------------------------------------------------------" << endl;
			ss_enh_out	<< endl 
						<< setw(outputFileWidth) << "element " << elem;
						
			ss_enh_out	<< endl << endl << "element volume:" << setw(outputFileWidth) << vol; 
						
			ss_enh_out	<< endl << endl << "detA_min:" << setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdetAmin]
						<< setw(outputFileWidth) << "dissip_max:" << setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kdissip_max]; 
			ss_enh_out	<< endl << " normal_chosen: " << setw(outputFileWidth) << normal_chosen[0] 
						<< setw(outputFileWidth) << normal_chosen[1] <<  setw(outputFileWidth) << normal_chosen[2]; 
			ss_enh_out	<< endl << "slipdir_chosen: " << setw(outputFileWidth) << slipdir_chosen[0] 
						<< setw(outputFileWidth) << slipdir_chosen[1] <<  setw(outputFileWidth) << slipdir_chosen[2]; 
			ss_enh_out	<< endl << "tangent_chosen: " << setw(outputFileWidth) << tangent_chosen[0] 
						<< setw(outputFileWidth) << tangent_chosen[1] <<  setw(outputFileWidth) << tangent_chosen[2]; 
			ss_enh_out	<< endl << "mu_dir: " << setw(outputFileWidth) << mu_dir[0] 
						<< setw(outputFileWidth) << mu_dir[1] <<  setw(outputFileWidth) << mu_dir[2];
			
			ss_enh_out	<< endl << endl << "loc_flag" << setw(outputFileWidth) << "jump_displ" 
						<< setw(outputFileWidth) << "gamma_delta" <<  setw(outputFileWidth) << "Q"
						<< setw(outputFileWidth) << "P" <<  setw(outputFileWidth) << "q_St"
						<< setw(outputFileWidth) << "q_Sn" <<  setw(outputFileWidth) << "p_S"; 
			ss_enh_out	<< endl << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kLocFlag] << setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kJumpDispl] 
						<< setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kgamma_delta] <<  setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kQ]
						<< setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kP] <<  setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kq_St]
						<< setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kq_Sn] <<  setw(outputFileWidth) << fElementLocScalars[kNUM_SCALAR_TERMS*elem + kp_S]; 
			
			ss_enh_out	<< endl << endl << "cohesion" << setw(outputFileWidth) << "friction (rad)" 
						<< setw(outputFileWidth) << "dilation (rad)"; 
			ss_enh_out	<< endl << fElementLocInternalVars[kNUM_ISV_TERMS*elem + kCohesion] << setw(outputFileWidth) << fElementLocInternalVars[kNUM_ISV_TERMS*elem + kFriction] 
						<< setw(outputFileWidth) << fElementLocInternalVars[kNUM_ISV_TERMS*elem + kDilation] << endl;		
		}
	}
	
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
		
		/* get De matrix */
		fDe.SetToScaled(scale, fCurrMaterial->ce_ijkl());
		if (print) cout << "\nelas_modulus: "<<fCurrMaterial->ce_ijkl();
							
		/* multiply b(transpose) * db, taking account of symmetry, */
		/* and accumulate in elstif */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
		
		/*
		check for localization if element has not localized, or has localized
		but not yet traced
		*/
		if ( fabs(loc_flag - 2.0) < smallnum ) 
		{
			//modify stiffness matrix
			
			const int ip = fShapes->CurrIP()+1;
			int array_location = fShapes->CurrIP()*NumSD();
			fElementLocGradEnhIP[CurrIP(),0] = grad_enh[array_location];
			fElementLocGradEnhIP[CurrIP(),1] = grad_enh[array_location+1];
			fElementLocGradEnhIP[CurrIP(),2] = grad_enh[array_location+2];
			if (fDeBug)
			{
				ss_enh_out	<< endl << "GradEnh at IP" << ip
							<< setw(outputFileWidth) << fElementLocGradEnhIP[CurrIP(),0] 
							<< setw(outputFileWidth) << fElementLocGradEnhIP[CurrIP(),1] 
							<< setw(outputFileWidth) << fElementLocGradEnhIP[CurrIP(),2];
			}
		}
		
	} // while next IP
	
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

