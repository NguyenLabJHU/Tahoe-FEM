/* $Id: FSSolidMatList2DT.cpp,v 1.4 2004-08-01 20:42:03 paklein Exp $ */
#include "FSSolidMatList2DT.h"
#include "FSMatSupportT.h"


#include "SolidMaterialsConfig.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#endif

#include "FDHookeanMat2DT.h"
#include "FDKStV2D.h"
#include "FDCubic2DT.h"
#include "SimoIso2D.h"
#include "QuadLog2D.h"
#include "QuadLogOgden2DT.h"

#ifdef CAUCHY_BORN_MATERIAL
#include "EAMFCC2D.h"
#include "LJTr2D.h"
#include "Hex2D.h"
#endif

#ifdef MODCBSW_MATERIAL
#include "ModCB2DT.h"
#endif

#ifdef VIB_MATERIAL
#include "IsoVIB2D.h"
#include "J2IsoVIB2DLinHardT.h"
#include "VIB2D.h"
#include "D2VIB2D_a.h"
#include "OgdenIsoVIB2D.h"
#endif

#ifdef PLASTICITY_MACRO_MATERIAL
#include "HyperEVP2D.h"
#include "BCJHypo2D.h"
#include "BCJHypoIsoDamageKE2D.h"
#include "BCJHypoIsoDamageYC2D.h"
#endif

#ifdef PLASTICITY_CRYSTAL_MATERIAL
#include "LocalCrystalPlast2D.h"
#include "GradCrystalPlast2D.h"
#include "LocalCrystalPlastFp2D.h"
#include "GradCrystalPlastFp2D.h"
#endif

#ifdef PLASTICITY_J2_MATERIAL
#include "J2Simo2D.h"
#include "J2QL2DLinHardT.h"
#endif

#ifdef VISCOELASTIC_MATERIALS_DEV
#include "RGVIB2D.h"
#include "RGSplit3D.h"
#include "FDSV_KStV2D.h"
#endif

#ifdef VISCOELASTICITY
#include "RGSplitT.h"
#endif

#ifdef ELASTIC_OGDEN_MATERIAL_DEV
#include "OgdenMaterialT.h"
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
#include "ABAQUS_BCJ.h"
#include "ABAQUS_BCJ_ISO.h"
#include "ABAQUS_VUMAT_BCJ.h"
#endif
#endif

#ifdef THERMO_VISCO_PLASTIC_MATERIAL
#include "tevp2D.h"
#include "povirk2D.h"
#endif

using namespace Tahoe;

/* constructor */
FSSolidMatList2DT::FSSolidMatList2DT(int length, const FSMatSupportT& support):
	SolidMatListT(length, support),
	fFSMatSupport(&support)
{
#pragma message("check spatial dimension of material support")
	SetName("large_strain_material_2D");
}

FSSolidMatList2DT::FSSolidMatList2DT(void):
	fFSMatSupport(NULL)
{
	SetName("large_strain_material_2D");
}

/* return true if the list contains plane stress models */
bool FSSolidMatList2DT::HasPlaneStress(void) const
{
	/* check materials */
	for (int i = 0; i < Length(); i++)
	{
		/* get pointer to Material2DT */
		const ContinuumMaterialT* cont_mat = fArray[i];
		const SolidMaterialT* sol_mat = TB_DYNAMIC_CAST(const SolidMaterialT*, cont_mat);
		
		/* assume materials that don't have Material2DT are plane strain */
		if (sol_mat && sol_mat->Constraint() == SolidMaterialT::kPlaneStress) 
			return true;
	}
	return false;
}

/* information about subordinate parameter lists */
void FSSolidMatList2DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("fs_material_list_2D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void FSSolidMatList2DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "fs_material_list_2D")
	{
		order = ParameterListT::Choice;
	
		sub_lists.AddSub("large_strain_Hookean_2D");
		sub_lists.AddSub("large_strain_cubic_2D");
		sub_lists.AddSub("large_strain_StVenant_2D");
		sub_lists.AddSub("Simo_isotropic_2D");
		sub_lists.AddSub("quad_log_2D");
		sub_lists.AddSub("quad_log_Ogden_2D");
#ifdef PLASTICITY_J2_MATERIAL
		sub_lists.AddSub("Simo_J2_2D");
		sub_lists.AddSub("quad_log_J2_2D");
#endif

#ifdef CAUCHY_BORN_MATERIAL
		sub_lists.AddSub("LJ_triangular_2D");
		sub_lists.AddSub("hex_2D");
		sub_lists.AddSub("FCC_EAM_2D");
#endif

#ifdef MODCBSW_MATERIAL
		sub_lists.AddSub("Cauchy-Born_diamond_2D");
#endif

#ifdef VIB_MATERIAL
		sub_lists.AddSub("VIB_2D");
		sub_lists.AddSub("isotropic_VIB_2D");
		sub_lists.AddSub("Ogden_isotropic_VIB_2D");
#endif

#ifdef VISCOELASTICITY
		sub_lists.AddSub("Reese-Govindjee_split");
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
		sub_lists.AddSub("ABAQUS_UMAT_BCJ");
		sub_lists.AddSub("ABAQUS_VUMAT_BCJ");
		sub_lists.AddSub("ABAQUS_UMAT_BCJ_iso-damage");
#endif
#endif
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSSolidMatList2DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	FSSolidMatT* fs_solid_mat = NewFSSolidMat(name);
	if (fs_solid_mat)
		return fs_solid_mat;
	else /* inherited */
		return SolidMatListT::NewSub(name);
}

/* accept parameter list */
void FSSolidMatList2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SolidMatListT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	 * here we construct as many materials as are passed in */
	AutoArrayT<FSSolidMatT*> materials;
	const ArrayT<ParameterListT>& subs = list.Lists();
	for (int i = 0; i < subs.Length(); i++) {
		const ParameterListT& sub = subs[i];
		FSSolidMatT* mat = NewFSSolidMat(sub.Name());
		if (mat) {
			materials.Append(mat);
			mat->TakeParameterList(sub);

			/* set flags */
			if (mat->HasHistory()) fHasHistory = true;	
			if (mat->HasThermalStrain()) fHasThermal = true;
			if (mat->HasLocalization()) fHasLocalizers = true;
		}
	}

	/* transfer */
	Dimension(materials.Length());
	for (int i = 0; i < materials.Length(); i++)
		fArray[i] = materials[i];
}

/* construct the specified material or NULL if the request cannot be completed */
FSSolidMatT* FSSolidMatList2DT::NewFSSolidMat(const StringT& name) const
{
	FSSolidMatT* mat = NULL;

	if (name == "large_strain_Hookean_2D")
		mat = new FDHookeanMat2DT;
	else if (name == "large_strain_cubic_2D")
		mat = new FDCubic2DT;
	else if (name == "large_strain_StVenant_2D")
		mat = new FDKStV2D;
	else if (name == "Simo_isotropic_2D")
		mat = new SimoIso2D;
	else if (name == "quad_log_2D")
		mat = new QuadLog2D;
	else if (name == "quad_log_Ogden_2D")
		mat = new QuadLogOgden2DT;

#ifdef PLASTICITY_J2_MATERIAL
	else if (name == "Simo_J2_2D")
		mat = new J2Simo2D;
	else if (name == "quad_log_J2_2D")
		mat = new J2QL2DLinHardT;
#endif

#ifdef CAUCHY_BORN_MATERIAL
	else if (name == "LJ_triangular_2D")
		mat = new LJTr2D;
	else if (name == "hex_2D")
		mat = new Hex2D;
	else if (name == "FCC_EAM_2D")
		mat = new EAMFCC2D;
#endif

#ifdef MODCBSW_MATERIAL
	else if (name == "Cauchy-Born_diamond_2D")
		mat = new ModCB2DT;
#endif

#ifdef VIB_MATERIAL
	else if (name == "VIB_2D")
		mat = new VIB2D;
	else if (name == "isotropic_VIB_2D")
		mat = new IsoVIB2D;
	else if (name == "Ogden_isotropic_VIB_2D")
		mat = new OgdenIsoVIB2D;
#endif

#ifdef VISCOELASTICITY
	else if (name == "Reese-Govindjee_split")
		mat= new RGSplitT;
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
	else if (name == "ABAQUS_UMAT_BCJ")
		mat= new ABAQUS_BCJ;
	else if (name == "ABAQUS_VUMAT_BCJ")
		mat= new ABAQUS_VUMAT_BCJ;
	else if (name == "ABAQUS_UMAT_BCJ_iso-damage")
		mat= new ABAQUS_BCJ_ISO;
#endif
#endif

	/* set support */
	if (mat) mat->SetFSMatSupport(fFSMatSupport);

	return mat;

}
