/* $Id: FSSolidMatList3DT.cpp,v 1.2 2004-07-15 08:28:28 paklein Exp $ */
/* created: paklein (02/14/1997) */
#include "FSSolidMatList3DT.h"

#include "SolidMaterialsConfig.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#endif

#include "FDKStV.h"
#include "FDCubicT.h"
#include "QuadLog3D.h"
#include "SimoIso3D.h"
#include "QuadLogOgden3DT.h"

#ifdef CAUCHY_BORN_MATERIAL
#include "EAMFCC3DMatT.h"
#include "FCC3D.h"
#endif

#ifdef MODCBSW_MATERIAL
#include "ModCB3DT.h"
#endif

#ifdef VIB_MATERIAL
#include "VIB3D.h"
#include "IsoVIB3D.h"
#include "J2IsoVIB3DLinHardT.h"
#include "OgdenIsoVIB3D.h"
#endif

#ifdef VISCOELASTICITY
#include "RGSplitT.h"
#endif

#ifdef VISCOELASTIC_MATERIALS_DEV
#include "FDSV_KStV3D.h"
#include "RGSplit3D.h"
#endif

#ifdef ELASTIC_OGDEN_MATERIAL_DEV
#include "OgdenMaterialT.h"
#endif

#ifdef PLASTICITY_CRYSTAL_MATERIAL
#include "LocalCrystalPlast.h"
#include "LocalCrystalPlast_C.h"
#include "GradCrystalPlast.h"
#include "LocalCrystalPlastFp.h"
#include "LocalCrystalPlastFp_C.h"
#include "GradCrystalPlastFp.h"
#endif

#ifdef PLASTICITY_MACRO_MATERIAL
#include "HyperEVP3D.h"
#include "BCJHypo3D.h"
#include "BCJHypoIsoDamageKE3D.h"
#include "BCJHypoIsoDamageYC3D.h"
#endif

#ifdef PLASTICITY_J2_MATERIAL
#include "J2Simo3D.h"
#include "J2QLLinHardT.h"
#endif

#ifdef ELASTICITY_CRYSTAL_MATERIAL
#include "FDCrystalElast.h"
#endif

#ifdef THERMO_VISCO_PLASTIC_MATERIAL
#include "tevp3D.h"
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
#include "ABAQUS_BCJ.h"
#include "ABAQUS_BCJ_ISO.h"
#include "ABAQUS_VUMAT_BCJ.h"
#endif
#endif

#ifdef SIERRA_MATERIAL
#include "SIERRA_HypoElasticT.h"
#ifdef __FOSSUM__
#include "SIERRA_Isotropic_Geomaterial.h"
#endif /* __FOSSUM__ */
#endif /* SIERRA_MATERIAL */

using namespace Tahoe;

/* constructors */
FSSolidMatList3DT::FSSolidMatList3DT(int length, const FSMatSupportT& support):
	SolidMatListT(length, support),
	fFSMatSupport(&support)
{
#pragma message("check spatial dimension of material support")

	SetName("large_strain_material_3D");
}

FSSolidMatList3DT::FSSolidMatList3DT(void):
	fFSMatSupport(NULL)
{
	SetName("large_strain_material_3D");
}

/* information about subordinate parameter lists */
void FSSolidMatList3DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("fs_material_list_3D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void FSSolidMatList3DT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const
{
	if (name == "fs_material_list_3D")
	{
		order = ParameterListT::Choice;
	
		sub_lists.AddSub("large_strain_cubic");
		sub_lists.AddSub("large_strain_StVenant");
		sub_lists.AddSub("Simo_isotropic");
		sub_lists.AddSub("quad_log");
		sub_lists.AddSub("quad_log_Ogden");

#ifdef PLASTICITY_J2_MATERIAL
		sub_lists.AddSub("Simo_J2");
		sub_lists.AddSub("quad_log_J2");
#endif

#ifdef CAUCHY_BORN_MATERIAL
		sub_lists.AddSub("FCC_3D");
		sub_lists.AddSub("FCC_EAM");		
#endif

#ifdef MODCBSW_MATERIAL
		sub_lists.AddSub("Cauchy-Born_diamond");
#endif

#ifdef VIB_MATERIAL
		sub_lists.AddSub("VIB");
		sub_lists.AddSub("isotropic_VIB");
		sub_lists.AddSub("Ogden_isotropic_VIB");
#endif

#ifdef VISCOELASTICITY
	sub_lists.AddSub("Reese-Govindjee_split");
#endif
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSSolidMatList3DT::NewSub(const StringT& name) const
{
	/* try to construct material */
	FSSolidMatT* fs_solid_mat = NewFSSolidMat(name);
	if (fs_solid_mat)
		return fs_solid_mat;
	else /* inherited */
		return SolidMatListT::NewSub(name);
}

/* accept parameter list */
void FSSolidMatList3DT::TakeParameterList(const ParameterListT& list)
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
FSSolidMatT* FSSolidMatList3DT::NewFSSolidMat(const StringT& name) const
{
	FSSolidMatT* mat = NULL;

	if (name == "large_strain_cubic")
		mat = new FDCubicT;
	else if (name == "large_strain_StVenant")
		mat = new FDKStV;
	else if (name == "Simo_isotropic")
		mat = new SimoIso3D;
	else if (name == "quad_log")
		mat = new QuadLog3D;
	else if (name == "quad_log_Ogden")
		mat = new QuadLogOgden3DT;

#ifdef PLASTICITY_J2_MATERIAL
	else if (name == "Simo_J2")
		mat = new J2Simo3D;
	else if (name == "quad_log_J2")
		mat = new J2QLLinHardT;
#endif

#ifdef CAUCHY_BORN_MATERIAL
	else if (name == "FCC_3D")
		mat = new FCC3D;
	else if (name == "FCC_EAM")
		mat = new EAMFCC3DMatT;
#endif

#ifdef MODCBSW_MATERIAL
	else if (name == "Cauchy-Born_diamond")
		mat= new ModCB3DT;
#endif

#ifdef VIB_MATERIAL
	else if (name == "VIB")
		mat= new VIB3D;
	else if (name == "isotropic_VIB")
		mat= new IsoVIB3D;
	else if (name == "Ogden_isotropic_VIB")
		mat= new OgdenIsoVIB3D;
#endif

#ifdef VISCOELASTICITY
	else if (name == "Reese-Govindjee_split")
		mat= new RGSplitT;
#endif

	/* set support */
	if (mat) mat->SetFSMatSupport(fFSMatSupport);

	return mat;

}
