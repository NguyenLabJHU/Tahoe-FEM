/* $Id: ElementListT.cpp,v 1.97 2004-07-27 20:08:50 raregue Exp $ */
/* created: paklein (04/20/1998) */
#include "ElementListT.h"
#include "ElementsConfig.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentElementsConfig.h"
#endif

#include <iostream.h>
#include "ifstreamT.h"
#include "StringT.h"
#include "ElementT.h"
#include "ElementSupportT.h"
#include "GeometryT.h"

#include "ElementBaseT.h"

#ifdef ADHESION_ELEMENT
#include "AdhesionT.h"
#endif

#ifdef COHESIVE_SURFACE_ELEMENT
#include "CSEIsoT.h"
#include "CSEAnisoT.h"
#include "CSESymAnisoT.h"
#include "MeshFreeCSEAnisoT.h"
#include "ThermalSurfaceT.h"
#ifdef COHESIVE_SURFACE_ELEMENT_DEV
#include "RigidCSEAnisoT.h"
#include "NodalRigidCSEAnisoT.h"
#endif
#endif

#ifdef CONTINUUM_ELEMENT
#include "ViscousDragT.h"
#include "SmallStrainT.h"
#include "SmallStrainAxiT.h"
#include "UpdatedLagrangianT.h"
#include "UpdatedLagrangianAxiT.h"
#include "UpLagAdaptiveT.h"
#include "TotalLagrangianT.h"
#include "TotalLagrangianAxiT.h"
#include "LocalizerT.h"
#include "SimoFiniteStrainT.h"
#include "SimoQ1P0.h"
#include "SimoQ1P0_inv.h"
#include "SimoQ1P0Axi.h"
#include "SimoQ1P0Axi_inv.h"
#include "DiffusionElementT.h"
#include "NLDiffusionElementT.h"
#include "MeshFreeSSSolidT.h"
#include "MeshFreeFSSolidT.h"
#include "MeshFreeFSSolidAxiT.h"
#include "D2MeshFreeFSSolidT.h"
#include "SS_SCNIMFT.h"
#include "FS_SCNIMFT.h"
#include "UpLagr_ExternalFieldT.h"
#ifdef SIMPLE_SOLID_DEV
#include "TotalLagrangianFlatT.h"
#endif
#endif

#ifdef BRIDGING_ELEMENT
#include "BridgingScaleT.h"
#include "MeshfreeBridgingT.h"
#endif

#ifdef CONTACT_ELEMENT
#include "PenaltyContact2DT.h"
#include "PenaltyContact3DT.h"
#include "AugLagContact2DT.h"
#include "AugLagContact3DT.h"
#include "ACME_Contact3DT.h"
#include "PenaltyContactDrag2DT.h"
#include "PenaltyContactDrag3DT.h"
#ifdef CONTINUUM_ELEMENT /* need meshfree code */
#include "MFPenaltyContact2DT.h"
#endif
#endif

#ifdef PARTICLE_ELEMENT
#include "ParticlePairT.h"
#include "EAMT.h"
#endif

#ifdef SPRING_ELEMENT
#include "SWDiamondT.h"
#include "MixedSWDiamondT.h"
#include "VirtualSWDC.h"
#include "RodT.h"
#include "UnConnectedRodT.h"
#include "VirtualRodT.h"
#endif

#ifdef CONTACT_ELEMENT_DEV
#include "MultiplierContactElement3DT.h"
#include "MultiplierContactElement2DT.h"
#include "PenaltyContactElement2DT.h"
#include "PenaltyContactElement3DT.h"
#endif

#ifdef BEM_ELEMENT_DEV
#include "BEMelement.h"
#endif

#ifdef MULTISCALE_ELEMENT_DEV
#include "StaggeredMultiScaleT.h"
#endif

#ifdef MULTISCALE_APS_DEV
#include "APS_AssemblyT.h"
#endif

#ifdef MULTISCALE_APS_V_DEV
#include "APS_V_AssemblyT.h"
#endif

#ifdef MESHFREE_GRAD_PLAST_DEV
#include "MFGP_AssemblyT.h"
#endif

#ifdef ENHANCED_STRAIN_LOC_DEV
#include "SmallStrainEnhLocT.h"
#endif

#ifdef GRAD_SMALL_STRAIN_DEV
#include "GradSmallStrainT.h"
#endif

#ifdef SOLID_ELEMENT_DEV
#ifdef MATERIAL_FORCE_ELEMENT_DEV
#include "SmallStrainQ2P1.h"
#include "UpLagMF.h"
#include "SmallStrainMF.h"
#include "SSMF.h"
#include "SSQ2P1MF.h"
#include "SmallStrainQ1P0.h"
#include "SSQ1P0MF.h"
#endif /* MATERIAL_FORCE_ELEMENT_DEV */

#ifdef SPLIT_INTEGRATION_DEV
#include "SplitIntegrationT.h"
#endif
#endif

using namespace Tahoe;

/* constructors */
ElementListT::ElementListT(FEManagerT& fe):
	ParameterInterfaceT("element_list"),
	fHasContact(false)
{
	/* initialize element support */
	fSupport.SetFEManager(&fe);
}

/* destructor */
ElementListT::~ElementListT(void)
{
	/* activate all */
	if (fAllElementGroups.Length() > 0) {
		ArrayT<bool> mask(fAllElementGroups.Length());
		mask = true;
		SetActiveElementGroupMask(mask);
	}
}

/* returns true of ALL element groups have interpolant DOF's */
bool ElementListT::InterpolantDOFs(void) const
{
	bool are_interp = true;
	for (int i = 0; i < Length(); i++)
		if (!fArray[i]->InterpolantDOFs()) are_interp = false;

	return are_interp;
}

/* change the number of active element groups */
void ElementListT::SetActiveElementGroupMask(const ArrayT<bool>& mask)
{
	/* first time */
	if (fAllElementGroups.Length() == 0) 
	{
		/* cache all pointers */
		fAllElementGroups.Dimension(Length());
		for (int i = 0; i < fAllElementGroups.Length(); i++)
		{
			ElementBaseT* element = (*this)[i];
			fAllElementGroups[i] = element;
		}
	}

	/* check */
	if (mask.Length() != fAllElementGroups.Length())
		ExceptionT::SizeMismatch("ElementListT::SetActiveElementGroupMask",
			"expecting mask length %d not %d", fAllElementGroups.Length(), mask.Length());

	/* reset active element groups */
	int num_active = 0;
	for (int i = 0; i < mask.Length(); i++)
		if (mask[i])
			num_active++;
			
	/* cast this to an ArrayT */
	ArrayT<ElementBaseT*>& element_list = *this;
	element_list.Dimension(num_active);
	num_active = 0;
	for (int i = 0; i < mask.Length(); i++)
		if (mask[i])
			element_list[num_active++] = fAllElementGroups[i];
}

/* information about subordinate parameter lists */
void ElementListT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* the element groups - an array of choices */
	sub_list.AddSub("element_groups", ParameterListT::OnePlus, true);
}

/* return the description of the given inline subordinate parameter list */
void ElementListT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "element_groups")
	{
		order = ParameterListT::Choice;

#ifdef COHESIVE_SURFACE_ELEMENT
		sub_lists.AddSub("isotropic_CSE");
		sub_lists.AddSub("anisotropic_CSE");
		sub_lists.AddSub("anisotropic_symmetry_CSE");
		sub_lists.AddSub("thermal_CSE");
#endif

#ifdef ADHESION_ELEMENT
		sub_lists.AddSub("adhesion");
#endif

#ifdef CONTACT_ELEMENT
		sub_lists.AddSub("contact_2D_penalty");
		sub_lists.AddSub("contact_3D_penalty");

		sub_lists.AddSub("contact_2D_multiplier");
		sub_lists.AddSub("contact_3D_multiplier");

		sub_lists.AddSub("contact_drag_2D_penalty");
		sub_lists.AddSub("contact_drag_3D_penalty");
#endif

#ifdef CONTACT_ELEMENT
#ifdef CONTINUUM_ELEMENT /* need meshfree code */
		sub_lists.AddSub("meshfree_contact_2D_penalty");
#endif
#endif

#ifdef PARTICLE_ELEMENT
		sub_lists.AddSub("particle_pair");
		sub_lists.AddSub("particle_EAM");
#endif

#ifdef CONTINUUM_ELEMENT
		sub_lists.AddSub("diffusion");
		sub_lists.AddSub("nonlinear_diffusion");
		sub_lists.AddSub("small_strain");
		sub_lists.AddSub("updated_lagrangian");
		sub_lists.AddSub("updated_lagrangian_Q1P0");
		sub_lists.AddSub("updated_lagrangian_Q1P0_inv");
		sub_lists.AddSub("total_lagrangian");
		sub_lists.AddSub("small_strain_meshfree");
		sub_lists.AddSub("large_strain_meshfree");
		sub_lists.AddSub("small_strain_axi");
		sub_lists.AddSub("updated_lagrangian_axi");
		sub_lists.AddSub("total_lagrangian_axi");
		sub_lists.AddSub("updated_lagrangian_Q1P0_axi");
		sub_lists.AddSub("updated_lagrangian_Q1P0_inv_axi");
		sub_lists.AddSub("large_strain_meshfree_axi");

#ifdef BRIDGING_ELEMENT
		sub_lists.AddSub("bridging");
		sub_lists.AddSub("meshfree_bridging");
#endif
#endif /* CONTINUUM_ELEMENT */

#ifdef GRAD_SMALL_STRAIN_DEV
		sub_lists.AddSub("grad_small_strain");
#endif

/*
#ifdef MULTISCALE_ELEMENT_DEV
		sub_lists.AddSub("variational_multiscale");
#endif
*/

#ifdef MULTISCALE_APS_DEV
		sub_lists.AddSub("antiplane_shear_grad_plast");
#endif

/*
#ifdef MULTISCALE_APS_V_DEV
		sub_lists.AddSub("antiplane_shear_grad_plast_V");
#endif

#ifdef MESHFREE_GRAD_PLAST_DEV
		sub_lists.AddSub("meshfree_grad_plast");
#endif
*/

#ifdef ENHANCED_STRAIN_LOC_DEV
		sub_lists.AddSub("ss_enhanced_strain_loc");
#endif

	}
	else /* inherited */
		ParameterInterfaceT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ElementListT::NewSub(const StringT& name) const
{
	/* try to construct element */
	ElementBaseT* element = NewElement(name);
	if (element)
		return element;
	else /* inherited */	
		return ParameterInterfaceT::NewSub(name);
}

/* accept parameter list */
void ElementListT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* dimension */
	const ArrayT<ParameterListT>& subs = list.Lists();
	Dimension(subs.Length());
	for (int i = 0; i < Length(); i++) {

		/* construct element */
		ElementBaseT* element = NewElement(subs[i].Name());
		if (!element)
			ExceptionT::GeneralFail("ElementListT::TakeParameterList", "could not construct \"%s\"");
		
		/* initialize */
		element->TakeParameterList(subs[i]);
		
		/* store */
		fArray[i] = element;
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* return a pointer to a new element group or NULL if the request cannot be completed */
ElementBaseT* ElementListT::NewElement(const StringT& name) const
{
	if (false) /* dummy */
		return NULL;

#ifdef COHESIVE_SURFACE_ELEMENT	
	else if (name == "isotropic_CSE")
		return new CSEIsoT(fSupport);
		
	else if (name == "anisotropic_CSE")
		return new CSEAnisoT(fSupport);

	else if (name == "anisotropic_symmetry_CSE")
		return new CSESymAnisoT(fSupport);

	else if (name == "thermal_CSE")
		return new ThermalSurfaceT(fSupport);
#endif

#ifdef ADHESION_ELEMENT
	else if (name == "adhesion")
		return new AdhesionT(fSupport);
#endif

#ifdef CONTACT_ELEMENT
	else if (name == "contact_2D_penalty")
		return new PenaltyContact2DT(fSupport);
	else if (name == "contact_3D_penalty")
		return new PenaltyContact3DT(fSupport);

	else if (name == "contact_2D_multiplier")
		return new AugLagContact2DT(fSupport);
	else if (name == "contact_3D_multiplier")
		return new AugLagContact3DT(fSupport);

	else if (name == "contact_drag_2D_penalty")
		return new PenaltyContactDrag2DT(fSupport);
	else if (name == "contact_drag_3D_penalty")
		return new PenaltyContactDrag3DT(fSupport);
#endif

#ifdef CONTACT_ELEMENT
#ifdef CONTINUUM_ELEMENT /* need meshfree code */
	else if (name == "meshfree_contact_2D_penalty")
		return new MFPenaltyContact2DT(fSupport);
#endif
#endif

#ifdef PARTICLE_ELEMENT
	else if (name == "particle_pair")
		return new ParticlePairT(fSupport);
	else if (name == "particle_EAM")
		return new EAMT(fSupport);
#endif

#ifdef CONTINUUM_ELEMENT
	else if (name == "diffusion")
		return new DiffusionElementT(fSupport);
	else if (name == "nonlinear_diffusion")
		return new NLDiffusionElementT(fSupport);
	else if (name == "small_strain")
		return new SmallStrainT(fSupport);
	else if (name == "updated_lagrangian")
		return new UpdatedLagrangianT(fSupport);
	else if (name == "updated_lagrangian_Q1P0")
		return new SimoQ1P0(fSupport);
	else if (name == "updated_lagrangian_Q1P0_inv")
		return new SimoQ1P0_inv(fSupport);
	else if (name == "total_lagrangian")
		return new TotalLagrangianT(fSupport);
	else if (name == "small_strain_meshfree")
		return new MeshFreeSSSolidT(fSupport);
	else if (name == "large_strain_meshfree")
		return new MeshFreeFSSolidT(fSupport);
	else if (name == "small_strain_axi")
		return new SmallStrainAxiT(fSupport);
	else if (name == "updated_lagrangian_axi")
		return new UpdatedLagrangianAxiT(fSupport);
	else if (name == "total_lagrangian_axi")
		return new TotalLagrangianAxiT(fSupport);
	else if (name == "updated_lagrangian_Q1P0_axi")
		return new SimoQ1P0Axi(fSupport);
	else if (name == "updated_lagrangian_Q1P0_inv_axi")
		return new SimoQ1P0Axi_inv(fSupport);
	else if (name == "large_strain_meshfree_axi")
		return new MeshFreeFSSolidAxiT(fSupport);

#ifdef BRIDGING_ELEMENT
	else if (name == "bridging")
		return new BridgingScaleT(fSupport);
	else if (name == "meshfree_bridging")
		return new MeshfreeBridgingT(fSupport);
#endif
#endif /* CONTINUUM_ELEMENT */

#ifdef GRAD_SMALL_STRAIN_DEV
	else if (name == "grad_small_strain")
		return new GradSmallStrainT(fSupport);
#endif	

/*
#ifdef MULTISCALE_ELEMENT_DEV
	else if (name == "variational_multiscale")
		return new StaggeredMultiScaleT(fSupport);
#endif
*/

#ifdef MULTISCALE_APS_DEV
	else if (name == "antiplane_shear_grad_plast")
		return new APS_AssemblyT(fSupport);
#endif

/*
#ifdef MULTISCALE_APS_V_DEV
	else if (name == "antiplane_shear_grad_plast_V")
		return new APS_V_AssemblyT(fSupport);
#endif

#ifdef MESHFREE_GRAD_PLAST_DEV
	else if (name == "meshfree_grad_plast")
		return new MFGP_AssemblyT(fSupport);
#endif
*/

/*
#ifdef ENHANCED_STRAIN_LOC_DEV
	else if (name == "ss_enhanced_strain_loc")
		return new SmallStrainEnhLocT(fSupport);
#endif
*/

	/* default */	
	else
		return NULL;
}
