/* $Id: DiffusionMatListT.cpp,v 1.10 2004-07-15 08:26:22 paklein Exp $ */
/* created: paklein (02/14/1997) */
#include "DiffusionMatListT.h"
#include "DiffusionMatSupportT.h"

/* diffusion materials */
#include "DiffusionMaterialT.h"
#include "NLDiffusionMaterialT.h"

using namespace Tahoe;

/* constructors */
DiffusionMatListT::	DiffusionMatListT(int length, const DiffusionMatSupportT& support):
	MaterialListT(length),
	fDiffusionMatSupport(&support)
{
	SetName("diffusion_material");
}

DiffusionMatListT::	DiffusionMatListT(void):
	fDiffusionMatSupport(NULL)
{
	SetName("diffusion_material");
}

/* information about subordinate parameter lists */
void DiffusionMatListT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	MaterialListT::DefineSubs(sub_list);

	/* an array of choices */
	sub_list.AddSub("diffusion_material_list", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void DiffusionMatListT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	/* list of choice of materials */
	if (name == "diffusion_material_list")
	{
		order = ParameterListT::Choice;
	
		/* diffusion materials */
		sub_lists.AddSub("linear_diffusion_material");
		sub_lists.AddSub("nonlinear_diffusion_material");
	}	
	else /* inherited */
		ParameterInterfaceT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* DiffusionMatListT::NewSub(const StringT& name) const
{
	/* try to construct material */
	DiffusionMaterialT* material = NewDiffusionMaterial(name);
	if (material)
		return material;
	else /* inherited */
		return MaterialListT::NewSub(name);
}

/* accept parameter list */
void DiffusionMatListT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	MaterialListT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	 * here we construct as many materials as are passed in */
	AutoArrayT<DiffusionMaterialT*> materials;
	const ArrayT<ParameterListT>& subs = list.Lists();
	for (int i = 0; i < subs.Length(); i++) {
		const ParameterListT& sub = subs[i];
		DiffusionMaterialT* mat = NewDiffusionMaterial(sub.Name());
		if (mat) {
			materials.Append(mat);
			mat->TakeParameterList(sub);
		}
	}

	/* transfer */
	Dimension(materials.Length());
	for (int i = 0; i < materials.Length(); i++)
		fArray[i] = materials[i];
	
}

/* construct the specified material or NULL if the request cannot be completed */
DiffusionMaterialT* DiffusionMatListT::NewDiffusionMaterial(const StringT& name) const
{
	DiffusionMaterialT* mat = NULL;

	if (name == "linear_diffusion_material")
		mat = new DiffusionMaterialT;	
	else if (name == "nonlinear_diffusion_material")
		mat = new NLDiffusionMaterialT;

	/* set support */
	if (mat) mat->SetDiffusionMatSupport(fDiffusionMatSupport);

	return mat;
}
