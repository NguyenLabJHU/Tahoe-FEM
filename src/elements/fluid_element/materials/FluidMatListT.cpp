/* $Header: /home/regueiro/tahoe_cloudforge_repo_snapshots/development/src/elements/fluid_element/materials/FluidMatListT.cpp,v 1.4 2006-08-18 01:23:44 a-kopacz Exp $ */
/* created: tdnguye (07/12/2006) */
#include "FluidMatListT.h"
#include "FluidMatSupportT.h"

/* fluid materials */
#include "FluidMaterialT.h"

using namespace Tahoe;

/* constructors */
FluidMatListT::	FluidMatListT(int length, const FluidMatSupportT& support):
	MaterialListT(length),
	fFluidMatSupport(&support)
{
	SetName("fluid_material");
}

FluidMatListT::	FluidMatListT(void):
	fFluidMatSupport(NULL)
{
	SetName("fluid_material");
}

/* information about subordinate parameter lists */
void FluidMatListT::DefineSubs(SubListT& sub_list) const
{
	//WriteCallLocation("DefineSubs"); //DEBUG

	/* inherited */
	MaterialListT::DefineSubs(sub_list);

	/* an array of choices */
	sub_list.AddSub("fluid_material_list", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void FluidMatListT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	//WriteCallLocation("DefineInlineSub"); //DEBUG

	/* list of choice of materials */
	if (name == "fluid_material_list")
	{
		order = ParameterListT::Choice;
		/* fluid materials */
		sub_lists.AddSub("linear_fluid_material");
	}	
	else /* inherited */
		ParameterInterfaceT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FluidMatListT::NewSub(const StringT& name) const
{
	//WriteCallLocation("NewSub"); //DEBUG

	/* try to construct material */
	FluidMaterialT* material = NewFluidMaterial(name);
	if (material)
		return material;
	else /* inherited */
		return MaterialListT::NewSub(name);
}

/* accept parameter list */
void FluidMatListT::TakeParameterList(const ParameterListT& list)
{
	//WriteCallLocation("TakeParameterList"); //DEBUG

	/* inherited */
	MaterialListT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	* here we construct as many materials as are passed in */
	AutoArrayT<FluidMaterialT*> materials;
	const ArrayT<ParameterListT>& subs = list.Lists();
	for (int i = 0; i < subs.Length(); i++) 
	{
		const ParameterListT& sub = subs[i];
		FluidMaterialT* mat = NewFluidMaterial(sub.Name());
		if (mat) 
		{
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
FluidMaterialT* FluidMatListT::NewFluidMaterial(const StringT& name) const
{
	//WriteCallLocation("NewFluidMaterial"); //DEBUG

	FluidMaterialT* mat = NULL;

	if (name == "linear_fluid_material")
		mat = new FluidMaterialT;

	/* set support */
	if (mat) mat->SetFluidMatSupport(fFluidMatSupport);

	return mat;
}

/** FOR DEBUGGING PURPOSES ONLY */
void FluidMatListT::WriteCallLocation( char* loc ) const
{
	cout << "\n Inside of FluidMatListT::" << loc << endl;
}
