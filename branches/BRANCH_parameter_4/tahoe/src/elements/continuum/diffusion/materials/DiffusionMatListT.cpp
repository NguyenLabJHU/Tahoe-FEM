/* $Id: DiffusionMatListT.cpp,v 1.9.2.2 2004-07-07 15:28:04 paklein Exp $ */
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

/* read material data from the input stream */
void DiffusionMatListT::ReadMaterialData(ifstreamT& in)
{
	const char caller[] = "DiffusionMatListT::ReadMaterialData";

ExceptionT::Stop(caller);
#if 0
	/* read material data */
	for (int i = 0; i < fLength; i++)
	{
		int matnum, matcode;
		in >> matnum; matnum--;
		in >> matcode;
		
		/* checks */
		if (matnum < 0  || matnum >= fLength) ExceptionT::BadInputValue(caller);

		/* repeated material number */
		if (fArray[matnum] != NULL)
			ExceptionT::BadInputValue(caller, "repeated material number: %d", matnum+1);
		
		/* add to the list of materials */
		switch (matcode)
		{
			case kLinear:
			{
				fArray[matnum] = new DiffusionMaterialT(in, *fDiffusionMatSupport);
				break;
			}
			case kNonLinear:
			{
				fArray[matnum] = new NLDiffusionMaterialT(in, *fDiffusionMatSupport);
				break;
			}
			default:
				ExceptionT::BadInputValue(caller, "unknown material code: %d", matcode);
		}

		/* verify construction */
		if (!fArray[matnum]) ExceptionT::OutOfMemory(caller);
	}
#endif
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
void DiffusionMatListT::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
	SubListT& sub_sub_list) const
{
	/* list of choice of materials */
	if (sub == "diffusion_material_list")
	{
		order = ParameterListT::Choice;
	
		/* diffusion materials */
		sub_sub_list.AddSub("linear_diffusion_material");
		sub_sub_list.AddSub("nonlinear_diffusion_material");
	}	
	else /* inherited */
		ParameterInterfaceT::DefineInlineSub(sub, order, sub_sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* DiffusionMatListT::NewSub(const StringT& list_name) const
{
	/* try to construct material */
	DiffusionMaterialT* material = NewDiffusionMaterial(list_name);
	if (material)
		return material;
	else /* inherited */
		return MaterialListT::NewSub(list_name);
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
DiffusionMaterialT* DiffusionMatListT::NewDiffusionMaterial(const StringT& list_name) const
{
	DiffusionMaterialT* mat = NULL;

	if (list_name == "linear_diffusion_material")
		mat = new DiffusionMaterialT;	
	else if (list_name == "nonlinear_diffusion_material")
		mat = new NLDiffusionMaterialT;

	/* set support */
	if (mat) mat->SetDiffusionMatSupport(fDiffusionMatSupport);

	return mat;
}
