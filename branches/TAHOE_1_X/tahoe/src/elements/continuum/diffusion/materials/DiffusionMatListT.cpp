/* $Id: DiffusionMatListT.cpp,v 1.9 2004-06-17 07:40:04 paklein Exp $ */
/* created: paklein (02/14/1997) */
#include "DiffusionMatListT.h"
#include "DiffusionMatSupportT.h"
#include "ifstreamT.h"

/* diffusion materials */
#include "DiffusionMaterialT.h"
#include "NLDiffusionMaterialT.h"

using namespace Tahoe;

/* constructors */
DiffusionMatListT::	DiffusionMatListT(int length, const DiffusionMatSupportT& support):
	MaterialListT(length),
	fDiffusionMatSupport(&support)
{
	SetName("diffusion_materials");
}

DiffusionMatListT::	DiffusionMatListT(void):
	fDiffusionMatSupport(NULL)
{
	SetName("diffusion_materials");
}

/* read material data from the input stream */
void DiffusionMatListT::ReadMaterialData(ifstreamT& in)
{
	const char caller[] = "DiffusionMatListT::ReadMaterialData";

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
}

/* information about subordinate parameter lists */
void DiffusionMatListT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	MaterialListT::DefineSubs(sub_list);

	/* an array of choices */
	sub_list.AddSub("diffusion_material_list", ParameterListT::OnePlus, true);

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
		sub_sub_list.AddSub("linear_diffusion");
		sub_sub_list.AddSub("nonlinear_diffusion");
	}	
	else /* inherited */
		ParameterInterfaceT::DefineInlineSub(sub, order, sub_sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* DiffusionMatListT::NewSub(const StringT& list_name) const
{
	if (list_name == "linear_diffusion")
		return new DiffusionMaterialT;	
	else if (list_name == "nonlinear_diffusion")
		return new NLDiffusionMaterialT;
	else /* inherited */
		return MaterialListT::NewSub(list_name);
}

