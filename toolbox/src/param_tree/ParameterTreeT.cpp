/* $Id: ParameterTreeT.cpp,v 1.1.2.3 2003-04-28 17:05:30 paklein Exp $ */
#include "ParameterTreeT.h"
#include "ParameterInterfaceT.h"

using namespace Tahoe;

ParameterTreeT::ParameterTreeT(void) {}

ParameterTreeT::~ParameterTreeT(void)
{
	/* free parameter list data */
	for (int i = 0; i < fBranches.Length(); i++)
		delete fBranches[i];

	/* free instantiations */
	for (int i = 0; i < fDeleteMe.Length(); i++)
		delete fDeleteMe[i];
}

/* add a branch to the tree with the given root */
void ParameterTreeT::BuildDescription(ParameterInterfaceT& root)
{
	/* branch exists */
	if (fDictionary.HasKey(root.Name())) return;

	/* collect parameters into a new list */
	ParameterListT* root_list = new ParameterListT(root.Name());
	BuildBranch(root, *root_list);
	
	/* store branch */
	fBranches.Append(root_list);
}

/* add only parts of the branch contained in the guide parameter list */
void ParameterTreeT::BuildDescription(ParameterInterfaceT& root, const ParameterListT& guide)
{
	/* branch exists */
	if (fDictionary.HasKey(root.Name())) return;

	/* collect parameters into a new list */
	ParameterListT* root_list = new ParameterListT(root.Name());
	BuildBranch(root, guide, *root_list);
	
	/* store branch */
	fBranches.Append(root_list);
}

/* return the parameter list for the given branch or NULL if it does not exist */
const ParameterListT* ParameterTreeT::Branch(const StringT& name) const
{
	for (int i = 0; i < fBranches.Length(); i++)
		if (fBranches[i]->Name() == name)
			return fBranches[i];
		
	return NULL;
}

/*************************************************************************
 * Private
 *************************************************************************/

/* build the branch */
void ParameterTreeT::BuildBranch(ParameterInterfaceT& source, ParameterListT& params)
{
	const char caller[] = "ParameterTreeT::BuildBranch";
	
	/* contents of subordinate lists */
	ArrayT<StringT> sub_lists;
	ArrayT<ParameterListT::OccurrenceT> occur;
	ArrayT<bool> is_inline;

	/* defining self */
	if (source.Name() == params.Name()) {
		
		/* collect parameters */
		source.DefineParameters(params);
			
		/* add list to the dictionary */
		if (!fDictionary.Insert(source.Name(), &source))
			ExceptionT::GeneralFail(caller, "dictionary already contains \"%s\"",
				source.Name().Pointer());

		/* collect information on subordinate lists */
		source.SubNames(sub_lists, occur, is_inline);
	}
	else /* defining inlined list */
	{
		/* list must be inlined */
		if (!params.Inline())
			ExceptionT::GeneralFail(caller, "source \"%s\" cannot build list for non-inlined \"%s\"",
				source.Name().Pointer(), params.Name().Pointer());
		
		/* collect information on inlined list */
		ParameterListT::ListOrderT order;
		source.DefineInlineSub(params.Name(), order, sub_lists, occur, is_inline);
		params.SetListOrder(order);
	}

	/* define sublists */
	for (int i = 0; i < sub_lists.Length(); i++) {

		/* list is inline */
		if (is_inline[i])
		{
			/* source builds inline list */
			ParameterListT inline_params(sub_lists[i]);
			inline_params.SetInline(true);
			BuildBranch(source, inline_params);
			
			/* add list */
			if (!params.AddList(inline_params, occur[i]))
				ExceptionT::GeneralFail(caller, "could not add inline \"%s\" to list \"%s\"",
					inline_params.Name().Pointer(), params.Name().Pointer());
		}
		/* already in the tree */
		else if (fDictionary.HasKey(sub_lists[i]))
		{
			/* add as reference */
			if (!params.AddReference(sub_lists[i], occur[i]))
				ExceptionT::GeneralFail(caller, "could not add reference \"%s\" to list \"%s\"",
					sub_lists[i].Pointer(), params.Name().Pointer());
		}
		else
		{
			/* sublist */
			ParameterInterfaceT* sub = source.NewSub(sub_lists[i]);
			if (!sub)
				ExceptionT::GeneralFail(caller, "source \"%s\" did not return sublist \"%s\"",
					params.Name().Pointer(), sub_lists[i].Pointer());
					
			/* build the sublist */
			ParameterListT sub_params(sub->Name());
			BuildBranch(*sub, sub_params);
			
			/* add list */
			if (!params.AddList(sub_params, occur[i]))
				ExceptionT::GeneralFail(caller, "could not add sublist \"%s\" to list \"%s\"",
					sub_params.Name().Pointer(), params.Name().Pointer());
	
			/* add to list of items to delete */
			fDeleteMe.Append(sub);
		}
	}	
}

/* build the branch */
void ParameterTreeT::BuildBranch(ParameterInterfaceT& source, const ParameterListT& guide, 
	ParameterListT& params)
{
	const char caller[] = "ParameterTreeT::BuildBranch";

	/* collect parameters */
	source.DefineParameters(params);
	
	/* add list to the dictionary */
	if (!fDictionary.Insert(source.Name(), &source))
		ExceptionT::GeneralFail(caller, "dictionary already contains \"%s\"",
			source.Name().Pointer());

	/* construct sublists */
	const ArrayT<ParameterListT>& sub_lists = guide.Lists();
	const ArrayT<ParameterListT::OccurrenceT>& occur = guide.ListOccurrences();
	for (int i = 0; i < sub_lists.Length(); i++) {

		/* fetch sublist */
		ParameterInterfaceT* sub = source.NewSub(sub_lists[i].Name());
		if (!sub)
			ExceptionT::GeneralFail(caller, "source \"%s\" did not return sublist \"%s\"",
				params.Name().Pointer(), sub_lists[i].Name().Pointer());
					
		/* build the sublist */
		ParameterListT sub_params;
		BuildBranch(*sub, sub_lists[i], sub_params);
			
		/* add list */
		if (!params.AddList(sub_params, occur[i]))
			ExceptionT::GeneralFail(caller, "could not add sublist \"%s\" to list \"%s\"",
				sub_params.Name().Pointer(), params.Name().Pointer());
				
		/* add to list of items to delete */
		fDeleteMe.Append(sub);
	}
}
