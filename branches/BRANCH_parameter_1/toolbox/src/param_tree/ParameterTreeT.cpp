/* $Id: ParameterTreeT.cpp,v 1.1 2003-04-26 19:08:43 paklein Exp $ */
#include "ParameterTreeT.h"
#include "ParameterInterfaceT.h"

using namespace Tahoe;

ParameterTreeT::ParameterTreeT(void) {}

ParameterTreeT::~ParameterTreeT(void)
{
	for (int i = 0; i < fBranches.Length(); i++)
		delete fBranches[i];
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

	/* collect parameters */
	source.DefineParameters(params);
	
	/* add list to the dictionary */
	if (!fDictionary.Insert(source.Name(), &source))
		ExceptionT::GeneralFail(caller, "dictionary already contains \"%s\"",
			source.Name().Pointer());

	/* sublists */
	ArrayT<StringT> sub_lists;
	ArrayT<ParameterListT::OccurrenceT> occur;
	source.SubListNames(sub_lists, occur);
	for (int i = 0; i < sub_lists.Length(); i++) {

		/* already in the tree */
		if (fDictionary.HasKey(sub_lists[i]))
		{
			/* add as reference */
			if (!params.AddReference(sub_lists[i], occur[i]))
				ExceptionT::GeneralFail(caller, "could not add reference \"%s\" to list \"%s\"",
					sub_lists[i].Pointer(), params.Name().Pointer());
		}
		else
		{
			/* sublist */
			ParameterInterfaceT* sub = source.SubList(sub_lists[i]);
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
		ParameterInterfaceT* sub = source.SubList(sub_lists[i].Name());
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
	}
}
