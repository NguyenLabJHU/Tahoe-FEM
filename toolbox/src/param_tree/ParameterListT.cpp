/* $Id: ParameterListT.cpp,v 1.21 2004-09-14 17:58:17 paklein Exp $ */
#include "ParameterListT.h"
#include "ParameterInterfaceT.h"
#include <string.h>

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<ParameterListT>::fByteCopy = false;
DEFINE_TEMPLATE_STATIC const bool ArrayT<ParameterListT*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<ParameterListT::OccurrenceT>::fByteCopy = false;
}

/* constructor */
ParameterListT::ParameterListT(const char* name):
	fName(name),
	fListOrder(Sequence),
	fInline(false),
	fDuplicateListNames(true)
{

}

/* default constructor */
ParameterListT::ParameterListT(void):
	fListOrder(Sequence),
	fInline(false),
	fDuplicateListNames(true)
{

}

/* set/change the list type */
void ParameterListT::SetInline(bool is_inline)
{
	if (is_inline && fParameters.Length() > 0)
		ExceptionT::GeneralFail("ParameterListT::SetInline", 
			"lists with parameters cannot inlined");
	fInline = is_inline;
}

/* set/change the list order */
void ParameterListT::SetListOrder(ListOrderT list_order)
{
	if (list_order == Choice) {
	
		/* all list occurrences must be Once */
		for (int i = 0; i < fParameterListsOccur.Length(); i++)
			if (fParameterListsOccur[i] != Once)
			ExceptionT::GeneralFail("ParameterListT::SetListOrder", 
				"for list order \"Choice\" all lists must occur \"Once\"");
	}

	fListOrder = list_order;
}

/* number of nested parameter lists with the given name */
int ParameterListT::NumLists(const char* name) const
{
	int count = 0;
	for (int i = 0; i < fParameterLists.Length(); i++)
		if (fParameterLists[i].Name() == name)
			count ++;
	return count;
}

/* add parameter */
bool ParameterListT::AddParameter(const ParameterT& param, OccurrenceT occur)
{
	if (fInline)
		ExceptionT::GeneralFail("ParameterListT::AddParameter", 
			"inlined list \"%s\" cannot have parameters",
			Name().Pointer());

	/* "description" is reserved */
	if (param.Name() == "description") {
		cout << "\n ParameterListT::AddParameter: parameter name \"description\" is reserved" << endl;
		return false;
	}

	/* scan name */
	for (int i = 0; i < fParameters.Length(); i++)
		if (fParameters[i].Name() == param.Name())
			return false;
	
	/* add if no matches */
	fParameters.Append(param);
	fParametersOccur.Append(occur);
	return true;
}

/* remove parameter */
bool ParameterListT::RemoveParameter(const char* name)
{
	/* scan name */
	int index = -1;
	for (int i = 0; index == -1 && i < fParameters.Length(); i++)
		if (fParameters[i].Name() == name)
			index = i;
	
	/* remove if found */
	if (index != -1) {
		fParameters.DeleteAt(index);
		fParametersOccur.DeleteAt(index);
		return true;
	}
	else
		return false;
}

/* add a parameter list */
bool ParameterListT::AddList(const ParameterListT& param_list, OccurrenceT occur)
{	
	/* check occurrence */
	if (fListOrder == Choice && occur != Once)
		ExceptionT::GeneralFail("ParameterListT::AddList", 
			"for list order \"Choice\" all lists must occur \"Once\"");

	/* scan name */
	if (!fDuplicateListNames)
		for (int i = 0; i < fParameterLists.Length(); i++)
			if (fParameterLists[i].Name() == param_list.Name())
				return false;

	/* add to list */
	fParameterLists.Append(param_list);
	fParameterListsOccur.Append(occur);
	return true;
}

/* remove a parameter list */
bool ParameterListT::RemoveList(const char* name)
{	
	/* scan name */
	int index = -1;
	for (int i = 0; index == -1 && i < fParameterLists.Length(); i++)
		if (fParameterLists[i].Name() == name)
			index = i;

	/* remove from to list */
	if (index != -1) {
		fParameterLists.DeleteAt(index);
		fParameterListsOccur.DeleteAt(index);
		return true;
	}
	else
		return false;
}

/* search for list by name */
const ParameterListT* ParameterListT::FindList(const char* search_name, int instance) const
{
	/* search list */
	int count = 0;
	for (int i = 0; i < fParameterLists.Length(); i++) {
	
		/* candidate list */
		const ParameterListT& list = fParameterLists[i];
		const char* list_name = list.Name();

		/* look for match */
		if (strstr(list_name, search_name) != NULL)
			if (count++ == instance)
				return &list;
	}
	
	/* no match */
	return NULL;
}

/* return the list associated a choice */
const ParameterListT* ParameterListT::ListChoice(const ParameterInterfaceT& source, 
	const char* choice_name, int instance) const
{
	const char caller[] = "ParameterListT::ListChoice";

	/* check choice with an interface - sub's with interfaces still may or may not
	 * have been declared "inline". We could check this with a call to DefineSubs.
	 * Instead, we look for a named (nonlined) choice first and fall through to the
	 * search for an inlined choice if not found */
	ParameterInterfaceT* choice = source.NewSub(choice_name);
	if (choice) {
	
		/* verify the list is a choice */
		if (choice->ListOrder() != ParameterListT::Choice)
			ExceptionT::GeneralFail(caller, "\"%s\" in \"%s\" is not a choice",
				choice_name, source.Name().Pointer());
		delete choice;
	
		/* search */
		int count = 0;
		for (int i = 0; i < fParameterLists.Length(); i++)
		{
			const StringT& name = fParameterLists[i].Name();
		
			/* choice is wrapped with choice name */
			if (name == choice_name)
				if (count++ == instance) {

					/* lists within */
					const ArrayT<ParameterListT>& lists = fParameterLists[i].Lists();
					if (lists.Length() != 1)
						ExceptionT::GeneralFail(caller, "expecting only one choice in \"%s\"",
							name.Pointer());
					
					/* return */
					return lists.Pointer();
				}
		}
	}
			
	/* check choice as inline sub */
	ParameterListT::ListOrderT order;
	SubListT sub_sub_list;
	source.DefineInlineSub(choice_name, order, sub_sub_list);
	if (sub_sub_list.Length() > 0) {
	
		/* must be a choice */
		if (order != ParameterListT::Choice)
			ExceptionT::GeneralFail(caller, "\"%s\" in \"%s\" is not a choice",
				choice_name, source.Name().Pointer());
	
		/* search */
		int count = 0;
		for (int i = 0; i < fParameterLists.Length(); i++)
		{
			const StringT& name = fParameterLists[i].Name();
		
			/* run through choices */
			for (int j = 0; j < sub_sub_list.Length(); j++)
				if (name == sub_sub_list[j].Name())
					if (count++ == instance)
						return fParameterLists.Pointer(i);
		}
	}

	/* failed */
	return NULL;
}

/* return the non-const pointer to the given parameter or NULL if the list is not found */
const ParameterT* ParameterListT::Parameter(const char* name) const
{
	/* search list */
	for (int i = 0; i < fParameters.Length(); i++)
		if (fParameters[i].Name() == name)
			return fParameters.Pointer(i);

	/* fail */
	return NULL;
}

void ParameterListT::GetParameter(const char* name, int& a) const
{
	const char caller[] = "ParameterListT::GetParameter";
	const ParameterT* parameter = Parameter(name);
	if (!parameter)
		ExceptionT::GeneralFail(caller, "parameter \"%s\" not found in \"%s\"", 
			name, Name().Pointer());
	a = *parameter;
}

void ParameterListT::GetParameter(const char* name, double& a) const
{
	const char caller[] = "ParameterListT::GetParameter";
	const ParameterT* parameter = Parameter(name);
	if (!parameter)
		ExceptionT::GeneralFail(caller, "parameter \"%s\" not found in \"%s\"", 
			name, Name().Pointer());
	a = *parameter;
}

void ParameterListT::GetParameter(const char* name, StringT& a) const
{
	const char caller[] = "ParameterListT::GetParameter";
	const ParameterT* parameter = Parameter(name);
	if (!parameter)
		ExceptionT::GeneralFail(caller, "parameter \"%s\" not found in \"%s\"", 
			name, Name().Pointer());
	a = *parameter;
}

void ParameterListT::GetParameter(const char* name, bool& a) const
{
	const char caller[] = "ParameterListT::GetParameter";
	const ParameterT* parameter = Parameter(name);
	if (!parameter)
		ExceptionT::GeneralFail(caller, "parameter \"%s\" not found in \"%s\"", 
			name, Name().Pointer());
	a = *parameter;
}

/* return the given parameter */
const ParameterT& ParameterListT::GetParameter(const char* name) const
{
	const ParameterT* parameter = Parameter(name);
	if (!parameter)
		ExceptionT::GeneralFail("ParameterListT::GetParameter", 
			"parameter \"%s\" not found in \"%s\"", 
			name, Name().Pointer());
	return *parameter;
}

ParameterT& ParameterListT::GetParameter(const char* name)
{
	ParameterT* parameter = Parameter(name);
	if (!parameter)
		ExceptionT::GeneralFail("ParameterListT::GetParameter", 
			"parameter \"%s\" not found in \"%s\"", 
			name, Name().Pointer());
	return *parameter;
}

/* search for parameter by name */
const ParameterT* ParameterListT::FindParameter(const char* search_name, int instance) const
{
	/* search list */
	int count = 0;
	for (int i = 0; i < fParameters.Length(); i++) {
	
		/* candidate list */
		const ParameterT& param = fParameters[i];
		const char* name = param.Name();

		/* look for match */
		if (strstr(name, search_name) != NULL)
			if (count++ == instance)
				return &param;
	}
	
	/* no match */
	return NULL;
}

/**********************************************************************
 * Private
 **********************************************************************/

void ParameterListT::Clear(void)
{
	fName.Clear();
	fDescription.Clear();
	fParameters.Dimension(0);
	fParametersOccur.Dimension(0);
	fParameterLists.Dimension(0);
	fParameterListsOccur.Dimension(0);
}
