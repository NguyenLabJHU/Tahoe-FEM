/* $Id: ParameterListT.cpp,v 1.6.2.3 2003-05-03 09:06:52 paklein Exp $ */
#include "ParameterListT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
const bool ArrayT<ParameterListT>::fByteCopy = false;
const bool ArrayT<ParameterListT*>::fByteCopy = true;
const bool ArrayT<ParameterListT::OccurrenceT>::fByteCopy = false;
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
			"inlined lists cannot have parameters");

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

/* add a reference */
bool ParameterListT::AddReference(const char* ref, OccurrenceT occur)
{
	/* scan name */
	for (int i = 0; i < fReferences.Length(); i++)
		if (fReferences[i] == ref)
			return false;
	
	/* add if no matches */
	fReferences.Append(ref);
	fReferencesOccur.Append(occur);
	return true;
}

/* return the pointer to the given list or NULL if the list is not found */
const ParameterListT* ParameterListT::List(const char* name, int instance) const
{
	/* search list */
	int count = 0;
	for (int i = 0; i < fParameterLists.Length(); i++)
		if (fParameterLists[i].Name() == name)
			if (++count == instance) 
				return fParameterLists.Pointer(i);

	/* fail */
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

/* create a validated parameter list */
void ParameterListT::Validate(const ParameterListT& source, const ParameterListT& description)
{
	const char caller[] = "ParameterListT::Validate";
	try { /* catch all exceptions */

	/* empty current values */
	Clear();

	/* check name */
	if (source.Name() != description.Name())
		ExceptionT::BadInputValue(caller, "source list \"%s\" should be \"%s\"",
			source.Name().Pointer(), description.Name().Pointer());

	/* set name and description */
	SetName(source.Name());
	SetDescription(source.Description());

	/* parameter lists */
	const ArrayT<ParameterT>& source_parameters = source.Parameters();
	const ArrayT<ParameterT>& parameters = description.Parameters();
	const ArrayT<ParameterListT::OccurrenceT>& occurrences = description.ParameterOccurrences();

	/* validate parameters */
	for (int i = 0; i < parameters.Length(); i++)
	{
		const ParameterT& parameter = parameters[i];
		OccurrenceT occurrence = occurrences[i];

		/* search entire list */
		int count = 0;
		for (int j = 0; j < source_parameters.Length(); j++) {
				
			/* name match */
			if (source_parameters[j].Name() == parameter.Name()) {
			
				/* check occurrence */
				if (count > 0 && (occurrence == Once || occurrence == ZeroOrOnce))
					ExceptionT::BadInputValue(caller, "parameter \"%s\" cannot appear in \"%s\" more than once",
						parameter.Name().Pointer(), source.Name().Pointer());
									
				/* get value */
				ParameterT new_parameter(parameter.Type(), parameter.Name());
				new_parameter.SetDescription(parameter.Description());
				if (source_parameters[j].Type() == parameter.Type())
					new_parameter = source_parameters[j];
				else if (source_parameters[j].Type() == ParameterT::String)
				{
					const StringT& value_str = source_parameters[j];
					new_parameter.FromString(value_str);	
				}
				else
					ExceptionT::BadInputValue(caller, "source for \"%s\" must have type %d or %d: %d", 
						parameter.Name().Pointer(), ParameterT::String, parameter.Type());

				/* increment count */
				count++;
				
				/* check limits */
				if (parameter.InBounds(new_parameter))
					/* add it */
					AddParameter(new_parameter);
				else {
				
					/* again, verbose */
					parameter.InBounds(new_parameter, true);
					ExceptionT::BadInputValue(caller, "improper value for parameter \"%s\" in \"%s\"",
						parameter.Name().Pointer(), source.Name().Pointer());
				}
			}
		}
				
		/* look for default value */
		if (count == 0 && (occurrence == Once || occurrence == OnePlus)) {
			const ValueT* default_value = parameter.Default();
			if (default_value) {
				ParameterT new_parameter(parameter.Type(), parameter.Name());
				new_parameter.SetDescription(parameter.Description());
				new_parameter = *default_value;

				/* check limits */
				if (parameter.InBounds(new_parameter))
					/* add it */
					AddParameter(new_parameter);
				else {
					
					/* again, verbose */
					parameter.InBounds(new_parameter, true);
					ExceptionT::BadInputValue(caller, "improper value for parameter \"%s\" in \"%s\"",
						parameter.Name().Pointer(), source.Name().Pointer());
				}
			}
			else
				ExceptionT::BadInputValue(caller, 
					"required value \"%s\" is missing from \"%s\" and has no default value",
					parameter.Name().Pointer(), source.Name().Pointer());
		}
	}

	/* process sub-lists */
	const ArrayT<ParameterListT>& lists = description.Lists();
	const ArrayT<ParameterListT::OccurrenceT>& lists_occur = description.ListOccurrences();
	const ArrayT<ParameterListT>& source_lists = source.Lists();
	switch (description.ListOrder())
	{
		/* ordered sequence */
		case ParameterListT::Sequence:
		{
			int source_dex = 0;
			int descript_dex = 0;
			for (descript_dex = 0; descript_dex < lists.Length() && source_dex < source_lists.Length(); descript_dex++)
			{
				const ParameterListT& list = lists[descript_dex];
				ParameterListT::OccurrenceT occur = lists_occur[descript_dex];
				switch (occur)
				{
					case ParameterListT::Once:
					case ParameterListT::ZeroOrOnce:
					{
						/* look for one */
						int count = 0;
						if (source_lists[source_dex].Name() == list.Name()) 
						{
							count++;
							
							/* create validated sub-list */
							ParameterListT new_list;
							new_list.Validate(source_lists[source_dex], list);

							/* add it */
							AddList(new_list);
						
							/* next source item */
							source_dex++;
						}

						/* check */
						if (occur == ParameterListT::Once && count != 1)
							ExceptionT::BadInputValue(caller, "item %d in \"%s\" must be \"%s\"",
								NumLists() + 1, source.Name().Pointer(), list.Name().Pointer());
						
						break;
					}
					case ParameterListT::OnePlus:
					case ParameterListT::Any:
					{
						/* accept many */
						int count = 0;
						while (source_dex < source_lists.Length() && source_lists[source_dex].Name() == list.Name())
						{
							count++;

							/* create validated sub-list */
							ParameterListT new_list;
							new_list.Validate(source_lists[source_dex], list);

							/* add it */
							AddList(new_list);
						
							/* next source item */
							source_dex++;
						}
						
						/* must have at least one */
						if (lists_occur[descript_dex] == ParameterListT::OnePlus && count < 1)
							ExceptionT::BadInputValue(caller, "list \"%s\" must contain at least one \"%s\" beginning at position %d",
								source.Name().Pointer(), list.Name().Pointer(), NumLists() + 1);		
						break;
					}
					default:
						ExceptionT::GeneralFail(caller, "unrecognized occurrence %d", occur);
				}
			}
			
			/* check any remaining */
			for (; descript_dex < lists.Length(); descript_dex++)
				if (lists_occur[descript_dex] != ParameterListT::ZeroOrOnce &&
				    lists_occur[descript_dex] != ParameterListT::Any)
					ExceptionT::GeneralFail(caller, "\"%s\" required by \"%s\"",
						lists[descript_dex].Name().Pointer(), source.Name().Pointer());
			
			break;
		}
#if 0
		/* unordered list */
		case ParameterListT::Choice:
		{
			//not implemented
			break;
		}
#endif
		default:
			ExceptionT::GeneralFail(caller, "unrecognized list order %d", description.ListOrder());
	}

	} /* catch all exceptions */ 

	/* rethrow */
	catch (ExceptionT::CodeT error) {
		ExceptionT::BadInputValue(caller, "exception %d validating list \"%s\"", 
			error, source.Name().Pointer());
	}
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
	fReferences.Dimension(0);
	fReferencesOccur.Dimension(0);
}
