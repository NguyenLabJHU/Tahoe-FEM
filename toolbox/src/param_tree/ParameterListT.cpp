/* $Id: ParameterListT.cpp,v 1.4 2003-04-22 18:32:16 paklein Exp $ */
#include "ParameterListT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
const bool ArrayT<ParameterListT>::fByteCopy = false;
const bool ArrayT<ParameterListT*>::fByteCopy = true;
const bool ArrayT<ParameterListT::OccurrenceT>::fByteCopy = false;
}

/* add parameter */
bool ParameterListT::AddParameter(const ParameterT& param, OccurrenceT occur)
{
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
	/* scan name */
	for (int i = 0; i < fParameterLists.Length(); i++)
		if (fParameterLists[i].Name() == param_list.Name())
			return false;
	
	/* add if no matches */
	fParameterLists.Append(param_list);
	fParameterListsOccur.Append(occur);
	return true;
}

/* add a reference */
bool ParameterListT::AddReference(const StringT& ref, OccurrenceT occur)
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
const ParameterListT* ParameterListT::List(const StringT& name) const
{
	/* search list */
	for (int i = 0; i < fParameterLists.Length(); i++)
		if (fParameterLists[i].Name() == name)
			return fParameterLists.Pointer(i);

	/* fail */
	return NULL;
}

/* return the non-const pointer to the given parameter or NULL if the list is not found */
const ParameterT* ParameterListT::Parameter(const StringT& name) const
{
	/* search list */
	for (int i = 0; i < fParameters.Length(); i++)
		if (fParameters[i].Name() == name)
			return fParameters.Pointer(i);

	/* fail */
	return NULL;
}

void ParameterListT::GetParameter(const StringT& name, int& a) const
{
	const char caller[] = "ParameterListT::GetParameter";
	const ParameterT* parameter = Parameter(name);
	if (!parameter)
		ExceptionT::GeneralFail(caller, "parameter \"%s\" not found", name.Pointer());
	a = *parameter;
}

void ParameterListT::GetParameter(const StringT& name, double& a) const
{
	const char caller[] = "ParameterListT::GetParameter";
	const ParameterT* parameter = Parameter(name);
	if (!parameter)
		ExceptionT::GeneralFail(caller, "parameter \"%s\" not found", name.Pointer());
	a = *parameter;
}

void ParameterListT::GetParameter(const StringT& name, StringT& a) const
{
	const char caller[] = "ParameterListT::GetParameter";
	const ParameterT* parameter = Parameter(name);
	if (!parameter)
		ExceptionT::GeneralFail(caller, "parameter \"%s\" not found", name.Pointer());
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

	/* set name */
	SetName(source.Name());

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
		bool add_parameter = false;
		ParameterT new_parameter(parameter.Type(), parameter.Name());
		int count = 0;
		for (int j = 0; j < source_parameters.Length(); j++) {
				
			/* name match */
			if (source_parameters[j].Name() == parameter.Name()) {
			
				/* check occurrence */
				if (count > 0 && (occurrence == Once || occurrence == ZeroOrOnce))
					ExceptionT::BadInputValue(caller, "parameter \"%s\" cannot appear in \"%s\" more than once",
						parameter.Name().Pointer(), source.Name().Pointer());
									
				/* get value */
				if (source_parameters[j].Type() == parameter.Type())
					new_parameter = source_parameters[j];
				else if (source_parameters[j].Type() == ParameterT::String)
					new_parameter.FromString(source_parameters[j]);
				else
					ExceptionT::BadInputValue(caller, "source for \"%s\" must have type %d or %d: %d", 
						parameter.Name().Pointer(), ParameterT::String, parameter.Type());

				/* increment count */
				count++;
				
				/* add it */
				add_parameter = true;
			}
		}
				
		/* look for default value */
		if (count == 0 && (occurrence == Once || occurrence == OnePlus)) {
			const ValueT* default_value = parameter.Default();
			if (default_value) {
				new_parameter = *default_value;
				add_parameter = true;
			}
			else
				ExceptionT::BadInputValue(caller, 
					"required value \"%s\" is missing from \"%s\" and has no default value",
					parameter.Name().Pointer(), source.Name().Pointer());
		}
		
		/* add it */
		if (add_parameter) {
			
			/* check limits */
			const ArrayT<LimitT>& limits = parameter.Limits();

			/* add to the list */
			AddParameter(new_parameter);
		}
	}

	/* process sub-lists */
	const ArrayT<ParameterListT>& lists = description.Lists();
	const ArrayT<ParameterListT::OccurrenceT>& lists_occur = description.ListOccurrences();
	const ArrayT<ParameterListT>& source_lists = source.Lists();

	/* validate parameters */
	for (int i = 0; i < lists.Length(); i++)
	{
		const ParameterListT& list = lists[i];
		OccurrenceT occurrence = lists_occur[i];

		/* search through sublists */
		int count = 0;
		for (int j = 0; j < source_lists.Length(); j++) {
				
			/* name match */
			if (source_lists[j].Name() == list.Name()) {
			
				/* check occurrence */
				if (count > 0 && (occurrence == Once || occurrence == ZeroOrOnce))
					ExceptionT::BadInputValue(caller, "parameter list \"%s\" cannot appear in \"%s\" more than once",
						list.Name().Pointer(), source.Name().Pointer());
									
				/* create validated sub-list */
				ParameterListT new_list;
				new_list.Validate(source_lists[j], list);

				/* add it */
				AddList(new_list);

				/* increment count */
				count++;
			}
		}
		
		/* check for required lists */
		if (count == 0 && (occurrence == Once || occurrence == OnePlus))
			ExceptionT::BadInputValue(caller,
				"required parameter list \"%s\" is missing from \"%s\"",
				list.Name().Pointer(), source.Name().Pointer());		
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
