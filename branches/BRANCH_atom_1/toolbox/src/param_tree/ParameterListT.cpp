/* $Id: ParameterListT.cpp,v 1.3 2002-11-18 09:59:03 paklein Exp $ */
#include "ParameterListT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
const bool ArrayT<ParameterListT>::fByteCopy = false;
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
