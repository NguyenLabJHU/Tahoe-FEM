/* $Id: ParameterListT.cpp,v 1.1 2002-09-06 05:51:28 paklein Exp $ */
#include "ParameterListT.h"

using namespace Tahoe;

/* add parameter */
bool ParameterListT::AddParameter(const ParameterT& param, OccurrenceT occur)
{
	/* scan name */
	for (int i = 0; i < fParameters.Length(); i++)
		if (fParameters[i].Name() == param.Name())
			return false;
	
	/* add if no matches */
	fParameters.Append(param);
	fParametersOccur.Append(occur);
	return true;
}
