/* $Id: iConsoleObjectT.cpp,v 1.6 2002-07-05 17:16:01 paklein Exp $ */
/* created: paklein (12/21/2000) */

#include "iConsoleObjectT.h"
#include <ctype.h>

/* array behavior */

using namespace Tahoe;

namespace Tahoe {
const bool ArrayT<iConsoleObjectT*>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
iConsoleObjectT::iConsoleObjectT(void):
	fSuper(NULL),
	fSubs(0)
{
	StringT name("<none>");
	iSetName(name);
}

/* subs control - return true if changed */
bool iConsoleObjectT::iAddSub(iConsoleObjectT& sub)
{
	if (sub.fSuper == NULL && fSubs.AppendUnique(&sub) == 1)
	{
		sub.fSuper = this;
		return true;
	}
	else
		return false;
}

bool iConsoleObjectT::iDeleteSub(iConsoleObjectT& sub)
{
	for (int i = 0; i < fSubs.Length(); i++)
		if (fSubs[i] == &sub)
		{
			fSubs.DeleteAt(i);
			sub.fSuper = NULL;
			return true;
		}
	return false;
}
	
/************************************************************************
* Protected
************************************************************************/

/* set name string */
void iConsoleObjectT::iSetName(const StringT& name)
{
	fName = name;
	char* str = fName;
	for (int i = 0; i < strlen(str); i++)
		if (isspace(str[i]))
			str[i] = '_';
}
