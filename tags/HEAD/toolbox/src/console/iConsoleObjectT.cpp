/* $Id: iConsoleObjectT.cpp,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (12/21/2000)                                          */
/* iConsoleObjectT.cpp                                                    */

#include "iConsoleObjectT.h"
#include <ctype.h>

/* array behavior */
const bool ArrayT<iConsoleObjectT*>::fByteCopy = true;

/* constructor */
iConsoleObjectT::iConsoleObjectT(void):
	fSuper(NULL),
	fSubs(0, true)
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
