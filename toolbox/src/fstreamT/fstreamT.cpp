/* $Id: fstreamT.cpp,v 1.1 2002-01-05 06:55:00 paklein Exp $ */

#include "fstreamT.h"

#include "Environment.h"
#include "ExceptionCodes.h"

/* temporary */
void fstreamT::FixPath(const char* path_old, StringT& path)
{
#if defined(__MWERKS__) && !defined(__MACH__)
if (__MWERKS__ < 0x2402) /* old versions are OK */
	path = path_old;
else if(__MWERKS__ <= 0x2406)
{
	/* copy */
	path = path_old;
	
	/* advance left pointer */
	char* pL = path;
	while (*(pL + 1) == ':') pL++;

	/* simplify */
	bool changed;
	do {
	
		/* scan for :: */
		changed = false;
		for (char* p = pL + 1; !changed && *p != '\0'; p++)
			if (*p == ':')
			{ 
				/* cut */
				if (*(p+1) == ':')
				{
					char* p0 = path;
					path.Delete(pL - p0 + 1, p - p0 + 1);
					changed = true;
				}
				else /* advance left pointer */
					pL = p;
			} 
	} while (changed);
}
else /* stop */
{
	cout << "ifstreamT::FixPath: __MWERKS__ <= 0x2406. Still need fix?" << endl;
	throw eStop;
}
#else
	path = path_old;
#endif
}
