/* $Id: fstreamT.cpp,v 1.2 2002-01-07 20:40:13 paklein Exp $ */

#include "fstreamT.h"

#include "Environment.h"
#include "ExceptionCodes.h"

/* temporary */
void fstreamT::FixPath(const char* path_old, StringT& path)
{
	/* workaround for CW7.x bug */
	if (need_MW_workaround())
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
	else
		path = path_old;
}

bool fstreamT::need_MW_workaround(void)
{
#if defined(__MWERKS__) && !defined(__MACH__)
if (__MWERKS__ < 0x2402) /* old versions are OK */
	return false;
else if(__MWERKS__ <= 0x2406)
	return true;
else /* stop */
{
	cout << "ifstreamT::need_MW_workaround: __MWERKS__ > 0x2406. Still need fix?" << endl;
	throw eStop;
}
#else
	return false;
#endif
}

