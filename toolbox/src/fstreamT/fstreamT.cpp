/* $Id: fstreamT.cpp,v 1.11 2005-01-31 07:07:51 paklein Exp $ */

#include "fstreamT.h"
#include "Environment.h"
#include "ExceptionT.h"
#include "ifstreamT.h"

using namespace Tahoe;

/* initialize static variables */
StringT fstreamT::fRoot;

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
if (__MWERKS__ < 0x2402 || __MWERKS__ > 0x3200) /* old versions are OK */
	return false;
else if(__MWERKS__ <= 0x2407)
	return true;
else /* stop */
{
	cout << "ifstreamT::need_MW_workaround: __MWERKS__ > 0x2407. Still need fix?" << endl;
	return true;	
}
#else
	return false;
#endif
}

/* set file name string - does not change stream */
void fstreamT::set_filename(const char* name)
{
	/* store file name */
	fFileName = name;
	fFileName.ToNativePathName();
}

/* returns true if the given file is found */
bool fstreamT::Exists(const char* path)
{
	ifstreamT file(path);
	return file.is_open() == 1;
}

/* clear input stream to the end of the line */
int fstreamT::ClearLine(istream& in)
{
	/* clear to end of line */
	int ct = 0;
	char a = in.get();
	while (in.good() && a != '\n') {
		ct++;
		a = in.get();
	}
	return ct;
}

/* set the file root, pass NULL to clear */
void fstreamT::SetRoot(const char* root)
{
	if (root)
		fRoot = root;
	else
		fRoot.Clear();
}
