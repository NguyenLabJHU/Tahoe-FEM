/* $Id: fstreamT.h,v 1.7 2002-12-02 09:37:02 paklein Exp $ */
/* created: paklein (12/30/2000) */

#ifndef _FSTREAM_T_H_
#define _FSTREAM_T_H_

#include "ifstreamT.h"
#include "ofstreamT.h"

namespace Tahoe {

/** interface for stream utilities */
class fstreamT
{
  public:

	/** C-C++(2.4.6)/MSL doesn't to path's right. temporarily needed to
	 * remove "ups" and "downs" from pathnames. Used only for Mac */
	static void FixPath(const char* path_old, StringT& path);

	/** work around check. This function returns false if the MW
	 * compiler is older than those released CW7.x, true if it is
	 * a compiler version for which the workaround is needed and
	 * throws eStop if the compiler version is newer. */
	static bool need_MW_workaround(void);
	
	/** returns true if the given file is found */
	static bool Exists(const char* path);

	/** clear input stream to the end of the line. A fix needed because GCC 3.1 has
	 * a bug in getline which does not return when it hits the line delimiter '\n'
	 * \return the number of characters read from the stream including the newline */
	static int ClearLine(istream& in);
};

} // namespace Tahoe 
#endif /* _FSTREAM_T_H_ */
