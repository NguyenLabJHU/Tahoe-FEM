/* $Id: fstreamT.h,v 1.4.4.1 2002-06-27 18:01:08 cjkimme Exp $ */
/* created: paklein (12/30/2000) */

#ifndef _FSTREAM_T_H_
#define _FSTREAM_T_H_

#include "ifstreamT.h"
#include "ofstreamT.h"

/** interface for stream utilities */

namespace Tahoe {

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
};

} // namespace Tahoe 
#endif /* _FSTREAM_T_H_ */
