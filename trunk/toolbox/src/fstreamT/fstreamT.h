/* $Id: fstreamT.h,v 1.2 2002-01-05 06:55:01 paklein Exp $ */
/* created: paklein (12/30/2000) */

#ifndef _FSTREAM_T_H_
#define _FSTREAM_T_H_

#include "ifstreamT.h"
#include "ofstreamT.h"

/** interface for stream utilities */
class fstreamT
{
  public:

	/** C-C++(2.4.6)/MSL doesn't to path's right. temporarily needed to
	 * remove "ups" and "downs" from pathnames. Used only for Mac */
	static void FixPath(const char* path_old, StringT& path);
};

#endif /* _FSTREAM_T_H_ */
