/* Metrowerks Standard Library
 * Copyright © 1995-2002 Metrowerks Corporation.  All rights reserved.
 *
 *  $Date: 2002-10-19 01:21:28 $ 
 *  $Revision: 1.2 $ 
 */

/*	This file contains the code that gets precompiled.	*/

#if __PALMOS_TRAPS__

	#include <ansi_prefix.Palm_OS.h>

#elif macintosh
	/*	option 1
	 *	to have ansi_prefix not include MacHeaders set MSL_USE_PRECOMPILED_HEADERS to something other than 0 or 1
	 *	uncommenting the following line will achieve this.
	 */

	/* #define MSL_USE_PRECOMPILED_HEADERS 2 */

	/* 	option 2
	 *	to generate MacHeaders as part of the MSL pch leave everything as is.
	 */

	#ifndef MSL_USE_PRECOMPILED_HEADERS
	#define MSL_USE_PRECOMPILED_HEADERS 1
	#endif
	
	#include <ansi_prefix.mac.h>

#elif __MACH__

	#if !_MSL_USING_MW_C_HEADERS
		#include <ansi_prefix.mach.h>
	#endif

#elif __INTEL__ && !__BEOS__

	#include <ansi_prefix.win32.h>

#elif __PPC_EABI__

	/* prefix file is in Language preference panel */
	
#else
	#error "OS currently unsupported"
#endif

#include <iosfwd>
// Support
#include <exception>
#include <new>
#include <limits>
#include <typeinfo>
// Diagnostics
#include <stdexcept>
// Iterators
#include <iterator>
// Utilities
#include <functional>
#include <memory>
#include <utility>
// Algorithms
#include <algorithm>
// Strings
#include <string>
// Containers
#include <bitset>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <queue>
#include <stack>
// Localization
#include <locale>
// Input/Output
#include <ios>
#include <streambuf>
#include <istream>
#include <ostream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
// Numerics
//#include <complex>
//#include <numeric>
//#include <valarray>
// Extensions
//#include <hash_map>
//#include <hash_set>
//#include <slist>
