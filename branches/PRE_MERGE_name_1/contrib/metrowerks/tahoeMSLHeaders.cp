/* $Id: tahoeMSLHeaders.cp,v 1.1 2002-04-30 23:11:28 paklein Exp $ */
/* This file contains the code that gets precompiled.	   */
/* It is basically the same as the file "MSLHeaders++.cp"  */
/* included with the CodeWarrior distribution, except that */
/* headers conflicting with tahoe/toolbox definitions have */
/* been omitted */

#if macintosh
	/*	option 1
	 *	to have ansi_prefix not include MacHeaders set MSL_USE_PRECOMPILED_HEADERS to something other than 0 or 1
	 *	uncommenting the following line will achieve this.
	 */

	/*	#define MSL_USE_PRECOMPILED_HEADERS 2 */

	/* 	option 2
	 *	to generate MacHeaders as part of the MSL pch leave everything as is.
	 */

	#ifndef MSL_USE_PRECOMPILED_HEADERS
	#define MSL_USE_PRECOMPILED_HEADERS 1	
	#endif
	
	#include <ansi_prefix.mac.h>

#elif __MACH__

	#include <ansi_prefix.mach.h>

#elif __INTEL__ && !__BEOS__

	#include <ansi_prefix.win32.h>

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
#include <hash_map>
#include <hash_set>
#include <slist>
