
/*****************************************************
*													 *
*	Special	macros required	to correctly build the	 *
*	the	CodeWarrior	MacOS version of libf2c.		 *
*													 *
*****************************************************/


#pragma once


#ifndef LIBF2C_CW_H

#define LIBF2C_CW_H


#if macintosh

	#define NON_UNIX_STDIO   	/* Force ANSI standard IO calls */
	#define _POSIX_SOURCE    	/* Force use of tmpnam() instead of mktemp() */
	#define USE_CLOCK			/* Force use of ANSI clock() routines */

	#define IMT_F2C_BUG_FIX		/* Trigger bug fixes that IMT implemented */
	#define CW_F2C_MAC			/* Trigger CodeWarrior MacOS specific library modifications */

//	#define Allow_TYQUAD 		/* Include Quad Integer support */

#elif __INTEL__

	#define NON_UNIX_STDIO   	/* Force ANSI standard IO calls */
	#define _POSIX_SOURCE    	/* Force use of tmpnam() instead of mktemp() */
	#define USE_CLOCK			/* Force use of ANSI clock() routines */

	#define IMT_F2C_BUG_FIX		/* Trigger bug fixes that IMT implemented */
	#define CW_F2C_WIN32		/* Trigger CodeWarrior Win32 specific library modifications */

	#define Allow_TYQUAD 		/* Include Quad Integer support */

	#define fileno	_fileno		/* Overcome bug in CW Pro 5 MSL release */

#endif	/* Target CPU */


#endif /* LIBF2C_CW_H */




