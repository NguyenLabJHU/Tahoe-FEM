/* $Id: fortran_names.h,v 1.7 2004-11-12 20:58:32 paklein Exp $ */
#ifndef _FORTRAN_NAMES_H_
#define _FORTRAN_NAMES_H_
/*
 * Fortran <--> C/C++ interfacing stuff
 * revised from fortran.h from: 
 *    http://www.aei.mpg.de/~jthorn/c2f.html#fortran.h
 */

/*
 * Names of Fortran routines are often altered by the compiler/loader.  The
 * following macro should be used to call a Fortran routine from C code, i.e.
 *	call sgefa(...)			-- Fortran code
 *	FORTRAN_NAME(sgefa)(...);	-- C code to do the same thing
 *
 * Unfortunately, the "alterations" are generally at the token level, and this
 * can't be done portably in pre-ANSI C.  In ANSI C, the preprocessor "token
 * splicing" facility is designed to handle just this sort of thing, but in
 * pre-ANSI C we have to use rather ugly system-dependent hacks of the sort
 * exemplified below.
 *
 * C code should reference Fortran names in lower case.
 */

/* SGI - single trailing underscore */
#if defined(__SGI__)
#ifdef __STDC__
#define FORTRAN_NAME(n_)	n_ ## _
#else
#define FORTRAN_NAME(n_)	n_/**/_
#endif

/* DEC - single trailing underscore */
#elif defined(__DEC__)
#define FORTRAN_NAME(n_)	n_ ## _

/* AIX - no mangling */
#elif defined(__AIX__)
#define FORTRAN_NAME(n_)	n_

/* Intel compilers */
#elif defined(__INTEL_CC__)
#define FORTRAN_NAME(n_)        n_ ## _

/* Sun not GNU - single trailing underscore */
#elif defined(__SUN__) && !defined(__GNUC__)
#ifdef __STDC__
#define FORTRAN_NAME(n_)	n_ ## _
#else
#define FORTRAN_NAME(n_)	n_/**/_
#endif

/* GNU - double trailing underscore */
#elif defined(__GNUC__) /* GCC predefined macro */
#ifdef __STDC__
#ifdef NO_SECOND_UNDERSCORE
#define FORTRAN_NAME(n_)	n_ ## _
#else
#define FORTRAN_NAME(n_)	n_ ## __
#endif
#else
#ifdef NO_SECOND_UNDERSCORE
#define FORTRAN_NAME(n_)	n_/**/_
#else
#define FORTRAN_NAME(n_)	n_/**/__
#endif
#endif

/* Metrowerks - assume GNU linkage */
#elif defined(__MWERKS__)
#define FORTRAN_NAME(n_)	n_ ## __

#else
#error "don't know Fortran function/subroutine naming convention for this system!"
#endif

/*****************************************************************************/

#endif	/* _FORTTRAN_NAMES_H_ */
