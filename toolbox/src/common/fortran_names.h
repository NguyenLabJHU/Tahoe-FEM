/* $Id: fortran_names.h,v 1.1 2003-03-09 04:27:57 paklein Exp $ */
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
#elif defined(__ALPHA__)
#ifdef __STDC__
#define FORTRAN_NAME(n_)	n_ ## _
#else
#define FORTRAN_NAME(n_)	n_/**/_
#endif

/* Sun not GNU - single trailing underscore */
#elif defined(__SUN__) && !defined(__GNU__)
#ifdef __STDC__
#define FORTRAN_NAME(n_)	n_ ## _
#else
#define FORTRAN_NAME(n_)	n_/**/_
#endif

/* GNU - double trailing underscore */
#elif defined(__GNU__)
#ifdef __STDC__
#define FORTRAN_NAME(n_)	n_ ## __
#else
#define FORTRAN_NAME(n_)	n_/**/__
#endif

#else
#error "don't know Fortran function/subroutine naming convention for this system!"
#endif

/*****************************************************************************/

#endif	/* _FORTTRAN_NAMES_H_ */
