/* $Id: Environment.h,v 1.1.1.1 2001-01-25 20:56:28 paklein Exp $ */
/* created: paklein (02/10/1997)                                          */
/* Environment.h                                                          */
/* defining environment-specific preprocessor symbols and options         */

#ifndef _ENVIRONMENT_H_
#define _ENVIRONMENT_H_

/*************************************************************************/
/**************************** platform names *****************************/
/*************************************************************************/

/* CodeWarrior */
#ifdef __MWERKS__
#ifdef __INTEL__
#define _WINNT_
#else
#define _MACOS_
#endif	
#endif

/* Visual C++ */
#ifdef _MSC_VER
#define _WINNT_
#pragma warning(disable:4068) //disable unknown MWERKS pragma warnings
#endif

/* something UNIX */
#ifndef _MACOS_
#ifndef _WINNT_
#define _UNIX__
#endif
#endif

/*************************************************************************/
/************************* error checking code ***************************/
/*************************************************************************/

#ifndef __MWERKS__ /* for compilation outside CodeWarrior */
/* compiler options */
#ifdef NDEBUG
#define extended_errorcheck	0	/* no error checking */
#else
#define extended_errorcheck	1	/* peform error checking */
#endif /* NDEBUG */
#define __option(x) x
#endif /* __MWERKS__ */

/*************************************************************************/
/*********************** language support options ************************/
/*************************************************************************/

/* namespaces with */
/*              CWPro > 3                           CWPro >= 5.3 ? */
#if defined(MSIPL_USING_NAMESPACE) || defined(_MSL_USING_NAMESPACE)
#ifdef __cplusplus
using namespace std;
#endif
#endif

/* using Metrowerks Standard Library */
#ifdef __mslGlobals_h
#define _MW_MSL_
#endif

/* compiler supports RTTI */
#ifdef __SUNPRO_CC
/* 	#define __NO_RTTI__ -> v4.2, not needed for v5.0? */
using namespace std;
#endif
/* NOTE: v4.2 support RTTI, but could not get it to work if the pointer */
/*       I was casting was the pointer to a purely virtual base class   */
/*       type. PAK (01/29/1999)                                           */

/* failure of new throws bad_alloc */
#if defined(__DEC__) || defined(__ALASKA__) || defined(_MW_MSL_)
#define __NEW_THROWS__
#endif

#endif /* _ENVIRONMENT_H_ */
