/* $Id: Environment.h,v 1.12 2005-04-30 21:14:46 paklein Exp $ */
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
#if (__MWERKS__ > 0) && (__MWERKS__ < 0x7000) /* later version predefine */
using namespace std;
#endif
#endif
#endif

/* using Metrowerks Standard Library */
#ifdef __MWERKS__
#if (__MWERKS__ == 1) /* older versions just define as 0/1 */
#ifdef __mslGlobals_h
#define _MW_MSL_
#endif /* __mslGlobals_h */
#else /* later versions all use MSL */
#define _MW_MSL_
#endif
#endif /* __MWERKS__ */

/* compiler supports RTTI */
#ifdef __SUNPRO_CC
/* 	#define __NO_RTTI__ -> v4.2, not needed for v5.0? */
using namespace std;
#undef _MW_MSL_
#endif
/* NOTE: v4.2 support RTTI, but could not get it to work if the pointer */
/*       I was casting was the pointer to a purely virtual base class   */
/*       type. PAK (01/29/1999)                                           */

/* failure of new throws bad_alloc */
#if defined(__DEC__) || defined(__ALASKA__) || defined(_MW_MSL_)
#define __NEW_THROWS__
#endif

/* GNU */
#ifdef __GNUC__

/* version 3.x */
#if (__GNUC__ == 3) /* GCC predefined macro */
#define __GCC_3__
using namespace std;
#endif

/* version 4.x */
#if (__GNUC__ == 4) /* GCC predefined macro */
#define __GCC_4__
using namespace std;
#ifndef NEW_STATIC_TEMPLATE_SYNTAX
#define NEW_STATIC_TEMPLATE_SYNTAX
#endif
#endif

#endif /* __GNU__ */

/* IBM XL C/C++ for OS X */
#if defined(__DARWIN__) && defined(__XL__)
#define __GCC_3__   /* uses GCC 3.x headers */
using namespace std;
#ifndef NEW_STATIC_TEMPLATE_SYNTAX
#define NEW_STATIC_TEMPLATE_SYNTAX
#endif
#endif

/* explicit definitions of static template data */
#if defined(NEW_STATIC_TEMPLATE_SYNTAX)
#define DEFINE_TEMPLATE_STATIC template<>
#else
#define DEFINE_TEMPLATE_STATIC /* */
#endif

/* wrapper for dynamic casts */
#define TB_DYNAMIC_CAST(type_name, arg) TB_CAST_MACRO(type_name, arg)
#ifndef __NO_RTTI__
#define TB_CAST_MACRO(type_name, arg) dynamic_cast<type_name>(arg)
#else
#define TB_CAST_MACRO(type_name, arg) (type_name) arg
#endif

#endif /* _ENVIRONMENT_H_ */
