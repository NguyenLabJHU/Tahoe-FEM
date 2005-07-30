/*====================================================================
 * 
 * File: xerbla.c 
 *
 * non-standard ANSI utilities
 *
 *====================================================================*/

#include <stdio.h>
#include <stdlib.h>

#include "f2c.h"

/*  XERBLA  is an error handler for the LAPACK routines. */
/*  It is called by an LAPACK routine if an input parameter has an */
/*  invalid value.  A message is printed and execution stops. */
/* int xerbla(char* srname, long int* info) */
integer xerbla(char* srname, integer* info)
{
	printf("** On entry to %s parameter number %i had an illegal value**\n",
		srname,info);
	fflush(stdout);	

	/* terminate execution - raise(SIGABRT) to catch */	
	exit(-1);
	
	return 0;
}
