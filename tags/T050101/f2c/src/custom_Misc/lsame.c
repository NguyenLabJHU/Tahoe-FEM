/*====================================================================
 * 
 * File: lsame.c 
 *
 *====================================================================*/

#include "f2c.h"
#undef abs /* f2c definition conflicts with stdlib */

#include <stdlib.h>
#include <ctype.h>

/* returns 1 of a is same as b regardless of case */
integer lsame(char *a, char *b)
{
	char A = toupper(*a);
	char B = toupper(*b);

	return (A==B);
}
