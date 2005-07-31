/*====================================================================
 * 
 * File: utilsx.c 
 *
 * non-standard ANSI utilities
 *
 *====================================================================*/

#include "utilsx.h"

#include <string.h>
#include <stdlib.h>

/* returns copy of str - dynamically allocated */
char* az_strdup(const char* str)
{
	/* allocate memory */
	char* new_str = (char*) malloc(sizeof(char)*(strlen(str) + 1));

	/* copy */
	return strcpy(new_str, str);
}
