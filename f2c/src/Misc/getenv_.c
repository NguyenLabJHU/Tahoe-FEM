#include "f2c.h"

/*
 * getenv - f77 subroutine to return environment variables
 *
 * called by:
 *	call getenv (ENV_NAME, char_var)
 * where:
 *	ENV_NAME is the name of an environment variable
 *	char_var is a character variable which will receive
 *		the current value of ENV_NAME, or all blanks
 *		if ENV_NAME is not defined
 */

#ifdef KR_headers
VOID getenv_(fname, value, flen, vlen) char *value, *fname; ftnlen vlen, flen;
#else
void getenv_(char *fname, char *value, ftnlen flen, ftnlen vlen)
#endif
{

/* IMT 13Aug97  Macintosh doesn't have any environment variables */

#if defined(SPM_F2C) || defined(TPM_F2C) || defined(CW_F2C_MAC) || defined(MPW_CW_F2C)

	while( --vlen >= 0 )				/* On Mac just fill it with blanks */
		*value++ = ' ';
		
#else	

extern char **environ;
register char *ep, *fp, *flast;
register char **env = environ;

flast = fname + flen;
for(fp = fname ; fp < flast ; ++fp)
	if(*fp == ' ')
		{
		flast = fp;
		break;
		}

while (ep = *env++)
	{
	for(fp = fname; fp<flast ; )
		if(*fp++ != *ep++)
			goto endloop;

	if(*ep++ == '=') {	/* copy right hand side */
		while( *ep && --vlen>=0 )
			*value++ = *ep++;

		goto blank;
		}
endloop: ;
	}

blank:
	while( --vlen >= 0 )
		*value++ = ' ';
		
#endif	/* TPM_F2C, SPM_F2C, CW_F2C_MAC and MPW_CW_F2C */
}
