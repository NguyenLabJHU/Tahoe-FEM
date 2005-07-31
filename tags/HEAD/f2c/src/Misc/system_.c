/* f77 interface to system routine */

#include "f2c.h"

#ifdef KR_headers
extern char *F77_aloc();

 integer
system_(s, n) register char *s; ftnlen n;
#else
#undef abs
#undef min
#undef max
#include "stdlib.h"
extern char *F77_aloc(ftnlen, char*);

 integer
system_(register char *s, ftnlen n)
#endif
{

/* IMT 13Aug97  Macintosh doesn't do system commands */

#if defined(SPM_F2C) || defined(TPM_F2C) || defined(CW_F2C_MAC) || defined(MPW_CW_F2C)

	/* 
		One could probably implement an AE version of system command 
		processing.  I'm too lazy to do it.
	*/
	
	return 0;
			
#else	

	char buff0[256], *buff;
	register char *bp, *blast;
	integer rv;

	buff = bp = n < sizeof(buff0)
			? buff0 : F77_aloc(n+1, "system_");
	blast = bp + n;

	while(bp < blast && *s)
		*bp++ = *s++;
	*bp = 0;
	rv = system(buff);
	if (buff != buff0)
		free(buff);
	return rv;
		
#endif	/* TPM_F2C, SPM_F2C, CW_F2C_MAC and MPW_CW_F2C */
}
