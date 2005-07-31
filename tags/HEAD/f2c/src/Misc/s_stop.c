
/* IMT 10Sep95  Declare jump buffer used to recover from exception exits & aborts */
#if defined(TPM_F2C) || defined(SPM_F2C) || defined(CW_F2C_MAC) || defined(CW_F2C_WIN32) 
#include <setjmp.h>
extern jmp_buf gRecoverToConsole;
#endif /* Macintosh C compilers and CW Win32 */


#include "stdio.h"
#include "f2c.h"

#ifdef KR_headers
extern void f_exit();
VOID s_stop(s, n) char *s; ftnlen n;
#else
#undef abs
#undef min
#undef max
#include "stdlib.h"
#ifdef __cplusplus
extern "C" {
#endif
void f_exit(void);

int s_stop(char *s, ftnlen n)
#endif
{
int i;

if(n > 0)
	{
	fprintf(stderr, "STOP ");
	for(i = 0; i<n ; ++i)
		putc(*s++, stderr);
	fprintf(stderr, " statement executed\n");
	}
#ifdef NO_ONEXIT
f_exit();
#endif

/* PAK 09May00 just use standard exit - missing definition of gRecoverToConsole */
/* IMT 10Sep95  Use jump buffer used to recover to console instead of exit() */
/* #if defined(TPM_F2C) || defined(SPM_F2C) || defined(CW_F2C_MAC) || defined(CW_F2C_WIN32) */
/* longjmp( gRecoverToConsole, 1 ); */
/* #else */
exit(0);
/* #endif  */
/* Macintosh compilers and CW Win32 */

	return 0; /* NOT REACHED */
#ifdef __cplusplus

}
#endif
}
