/* This gives the effect of

	subroutine exit(rc)
	integer*4 rc
	stop
	end

 * with the added side effect of supplying rc as the program's exit code.
 */


/* IMT 10Sep95  Declare jump buffer used to recover from exception exits & aborts */
#if defined(TPM_F2C) || defined(SPM_F2C) || defined(CW_F2C_MAC) || defined(CW_F2C_WIN32)
#include <setjmp.h>
extern jmp_buf gRecoverToConsole;
#endif /* Macintosh C compilers and CW Win32 */


#include "f2c.h"
#undef abs
#undef min
#undef max
#ifndef KR_headers
#include "stdlib.h"
#ifdef __cplusplus
extern "C" {
#endif
extern void f_exit(void);
#endif

 void
#ifdef KR_headers
exit_(rc) integer *rc;
#else
exit_(integer *rc)
#endif
{
#ifdef NO_ONEXIT
	f_exit();
#endif


/* IMT 15Aug97  USe jump buffer used to recover to console instead of exit() */
#if defined(TPM_F2C) || defined(SPM_F2C) || defined(CW_F2C_MAC) || defined(CW_F2C_WIN32)
longjmp( gRecoverToConsole, (*rc ? *rc : 1) );
#else
	exit(*rc);
#endif /* Macintosh compilers and CW Win32 */

}
	
#ifdef __cplusplus
}
#endif
