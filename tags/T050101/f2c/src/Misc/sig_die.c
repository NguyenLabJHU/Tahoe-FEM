#include "stdio.h"
#include "signal.h"


#if defined(TPM_F2C) || defined(SPM_F2C) || defined(CW_F2C_MAC) || defined(CW_F2C_WIN32)
#include <setjmp.h>
extern jmp_buf gRecoverToConsole;					/* IMT 13Aug97  For recovery to console */
#endif 


#ifndef SIGIOT
#ifdef SIGABRT
#define SIGIOT SIGABRT
#endif
#endif

#ifdef KR_headers
void sig_die(s, kill) register char *s; int kill;
#else
#include "stdlib.h"
#ifdef __cplusplus
extern "C" {
#endif
 extern void f_exit(void);

void sig_die(register char *s, int kill)
#endif
{
	/* print error message, then clear buffers */
	fprintf(stderr, "%s\n", s);

/* PAK 09May00 just use standard exit - missing definition of gRecoverToConsole */
/* IMT 13Aug97  MacOS recovers to console */
/* #if defined(TPM_F2C) || defined(SPM_F2C) || defined(CW_F2C_MAC) || defined(CW_F2C_WIN32) */
/*	fflush(stderr); */
/*	fflush(stdout); */
/*	longjmp( gRecoverToConsole, 1 ); */
/* #else */


	if(kill)
		{
		fflush(stderr);
		f_exit();
		fflush(stderr);
		/* now get a core */
#ifdef SIGIOT
		signal(SIGIOT, SIG_DFL);
#endif
		abort();
		}
	else {
#ifdef NO_ONEXIT
		f_exit();
#endif
		exit(1);
		}

/* TPM, SPM, and CW */
/* #endif */

	}
#ifdef __cplusplus
}
#endif
