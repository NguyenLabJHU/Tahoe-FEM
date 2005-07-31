#ifdef KR_headers
extern FILE *fdopen();
#else
#ifdef MSDOS
#include "io.h"
#ifndef WATCOM
#define close _close
#define creat _creat
#define open _open
#define read _read
#define write _write
#endif /*WATCOM*/
#endif /*MSDOS*/
#ifdef __cplusplus
extern "C" {
#endif
#ifndef MSDOS
#ifdef OPEN_DECL
extern int creat(const char*,int), open(const char*,int);
#endif

/* IMT 14Aug97 Use CW provided unix prototypes (proto for read() clashes on 3rd argument */
/* IMT 01Oct98 Added CW_F2C_WIN32 */
#if defined(CW_F2C_MAC) || defined(CW_F2C_WIN32) || defined(MPW_CW_F2C) 	/* [ */

/* IMT 01Oct98 Defining VOID breaks some of the Win32 headers included via <unix.h> */
#if defined(CW_F2C_WIN32) 	/* [ */
	#undef VOID
#endif	/* ] CW_F2C_WIN32 */

#include <unix.h>

/* IMT 01Oct98 Make sure <unix.h> really did define VOID */
#if defined(CW_F2C_WIN32) && !defined(VOID)		/* [ */
	#define VOID void
#endif	/* ] CW_F2C_WIN32 && VOID */

#else	/* ] [ */
extern int close(int);
extern int read(int,void*,size_t), write(int,void*,size_t);
extern int unlink(const char*);
#endif	/* ] */


#ifndef _POSIX_SOURCE
#ifndef NON_UNIX_STDIO
extern FILE *fdopen(int, const char*);
#endif
#endif
#endif /*KR_HEADERS*/

extern char *mktemp(char*);

#ifdef __cplusplus
	}
#endif
#endif

#include "fcntl.h"

#ifndef O_WRONLY
#define O_RDONLY 0
#define O_WRONLY 1
#endif
