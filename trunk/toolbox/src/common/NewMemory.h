/* $Id: NewMemory.h,v 1.1 2001-07-22 20:02:16 paklein Exp $ */

/* memory allocation function to handle platform dependent
 * differences associated with memory allocation failures and
 * system checks for the amount of free memory */

#ifndef _NEW_MEMORY_H_
#define _NEW_MEMORY_H_

/* environment */
#include "Environment.h"

/* new that throws bad_alloc on fail instead of returning NULL */
#ifdef __NEW_THROWS__
#include <new.h>
#endif

/* certain CPLANT's don't like to hit memory limits */
#if defined(__DELMAR__) || defined(__ASILOMAR__)
#include <sys/sysinfo.h>
#define _SPACE_CHECK_
const int  min_space_check_size = 100;       /* skip check for very small arrays */
const long min_free_memory      = 5*1000000; /* ~5MB */
#endif

#if _SPACE_CHECK_
inline bool HasFreeMemory(long size)
{
	/* get sys info structure */
	struct sysinfo s_info;
	sysinfo(&s_info);

	/* check free memory */
	return size < s_info.freeram;
};
#endif

template <class TYPE>
TYPE* New(int size)
{
	TYPE* p = NULL;
	if (size > 0)
	{
#ifdef _SPACE_CHECK_
		if (size >= min_space_check_size)
			if (!HasFreeMemory(size*sizeof(TYPE))) return NULL;
#endif
	
#ifdef __NEW_THROWS__
		try { p = new TYPE[size]; }
		catch (bad_alloc) { p = NULL; }
#else
		p = new TYPE[size];
#endif
	}
	return p;
};

/* sysinfo structure for DEC Alpha-Linux */
#if 0
struct sysinfo {
  long uptime;                    /* Seconds since boot */
  unsigned long loads[3];         /* 1, 5, and 15 minute load averages */
  unsigned long totalram;         /* Total usable main memory size */
  unsigned long freeram;          /* Available memory size */
  unsigned long sharedram;        /* Amount of shared memory */
  unsigned long bufferram;        /* Memory used by buffers */
  unsigned long totalswap;        /* Total swap space size */
  unsigned long freeswap;         /* swap space still available */
  unsigned short procs;           /* Number of current processes */
  char _f[22];                    /* Pads structure to 64 bytes */
};
#endif

#endif /* _NEW_MEMORY_H_ */
