/* $Id: NewMemory.h,v 1.6 2002-07-02 19:56:55 cjkimme Exp $ */

/* memory allocation function to handle platform dependent
 * differences associated with memory allocation failures and
 * system checks for the amount of free memory */

#ifndef _NEW_MEMORY_H_
#define _NEW_MEMORY_H_

#include <iostream.h>

/* environment */
#include "Environment.h"
#include "ExceptionCodes.h"

/* new that throws bad_alloc on fail instead of returning NULL */
#ifdef __NEW_THROWS__
#include <new.h>
#endif

/* certain CPLANT's don't like to hit memory limits */
#if defined(__DELMAR__) || defined(__ASILOMAR__)
#define _SPACE_CHECK_

/* memory check methods */
#define MEM_CHECK_SYSINFO      0
#define MEM_CHECK_PROC_MEMINFO 1

#define MEM_CHECK_METHOD MEM_CHECK_PROC_MEMINFO

/* size parameters */

namespace Tahoe {

const int min_space_check_size = 100; /* skip check for very small arrays */
const unsigned long min_free_memory = 5*1000000;   /*  ~5MB */
const unsigned long max_used_memory = 210*1000000; /* 210MB */

/* headers needed for memory check implementation */
#if (MEM_CHECK_METHOD == MEM_CHECK_SYSINFO)
#include <sys/sysinfo.h>
#else /* (MEM_CHECK_METHOD == MEM_CHECK_PROC_MEMINFO) */
#include <fstream.h>
#endif /* (MEM_CHECK_METHOD == MEM_CHECK_SYSINFO) */

#endif /* defined(__DELMAR__) || defined(__ASILOMAR__) */

#ifdef _SPACE_CHECK_
#if (MEM_CHECK_METHOD == MEM_CHECK_SYSINFO)
/* assess free memory based on the return from sysinfo. The
 * free memory returned here does not include memory used for
 * disk cache */
inline bool HasFreeMemory(unsigned long size)
{
	/* get sys info structure */
	struct sysinfo s_info;
	sysinfo(&s_info);

	/* check free memory */
	return s_info.freeram - size > min_free_memory;
};
#else
/* assess free memory based on the contents of /proc/meminfo:

        total:    used:    free:  shared: buffers:  cached:
Mem:  260939776 55083008 205856768 23363584   516096 23576576

* the total free memory is taken as: free + buffers + cached */
inline bool HasFreeMemory(unsigned long size)
{
	/* buffer */
	char buffer[255];
	
	/* meminfo file */
	ifstream in("/proc/meminfo");

	/* clear 1st line and 1st word on second line */
	in.getline(buffer, 254);
	in >> buffer;
	unsigned long total, used, free, shared, buffers, cached;
	total = used = free = shared = buffers = cached = 0;
	in >> total >> used >> free >> shared >> buffers >> cached;
	if (free < 1 || buffers < 1 || cached < 1)
	{
		cout << "::HasFreeMemory: error reading \"/proc/meminfo\"" << endl;
		throw eGeneralFail;
	}

//DEBUG
#if 1
cout << "  max: " << max_used_memory << '\n';
cout << " used: " << used << '\n';
cout << " free: " << free << '\n';
cout << " buff: " << buffers << '\n';
cout << " cach: " << cached << '\n';
cout << " used - cached: " << used - cached << endl;
#endif	


//NOTE - cannot subtract unsigned long's, will result in add!

	/* check free memory */
//bool has_free = (cached > used) ? true : used > cached + max_used_memory;
bool has_free = free + buffers + cached > min_free_memory + size;
//bool has_free = (used + size) < max_used_memory;
//DEBUG
#if 1
if (!has_free)
  cout << "FAIL: allocation would violate min_free_memory" << endl;
#endif
	return has_free;
};
#endif /* defined(__DELMAR__) || defined(__ASILOMAR__) */
#endif /* _SPACE_CHECK_ */

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

} // namespace Tahoe 
#endif /* _NEW_MEMORY_H_ */
