/* $Id: CommunicatorT.cpp,v 1.17 2004-07-27 17:48:22 paklein Exp $ */
#include "CommunicatorT.h"
#include "ExceptionT.h"
#include <iostream.h>
#include "nArrayT.h"
#include "ofstreamT.h"

/* to handle the variable number of arguments in Throw() */
#if defined(__SGI__) || defined(__DELMAR__) || defined(__PGI__) || defined(__JANUS__) || defined(__MWERKS__)
#include <stdio.h>
#include <stdarg.h>
#else
#include <cstdio>
#include <cstdarg>
#endif

/* error checking - these should be set by the compiler macros 
#undef CHECK_MPI_STATUS
#undef CHECK_MPI_RETURN
*/

/* buffer for vsprintf */
static char message_buffer[255];

using namespace Tahoe;

/* initialize static variables */
int CommunicatorT::fCount = 0;
int* CommunicatorT::fargc = NULL;
char*** CommunicatorT::fargv = NULL;

/* create communicator including all processes */
CommunicatorT::CommunicatorT(void):
	fComm(MPI_COMM_NULL),
	fSize(1),
	fRank(0),
	fLogLevel(kUrgent),
	fLog(&cout),
	fLastTime(0)
{
	/* check MPI environment */
	Init();

#ifdef __TAHOE_MPI__
	fComm = MPI_COMM_WORLD;

	int ret = MPI_Comm_size(fComm, &fSize);
#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::CommunicatorT", "MPI_Comm_size failed");
#endif

	ret = MPI_Comm_rank(fComm, &fRank);
#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::CommunicatorT", "MPI_Comm_rank failed");
#endif

#endif
}

/* copy constructor */
CommunicatorT::CommunicatorT(const CommunicatorT& source):
	fComm(MPI_COMM_NULL),
	fSize(1),
	fRank(0),
	fLogLevel(source.fLogLevel),
	fLog(source.fLog),
	fLastTime(source.fLastTime)
{
	/* check MPI environment */
	Init();

#ifdef __TAHOE_MPI__
	/* ducplicate the communicator */
	int ret = MPI_Comm_dup(source.fComm, &fComm);
#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::CommunicatorT", "MPI_Comm_dup failed");
#endif

	/* get rank and size */
	if (fComm != MPI_COMM_NULL) {

		int ret = MPI_SUCCESS;
		ret = MPI_Comm_size(fComm, &fSize);
		ret = MPI_Comm_rank(fComm, &fRank);

#ifdef CHECK_MPI_RETURN
		if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::CommunicatorT", "MPI_Comm_size failed");
		if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::CommunicatorT", "MPI_Comm_rank failed");
#endif
	}
#endif
}

/* split communicator */
CommunicatorT::CommunicatorT(const CommunicatorT& source, int color, int rank_key):
	fComm(MPI_COMM_NULL),
	fSize(1),
	fRank(0),
	fLogLevel(source.fLogLevel),
	fLog(source.fLog),
	fLastTime(source.fLastTime)	
{
	/* check MPI environment */
	Init();

#ifdef __TAHOE_MPI__
	/* split the communicator */
	int ret = MPI_Comm_split(source.fComm, color, rank_key, &fComm);
#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::CommunicatorT", "MPI_Comm_split failed");
#endif

	/* get rank and size */
	if (fComm != MPI_COMM_NULL) {

		int ret = MPI_SUCCESS;
		ret = MPI_Comm_size(fComm, &fSize);
		ret = MPI_Comm_rank(fComm, &fRank);

#ifdef CHECK_MPI_RETURN
		if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::CommunicatorT", "MPI_Comm_size failed");
		if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::CommunicatorT", "MPI_Comm_rank failed");
#endif
	}
#else
#pragma unused(color)
#pragma unused(rank_key)
#endif
}

/* destructor */
CommunicatorT::~CommunicatorT(void)
{
	/* check MPI environment */
	Finalize();
}

void CommunicatorT::Log(LogLevelT priority, const char* caller) const
{
	if (priority >= LogLevel())
	{
		/* log message */
		doLog(caller, NULL);

		/* throw exception */
		if (priority == kFail) ExceptionT::MPIFail(caller);
	}
}

void CommunicatorT::Log(LogLevelT priority, const char* caller, const char* fmt, ...) const
{
	if (priority >= LogLevel())
	{
		/* log message */
		va_list argp;
		va_start(argp, fmt);
		vsprintf(message_buffer, fmt, argp);
		va_end(argp);
		doLog(caller, message_buffer);

		/* throw exception */
		if (priority == kFail) ExceptionT::MPIFail(caller, message_buffer);
	}
}

/* return elapsed time */
double CommunicatorT::Time(void)
{
#ifdef __TAHOE_MPI__
	return MPI_Wtime();
#else
	time_t t;
	time(&t);
	return double(t)/CLOCKS_PER_SEC;
#endif
}

/* (re-)set the logging stream */
void CommunicatorT::SetLog(ostream& log)
{
	Log(kUrgent, "SetLog", "new log stream");
	fLog->flush();
	fLog = &log;
}

void CommunicatorT::SetLog(ofstreamT& log)
{
	Log(kUrgent, "SetLog", "new log stream: %s", log.filename());
	fLog->flush();
	fLog = &log;
}

/* maximum over single integers */
int CommunicatorT::Max(int a) const
{
	Log(kModerate, "Max", "in = %d", a);

#ifdef __TAHOE_MPI__
	if (Size() > 1)
	{
		int b = a;
		int ret = MPI_Allreduce(&b, &a, 1, MPI_INT, MPI_MAX, fComm);

#ifdef CHECK_MPI_RETURN
		if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::Max", "MPI_Allreduce failed");
#endif
	}
#endif

	Log(kModerate, "Max", "out = %d", a);

	return a;
}

/* minimum over single integers */
int CommunicatorT::Min(int a) const
{
	Log(kModerate, "Min", "in = %d", a);

#ifdef __TAHOE_MPI__
	if (Size() > 1)
	{
		int b = a;
		int ret = MPI_Allreduce(&b, &a, 1, MPI_INT, MPI_MIN, fComm);

#ifdef CHECK_MPI_RETURN
		if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::Min", "MPI_Allreduce failed");
#endif
	}
#endif

	Log(kModerate, "Min", "out = %d", a);
	return a;
}

/* minimum over single integers */
int CommunicatorT::Sum(int a) const
{
	Log(kModerate, "Sum", "in = %d", a);

#ifdef __TAHOE_MPI__
	if (Size() > 1)
	{
		int b = a;
		int ret = MPI_Allreduce(&b, &a, 1, MPI_INT, MPI_SUM, fComm);
		
#ifdef CHECK_MPI_RETURN
		if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::Sum", "MPI_Allreduce failed");
#endif
	}
#endif

	Log(kModerate, "Sum", "out = %d", a);
	return a;
}

/* element-by-element sum of a vector */
void CommunicatorT::Sum(const nArrayT<double>& my_values, nArrayT<double>& sum) const
{
	const char caller[] = "CommunicatorT::Sum";

	/* dimension check */
	if (my_values.Length() != sum.Length())
		Log(kFail, caller, "sum %d insufficient size to hold %d", my_values.Length(), sum.Length());

	/* log input */
	if (LogLevel() == kLow)
		Log() << setw(10) << "in:\n" << my_values.wrap(5) << '\n';

	if (Size() == 1)
		sum = my_values;
#ifdef __TAHOE_MPI__
	else if (Size() > 1)
	{
		int ret = MPI_Allreduce((void*) my_values.Pointer(), sum.Pointer(), my_values.Length(), MPI_DOUBLE, MPI_SUM, fComm);
		
#ifdef CHECK_MPI_RETURN
		if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Allreduce failed");
#endif
	}
#endif

	/* log sum */
	if (LogLevel() == kLow)
		Log() << setw(10) << "sum:\n" << sum.wrap(5) << '\n';
}

/* maximum over single doubles */
double CommunicatorT::Max(double a) const
{
	Log(kModerate, "Max", "in = %e", a);

#ifdef __TAHOE_MPI__
	if (Size() > 1)
	{
		double b = a;
		int ret = MPI_Allreduce(&b, &a, 1, MPI_DOUBLE, MPI_MAX, fComm);
		
#ifdef CHECK_MPI_RETURN
		if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::Max", "MPI_Allreduce failed");
#endif
	}
#endif

	Log(kModerate, "Max", "out = %e", a);
	return a;
}

/* minimum over single doubles */
double CommunicatorT::Min(double a) const
{
	Log(kModerate, "Min", "in = %e", a);

#ifdef __TAHOE_MPI__
	if (Size() > 1)
	{
		double b = a;
		int ret = MPI_Allreduce(&b, &a, 1, MPI_DOUBLE, MPI_MIN, fComm);
		
#ifdef CHECK_MPI_RETURN
		if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::Min", "MPI_Allreduce failed");
#endif
	}
#endif

	Log(kModerate, "Min", "out = %e", a);
	return a;
}

/* minimum over single integers */
double CommunicatorT::Sum(double a) const
{
	Log(kModerate, "Sum", "in = %e", a);

#ifdef __TAHOE_MPI__
	if (Size() > 1)
	{
		double b = a;
		int ret = MPI_Allreduce(&b, &a, 1, MPI_DOUBLE, MPI_SUM, fComm);

#ifdef CHECK_MPI_RETURN
		if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::Sum", "MPI_Allreduce failed");
#endif
	}
#endif

	Log(kModerate, "Sum", "out = %e", a);
	return a;
}

/* gather single integer. Called by destination process. */
void CommunicatorT::Gather(int a, nArrayT<int>& gather) const
{
	Log(kModerate, "Gather", "in = %d, destination = %d", a, Rank());

	/* check */
	if (gather.Length() != Size()) ExceptionT::SizeMismatch("CommunicatorT::Gather");

	/* this */
	gather[Rank()] = a;

#ifdef __TAHOE_MPI__
	int ret = MPI_Gather(&a, 1, MPI_INT, gather.Pointer(), 1, MPI_INT, Rank(), fComm);

#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::Gather", "MPI_Gather failed");
#endif

#endif

	/* verbose */
	if (LogLevel() == kLow)
	{
		Log(kLow, "Gather", "destination");
		Log() << gather.wrap(10) << '\n';
	}
}

/* gather single integer. Called by sending processes. */
void CommunicatorT::Gather(int a, int destination) const
{
	const char caller[] = "CommunicatorT::Gather";
	Log(kModerate, caller, "in = %d, destination = %d", a, destination);

	/* check */
	if (destination == Rank())
		ExceptionT::GeneralFail(caller, 
			"destination %d must be different from rank %d", destination, Rank());

#ifdef __TAHOE_MPI__
	int* tmp = NULL;
	int ret = MPI_Gather(&a, 1, MPI_INT, tmp, 1, MPI_INT, destination, fComm);

#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Gather failed");
#endif

#endif
}

/* gather an array of integers. Called by destination process. */
void CommunicatorT::Gather(const nArrayT<int>& my_values, nArrayT<int>& gather, 
	const nArrayT<int>& counts, const nArrayT<int>& displacements) const
{
	const char caller[] = "CommunicatorT::Gather";
	Log(kModerate, caller, "sending %d to %d (self)", my_values.Length(), Rank());
	if (LogLevel() == kLow) Log() << my_values.wrap(10) << '\n';

#if __option(extended_errorcheck)
	if (displacements.Length() != Size()) ExceptionT::SizeMismatch(caller);
	if (counts[Rank()] != my_values.Length()) ExceptionT::SizeMismatch(caller);
	if (counts.Length() != displacements.Length()) ExceptionT::SizeMismatch(caller);
	if (displacements.Last() + counts.Last() > gather.Length()) ExceptionT::SizeMismatch(caller);
		/* assume counts and displacements are OK interms of overlap and assume the
		 * displacements are monotonically increasing */
#endif

#ifdef __TAHOE_MPI__
	int ret = MPI_Gatherv((void*) my_values.Pointer(), my_values.Length(), MPI_INT, 
		gather.Pointer(), (int*) counts.Pointer(), (int*) displacements.Pointer(), MPI_INT, 
		Rank(), fComm);

#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Gatherv failed");
#endif

#else
	/* write my into gather */
	int len = my_values.Length();
	int offset = displacements[Rank()];
	gather.CopyPart(offset, my_values, 0, len);
#endif

	if (LogLevel() == kLow) {
		Log(kLow, caller);
		Log() << gather.wrap(10) << '\n';
	}
}

/* gather an array of integers. Called by sending processes. */
void CommunicatorT::Gather(const nArrayT<int>& my_values, int destination) const
{
	const char caller[] = "CommunicatorT::Gather";
	Log(kModerate, caller, "sending %d to %d", my_values.Length(), destination);
	if (LogLevel() == kLow) Log() << my_values.wrap(10) << '\n';

	/* check */
	if (destination < 0 || destination >= Size()) ExceptionT::SizeMismatch(caller);

#ifdef __TAHOE_MPI__
	int ret = MPI_Gatherv((void*) my_values.Pointer(), my_values.Length(), MPI_INT, 
		NULL, NULL, NULL, MPI_INT, destination, fComm);
		
#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Gatherv failed");
#endif

#endif
}

/* synchronize all processes */
void CommunicatorT::Barrier(void) const
{
	Log(kLow, "Barrier");

#ifdef __TAHOE_MPI__
	int ret = MPI_Barrier(fComm);
#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::Barrier", "MPI_Barrier failed");
#endif
#endif

	Log(kLow, "Barrier");
}

/* gather single integer to all processes. */
void CommunicatorT::AllGather(int a, nArrayT<int>& gather) const
{
	const char caller[] = "CommunicatorT::AllGather";
	Log(kModerate, caller, "in = %d", a);

	/* check */
	if (gather.Length() != Size())
		ExceptionT::SizeMismatch(caller, "buffer length %d does not match size %d", 
			gather.Length(), Size());

	/* this */
	gather[Rank()] = a;

#ifdef __TAHOE_MPI__
	int ret = MPI_Allgather(&a, 1, MPI_INT, gather.Pointer(), 1, MPI_INT, fComm);
#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Allgather failed");
#endif
#endif

//NOTE: need to implement gather with sends and receives for CPLANT and other
//      platforms where Allgather seems to be troublesome.

	/* verbose */
	if (LogLevel() == kLow)
	{
		Log(kLow, caller);
		Log() << gather.wrap(10) << '\n';
	}
}

/* gather multiple values from all */
void CommunicatorT::AllGather(const nArrayT<double>& my, nArrayT<double>& gather) const
{
	const char caller[] = "CommunicatorT::AllGather";
	Log(kModerate, caller, "sending %d", my.Length());
	if (LogLevel() == kLow) Log() << my.wrap(5) << '\n';

	/* check */
	if (my.Length()*Size() != gather.Length()) 
		ExceptionT::SizeMismatch(caller, 
			"buffer length %d does not match size %d", gather.Length(), my.Length()*Size());

#ifdef __TAHOE_MPI__
	int len = my.Length();
	int ret = MPI_Allgather((void*) my.Pointer(), len, MPI_DOUBLE, gather.Pointer(), len, MPI_DOUBLE, fComm);

#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Allgather failed");
#endif

#else
	/* write my into gather */
	int len = my.Length();
	gather.CopyPart(Rank()*len, my, 0, len);
#endif

//NOTE: need to implement gather with sends and receives for CPLANT and other
//      platforms where Allgather seems to be troublesome.

	Log(kModerate, caller, "gathered %d", gather.Length());
	if (LogLevel() == kLow) Log() << gather.wrap(5) << '\n';
}

/* gather multiple values from all */
void CommunicatorT::AllGather(const nArrayT<int>& my, nArrayT<int>& gather) const
{
	const char caller[] = "CommunicatorT::AllGather";
	Log(kModerate, caller, "sending %d", my.Length());
	if (LogLevel() == kLow) Log() << my.wrap(5) << '\n';

	/* check */
	if (my.Length()*Size() != gather.Length()) 
		ExceptionT::SizeMismatch(caller, 
			"buffer length %d does not match size %d", gather.Length(), my.Length()*Size());

#ifdef __TAHOE_MPI__
	int len = my.Length();
	int ret = MPI_Allgather((void*) my.Pointer(), len, MPI_INT, gather.Pointer(), len, MPI_INT, fComm);

#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Allgather failed");
#endif

#else
	/* write my into gather */
	int len = my.Length();
	gather.CopyPart(Rank()*len, my, 0, len);
#endif

//NOTE: need to implement gather with sends and receives for CPLANT and other
//      platforms where Allgather seems to be troublesome.

	Log(kModerate, caller, "gathered %d", gather.Length());
	if (LogLevel() == kLow) Log() << gather.wrap(5) << '\n';
}

/* gather multiple values from all */
void CommunicatorT::AllGather(const nArrayT<double>& my, nArrayT<double>& gather,
	const nArrayT<int>& counts, const nArrayT<int>& displacements) const
{
	const char caller[] = "CommunicatorT::Gather";
	Log(kModerate, caller, "sending %d", my.Length());
	if (LogLevel() == kLow)
	{
		Log() << "counts = \n" << counts.wrap(5) << '\n';
		Log() << "displacements = \n" << displacements.wrap(5) << '\n';
		Log() << "sending = \n" << my.wrap(5) << '\n';		
	}

#if __option(extended_errorcheck)
	if (counts[Rank()] != my.Length()) ExceptionT::SizeMismatch(caller);
	if (counts.Length() != displacements.Length()) ExceptionT::SizeMismatch(caller);
	if (displacements.Last() + counts.Last() > gather.Length()) ExceptionT::SizeMismatch(caller);
		/* assume counts and displacements are OK interms of overlap and assume the
		 * displacements are monotonically increasing */
#endif

#ifdef __TAHOE_MPI__
	int len = my.Length();
	int ret = MPI_Allgatherv((void*) my.Pointer(), len, MPI_DOUBLE, gather.Pointer(), 
				(int*) counts.Pointer(), (int*) displacements.Pointer(), MPI_DOUBLE, fComm);

#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Allgatherv failed");
#endif

#else
	/* write my into gather */
	int len = my.Length();
	int offset = displacements[Rank()];
	gather.CopyPart(offset, my, 0, len);
#endif

//NOTE: need to implement gather with sends and receives for CPLANT and other
//      platforms where Allgather seems to be troublesome.

	Log(kModerate, caller, "gathered %d", gather.Length());
	if (LogLevel() == kLow) Log() << gather.wrap(5) << '\n';
}

/* gather multiple values from all */
void CommunicatorT::AllGather(const nArrayT<int>& my, nArrayT<int>& gather,
	const nArrayT<int>& counts, const nArrayT<int>& displacements) const
{
	const char caller[] = "CommunicatorT::Gather";
	Log(kModerate, caller, "sending %d", my.Length());
	if (LogLevel() == kLow)
	{
		Log() << "counts = \n" << counts.wrap(5) << '\n';
		Log() << "displacements = \n" << displacements.wrap(5) << '\n';
		Log() << "sending = \n" << my.wrap(5) << '\n';		
	}

#if __option(extended_errorcheck)
	if (counts[Rank()] != my.Length()) ExceptionT::SizeMismatch(caller);
	if (counts.Length() != displacements.Length()) ExceptionT::SizeMismatch(caller);
	if (displacements.Last() + counts.Last() > gather.Length()) ExceptionT::SizeMismatch(caller);
		/* assume counts and displacements are OK interms of overlap and assume the
		 * displacements are monotonically increasing */
#endif

#ifdef __TAHOE_MPI__
	int len = my.Length();
	int ret = MPI_Allgatherv((void*) my.Pointer(), len, MPI_INT, gather.Pointer(), 
				(int*) counts.Pointer(), (int*) displacements.Pointer(), MPI_INT, fComm);

#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Allgatherv failed");
#endif

#else
	/* write my into gather */
	int len = my.Length();
	int offset = displacements[Rank()];
	gather.CopyPart(offset, my, 0, len);
#endif

//NOTE: need to implement gather with sends and receives for CPLANT and other
//      platforms where Allgather seems to be troublesome.

	Log(kModerate, caller, "gathered %d", gather.Length());
	if (LogLevel() == kLow) Log() << gather.wrap(5) << '\n';
}

/* broadcast character array */
void CommunicatorT::Broadcast(int source, ArrayT<char>& data)
{
	if (source == Rank())
	{
		Log(kModerate, "Broadcast", "sending %d", data.Length());
		if (LogLevel() == kLow) Log() << setw(10) << "data: " << data.Pointer() << '\n';
	}

#ifdef __TAHOE_MPI__
	int ret = MPI_Bcast(data.Pointer(), data.Length(), MPI_CHAR, source, fComm);

#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::Broadcast", "MPI_Bcast failed");
#endif

#endif

	if (source != Rank())
	{
		Log(kModerate, "Broadcast", "received %d", data.Length());
		if (LogLevel() == kLow) Log() << setw(10) << "data: " << data.Pointer() << '\n';
	}
}

/* free any uncompleted requests */
void CommunicatorT::FreeRequests(ArrayT<MPI_Request>& requests) const
{
#ifdef __TAHOE_MPI__
	const char caller[] = "CommunicatorT::FreeRequests";

	/* free any uncompleted receive requests */
	for (int i = 0; i < requests.Length(); i++)
		if (requests[i] != MPI_REQUEST_NULL)
		{
			Log(kLow, "caller", "cancelling request %d/%d", i+1, requests.Length());
		
			/* cancel request */
			MPI_Cancel(&requests[i]);
			MPI_Status status;
			MPI_Wait(&requests[i], &status);
			int flag;
			MPI_Test_cancelled(&status, &flag);
			if (flag )
				Log(kLow, "caller", "cancelling request %d/%d: DONE", i+1, requests.Length());
			else	
				Log(kLow, "caller", "cancelling request %d/%d: FAIL", i+1, requests.Length());
		}
#else
#pragma unused(requests)
#endif
}

/* write status information */
void CommunicatorT::WriteStatus(ostream& out, const char* caller, 
	const MPI_Status& status) const
{
#ifdef __TAHOE_MPI__
	out << "\n " << caller << ": MPI status returned with error\n"
        <<   "     status.MPI_SOURCE: " << status.MPI_SOURCE << '\n'
        <<   "        status.MPI_TAG: " << status.MPI_TAG << '\n'
        <<   "      status.MPI_ERROR: " << status.MPI_ERROR << endl;
#else
#pragma unused(out)
#pragma unused(caller)
#pragma unused(status)
#endif
}

/* post non-blocking send */
void CommunicatorT::PostSend(const nArrayT<double>& data, int destination, int tag, 
	MPI_Request& request) const
{
	const char caller[] = "CommunicatorT::PostSend";

	Log(kModerate, caller, "posting send of %d to %d with tag %d", 
		data.Length(), destination, tag);
	if (LogLevel() == kLow)
		Log() << setw(10) << "data:\n" << data.wrap(5) << '\n';

	if (destination < 0 || destination >= Size())
		Log(kFail, caller, "destination ! (0 <= %d <= %d)", destination, Size());

#ifdef __TAHOE_MPI__
	/* post send */
	int ret = MPI_Isend((void*) data.Pointer(), data.Length(), MPI_DOUBLE, destination, tag, fComm, &request);

#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Isend failed");
#endif

#else
#pragma unused(request)
#endif
}

void CommunicatorT::PostSend(const nArrayT<int>& data, int destination, int tag, 
	MPI_Request& request) const
{
	const char caller[] = "CommunicatorT::PostSend";

	Log(kModerate, caller, "posting send of %d to %d with tag %d", 
		data.Length(), destination, tag);
	if (LogLevel() == kLow)
		Log() << setw(10) << "data:\n" << data.wrap(5) << '\n';

	if (destination < 0 || destination >= Size())
		Log(kFail, caller, "destination ! (0 <= %d <= %d)", destination, Size());

#ifdef __TAHOE_MPI__
	/* post send */
	int ret = MPI_Isend((void*) data.Pointer(), data.Length(), MPI_INT, destination, tag, fComm, &request);
	
#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Isend failed");
#endif

#else
#pragma unused(request)
#endif
}

/* post non-blocking send */
void CommunicatorT::PostReceive(nArrayT<double>& data, int source, 
	int tag, MPI_Request& request) const
{
	const char caller[] = "CommunicatorT::PostReceive";

	Log(kModerate, caller, "posting receive of %d from %d with tag %d", 
		data.Length(), source, tag);

	if (source < 0 || source >= Size())
		Log(kFail, caller, "source ! (0 <= %d <= %d)", source, Size());

#ifdef __TAHOE_MPI__
	/* post receive */
	int ret = MPI_Irecv(data.Pointer(), data.Length(), MPI_DOUBLE, source, tag, fComm, &request);

#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Irecv failed");
#endif

#else
#pragma unused(request)
#endif
}

void CommunicatorT::PostReceive(nArrayT<int>& data, int source, 
	int tag, MPI_Request& request) const
{
	const char caller[] = "CommunicatorT::PostReceive";

	Log(kModerate, caller, "posting receive of %d from %d with tag %d", 
		data.Length(), source, tag);

	if (source < 0 || source >= Size())
		Log(kFail, caller, "source ! (0 <= %d <= %d)", source, Size());

#ifdef __TAHOE_MPI__
	/* post receive */
	int ret = MPI_Irecv(data.Pointer(), data.Length(), MPI_INT, source, tag, fComm, &request);

#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Irecv failed");
#endif

#else
#pragma unused(request)
#endif
}

/* return the index the next receive */
void CommunicatorT::WaitReceive(const ArrayT<MPI_Request>& requests, int& index, int& source) const
{
	const char caller[] = "CommunicatorT::WaitReceive";
	Log(kModerate, caller, "waiting for 1 of %d", requests.Length());

	index = source = -1;

#ifdef __TAHOE_MPI__
	/* grab completed receive */
	MPI_Status status;
	int ret = MPI_Waitany(requests.Length(), (MPI_Request*) requests.Pointer(), &index, &status);

#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Waitany failed");
#endif

#ifdef CHECK_MPI_STATUS
	if (status.MPI_ERROR != MPI_SUCCESS) Log(kFail, caller, "bad status: %d", status.MPI_ERROR);
#endif
	
	source = status.MPI_SOURCE;
#endif

	Log(kModerate, caller, "received request at index %d from %d", index, source);
}

/* block until all sends posted with CommunicatorT::PostSend have completed */
void CommunicatorT::WaitSends(const ArrayT<MPI_Request>& requests)
{
	const char caller[] = "CommunicatorT::WaitSends";
	Log(kModerate, caller, "waiting for 1 of %d", requests.Length());

	/* complete all sends */
	for (int i = 0; i < requests.Length(); i++)
	{
		int index = -1;

#ifdef __TAHOE_MPI__
		/* grab completed receive */
		MPI_Status status;
		int ret = MPI_Waitany(requests.Length(), (MPI_Request*) requests.Pointer(), &index, &status);

#ifdef CHECK_MPI_RETURN
		if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Waitany failed");
#endif

#ifdef CHECK_MPI_STATUS
		if (status.MPI_ERROR != MPI_SUCCESS) {
			WriteStatus(Log(), caller, status);
			Log(kFail, caller, "bad status: %d", status.MPI_ERROR);
		}
#endif

#endif

		Log(kLow, caller, "completing send at index %d", index);
	}
}

/* post blocking send */
void CommunicatorT::Send(const nArrayT<double>& data, int destination, int tag) const
{
	const char caller[] = "CommunicatorT::Send";

	Log(kModerate, caller, "posting send of %d to %d with tag %d", 
		data.Length(), destination, tag);
	if (LogLevel() == kLow)
		Log() << setw(10) << "data:\n" << data.wrap(5) << '\n';

	if (destination < 0 || destination >= Size())
		Log(kFail, caller, "destination ! (0 <= %d <= %d)", destination, Size());

#ifdef __TAHOE_MPI__
	/* post send */
	int ret = MPI_Send((void*) data.Pointer(), data.Length(), MPI_DOUBLE, destination, tag, fComm);

#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Send failed");
#endif

#endif
}

void CommunicatorT::Send(const nArrayT<int>& data, int destination, int tag) const
{
	const char caller[] = "CommunicatorT::Send";

	Log(kModerate, caller, "posting send of %d to %d with tag %d", 
		data.Length(), destination, tag);
	if (LogLevel() == kLow)
		Log() << setw(10) << "data:\n" << data.wrap(5) << '\n';

	if (destination < 0 || destination >= Size())
		Log(kFail, caller, "destination ! (0 <= %d <= %d)", destination, Size());

#ifdef __TAHOE_MPI__
	/* post send */
	int ret = MPI_Send((void*) data.Pointer(), data.Length(), MPI_INT, destination, tag, fComm);

#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Send failed");
#endif

#endif
}

/* post blocking receive */
void CommunicatorT::Receive(nArrayT<double>& data, int source, int tag) const
{
	const char caller[] = "CommunicatorT::Receive";

	Log(kModerate, caller, "posting receive for %d from %d with tag %d", 
		data.Length(), source, tag);

	if (source < 0 || source >= Size())
		Log(kFail, caller, "source ! (0 <= %d <= %d)", source, Size());

#ifdef __TAHOE_MPI__
	/* post receive */
	MPI_Status status;
	int ret = MPI_Recv(data.Pointer(), data.Length(), MPI_DOUBLE, source, tag, fComm, &status);

#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Recv failed");
#endif

#ifdef CHECK_MPI_STATUS
	if (status.MPI_ERROR != MPI_SUCCESS) {
		WriteStatus(Log(), caller, status);
		Log(kFail, caller, "bad status: %d", status.MPI_ERROR);
	}
#endif

#endif

	if (LogLevel() == kLow)
		Log() << setw(10) << "data:\n" << data.wrap(5) << '\n';
}

void CommunicatorT::Receive(nArrayT<int>& data, int source, int tag) const
{
	const char caller[] = "CommunicatorT::Receive";

	Log(kModerate, caller, "posting receive for %d from %d with tag %d", 
		data.Length(), source, tag);

	if (source < 0 || source >= Size())
		Log(kFail, caller, "source ! (0 <= %d <= %d)", source, Size());

#ifdef __TAHOE_MPI__
	/* post receive */
	MPI_Status status;
	int ret = MPI_Recv(data.Pointer(), data.Length(), MPI_INT, source, tag, fComm, &status);

#ifdef CHECK_MPI_RETURN
	if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Recv failed");
#endif

#ifdef CHECK_MPI_STATUS
	if (status.MPI_ERROR != MPI_SUCCESS) {
		WriteStatus(Log(), caller, status);
		Log(kFail, caller, "bad status: %d", status.MPI_ERROR);
	}
#endif

#endif

	if (LogLevel() == kLow)
		Log() << setw(10) << "data:\n" << data.wrap(5) << '\n';
}

/*************************************************************************
 * Private
 *************************************************************************/

void CommunicatorT::Init(void)
{
	const char caller[] = "CommunicatorT::Init";

	/* environment was shut down */
	if (fCount == -1)
		ExceptionT::MPIFail(caller, "cannot restart MPI environment");

	/* communicator count */
	fCount++;

#ifdef __TAHOE_MPI__
	if (fCount == 1)
	{
		int *argc_ = fargc, argc = 1;
		char ***argv_ = fargv, **argv = NULL;

		/* need dummy arguments */
		if (!fargv) {
			argc_ = &argc;
			argv_ = &argv;
			argv = new char*[1];
			argv[0] = new char[2];
			argv[0][0] = 'a';
			argv[0][1] = '\0';
		}

		/* initialize MPI environment */
		int ret = MPI_Init(argc_, argv_);
#ifdef CHECK_MPI_RETURN
		if (ret != MPI_SUCCESS) Log(kFail, caller, "MPI_Init failed");
#endif
		
		/* free */
		if (!fargv) {
			for (int i = 0; i < argc; i++)
				delete[] argv[i];
			delete[] argv;
		}
	}
#endif	

	Log(kModerate, caller, "communicator count = %d", fCount);
}

void CommunicatorT::Finalize(void)
{
	/* environment was shut down */
	if (fCount == -1)
		ExceptionT::MPIFail("CommunicatorT::Finalize", "MPI environment already down");

	/* communicator count */
	fCount--;

#ifdef __TAHOE_MPI__
	if (fCount == 0) /* last communicator */
	{
		/* shut down MPI environment */
		int ret = MPI_Finalize();
#ifdef CHECK_MPI_RETURN
		if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::Finalize", "MPI_Finalized failed");
#endif
	}

	/* free communicator */
	if (fComm != MPI_COMM_WORLD && fComm != MPI_COMM_NULL)
	{
		int ret = MPI_Comm_free(&fComm);
#ifdef CHECK_MPI_RETURN
		if (ret != MPI_SUCCESS) Log(kFail, "CommunicatorT::Finalize", "MPI_Comm_free failed");
#endif
	}
#endif

	/* close */
	Log(kModerate, "Finalize", "communicator count = %d", fCount);
	if (fCount == 0) {
		fCount = -1;
		Log(kModerate, "Finalize", "communicator count = %d", fCount);
	}
}

void CommunicatorT::doLog(const char* caller, const char* message) const
{
	ostream& out = *fLog;

	/* current clock */
	time_t t;
	time(&t);

	/* cast away const-ness */
	CommunicatorT* non_const_this = (CommunicatorT*) this;

	/* elapsed time */
#ifdef __TAHOE_MPI__
	double new_time = MPI_Wtime();
	double elapsed_time = new_time - fLastTime;
	non_const_this->fLastTime = new_time;
#else
	double elapsed_time = double(t - fLastTime)/CLOCKS_PER_SEC;
	non_const_this->fLastTime = t;
#endif	    	    

	/* write info */
	out << '\n' << setw(10) << "CommunicatorT::doLog: rank: " << Rank() << '\n';
	if (caller) out << setw(10) << " caller: " << caller << '\n';		
	out << setw(10) << "  clock: " << ctime(&t);
	out << setw(10) << "elapsed: " << elapsed_time << " sec" << '\n';
	if (message) out << setw(10) << "message: " << message << '\n';

	/* flush stream */
	if (fLogLevel != kLow) out.flush();
}
