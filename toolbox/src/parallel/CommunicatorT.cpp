/* $Id: CommunicatorT.cpp,v 1.7 2002-12-05 08:25:19 paklein Exp $ */
#include "CommunicatorT.h"
#include "ExceptionT.h"
#include <iostream.h>
#include <time.h>
#include "ArrayT.h"

using namespace Tahoe;

/* initialize static variables */
int CommunicatorT::fCount = 0;
int* CommunicatorT::fargc = NULL;
char*** CommunicatorT::fargv = NULL;

/* create communicator including all processes */
CommunicatorT::CommunicatorT(void):
	fLogLevel(kSilent),
	fLog(&cout)
{
	/* check MPI environment */
	Init();

#ifdef __TAHOE_MPI__
	fComm = MPI_COMM_WORLD;
	if (MPI_Comm_size(fComm, &fSize) != MPI_SUCCESS) {
		Log("CommunicatorT::CommunicatorT", "MPI_Comm_size failed", true);
		throw ExceptionT::kMPIFail;
	}
	if (MPI_Comm_rank(fComm, &fRank) != MPI_SUCCESS) {
		Log("CommunicatorT::CommunicatorT", "MPI_Comm_rank failed", true);
		throw ExceptionT::kMPIFail;
	}
#else
	fComm = 0;
	fSize = 1;
	fRank = 0;
#endif
}

/* copy constructor */
CommunicatorT::CommunicatorT(const CommunicatorT& source):
	fComm(source.fComm),
	fSize(source.fSize),
	fRank(source.fRank),
	fLogLevel(source.fLogLevel),
	fLog(source.fLog)
{
	/* check MPI environment */
	Init();
}

/* destructor */
CommunicatorT::~CommunicatorT(void)
{
	/* check MPI environment */
	Finalize();
}

/* write message to log */
void CommunicatorT::Log(const char* caller, const char* message, bool force) const
{
	if (fLogLevel != kSilent || force)
		LogHead(caller) << message << endl;
}

/* (re-)set the logging stream */
void CommunicatorT::SetLog(ostream& log)
{
	Log("CommunicatorT::SetLog", "closing log stream");
	fLog = &log;
	Log("CommunicatorT::SetLog", "opening log stream");
}

/* maximum over single integers */
int CommunicatorT::Max(int a) const
{
	if (fLogLevel != kSilent)
		*fLog << '\n' << Rank() << ": CommunicatorT::Max: in = " << a << endl;

#ifdef __TAHOE_MPI__
	if (Size() > 1)
	{
		int b = a;
		if (MPI_Allreduce(&b, &a, 1, MPI_INT, MPI_MAX, fComm) != MPI_SUCCESS) {
			Log("CommunicatorT::Max", "MPI_Allreduce failed", true);
			throw ExceptionT::kMPIFail;
		}
	}
#endif

	if (fLogLevel != kSilent)
		*fLog << Rank() << ": CommunicatorT::Max: out = " << a << endl;

	return a;
}

/* minimum over single integers */
int CommunicatorT::Min(int a) const
{
	if (fLogLevel != kSilent)
		*fLog << '\n' << Rank() << ": CommunicatorT::Min: in = " << a << endl;

#ifdef __TAHOE_MPI__
	if (Size() > 1)
	{
		int b = a;
		if (MPI_Allreduce(&b, &a, 1, MPI_INT, MPI_MIN, fComm) != MPI_SUCCESS) {
			Log("CommunicatorT::Min", "MPI_Allreduce failed", true);
			throw ExceptionT::kMPIFail;
		}
	}
#endif

	if (fLogLevel != kSilent)
		*fLog << Rank() << ": CommunicatorT::Min: out = " << a << endl;

	return a;
}

/* minimum over single integers */
int CommunicatorT::Sum(int a) const
{
	if (fLogLevel != kSilent)
		*fLog << '\n' << Rank() << ": CommunicatorT::Sum: in = " << a << endl;

#ifdef __TAHOE_MPI__
	if (Size() > 1)
	{
		int b = a;
		if (MPI_Allreduce(&b, &a, 1, MPI_INT, MPI_SUM, fComm) != MPI_SUCCESS) {
			Log("CommunicatorT::Sum", "MPI_Allreduce failed", true);
			throw ExceptionT::kMPIFail;
		}
	}
#endif

	if (fLogLevel != kSilent)
		*fLog << Rank() << ": CommunicatorT::Sum: out = " << a << endl;

	return a;
}

/* maximum over single doubles */
double CommunicatorT::Max(double a) const
{
	if (fLogLevel != kSilent)
		*fLog << '\n' << Rank() << ": CommunicatorT::Max: in = " << a << endl;

#ifdef __TAHOE_MPI__
	if (Size() > 1)
	{
		double b = a;
		if (MPI_Allreduce(&b, &a, 1, MPI_DOUBLE, MPI_MAX, fComm) != MPI_SUCCESS) {
			Log("CommunicatorT::Max", "MPI_Allreduce failed", true);
			throw ExceptionT::kMPIFail;
		}
	}
#endif
	if (fLogLevel != kSilent)
		*fLog << Rank() << ": CommunicatorT::Max: out = " << a << endl;

	return a;
}

/* minimum over single doubles */
double CommunicatorT::Min(double a) const
{
	if (fLogLevel != kSilent)
		*fLog << '\n' << Rank() << ": CommunicatorT::Min: in = " << a << endl;

#ifdef __TAHOE_MPI__
	if (Size() > 1)
	{
		double b = a;
		if (MPI_Allreduce(&b, &a, 1, MPI_DOUBLE, MPI_MIN, fComm) != MPI_SUCCESS) {
			Log("CommunicatorT::Min", "MPI_Allreduce failed", true);
			throw ExceptionT::kMPIFail;
		}
	}
#endif

	if (fLogLevel != kSilent)
		*fLog << Rank() << ": CommunicatorT::Min: out = " << a << endl;

	return a;
}

/* minimum over single integers */
double CommunicatorT::Sum(double a) const
{
	if (fLogLevel != kSilent)
		*fLog << '\n' << Rank() << ": CommunicatorT::Sum: in = " << a << endl;

#ifdef __TAHOE_MPI__
	if (Size() > 1)
	{
		double b = a;
		if (MPI_Allreduce(&b, &a, 1, MPI_DOUBLE, MPI_SUM, fComm) != MPI_SUCCESS) {
			Log("CommunicatorT::Sum", "MPI_Allreduce failed", true);
			throw ExceptionT::kMPIFail;
		}
	}
#endif

	if (fLogLevel != kSilent)
		*fLog << Rank() << ": CommunicatorT::Sum: out = " << a << endl;

	return a;
}

/* gather single integer. Called by destination process. */
void CommunicatorT::Gather(int a, ArrayT<int>& gather) const
{
	/* check */
	if (gather.Length() != Size()) ExceptionT::SizeMismatch("CommunicatorT::Gather");

	/* this */
	gather[Rank()] = a;

#ifdef __TAHOE_MPI__
	if (MPI_Gather(&a, 1, MPI_INT, gather.Pointer(), 1, MPI_INT, Rank(), fComm) != MPI_SUCCESS) {
		Log("CommunicatorT::Gather", "MPI_Gather failed", true);
		throw ExceptionT::kMPIFail;
	}
#endif
}

/* gather single integer. Called by sending processes. */
void CommunicatorT::Gather(int a, int destination) const
{
	/* check */
	if (destination == Rank()) {
		cout << "\n CommunicatorT::Gather: destination must be different from rank: " 
		     << destination << endl;
		throw ExceptionT::kGeneralFail;
	}

#ifdef __TAHOE_MPI__
	int* tmp = NULL;
	if (MPI_Gather(&a, 1, MPI_INT, tmp, 1, MPI_INT, destination, fComm) != MPI_SUCCESS) {
		Log("CommunicatorT::Gather", "MPI_Gather failed", true);
		throw ExceptionT::kMPIFail;
	}
#else
#pragma unused(a)
#endif
}

/* synchronize all processes */
void CommunicatorT::Barrier(void) const
{
#ifdef __TAHOE_MPI__
	if (MPI_Barrier(fComm) != MPI_SUCCESS) {
		Log("CommunicatorT::Barrier", "MPI_Barrier failed", true);
		throw ExceptionT::kMPIFail;
	}
#endif
}

/* gather single integer to all processes. */
void CommunicatorT::AllGather(int a, ArrayT<int>& gather) const
{
	/* check */
	if (gather.Length() != Size()) ExceptionT::SizeMismatch("CommunicatorT::AllGather");

	/* this */
	gather[Rank()] = a;

#ifdef __TAHOE_MPI__
	if (MPI_Allgather(&a, 1, MPI_INT, gather.Pointer(), 1, MPI_INT, fComm) != MPI_SUCCESS) {
		Log("CommunicatorT::Gather", "MPI_Allgather failed", true);
		throw ExceptionT::kMPIFail;
	}
#endif

//NOTE: need to implement gather with sends and receives for CPLANT and other
//      platforms where Allgather seems to be troublesome.
}

/* gather multiple values from all */
void CommunicatorT::AllGather(const ArrayT<double>& my, ArrayT<double>& gather) const
{
	/* check */
	if (my.Length()*Size() != gather.Length()) ExceptionT::SizeMismatch("CommunicatorT::AllGather");

#ifdef __TAHOE_MPI__
	int len = my.Length();
	if (MPI_Allgather(my.Pointer(), len, MPI_DOUBLE, gather.Pointer(), len, MPI_DOUBLE, fComm) != MPI_SUCCESS) {
		Log("CommunicatorT::Gather", "MPI_Allgather failed", true);
		throw ExceptionT::kMPIFail;
	}
#else
	/* write my into gather */
	int len = my.Length();
	gather.CopyPart(Rank()*len, my, 0, len);
#endif

//NOTE: need to implement gather with sends and receives for CPLANT and other
//      platforms where Allgather seems to be troublesome.
}

/* gather multiple values from all */
void CommunicatorT::AllGather(const ArrayT<int>& my, ArrayT<int>& gather) const
{
	/* check */
	if (my.Length()*Size() != gather.Length()) ExceptionT::SizeMismatch("CommunicatorT::AllGather");

#ifdef __TAHOE_MPI__
	int len = my.Length();
	if (MPI_Allgather(my.Pointer(), len, MPI_INT, gather.Pointer(), len, MPI_INT, fComm) != MPI_SUCCESS) {
		Log("CommunicatorT::Gather", "MPI_Allgather failed", true);
		throw ExceptionT::kMPIFail;
	}
#else
	/* write my into gather */
	int len = my.Length();
	gather.CopyPart(Rank()*len, my, 0, len);
#endif

//NOTE: need to implement gather with sends and receives for CPLANT and other
//      platforms where Allgather seems to be troublesome.
}

/* gather multiple values from all */
void CommunicatorT::AllGather(const ArrayT<int>& counts, const ArrayT<int>& displacements, 
	const ArrayT<double>& my, ArrayT<double>& gather) const
{
#if __option(extended_errorcheck)
	if (counts[Rank()] != my.Length()) ExceptionT::SizeMismatch();
	if (counts.Length() != displacements.Length()) ExceptionT::SizeMismatch();
	if (gather.Length() > 1 && 
		(displacements.Last() + counts.Last() >= gather.Length())) ExceptionT::SizeMismatch();
		/* assume counts and displacements are OK interms of overlap and assume the
		 * displacements are monotonically increasing */
#endif

#ifdef __TAHOE_MPI__
	int len = my.Length();
	if (MPI_Allgatherv(my.Pointer(), len, MPI_DOUBLE, 
		gather.Pointer(), counts.Pointer(), displacements.Pointer(), MPI_DOUBLE, fComm) != MPI_SUCCESS) {
		Log("CommunicatorT::Gather", "MPI_Allgatherv failed", true);
		throw ExceptionT::kMPIFail;
	}
#else
	/* write my into gather */
	int len = my.Length();
	int offset = displacements[Rank()];
	gather.CopyPart(offset, my, 0, len);
#endif

//NOTE: need to implement gather with sends and receives for CPLANT and other
//      platforms where Allgather seems to be troublesome.
}

/* gather multiple values from all */
void CommunicatorT::AllGather(const ArrayT<int>& counts, const ArrayT<int>& displacements, 
	const ArrayT<int>& my, ArrayT<int>& gather) const
{
#if __option(extended_errorcheck)
	if (counts[Rank()] != my.Length()) ExceptionT::SizeMismatch();
	if (counts.Length() != displacements.Length()) ExceptionT::SizeMismatch();
	if (gather.Length() > 1 && 
		(displacements.Last() + counts.Last() >= gather.Length())) ExceptionT::SizeMismatch();
		/* assume counts and displacements are OK interms of overlap and assume the
		 * displacements are monotonically increasing */
#endif

#ifdef __TAHOE_MPI__
	int len = my.Length();
	if (MPI_Allgatherv(my.Pointer(), len, MPI_INT, 
		gather.Pointer(), counts.Pointer(), displacements.Pointer(), MPI_INT, fComm) != MPI_SUCCESS) {
		Log("CommunicatorT::Gather", "MPI_Allgatherv failed", true);
		throw ExceptionT::kMPIFail;
	}
#else
	/* write my into gather */
	int len = my.Length();
	int offset = displacements[Rank()];
	gather.CopyPart(offset, my, 0, len);
#endif

//NOTE: need to implement gather with sends and receives for CPLANT and other
//      platforms where Allgather seems to be troublesome.
}

/* broadcast character array */
void CommunicatorT::Broadcast(ArrayT<char>& data)
{
#ifdef __TAHOE_MPI__
	if (MPI_Bcast(data.Pointer(), data.Length(), MPI_CHAR, 0, fComm) != MPI_SUCCESS) {
		Log("CommunicatorT::Broadcast", "MPI_Bcast failed", true);
		throw ExceptionT::kMPIFail;
	}
#else
#pragma unused(data)
#endif
}

/*************************************************************************
* Protected
*************************************************************************/

/* returns the time string */
const char* CommunicatorT::WallTime(void) const
{
	time_t t;
	time(&t);
	return ctime(&t);
}

/* write log header */
ostream& CommunicatorT::LogHead(const char* caller) const
{
	*fLog << Rank() << ": "
	     << WallTime() << '\t'
	     << caller << ": ";
	return *fLog;
}

/*************************************************************************
* Private
*************************************************************************/

void CommunicatorT::Init(void)
{
	/* environment was shut down */
	if (fCount == -1) {
		cout << "\n CommunicatorT::Init: cannot restart MPI environment" << endl;
		throw ExceptionT::kMPIFail;
	}

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
		if (MPI_Init(argc_, argv_) != MPI_SUCCESS) {
			Log("CommunicatorT::Init", "MPI_Init failed", true);
			throw ExceptionT::kMPIFail;
		}
		
		/* free */
		if (!fargv) {
			for (int i = 0; i < argc; i++)
				delete[] argv[i];
			delete[] argv;
		}
	}
#endif	
}

void CommunicatorT::Finalize(void)
{
	/* environment was shut down */
	if (fCount == -1) {
		cout << "\n CommunicatorT::Finalize: MPI environment already down" << endl;
		throw ExceptionT::kMPIFail;
	}

	/* communicator count */
	fCount--;

#ifdef __TAHOE_MPI__
	if (fCount == 0)
	{
		/* shut down MPI environment */
		if (MPI_Finalize() != MPI_SUCCESS) {
			Log("CommunicatorT::Finalize", "MPI_Finalized failed", true);
			throw ExceptionT::kMPIFail;
		}
	}
#endif

	/* close */
	if (fCount == 0) fCount = -1;
}
