/* $Id: CommunicatorT.cpp,v 1.1 2002-02-12 02:04:07 paklein Exp $ */

#include "CommunicatorT.h"
#include "ExceptionCodes.h"
#include <iostream.h>
#include <time.h>
#include "ArrayT.h"

/* create communicator including all processes */
CommunicatorT::CommunicatorT(void):
	fLogLevel(kSilent),
	fLog(&cout)
{
#ifdef __MPI__
	fComm = MPI_COMM_WORLD;
	if (MPI_Comm_size(fComm, &fSize) != MPI_SUCCESS) {
		Log("CommunicatorT::CommunicatorT", "MPI_Comm_size failed", true);
		throw eMPIFail;
	}
	if (MPI_Comm_rank(fComm, &fRank) != MPI_SUCCESS) {
		Log("CommunicatorT::CommunicatorT", "MPI_Comm_rank failed", true);
		throw eMPIFail;
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
#ifdef __MPI__
	if (Size() > 1)
	{
		int b = a;
		if (MPI_Allreduce(&b, &a, 1, MPI_INT, MPI_MAX, fComm) != MPI_SUCCESS) {
			Log("CommunicatorT::Max", "MPI_Allreduce failed", true);
			throw eMPIFail;
		}
	}
#endif
	return a;
}

/* minimum over single integers */
int CommunicatorT::Min(int a) const
{
#ifdef __MPI__
	if (Size() > 1)
	{
		int b = a;
		if (MPI_Allreduce(&b, &a, 1, MPI_INT, MPI_MIN, fComm) != MPI_SUCCESS) {
			Log("CommunicatorT::Min", "MPI_Allreduce failed", true);
			throw eMPIFail;
		}
	}
#endif
	return a;
}

/* minimum over single integers */
int CommunicatorT::Sum(int a) const
{
#ifdef __MPI__
	if (Size() > 1)
	{
		int b = a;
		if (MPI_Allreduce(&b, &a, 1, MPI_INT, MPI_SUM, fComm) != MPI_SUCCESS) {
			Log("CommunicatorT::Sum", "MPI_Allreduce failed", true);
			throw eMPIFail;
		}
	}
#endif
	return a;
}

/* maximum over single doubles */
double CommunicatorT::Max(double a) const
{
#ifdef __MPI__
	if (Size() > 1)
	{
		double b = a;
		if (MPI_Allreduce(&b, &a, 1, MPI_DOUBLE, MPI_MAX, fComm) != MPI_SUCCESS) {
			Log("CommunicatorT::Max", "MPI_Allreduce failed", true);
			throw eMPIFail;
		}
	}
#endif
	return a;
}

/* minimum over single doubles */
double CommunicatorT::Min(double a) const
{
#ifdef __MPI__
	if (Size() > 1)
	{
		double b = a;
		if (MPI_Allreduce(&b, &a, 1, MPI_DOUBLE, MPI_MIN, fComm) != MPI_SUCCESS) {
			Log("CommunicatorT::Min", "MPI_Allreduce failed", true);
			throw eMPIFail;
		}
	}
#endif
	return a;
}

/* minimum over single integers */
double CommunicatorT::Sum(double a) const
{
#ifdef __MPI__
	if (Size() > 1)
	{
		double b = a;
		if (MPI_Allreduce(&b, &a, 1, MPI_DOUBLE, MPI_SUM, fComm) != MPI_SUCCESS) {
			Log("CommunicatorT::Sum", "MPI_Allreduce failed", true);
			throw eMPIFail;
		}
	}
#endif
	return a;
}

/* gather single integer. Called by destination process. */
void CommunicatorT::Gather(int a, ArrayT<int>& gather) const
{
	/* allocate */
	gather.Allocate(Size());
	gather[Rank()] = a;

#ifdef __MPI__
	if (MPI_Gather(&a, 1, MPI_INT, gather.Pointer(), 1, MPI_INT, Rank(), fComm) != MPI_SUCCESS) {
		Log("CommunicatorT::Gather", "MPI_Gather failed", true);
		throw eMPIFail;
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
		throw eGeneralFail;
	}

#ifdef __MPI__
	int* tmp = NULL;
	if (MPI_Gather(&a, 1, MPI_INT, tmp, 1, MPI_INT, destination, fComm) != MPI_SUCCESS) {
		Log("CommunicatorT::Gather", "MPI_Gather failed", true);
		throw eMPIFail;
	}
#else
#pragma unused(a)
#endif
}

/* synchronize all processes */
void CommunicatorT::Barrier(void) const
{
#ifdef __MPI__
	if (MPI_Barrier(fComm) != MPI_SUCCESS) {
		Log("CommunicatorT::Barrier", "MPI_Barrier failed", true);
		throw eMPIFail;
	}
#endif
}

/* gather single integer to all processes. */
void CommunicatorT::AllGather(int a, ArrayT<int>& gather) const
{
	/* allocate */
	gather.Allocate(Size());
	gather[Rank()] = a;

#ifdef __MPI__
	if (MPI_Allgather(&a, 1, MPI_INT, gather.Pointer(), 1, MPI_INT, fComm) != MPI_SUCCESS) {
		Log("CommunicatorT::Gather", "MPI_Allgather failed", true);
		throw eMPIFail;
	}
#endif

//NOTE: need to implement gather with sends and receives for CPLANT and other
//      platforms where Allgather seems to be troublesome.
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
	     << WallTime() << ": "
	     << caller << ": ";
	return *fLog;
}
