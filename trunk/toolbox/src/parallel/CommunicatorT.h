/* $Id: CommunicatorT.h,v 1.12 2004-03-16 05:35:51 paklein Exp $ */
#ifndef _COMMUNICATOR_T_H_
#define _COMMUNICATOR_T_H_

#include "Environment.h"
#include <time.h>

#ifdef __TAHOE_MPI__
#include "mpi.h"
#else
#define MPI_COMM_NULL  0
#define MPI_UNDEFINED -1
#define MPI_Comm    long
#define MPI_Request long
#define MPI_Status  long
#define MPI_Op      long
#endif

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
template <class TYPE> class nArrayT;
template <class TYPE> class ArrayT;
class ofstreamT;

/** interface to handle process to process communication. If compiled
 * with -D__TAHOE_MPI__, will implemented process-to-process communication
 * using MPI. Otherwise, will assume single process. */
class CommunicatorT
{
  public:

	/** class constants */
	enum ConstantsT {
		kNoColor = MPI_UNDEFINED
	};

	/** \name constructor */
	/*@{*/
	/** create communicator including all processes */
	CommunicatorT(void); 
	
	/** copy constructor */
	CommunicatorT(const CommunicatorT& source); 

	/** split communicator.
	 * Construct a communicator as subset of the source based on color. Each
	 * process returns a CommunicatorT which includes the other processes that
	 * have called the constructor with the same color.
	 * \param source communicator to split. All processes in source must be
	 *        synchronized in calling this constructor.
	 * \param color nonnegative key to determine subset of processes belong to this 
	 *        group. Pass CommunicatorT::kNoColor to be left out of any new groups.
	 * \param rank_key used to determine the rank in the new communicator by sorting
	 *        keys in ascending order, breaking ties in a consistent way */
	CommunicatorT(const CommunicatorT& source, int color, int rank_key = 0); 
	/*@}*/

	/** destructor */
	~CommunicatorT(void);
	
	/** (optionally) set argc and argv to pass command line arguments */
	static void SetArgv(int* argc, char*** argv) {
			fargc = argc;
			fargv = argv;
		};

	/** communicator size */
	int Size(void) const { return fSize; };

	/** rank of this process */
	int Rank(void) const { return fRank; };

	/** return the raw MPI_Comm */
	const MPI_Comm& Comm(void) const { return fComm; };

	/** type conversion operator */
	operator const MPI_Comm() const { return Comm(); };
	
	/** returns true if an MP environment is active */
	static bool ActiveMP(void);

	/** \name logging */
	/*@{*/
	/** logging level */
	enum LogLevelT {
           kLow = 0, /**< log everything including message data */ 
      kModerate = 1, /**< log destination/source */
        kUrgent = 2, /**< serious error messages only */
          kFail = 3  /**< log message and then throw an ExceptionT::kMPIFail */
           };

	/** write log message
	 * \param priority message suppressed if the priority doesn't 
	 *        match the current CommunicatorT::LogLevel
	 * \param caller the calling subroutine */
	void Log(LogLevelT priority, const char* caller) const;

	/** write log message
	 * \param priority message suppressed if the priority doesn't 
	 *        match the current CommunicatorT::LogLevel
	 * \param caller the calling subroutine
	 * \param fmt formatted message string */
	void Log(LogLevelT priority, const char* caller, const char* fmt, ...) const;

	/** logging stream */
	ostream& Log(void) const { return *fLog; };
	
	void SetLog(ostream& log);
	void SetLog(ofstreamT& log);
	
	/** logging level */
	LogLevelT LogLevel(void) const { return fLogLevel; };

	/** (re-)set logging level */
	void SetLogLevel(LogLevelT log_level) { fLogLevel = log_level; };

	/** return elapsed time */
	static double Time(void);
	/*@}*/
	
	/** maximum over single integers returned to all */
	int Max(int a) const;

	/** minimum over single integers returned to all */
	int Min(int a) const;

	/** sum of single integers returned to all */
	int Sum(int a) const;

	/** maximum over single doubles returned to all */
	double Max(double a) const;

	/** minimum over single doubles returned to all */
	double Min(double a) const;

	/** sum of single doubles returned to all */
	double Sum(double a) const;

	/** element-by-element sum of a vector */
	void Sum(const nArrayT<double>& my_values, nArrayT<double>& sum) const;

	/** \name gather from all to one
	 * Gather single integer from all processes to one. */
	/*@{*/	
	/** Called by destination process. 
	 * \param gather destination for gathered values. Must be dimensioned
	 *        to length CommunicatorT::Size before the call. */
	void Gather(int a, nArrayT<int>& gather) const;

	/** gather single integer. Called by sending processes. */
	void Gather(int a, int destination) const;
	/*@}*/

	/** \name gather an array of integers from all to one
	 * Gather an array of integers from all processes to one. */
	/*@{*/	
	/** Called by destination process. 
	 * \param gather destination for gathered values. Must be dimensioned
	 *        to length CommunicatorT::Size before the call. */
	void Gather(const nArrayT<int>& my_values, nArrayT<int>& gather, 
		const nArrayT<int>& counts, const nArrayT<int>& displacements) const;

	/** gather an array of integers. Called by sending processes. Number of
	 * values sent from this process must already been known by the receiving
	 * process to call the matching Gather method. */
	void Gather(const nArrayT<int>& my_values, int destination) const;
	/*@}*/

	/** gather single integer to all processes. 
	 * \param gather destination for gathered values. Must be dimensioned
	 *        to length CommunicatorT::Size before the call. */
	void AllGather(int a, nArrayT<int>& gather) const;

	/** \name gather multiple values from all. 
	 * All processes sending the same amount of information to all. */
	/*@{*/
	/** gather double's */
	void AllGather(const nArrayT<double>& my, nArrayT<double>& gather) const;

	/** gather int's */
	void AllGather(const nArrayT<int>& my, nArrayT<int>& gather) const;
	/*@}*/

	/** \name gather multiple values from all. 
	 * Arbitrary data size from each process 
	 * \param my data sent from this processor
	 * \param gather the destination for the gathered data
	 * \param counts amount of data sent from each process
	 * \param displacements offset in the destination array where the data
	 *        from the corresponding rank should be written */
	/*@{*/
	/** gather double's */
	void AllGather(const nArrayT<double>& my, nArrayT<double>& gather,
		const nArrayT<int>& counts, const nArrayT<int>& displacements) const;

	/** gather int's */
	void AllGather(const nArrayT<int>& my, nArrayT<int>& gather,
		const nArrayT<int>& counts, const nArrayT<int>& displacements) const;
	/*@}*/

	/** broadcast character array */
	void Broadcast(int source, ArrayT<char>& data);
	
	/** synchronize all processes */
	void Barrier(void) const;

	/** \name non-blocking, point-to-point communications */
	/*@{*/
	/** post non-blocking send for double data */
	void PostSend(const nArrayT<double>& data, int destination, int tag, MPI_Request& request) const;

	/** post non-blocking send for int data */
	void PostSend(const nArrayT<int>& data, int destination, int tag, MPI_Request& request) const;

	/** post non-blocking receive double data */
	void PostReceive(nArrayT<double>& data, int source, int tag, MPI_Request& request) const;

	/** post non-blocking receive int data */
	void PostReceive(nArrayT<int>& data, int source, int tag, MPI_Request& request) const;
	
	/** return the index the next receive. Call blocks until the next message is received.
	 * \param requests list of requests previously posted with CommunicatorT::PostReceive 
	 * \return index in requests of the next received message*/
	void WaitReceive(const ArrayT<MPI_Request>& requests, int& index, int& source) const;
	
	/** block until all sends posted with CommunicatorT::PostSend have completed */
	void WaitSends(const ArrayT<MPI_Request>& requests);
	/*@}*/

	/** write status information */
	void WriteStatus(ostream& out, const char* caller, const MPI_Status& status) const;

	/** free any uncompleted requests. This often doesn't return successfully
	 * because one or more of the other processes has quit */
	void FreeRequests(ArrayT<MPI_Request>& requests) const;

	/** \name blocking, point-to-point communications */
	/*@{*/
	/** post blocking send for an array of doubles */
	void Send(const nArrayT<double>& data, int destination, int tag) const;

	/** post blocking send for an array of ints */
	void Send(const nArrayT<int>& data, int destination, int tag) const;

	/** post blocking receive for an array of doubles */
	void Receive(nArrayT<double>& data, int source, int tag) const;

	/** post blocking receive for an array ints */
	void Receive(nArrayT<int>& data, int source, int tag) const;
	/*@}*/

  private:
  
  	/** \name initialize/finalize MPI environment */
  	/*@{*/
	void Init(void);
	void Finalize(void);
  	/*@}*/

	/** write log message */
	void doLog(const char* caller, const char* message) const;
  	
  private:
  
  	/** \name MPI parameters */
  	/*@{*/
  	MPI_Comm fComm; /**< MPI communicator */
  	int fSize;      /**< communicator size */
  	int fRank;      /**< rank of this process */
  	/*@}*/
  	
  	/** \name logging */
  	/*@{*/
  	LogLevelT fLogLevel;
  	ostream*  fLog;
  	/*@}*/
  	
	/** count of communicators */
	static int fCount;
	static int* fargc;
	static char*** fargv;

	/** running time stamp needed for calculating intervals
	 * between calls to CommunicatorT:doLog. Uses MPI timing
	 * if available */
#ifdef __TAHOE_MPI__
	double fLastTime;
#else
	clock_t fLastTime;
#endif
};

/* returns true if an MP environment is active */
inline bool CommunicatorT::ActiveMP(void)
{
#ifdef __TAHOE_MPI__
	return true;
#else
	return false;
#endif
}

} // namespace Tahoe 
#endif /* _COMMUNICATOR_T_H_ */
