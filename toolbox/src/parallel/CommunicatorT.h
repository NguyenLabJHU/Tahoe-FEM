/* $Id: CommunicatorT.h,v 1.8.2.1 2002-12-10 17:03:41 paklein Exp $ */
#ifndef _COMMUNICATOR_T_H_
#define _COMMUNICATOR_T_H_

#include "Environment.h"

#ifdef __TAHOE_MPI__
#include "mpi.h"
#else
#define MPI_Comm long
#define MPI_Request long
#endif

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
template <class TYPE> class nArrayT;
template <class TYPE> class ArrayT;

/** interface to handle process to process communication. If compiled
 * with -D__TAHOE_MPI__, will implemented process-to-process communication
 * using MPI. Otherwise, will assume single process. */
class CommunicatorT
{
  public:

	/** \name constructor */
	/*@{*/
	/** create communicator including all processes */
	CommunicatorT(void); 
	
	/** copy constructor */
	CommunicatorT(const CommunicatorT& source); 
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
	
	/** (re-)set the logging stream */
	void SetLog(ostream& log);
	
	/** logging level */
	LogLevelT LogLevel(void) const { return fLogLevel; };

	/** (re-)set logging level */
	void SetLogLevel(LogLevelT log_level) { fLogLevel = log_level; };
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

	/** \name gather from all to one */
	/*@{*/	
	/** gather single integer. Called by destination process. 
	 * \param gather destination for gathered values. Must be dimensioned
	 *        to length CommunicatorT::Size before the call. */
	void Gather(int a, nArrayT<int>& gather) const;

	/** gather single integer. Called by sending processes. */
	void Gather(int a, int destination) const;
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
	 * \param counts amount of data sent from each process
	 * \param displacements offset in the destination array where the data
	 *        from the corresponding rank should be written
	 * \param my data sent from this processor
	 * \param gather the destination for the gathered data */
	/*@{*/
	/** gather double's */
	void AllGather(const nArrayT<int>& counts, const nArrayT<int>& displacements, 
		const nArrayT<double>& my, nArrayT<double>& gather) const;

	/** gather int's */
	void AllGather(const nArrayT<int>& counts, const nArrayT<int>& displacements, 
		const nArrayT<int>& my, nArrayT<int>& gather) const;
	/*@}*/

	/** broadcast character array */
	void Broadcast(int source, ArrayT<char>& data);
	
	/** synchronize all processes */
	void Barrier(void) const;

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

#if 0
/** information for a non-blocking send */
class SendPacketT
{
  public:
  
  	/** message type enumeration */
  	int TypeT {Integer, Double, String};


  private:
	
	/** message tag */ 
  	int fMessageTag;
  	
  	/** message data type */
  	TypeT fType;

	/** requests */
	ArrayT<MPI_Request> fRequest;
  	
  	/** \name pointers to message buffers */
	/*@{*/
  	ArrayT<ArrayT<int>* >    fIntArray;
  	ArrayT<ArrayT<double>* > fDoubleArray;
  	ArrayT<ArrayT<char>* >   fCharArray;
	/*@}*/  	
};
#endif

} // namespace Tahoe 
#endif /* _COMMUNICATOR_T_H_ */
