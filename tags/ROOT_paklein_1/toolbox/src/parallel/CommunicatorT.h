/* $Id: CommunicatorT.h,v 1.5 2002-08-21 07:19:29 paklein Exp $ */
#ifndef _COMMUNICATOR_T_H_
#define _COMMUNICATOR_T_H_

#include "Environment.h"

#ifdef __MPI__
#include "mpi.h"
#else
#define MPI_Comm long
#define MPI_Request long
#endif

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
template <class TYPE> class ArrayT;

/** interface to handle process to process communication. If compiled
 * with -D__MPI__, will implemented process-to-process communication
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
	enum LogLevelT {kSilent = 0, /**< serious error messages only */
	               kAddress = 1, /**< log destination/source */
	                   kAll = 2  /**< log everything including message data */ };

	/** write message to log 
	 * \param caller the calling subroutine 
	 * \param message the log message
	 * \param force write message regardless of the logging level */
	void Log(const char* caller, const char* message, bool force = false) const;

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
	
	/** gather single integer. Called by destination process. */
	void Gather(int a, ArrayT<int>& gather) const;

	/** gather single integer. Called by sending processes. */
	void Gather(int a, int destination) const;

	/** gather single integer to all processes. */
	void AllGather(int a, ArrayT<int>& gather) const;

	/** broadcast character array */
	void Broadcast(ArrayT<char>& data);
	
	/** synchronize all processes */
	void Barrier(void) const;

  protected:
  
  	/** time as a string */
	const char* WallTime(void) const;

	/** write log header */
	ostream& LogHead(const char* caller) const;

  private:
  
  	/** \name initialize/finalize MPI environment */
  	/*@{*/
	void Init(void);
	void Finalize(void);
  	/*@}*/
  	
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
#ifdef __MPI__
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