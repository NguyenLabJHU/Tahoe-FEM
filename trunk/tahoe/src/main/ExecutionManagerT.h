/* $Id: ExecutionManagerT.h,v 1.7 2004-07-15 08:31:03 paklein Exp $ */
/* created: paklein (08/27/1997) */

#ifndef _EXECMAN_T_H_
#define _EXECMAN_T_H_

#include "Environment.h"

/* direct members */
#include "AutoArrayT.h"
#include "StringT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class StringT;
class CommunicatorT;

/** runs tree of input file driven jobs. Derived classes \a must overload 
 * ExecutionManagerT::RunJob() */
class ExecutionManagerT
{
public:

	/** constructor.
	 * \param argc number of command line arguments passed in
	 * \param argv list of command line arguments
	 * \param job_char first-in-file character signaling a job file
	 * \param batch_char first-in-file character signaling a batch file 
	 * \param comm MP environment
	 * \param jobcharputback set to 1 if job_char should be returned to the 
	 *        input stream before it is passed to the analysis object. */
	ExecutionManagerT(int argc, char* argv[], char job_char, char batch_char,
		CommunicatorT& comm, int jobcharputback = 1);

	/** destructor */
	virtual ~ExecutionManagerT(void);

	/** prompt for and execute input files until "quit" */
	virtual void Run(void);

protected:

	/** MUST be overloaded */
	virtual void RunJob(ifstreamT& in, ostream& status) = 0;

	/* returns the index of the requested option */
	bool CommandLineOption(const char* str) const;
	bool CommandLineOption(const char* str, int& index) const;

	/** add the command line option to the list. \returns true if the option was
	 * added, false otherwise. */
	virtual bool AddCommandLineOption(const char* str);

	/** remove the command line option to the list. \returns true if the option was
	 * removed, false otherwise. */
	virtual bool RemoveCommandLineOption(const char* str);

	/** Recursive dispatch */
	virtual void JobOrBatch(ifstreamT& in, ostream& status);

private:
	
	/* Batch file processing */
	void RunBatch(ifstreamT& in, ostream& status);

	/** prompt user for input */
	void Prompt(const char* prompt, const char* default_input, StringT& line) const;
		
protected:

	/* filetype character codes */
	char fJobChar;
	char fBatchChar;
	
	/** MP environment */
	CommunicatorT& fComm;

	/* put back flag */
	int fJobCharPutBack;

	/* command line arguments */
	AutoArrayT<StringT> fCommandLineOptions;
	
private:  	

	/* batch file recursion depth - safety */
	int fRecursionDepth;
};

/* inlines */
inline bool ExecutionManagerT::CommandLineOption(const char* str) const
{
	int index;
	return CommandLineOption(str, index);
}

} // namespace Tahoe 
#endif /* _EXECMAN_T_H_ */
