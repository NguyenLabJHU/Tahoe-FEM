/* $Id: ExecutionManagerT.h,v 1.2 2002-01-03 19:10:28 paklein Exp $ */
/* created: paklein (08/27/1997) */

#ifndef _EXECMAN_T_H_
#define _EXECMAN_T_H_

#include "Environment.h"

/* direct members */
#include "AutoArrayT.h"
#include "StringT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class StringT;

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
	 * \param jobcharputback set to 1 if job_char should be returned to the 
	 *        input stream before it is passed to the analysis object. */
	ExecutionManagerT(int argc, char* argv[], char job_char, char batch_char,
		int jobcharputback = 1);

	/** destructor */
	virtual ~ExecutionManagerT(void);

	/** prompt for and execute input files until "quit" */
	virtual void Run(void);

protected:

	/* run function */
	void Run_serial(void);
	void Run_parallel(void);

	/* MUST be overloaded */
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

private:

	/* Recursive dispatch */
	void JobOrBatch(ifstreamT& in, ostream& status);
	
	/* Batch file processing */
	void RunBatch(ifstreamT& in, ostream& status);
	
	/** open stream with prompt. Input names starting with '-' will be
	 * treated as options that will be handled with ExecutionManagerT::AddCommandLineOption.
	 * \param prompt prompt that will appear on command line
	 * \param skipname signal to exit loop to find string 
	 * \param defaultname name returned if an empty string is given 
	 * \param in stream to open. */
	int OpenWithPrompt(const char* prompt, const char* skipname,
		const char* defaultname, ifstreamT& in);
		
protected:

	/* filetype character codes */
	char fJobChar;
	char fBatchChar;

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

#endif /* _EXECMAN_T_H_ */