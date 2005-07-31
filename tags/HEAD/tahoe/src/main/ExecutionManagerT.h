/* $Id: ExecutionManagerT.h,v 1.1.1.1 2001-01-29 08:20:21 paklein Exp $ */
/* created: paklein (08/27/1997)                                          */
/* Manages input file driven jobs.                                        */
/* MUST overload private:RunJob().                                        */

#ifndef _EXECMAN_T_H_
#define _EXECMAN_T_H_

#include "Environment.h"

/* direct members */
#include "ArrayT.h"
#include "StringT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class StringT;

class ExecutionManagerT
{
public:

	/* Constructors */
	ExecutionManagerT(int argc, char* argv[], char job_char, char batch_char,
		int jobcharputback = 1);

	/* Destructor */
	virtual ~ExecutionManagerT(void);

	/* Prompt input files until "quit" */
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
	void AddCommandLineOption(const char* str);

private:

	/* Recursive dispatch */
	void JobOrBatch(ifstreamT& in, ostream& status);
	
	/* Batch file processing */
	void RunBatch(ifstreamT& in, ostream& status);
	
protected:

	/* filetype character codes */
	char fJobChar;
	char fBatchChar;

	/* put back flag */
	int fJobCharPutBack;

	/* command line arguments */
	ArrayT<StringT> fCommandLineOptions;
	
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
