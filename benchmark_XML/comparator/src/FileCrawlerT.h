/* $Id: FileCrawlerT.h,v 1.1 2001-06-11 02:10:12 paklein Exp $ */

#ifndef _FILE_CRAWLER_T_H_
#define _FILE_CRAWLER_T_H_

#include "Environment.h"

/* direct members */
#include "ArrayT.h"
#include "StringT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class StringT;

class FileCrawlerT
{
public:

	/* constructor */
	FileCrawlerT(int argc, char* argv[], char job_char, char batch_char,
		int jobcharputback = 1);

	/* destructor */
	virtual ~FileCrawlerT(void);

	/* prompt input files until "quit" */
	virtual void Run(void);

protected:

	/* MUST be overloaded */
	virtual void RunJob(ifstreamT& in, ostream& status) = 0;

	/* batch file processing */
	virtual void RunBatch(ifstreamT& in, ostream& status);

	/* returns the index of the requested option */
	bool CommandLineOption(const char* str) const;
	bool CommandLineOption(const char* str, int& index) const;
	void AddCommandLineOption(const char* str);

private:

	/* recursive dispatch */
	virtual void JobOrBatch(ifstreamT& in, ostream& status);
	
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
inline bool FileCrawlerT::CommandLineOption(const char* str) const
{
	int index;
	return CommandLineOption(str, index);
}

#endif /* _FILE_CRAWLER_T_H_ */
