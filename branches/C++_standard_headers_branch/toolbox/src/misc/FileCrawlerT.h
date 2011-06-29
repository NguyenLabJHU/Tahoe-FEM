/* $Id: FileCrawlerT.h,v 1.6 2004-02-26 08:56:02 paklein Exp $ */

#ifndef _FILE_CRAWLER_T_H_
#define _FILE_CRAWLER_T_H_

#include "Environment.h"

/* direct members */
#include "ArrayT.h"
#include "StringT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class StringT;

/** run one file or a tree of files. Gives console prompt
 * for a file name which can be either a "job" file or a "batch"
 * file depending on the leading character. Lines beginning with
 * '-' are assumed to be application specific commands, which
 * are passed to derived classes for handling. */
class FileCrawlerT
{
public:

	/** constructor */
	FileCrawlerT(int argc, char* argv[], char job_char, char batch_char);

	/** destructor */
	virtual ~FileCrawlerT(void);

	/** prompt input files until "quit" */
	virtual void Run(void);

protected:

	/** application-specific job execution */
	virtual void RunJob(ifstreamT& in, ostream& status) = 0;

	/** handle batch file command */
	virtual void BatchFileCommand(const StringT& command, ifstreamT& batch) = 0;

	/** batch file processing */
	virtual void RunBatch(ifstreamT& in, ostream& status);

	/** returns the index of the requested option */
	bool CommandLineOption(const char* str) const;
	bool CommandLineOption(const char* str, int& index) const;
	void AddCommandLineOption(const char* str);

	/** recursive dispatch */
	virtual void JobOrBatch(ifstreamT& in, ostream& status);
	
protected:

	/* filetype character codes */
	char fJobChar;
	char fBatchChar;

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

} // namespace Tahoe 
#endif /* _FILE_CRAWLER_T_H_ */
