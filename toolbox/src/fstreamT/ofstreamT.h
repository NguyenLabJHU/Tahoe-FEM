/* $Id: ofstreamT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (12/30/2000)                                          */

#ifndef _OFSTREAM_T_H_
#define _OFSTREAM_T_H_

/* base class */
#include "ios_fwd_decl.h"
#include <fstream.h>
#include <stddef.h>

class ofstreamT: public ofstream
{
public:

	/* constructors */
	ofstreamT(void);
	ofstreamT(const char* file_name, bool append = false);

	/* destructor */
	~ofstreamT(void);

	/* open stream */
	void open(const char* file_name);
	void open_append(const char* file_name);
	int is_open(void);
	
	/* close stream */
	void close(void);

	/* return the filename - NULL if no file is open */
	const char* filename(void) const;

	/* set stream formats */
	static void format_stream(ostream& out);

private:

	/* copy the string to fFileName */
	void CopyName(const char* filename);

private:

	/* the filename */
	char* fFileName;

	/* NULL file name */
	static const char fNULLFileName;
};

/* inlines */
inline const char* ofstreamT::filename(void) const
{
	return (fFileName != NULL) ? fFileName : &fNULLFileName;
}

#endif /* _OFSTREAM_T_H_ */
