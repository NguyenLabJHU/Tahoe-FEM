/* $Id: ofstreamT.h,v 1.3 2002-07-02 19:57:13 cjkimme Exp $ */
/* created: paklein (12/30/2000)                                          */

#ifndef _OFSTREAM_T_H_
#define _OFSTREAM_T_H_

/* base class */
#include "ios_fwd_decl.h"
#include <fstream.h>
#include <stddef.h>

/* direct members */
#include "StringT.h"


namespace Tahoe {

class ofstreamT: public ofstream
{
public:

	/* constructors */
	ofstreamT(void);
	ofstreamT(const char* file_name, bool append = false);

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

	/* the filename */
	StringT fFileName;
};

/* inlines */
inline const char* ofstreamT::filename(void) const { return fFileName; }

} // namespace Tahoe 
#endif /* _OFSTREAM_T_H_ */
