/* $Id: ofstreamT.h,v 1.5 2004-06-17 06:37:52 paklein Exp $ */
/* created: paklein (12/30/2000) */
#ifndef _OFSTREAM_T_H_
#define _OFSTREAM_T_H_

/* base class */
#include "fstreamT.h"
#include "ios_fwd_decl.h"
#include <fstream.h>
#include <stddef.h>

namespace Tahoe {

class ofstreamT: public ofstream, public fstreamT
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

	/* set stream formats */
	static void format_stream(ostream& out);

private:

	/** copy constructor not allowed */
	ofstreamT(const ofstreamT&);
};

} /* namespace Tahoe */

#endif /* _OFSTREAM_T_H_ */
