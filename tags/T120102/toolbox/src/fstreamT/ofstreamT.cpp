/* $Id: ofstreamT.cpp,v 1.5 2002-10-20 22:39:00 paklein Exp $ */
/* created: paklein (12/30/2000) */

#include "ofstreamT.h"

/* ANSI */
#include <iostream.h>
#include <string.h>

#include "fstreamT.h"
#include "Environment.h"
#include "ExceptionT.h"

#include "StringT.h"

/* parameter */

using namespace Tahoe;

const int kLineLength = 255;

/* constructors */
ofstreamT::ofstreamT(void)
{
	format_stream(*this);
}	

ofstreamT::ofstreamT(const char* file_name, bool append)
{
	format_stream(*this);
	if (append)
		open_append(file_name);
	else
		open(file_name);
}

/* open stream */
void ofstreamT::open(const char* stream)
{
	/* close stream if already open */
	if (is_open()) close();

	/* copy file name */
	fFileName = stream;
	fFileName.ToNativePathName();

	/* ANSI */
	ofstream::open(fFileName);
}

void ofstreamT::open_append(const char* stream)
{
	/* close stream if already open */
	if (is_open()) close();

	/* copy file name */
	fFileName = stream;
	fFileName.ToNativePathName();

	/* ANSI */
	ofstream::open(stream, ios::app);
}

int ofstreamT::is_open(void)
{
#ifdef __MWERKS__
	return ofstream::is_open();
#else
// is_open is only defined for filebuf not ostream or istream,
// and isn't defined as const
ofstream* non_const_ofstr = (ofstream*) this;
filebuf* fbuf = non_const_ofstr->rdbuf();
return fbuf->is_open();
#endif
}

/* close stream */
void ofstreamT::close(void)
{
	if (fstreamT::need_MW_workaround())
	{
		/* ANSI */
		if (is_open()) ofstream::close();

		/* clear name */
		fFileName.Clear();
	}
	else /* original code */
	{
	/* ANSI */
	ofstream::close();

	/* clear name */
	fFileName.Clear();
	}
}

/* set stream formats */
void ofstreamT::format_stream(ostream& out)
{
	out.precision(kPrecision);
	out.setf(ios::showpoint);
	out.setf(ios::right, ios::adjustfield);
	out.setf(ios::scientific, ios::floatfield);
}
