/* $Id: ofstreamT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (12/30/2000)                                          */
/* interface                                                              */

#include "ofstreamT.h"

/* ANSI */
#include <iostream.h>
#include <string.h>

#include "Environment.h"
#include "ExceptionCodes.h"

#include "StringT.h"

/* parameter */
const int kLineLength = 255;

/* static variables */
const char ofstreamT::fNULLFileName = '\0';

/* constructors */
ofstreamT::ofstreamT(void):
	fFileName(NULL)
{
	format_stream(*this);
}	

ofstreamT::ofstreamT(const char* file_name, bool append):
	fFileName(NULL)
{
	format_stream(*this);
	if (append)
		open_append(file_name);
	else
		open(file_name);
}

/* destructor */
ofstreamT::~ofstreamT(void)
{
	delete[] fFileName;
	fFileName = NULL;
}

/* open stream */
void ofstreamT::open(const char* stream)
{
	/* close stream if already open */
	if (is_open()) close();

	/* ANSI */
	ofstream::open(stream);
	
	/* store file name */
	CopyName(stream);
}

void ofstreamT::open_append(const char* stream)
{
	/* close stream if already open */
	if (is_open()) close();

	/* ANSI */
	ofstream::open(stream, ios::app);
	
	/* store file name */
	CopyName(stream);
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
	/* ANSI */
	ofstream::close();

	/* clear name */
	delete[] fFileName;
	fFileName = NULL;
}

/* set stream formats */
void ofstreamT::format_stream(ostream& out)
{
	out.precision(kPrecision);
	out.setf(ios::showpoint);
	out.setf(ios::right, ios::adjustfield);
	out.setf(ios::scientific, ios::floatfield);
}

/*************************************************************************
* Private
*************************************************************************/

/* copy the string to fFileName */
void ofstreamT::CopyName(const char* filename)
{
	/* no copies to self */
	if (filename == fFileName) return;

	/* free existing memory */
	delete[] fFileName;
	
	/* check file name */
	if (strlen(filename) == 0)
	{
		cout << "\n ofstreamT::CopyName: zero length filename" << endl;
		throw eGeneralFail;
	}
	
	/* allocate new memory */
	fFileName = new char[strlen(filename) + 1];
	if (!fFileName)
	{
		cout << "\n ofstreamT::CopyName: out of memory" << endl;
		throw eOutOfMemory;
	}
	
	/* copy in */
	memcpy(fFileName, filename, sizeof(char)*(strlen(filename) + 1));
}
