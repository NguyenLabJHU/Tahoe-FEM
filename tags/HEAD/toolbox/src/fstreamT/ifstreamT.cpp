/* $Id: ifstreamT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (03/03/1999)                                          */
/* interface                                                              */

#include "ifstreamT.h"

/* ANSI */
#include <iostream.h>
#include <string.h>
#include <ctype.h>

#include "Environment.h"
#include "ExceptionCodes.h"

#include "StringT.h"

/* parameter */
const int kLineLength = 255;

/* static variables */
const char ifstreamT::fNULLFileName = '\0';
const bool ArrayT<ifstreamT*>::fByteCopy = true; // array behavior

/* constructors */
ifstreamT::ifstreamT(void):
	fSkipComments(0),
	fMarker('0'),
	fFileName(NULL)
{

}	

ifstreamT::ifstreamT(const char* file_name):
	fSkipComments(0),
	fMarker('0'),
	fFileName(NULL)
{
	open(file_name);
}

ifstreamT::ifstreamT(char marker):
	fSkipComments(1),
	fMarker(marker),
	fFileName(NULL)
{

}	

ifstreamT::ifstreamT(char marker, const char* file_name):
	fSkipComments(1),
	fMarker(marker),
	fFileName(NULL)
{
	open(file_name);
}

/* destructor */
ifstreamT::~ifstreamT(void)
{
	delete[] fFileName;
	fFileName = NULL;
}

/* open stream */
void ifstreamT::open(const char* file_name)
{
	/* close stream if already open */
	if (is_open()) close();

	/* ANSI */
	ifstream::open(file_name);
	
	/* store file name */
	CopyName(file_name);
}

int ifstreamT::open(const char* prompt, const char* skipname,
	const char* defaultname)
{
	/* close stream if already open */
	if (is_open()) close();

	/* attempt open */
	return OpenWithPrompt(prompt, skipname, defaultname);
}

int ifstreamT::is_open(void)
{
#ifdef __MWERKS__
	return ifstream::is_open();
#else
// is_open is only defined for filebuf not ostream or istream,
// and isn't defined as const
ifstream* non_const_ifstr = (ifstream*) this;
filebuf* fbuf = non_const_ifstr->rdbuf();
return fbuf->is_open();
#endif
}

/* close stream */
void ifstreamT::close(void)
{
	/* ANSI */
	ifstream::close();

	/* clear name */
	delete[] fFileName;
	fFileName = NULL;
}

/* return the next character (skipping whitespace and comments)
* without removing it from the stream */
char ifstreamT::next_char(void)
{
	/* advance */
	do_skip_comments();

	/* get next character */
	char c;	
	get(c);
	while (good() && isspace(c)) get(c);	

	/* restore last character */
	if (good()) putback(c);

	return c;
}

/* set file name string - does not change stream */
void ifstreamT::set_filename(const char* name)
{
	/* store file name */
	CopyName(name);
}

/* adjusting stream position, returns the number of lines rewound */
int ifstreamT::rewind(int num_lines)
{	
	int line_count = 0;
	char c;

#ifdef __MWERKS__
#ifdef _MW_MSL_ // for CWPro 5
	streampos pos = tellg();
	while (pos >= 0 && line_count < num_lines)
	{
		pos -= 1;
		seekg(pos);
		c = peek();
		if (c == '\n') line_count++;
	}
#else // not MSL
	/* must work directly with stream buffer */
	filebuf* buf = rdbuf();
	streampos pos = buf->pubseekoff(0, ios::cur);
	while (pos.offset() >= 0 && line_count < num_lines)
	{
		pos = buf->pubseekoff(-1, ios::cur);
		c = peek();
		if (c == '\n') line_count++;
	}
#endif // _MW_MSL_
#else  // not CodeWarrior
#ifdef __SUNPRO_CC
	streampos pos = tellg();
	while (pos >= 0 && line_count < num_lines)
	{
		pos -= 1;
		seekg(pos);
		c = peek();
		if (c == '\n') line_count++;
	}
#else // not Sun
	streampos pos = tellg();
	while (pos >= 0 && line_count < num_lines)
	{
		seekg(pos--);
		c = peek();
		if (c == '\n') line_count++;
	}
#endif // __SUNPRO_CC
#endif // __MWERKS__
	
	/* advance passed newline */
	if (c == '\n') get(c);
	
	return line_count;
}

/* advance to the end of the line */
void ifstreamT::clear_line(void)
{
	char line[kLineLength];
	getline(line, kLineLength-1);
}

/* advances passed comments */
void ifstreamT::do_skip_comments(void)
{
	if (!fSkipComments || !good()) return;

	char c;
	get(c);
	while (good() && (isspace(c) || c == fMarker))
	{
		/* remove comment line */
		if (c == fMarker) clear_line();
		
		get(c);	
	}

	/* don't die while skipping comments */
	if (good())
		putback(c);
	else
		clear();
}

/* extraction of streams */
ifstreamT& ifstreamT::operator>>(bool& a)
{
	/* read int */
	int i;
	(*this) >> i;
	
	if (i == 0)
		a = false;
	else if (i == 1)
		a = true;
	else
	{
		cout << "\n ifstreamT::operator>>bool&: expecting 0 or 1 from stream" << endl;
		throw eBadInputValue;
	}
	
	return *this;
}

/*************************************************************************
* Private
*************************************************************************/

/* copy the string to fFileName */
void ifstreamT::CopyName(const char* filename)
{
	/* no copies to self */
	if (filename == fFileName) return;

	/* free existing memory */
	delete[] fFileName;
	
	/* check file name */
	if (strlen(filename) == 0)
	{
		cout << "\n ifstreamT::CopyName: zero length filename" << endl;
		throw eGeneralFail;
	}
	
	/* allocate new memory */
	fFileName = new char[strlen(filename) + 1];
	if (!fFileName)
	{
		cout << "\n ifstreamT::CopyName: out of memory" << endl;
		throw eOutOfMemory;
	}
	
	/* copy in */
	memcpy(fFileName, filename, sizeof(char)*(strlen(filename) + 1));
}

/* open stream with prompt - return 1 if successful */
int ifstreamT::OpenWithPrompt(const char* prompt, const char* skipname,
	const char* defaultname)
{
	StringT newfilename;
	int maxtry = 20;
	int count  = 0;
	while (1)
	{
		cout << '\n' << prompt << "\n(\"" << skipname << "\" to exit";

		/* default */
		if (defaultname != NULL && strlen(defaultname) > 0)
		{
			cout << ", <RETURN> for \"" << defaultname << "\"): ";
#if (defined __SGI__ && defined __MPI__)
			cout << '\n';
#endif
			
			/* new filename */
			char test = cin.peek();
			if (test != '\n')
			{
				/* take first word */
				cin >> newfilename;
			}
			else
			{
				/* copy default */
				newfilename = defaultname;
			}				
		}
		else
		{
			cout << "): ";
#if (defined __SGI__ && defined __MPI__)
			cout << '\n';
#endif					
			cin >> newfilename;
		}
		
		/* clear to end of line */
		char line[255];
		cin.getline(line, 254);

		/* check exit */
		if (strcmp(newfilename, skipname) == 0)
			return 0;
		/* attempt open file */
		else
		{
			/* convert to native path */
			newfilename.ToNativePathName();
		
			/* attempt open */
			open(newfilename);
			
			/* check file open */
			if (is_open())
			{
				/* store file name */
				CopyName(newfilename);
			
				return 1;	
			}
			else
			{
				cout << "\nError: filename: " << newfilename << " not found\n";
				clear();
			}
			
			/* could not find file */
			if (++count == maxtry)
			{
				cout << "\n StringT::OpenInputStream: could not find file after ";
				cout << maxtry << " iterations" << endl;
				throw eGeneralFail;
			}
		}
	}	
}
