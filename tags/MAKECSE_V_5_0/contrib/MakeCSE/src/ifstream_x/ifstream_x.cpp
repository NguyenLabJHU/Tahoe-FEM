/*
 * File: ifstream_x.cpp
 */

/*
 * created      : PAK (03/03/99)
 * last modified: PAK (03/03/99)
 */

/* interface */
#include "ifstream_x.h"

/* ANSI */
#include <iostream.h>
#include <string.h>
#include <ctype.h>

#include "Environment.h"
#include "ExceptionCodes.h"

#include "StringT.h"

/* parameter */
const int kLineLength = 255;

/* constructors */
ifstream_x::ifstream_x(void):
	fSkipComments(0),
	fMarker('0'),
	fFileName(NULL)
{

}	

ifstream_x::ifstream_x(const char* stream):
	fifstream(stream),
	fSkipComments(0),
	fFileName(NULL)
{
	CopyName(stream);
}

ifstream_x::ifstream_x(char marker):
	fSkipComments(1),
	fMarker(marker),
	fFileName(NULL)
{

}	

ifstream_x::ifstream_x(char marker, const char* stream):
	fifstream(stream),
	fSkipComments(1),
	fMarker(marker),
	fFileName(NULL)
{
	CopyName(stream);
}

/* destructor */
ifstream_x::~ifstream_x(void)
{
	delete[] fFileName;
	fFileName = NULL;
}

/* open stream */
void ifstream_x::open(const char* stream)
{
	/* close stream if already open */
	if (is_open()) close();

	/* ANSI */
	fifstream.open(stream);
	
	/* store file name */
	CopyName(stream);
}

int ifstream_x::open(const char* prompt, const char* skipname,
	const char* defaultname)
{
	/* close stream if already open */
	if (is_open()) close();

	/* attempt open */
	return OpenWithPrompt(prompt, skipname, defaultname);
}

int ifstream_x::is_open(void)
{ 
#ifdef __MWERKS__
	return fifstream.is_open(); 
#else
// is_open is only defined for filebuf not ostream or istream, 
// and isn't defined as const
  ifstream* non_const_ifstr = (ifstream*) &fifstream;
  filebuf* fbuf = non_const_ifstr->rdbuf();
  return fbuf->is_open();
#endif
}

/* close stream */
void ifstream_x::close(void)
{
	/* ANSI */
	fifstream.close();

	/* clear name */
	delete[] fFileName;
	fFileName = NULL;
}

/* return the next character (skipping whitespace and comments)
 * without removing it from the stream */
char ifstream_x::next_char(void)
{
	/* advance */
	do_skip_comments();

	/* get next character */
	char c;	
	fifstream.get(c);
	while (fifstream.good() && isspace(c)) fifstream.get(c);	

	/* restore last character */
	if (fifstream.good()) fifstream.putback(c);

	return c;
}

/* force name */
void ifstream_x::set_filename(const char* name)
{
	/* store file name */
	CopyName(name);
}

/* adjusting stream position, returns the number of lines rewound */
int ifstream_x::rewind(int num_lines)
{	
	int line_count = 0;
	char c;

#ifdef __MWERKS__
#ifdef _MW_MSL_ // for CWPro 5
	ifstream& ifstr = *this;
	streampos pos = ifstr.tellg();
	while (pos >= 0 && line_count < num_lines)
	{
		pos -= 1;
		ifstr.seekg(pos);
		c = ifstr.peek();
		if (c == '\n') line_count++;
	}
#else // not MSL
	/* must work directly with stream buffer */
	filebuf* buf = fifstream.rdbuf();
	streampos pos = buf->pubseekoff(0, ios::cur);
	while (pos.offset() >= 0 && line_count < num_lines)
	{
		pos = buf->pubseekoff(-1, ios::cur);
		c = fifstream.peek();
		if (c == '\n') line_count++;
	}
#endif // _MW_MSL_
#else  // not CodeWarrior 
#ifdef __SUNPRO_CC
	ifstream& ifstr = *this;
	streampos pos = ifstr.tellg();
	while (pos >= 0 && line_count < num_lines)
	{
		pos -= 1;
		ifstr.seekg(pos);
		c = ifstr.peek();
		if (c == '\n') line_count++;
	}
#else // not Sun
	ifstream& ifstr = *this;
	streampos pos = ifstr.tellg();
	while (pos >= 0 && line_count < num_lines)
	{
		ifstr.seekg(pos--);
		c = ifstr.peek();
		if (c == '\n') line_count++;
	}
#endif // __SUNPRO_CC
#endif // __MWERKS__
	
	/* advance passed newline */
	if (c == '\n') get(c);
	
	return line_count;
}

/* advance to the end of the line */
void ifstream_x::clear_line(void)
{
	char line[kLineLength];
	fifstream.getline(line, kLineLength-1);
}

/* advances passed comments */
void ifstream_x::do_skip_comments(void)
{
	if (!fSkipComments || !fifstream.good()) return;

	char c;
	fifstream.get(c);
	while (fifstream.good() && (isspace(c) || c == fMarker))
	{
		/* remove comment line */
		if (c == fMarker) clear_line();
		
		fifstream.get(c);	
	}

	/* don't die while skipping comments */
	if (fifstream.good()) 
		fifstream.putback(c);
	else
		fifstream.clear();
}

/* extraction of streams */
ifstream_x& ifstream_x::operator>>(bool& a)
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
		cout << "\n ifstream_x::operator>>bool&: expecting 0 or 1 from stream" << endl;
		throw eBadInputValue;
	}
	
	return *this;
}

/*************************************************************************
 * Private
 *************************************************************************/

/* copy the string to fFileName */
void ifstream_x::CopyName(const char* filename)
{
	/* no copies to self */
	if (filename == fFileName) return;

	/* free existing memory */
	delete[] fFileName;
	
	/* check file name */
	if (strlen(filename) == 0)
	{
		cout << "\n ifstream_x::CopyName: zero length filename" << endl;
		throw eGeneralFail;
	}
	
	/* allocate new memory */
	fFileName = new char[strlen(filename) + 1];
	if (!fFileName) 
	{
		cout << "\n ifstream_x::CopyName: out of memory" << endl;
		throw eOutOfMemory;
	}
	
	/* copy in */
	memcpy(fFileName, filename, sizeof(char)*(strlen(filename) + 1));
}

/* open stream with prompt - return 1 if successful */
int ifstream_x::OpenWithPrompt(const char* prompt, const char* skipname,
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

			char test;
			cin.get(test); //clear newline from last time
			if (test == '\n')
				cin.get(test);
			
			/* new filename */
			if (test != '\n')
			{
				cin.putback(test);
				cin >> newfilename;
			}
			else
			{
				/* restore newline */
				cin.putback(test); 
				
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

		/* check exit */
		if (strcmp(newfilename, skipname) == 0)
			return 0;
		/* attempt open file */
		else
		{
			/* convert to native path */
			newfilename.ToNativePathName();
		
			/* attempt open */
			fifstream.open(newfilename);
			
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
				fifstream.clear();
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
