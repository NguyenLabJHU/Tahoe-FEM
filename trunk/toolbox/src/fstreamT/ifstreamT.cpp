/* $Id: ifstreamT.cpp,v 1.27 2004-11-12 20:58:37 paklein Exp $ */
/* created: paklein (03/03/1999) */
#include "ifstreamT.h"

/* ANSI */
#include <iostream.h>
#include <string.h>
#include <ctype.h>

#include "fstreamT.h"

using namespace Tahoe;

/* parameter */
const int kLineLength = 255;

/* static variables */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<ifstreamT*>::fByteCopy = true; // array behavior
} /* namespace Tahoe */

/* constructors */
ifstreamT::ifstreamT(void):
	fSkipComments(0),
	fMarker('0')
{

}	

ifstreamT::ifstreamT(const char* file_name):
	fSkipComments(0),
	fMarker('0')
{
	open(file_name);
}

ifstreamT::ifstreamT(char marker):
	fSkipComments(1),
	fMarker(marker)
{

}	

ifstreamT::ifstreamT(char marker, const char* file_name):
	fSkipComments(1),
	fMarker(marker)
{
	open(file_name);
}

/* open stream */
void ifstreamT::open(const char* file_name)
{
	/* close stream if already open */
	if (is_open()) close();

	/* translate name */
	fFileName = file_name;
	fFileName.ToNativePathName();
	
	//TEMP - problem with CW7
	StringT old = fFileName;
	fstreamT::FixPath(old, fFileName);
	//TEMP

	/* ANSI */
	ifstream::open(fFileName);
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
	if (fstreamT::need_MW_workaround())
	{
		/* ANSI */
		if (is_open()) ifstream::close();

		/* clear name */
		fFileName.Clear();
	}
	else /* original code */
	{
	/* ANSI */
	ifstream::close();

	/* clear name */
	fFileName.Clear();
	}
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
	if (good()) ifstream::putback(c);

	return c;
}

/* put a character back in the stream */
istream& ifstreamT::putback(char a)
{
	/* do not allow putback with skip comments activated since
	 * the stream has probably been advanced */
	if (fSkipComments)
	{
		cout << "\n ifstreamT::putback: not allowed while comment skipping is enabled" << endl;
		throw ExceptionT::kGeneralFail;
	}
	
	/* inherited */
	return ifstream::putback(a);
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
#if defined(__GCC_3__) || defined( __SUNPRO_CC) || (defined(__GNUC__) && defined(__PGI__)) || (defined(__DEC__) && defined (__USE_STD_IOSTREAM)) || (defined(__JANUS__) && defined(__ROGUE_STL__))
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
	
	int skip = fSkipComments;
	fSkipComments = 0;
	getline(line, kLineLength-1);
	fSkipComments = skip;
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
		ifstream::putback(c);
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
		throw ExceptionT::kBadInputValue;
	}
	
	return *this;
}

/* stream search */
bool ifstreamT::FindString(const char* key, StringT& line)
{
	/* clear return string */
	line.Clear();

	/* read line-by-line */
	bool found = false;
	while (!found && good())
	{
		char buffer[kLineLength];
		getline(buffer, kLineLength-1);
		
		/* long line */
		if (!good()) {		
			/* correct anything but eof */
			if (!eof()) clear();
		}

		/* found string */
		if (strstr(buffer, key) != NULL)
		{
			line = buffer;
			found = true;
		}
	}
	return found;
}

/*************************************************************************
* Private
*************************************************************************/

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
#ifdef __SGI__
			cout.flush();
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
#ifdef __SGI__
			cout.flush();
#endif					
			cin >> newfilename;
		}
		
		/* clear to end of line */
		fstreamT::ClearLine(cin);

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
				fFileName = newfilename;
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
				throw ExceptionT::kGeneralFail;
			}
		}
	}	
}
