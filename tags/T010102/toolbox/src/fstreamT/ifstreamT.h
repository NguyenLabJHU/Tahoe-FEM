/* $Id: ifstreamT.h,v 1.6 2001-12-16 23:50:57 paklein Exp $ */
/* created: paklein (03/03/1999) */

#ifndef _IFSTREAM_T_H_
#define _IFSTREAM_T_H_

#include "Environment.h"

/* base class */
#include "ios_fwd_decl.h"
#include <fstream.h>
#include <stddef.h>

/* direct members */
#include "StringT.h"

/** input file stream with extended capabilities */
class ifstreamT: public ifstream
{
public:

	/* constructors */
	ifstreamT(void);
	ifstreamT(const char* file_name);
	ifstreamT(char marker);
	ifstreamT(char marker, const char* file_name);

	/* open stream */
	void open(const char* file_name);
	int open(const char* prompt, const char* skipname,
		const char* defaultname = NULL);
	int is_open(void);
	
	/* close stream */
	void close(void);

	/* comment marker */
	void set_marker(char marker);
	void clear_marker(void);
	char comment_marker(void) const;
	int skip_comments(void) const;
	
	/* put a character back in the stream */
	istream& putback(char a);
	
	/* return the next character (skipping whitespace and comments)
	 * without removing it from the stream */
	char next_char(void);
	
	/* return the filename - NULL if no file is open */
	const char* filename(void) const;
	
	/* set file name string - does not change stream */
	void set_filename(const char* name);

	/* adjusting stream position - returns actual number of rewound lines */
	int rewind(int num_lines = 1);

	/* advance to the end of the line (or next 255 characters) */
	void clear_line(void);

	/* advances passed comments */
	void do_skip_comments(void);

	/* extraction of streams */
	ifstreamT& operator>>(bool& a);

	/** stream search. read lines from the stream looking for key
	 * \param key search string
	 * \param line entire line from string containing key
	 * \return true if key found, else false */
	bool FindString(const char* key, StringT& line);

	/** C-C++(2.4.6)/MSL doesn't to path's right. temporarily needed to
	 * remove "ups" and "downs" from pathnames. Used only for Mac */
	void FixPath(const char* path_old, StringT& path) const; 

private:

	/* open stream with prompt - return 1 if successful */
	int OpenWithPrompt(const char* prompt, const char* skipname,
		const char* defaultname);
	
private:

	/* comment marker */
	int  fSkipComments;
	char fMarker;
	
	/* the filename */
	StringT fFileName;
};

/* inlines */
inline char ifstreamT::comment_marker(void) const { return fMarker; }
inline int ifstreamT::skip_comments(void) const { return fSkipComments; }
inline const char* ifstreamT::filename(void) const { return fFileName; }

inline void ifstreamT::set_marker(char marker)
{
	fSkipComments = 1;
	fMarker = marker;
}

inline void ifstreamT::clear_marker(void) { fSkipComments = 0; }

/* extraction operator */
template <class TYPE> ifstreamT& operator>>(ifstreamT& str, TYPE& data);
template <class TYPE>
ifstreamT& operator>>(ifstreamT& str, TYPE& data)
{
	/* advance */
	str.do_skip_comments();
		
	/* ANSI */
	ifstream& ifstr = str;
	ifstr >> data;

	/* advance */
	str.do_skip_comments();
	
	return str;
};

#endif /* _IFSTREAM_X_H_ */
