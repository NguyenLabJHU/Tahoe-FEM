/* $Id: ifstreamT.h,v 1.14 2004-06-17 06:37:52 paklein Exp $ */
/* created: paklein (03/03/1999) */
#ifndef _IFSTREAM_T_H_
#define _IFSTREAM_T_H_

/* base class */
#include "fstreamT.h"
#include <fstream.h>
#include <stddef.h>

#include "ios_fwd_decl.h"

namespace Tahoe {

/** input file stream with extended capabilities */
class ifstreamT: public ifstream, public fstreamT
{
public:

	/** \name constructors */
	/*@{*/
	ifstreamT(void);
	ifstreamT(const char* file_name);
	ifstreamT(char marker);
	ifstreamT(char marker, const char* file_name);
	/*@}*/

	/** \name opening stream */
	/*@{*/
	void open(const char* file_name);
	int open(const char* prompt, const char* skipname,
		const char* defaultname = NULL);
	int is_open(void);
	/*@}*/
	
	/** close stream */
	void close(void);

	/** \name comment marker */
	/*@{*/
	void set_marker(char marker);
	void clear_marker(void);
	char comment_marker(void) const;
	int skip_comments(void) const;
	/*@}*/
	
	/** put a character back in the stream */
	istream& putback(char a);
	
	/** return the next character (skipping whitespace and comments)
	 * without removing it from the stream */
	char next_char(void);
	
	/** get the next line from stream. Ignores comment lines */
	ifstreamT& getline(char* s, int n, char delimiter = '\n');

	/** adjusting stream position
	 * \return the actual number of rewound lines */
	int rewind(int num_lines = 1);

	/** advance to the end of the line (or next 255 characters) */
	void clear_line(void);

	/** advances passed comments */
	void do_skip_comments(void);

	/** extraction of streams */
	ifstreamT& operator>>(bool& a);

	/** stream search. read lines from the stream looking for key
	 * \param key search string
	 * \param line entire line from string containing key
	 * \return true if key found, else false */
	bool FindString(const char* key, StringT& line);

private:

	/** open stream with prompt
	 * \return 1 if successful */
	int OpenWithPrompt(const char* prompt, const char* skipname,
		const char* defaultname);
	
private:

	/** \name comment marker */
	/*@{*/
	int  fSkipComments;
	char fMarker;
	/*@}*/
};

/* inlines */
inline char ifstreamT::comment_marker(void) const { return fMarker; }
inline int ifstreamT::skip_comments(void) const { return fSkipComments; }

inline void ifstreamT::set_marker(char marker)
{
	fSkipComments = 1;
	fMarker = marker;
}

inline void ifstreamT::clear_marker(void) { fSkipComments = 0; }

/* get the next line from stream. Ignores comment lines */
inline ifstreamT& ifstreamT::getline(char* s, int n, char delimiter)
{
	/* advance */
	do_skip_comments();
		
	/* ANSI */
	ifstream::getline(s, n, delimiter);

	/* advance */
	do_skip_comments();

	return *this;
}

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
}

} // namespace Tahoe 
#endif /* _IFSTREAM_X_H_ */
