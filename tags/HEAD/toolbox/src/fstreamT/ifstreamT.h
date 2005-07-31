/* $Id: ifstreamT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (03/03/1999)                                          */

#ifndef _IFSTREAM_T_H_
#define _IFSTREAM_T_H_

/* base class */
#include "ios_fwd_decl.h"
#include <fstream.h>
#include <stddef.h>

class ifstreamT: public ifstream
{
public:

	/* constructors */
	ifstreamT(void);
	ifstreamT(const char* file_name);
	ifstreamT(char marker);
	ifstreamT(char marker, const char* file_name);

	/* destructor */
	~ifstreamT(void);

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

private:

	/* copy the string to fFileName */
	void CopyName(const char* filename);

	/* open stream with prompt - return 1 if successful */
	int OpenWithPrompt(const char* prompt, const char* skipname,
		const char* defaultname);
	
private:

	/* comment marker */
	int  fSkipComments;
	char fMarker;
	
	/* the filename */
	char* fFileName;
	
	/* NULL file name */
	static const char fNULLFileName;
};

/* inlines */
inline char ifstreamT::comment_marker(void) const { return fMarker; }
inline int ifstreamT::skip_comments(void) const { return fSkipComments; }
inline const char* ifstreamT::filename(void) const
{
	return (fFileName != NULL) ? fFileName : &fNULLFileName;
}

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
