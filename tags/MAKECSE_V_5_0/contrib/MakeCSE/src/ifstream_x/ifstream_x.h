/*
 * File: ifstream_x.h
 *
 * Input file stream which skips to the end of the line when the
 * comment marker is read.
 */

/*
 * created      : PAK (03/03/99)
 * last modified: PAK (03/03/99)
 */

#ifndef _IFSTREAM_X_H_
#define _IFSTREAM_X_H_

/* base class */
#include <fstream.h>
#include <stddef.h>

class ifstream_x
{
  public:

	/* constructors */
	ifstream_x(void);
	ifstream_x(const char* stream);
	ifstream_x(char marker);
	ifstream_x(char marker, const char* stream);

	/* destructor */
	virtual ~ifstream_x(void);

	/* type conversion operator */
	operator ifstream&();

	/* open stream */
	void open(const char* stream);
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
	
	/* return the filename - will be NULL if no file is open */
	const char* filename(void) const;
	void set_filename(const char* name); // force name

	/* overloaded istream functions */
	ifstream_x& putback(char c);
	ifstream_x& get(char& c);
	int good(void) const;
	void clear(void);

	/* adjusting stream position - returns actual number of rewound lines */
	int rewind(int num_lines = 1);

// Add next:
// (1) ofstream_x class with autoformatting
// (1.1) also has returns d_width and i_width, and other formatting
// (2) pull all file functions out of the FEManagerT
// (3) changes names to ifstreamT and ofstreamT
// (4) include both in single fstreamT .h/.cpp

	/* advance to the end of the line (or next 255 characters) */
	void clear_line(void);

	/* advances passed comments */
	void do_skip_comments(void);

	/* extraction of streams */
	ifstream_x& operator>>(bool& a);

  private:

	/* copy the string to fFileName */
	void CopyName(const char* filename);

	/* open stream with prompt - return 1 if successful */
	int OpenWithPrompt(const char* prompt, const char* skipname, 
		const char* defaultname);
	
  private:
  
  	/* the stream */
  	ifstream fifstream;
  
  	/* comment marker */
  	int  fSkipComments;
  	char fMarker;
  	
  	/* the filename */
  	char* fFileName;
};

/* inlines */
inline char ifstream_x::comment_marker(void) const { return fMarker; }
inline int ifstream_x::skip_comments(void) const { return fSkipComments; }
inline const char* ifstream_x::filename(void) const { return fFileName; }

inline void ifstream_x::set_marker(char marker)
{
	fSkipComments = 1;
	fMarker = marker;
}

inline void ifstream_x::clear_marker(void) { fSkipComments = 0; }
inline ifstream_x::operator ifstream&() { return  fifstream; }

/* overloaded istream functions */
inline ifstream_x& ifstream_x::putback(char c)
{
	fifstream.putback(c);
	return *this;
}

inline ifstream_x& ifstream_x::get(char& c)
{
	fifstream.get(c);
	return *this;
}

inline int ifstream_x::good(void) const { return fifstream.good(); }
inline void ifstream_x::clear(void) { fifstream.clear(); }

/* extraction operator */
template <class TYPE> ifstream_x& operator>>(ifstream_x& str, TYPE& data);
template <class TYPE> 
ifstream_x& operator>>(ifstream_x& str, TYPE& data)
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
