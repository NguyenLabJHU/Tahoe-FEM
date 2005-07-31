/* $Id: StringT.h,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (08/01/1996)                                          */

#ifndef _STRING_T_H_
#define _STRING_T_H_

/* Environmental */
#include "Environment.h"

/* base class */
#include "ArrayT.h"

/* forward declarations */
#include "ios_fwd_decl.h"

class StringT: public ArrayT<char>
{
public:

	/* constructors */
	StringT(void);
	StringT(const StringT& string);
	StringT(const char* string);
	StringT(int length);
	
	/* type conversion operator - allows use of StringT in all const char*
	 * ANSI C functions. */
	operator const char*() const;
	operator char*() const;
	
	/* input initializer */
	friend istream& operator>>(istream& in, StringT& string);
	friend ostream& operator<<(ostream& out, const StringT& string);

	/* assignment operator */
	StringT& operator=(const char* string);
	StringT& operator=(const StringT& string);

	/* copy what fits into the current string length. returns new string length */
	int CopyIn(const char* string);

	/* in/equality operators */
	int operator==(const StringT& rhs) const;
	int operator==(const char* string) const;
	friend int operator==(const char* str_lhs, const StringT& str_rhs);

	int operator!=(const StringT& rhs) const;
	int operator!=(const char* string) const;
	friend int operator!=(const char* str_lhs, const StringT& str_rhs);

	/* convert all to uppercase */
	const StringT& ToUpper(void);

	/* read a line from the input stream, where a line is the next
	 * kFileNameLength characters or fewer characters terminated
	 * by a newline. */
	void GetLineFromStream(istream& in);

	/* drop the last ".xxx" extension to the string */
	StringT& Root(void);
	StringT& Root(const char* s);
	
	/* returns the last ".xxx" extension to the string */
	StringT& Suffix(void);
	StringT& Suffix(const char* s);

	/* append characters to the string */
	StringT& Append(const char* s);
	StringT& Append(const char* s1, const char* s2);
	StringT& Append(const char* s1, const char* s2, const char* s3);
	StringT& Append(char c);
	
	/* append an integer - width specifies the minimum number of digits
	 * that will be appended, padded by zeroes if number has fewer
	 * digits than width */
	StringT& Append(int number, int width = 0);
	StringT& Append(const char* s, int number, int width = 0);

	/* insert characters at the beginning of the string */
	StringT& Prepend(const char* s);
	StringT& Prepend(const char* s1, const char* s2);
	
	/* drop n characters from the string from the start (n > 0) or
	 * from the end (n < 0) */
	StringT& Drop(int n);
	
	/* take n characters from the source from the start (n > 0) or
	 * from the end (n < 0) */
	StringT& Take(const StringT& source, int n);
	StringT& Take(const StringT& source, int start, int end);
	
	/* copy the first word from the source and return number of characters scanned */
	StringT& FirstWord(const StringT& source, int& count, bool C_word_only);

	/* drop leading white space */
	StringT& DropLeadingSpace(void);
	StringT& DropTrailingSpace(void);
	
	/* return a string with the extension and suffix tacked into
	 * the root of this. The default extension for this is ".in" */
	StringT& DefaultName(const StringT& sourcename);
	StringT& DefaultName(const StringT& sourcename, const char* extout, int suffix = -1);
	StringT& DefaultName(const StringT& sourcename, const char* extint, const char* extout,
		int suffix);

	/* convert string to native, relative file path */
	void ToNativePathName(void);
	void ToMacOSPath(void);			
	void ToWinNTPath(void);			
	void ToUNIXPath(void);			

	/* print ASCII codes */
	void PrintCodes(ostream& out) const;

private:
	
	/* to NT or UNIX path - from and to are the delimiting characters */
	void ToNTorUNIX(char from, char to);

	/* deep copy of string */
	void CopyString(const char* string);
	
	/* returns the character string corresponding to the number */
	void IntegerToString(int number, char* string) const;
};

/* inlines */

/* constructors */
inline StringT::StringT(void) { operator=("\0"); }
inline StringT::StringT(const StringT& string): ArrayT<char>(string) { }
inline StringT::StringT(const char* string) { operator=(string); }

/* type conversion operator - allows use of StringT in all const char*
* ANSI C functions */
inline StringT::operator const char*() const { return Pointer(); }
inline StringT::operator char*() const { return Pointer(); }

/* assignment operator */
inline StringT& StringT::operator=(const StringT& string)
{
	return operator=(string.Pointer());
}

/* equality operator */
inline int StringT::operator==(const StringT& rhs) const
{
	return operator==(rhs.Pointer());
}

inline int StringT::operator!=(const StringT& rhs) const
{
	return operator!=(rhs.Pointer());
}

/* redirects */
inline StringT& StringT::DefaultName(const StringT& sourcename)
{
	return DefaultName(sourcename, ".in", ".out", -1);
}

inline StringT& StringT::DefaultName(const StringT& sourcename,
	const char* extout, int suffix)
{
	return DefaultName(sourcename, ".in", extout, suffix);
}

#endif /* _STRING_T_H_ */
