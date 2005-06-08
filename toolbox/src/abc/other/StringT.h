/* $Id: StringT.h,v 1.25 2005-06-08 17:20:20 paklein Exp $ */
/* created: paklein (08/01/1996) */
#ifndef _STRING_T_H_
#define _STRING_T_H_

/* Environmental */
#include "Environment.h"

/* base class */
#include "ArrayT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** string class */
class StringT: public ArrayT<char>
{
public:

	/** \name constructors */
	/*@{*/
	StringT(void);
	StringT(const StringT& string);
	StringT(const char* string);
	explicit StringT(int length);
	/*@}*/
	
	/** \name type conversion operator
	 * allows use of StringT in all const char* ANSI C functions. */
	/*@{*/
	operator char*() { return Pointer(); };
	operator const char*() const { return Pointer(); };
	/*@}*/
	
	/** input initializer */
	friend istream& operator>>(istream& in, StringT& string);
	friend ostream& operator<<(ostream& out, const StringT& string);

	/** \name assignment operators
	 * There's no operator=(char) because it leads to ambiguous conversion
	 * of int's and pointers. */
	/*@{*/
	void Fill(char a);
	StringT& operator=(const char* string);
	StringT& operator=(const StringT& string);
	/*@}*/

	/** copy what fits into the current string length. returns new string length */
	int CopyIn(const char* string);

	/** make string empty */
	void Clear(void);

	/** return the string length. the string length is the number of characters
	 * before the first '\0' character, i.e., the ANSI C strlen() function */
	int StringLength(void) const;

	/** \name equality operators */
	/*@{*/
	int operator==(const StringT& rhs) const;
	int operator==(const char* string) const;
	int operator==(char* string) const { return operator==((const char *) string); };
	friend int operator==(const char* str_lhs, const StringT& str_rhs);
	/*@}*/

	/** string match. Return a pointer to the first occurrence of the search string in
	 * this or NULL if the string is not found */
	const char* StringMatch(const char* search) const;

	/** \name comparisons */
	/*@{*/
	/** return true of {*this, rhs} are in alphabetical order, 0 otherwise. Note that
	 * the order is case sensitive [A-Z] < [a-z], i.e., capital letters come before
	 * lower case. */
	int operator<(const StringT& rhs) const;
	
	/** return true of {*this, rhs} are in reverse alphabetical order, 0 otherwise. Note that
	 * the order is case sensitive [A-Z] < [a-z], i.e., capital letters come before
	 * lower case. */
	int operator>(const StringT& rhs) const;
	/*@}*/

	/** \name inequality operators */
	/*@{*/
	int operator!=(const StringT& rhs) const;
	int operator!=(const char* string) const;
	int operator!=(char* string) const { return operator!=((const char*) string); };
	friend int operator!=(const char* str_lhs, const StringT& str_rhs);
	/*@}*/

	/** \name reformatting */
	/*@{*/
	/** convert all to uppercase */
	const StringT& ToUpper(void);

	/** convert all to lowercase */
	const StringT& ToLower(void);

	/** replace one character with another */
	void Replace(char find, char replace);
	/*@}*/

	/** read a line from the input stream, where a line is the next
	 * kFileNameLength characters or fewer characters terminated
	 * by a newline. */
	/*@{*/
	void GetLineFromStream(istream& in);
	void GetLineFromStream(ifstreamT& in);
	/*@}*/

	/** \name drop the last ".xxx" extension to the string */
	/*@{*/
	StringT& Root(char marker = '.');
	StringT& Root(const char* s, char marker = '.');
	/*@}*/
	
	/** \name returns the last ".xxx" extension to the string */
	/*@{*/
	StringT& Suffix(char marker = '.');
	StringT& Suffix(const char* s, char marker = '.');
	/*@}*/

	/** returns the path part of the full path to a file - drops the file
	 * from the full path to a file, keeping the directory separator */
	/*@{*/
	StringT& FilePath(void);
	StringT& FilePath(const char* s);
	/*@}*/

	/** \name append characters to the string */
	/*@{*/
	StringT& Append(const char* s);
	StringT& Append(const char* s1, const char* s2);
	StringT& Append(const char* s1, const char* s2, const char* s3);
	StringT& Append(char c);
	/*@}*/
	
	/** \name append an integer
	 * width specifies the minimum number of digits
	 * that will be appended, padded by zeroes if number has fewer
	 * digits than width */
	/*@{*/
	StringT& Append(int number, int width = 0);
	StringT& Append(const char* s, int number, int width = 0);
	/*@}*/

	/** append a floating point number */
	StringT& Append(double number, int precision = 2);

	/** insert characters at the beginning of the string */
	/*@{*/
	StringT& Prepend(const char* s);
	StringT& Prepend(const char* s1, const char* s2);
	/*@}*/
	
	/** drop n characters from the string from the start (n > 0) or
	 * from the end (n < 0) */
	StringT& Drop(int n);

	/** delete characters from the string from start to end, inclusive. */
	StringT& Delete(int start, int end);

	/** delete the specified character from the string. */
	StringT& Delete(int position) { return Delete(position, position); };
	
	/** take n characters from the source from the start (n > 0) or
	 * from the end (n < 0) */
	StringT& Take(const StringT& source, int n);

	/** copy a section of the source string */
	StringT& Take(const StringT& source, int start, int end);
	
	/** copy the first word from the source and return number of characters scanned.
	 * There are multiple rules available for determining word boundaries.
	 * \param source string from which to extract the first word
	 * \param count returns with the number of characters scanned from source. This
	 *        includes any leading white space, or non-alphanumeric characters that
	 *        were skipped. The number can be used to with StringT::Drop to remove 
	 *        the data associated with the first word from source.
	 * \param C_word_only if true, uses C variable naming rules, something like [a-zA-Z0-9_]+,
	 *        to determine word boundaries. If false, uses white space to determine word
	 *        boundaries or returns an entire quoted string, containing any collection of
	 *        characters. */
	StringT& FirstWord(const StringT& source, int& count, bool C_word_only);

	/** drop leading white space */
	StringT& DropLeadingSpace(void);

	/** drop trailing white space */
	StringT& DropTrailingSpace(void);

	/** \name convert string to native, relative file path */
	/*@{*/
	void ToNativePathName(void);
	void ToMacOSPath(void);			
	void ToWinNTPath(void);			
	void ToUNIXPath(void);
	static char DirectorySeparator(void);	
	/*@}*/

	/** print ASCII codes */
	void PrintCodes(ostream& out) const;

	/** version number comparison - returns 0 if the versions numbers are
	 * the same, -1 if v1 is older than v2, 1 if v1 is newer than v2 */
	static int versioncmp(const char* v1, const char* v2);
	
	/** extract double. perform type conversion on the tail of the string
	 * \param key last character before the start of the tail
	 * \param value conversion of tail, returns 0.0 if key not found
	 * \return true if key found, else returns false */
	bool Tail(char key, double& value) const;

	/** extract integer. perform type conversion on the tail of the string
	 * \param key last character before the start of the tail
	 * \param value conversion of tail, returns 0 if key not found
	 * \return true if key found, else returns false */
	bool Tail(char key, int& value) const;

	/** extract StringT. perform type conversion on the tail of the string
	 * \param key last character before the start of the tail
	 * \param value conversion of tail, returns 0 if key not found. String
	 *        will not contain any leading or trailing white space
	 * \return true if key found, else returns false */
	bool Tail(char key, StringT& value) const;

	/** \name scan for character
	 * return the position of the first or last occurrence of the character
	 * in the string. Returns -1 if the character is not present */
	/*@{*/
	int FirstPositionOf(char a) const;
	int LastPositionOf(char a) const;
	/*@}*/

	/** \name check for valid XML name */
	/*@{*/
	static bool IsXMLName(const char* s);
	bool IsXMLName(void) const { return IsXMLName(*this); };
	/*@}*/

private:
	
	/** to NT or UNIX path - from and to are the delimiting characters */
	void ToNTorUNIX(char from, char to);

	/** deep copy of string */
	void CopyString(const char* string);
	
	/** returns the character string corresponding to the number */
	void IntegerToString(int number, char* string) const;

private:

	/** character for empty strings */
	static char fEmpty;
};

/* inlines */

/* constructors */
inline StringT::StringT(void) { Alias(1, &fEmpty); }
inline StringT::StringT(const StringT& string): ArrayT<char>(string) { }
inline StringT::StringT(const char* string) { operator=(string); }

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

/* string length */
inline int StringT::StringLength(void) const { return strlen(*this); }

inline char StringT::DirectorySeparator(void)
{
#ifdef _MACOS_
#ifdef __MACH__
	return '/';
#else
	return ':';
#endif
#elif defined(_WINNT_)
	return '\\';
#endif

	/* UNIX by default */
	return '/';
}
}//namespace Tahoe
#endif /* _STRING_T_H_ */
