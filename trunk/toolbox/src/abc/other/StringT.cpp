/* $Id: StringT.cpp,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (08/01/1996)                                          */

#include "StringT.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>

#include "Constants.h"
#include "ExceptionCodes.h"

/* array behavior */
const bool ArrayT<StringT>::fByteCopy = false;
const bool ArrayT<StringT*>::fByteCopy = true;

/* line length */
const int kLineLength = 254;

/* constructor */
StringT::StringT(int length):
	ArrayT<char>(length)
{
	/* zero word length */
	if (length > 0)
		(*this)[0] = '\0';
	else
		*this = "\0";
}

/* input initializer */
istream& operator>>(istream& in, StringT& stringT)
{
	char string[kLineLength];
	in >> string;
	stringT = string;
	
	return(in);
}

ostream& operator<<(ostream& out, const StringT& string)
{
	out << (const char*) string.Pointer();

	return(out);
}

/* assignment operator */
StringT& StringT::operator=(const char* string)
{
	if (string != Pointer())
	{
		/* allocate memory */
		int length = strlen(string) + 1;
		Allocate(length);
	
		/* byte copy */
		memcpy(Pointer(), string, sizeof(char)*length);
	}

	return(*this);
}

/* copy what fits without resizing. new string length */
int StringT::CopyIn(const char* string)
{
	if (string == Pointer())
	{
		StringT tmp(string);
		return CopyIn(tmp);
	}
	else
	{
		int copy_length = strlen(string);
		if (copy_length >= fLength) copy_length = fLength - 1;
		if (copy_length > 0)
		{
			/* byte copy */
			memcpy(Pointer(), string, sizeof(char)*copy_length);
			fArray[copy_length] = '\0';
		}
		else if (copy_length < 0)
			copy_length = 0;
			
		return copy_length;
	}	
}

/* in/equality operators */
int StringT::operator==(const char* string) const
{
	return (strcmp(*this,string) == 0);
}

int operator==(const char* str_lhs, const StringT& str_rhs)
{
	return (str_rhs == str_lhs);
}

int StringT::operator!=(const char* string) const
{
	return (strcmp(*this,string) != 0);
}

int operator!=(const char* str_lhs, const StringT& str_rhs)
{
	return (str_rhs != str_lhs);
}

/* convert all to uppercase */
const StringT& StringT::ToUpper(void)
{
int num_chars = strlen(*this);
	char* p = Pointer();
	for (int i = 0; i < num_chars; i++)
	{
	    *p = toupper(*p);
		p++;
	}

	return *this;
}	

/* read a line from the input stream, where a line is the next
* kLineLength characters or fewer characters terminated
* by a newline */
void StringT::GetLineFromStream(istream& in)
{
	char string[kLineLength];
	in.getline(string, kLineLength-1);
	
	operator=(string);
}

/* append characters to the string */
StringT& StringT::Append(const char* s)
{
	if (s == Pointer())
	{
		StringT tmp(s);
		return Append(tmp);
	}
	else
	{
		int old_len = strlen(*this);
		int   s_len = strlen(s);
		int new_len = old_len + s_len + 1;

		/* reallocate (copy contents) */
		Resize(new_len, true);
	
		/* append new string (with '\0') */
		memcpy(Pointer(old_len), s, sizeof(char)*(s_len + 1));
		return *this;
	}
}

/* drop the last ".xxx" extension to the string */
StringT& StringT::Root(void)
{
	/* find last "." */
	int len = strlen(*this);
	char* p = fArray + len - 1;
	for (int i = 1; i < len && *p != '.'; i++) p--;
	
	/* drop tail */
	if (*p == '.')
	{
		int new_len = p - fArray + 1;
		Resize(new_len, true);
		fArray[new_len - 1] = '\0';
	}
	return *this;
}

StringT& StringT::Root(const char* s)
{
	/* check for self */
	if (s == Pointer())
		return Root();
	else
	{	
		/* find last "." */
		int len = strlen(s);
		const char* p = s + len - 1;
		for (int i = 0; i < len && *p != '.'; i++) p--;

		/* drop tail */
		if (*p == '.')
		{
			int new_len = p - s + 1;
			Allocate(new_len);
			memcpy(fArray, s, sizeof(char)*new_len);
			fArray[new_len - 1] = '\0';
		}	
		return *this;
	}
}

/* returns the last ".xxx" extension to the string */
StringT& StringT::Suffix(void)
{
	/* find last "." */
	int len = strlen(*this);
	char* p = fArray + len - 1;
	int i;
	for (i = 1; i < len && *p != '.'; i++) p--;
	
	/* take tail */
	if (*p == '.')
	{
		/* copy down */
		memmove(fArray, p, i + 1);

		/* shrink */	
		int new_len = strlen(fArray) + 1;
		Resize(new_len, true);
		
		return *this;
	}
	return *this = "\0";
}

StringT& StringT::Suffix(const char* s)
{
	/* check for self */
	if (s == Pointer())
		return Suffix();
	else
	{	
		/* find last "." */
		int len = strlen(s);
		const char* p = s + len - 1;
		for (int i = 0; i < len && *p != '.'; i++) p--;

		/* drop tail */
		if (*p == '.')
			return *this = p;
		else
			return *this = "\0";
	}
}

StringT& StringT::Append(const char* s1, const char* s2)
{
	int old_len = strlen(*this);
	int  s1_len = strlen(s1);
	int  s2_len = strlen(s2);
	int new_len = old_len + s1_len + s2_len + 1;

	/* reallocate (copy contents) */
	Resize(new_len, true);
	
	/* append new strings (terminate with '\0') */
	char* pstr = Pointer(old_len);
	memcpy(pstr, s1, sizeof(char)*s1_len);
	
	pstr += s1_len;
	memcpy(pstr, s2, sizeof(char)*(s2_len+1));
	
	return *this;
}
	
StringT& StringT::Append(const char* s1, const char* s2, const char* s3)
{
	int old_len = strlen(*this);
	int  s1_len = strlen(s1);
	int  s2_len = strlen(s2);
	int  s3_len = strlen(s3);
	int new_len = old_len + s1_len + s2_len + s3_len + 1;

	/* reallocate (copy contents) */
	Resize(new_len, true);
	
	/* append new strings (terminate with '\0') */
	char* pstr = Pointer(old_len);
	memcpy(pstr, s1, sizeof(char)*s1_len);
	
	pstr += s1_len;
	memcpy(pstr, s2, sizeof(char)*s2_len);

	pstr += s2_len;
	memcpy(pstr, s3, sizeof(char)*(s3_len+1));
	
	return *this;
}

StringT& StringT::Append(char c)
{
	char tmp[2] = {'\0','\0'};
	tmp[0] = c;
	return Append(tmp);
}

/* append an integer - width specifies the minimum number of digits
* that will be appended, padded by zeroes if number has fewer
* digits than width */
StringT& StringT::Append(int number, int width)
{
	/* convert to string */
	char num_str[] = "000000000";
	IntegerToString(number, num_str);

	if (strlen(num_str) < width)
	{
		/* set padding */
		char pad_str[] = "0000000000";
		int pad_len = width - strlen(num_str);
		if (pad_len > strlen(pad_str))
		{
			cout << "\n StringT::Append: padding limit: ";
			cout << strlen(pad_str) << endl;
			throw eGeneralFail;
		}		
		else
			pad_str[pad_len] = '\0';
		
		return Append(pad_str, num_str);	
	}
	else
		return Append(num_str);
}

StringT& StringT::Append(const char* s, int number, int width)
{
	/* convert to string */
	char num_str[] = "000000000";
	IntegerToString(number, num_str);

	if (strlen(num_str) < width)
	{
		/* set padding */
		char pad_str[] = "0000000000";
		int pad_len = width - strlen(num_str);
		if (pad_len > strlen(pad_str))
		{
			cout << "\n StringT::Append: padding limit: ";
			cout << strlen(pad_str) << endl;
			throw eGeneralFail;
		}		
		else
			pad_str[pad_len] = '\0';
		
		return Append(s, pad_str, num_str);	
	}
	else
		return Append(s, num_str);
}

/* insert characters at the beginning of the string */
StringT& StringT::Prepend(const char* s)
{
	if (s == Pointer())
	{
		StringT tmp(s);
		return Prepend(s);
	}
	else
	{
		int s_len = strlen(s);
		if (s_len > 0)
		{
			/* lengths */
			int old_len = strlen(*this);
			Resize(old_len + s_len + 1, true);
	
			/* shift back */
			char* str = *this;
			memmove(str + s_len, str, old_len + 1);
		
			/* copy in */
			memcpy(str, s, s_len);
		}
		return *this;
	}
}

StringT& StringT::Prepend(const char* s1, const char* s2)
{
	int s1_len = strlen(s1);
	int s2_len = strlen(s2);
	if (s1_len > 0 || s2_len > 0)
	{
		/* lengths */
		int old_len = strlen(*this);
		Resize(old_len + s1_len + s2_len + 1, true);

		/* shift back */
		char* str = *this;
		memmove(str + s1_len + s2_len, str, old_len + 1);
		
		/* copy s1 in*/
		memcpy(str, s1, s1_len);

		/* copy s2 in*/
		str += s1_len;
		memcpy(str, s2, s2_len);
	}
	return *this;
}

/* drop n characters from the string from the start (n > 0) or
* from the end (n < 0) */
StringT& StringT::Drop(int n)
{
	/* check */
	if ((int) fabs(n) > strlen(*this)) throw eOutOfRange;
	//NOTE - SUNWspro 5.0 doesn't like int(fabs(n))
	
	if (n > 0)
	{
		char* str = *this;
		memmove(str, str + n, strlen(str) + 1 - n);
	}
	else if (n < 0)
		(*this)[strlen(*this) + n] = '\0';

	return *this;
}

/* take n characters from the source from the start (n > 0) or
* from the end (n < 0) */
StringT& StringT::Take(const StringT& source, int n)
{
	/* self */
	if (source.Pointer() == Pointer())
	{
		StringT tmp(source);
		return Take(tmp, n);
	}
	else
	{
		/* check */
		int size = (n < 0) ? -n : n;
		if (n > strlen(source)) throw eOutOfRange;

		/* allocate */
		Allocate(size + 1);
		if (n > 0)
		{
			memcpy(Pointer(), source.Pointer(), size);
			(*this)[size] = '\0';
		}
		else if (n < 0)
		{
			const char* src = source.Pointer(strlen(source) - size);
			memcpy(Pointer(), src, size + 1);
		}

		return *this;
	}
}

StringT& StringT::Take(const StringT& source, int start, int end)
{
	if (source.Pointer() == Pointer())
	{
		StringT tmp(source);
		return Take(tmp, start, end);
	}
	else
	{
		/* checks */
		if (start < 0 ||
		    end > strlen(source) ||
		    end < start) throw eOutOfRange;

		/* allocate */
		int n = end - start + 1;
		Allocate(n + 1);
		if (n > 0) memcpy(Pointer(), source.Pointer(start), n);
		(*this)[n] = '\0';
		return *this;
	}
}

/* copy the first word from the source */
StringT& StringT::FirstWord(const StringT& source, int& count, bool C_word_only)
{
	if (source.Pointer() == Pointer())
	{
		StringT tmp(source);
		return FirstWord(tmp, count, C_word_only);
	}
	else
	{
		/* skip leading white space */
		int length = strlen(source);
		count = 0;
		char* str = source;
		while (count < length && isspace(*str))
		{
			str++;
			count++;
		}

		/* count word length */
		char* start = str;
		int word_count = 0;
		if (C_word_only)
		{
			while (count < length && (isalnum(*str) || *str == '_'))
			{
				str++;
				count++;
				word_count++;
			}
		}
		else
		{
			while (count < length && !isspace(*str))
			{
				str++;
				count++;
				word_count++;
			}
		}
		
		/* count trailing whitespace */
		while (count < length && isspace(*str))
		{
			str++;
			count++;
		}

		/* copy word */
		Allocate(word_count + 1);
		if (word_count > 0) memcpy(*this, start, word_count);
		(*this)[word_count] = '\0';
	
		return *this;
	}
}

/* drop leading white space */
StringT& StringT::DropLeadingSpace(void)
{
	/* count leading white space */
	char* str = Pointer();
	int count = 0;
	int length = strlen(str);
	while (count < length && isspace(*str++)) count++;

	/* remove leading white space */
	if (count > 0) Drop(count);
	return *this;
}

StringT& StringT::DropTrailingSpace(void)
{
	char* str = Pointer();
	int length = strlen(str);
	int dex = length - 1;
	while (dex > -1 && isspace(str[dex])) dex--;
		
	/* drop */	
	if (dex != length - 1) str[dex + 1] = '\0';
	return *this;
}

/* return a string with the extension and suffix tacked into
* the root of this.  The default extension for this is ".in" */
StringT& StringT::DefaultName(const StringT& sourcename, const char* extint,
	const char* extout, int suffix)
{
	char* source = sourcename.Pointer();
	char  defname[kLineLength];
	
	char  num[]     = "01234567890123456789";
	char  teststr[] = "01234567890123456789";
	
	int lext  = strlen(extint);
	int lroot = strlen(source) - lext;

	strncpy(teststr, &source[lroot], lext+1);
	
	if ( strcmp(extint,teststr) == 0 )
	{
		strncpy(defname, source, lroot);
		defname[lroot] = '\0';
	}
	else
		strcpy(defname, source);
	
	strcat(defname,extout);
	
	if (suffix > -1)
	{
		IntegerToString(suffix,num);
		strcat(defname, num);
	}
	
	/* drop in the result */
	operator=(defname);
	
	return *this;
}

/* convert string to native, relative file path */
void StringT::ToNativePathName(void)
{
#ifdef _MACOS_
	ToMacOSPath();
	return;
#endif

#ifdef _WINNT_
	ToWinNTPath();
	return;
#endif

#ifdef _UNIX__
	ToUNIXPath();
	return;
#endif

	/* fall through */
	cout << "\n StringT::ToNativePathName: unknown platform" << endl;
	throw eGeneralFail;
}

/* path string translators */
void StringT::ToMacOSPath(void)
{
	if (fArray[0] == ':')
		return;
	else
	{
		/* length check */
		if (strlen(fArray) >= kLineLength-1)
		{	
			cout << "\n StringT::ToMacOSPath: path too long" << endl;
			throw eGeneralFail;
		}

		char temp[kLineLength];
		
		int len = strlen(fArray);
		int i = 0;

		temp[i++] = ':';
		for (int j = 0; j < len; j++)
		{
			if (fArray[j] == '\\' || fArray[j] == '/')
				temp[i++] = ':';
			else if (fArray[j] == '.')
			{
				if (fArray[j+1] == '.')
				{
					temp[i++] = ':';
					j += 2;
				}	
				else
					temp[i++] = fArray[j];
			}
			else
				temp[i++] = fArray[j];
		}

		/* append term null character */
		temp[i] = '\0';

		/* copy in */
		*this = temp;
	}
}

void StringT::ToWinNTPath(void)
{
	if (fArray[0] == '\\')
		return;
	else
	{
		/* length check */
		if (strlen(fArray) >= kLineLength-1)
		{	
			cout << "\n StringT::ToWinNTPath: path too long" << endl;
			throw eGeneralFail;
		}

		/* driver */
		ToNTorUNIX('/', '\\');
	}
}

void StringT::ToUNIXPath(void)
{
	if (fArray[0] == '/')
		return;
	else
	{
		/* length check */
		if (strlen(fArray) >= kLineLength-1)
		{	
			cout << "\n StringT::ToUNIXPath: path too long" << endl;
			throw eGeneralFail;
		}

		/* driver */
		ToNTorUNIX('\\', '/');
	}
}

/* print ASCII codes */
void StringT::PrintCodes(ostream& out) const
{
	for (int i = 0; i < Length(); i++)
		out << setw(5) << int((*this)[i]);
	out << '\n';
}

/**********************************************************************
* Private
**********************************************************************/

/* returns the character string corresponding to the number */
void StringT::IntegerToString(int number, char* string) const
{
	if (number == 0)
	{
		/* check space */
		if (strlen(string) < 2) throw eGeneralFail;
		
		string[0] = '0';
		string[1] = '\0';
	}
	else
	{
		int lens = (number > 0) ? int(log10(double(number)) + 2) : 1;
			// extra space for '\0';
	
		/* check that string has enough space */
		if (strlen(string) + 1 < lens) throw eSizeMismatch;
		
		/* set all bytes to 0! */
		memset(string, '\0', sizeof(char)*lens);
			
		for (int i = 0; i < lens-1; i++)
		{
			int power = int(pow(10.0,lens-2-i));
			int digit = number/power;
	
			string[i] = '0' + digit;
			
			number -= digit*power;
		}
	}
}

void StringT::ToNTorUNIX(char from, char to)
{
	/* check */
	if ((from != '\\' && from != '/') ||
	     ( to != '\\' &&   to != '/')) throw eGeneralFail;

	char temp[kLineLength];

	int len = strlen(fArray);
	int i = 0;

	for (int j = 0; j < len; j++)
	{
		if (fArray[j] == from)
			temp[i++] = to;
		else if (fArray[j] == ':')	
		{
			if (j > 0) temp[i++] = to;
		
			if (fArray[j+1] == ':')
			{
				while (fArray[j+1] == ':')
				{
					temp[i++] = '.';
					temp[i++] = '.';
					temp[i++] = to;
					j++;						
				}
			}
		}
		else
			temp[i++] = fArray[j];
	}

	/* append term null character */
	temp[i] = '\0';
	
	/* copy in */
	*this = temp;
}
