/* $Id: StringT.cpp,v 1.24 2002-09-03 07:06:06 paklein Exp $ */
/* created: paklein (08/01/1996) */

#include "StringT.h"
#include "ifstreamT.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#ifdef _MSC_VER
#include <strstrea.h>
#else
#include <strstream.h>
#endif

/* array behavior */

using namespace Tahoe;

namespace Tahoe {
const bool ArrayT<StringT>::fByteCopy = false;
const bool ArrayT<StringT*>::fByteCopy = true;
const bool ArrayT<const StringT*>::fByteCopy = true;
} /* namespace Tahoe */

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

namespace Tahoe {

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

} // namespace Tahoe

/* assignment operator */
StringT& StringT::operator=(const char* string)
{
	if (string == NULL)
		Clear();
	else if (string != Pointer())
	{
		/* allocate memory */
		int length = strlen(string) + 1;
		Dimension(length);
	
		/* byte copy */
		memcpy(Pointer(), string, sizeof(char)*length);
	}

	return(*this);
}

/* make string empty */
void StringT::Clear(void)
{
	/* should have at least one character */
	if (Length() < 1) throw eGeneralFail;
	
	/* zero string length */
	fArray[0] = '\0';
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

namespace Tahoe {
int operator==(const char* str_lhs, const StringT& str_rhs)
{
	return (str_rhs == str_lhs);
}
} /* namespace Tahoe */

int StringT::operator!=(const char* string) const
{
	return (strcmp(*this,string) != 0);
}

namespace Tahoe {
int operator!=(const char* str_lhs, const StringT& str_rhs)
{
	return (str_rhs != str_lhs);
}
} /* namespace Tahoe */

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

void StringT::GetLineFromStream(ifstreamT& in)
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

/* append a floating point number */
StringT& StringT::Append(double number, int precision)
{
	char buffer[51];
	memset(buffer, '\0', 51);
	ostrstream out(buffer, 50);
	out.precision(precision);
	out.setf(ios::showpoint);
	out.setf(ios::right, ios::adjustfield);
	out.setf(ios::scientific, ios::floatfield);

	out << number;

	Append(buffer);
	return *this;
}

/* drop the last ".xxx" extension to the string */
StringT& StringT::Root(char marker)
{
	/* find last instance or marker */
	int len = strlen(*this);
	char* p = fArray + len - 1;
	for (int i = 1; i < len && *p != marker; i++) p--;
	
	/* drop tail */
	if (*p == marker)
	{
		int new_len = p - fArray + 1;
		Resize(new_len, true);
		fArray[new_len - 1] = '\0';
	}
	return *this;
}

StringT& StringT::Root(const char* s, char marker)
{
	/* check for self */
	if (s == Pointer())
		return Root(marker);
	else
	{	
		/* find last marker */
		int len = strlen(s);
		const char* p = s + len - 1;
		for (int i = 0; i < len && *p != marker; i++) p--;

		/* drop tail */
		if (*p == marker)
		{
			int new_len = p - s + 1;
			Dimension(new_len);
			memcpy(fArray, s, sizeof(char)*new_len);
			fArray[new_len - 1] = '\0';
		} /* keep whole string */
		else
			*this = s;
		return *this;
	}
}

/* returns the last ".xxx" extension to the string */
StringT& StringT::Suffix(char marker)
{
	/* find last instance of marker */
	int len = strlen(*this);
	char* p = fArray + len - 1;
	int i;
	for (i = 1; i < len && *p != marker; i++) p--;
	
	/* take tail */
	if (*p == marker)
	{
		/* copy down */
		memmove(fArray, p, i + 1);

		/* shrink */	
		int new_len = strlen(fArray) + 1;
		Resize(new_len, true);
	} 
	else /* empty */
		*this = "\0";
	return *this;
}

StringT& StringT::Suffix(const char* s, char marker)
{
	/* check for self */
	if (s == Pointer())
		return Suffix(marker);
	else
	{	
		/* find last instance of marker */
		int len = strlen(s);
		const char* p = s + len - 1;
		for (int i = 0; i < len && *p != marker; i++) p--;

		/* drop tail */
		if (*p == marker)
			return *this = p;
		else
			return *this = "\0";
	}
}

/* returns the path part of the full path to a file - drops the file
 * from the full path to a file, keeping the directory separator */
StringT& StringT::FilePath(void)
{
	/* convert to native file path */
	ToNativePathName();

	/* directory separator */
	char separator = DirectorySeparator();

	/* find last instance of separator */
	int len = strlen(*this);
	char* p = fArray + len - 1;
	int i;
	for (i = 1; i < len && *p != separator; i++) p--;

	/* take tail including separator (unless MacOS) */
	if (*p == separator)
	{
#if defined(_MACOS_) && !defined(__MACH__)
		int offset = 1;
#else
		int offset = 2;		
#endif
		int new_len = p - fArray + offset;
		Resize(new_len, true);
		fArray[new_len - 1] = '\0';
	}
	/* no path */
	else
		fArray[0] = '\0';

	return *this;
}

StringT& StringT::FilePath(const char* s)
{
	/* copy */
	if (s != Pointer()) *this = s;
	
	/* drop file name */
	return FilePath();
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
	if ((int) fabs(double(n)) > strlen(*this)) throw eOutOfRange;
	//NOTE - SUNWspro 5.0 doesn't like int(fabs(n))
	
	if (n > 0)
	{
		char* str = *this;
		memmove(str, str + n, strlen(str) + 1 - n);
	}
	else if (n < 0)
		(*this)[int(strlen(*this)) + n] = '\0';

	return *this;
}

/* delete characters from the string from start to end inclusive */
StringT& StringT::Delete(int start, int end)
{
#if __option(extended_errorcheck)
	if (end < start || start < 0 || end >= strlen(*this))
		throw eOutOfRange;
#endif

	char* str = *this;
	memmove(str + start, str + end + 1, strlen(str) - end);

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
		Dimension(size + 1);
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
		Dimension(n + 1);
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
			/* return quoted strings as one word */
			if (*str == '"')
			{
				/* skip forward */
				str++;
				start++;
				count++;
				
				/* bound whole quoted string */
				while (count < length && *str != '"')
				{
					str++;
					count++;
					word_count++;
				}
				
				/* skip trailing quote (if present) */
				if (*str == '"')
				{
					str++;
					count++;
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
		}
		
		/* count trailing whitespace */
		while (count < length && isspace(*str))
		{
			str++;
			count++;
		}

		/* copy word */
		Dimension(word_count + 1);
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

/* convert string to native, relative file path */
void StringT::ToNativePathName(void)
{
#ifdef _MACOS_
#ifdef __MACH__
	ToUNIXPath();
#else
	ToMacOSPath();
#endif
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
		int j_start = 0;

		/* relative path */
		if (fArray[0] != '/' && fArray[0] != '\\')
		{
			temp[i++] = ':';
			
			/* forward */
			if (fArray[0] == '.' && fArray[1] == '/')
				j_start = 2;
		}
		/* forward */
		else
			while (fArray[j_start] == '/' || fArray[j_start] == '\\')
				j_start++;

		for (int j = j_start; j < len; j++)
		{
			if (fArray[j] == '\\' || fArray[j] == '/')
			{
				temp[i++] = ':';
				
				/* forward */
				while (fArray[j+1] == '/' || fArray[j+1] == '\\') j++;
			}
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
	if (fArray[0] == '\\') /* native path starting with root */
		return;
	else if (fArray[1] == ':' && 
	         (fArray[2] == '\\' || fArray[2] == '/')) /* native path starting with drive */
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

/* version number comparison - returns 0 if the versions numbers are
 * the same, -1 if v1 is older than v2, 1 if v1 is newer than v2 */
int StringT::versioncmp(const char* v1_, const char* v2_)
{
	/* skip leading char's */
	while (!isdigit(*v1_) && !isdigit(*v2_))
	{
		if (*v1_ != *v2_)
		{
			cout << "\n StringT::versioncmp: incompatible version numbers:\n" 
			     << '\t' << v1_ << '\n'
			     << '\t' << v2_ << endl;
			throw eGeneralFail;
		}
		v1_++; v2_++;
	}

	/* check */
	int l1 = strlen(v1_);
	int l2 = strlen(v2_);	
	if (l1 > 50 || l2 > 50)
	{
		cout << "StringT::versioncmp: exceeded maximum version string length: 50" << endl;
		throw eGeneralFail;
	}

	/* copy in and pad with trailing space:
	 * SGI, GNU, ?: hit EOF after reading last value
	 *  DEC, CW, ?: hit EOF when reading last value */
	char v1[52];
	memcpy(v1, v1_, l1);
	v1[l1] = '#'; // my EOF
	v1[l1+1] = '\0';
	char v2[52];
	memcpy(v2, v2_, l2);
	v2[l2] = '#'; // my EOF
	v2[l2+1] = '\0';

	istrstream s1(v1), s2(v2);
	while (true)
	{
		int i1 = -9999, i2 = -9999;
		s1 >> i1;
		s2 >> i2;
		
		/* error reading */
		if (i1 == -9999)
		{
			cout << "\n StringT::versioncmp: error reading version number: "
			     << v1 << endl;
			throw eGeneralFail;
		}
		if (i2 == -9999)
		{
			cout << "\n StringT::versioncmp: error reading version number: "
			     << v2 << endl;
			throw eGeneralFail;
		}
		
		/* resolve */
		if (i1 > i2)
			return 1;
		else if (i2 > i1)
			return -1;
		else /* the same so far */
		{
			/* get next char */
			char a1, a2;
			s1.get(a1);
			s2.get(a2);
			
			/* same */
			if (a1 == '#' && a2 == '#')
				return 0;
			else if (a1 == '#')
				return -1;
			else if (a2 == '#')
				return 1;
			/* check */
			else if (a1 != '.' || a2 != '.')
			{
				cout << "\n StringT::versioncmp: illegal character" << endl;
				throw eGeneralFail;
			}
		}
	}
}

/* extract double */
bool StringT::Tail(char key, double& value) const
{
	char* p = Pointer();
	while (*p != '\0' && *p != key) p++;

	value = 0.0;
	if (*p == key)
	{
		istrstream in(p + 1);
		in >> value;
		return true;
	}
	return false;
}

/* extract integer */
bool StringT::Tail(char key, int& value) const
{
	char* p = Pointer();
	while (*p != '\0' && *p != key) p++;

	value = 0;
	if (*p == key)
	{
		istrstream in(p + 1);
		in >> value;
		return true;
		
	}
	return false;
}

/* extract string */
bool StringT::Tail(char key, StringT& value) const
{
	char* p = Pointer();
	while (*p != '\0' && *p != key) p++;

	value = 0;
	if (*p == key)
	{
		istrstream in(p + 1);
		in >> value;
		return true;
		
	}
	return false;
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