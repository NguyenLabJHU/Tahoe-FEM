/*
 * File: iArrayT.cpp
 *
 */

/*
 * created      : PAK (08/10/96)
 * last modified: PAK (08/20/97)
 */

#include <iostream.h>
#include <iomanip.h>

#include "Constants.h"
#include "iArrayT.h"

/* constructor */ 
iArrayT::iArrayT(void) { }
iArrayT::iArrayT(int length): nArrayT<int>(length) { }
iArrayT::iArrayT(int length, int* p): nArrayT<int>(length,p) { }
iArrayT::iArrayT(const iArrayT& source): nArrayT<int>(source) { }

/* flagging operations */
void iArrayT::ChangeValue(int from, int to)
{
	int* p = Pointer();

	for (int i = 0; i < Length(); i++)
	{
		int& ptemp = *p++;
	
		if (ptemp == from)
			ptemp = to;
	}
}

void iArrayT::PrintValued(ostream& out, int value, int wrapat) const
{
	int* p = Pointer();
	int  count = 0;

	for (int i = 0; i < Length(); i++)
		if (*p++ == value)
		{
			out << setw(kIntWidth) << i + 1;
			if (++count == wrapat)
			{
				out << '\n';
				count = 0;
			}	
		}
		
	/* final wrap */		
	if (count != 0)	out << '\n';
}

int iArrayT::Count(int value) const
{
	int* p = Pointer();
	int  count = 0;

	for (int i = 0; i < Length(); i++)
		if (*p++ == value)
			count++;

	return count;			
}

int iArrayT::HasValue(int value) const
{
	int* p = Pointer();

	for (int i = 0; i < Length(); i++)
		if (*p++ == value)
			return 1;

	return 0;			
}

int iArrayT::HasValue(int value, int& index) const
{
	index  = -1;
	int* p = Pointer();

	for (int i = 0; i < Length(); i++)
		if (*p++ == value)
		{
			index = i;
			return 1;
		}

	return 0;			
}

/* set array value to its position in the array */
void iArrayT::SetValueToPosition(void)
{
	int* p = Pointer();
	for (int i = 0; i < Length(); i++)
		*p++ = i;
}
