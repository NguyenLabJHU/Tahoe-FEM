/* $Id: iArrayT.cpp,v 1.3 2001-07-31 18:32:12 sawimme Exp $ */
/* created: paklein (08/10/1996)                                          */

#include "iArrayT.h"
#include <iostream.h>
#include <iomanip.h>
#include "Constants.h"

/* array behavior */
const bool ArrayT<iArrayT*>::fByteCopy = true;
const bool ArrayT<const iArrayT*>::fByteCopy = true;
const bool ArrayT<iArrayT>::fByteCopy = false; 

/* constructor */
iArrayT::iArrayT(void) { }
iArrayT::iArrayT(int length): nArrayT<int>(length) { }
iArrayT::iArrayT(int length, int* p): nArrayT<int>(length,p) { }
iArrayT::iArrayT(const iArrayT& source): nArrayT<int>(source) { }

/* flagging operations */
int iArrayT::ChangeValue(int from, int to)
{
	int count = 0;
	int* p = Pointer();
	for (int i = 0; i < Length(); i++)
	{
		int& ptemp = *p++;
		if (ptemp == from)
		{
			ptemp = to;
			count++;
		}
	}
	return count;
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
