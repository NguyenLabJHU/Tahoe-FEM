/* $Id: iAutoArrayT.cpp,v 1.7 2003-10-27 19:50:33 paklein Exp $ */
/* created: paklein (02/08/1999) */
#include "iAutoArrayT.h"

#include <iostream.h>
#include <iomanip.h>

#include "toolboxConstants.h"

using namespace Tahoe;

namespace Tahoe {
template<> const bool ArrayT<iAutoArrayT*>::fByteCopy = true; 
template<> const bool ArrayT<iAutoArrayT>::fByteCopy = false; 
} /* namespace Tahoe */

/* max and min */
int iAutoArrayT::Max(void) const
{
	int* pthis = Pointer();
	int  max   = *pthis++;
	for (int i = 1; i < Length(); i++)
	{
		if (*pthis > max) max = *pthis;
		pthis++;
	}
	return max;
}

int iAutoArrayT::Min(void) const
{
	int* pthis = Pointer();
	int  min   = *pthis++;
	for (int i = 1; i < Length(); i++)
	{
		if (*pthis < min) min = *pthis;
		pthis++;
	}
	return min;
}

void iAutoArrayT::MinMax(int& min, int& max) const
{
	int* pthis = Pointer();
	min = max = *pthis++;
	for (int i = 1; i < Length(); i++)
	{
		if (*pthis < min)
			min = *pthis;
		else if (*pthis > max)
			max = *pthis;
		pthis++;
	}
}

/* flagging operations */
void iAutoArrayT::ChangeValue(int from, int to)
{
	int* pthis = Pointer();
	for (int i = 0; i < Length(); i++)
	{
		if (*pthis == from) *pthis = to;
		pthis++;
	}
}

int	iAutoArrayT::Count(int value) const
{
	int count = 0;
	int* pthis = Pointer();
	for (int i = 0; i < Length(); i++)
		if (*pthis++ == value) count++;
		
	return count;
}

namespace Tahoe {
ostream& operator<<(ostream& out, const iAutoArrayT& array)
{
	int* p = array.Pointer();
	for (int i = 0; i < array.Length(); i++)
		out << setw(kIntWidth) << *p++ << '\n';

	return out;
};
} /* namespace Tahoe */
