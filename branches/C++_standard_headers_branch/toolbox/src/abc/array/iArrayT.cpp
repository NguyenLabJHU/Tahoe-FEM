/* $Id: iArrayT.cpp,v 1.18 2003-11-21 22:41:30 paklein Exp $ */
/* created: paklein (08/10/1996) */
#include "iArrayT.h"
#include <iostream.h>
#include <iomanip.h>
#include "toolboxConstants.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<iArrayT*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<const iArrayT*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<iArrayT>::fByteCopy = false; 
} /* namespace Tahoe */

/* constructor */
iArrayT::iArrayT(void) { }
iArrayT::iArrayT(int length): nArrayT<int>(length) { }
iArrayT::iArrayT(int length, const int* p): nArrayT<int>(length,p) { }
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
	const int* p = Pointer();
	int  count = 0;
	for (int i = 0; i < Length(); i++)
		if (*p++ == value)
			count++;
	return count;			
}

int iArrayT::HasValue(int value) const
{
	const int* p = Pointer();
	for (int i = 0; i < Length(); i++)
		if (*p++ == value)
			return 1;
	return 0;			
}

int iArrayT::HasValue(int value, int& index) const
{
	index  = -1;
	const int* p = Pointer();
	for (int i = 0; i < Length(); i++)
		if (*p++ == value)
		{
			index = i;
			return 1;
		}
	return 0;			
}

/* determine union of the given array */
iArrayT& iArrayT::Union(const nArrayT<int>& source)
{
	/* quick exit */
	if (source.Length() == 0)
		Dimension(0);
	else
	{
		/* range of data */
		int min, max;
		source.MinMax(min, max);
		int range = max - min + 1;

		/* local map */
		iArrayT node_map(range);
		node_map = 0;

		/* mark nodes used */
		const int* p = source.Pointer();
		for (int i = 0; i < source.Length(); i++)
			node_map[*p++ - min] = 1;

		/* collect list */
		Dimension(node_map.Count(1));
		int dex = 0;
		p = node_map.Pointer();
		int* pthis = Pointer();
		for (int j = 0; j < node_map.Length(); j++)
			if (*p++ == 1) pthis[dex++] = j + min;
	}
	return *this;
}

/* determine the union of the given arrays */
iArrayT& iArrayT::Union(const ArrayT<const nArrayT<int>*>& source)
{
	/* quick exit */
	if (source.Length() == 0)
		Dimension(0);
	else
	{
		/* verify list and skip empties */
		iArrayT empty(source.Length());
		empty = 1;
		for (int i = 0; i < empty.Length(); i++)
		{
			const nArrayT<int>* a = source[i];
			if (!a) ExceptionT::GeneralFail("iArrayT::Union", "source array %d is NULL", i);
			if (a->Length() > 0) empty[i] = 0;
		}
	
		/* range of data */
		int count = empty.Count(0);
		if (count == 0)
			Dimension(0);
		else
		{
			iArrayT mins(count), maxs(count);
			int dex = 0;
			for (int i = 0; i < source.Length(); i++)
				if (!empty[i])
				{
					source[i]->MinMax(mins[dex], maxs[dex]);				
					dex++;
				}
			int min = mins.Min();
			int max = maxs.Max();
		
			/* node map */
			int range = max - min + 1;
			iArrayT node_map(range);
			node_map = 0;
	
			/* mark nodes used */
			for (int j = 0; j < source.Length(); j++)
				if (!empty[j]) 
				{
					const nArrayT<int>& src = *source[j];
					const int* p = src.Pointer();
					for (int i = 0; i < src.Length(); i++)
						node_map[*p++ - min] = 1;
				}

			/* collect list */
			Dimension(node_map.Count(1));
			dex = 0;
			int* p = node_map.Pointer();
			int* pthis = Pointer();
			for (int j = 0; j < node_map.Length(); j++)
				if (*p++ == 1) pthis[dex++] = j + min;
		}
	}

	return *this;
}
