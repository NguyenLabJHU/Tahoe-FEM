/* $Id: MSRBuilderT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (07/30/1998)                                          */
/* class to generate MSR matrix structure data                            */

#include "MSRBuilderT.h"
#include "Constants.h"
#include "ExceptionCodes.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "RaggedArray2DT.h"

/* constructor */
MSRBuilderT::MSRBuilderT(bool upper_only): fUpperOnly(upper_only) { }

/* return the MSR database and pointers to the start of each row */
void MSRBuilderT::SetMSRData(const iArrayT& activerows, iArrayT& MSRdata)
{
	/* set the graph data */
	MakeGraph(activerows, true, fUpperOnly);

	/* within range and in ascending order */
	CheckActiveSet(activerows);
	int numinactive = NumNodes() - activerows.Length();
	
	/* graph data */
	int row_shift;
	const RaggedArray2DT<int>& edgelist = EdgeList(row_shift);
	
	/* allocate space for MSR structure data */
	MSRdata.Allocate(edgelist.Length() + numinactive + 1);
	
	/* generate compressed MSR structure data */
	GenerateMSR(row_shift, edgelist, activerows, MSRdata);
}

/* active equations must be within range and in ascending order */
void MSRBuilderT::CheckActiveSet(const iArrayT& activerows) const
{
	/* check range of active rows */
	int min, max;
	activerows.MinMax(min,max);
	int num_nodes = NumNodes();
	if (min - fShift < 0 || max - fShift > num_nodes - 1)
	{
		cout << "\n MSRBuilderT::SetMSRData: active equations are out of range:";
		cout << "    {min,max} = {" << min << "," << max << "}"<< endl;
		throw eOutOfRange;
	}

	/* must be in ascending order */
	int* pactive = activerows.Pointer() + 1;
	for (int i = 1; i < activerows.Length(); i++)
	{
		if (*(pactive-1) >= *pactive)
		{
			cout << "\n MSRBuilderT::SetMSRData: active rows must be unique and";
			cout << " in ascending order." << endl;
			throw eGeneralFail;
		}
		pactive++;
	}
}

/* generate MSR structure */
void MSRBuilderT::GenerateMSR(int row_shift, const RaggedArray2DT<int>& edgelist,
	const iArrayT& activeeqs, iArrayT& MSRdata) const
{	
	/* dimensions */
	int numactive = activeeqs.Length();
	int MSRcolumnoffset = numactive + 1;

	int* pstarts = MSRdata.Pointer();
	int* pcol    = MSRdata.Pointer(MSRcolumnoffset);
	int* pactive = activeeqs.Pointer();
	
	*pstarts = MSRcolumnoffset;
	pstarts++; /* one column ahead */
	for (int i = 0; i < numactive; i++)
	{
		int  dex   = *pactive - row_shift;
		int  count = edgelist.MinorDim(dex);
		int* pdata = edgelist(dex);
		int  offdiagsize = count - 1;
		
		/* column start offsets */
		*pstarts = *(pstarts - 1) + offdiagsize;
		
		/* copy in off-diagonal data (skip self in 1st slot) */
		memcpy(pcol, pdata+1, offdiagsize*sizeof(int));
		
		/* sort column data */
		SortAscending(pcol, offdiagsize);
		
#if __option(extended_errorcheck)
		if (fUpperOnly)
			for (int j = 0; j < count; j++)
				if (pdata[j] < *pactive - 1) //OFFSET
				{
					cout << "\n MSRBuilderT::GenerateMSR: check of upper only fails on row: "
					     << *pactive << '\n';
					iArrayT tmp(count, pdata);
					cout << tmp.wrap(8) << endl;
					throw eGeneralFail;
				}
#endif

		/* next */
		pcol += offdiagsize;
		pactive++;
		pstarts++;
	}
}

/* write MSR data to output stream */
void MSRBuilderT::WriteMSRData(ostream& out, const iArrayT& activerows,
	const iArrayT& MSRdata) const
{
	iArrayT rowdata;
	int wrap = 12;
	out << " number of active rows: " << activerows.Length() << endl;
	for (int i = 0; i < activerows.Length(); i++)
	{
		out << " row: " << activerows[i] << '\n';
		int row_length = MSRdata[i+1] - MSRdata[i];
		if (row_length > 0)
		{
			rowdata.Set(MSRdata[i+1] - MSRdata[i], MSRdata.Pointer(MSRdata[i]));
			out << rowdata.wrap_tight(wrap) << '\n';
		}
	}
	out << '\n';
}

/* This routine was taken from Knuth: Sorting and Searching. It puts the input
* data list into a heap and then sorts it. */
void MSRBuilderT::SortAscending(int* list, int N) const
{
	int    l, r, RR, K, j, i, flag;

	if (N <= 1) return;

	l   = N / 2 + 1;
	r   = N - 1;
	l   = l - 1;
	RR  = list[l - 1];
	K   = list[l - 1];

	while (r != 0) {
		j = l;
		flag = 1;

		while (flag == 1) {
	i = j;
	j = j + j;

			if (j > r + 1)
				flag = 0;
			else {
				if (j < r + 1)
					if (list[j] > list[j - 1]) j = j + 1;

				if (list[j - 1] > K) {
					list[ i - 1] = list[ j - 1];
				}
				else {
					flag = 0;
				}
			}
		}

		list[ i - 1] = RR;

		if (l == 1) {
			RR  = list [r];

			K = list[r];
			list[r ] = list[0];
			r = r - 1;
		}
		else {
			l   = l - 1;
			RR  = list[ l - 1];
			K   = list[l - 1];
		}
	}

	list[ 0] = RR;
}
