/* $Id: MSRBuilderT.cpp,v 1.8 2004-03-16 19:26:32 paklein Exp $ */
/* created: paklein (07/30/1998) */
#include "MSRBuilderT.h"
#include "toolboxConstants.h"
#include "ExceptionT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "RaggedArray2DT.h"

using namespace Tahoe;

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
	MSRdata.Dimension(edgelist.Length() + numinactive + 1);
	
	/* generate compressed MSR structure data */
	GenerateMSR(row_shift, edgelist, activerows, MSRdata);
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
			rowdata.Alias(MSRdata[i+1] - MSRdata[i], MSRdata.Pointer(MSRdata[i]));
			out << rowdata.wrap_tight(wrap) << '\n';
		}
	}
	out << '\n';
}

/* return the PSPASES data structure */
void MSRBuilderT::SetPSPASESData(const iArrayT& activerows, iArray2DT& aptrs, iArrayT& ainds)
{
	/* set the graph data - data for each row in activerows needs to be sequential */
	MakeGraph(activerows, true, fUpperOnly);

	/* within range and in ascending order */
	CheckActiveSet(activerows);
	int numinactive = NumNodes() - activerows.Length();
	
	/* graph data */
	int row_shift;
	const RaggedArray2DT<int>& edgelist = EdgeList(row_shift);
	
	/* dimension */
	aptrs.Dimension(activerows.Length(), 2);
	int ainds_dim = 0;
	for (int i = 0; i < activerows.Length(); i++)
		ainds_dim += edgelist.MinorDim(activerows[i] - row_shift);
	ainds.Dimension(ainds_dim);
	
	/* generate compressed MSR structure data */
	GeneratePSPASES(row_shift, edgelist, activerows, aptrs, ainds);
}

/* return the SuperLU data structure */
void MSRBuilderT::SetSuperLUData(const iArrayT& activerows, iArrayT& rowptr, iArrayT& colind)
{
	/* set the graph data */
	MakeGraph(activerows, true, fUpperOnly);

	/* within range and in ascending order */
	CheckActiveSet(activerows);
	int numinactive = NumNodes() - activerows.Length();
	
	/* graph data */
	int row_shift;
	const RaggedArray2DT<int>& edgelist = EdgeList(row_shift);

	/* dimension */
	rowptr.Dimension(activerows.Length() + 1); /* last entry points just beyond last row */
	int colind_dim = 0;
	for (int i = 0; i < activerows.Length(); i++)
		colind_dim += edgelist.MinorDim(activerows[i] - row_shift);
	colind.Dimension(colind_dim);

	/* generate compressed MSR structure data */
	GenerateSuperLU(row_shift, edgelist, activerows, rowptr, colind);
}

/************************************************************************
 * Private
 ************************************************************************/

/* active equations must be within range and in ascending order */
void MSRBuilderT::CheckActiveSet(const iArrayT& activerows) const
{
	const char caller[] = "MSRBuilderT::CheckActiveSet";

	/* check range of active rows */
	int min = 0, max = 0;
	if (activerows.Length() > 0) activerows.MinMax(min,max);
	int num_nodes = NumNodes();
	if (activerows.Length() > 0 && (min - fShift < 0 || max - fShift > num_nodes - 1))
		ExceptionT::OutOfRange(caller, "active equations out of range: {min, max} = {%d,%d}", min, max);

	/* must be in ascending order */
	const int* pactive = activerows.Pointer() + 1;
	for (int i = 1; i < activerows.Length(); i++)
	{
		if (*(pactive-1) >= *pactive)
			ExceptionT::GeneralFail(caller, "active rows not unique and ascending");

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
	const int* pactive = activeeqs.Pointer();
	
	*pstarts = MSRcolumnoffset;
	pstarts++; /* one column ahead */
	for (int i = 0; i < numactive; i++)
	{
		int  dex   = *pactive - row_shift;
		int  count = edgelist.MinorDim(dex);
		const int* pdata = edgelist(dex);
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
					throw ExceptionT::kGeneralFail;
				}
#endif

		/* next */
		pcol += offdiagsize;
		pactive++;
		pstarts++;
	}
}

/* generate PSPASES structure */
void MSRBuilderT::GeneratePSPASES(int row_shift, const RaggedArray2DT<int>& edgelist, const iArrayT& activeeqs, 
	iArray2DT& aptrs, iArrayT& ainds)
{
	/* dimensions */
	int numactive = activeeqs.Length();

	int index = 1;	
	const int* pactive = activeeqs.Pointer();	
	int* painds = ainds.Pointer();
	for (int i = 0; i < numactive; i++)
	{
		/* row information */	
		int  dex   = *pactive - row_shift;
		int  count = edgelist.MinorDim(dex);
		const int* pdata = edgelist(dex);

		/* set apts data */
		aptrs(i,0) = index;
		aptrs(i,1) = count;
			
		/* copy and sort column data */
		memcpy(painds, pdata, count*sizeof(int));
		SortAscending(painds, count);

		/* next row */
		index += count;
		painds += count;
		pactive++;
	}
}

/* generate SuperLU structure */
void MSRBuilderT::GenerateSuperLU(int row_shift, const RaggedArray2DT<int>& edgelist, const iArrayT& activeeqs,
	iArrayT& rowptr, iArrayT& colind)
{
	/* dimensions */
	int numactive = activeeqs.Length();

	/* save, sort, and count */
	const int* pactive = activeeqs.Pointer();	
	int* pcolind = colind.Pointer();
	rowptr[0] = 0;
	for (int i = 0; i < numactive; i++)
	{
		/* row information */	
		int dex = *pactive - row_shift;
		int count = edgelist.MinorDim(dex);
		const int* pdata = edgelist(dex);
			
		/* copy and sort column data */
		memcpy(pcolind, pdata, count*sizeof(int));
		SortAscending(pcolind, count);

		/* set rowptr data */
		rowptr[i+1] = rowptr[i] + count;

		/* next row */
		pcolind += count;
		pactive++;
	}
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
