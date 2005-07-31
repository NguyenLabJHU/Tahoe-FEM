/* $Id: MSRBuilderT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (07/30/1998)                                          */
/* class to generate MSR matrix structure data                            */
/* INPUT: equation topology                                               */
/* Equation topology can be described in sets which as added              */
/* using AddEquationSet(). Active equations are specified with            */
/* any numbers > 0. Equation numbers less than 1 are assumed to           */
/* be inactive and are ignored.                                           */
/* OUTPUT: MSR data                                                       */
/* SetMSRData() generates the MSR structure data for the active rows      */
/* specified in the input. Active rows are numbered {0,...,max}. The      */
/* MSR data generated in also numbered with zero offset.                  */

#ifndef _MSR_BUILDER_T_H_
#define _MSR_BUILDER_T_H_

/* base class */
#include "GraphT.h"

class MSRBuilderT: public GraphT
{
public:

	/* constructor */
	MSRBuilderT(bool upper_only);
	
	/* return the MSR database and pointers to the start of
	 * each row */
	void SetMSRData(const iArrayT& activerows, iArrayT& MSRdata);

	/* write MSR data to output stream */
	void WriteMSRData(ostream& out, const iArrayT& activerows,
		const iArrayT& MSRdata) const;
	
private:

	/* active equations must be within range and in ascending order */
	void CheckActiveSet(const iArrayT& activeeqs) const;

	/* generate MSR structure */
	void GenerateMSR(int row_shift, const RaggedArray2DT<int>& edgelist,
		const iArrayT& activerows, iArrayT& MSRdata) const;

	/* This routine was taken from Knuth: Sorting and Searching. It puts the input
	 * data list into a heap and then sorts it. */
	void SortAscending(int* list, int N) const;
	
private:
	
	/* upper triangle only */
	bool fUpperOnly;
};

#endif /* _MSR_BUILDER_T_H_ */
