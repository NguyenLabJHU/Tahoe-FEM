/* $Id: MSRBuilderT.h,v 1.6 2004-03-16 19:26:33 paklein Exp $ */
/* created: paklein (07/30/1998) */
#ifndef _MSR_BUILDER_T_H_
#define _MSR_BUILDER_T_H_

/* base class */
#include "GraphT.h"

namespace Tahoe {

/** class to generate MSR matrix structure data.
 * INPUT: equation topology
 * Equation topology can be described in sets which as added
 * using MSRBuilderT::AddGroup. Active equations are specified with
 * any numbers > 0. Equation numbers less than 1 are assumed to
 * be inactive and are ignored.
 * OUTPUT: MSR data
 * MSRBuilderT::SetMSRData generates the MSR structure data for the active rows
 * specified in the input. Active rows are numbered {0,...,max}. The
 * MSR data generated in also numbered with zero offset. */
class MSRBuilderT: public GraphT
{
public:

	/** constructor */
	MSRBuilderT(bool upper_only);

	/** \name MSR structure */
	/*@{*/
	/** return the MSR data structure */
	void SetMSRData(const iArrayT& activerows, iArrayT& MSRdata);

	/** write MSR data to output stream */
	void WriteMSRData(ostream& out, const iArrayT& activerows,
		const iArrayT& MSRdata) const;
	/*@}*/

	/** \name PSPASES (compressed row) structure */
	/*@{*/
	/** return the PSPASES data structure */
	void SetPSPASESData(const iArrayT& activerows, iArray2DT& aptrs, iArrayT& ainds);
	/*@}*/

	/** \name distributed SuperLU (compressed row) structure */
	/*@{*/
	/** return the distributed SuperLU data structure */
	void SetSuperLUData(const iArrayT& activerows, iArrayT& rowptr, iArrayT& colind);
	/*@}*/

private:

	/** active equations must be within range and in ascending order */
	void CheckActiveSet(const iArrayT& activeeqs) const;

	/** generate MSR structure */
	void GenerateMSR(int row_shift, const RaggedArray2DT<int>& edgelist,
		const iArrayT& activerows, iArrayT& MSRdata) const;

	/** generate PSPASES structure */
	void GeneratePSPASES(int row_shift, const RaggedArray2DT<int>& edgelist, const iArrayT& activerows, 
		iArray2DT& aptrs, iArrayT& ainds);

	/** generate SuperLU structure */
	void GenerateSuperLU(int row_shift, const RaggedArray2DT<int>& edgelist, const iArrayT& activerows,
		iArrayT& rowptr, iArrayT& colind);

	/** This routine was taken from Knuth: Sorting and Searching. It puts the input
	 * data list into a heap and then sorts it. */
	void SortAscending(int* list, int N) const;
	
private:
	
	/** true if upper triangle only. Passed in MSRBuilderT::MSRBuilderT. */
	bool fUpperOnly;
};

} // namespace Tahoe 
#endif /* _MSR_BUILDER_T_H_ */
