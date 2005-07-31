/* $Id: SPOOLESMatrixT.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (09/13/2000)                                          */

#ifndef _SPOOLES_MATRIX_T_H_
#define _SPOOLES_MATRIX_T_H_

/* base class */
#include "GlobalMatrixT.h"

/* library support options */
#ifdef __SPOOLES__

/* direct members */
#include "AutoArrayT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "dMatrixT.h"
#include "nVariMatrixT.h"

/* forward declarations */
class MSRBuilderT;

class SPOOLESMatrixT: public GlobalMatrixT
{
public:

	/* constuctor */
	SPOOLESMatrixT(ostream& out, int check_code, bool symmetric,
		bool pivoting);
	
	/* destructor */
	virtual ~SPOOLESMatrixT(void);

	/* set the internal matrix structure.
	 * NOTE: do not call Initialize() equation topology has been set
	 * with AddEquationSet() for all equation sets */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);

	/* add to structure - active equations only (# > 0)
	 * NOTE: structure not set until #CALL#. equation data
	 * must persist outside of Aztec until (at least) then */
	void AddEquationSet(const iArray2DT& eqnos);
	void AddEquationSet(const RaggedArray2DT<int>& eqnos);

	/* assemble the element contribution into the LHS matrix - assumes
	 * that elMat is square (n x n) and that eqnos is also length n.
	 * NOTE: assembly positions (equation numbers) = 1...fNumEQ */
	virtual void Assemble(const ElementMatrixT& elMat, const iArrayT& eqnos);

	/* set all matrix values to 0.0 */
	void Clear(void);

	/* write non-zero values to stream as {row,col,value} */
	void PrintNonZero(ostream& out) const;

	/* number scope and reordering */
	virtual EquationNumberScopeT EquationNumberScope(void) const;
	virtual bool RenumberEquations(void) const;

	/* assignment operator - not implemented */
	virtual GlobalMatrixT& operator=(const GlobalMatrixT& RHS);

protected:

	/* precondition matrix */
	virtual void Factorize(void);
	
	/* determine new search direction and put the results in result */
	virtual void BackSubstitute(dArrayT& result);

	/* rank check functions */
	virtual void PrintAllPivots(void) const;
	virtual void PrintZeroPivots(void) const;
	virtual void PrintLHS(void) const;

	/* copy MSR data to RCV */
	void GenerateRCV(iArrayT& r, iArrayT& c, dArrayT& v);

private:

	/* assemble row values into global structure using the global column
	 * indices given in coldex - status is 1 if successful, 0 otherwise */
	void AssembleRow(int row, int numvals, const int* col_dex,
		const double* rowvals, int& status);

	/* assemble diagonal values into the global structure
	 * status is 1 if successful, 0 otherwise */
	void AssembleDiagonals(int numvals, const int* rows, const double* vals,
		int& status);

	/* sets MSR data with column indices sorted in ascending order */
	void SetMSRData(void);

	/* allocate memory for quick find */
	void SetUpQuickFind(void);

	/* borrowed from Aztec v1.0 */
	int AZ_quick_find(int key, int list[], int length, int shift, int bins[]);
	int AZ_find_index(int key, int list[], int length);	
	void AZ_init_quick_find(int list[], int length, int *shift, int *bins);
	void AZ_sort(int list[], int N, int list2[], double list3[]);

protected:

	/* SPOOLES parameters */
	bool fSymmetric;
	bool fPivoting;	

	/* MSR database builder */
	MSRBuilderT* fMSRBuilder;

	/* matrix data in MSR format */
	iArrayT fupdate; // global indices updated on this processor
	iArrayT fbindx;  // MSR structure data
	dArrayT fval;    // matrix value array

	/* quick find data */
	int     fQF_shift;
	iArrayT fupdate_bin;
	iArrayT fsrow_dex; // space for sorted row indices
	dArrayT fsrow_val; // space for sorted row values

	/* for assembly operations */
	AutoArrayT<int> fRowEqnVec, fColEqnVec;
	AutoArrayT<int> fRowDexVec, fColDexVec;
	AutoArrayT<double> fValVec;
	dMatrixT             fActiveBlk;
	nVariMatrixT<double> fActiveBlkMan;
	iArrayT fActiveDex; // symmetric assembly
	
	/* solve count */
	int fSolveCount;
};

#endif /*__SPOOLES__ */
#endif /* _SPOOLES_MATRIX_T_H_ */
