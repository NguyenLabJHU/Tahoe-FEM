/* $Id: MSRMatrixT.h,v 1.6 2005-04-13 21:49:58 paklein Exp $ */
#ifndef _MSR_MATRIX_T_H_
#define _MSR_MATRIX_T_H_

/* base class */
#include "GlobalMatrixT.h"

/* direct members */
#include "AutoArrayT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "dMatrixT.h"
#include "nVariMatrixT.h"

namespace Tahoe {

/* forward declarations */
class MSRBuilderT;
class CommunicatorT;

/** Base class for matricies using the Modified Sparse Row (MSR) format. Methods beginning 
 * with names AZ_ are borrowed from Aztec 1.0. */
class MSRMatrixT: public GlobalMatrixT
{
public:

	/** constuctor */
	MSRMatrixT(ostream& out, int check_code, bool symmetric, const CommunicatorT& comm);

	/** copy constructor */
	MSRMatrixT(const MSRMatrixT& source);
	
	/** destructor */
	virtual ~MSRMatrixT(void);

	/** set the internal matrix structure. Initialization of the matrix requires two steps:
	 * -# set the matrix topology with MSRMatrixT::AddEquationSet
	 * -# call MSRMatrixT::Initialize after all equation sets have been added */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);

	/** \name set matrix structure
	 * active equations are eq > 0 */
	/*@{*/
	void AddEquationSet(const iArray2DT& eqnos);
	void AddEquationSet(const RaggedArray2DT<int>& eqnos);
	/*@}*/

	/** \name assembly operators
	 * active equations are eq > 0 */	
	/*@{*/
	virtual void Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos);
	virtual void Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
		const ArrayT<int>& col_eqnos);
	virtual void Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos);
	/*@}*/

	/** return the values along the diagonal of the matrix */
	virtual bool CopyDiagonal(dArrayT& diags) const;

	/** set all matrix values to 0.0 */
	void Clear(void);

	/** number scope */
	virtual EquationNumberScopeT EquationNumberScope(void) const { return kGlobal; };

	/** matrix prefers optimal ordering */
	virtual bool RenumberEquations(void) const { return false; };	

	/** assignment operator */
	MSRMatrixT& operator=(const MSRMatrixT& rhs);

protected:

	/** \name check functions */
	/*@{*/
	virtual void PrintAllPivots(void) const;
	virtual void PrintZeroPivots(void) const;
	virtual void PrintLHS(bool force = false) const;
	/*@}*/

	/** convert MSR data in MSRMatrixT::fMSRBuilder to RCV format
	 * \param drop_tol tolerance to drop values if absolute value is smaller
	 *        than the tolerance. Passing a negative value causes no
	 *        values to be dropped. */
	void GenerateRCV(iArrayT& r, iArrayT& c, dArrayT& v, double drop_tol) const;

private:

	/** assemble row values into global structure using the global column
	 * indices given in coldex - status is 1 if successful, 0 otherwise */
	void AssembleRow(int row, int numvals, const int* col_dex,
		const double* rowvals, int& status);

	/** assemble diagonal values into the global structure
	 * status is 1 if successful, 0 otherwise */
	void AssembleDiagonals(int numvals, const int* rows, const double* vals,
		int& status);

	/** sets MSR data with column indices sorted in ascending order */
	void SetMSRData(void);

	/** allocate memory for quick find */
	void SetUpQuickFind(void);

	/** \name borrowed from Aztec v1.0 */
	/*@{*/
	int AZ_quick_find(int key, int list[], int length, int shift, int bins[]);
	int AZ_find_index(int key, int list[], int length);	
	void AZ_init_quick_find(int list[], int length, int *shift, int *bins);
	void AZ_sort(int list[], int N, int list2[], double list3[]);
	/*@}*/

protected:

	/** true if only upper triangle of matrix is stored */
	bool fSymmetric;

	/** MSR database builder */
	MSRBuilderT* fMSRBuilder;

	/** \name matrix data in MSR format */
	/*@{*/
	iArrayT fupdate; /**< global indices updated on this processor */
	iArrayT fbindx;  /**< MSR structure data */
	dArrayT fval;    /**< matrix value array */
	/*@}*/

	/** \name quick find data */
	/*@{*/
	int     fQF_shift;
	iArrayT fupdate_bin;
	iArrayT fsrow_dex; /**< space for sorted row indices */
	dArrayT fsrow_val; /**< space for sorted row values */
	/*@}*/

	/** \name for assembly operations */
	/*@{*/
	AutoArrayT<int> fRowEqnVec, fColEqnVec;
	AutoArrayT<int> fRowDexVec, fColDexVec;
	AutoArrayT<double> fValVec;
	dMatrixT             fActiveBlk;
	nVariMatrixT<double> fActiveBlkMan;
	iArrayT fActiveDex; /**< symmetric assembly	only */
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _MSR_MATRIX_T_H_ */




