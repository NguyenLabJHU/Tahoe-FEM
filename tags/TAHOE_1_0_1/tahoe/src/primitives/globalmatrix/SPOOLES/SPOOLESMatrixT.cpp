/* $Id: SPOOLESMatrixT.cpp,v 1.12 2002-11-25 07:13:41 paklein Exp $ */
/* created: paklein (09/13/2000) */

#include "SPOOLESMatrixT.h"

/* library support options */
#ifdef __SPOOLES__

#include "MSRBuilderT.h"
#include "ElementMatrixT.h"
#include "SPOOLES.h"

#ifdef __SPOOLES_MT__
#include "SPOOLESMT.h"
#endif

/* message file name */

using namespace Tahoe;

const char SPOOLES_FILE[] = "SPOOLES.out";

/* constuctor */
SPOOLESMatrixT::SPOOLESMatrixT(ostream& out, int check_code,
	bool symmetric, bool pivoting):
	GlobalMatrixT(out, check_code),
	fSymmetric(symmetric),
	fPivoting(pivoting),
	fActiveBlkMan(0, fActiveBlk),
	fSolveCount(0)	
{
	fMSRBuilder = new MSRBuilderT(fSymmetric);
	if (!fMSRBuilder) throw ExceptionT::kOutOfMemory;
}

SPOOLESMatrixT::SPOOLESMatrixT(const SPOOLESMatrixT& source):
	GlobalMatrixT(source)
{
	cout << "\n SPOOLESMatrixT::SPOOLESMatrixT: not implemented" << endl;
	throw ExceptionT::kGeneralFail;
}

/* destructor */
SPOOLESMatrixT::~SPOOLESMatrixT(void)
{
	delete fMSRBuilder;
	fMSRBuilder = NULL;
}

/* set the internal matrix structure.
* NOTE: do not call Initialize() equation topology has been set
* with AddEquationSet() for all equation sets */
void SPOOLESMatrixT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* inherited */
	GlobalMatrixT::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* set matrix structure */
	SetMSRData();

	/* set up quick row finding */	
	SetUpQuickFind();
}

/* add to structure - active equations only (# > 0)
* NOTE: structure not set until #CALL#. equation data
* must persist outside of Aztec until (at least) then */
void SPOOLESMatrixT::AddEquationSet(const iArray2DT& eqnos)
{
	/* send to MSR data builder */
	fMSRBuilder->AddGroup(eqnos);
}

void SPOOLESMatrixT::AddEquationSet(const RaggedArray2DT<int>& eqnos)
{
	/* send to MSR data builder */
	fMSRBuilder->AddGroup(eqnos);
}

/* assemble the element contribution into the LHS matrix - assumes
* that elMat is square (n x n) and that eqnos is also length n.
*
* NOTE: assembly positions (equation numbers) = 1...fNumEQ */
void SPOOLESMatrixT::Assemble(const ElementMatrixT& elMat, const nArrayT<int>& eqnos)
{
#if __option (extended_errorcheck)
	if (elMat.Rows() != elMat.Cols()) throw ExceptionT::kGeneralFail;
	if (elMat.Rows() != eqnos.Length()) throw ExceptionT::kSizeMismatch;
#endif

	/* element matrix format */
	ElementMatrixT::FormatT format = elMat.Format();
	if (format == ElementMatrixT::kDiagonal)
	{
		/* extract values for active equation numbers */
		fRowDexVec.Dimension(0);	
		fValVec.Dimension(0);
		int end_update = fStartEQ + fLocNumEQ - 1;
		for (int i = 0; i < eqnos.Length(); i++)
		{
			int eq = eqnos[i];
//			if (eq > 0)
			if (eq >= fStartEQ && eq <= end_update)
			{
				fRowDexVec.Append(eq - 1); //OFFSET
				fValVec.Append(elMat(i,i));
			}
		}
	
		/* assemble */
		int status;
		AssembleDiagonals(fRowDexVec.Length(), fRowDexVec.Pointer(),
			fValVec.Pointer(), status);

		/* check completion */
		if (!status)
		{
			cout << "\n SPOOLESMatrixT::Assemble: error with equations:\n";
			cout << eqnos << endl;
			throw ExceptionT::kGeneralFail;
		}
	}
	else
	{
		/* assembly mode */
		if (fSymmetric)
		{
			/* sanity check */
			ElementMatrixT::FormatT format = elMat.Format();
			if (format != ElementMatrixT::kSymmetricUpper &&
			    format != ElementMatrixT::kDiagonal && // faster with dedicated function
			    format != ElementMatrixT::kSymmetric)  // faster with dedicated function
			{
				cout << "\n SPOOLESMatrixT::Assemble: unexpected matrix format:"
				     << format << endl;
				throw ExceptionT::kGeneralFail;
			}

			/* equation numbers -> active equation numbers */
			fRowDexVec.Dimension(0);
			fRowEqnVec.Dimension(0);
			for (int j = 0; j < eqnos.Length(); j++)
				if (eqnos[j] > 0)
				{
					fRowDexVec.Append(j);
					fRowEqnVec.Append(eqnos[j] - 1); //OFFSET
				}

			/* sort active rows by global equation number */
			fActiveDex.Alias(fRowDexVec);
			fActiveDex.SortAscending(fRowEqnVec);
			
			/* assemble (global) upper triangle */
			int end_update = fStartEQ + fLocNumEQ - 1;
			int num_active = fActiveDex.Length();
			fValVec.Dimension(num_active); // max size
			int status = 1;
			for (int i = 0; i < num_active && status; i++)
			{
				/* equation number */
				int eq_row = fRowEqnVec[i] + 1; //OFFSET
				if (eq_row >= fStartEQ && eq_row <= end_update)
				{
					/* collect (global) row values */
					int dex_i = fActiveDex[i];
					for (int j = i; j < num_active; j++)
					{
						int dex_j = fActiveDex[j];
						double value = (dex_j >= dex_i) ? // source is upper only
							elMat(dex_i, dex_j) :
							elMat(dex_j, dex_i);
						fValVec[j - i] = value;
					}
				
					/* assemble (global) row */
					AssembleRow(fRowEqnVec[i], num_active - i, fRowEqnVec.Pointer(i),
						fValVec.Pointer(), status);
				}
			}
		}
		else
		{
			/* equation numbers -> active element row numbers */
			fRowDexVec.Dimension(0);
			fColDexVec.Dimension(0);
			int end_update = fStartEQ + fLocNumEQ - 1;
			for (int j = 0; j < eqnos.Length(); j++)
			{
				int eq = eqnos[j];
				if (eq >= fStartEQ && eq <= end_update)
					fRowDexVec.Append(j);
				if (eq > 0)
					fColDexVec.Append(j);
			}

			/* fill element matrix */
			if (elMat.Format() == ElementMatrixT::kSymmetricUpper)
				elMat.CopySymmetric();

			/* copy active block */
			int num_rows = fRowDexVec.Length();
			int num_cols = fColDexVec.Length();
			fActiveBlkMan.SetDimensions(num_rows, num_cols);
			elMat.CopyBlock(fRowDexVec, fColDexVec, fActiveBlk);
		
			/* active equation numbers -> global row numbers */
			for (int r = 0; r < num_rows; r++)
				fRowDexVec[r] = eqnos[fRowDexVec[r]] - 1; //OFFSET

			/* active equation numbers -> global col numbers */
			for (int c = 0; c < num_cols; c++)
				fColDexVec[c] = eqnos[fColDexVec[c]] - 1; //OFFSET
	
			/* row-by-row assembly */
			fValVec.Dimension(num_cols);
			int status = 1;
			for (int i = 0; i < num_rows && status; i++)
			{
				/* copy row values */
				fActiveBlk.CopyRow(i, fValVec);
	
				/* assemble */
				AssembleRow(fRowDexVec[i], num_cols, fColDexVec.Pointer(),
					fValVec.Pointer(), status);
			}
	
			/* check completion */
			if (!status)
			{
				cout << "\n SPOOLESMatrixT::Assemble: error with equations:\n";
				cout << eqnos << endl;
				throw ExceptionT::kGeneralFail;
			}
		}
	}
}

void SPOOLESMatrixT::Assemble(const ElementMatrixT& elMat, const nArrayT<int>& row_eqnos,
	const nArrayT<int>& col_eqnos)
{
#if __option (extended_errorcheck)
	/* check dimensions */
	if (elMat.Rows() != row_eqnos.Length() ||
	    elMat.Cols() != col_eqnos.Length()) throw ExceptionT::kSizeMismatch;
#endif

	/* element matrix format */
	ElementMatrixT::FormatT format = elMat.Format();
	if (format == ElementMatrixT::kDiagonal)
	{
		cout << "\n SPOOLESMatrixT::Assemble(m, r, c): cannot assemble diagonal matrix" << endl;
		throw ExceptionT::kGeneralFail;
	}
	else
	{
		/* assembly mode */
		if (fSymmetric)
		{
			cout << "\n SPOOLESMatrixT::Assemble(m, r, c): cannot assemble symmetric matrix" << endl;
			throw ExceptionT::kGeneralFail;
		}
		else
		{
			/* equation number limit */
			int end_update = fStartEQ + fLocNumEQ - 1;

			/* equation numbers -> active element row numbers */
			fRowDexVec.Dimension(0);
			for (int j = 0; j < row_eqnos.Length(); j++)
			{
				int eq = row_eqnos[j];
				if (eq >= fStartEQ && eq <= end_update)
					fRowDexVec.Append(j);
			}

			/* equation numbers -> active element column numbers */
			fColDexVec.Dimension(0);
			for (int j = 0; j < col_eqnos.Length(); j++)
				if (col_eqnos[j] > 0)
					fColDexVec.Append(j);

			/* fill element matrix */
			if (elMat.Format() == ElementMatrixT::kSymmetricUpper)
				elMat.CopySymmetric();

			/* copy active block */
			int num_rows = fRowDexVec.Length();
			int num_cols = fColDexVec.Length();
			fActiveBlkMan.SetDimensions(num_rows, num_cols);
			elMat.CopyBlock(fRowDexVec, fColDexVec, fActiveBlk);
		
			/* active equation numbers -> global row numbers */
			for (int r = 0; r < num_rows; r++)
				fRowDexVec[r] = row_eqnos[fRowDexVec[r]] - 1; //OFFSET

			/* active equation numbers -> global col numbers */
			for (int c = 0; c < num_cols; c++)
				fColDexVec[c] = col_eqnos[fColDexVec[c]] - 1; //OFFSET
	
			/* row-by-row assembly */
			fValVec.Dimension(num_cols);
			int status = 1;
			for (int i = 0; i < num_rows && status; i++)
			{
				/* copy row values */
				fActiveBlk.CopyRow(i, fValVec);
	
				/* assemble */
				AssembleRow(fRowDexVec[i], num_cols, fColDexVec.Pointer(),
					fValVec.Pointer(), status);
			}
	
			/* check completion */
			if (!status)
			{
				cout << "\n SPOOLESMatrixT::Assemble: error with equations:\n";
				cout << " row:\n" << row_eqnos << endl;
				cout << " col:\n" << col_eqnos << endl;
				throw ExceptionT::kGeneralFail;
			}
		}
	}
}

void SPOOLESMatrixT::Assemble(const nArrayT<double>& diagonal_elMat, const nArrayT<int>& eqnos)
{
#pragma unused(diagonal_elMat)
#pragma unused(eqnos)

	ExceptionT::GeneralFail("SPOOLESMatrixT::Assemble", "not implemented");
}

/* set all matrix values to 0.0 */
void SPOOLESMatrixT::Clear(void) { fval = 0.0; }

/* number scope and reordering */
GlobalMatrixT::EquationNumberScopeT SPOOLESMatrixT::EquationNumberScope(void) const
{
	return kGlobal;
}

bool SPOOLESMatrixT::RenumberEquations(void) const { return false; }

/* assignment operator */
GlobalMatrixT& SPOOLESMatrixT::operator=(const SPOOLESMatrixT& rhs)
{
#pragma unused(rhs)

	cout << "\n SPOOLESMatrixT::operator= : not implemented" << endl;
	throw ExceptionT::kGeneralFail;
	return *this;
}

/* assignment operator */
GlobalMatrixT& SPOOLESMatrixT::operator=(const GlobalMatrixT& rhs)
{
#ifdef __NO_RTTI__
	cout << "\n SPOOLESMatrixT::operator= : requires RTTI" << endl;
	throw ExceptionT::kGeneralFail;
#endif

	const SPOOLESMatrixT* sp = dynamic_cast<const SPOOLESMatrixT*>(&rhs);
	if (!sp) {
		cout << "\n SPOOLESMatrixT::operator= : cast failed" << endl;
		throw ExceptionT::kGeneralFail;
	}
	return operator=(*sp);
}

/** return a clone of self */
GlobalMatrixT* SPOOLESMatrixT::Clone(void) const
{
	SPOOLESMatrixT* new_mat = new SPOOLESMatrixT(*this);
	return new_mat;
}

/*************************************************************************
* Protected
*************************************************************************/

/* precondition matrix */
void SPOOLESMatrixT::Factorize(void)
{
	/* preconditioning done during solve */
	fIsFactorized = 0; // undo GlobalMatrixT flag set
}
	
/* determine new search direction and put the results in result */
void SPOOLESMatrixT::BackSubstitute(dArrayT& result)
{
	/* check */
	if (fTotNumEQ != fLocNumEQ)
	{
		cout << "\n SPOOLESMatrixT::BackSubstitute: expecting total number of equations\n"
		     <<   "     " << fTotNumEQ
		     << " to be equal to the local number of equations " << fLocNumEQ << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* flag should not be set */
	if (fIsFactorized) throw ExceptionT::kGeneralFail;

	/* convert matrix to RCV */
	iArrayT r, c;
	dArrayT v;
	GenerateRCV(r, c, v, 1.0e-15);
	
	/* write matrix */
	if (fCheckCode == kPrintLHS) {
		int old_precision = fOut.precision();
		fOut.precision(12);
		int d_width = fOut.precision() + kDoubleExtra;
		fOut << "\n LHS matrix in {r,c,v} format:\n";
		fOut << " Number of values = " << r.Length() << '\n';
		for (int i = 0; i < r.Length(); i++)
			fOut << setw(kIntWidth) << r[i] + 1
			     << setw(kIntWidth) << c[i] + 1
			     << setw(d_width) << v[i] << '\n';
		fOut.flush();
		fOut.precision(old_precision);
	}

/* solving the system using either the serial or MT driver */
 
 int msglvl = 0; //  0: nothing
 //  1: scalar output (timing data) only
 // >1: verbose
 int matrix_type = SPOOLES_REAL;
 int symmetry_flag = (fSymmetric) ? SPOOLES_SYMMETRIC : SPOOLES_NONSYMMETRIC;
 int pivoting_flag = (fPivoting) ? SPOOLES_PIVOTING : SPOOLES_NO_PIVOTING;
 int seed = 1; /* any seed will do here */
 int OK;
 
#ifdef __SPOOLES_MT__
 /* call MT driver, __NUM_THREADS__ comes from make macro */
 //cout << "Using MT driver with " << __NUM_THREADS__ << " threads" << endl;

 OK = LU_MT_driver(msglvl, SPOOLES_FILE, matrix_type, symmetry_flag,
		   pivoting_flag, seed, result.Length(), result.Pointer(),
		   r.Length(), r.Pointer(), c.Pointer(), v.Pointer(), 
		   __NUM_THREADS__);
#endif
 
#ifndef __SPOOLES_MT__
 //cout << "Using SERIAL driver" << endl;
 OK = LU_serial_driver(msglvl, SPOOLES_FILE, matrix_type, symmetry_flag,
		       pivoting_flag, seed, result.Length(), result.Pointer(),
		       r.Length(), r.Pointer(), c.Pointer(), v.Pointer());
#endif
 
 if (OK != 1)
   {
     cout << "\n SPOOLESMatrixT::BackSubstitute: LU_MT_driver returned: "
	  << OK << endl;
     throw ExceptionT::kGeneralFail;
   }
}

/* rank check functions */
void SPOOLESMatrixT::PrintAllPivots(void) const
{
//not implemented
}

void SPOOLESMatrixT::PrintZeroPivots(void) const
{
//not implemented
}

void SPOOLESMatrixT::PrintLHS(void) const
{
//not implemented
}

/* copy MSR data to RCV */
void SPOOLESMatrixT::GenerateRCV(iArrayT& r, iArrayT& c, dArrayT& v, double drop_tol)
{
	/* count number of values */
	int num_vals = 0;
	double *pval = fval.Pointer();
	for (int i = 0; i < fLocNumEQ; i++) /* diagonals */
		if (fabs(*pval++) > drop_tol) num_vals++;
	pval++;
	for (int i = fLocNumEQ + 1; i < fval.Length(); i++) /* off-diagonals */
		if (fabs(*pval++) > drop_tol) num_vals++;

	/* overall dimension */
	r.Dimension(num_vals);
	c.Dimension(num_vals);
	v.Dimension(num_vals);

	/* start of off-diagonal data (MSR) */
	int* pcol = fbindx.Pointer(fLocNumEQ + 1);
	pval = fval.Pointer(fLocNumEQ + 1);

	/* output rows in ascending column order */
	int shift = fStartEQ - 1; //OFFSET
	int count = 0;
	for (int row = 0; row < fLocNumEQ; row++)
	{
		/* the diagonal */
		double& diag = fval[row];
		if (fabs(diag) > drop_tol) {
			r[count] = row + shift;
			c[count] = row + shift;
			v[count] = diag;
			count++;
		}

		int numvals = fbindx[row+1] - fbindx[row]; /* not incl. diagonal */
		for (int i = 0; i < numvals; i++)
		{
			double& val = *pval++;
			if (fabs(val) > drop_tol) {
				r[count] = row + shift;
				c[count] = *pcol++;
				v[count] = val;
				count++;
			} else pcol++;
		}
	}
	
	/* check */
	if (count != num_vals)
	{
		cout << "\n SPOOLESMatrixT::GenerateRCV: translation error:\n"
		     <<   "   expected number of values = " << num_vals << '\n'
		     <<   "            number of values = " << count << endl;
		throw ExceptionT::kGeneralFail;
	}
}

/*************************************************************************
* Private
*************************************************************************/

/* assemble row values into global structure using the global column
* indices given in coldex - status is 1 if successful, 0 otherwise */
void SPOOLESMatrixT::AssembleRow(int row, int numvals, const int* col_dex,
	const double* rowvals, int& status)
{
	/* find row index */
	int rowkey = AZ_quick_find(row, fupdate.Pointer(), fLocNumEQ, fQF_shift,
		fupdate_bin.Pointer());
	if (rowkey == -1)
		status = 0;
	else	
	{
		/* pointers */
		int*       bindx = fbindx.Pointer();
		double*      val = fval.Pointer();
		int*    srow_dex = fsrow_dex.Pointer();
		double* srow_val = fsrow_val.Pointer();
					
		/* copy */
		memcpy(srow_dex, col_dex, numvals*sizeof(int));
		memcpy(srow_val, rowvals, numvals*sizeof(double));
	
		/* sort copies */
		AZ_sort(srow_dex, numvals, NULL, srow_val);
			
		/* assemble */
		int*    pasmdex = srow_dex;
		double* pasmval = srow_val;

		int  rowlength = bindx[rowkey+1] - bindx[rowkey] + 1; // incl. diagonal
		int*      pdex = bindx + bindx[rowkey];
		double*   pval = val   + bindx[rowkey];
			
		int rowcount = 0;
		int valcount = 0;
		while (rowcount < rowlength && valcount < numvals)
		{
			/* diagonal value */
			if (*pasmdex == row)
			{
				/* assemble */
				val[rowkey] += *pasmval;
			
				/* increment assembled */
				pasmdex++;
				pasmval++;
				valcount++;
			}
			/* matching index */
			else if (*pasmdex == *pdex)
			{
				/* assemble */
				*pval += *pasmval;
					
				/* increment assembled */
				pasmdex++;
				pasmval++;
				valcount++;
			}
			/* no matching index */
			else
			{
				/* increment matrix */
				pdex++;
				pval++;
				rowcount++;
			}
		}
			
		/* check completion */
		if (rowcount == rowlength)
			status = 0;
		else
			status = 1;
	}
}

/* assemble diagonal values into the global structure
* status is 1 if successful, 0 otherwise */
void SPOOLESMatrixT::AssembleDiagonals(int numvals, const int* rows,
	const double* vals, int& status)
{
	/* assemble values */
	status = 1;
	for (int i = 0; i < numvals && status; i++)
	{
		int rowkey = AZ_quick_find(*rows, fupdate.Pointer(), fLocNumEQ,
			fQF_shift, fupdate_bin.Pointer());
		if (rowkey == -1)
			status = 0;
		else
		{
			/* assemble */
			fval[rowkey] += *vals;		
		
			/* next */
			rows++;
			vals++;
		}		
	}
}

/* sets MSR data with column indices sorted in ascending order */
void SPOOLESMatrixT::SetMSRData(void)
{
	/* set update vector - global numbering */
	fupdate.Dimension(fLocNumEQ);
	int* pupdate = fupdate.Pointer();
	int n_update = fStartEQ; //OFFSET
	for (int i = 0; i < fLocNumEQ; i++)
		(*pupdate++) = n_update++;

	/* set MSR structure data */
	fMSRBuilder->SetMSRData(fupdate, fbindx);

	/* shift from equation number to global row number */
	fupdate--; //OFFSET

//NOTE - equation numbers are already global (1,...)
#if 0 	
	/* offset to global equation numbers */
	int offset = fStartEQ - 1;
	if (offset > 0)
	{
		fupdate += offset;
		for (int i = fbindx[0]; i < fbindx.Length(); i++)
			fbindx[i] += offset;
	}
#endif
	
	/* allocate the matrix */
	fval.Dimension(fbindx.Length());

	/* clear equation lists */
	fMSRBuilder->ClearGroups();
}

/* allocate memory for quick find */
void SPOOLESMatrixT::SetUpQuickFind(void)
{
	/* quick find bin */
	fupdate_bin.Dimension(2 + (fLocNumEQ + 4)/4); /* oversize */

	/* initialize shift and bin */
	AZ_init_quick_find(fupdate.Pointer(), fLocNumEQ, &fQF_shift,
		fupdate_bin.Pointer());

	/* find max row length */
	int maxlength = 0;
	for (int i = 0; i < fLocNumEQ; i++)
	{
		int rowln = fbindx[i+1] - fbindx[i];
		maxlength = (rowln > maxlength) ? rowln : maxlength;
	}
	
	/* add space for the diagonal */
	maxlength += 1;

	/* allocate space for sorted row data */
	fsrow_dex.Dimension(maxlength);
	fsrow_val.Dimension(maxlength);
}

int SPOOLESMatrixT::AZ_quick_find(int key, int list[], int length, int shift,
	int bins[])
{
int i, loc, oldkey;
if (length == 0)            return -1;
if (key > list[length - 1]) return -1;
oldkey = key;
key   -= list[0];
if (key < 0) return -1;
loc = key >> shift;
i = AZ_find_index(oldkey, &(list[bins[loc]]), bins[loc + 1] - bins[loc]);
if (i == -1) return -1;
return (i + bins[loc]);

}

int SPOOLESMatrixT::AZ_find_index(int key, int list[], int length)
{
int start, end;
int mid;
if (length == 0) return -1;
start = 0;
end   = length - 1;
while (end - start > 1) {
mid = (start + end) / 2;
if (list[mid] < key) start = mid;
else end = mid;
}
if (list[start] == key) return start;
if (list[end] == key)   return end;
return -1;
}

void SPOOLESMatrixT::AZ_init_quick_find(int list[], int length, int *shift,
	int *bins)
{
register int i, j = 0;
int          range, temp;
if (length == 0) return;
range  = list[length - 1] - list[0];
*shift = 0;
while ((range >> (*shift)) > length / 4)
(*shift)++;
bins[j++] = 0;
for (i = 0; i < length; i++) {
temp = list[i] - list[0];
while ((temp >> (*shift)) >= j)
bins[j++] = i;
}
bins[j] = length;
}

void SPOOLESMatrixT::AZ_sort(int list[], int N, int list2[], double list3[])
{
int    l, r, RR, K, j, i, flag;
int    RR2;
double RR3;
if (N <= 1) return;
l   = N / 2 + 1;
r   = N - 1;
l   = l - 1;
RR  = list[l - 1];
K   = list[l - 1];
if ((list2 != NULL) && (list3 != NULL)) {
RR2 = list2[l - 1];
RR3 = list3[l - 1];
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
list2[i - 1] = list2[j - 1];
list3[i - 1] = list3[j - 1];
}
else {
flag = 0;
}
}
}
list[ i - 1] = RR;
list2[i - 1] = RR2;
list3[i - 1] = RR3;
if (l == 1) {
RR  = list [r];
RR2 = list2[r];
RR3 = list3[r];
K = list[r];
list[r ] = list[0];
list2[r] = list2[0];
list3[r] = list3[0];
r = r - 1;
}
else {
l   = l - 1;
RR  = list[ l - 1];
RR2 = list2[l - 1];
RR3 = list3[l - 1];
K   = list[l - 1];
}
}
list[ 0] = RR;
list2[0] = RR2;
list3[0] = RR3;
}
else if (list2 != NULL) {
RR2 = list2[l - 1];
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
list2[i - 1] = list2[j - 1];
}
else {
flag = 0;
}
}
}
list[ i - 1] = RR;
list2[i - 1] = RR2;
if (l == 1) {
RR  = list [r];
RR2 = list2[r];
K = list[r];
list[r ] = list[0];
list2[r] = list2[0];
r = r - 1;
}
else {
l   = l - 1;
RR  = list[ l - 1];
RR2 = list2[l - 1];
K   = list[l - 1];
}
}
list[ 0] = RR;
list2[0] = RR2;
}
else if (list3 != NULL) {
RR3 = list3[l - 1];
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
list3[i - 1] = list3[j - 1];
}
else {
flag = 0;
}
}
}
list[ i - 1] = RR;
list3[i - 1] = RR3;
if (l == 1) {
RR  = list [r];
RR3 = list3[r];

K = list[r];
list[r ] = list[0];
list3[r] = list3[0];
r = r - 1;
}
else {
l   = l - 1;
RR  = list[ l - 1];
RR3 = list3[l - 1];
K   = list[l - 1];
}
}
list[ 0] = RR;
list3[0] = RR3;
}
else {
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
}

#endif /* __SPOOLES__ */
