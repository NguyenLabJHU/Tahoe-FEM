/* $Id: SLUMatrix.cpp,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: rbridson (06/30/2000)                                         */
/* Implementation of interface to SuperLU solver library.                 */

#include "SLUMatrix.h"

/* library support */
#ifdef __SUPERLU__

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include "Constants.h"
#include "ExceptionCodes.h"

/* types that we use in these methods */
#include "dMatrixT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "RaggedArray2DT.h"
#include "ElementMatrixT.h"

/***************************************************************************
* Public
***************************************************************************/

/* Constructor */
SLUMatrix::SLUMatrix(ostream& out, int check_code):
	GlobalMatrixT(out, check_code)
{
	/* set up NULL structures */
	fMatrix.Stype = NC;      /* column-wise, no supernodes */
	fMatrix.Dtype = _DOUBLE; /* storing doubles */
	fMatrix.Mtype = GE;      /* general matrix */
	fMatrix.nrow = fLocNumEQ;
	fMatrix.ncol = fLocNumEQ;
	fMatrix.Store = malloc(sizeof(NCformat));
	if (!fMatrix.Store) throw eOutOfMemory;

	NCformat *A = (NCformat*)fMatrix.Store;
	A->nnz = 0;
	A->nzval = NULL;
	A->rowind = NULL;

	/* The only important thing to initialize in L and U are pointers */
	fLower.Store = NULL;
	fUpper.Store = NULL;
	fLUallocated = false;
}

/* Destructor */	
SLUMatrix::~SLUMatrix(void)
{
	Destroy_CompCol_Matrix (&fMatrix);

	if (fLUallocated)
	{
		Destroy_SuperNode_Matrix (&fLower);
		Destroy_CompCol_Matrix (&fUpper);
	}

	free(fPerm_c);
	free(fPerm_r);
	free(fEtree);

	/* SuperLU statistics */
	StatFree();
}

/* set the internal matrix structure.
* NOTE: do not call Initialize() until equation topology has been set
* with AddEquationSet() for all equation sets */
void SLUMatrix::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* inherited */
	GlobalMatrixT::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* check */
	if (tot_num_eq != loc_num_eq)
	{
		cout << "\n SLUMatrix::Initialize: expecting total number of equations\n"
		     <<   "     " << tot_num_eq
		     << " to be equal to the local number of equations " << loc_num_eq << endl;
		throw eGeneralFail;
	}

	/* A note on memory allocation: since SuperLU is a C library, */
	/* I use malloc/free instead of new/delete for the structures */
	/* that SuperLU accesses, just in case. */	
	fMatrix.nrow = fLocNumEQ;
	fMatrix.ncol = fLocNumEQ;

	NCformat *A = (NCformat*)fMatrix.Store;
	A->colptr = (int*) calloc(fLocNumEQ+1, sizeof(int));

	fPerm_c = (int*) malloc (fLocNumEQ*sizeof(int));
	if (!fPerm_c) throw eOutOfMemory;
	fIsColOrdered = false;
	fPerm_r = (int*) malloc (fLocNumEQ*sizeof(int));
	if (!fPerm_r) throw eOutOfMemory;
	fEtree = (int*) malloc (fLocNumEQ*sizeof(int));
	if (!fEtree) throw eOutOfMemory;

	/* SuperLU statistics */
	StatInit (kPanelSize, kSupernodeRelax);

	/* structure could be changing, so get rid of old factors etc. */
	if (fLUallocated)
	{
		Destroy_SuperNode_Matrix (&fLower);
		Destroy_CompCol_Matrix (&fUpper);
		fLUallocated = false;
	}
	fIsColOrdered = false;

	/* We now construct the sparsity pattern of A from the equation sets */
	/* check if A is already allocated */
	if (A->rowind)
	{
		free (A->rowind);
		free (A->nzval);
	}

	/* Begin by (over-)estimating the number of nonzeros per column */
	int *colLength = (int*) malloc (fLocNumEQ*sizeof(int));
	EstimateNNZ (colLength, A->nnz);

	/* Now allocate enough room for row indices (wait until later for */
	/* the nonzero values themselves) */
	A->rowind = (int*) malloc (A->nnz*sizeof(int));
	if (!A->rowind) throw eOutOfMemory;

	/* Using the upper bounds in colLength, set up provisional column */
	/* pointers */
	A->colptr[0] = 0;
	for (int i = 0; i < fLocNumEQ; ++i)
		A->colptr[i+1] = A->colptr[i] + colLength[i];

	/* and now go through all the elements, inserting all the equations */
	InsertEquations (A, colLength, A->nnz);

	/* Then we can compress A to eliminate spaces between columns */
	CompressColumns (A, colLength);
	free (colLength);  // no longer needed

	/* and finish by reallocating A->rowind and allocating A->nzval */
	A->rowind = (int*) realloc (A->rowind, A->nnz*sizeof(int));
	if (!A->rowind) throw eOutOfMemory;
	A->nzval = (void*) malloc (A->nnz*sizeof(double));
	if (!A->nzval) throw eOutOfMemory;

	/* output */
	fOut <<" Number of nonzeros in global matrix = "<< A->nnz <<"\n"<<endl;
	
	/* clear stored equation sets */
	fEqnos.Clear();
	fRaggedEqnos.Clear();	
}

/* set all matrix values to 0.0 */
void SLUMatrix::Clear(void)
{
	/* inherited */
	GlobalMatrixT::Clear();

	/* clear entries */
	NCformat *ncmat = (NCformat*)fMatrix.Store;
	memset (ncmat->nzval, 0, sizeof(double)*ncmat->colptr[fLocNumEQ]);

	/* if old factors are hanging around, get rid of them now */
	if (fLUallocated)
	{
		Destroy_SuperNode_Matrix (&fLower);
		Destroy_CompCol_Matrix (&fUpper);
		fLUallocated = false;
	}
}

/* add element group equations to the overall topology.
* NOTE: assembly positions (equation numbers) = 1...fLocNumEQ
* equations can be of fixed size (iArray2DT) or
* variable length (RaggedArray2DT) */
void SLUMatrix::AddEquationSet(const iArray2DT& eqset)
{
	fEqnos.AppendUnique (&eqset);
}

/* see AddEquationSet above */
void SLUMatrix::AddEquationSet(const RaggedArray2DT<int>& eqset)
{
	fRaggedEqnos.AppendUnique (&eqset);
}

/* assemble the element contribution into the LHS matrix - assumes
* that elMat is square (n x n) and that eqnos is also length n.
* NOTE: assembly positions (equation numbers) = 1...fLocNumEQ */
void SLUMatrix::Assemble(const ElementMatrixT& elMat, const iArrayT& eqnos)
{
	/* element matrix format */
	ElementMatrixT::FormatT format = elMat.Format();

	/* two cases: element matrix is diagonal, or it's not. */
	if (format == ElementMatrixT::kDiagonal)
	{
		/* less work to do! We only add diagonal entries */
		double *pelMat = elMat.Pointer();
		int inc = elMat.Rows() + 1; // how far apart diag entries are

		int nee = eqnos.Length();
		for (int eqdex = 0; eqdex < nee; ++eqdex)
		{
			int eqno = eqnos[eqdex] - 1;
			if (eqno > -1)   // active dof?
				(*this)(eqno,eqno) += *pelMat;
			pelMat += inc;
		}
	}
	else    /* otherwise there's a full matrix to deal with */
	{
		/* If it's symmetric and just a triangle is stored, */
		/* copy it over to get full storage */
		if (format == ElementMatrixT::kSymmetricUpper)
			elMat.CopySymmetric();

		int nee = eqnos.Length();  // number of equations for element
		for (int col = 0; col < nee; ++col)
		{
			int ceqno = eqnos[col] - 1;
			if (ceqno > -1)   // active dof?
			{
				for (int row = 0; row < nee; ++row)
				{
					int reqno = eqnos[row] - 1;
					if (reqno > -1) // active dof?
						(*this)(reqno,ceqno) += elMat(row,col);
				}
			}
		}
	}
}

/* assignment operator */
GlobalMatrixT& SLUMatrix::operator=(const GlobalMatrixT& RHS)
{
	// not implemented yet
	throw eGeneralFail;
	return GlobalMatrixT::operator=(RHS);
}

/* element accessor - READ ONLY */
double SLUMatrix::Element(int row, int col) const
{
	try
	{
		return (*this)(row,col);
	}
	catch (int code)
	{
		if (code == eOutOfRange)
			return 0.0;
		else
			throw code;
	}
	return 0.0;
}

/* number scope and reordering */
GlobalMatrixT::EquationNumberScopeT SLUMatrix::EquationNumberScope(void) const
{
	return kLocal;
}

bool SLUMatrix::RenumberEquations(void) const { return false; }

/***********************************************************************
* Protected
***********************************************************************/

/* output in sparse format */
ostream& operator<<(ostream& out, const SLUMatrix& matrix)
{
	NCformat *A = (NCformat*)(matrix.fMatrix.Store);
	int i, j;

	for (i = 0; i < matrix.fMatrix.ncol; ++i)
	{
//		out << "column " << i+1 << ":  [row] value\n";
		out << -1 << " " << 0.0 << "\n";
		for (j = A->colptr[i]; j < A->colptr[i+1]; ++j)
		{
//			out << "  [" << A->rowind[j]+1 << "]  ";
			out << A->rowind[j]+1 << " ";
			out << ((double*)A->nzval)[j] << "\n";
		}
	}

	return out;
}

/* decompose matrix into PLU */
void SLUMatrix::Factorize(void)
{
	int info;
	SuperMatrix AC;

	if (fLUallocated)
	{
		Destroy_SuperNode_Matrix (&fLower);
		Destroy_CompCol_Matrix (&fUpper);
		fLUallocated = false;
	}

	if (!fIsColOrdered)  /* if we don't have a column ordering yet */
	{
		/* figure out a good column ordering of A */
		OrderColumns();
		fIsColOrdered = true;
		/* compute etree and permute cols of A to get AC */
		sp_preorder ("N", &fMatrix, fPerm_c, fEtree, &AC);
	}
	else
	{	/* otherwise just use existing fPerm_c to permute cols of A */
		sp_preorder ("Y", &fMatrix, fPerm_c, fEtree, &AC);
	}

	/* do the factorization, figuring out the row ordering */
	dgstrf ("N", &AC, 1, 0.0, kSupernodeRelax, kPanelSize, fEtree,
	        NULL, 0, fPerm_r, fPerm_c, &fLower, &fUpper, &info);
	fLUallocated = true;

{
fOut << "Factorization: n = " << fLocNumEQ;
fOut << ", nnz(A) = " << ((NCformat*)fMatrix.Store)->nnz;
fOut << ", nnz(L) = " << ((SCformat*)fLower.Store)->nnz;
fOut << ", nnz(U) = " << ((NCformat*)fUpper.Store)->nnz;
fOut << endl;
}

	/* forget about the column-permuted version of A - but remember that */
	/* it is pointing to the same arrays as A, so don't free those */
	free(((NCPformat*)AC.Store)->colbeg);
	free(((NCPformat*)AC.Store)->colend);
	free(AC.Store);

	if (info) throw eGeneralFail;
}

/* solution driver */
void SLUMatrix::BackSubstitute(dArrayT& result)
{
	int info;
	SuperMatrix B;

	/* Put result (initially right-hand side) into supermatrix B */
	B.Stype = DN;
	B.Dtype = _DOUBLE;
	B.Mtype = GE;
	B.nrow = fLocNumEQ;
	B.ncol = 1;
	B.Store = malloc (sizeof(DNformat));
	if (!B.Store) throw eGeneralFail;
	((DNformat*)B.Store)->lda = fLocNumEQ;
((DNformat*)B.Store)->nzval = result.Pointer();

	/* solve */
	dgstrs ("N", &fLower, &fUpper, fPerm_r, fPerm_c, &B, &info);

	/* check for errors */
	if (info) throw eGeneralFail;

	/* B (hence result) now contains the solution. */

	/* free up the record we allocated */
	free (B.Store);
}

/* check functions */
void SLUMatrix::PrintAllPivots(void) const
{
// temp - not implemented yet. Maybe inappropriate
}

void SLUMatrix::PrintZeroPivots(void) const
{
// temp - not implemented yet. Maybe inappropriate
}

void SLUMatrix::PrintLHS(void) const
{
	if (fCheckCode != GlobalMatrixT::kPrintLHS)
		return;

	fOut << "\nLHS matrix:\n\n";
	fOut << (*this) << "\n\n";
}

/* element accessor - read and write, for assembly. Exception for */
/* access to unstored zero. row and col are */
double& SLUMatrix::operator()(int row, int col) const
{
	/* range checks */
	if (row < 0 || row >= fLocNumEQ) throw eGeneralFail;
	if (col < 0 || col >= fLocNumEQ) throw eGeneralFail;

	NCformat *A = (NCformat*)fMatrix.Store;

	/* look through column col for the given row index */
	for (int i = A->colptr[col]; i < A->colptr[col+1]; ++i)
	{
		if (row == A->rowind[i])
			return ((double*)A->nzval)[i];
	}
	/* otherwise this nonzero wasn't present */
	throw eOutOfRange;
}

/* (over)estimate the number of nonzeros, based on the equation sets */
void SLUMatrix::EstimateNNZ (int *colLength, int &totalnnz)
{
	/* The estimate is simple: forget about overlap of nodes between */
	/* elements, and just add up all the nonzeros of all the element */
	/* stiffness matrices. */
	
	totalnnz = 0;
	memset (colLength, 0, fLocNumEQ*sizeof(int));

	const iArray2DT* peq;
	fEqnos.Top();
	while (fEqnos.Next(peq))
	{
		int nel = peq->MajorDim();  /* number of elements */
		int nee = peq->MinorDim();  /* number of element equations */
		for (int j = 0; j < nel; ++j)
		{
			int *eleqnos = (*peq)(j);
			for (int k = 0; k < nee; ++k)
			{
				int eq = eleqnos[k]-1;
				if (eq > -1)
				{
					colLength[eq] += nee;
					totalnnz += nee;
				}
			}
		}
	}

	const RaggedArray2DT<int>* prageq;
	fRaggedEqnos.Top();
	while (fRaggedEqnos.Next(prageq))
	{
		int nel = prageq->MajorDim();  /* number of elements */
		for (int j = 0; j < nel; ++j)
			{
			int nee = prageq->MinorDim(j); /* no. element equations */
			int *eleqnos = (*prageq)(j);
			for (int k = 0; k < nee; ++k)
			{
				int eq = eleqnos[k]-1;
				if (eq > -1)
				{
					colLength[eq] += nee;
					totalnnz += nee;
				}
			}
		}
	}
}

/* insert all the element equations into A */
void SLUMatrix::InsertEquations (NCformat *A, int *colLength, int &nnz)
{
	/* Reset the lengths */
	memset (colLength, 0, fLocNumEQ*sizeof(int));
	nnz = 0;

	const iArray2DT* peq;
	fEqnos.Top();
	while (fEqnos.Next(peq))
	{
		int nel = peq->MajorDim();  /* number of elements */
		int nee = peq->MinorDim();  /* number of element equations */
		for (int j = 0; j < nel; ++j)
		{
			int *eleqnos = (*peq)(j);
			for (int k = 0; k < nee; ++k)
			{
				int eq = eleqnos[k]-1;
				if (eq > -1)
					InsertNZ (A, colLength, nnz, eq, nee, eleqnos);
			}
		}
	}

	const RaggedArray2DT<int>* prageq;
	fRaggedEqnos.Top();
	while (fRaggedEqnos.Next(prageq))
	{
		int nel = prageq->MajorDim();  /* number of elements */
		for (int j = 0; j < nel; ++j)
			{
			int nee = prageq->MinorDim(j); /* no. element equations */
			int *eleqnos = (*prageq)(j);
			for (int k = 0; k < nee; ++k)
			{
				int eq = eleqnos[k]-1;
				if (eq > -1)
					InsertNZ (A, colLength, nnz, eq, nee, eleqnos);
			}
		}
	}
}

/* Insert the list nzlist of nonzeros (with length nzlen) into column c of
* matrix A, using colLengths to keep track of column lengths, and keeping
* nnz up-to-date. The columns of A will have nonzeros in ascending order.
* The entries of nzlist are FORTRAN-indexed (starting at 1, so entries 0 or
* less are ignored) whereas A and c are C-indexed (starting at 0).
*/
void SLUMatrix::InsertNZ (NCformat *A, int *colLength, int &nnz, int c,
int nzlen, int *nzlist)
{
	for (int i = 0; i < nzlen; ++i)
	{
		int newnz = nzlist[i]-1;

		if (newnz < 0)
			continue;   /* ignore inactive dof's */

		/* insert newnz into column c */
		int j = 0;
		while (j < colLength[c])
		{
			int nzj = A->rowind[A->colptr[c]+j];
			if (newnz < nzj)
				break;
			else if (newnz == nzj)
				goto AlreadyThere; // column c already has it
			++j;
		}
		/* insert newnz before j */
		++colLength[c];  // new nonzero
		++nnz;
		/* move up list by one */
		for (int k = colLength[c]-1; k > j; --k)
			A->rowind[A->colptr[c]+k] = A->rowind[A->colptr[c]+k-1];
		/* put newnz in */
		A->rowind[A->colptr[c]+j] = newnz;

		AlreadyThere: {} // nothing to do if nzlist[i]-1 is there
	}
}

/* compress columns in A */
void SLUMatrix::CompressColumns (NCformat *A, const int *colLength)
{
	int writeto = 0,  /* where we are writing to in A->rowind */
	    colstart;     /* where a column will start in A->rowind */

	/* look over the columns */
	for (int i = 0; i < fLocNumEQ; ++i)
	{
		colstart = writeto;    /* where column i will start */
		/* copy i's list of nonzeros down to writeto */
		for (int j = A->colptr[i]; j < A->colptr[i]+colLength[i]; ++j)
		{
			A->rowind[writeto] = A->rowind[j];
			++writeto;
		}
		/* set column i's pointer to the new location in A->rowind */
		A->colptr[i] = colstart;
	}
	/* set the last pointer to writeto, which now should == nnz */
	A->colptr[fLocNumEQ] = writeto;
}

/* Figure out a good ordering of the columns in fPerm_c
* Because of its construction, the stiffness matrix has symmetric structure
* even if it doesn't have numerical symmetry. This means we have several
* options: any column ordering of course (e.g. colamd), standard orderings (if
* we don't think much pivoting will take place), or nested dissection with
* wide separators (and we don't need to explicitly form the structure of
* A'*A or A'+A).
*/
void SLUMatrix::OrderColumns (void)
{
NCformat *A = (NCformat*)fMatrix.Store;
get_colamd (fLocNumEQ, fLocNumEQ, A->nnz, A->colptr, A->rowind, fPerm_c);
}

#endif /* __SUPERLU__ */
