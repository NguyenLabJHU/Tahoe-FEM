/* $Id: SLUMatrix.cpp,v 1.10 2004-07-15 08:31:41 paklein Exp $ */
/* created: rbridson (06/30/2000) */
#include "SLUMatrix.h"

/* library support */
#ifdef __SUPERLU__

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include "toolboxConstants.h"
#include "ExceptionT.h"

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

using namespace Tahoe;

/* constructor */
SLUMatrix::SLUMatrix(ostream& out, int check_code):
	GlobalMatrixT(out, check_code),
	fperm_c(NULL),
	fperm_r(NULL),
	fetree(NULL)
{
	/* set up NULL structures */
	fA.Stype = SLU_NC; /* column-wise, no supernodes */
	fA.Dtype = SLU_D;  /* storing doubles */
	fA.Mtype = SLU_GE; /* general matrix */
	fA.nrow = 0;
	fA.ncol = 0;
	fA.Store = malloc(sizeof(NCformat));
	if (!fA.Store) ExceptionT::OutOfMemory("SLUMatrix::SLUMatrix");

	NCformat *A = (NCformat*) fA.Store;
	A->nnz = 0;
	A->nzval = NULL;
	A->rowind = NULL;
	A->colptr = NULL;

	/* The only important thing to initialize in L and U are pointers */
	fL.Store = NULL;
	fU.Store = NULL;
	fLUallocated = false;

    /* Set the default input options:
		options.Fact = DOFACT;
		options.Equil = YES;
    	options.ColPerm = COLAMD;
		options.DiagPivotThresh = 1.0;
    	options.Trans = NOTRANS;
    	options.IterRefine = NOREFINE;
    	options.SymmetricMode = NO;
    	options.PivotGrowth = NO;
    	options.ConditionNumber = NO;
    	options.PrintStat = YES; */
    set_default_options(&foptions);
}

/* Destructor */	
SLUMatrix::~SLUMatrix(void)
{
	/* free the matrix */
	Destroy_CompCol_Matrix(&fA);

	/* free upper and lower factors */
	if (fLUallocated) {
		Destroy_SuperNode_Matrix(&fL);
		Destroy_CompCol_Matrix(&fU);
	}

	free(fperm_c);
	free(fperm_r);
	free(fetree);
}

/* set the internal matrix structure.
* NOTE: do not call Initialize() until equation topology has been set
* with AddEquationSet() for all equation sets */
void SLUMatrix::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	const char caller[] = "SLUMatrix::Initialize";

	/* inherited */
	GlobalMatrixT::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* check */
	if (tot_num_eq != loc_num_eq)
		ExceptionT::GeneralFail(caller, "tot eq %d != loc eq %d", tot_num_eq, loc_num_eq);

	/* reset options */
	foptions.Fact = DOFACT;

	/* A note on memory allocation: since SuperLU is a C library, */
	/* I use malloc/free instead of new/delete for the structures */
	/* that SuperLU accesses, just in case. */	
	fA.nrow = fLocNumEQ;
	fA.ncol = fLocNumEQ;

	NCformat *A = (NCformat*) fA.Store;

	free(A->colptr);
	A->colptr = (int*) calloc(fLocNumEQ+1, sizeof(int));
	if (!A->colptr) ExceptionT::OutOfMemory(caller);

	free(fperm_c);
	fperm_c = (int*) malloc (fLocNumEQ*sizeof(int));
	if (!fperm_c) ExceptionT::OutOfMemory(caller);

	fIsColOrdered = false;

	free(fperm_r);
	fperm_r = (int*) malloc (fLocNumEQ*sizeof(int));
	if (!fperm_r) ExceptionT::OutOfMemory(caller);

	free(fetree);
	fetree = (int*) malloc (fLocNumEQ*sizeof(int));
	if (!fetree) ExceptionT::OutOfMemory(caller);

	/* structure could be changing, so get rid of old factors etc. */
	if (fLUallocated) {
		Destroy_SuperNode_Matrix(&fL);
		Destroy_CompCol_Matrix(&fU);
		fLUallocated = false;
	}
	fIsColOrdered = false;

	/* We now construct the sparsity pattern of A from the equation sets */
	/* check if A is already allocated */
	if (A->rowind) {
		free(A->rowind);
		free(A->nzval);
	}

	/* Begin by (over-)estimating the number of nonzeros per column */
	int *colLength = (int*) malloc (fLocNumEQ*sizeof(int));
	if (!colLength) ExceptionT::OutOfMemory(caller);
	EstimateNNZ(colLength, A->nnz);

	/* Now allocate enough room for row indices (wait until later for */
	/* the nonzero values themselves) */
	A->rowind = (int*) malloc (A->nnz*sizeof(int));
	if (!A->rowind) ExceptionT::OutOfMemory(caller);

	/* Using the upper bounds in colLength, set up provisional column */
	/* pointers */
	A->colptr[0] = 0;
	for (int i = 0; i < fLocNumEQ; ++i)
		A->colptr[i+1] = A->colptr[i] + colLength[i];

	/* and now go through all the elements, inserting all the equations */
	InsertEquations (A, colLength, A->nnz);

	/* Then we can compress A to eliminate spaces between columns */
	CompressColumns (A, colLength);
	free(colLength);  // no longer needed

	/* and finish by reallocating A->rowind and allocating A->nzval */
	A->rowind = (int*) realloc (A->rowind, A->nnz*sizeof(int));
	if (!A->rowind) ExceptionT::OutOfMemory(caller);

	A->nzval = (void*) malloc (A->nnz*sizeof(double));
	if (!A->nzval) ExceptionT::OutOfMemory(caller);

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
	NCformat *A = (NCformat*) fA.Store;
	memset(A->nzval, 0, sizeof(double)*A->colptr[fLocNumEQ]);

	/* if old factors are hanging around, get rid of them now */
	if (fLUallocated) {
		Destroy_SuperNode_Matrix (&fL);
		Destroy_CompCol_Matrix (&fU);
		fLUallocated = false;
	}

	/* set options flag for next call to solve */
	foptions.Fact = SamePattern;
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
void SLUMatrix::Assemble(const ElementMatrixT& elMat, const nArrayT<int>& eqnos)
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

void SLUMatrix::Assemble(const ElementMatrixT& elMat, const nArrayT<int>& row_eqnos,
	const nArrayT<int>& col_eqnos)
{
#pragma unused(elMat)
#pragma unused(row_eqnos)
#pragma unused(col_eqnos)
	ExceptionT::GeneralFail("SLUMatrix::Assemble", "non-square not implemented");
}

void SLUMatrix::Assemble(const nArrayT<double>& diagonal_elMat, const nArrayT<int>& eqnos)
{
#pragma unused(diagonal_elMat)
#pragma unused(eqnos)
	ExceptionT::GeneralFail("SLUMatrix::Assemble", "diagonal not implemented");
}

/* assignment operator */
const GlobalMatrixT& SLUMatrix::operator=(const GlobalMatrixT&)
{
	/* not implemented */
	ExceptionT::GeneralFail("SLUMatrix::operator=", "not implemented");
	return *this;
}

/* number scope and reordering */
GlobalMatrixT::EquationNumberScopeT SLUMatrix::EquationNumberScope(void) const
{
	return kLocal;
}

bool SLUMatrix::RenumberEquations(void) const { return false; }

/* return a clone of self */
GlobalMatrixT* SLUMatrix::Clone(void) const
{
	/* not implemented */
	ExceptionT::GeneralFail("SLUMatrix::operator=", "not implemented");
	return (GlobalMatrixT*) this;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* decompose matrix into PLU */
void SLUMatrix::Factorize(void) { /* nothing to do here */ }

/* solution driver */
void SLUMatrix::BackSubstitute(dArrayT& result)
{
	const char caller[] = "SLUMatrix::BackSubstitute";

	/* put result (initially right-hand side) into supermatrix B */
	SuperMatrix B;
	B.Stype = SLU_DN;
	B.Dtype = SLU_D;
	B.Mtype = SLU_GE;
	B.nrow = fLocNumEQ;
	B.ncol = 1;
	B.Store = malloc(sizeof(DNformat));
	if (!B.Store) ExceptionT::GeneralFail(caller);

	DNformat* BStore = (DNformat*) B.Store;
	BStore->lda = fLocNumEQ;
	BStore->nzval = result.Pointer();

	/* solve */
	int info;


    /* Initialize the statistics variables. */
	SuperLUStat_t stat;    
    StatInit(&stat);
    
    /* ------------------------------------------------------------
       WE SOLVE THE LINEAR SYSTEM FOR THE FIRST TIME: AX = B
       ------------------------------------------------------------*/
    dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
           &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
           &mem_usage, &stat, &info);

    printf("First system: dgssvx() returns info %d\n", info);

    if ( info == 0 || info == n+1 ) {

        /* This is how you could access the solution matrix. */
        double *sol = (double*) ((DNformat*) X.Store)->nzval; 

	if ( options.PivotGrowth ) printf("Recip. pivot growth = %e\n", rpg);
	if ( options.ConditionNumber )
	    printf("Recip. condition number = %e\n", rcond);
        Lstore = (SCformat *) L.Store;
        Ustore = (NCformat *) U.Store;
	printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
    	printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
    	printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
	       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
	       mem_usage.expansions);
	if ( options.IterRefine ) {
            printf("Iterative Refinement:\n");
	    printf("%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
	    for (i = 0; i < nrhs; ++i)
	      printf("%8d%8d%16e%16e\n", i+1, stat.RefineSteps, ferr[i], berr[i]);
	}
	fflush(stdout);

    } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: %d bytes\n", info - n);
    }

    if ( options.PrintStat ) StatPrint(&stat);
    StatFree(&stat);

	/* check for errors */
	//if (info) throw ExceptionT::kGeneralFail;

	/* B (hence result) now contains the solution. */

	/* free up the record we allocated */
	free(B.Store);
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

#endif /* __SUPERLU__ */
