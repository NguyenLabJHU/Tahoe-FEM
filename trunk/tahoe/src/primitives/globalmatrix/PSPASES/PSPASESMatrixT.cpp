/* $Id: PSPASESMatrixT.cpp,v 1.1 2004-03-14 00:09:40 paklein Exp $ */
/* created: paklein (09/13/2000) */
#include "PSPASESMatrixT.h"

/* library support options */
#ifdef __PSPASES__

#include "CommunicatorT.h"
#include "ElementMatrixT.h"
#include "MSRBuilderT.h"

using namespace Tahoe;

/* constuctor */
PSPASESMatrixT::PSPASESMatrixT(ostream& out, int check_code, CommunicatorT& comm):
	GlobalMatrixT(out, check_code),
	fComm(comm),
	fBuilder(NULL),
	faptrs_man(10, faptrs, 2),
	fIsSymFactorized(0)
{
	const char caller[] = "PSPASESMatrixT::PSPASESMatrixT";

	/* verify that number of processes is power of 2 */
	int size = fComm.Size();
	int power2 = 2;
	while (size != power2) {
		if (power2 > size)
			ExceptionT::GeneralFail(caller, "PSPACES requires nproc as 2^n with n >= 1");
		power2 *= 2;
	}

	fBuilder = new MSRBuilderT(false);
	if (!fBuilder) ExceptionT::OutOfMemory(caller);
}

PSPASESMatrixT::PSPASESMatrixT(const PSPASESMatrixT& source):
	GlobalMatrixT(source),
	fComm(source.fComm)
{
	ExceptionT::GeneralFail("PSPASESMatrixT::PSPASESMatrixT", "not implemented");
}

/* destructor */
PSPASESMatrixT::~PSPASESMatrixT(void)
{
	delete fBuilder;

	/* free storage (order matters) */
	int option_0 = 0;
	if (fIsFactorized) PSPACEC(&fNcomm, &option_0);
	if (fIsSymFactorized) PSPACEC(&fYcomm, &option_0);
}

/* add to structure */
void PSPASESMatrixT::AddEquationSet(const iArray2DT& eqnos) { fBuilder->AddGroup(eqnos); }
void PSPASESMatrixT::AddEquationSet(const RaggedArray2DT<int>& eqnos) { fBuilder->AddGroup(eqnos); }

/* set the internal matrix structure */
void PSPASESMatrixT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* free space */
	int option_0 = 0;
	if (fIsFactorized) PSPACEC(&fNcomm, &option_0);
	if (fIsSymFactorized) {
		PSPACEC(&fYcomm, &option_0);
		fIsSymFactorized;
	}

	/* inherited - initialize MSR data */
	GlobalMatrixT::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* exchange number of equations per processor */
	iArrayT num_eq(fComm.Size());
	fComm.AllGather(loc_num_eq, num_eq);
	
	/* set offsets array */
	frowdist.Dimension(num_eq.Length() + 1);
	frowdist[0] = 0;
	for (int i = 0; i < num_eq.Length(); i++)
		frowdist[i+1] = frowdist[i] + num_eq[i];

	/* set update vector - global numbering */
	iArrayT activerows(fLocNumEQ);
	int n_update = fStartEQ; //OFFSET
	for (int i = 0; i < fLocNumEQ; i++)
		activerows[i] = n_update++;

	/* set structure data */
	iArray2DT aptrs_tmp;
	iArrayT ainds_tmp;
	fBuilder->SetPSPASESData(activerows, aptrs_tmp, ainds_tmp);
	faptrs_man.SetMajorDimension(aptrs_tmp.MajorDim(), false);
	faptrs = aptrs_tmp;
	fainds = ainds_tmp;
	
	/* dimension other work space */
	favals.Dimension(fainds.Length());
	forder.Dimension(loc_num_eq);
	fsizes.Dimension(2*fComm.Size());
	fioptions_PSPACEO.Dimension(16);
	fioptions_PSPACEO = 0;
	fioptions_PSPACEY.Dimension(16);
	fioptions_PSPACEY = 0;
	fioptions_PSPACEY[0] = 64;
	fdoptions_PSPACEY.Dimension(16);
	fdoptions_PSPACEY = 0.0;
	fioptions_DPSPACET.Dimension(8);
	fioptions_DPSPACET = 0;
	fioptions_DPSPACET[0] = 1;

	/* clear equation lists */
	fBuilder->ClearGroups();
}

/* set all matrix values to 0.0 */
void PSPASESMatrixT::Clear(void) {
	dArrayT tmp;
	tmp.Alias(favals);
	tmp = 0.0;
	
	int option_0 = 0;
	if (fIsFactorized) {
		PSPACEC(&fNcomm, &option_0);	
		fIsFactorized = 0;
	}
}

/* assemble the element contribution */
void PSPASESMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos)
{
	const char caller[] = "PSPASESMatrixT::Assemble";

	/* element matrix format */
	ElementMatrixT::FormatT format = elMat.Format();

	/* two cases: element matrix is diagonal, or it's not. */
	if (format == ElementMatrixT::kDiagonal)
	{
		/* diagonal entries only */
		const double *pelMat = elMat.Pointer();
		int inc = elMat.Rows() + 1; /* offset between diag entries are */

		int nee = eqnos.Length();
		for (int eqdex = 0; eqdex < nee; ++eqdex) {
			int eqno = eqnos[eqdex] - 1;
			if (eqno > -1) /* active eqn */ { 
				double* a = (*this)(eqno, eqno);
				if (a)
					*a += *pelMat;
				else
					ExceptionT::OutOfRange(caller);
			}
			pelMat += inc;
		}
	}
	else if (format == ElementMatrixT::kNonSymmetric)
	{
		int nee = eqnos.Length();  // number of equations for element
		for (int col = 0; col < nee; ++col)
		{
			int ceqno = eqnos[col] - 1;
			if (ceqno > -1) /* active eqn */ {
				for (int row = 0; row <= col; ++row) {
					int reqno = eqnos[row] - 1;
					if (reqno > -1) /* active eqn */ {
					
						/* diagonal term */
						if (reqno == ceqno) {
							double* a = (*this)(reqno,ceqno);
							if (a)
								*a += elMat(row,col);
							else
								ExceptionT::OutOfRange(caller);
						}
						else /* off-diagonal term */ {
						
							/* symmetrize */
							double v = 0.5*(elMat(row,col) + elMat(col,row));
							double* a;

							a = (*this)(reqno,ceqno);
							if (a)
								*a += elMat(row,col);
							else
								ExceptionT::OutOfRange(caller);

							a = (*this)(ceqno,reqno);
							if (a)
								*a += elMat(row,col);
							else
								ExceptionT::OutOfRange(caller);
						}
					}
				}
			}
		}
	}
	else if (format == ElementMatrixT::kSymmetricUpper) /* still need to assemble whole matrix */
	{
		int nee = eqnos.Length();  // number of equations for element
		for (int col = 0; col < nee; ++col)
		{
			int ceqno = eqnos[col] - 1;
			if (ceqno > -1) /* active eqn */ {
				for (int row = 0; row <= col; ++row) {
					int reqno = eqnos[row] - 1;
					if (reqno > -1) /* active eqn */ {
						double* a = (*this)(reqno,ceqno);
						if (a)
							*a += elMat(row,col);
						else
							ExceptionT::OutOfRange(caller);
							
						/* transpose */
						if (reqno != ceqno) {
							double* a = (*this)(ceqno,reqno);
							if (a)
								*a += elMat(row,col);
							else
								ExceptionT::OutOfRange(caller);
						}
					}
				}
			}
		}
	}
	else
		ExceptionT::GeneralFail(caller, "unsupported element matrix format %d", format);
}

void PSPASESMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
	const ArrayT<int>& col_eqnos)
{
#pragma unused(elMat)
#pragma unused(row_eqnos)
#pragma unused(col_eqnos)
	ExceptionT::GeneralFail("SLUMatrix::Assemble", "non-square not implemented");
}

void PSPASESMatrixT::Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos)
{
#pragma unused(diagonal_elMat)
#pragma unused(eqnos)
	ExceptionT::GeneralFail("SLUMatrix::Assemble", "diagonal not implemented");
}

/* assignment operator */
GlobalMatrixT& PSPASESMatrixT::operator=(const PSPASESMatrixT& rhs)
{
#pragma unused(rhs)
	ExceptionT::GeneralFail("PSPASESMatrixT::operator=", "not implemented");
	return *this;
}

/* assignment operator */
GlobalMatrixT& PSPASESMatrixT::operator=(const GlobalMatrixT& rhs)
{
	const char caller[] = "PSPASESMatrixT::operator=";

#ifdef __NO_RTTI__
	ExceptionT::GeneralFail(caller, "requires RTTI");
#endif

	const PSPASESMatrixT* sp = TB_DYNAMIC_CAST(const PSPASESMatrixT*, &rhs);
	if (!sp)  ExceptionT::GeneralFail(caller, "cast const PSPASESMatrixT* failed");
	return operator=(*sp);
}

/** return a clone of self */
GlobalMatrixT* PSPASESMatrixT::Clone(void) const
{
	PSPASESMatrixT* new_mat = new PSPASESMatrixT(*this);
	return new_mat;
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* rank check functions */
void PSPASESMatrixT::PrintAllPivots(void) const
{
//not implemented
}

void PSPASESMatrixT::PrintZeroPivots(void) const
{
//not implemented
}

void PSPASESMatrixT::PrintLHS(bool force) const
{
#pragma unused(force)
	PSPASESMatrixT& A = const_cast<PSPASESMatrixT&>(*this);
	
	int d_width = OutputWidth(fOut, favals.Pointer());
	fOut << "\nLHS matrix:\n\n";
	for (int row = 0; row < fLocNumEQ; row++) {
		for (int col = 0; col < fTotNumEQ; col++) {
			double* pa = A(row + fStartEQ - 1, col);
			double a = (pa) ? *pa : 0.0;
			fOut << setw(d_width) << a;
		}
		fOut << '\n';
	}
	fOut << endl;
}

/* precondition matrix */
void PSPASESMatrixT::Factorize(void)
{
	/* MPI communicator */
	MPI_Comm comm = fComm.Comm();

	/* compute symbolic factorization */
	if (!fIsSymFactorized) {
	
		/* compute fill-reducing ordering */
		PSPACEO(frowdist.Pointer(), faptrs.Pointer(), fainds.Pointer(), 
			forder.Pointer(), fsizes.Pointer(), fioptions_PSPACEO.Pointer(), &comm);
       
       	/* compute shape of L by symbolic factorization */
		PSPACEY(frowdist.Pointer(), faptrs.Pointer(), fainds.Pointer(),
			forder.Pointer(), fsizes.Pointer(), fioptions_PSPACEY.Pointer(), fdoptions_PSPACEY.Pointer(), 
			&fYcomm, &comm);

		fIsSymFactorized = 1;
	}

	/* compute numerical factorization */
	if (!fIsFactorized) {

		/* compute numerical factorization */
		DPSPACEN(frowdist.Pointer(), faptrs.Pointer(), fainds.Pointer(),
			favals.Pointer(), &fYcomm, &fNcomm, &comm);

		fIsFactorized = 1;
	}
}
	
/* determine new search direction and put the results in result */
void PSPASESMatrixT::BackSubstitute(dArrayT& result)
{
	const char caller[] = "PSPASESMatrixT::BackSubstitute";

	/* check */
	if (fTotNumEQ != fLocNumEQ)
		ExceptionT::GeneralFail(caller, "total equations (%d) != local equations (%d)",
			fTotNumEQ, fLocNumEQ);

	/* flag should not be set */
	if (fIsFactorized) ExceptionT::GeneralFail(caller);

	/* MPI communicator */
	MPI_Comm comm = fComm.Comm();

	/* compute solution */
	int nrhs = 1;
	fb = result;
	DPSPACET(frowdist.Pointer(), &nrhs, fb.Pointer(), &fLocNumEQ, result.Pointer(),
		&fLocNumEQ, fioptions_DPSPACET.Pointer(), &fNcomm, &comm);

#if 0
	void CHECKB_AX(int *rowdista,int *aptrs,int *ainds,double *avals,
		int *rowdistb,int *pnrhs,double *b,int *pldb,double *x,
		int *pldx,double *perr,MPI_Comm *pcomm);
#endif
}

/* element accessor */
double* PSPASESMatrixT::operator()(int row, int col)
{
	const char caller[] = "PSPASESMatrixT::operator()";
	int row_loc = row - fStartEQ + 1;

#if __option(extended_errorcheck)
	/* range checks */
	if (row_loc < 0 || row_loc >= fLocNumEQ) ExceptionT::OutOfRange(caller);
	if (col < 0 || col >= fTotNumEQ) ExceptionT::OutOfRange(caller);
#endif

	/* equations are 1... */
	int* paptrs = faptrs(row_loc);
	int r_dex = paptrs[0] - 1; /* PSPASES uses 1... */
	int* painds = fainds.Pointer(r_dex);

	/* range */
	int min = 0;
	int max = paptrs[1] - 1; /* last col index in row */

	/* is last value */
	if (painds[max] == col)
		return favals.Pointer(r_dex + max);

	/* bisection */
	int c_dex = (max + min)/2;
	int c_hit = painds[c_dex];
	while (c_hit != col && max != min+1) {
	
		/* shift bounds */
		if (c_hit > col)
			max = c_dex;
		else /* painds[cdex] < col */
			min = c_dex;
	
		/* bisect */
		c_dex = (max + min)/2;
		c_hit = painds[c_dex];
	}
	
	/* found */
	if (c_hit == col)
		return favals.Pointer(r_dex + c_dex);
	else /* not found */
		return NULL;
}

#endif /* __PSPASES__ */
