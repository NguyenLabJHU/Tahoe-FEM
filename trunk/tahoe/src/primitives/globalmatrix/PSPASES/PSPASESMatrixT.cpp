/* $Id: PSPASESMatrixT.cpp,v 1.13 2004-11-19 08:30:12 paklein Exp $ */
/* created: paklein (09/13/2000) */
#include "PSPASESMatrixT.h"

/* library support options */
#ifdef __PSPASES__

#include "CommunicatorT.h"
#include "ElementMatrixT.h"
#include "MSRBuilderT.h"

#include "cpspases.h"
#ifdef __LINK_DUMMIES__
void PSPACEO(int*, int*, int*, int*, int*, int*, MPI_Comm*) {}
void PSPACEY(int*, int*, int*, int*, int*, int*, double*, long*, MPI_Comm*) {}
void DPSPACEN(int*, int*, int*, double*, long*, long*, MPI_Comm*) {}
void DPSPACEF(int*, int*, int*, double*, int*, double*, long*, MPI_Comm*) {}
void DPSPACET(int*, int*, double*, int*, double*, int*, int*, long*, MPI_Comm*) {}
void PSPACEC(long*, int*) {}
void CHECKB_AX(int*, int*, int*, double*, int*, int*, double*, int*, double*, int*, double*, MPI_Comm*) {}
#endif

using namespace Tahoe;

/* constuctor */
PSPASESMatrixT::PSPASESMatrixT(ostream& out, int check_code, CommunicatorT& comm):
	GlobalMatrixT(out, check_code),
	fComm(comm),
	fBuilder(NULL),
	faptrs_man(10, faptrs, 2),
	fIsSymFactorized(false),
	fIsNumFactorized(false)
{
	const char caller[] = "PSPASESMatrixT::PSPASESMatrixT";
	fBuilder = new MSRBuilderT(false);
	if (!fBuilder) ExceptionT::OutOfMemory(caller);
}

PSPASESMatrixT::PSPASESMatrixT(const PSPASESMatrixT& source):
	GlobalMatrixT(source),
	fComm(source.fComm),
	fIsSymFactorized(source.fIsSymFactorized),
	fIsNumFactorized(source.fIsNumFactorized)
{
	ExceptionT::GeneralFail("PSPASESMatrixT::PSPASESMatrixT", "not implemented");
}

/* destructor */
PSPASESMatrixT::~PSPASESMatrixT(void)
{
	delete fBuilder;

	/* free storage (order matters) */
	int option_0 = 0;
	if (fIsNumFactorized) PSPACEC(&fNcomm, &option_0);
	if (fIsSymFactorized) PSPACEC(&fYcomm, &option_0);
}

/* add to structure */
void PSPASESMatrixT::AddEquationSet(const iArray2DT& eqnos) { fBuilder->AddGroup(eqnos); }
void PSPASESMatrixT::AddEquationSet(const RaggedArray2DT<int>& eqnos) { fBuilder->AddGroup(eqnos); }

/* set the internal matrix structure */
void PSPASESMatrixT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* inherited - initialize MSR data */
	GlobalMatrixT::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* verify that number of processes is power of 2 */
	int size = fComm.Size();
	int power2 = 2;
	while (size != power2) {
		if (power2 > size)
			ExceptionT::GeneralFail("PSPASESMatrixT::Initialize", "PSPACES requires nproc as 2^n with n >= 1");
		power2 *= 2;
	}

	/* free space */
	int option_0 = 0;
	if (fIsNumFactorized) {
		PSPACEC(&fNcomm, &option_0);
		fIsNumFactorized = false;
	}
	if (fIsSymFactorized) {
		PSPACEC(&fYcomm, &option_0);
		fIsSymFactorized = false;
	}

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
	if (fIsNumFactorized) {
		PSPACEC(&fNcomm, &option_0);	
		fIsNumFactorized = 0;
	}
}

/* assemble the element contribution */
void PSPASESMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos)
{
	const char caller[] = "PSPASESMatrixT::Assemble";

	/* element matrix format */
	ElementMatrixT::FormatT format = elMat.Format();

	/* two cases: element matrix is diagonal, or it's not. */
	int end_update = fStartEQ + fLocNumEQ - 1;
	if (format == ElementMatrixT::kDiagonal)
	{
		/* diagonal entries only */
		const double *pelMat = elMat.Pointer();
		int inc = elMat.Rows() + 1; /* offset between diag entries are */
		int nee = eqnos.Length();
		for (int i = 0; i < nee; ++i) {
			int eq = eqnos[i];
			if (eq >= fStartEQ && eq <= end_update) /* active eqn */ {
				eq--;
				double* a = (*this)(eq,eq);
				if (a)
					*a += *pelMat;
				else
					ExceptionT::OutOfRange(caller);
			}
			pelMat += inc;
		}
	}
	else if (format == ElementMatrixT::kNonSymmetric || 
             format == ElementMatrixT::kSymmetric ||
             format == ElementMatrixT::kSymmetricUpper )
	{
		/* fill matrix */
		if (format != ElementMatrixT::kNonSymmetric)
			elMat.CopySymmetric();

		int nee = eqnos.Length();  // number of equations for element
		for (int col = 0; col < nee; ++col)
		{
			int ceqno = eqnos[col] - 1;
			if (ceqno > -1) /* active eqn */ {
				for (int row = 0; row < nee; ++row) {
					int reqno = eqnos[row];
					if (reqno >= fStartEQ && reqno <= end_update) /* active eqn */ {
						reqno--;
						double* a = (*this)(reqno,ceqno);
						if (a)
							*a += 0.5*(elMat(row,col) + elMat(col,row));
						else
							ExceptionT::OutOfRange(caller);
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
GlobalMatrixT& PSPASESMatrixT::operator=(const GlobalMatrixT& rhs)
{
	const char caller[] = "PSPASESMatrixT::operator=";

#ifdef __NO_RTTI__
	ExceptionT::GeneralFail(caller, "requires RTTI");
#endif

	const PSPASESMatrixT* sp = TB_DYNAMIC_CAST(const PSPASESMatrixT*, &rhs);
	if (!sp)  ExceptionT::GeneralFail(caller, "cast const PSPASESMatrixT* failed");
	return *this;
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
	if (!force && fCheckCode != GlobalMatrixT::kPrintLHS) return;

	fOut << "LHS: {r, c, v}: \n";
	for (int i = 0; i < fLocNumEQ; i++) {
	
		int index = faptrs(i,0) - 1;
		int count = faptrs(i,1);

		const int* col = fainds.Pointer(index);
		const double* val = favals.Pointer(index);
	
		for (int j = 0; j < count; j++)
			fOut << i+fStartEQ << " " << (*col++)+1 << " " << *val++ << '\n';
	}
	fOut << endl;
}

/* precondition matrix */
void PSPASESMatrixT::Factorize(void)
{
	/* quick exit */
	if (fIsNumFactorized) return;

	const char caller[] = "PSPASESMatrixT::Factorize";

	/* MPI communicator */
	MPI_Comm comm = fComm.Comm();

	/* compute symbolic factorization */
	if (!fIsSymFactorized) {

		/* make time stamp */
		fComm.Log(CommunicatorT::kUrgent, caller, "start symbolic factorization");
	
		/* compute fill-reducing ordering */
		PSPACEO(frowdist.Pointer(), faptrs.Pointer(), fainds.Pointer(), 
			forder.Pointer(), fsizes.Pointer(), fioptions_PSPACEO.Pointer(), &comm);
       
       	/* compute shape of L by symbolic factorization */
		PSPACEY(frowdist.Pointer(), faptrs.Pointer(), fainds.Pointer(),
			forder.Pointer(), fsizes.Pointer(), fioptions_PSPACEY.Pointer(), fdoptions_PSPACEY.Pointer(), 
			&fYcomm, &comm);

		/* make time stamp */
		fComm.Log(CommunicatorT::kUrgent, caller, "end symbolic factorization");

		fIsSymFactorized = true;
	}

	/* compute numerical factorization */
	if (!fIsNumFactorized) {

		/* make time stamp */
		fComm.Log(CommunicatorT::kUrgent, caller, "start numerical factorization");

		/* compute numerical factorization */
		DPSPACEN(frowdist.Pointer(), faptrs.Pointer(), fainds.Pointer(),
			favals.Pointer(), &fYcomm, &fNcomm, &comm);

		/* make time stamp */
		fComm.Log(CommunicatorT::kUrgent, caller, "end numerical factorization");

		fIsNumFactorized = true;
	}
}
	
/* determine new search direction and put the results in result */
void PSPASESMatrixT::BackSubstitute(dArrayT& result)
{
	const char caller[] = "PSPASESMatrixT::BackSubstitute";

	/* check */
	if (!fIsNumFactorized) ExceptionT::GeneralFail(caller, "matrix is not factorized");

	/* MPI communicator */
	MPI_Comm comm = fComm.Comm();

	/* compute solution */
	int nrhs = 1;
	fb = result;
	DPSPACET(frowdist.Pointer(), &nrhs, fb.Pointer(), &fLocNumEQ, result.Pointer(),
		&fLocNumEQ, fioptions_DPSPACET.Pointer(), &fNcomm, &comm);

#if __option(extended_errorcheck)
	double error = 0.0;
	CHECKB_AX(frowdist.Pointer(), faptrs.Pointer(), fainds.Pointer(),
		favals.Pointer(), frowdist.Pointer(), &nrhs, fb.Pointer(), &fLocNumEQ, result.Pointer(),
		&fLocNumEQ, &error, &comm);
	fOut << caller << ": max |B - AX| = " << error << endl;
#endif
}

#endif /* __PSPASES__ */
