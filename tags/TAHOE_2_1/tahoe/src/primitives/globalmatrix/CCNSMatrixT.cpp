/* $Id: CCNSMatrixT.cpp,v 1.21 2004-08-03 06:32:39 paklein Exp $ */
/* created: paklein (03/04/1998) */
#include "CCNSMatrixT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include <string.h>

#include "toolboxConstants.h"
#include "dMatrixT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "RaggedArray2DT.h"
#include "ElementMatrixT.h"

using namespace Tahoe;

/* constructor */
CCNSMatrixT::CCNSMatrixT(ostream& out, int check_code):
	GlobalMatrixT(out, check_code),
	famax(NULL),
	fNumberOfTerms(0),
	fMatrix(NULL),
	fu(NULL),
	fIsFactorized(false)
{

}

/* copy constructor */
CCNSMatrixT::CCNSMatrixT(const CCNSMatrixT& source):
	GlobalMatrixT(source),
	famax(NULL),
	fKU(NULL),
	fKL(NULL),
	fKD(NULL),
	fNumberOfTerms(0),
	fMatrix(NULL),
	fu(NULL),
	fIsFactorized(false)
{
	CCNSMatrixT::operator=(source);
}

CCNSMatrixT::~CCNSMatrixT(void)
{
	delete[] famax;
	delete[] fMatrix;
	delete[] fu;
}

/* set the internal matrix structure.
* NOTE: do not call Initialize() equation topology has been set
* with AddEquationSet() for all equation sets */
void CCNSMatrixT::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* inherited */
	GlobalMatrixT::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* check */
	if (tot_num_eq != loc_num_eq)
	{
		cout << "\n CCNSMatrixT::Initialize: expecting total number of equations\n"
		     <<   "     " << tot_num_eq
		     << " to be equal to the local number of equations " << loc_num_eq << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* allocate */	
	if (famax != NULL) delete[] famax;
	iArrayT i_memory;
	try {
		i_memory.Dimension(fLocNumEQ + 1);
		i_memory.ReleasePointer(&famax);
	}	
	catch (ExceptionT::CodeT error)
	{
		if (error == ExceptionT::kOutOfMemory)
		{
			cout << "\n CCNSMatrixT::Initialize: not enough memory" << endl;
			fOut << "\n CCNSMatrixT::Initialize: not enough memory" << endl;
		}
		throw error;
	}	

	/* compute matrix structure and return dimensions */
	int meanband;
	int band;
	ComputeSize(fNumberOfTerms, meanband, band);

	/* output */
	fOut << " Number of terms in global matrix. . . . . . . . = " << fNumberOfTerms << '\n';
	fOut << " Mean half bandwidth . . . . . . . . . . . . . . = " << meanband << '\n';
	fOut << " Bandwidth . . . . . . . . . . . . . . . . . . . = " << band     << '\n';


	/* allocate */	
	if (fu != NULL) delete[] fu;
	if (fMatrix != NULL) delete[] fMatrix;
	dArrayT d_memory;
	try {
		d_memory.Dimension(fLocNumEQ);
		d_memory.ReleasePointer(&fu);
		d_memory.Dimension(0); /* reset */

		d_memory.Dimension(fNumberOfTerms);
		d_memory.ReleasePointer(&fMatrix);
	}	
	catch (ExceptionT::CodeT error)
	{
		if (error == ExceptionT::kOutOfMemory)
		{
			cout << "\n CCNSMatrixT::Initialize: not enough memory" << endl;
			fOut << "\n CCNSMatrixT::Initialize: not enough memory" << endl;
		}
		throw error;
	}	
	
	/* set pointers */
	fKU = fMatrix;
	fKL = fKU + famax[fLocNumEQ];
	fKD = fKL + famax[fLocNumEQ];
	
	int computefilledin = 1;
	if (computefilledin)
	{
		int filledelements = NumberOfFilled();
		double percent_fill = (fNumberOfTerms != 0) ? 
			(100.0*filledelements)/fNumberOfTerms : 
			0.0;

		fOut << " Storage efficiency (% non-zero) . . . . . . . . = ";
		fOut << percent_fill << '\n';
	}
	/* flush stream */
	fOut << endl;

	/* clear stored equation sets */
	fEqnos.Clear();
	fRaggedEqnos.Clear();
	
	/* set flag */
}

/* set all matrix volues to 0.0 */
void CCNSMatrixT::Clear(void)
{
	/* inherited */
	GlobalMatrixT::Clear();

	/* byte set */
	memset(fMatrix, 0, sizeof(double)*fNumberOfTerms);
	fIsFactorized = false;
}

/* add element group equations to the overall topology.
* NOTE: assembly positions (equation numbers) = 1...fLocNumEQ
* equations can be of fixed size (iArray2DT) or
* variable length (RaggedArray2DT) */
void CCNSMatrixT::AddEquationSet(const iArray2DT& eqset)
{
	fEqnos.AppendUnique(&eqset);	
}

void CCNSMatrixT::AddEquationSet(const RaggedArray2DT<int>& eqset)
{
	fRaggedEqnos.AppendUnique(&eqset);	
}

/* assemble the element contribution into the LHS matrix - assumes
* that elMat is square (n x n) and that eqnos is also length n.
* NOTE: assembly positions (equation numbers) = 1...fLocNumEQ */
void CCNSMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos)
{
	/* element matrix format */
	ElementMatrixT::FormatT format = elMat.Format();

	if (format == ElementMatrixT::kDiagonal)
	{
		/* from diagonal only */
		const double* pelMat = elMat.Pointer();
		int inc = elMat.Rows() + 1;
		
		int size = eqnos.Length();
		int eqno;
		for (int eqdex = 0; eqdex < size; eqdex++)
		{
			eqno = eqnos[eqdex] - 1;
			if (eqno > -1)	/* active dof */
			{
#if __option (extended_errorcheck)
				if (eqno < 0 || eqno >= fLocNumEQ)
					ExceptionT::OutOfRange("CCNSMatrixT::Assemble", "bad equation number %d", eqno+1);
#endif
				/* assemble */
				fKD[eqno] += *pelMat;
			}
			
			pelMat += inc;
		}
	}
	else
	{   	
		/* copy to full symmetric */
		if (format == ElementMatrixT::kSymmetricUpper) elMat.CopySymmetric();

		/* assemble active degrees of freedom */
		int size = eqnos.Length();
		for (int col = 0; col < size; col++)
		{
			int ceqno = eqnos[col] - 1;	
			if (ceqno > -1)	
				for (int row = 0; row < size; row++)
				{
					int reqno = eqnos[row] - 1;
					if ( reqno > -1)
						(*this)(reqno,ceqno) += elMat(row,col);
				}
		}
	}
}

void CCNSMatrixT::Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
	const ArrayT<int>& col_eqnos)
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (elMat.Rows() != row_eqnos.Length() ||
	    elMat.Cols() != col_eqnos.Length()) throw ExceptionT::kSizeMismatch;
#endif

	/* element matrix format */
	ElementMatrixT::FormatT format = elMat.Format();

	if (format == ElementMatrixT::kDiagonal)
	{
		cout << "\n CCNSMatrixT::Assemble(m, r, c): cannot assemble diagonal matrix" << endl;
		throw ExceptionT::kGeneralFail;
	}
	else
	{   	
		/* copy to full symmetric */
		if (format == ElementMatrixT::kSymmetricUpper) elMat.CopySymmetric();

		/* assemble active degrees of freedom */
		int n_c = col_eqnos.Length();
		int n_r = row_eqnos.Length();
		for (int col = 0; col < n_c; col++)
		{
			int ceqno = col_eqnos[col] - 1;	
			if (ceqno > -1)	
				for (int row = 0; row < n_r; row++)
				{
					int reqno = row_eqnos[row] - 1;
					if ( reqno > -1)
						(*this)(reqno,ceqno) += elMat(row,col);
				}
		}
	}
}

void CCNSMatrixT::Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos)
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (diagonal_elMat.Length() != eqnos.Length()) ExceptionT::SizeMismatch("CCNSMatrixT::Assemble");
#endif

	for (int i = 0; i < eqnos.Length(); i++)
	{
		int eqno = eqnos[i] - 1;
		if (eqno > -1)
			fKD[eqno] += diagonal_elMat[i];
	}
}

/* returns 1 if the factorized matrix contains a negative
* pivot.  Matrix MUST be factorized.  Otherwise function
* returns 0 */
int CCNSMatrixT::HasNegativePivot(void) const
{
	if (!fIsFactorized) return 0;

	for (int i = 0; i < fLocNumEQ; i++)
		if (fKD[i] < 0.0)
			return 1;
	
	return 0;
}

/* number scope and reordering */
GlobalMatrixT::EquationNumberScopeT CCNSMatrixT::EquationNumberScope(void) const
{
	return kLocal;
}

bool CCNSMatrixT::RenumberEquations(void) const { return true; }

/* find the smallest and largest diagonal value */
void CCNSMatrixT::FindMinMaxPivot(double& min, double& max, double& abs_min, 
	double& abs_max) const
{
	if (fLocNumEQ == 0) min = max = abs_min = abs_max = 0.0;
	else
	{
		abs_min = abs_max = min = max = fKD[0];
		for (int i = 1; i < fLocNumEQ; i++)
		{
			double& diag = fKD[i];

			/* absolute */
			if (diag < min)
				min = diag;
			else if (diag > max)
				max = diag;
				
			/* magnitude */	
			if (fabs(diag) < fabs(abs_min))
				abs_min = diag;
			else if (fabs(diag) > fabs(abs_max))
				abs_max = diag;
		}
	}
}

/* assignment operator */
GlobalMatrixT& CCNSMatrixT::operator=(const CCNSMatrixT& rhs)
{
	fIsFactorized = rhs.fIsFactorized;

	cout << "\n CCNSMatrixT::operator= : not implemented" << endl;
	throw ExceptionT::kGeneralFail;
}

/* assignment operator */
GlobalMatrixT& CCNSMatrixT::operator=(const GlobalMatrixT& rhs)
{
#ifdef __NO_RTTI__
	cout << "\n CCNSMatrixT::operator= : requires RTTI" << endl;
	throw ExceptionT::kGeneralFail;
#endif

	const CCNSMatrixT* ccns = TB_DYNAMIC_CAST(const CCNSMatrixT*, &rhs);
	if (!ccns) {
		cout << "\n CCNSMatrixT::operator= : cast failed" << endl;
		throw ExceptionT::kGeneralFail;
	}
	return operator=(*ccns);
}

/* return a clone of self. Caller is responsible for disposing of the matrix */
GlobalMatrixT* CCNSMatrixT::Clone(void) const
{
	CCNSMatrixT* new_mat = new CCNSMatrixT(*this);
	return new_mat;
}

/* return the values along the diagonal of the matrix */
bool CCNSMatrixT::CopyDiagonal(dArrayT& diags) const
{
	/* cannot be factorized */
	if (fIsFactorized)
		return false;
	else
	{
		/* copy */
		dArrayT tmp(fLocNumEQ, fKD);
		diags = tmp;
		return true;		
	}
}

/**************************************************************************
* Protected
**************************************************************************/

/* element accessor */
double CCNSMatrixT::Element(int row, int col) const
{
	/* in the skyline */
	if (InSkyline(row,col))
		return (*this)(row,col);
	else
		return 0.0;
}

namespace Tahoe {

ostream& operator<<(ostream& out, const CCNSMatrixT& matrix)
{
	double* junk = NULL;
	int d_width = OutputWidth(out, junk);

	for (int row = 0; row < matrix.fLocNumEQ; row++)
	{
		for (int col = 0; col < matrix.fLocNumEQ; col++)
			out << setw(d_width) << matrix.Element(row,col);
	
		out << '\n';
	}

	return out;
}

}

/* solution routines */
void CCNSMatrixT::Factorize(void)
{
	/* quick exit */
	if (fIsFactorized)
		return;
	else /* factorization routine (GRF) */ {
		SolNonSymSysSkyLine(fKU, fKL, fKD, famax, fLocNumEQ);
		fIsFactorized = true;
	}
}

/* solves system. K has already been decomposed in LU form. */
void CCNSMatrixT::BackSubstitute(dArrayT& result)
{
	if (!fIsFactorized) ExceptionT::GeneralFail("CCNSMatrixT::BackSubstitute",
		"matrix is not factorized");

	/* solves L u'= F -> F := u' (GRF) */
	solvLT(fKL, result.Pointer(), famax, fLocNumEQ);
	
	/* solves U u = u' = F       (GRF) */
	solvUT(fKU, fKD, fu, result.Pointer(), famax, fLocNumEQ);

	/* copy solution back into result */
	memcpy(result.Pointer(),fu,sizeof(double)*fLocNumEQ);
}

/* check functions */
void CCNSMatrixT::PrintZeroPivots(void) const
{
	if (fCheckCode != GlobalMatrixT::kZeroPivots) return;
	int d_width = OutputWidth(fOut, fKD);

	/* pivot extrema */
	double min, max, abs_min, abs_max;
	FindMinMaxPivot(min, max, abs_min, abs_max);
	fOut << "\n Matrix pivots:\n"
	     <<   "     min = " << setw(d_width) << min << '\n' 
	     <<   "     max = " << setw(d_width) << max << '\n' 
	     <<   "   |min| = " << setw(d_width) << abs_min << '\n' 
	     <<   "   |max| = " << setw(d_width) << abs_max << '\n';

	/* write zero or negative pivots */
	int firstline = 1;
	for (int i = 0; i < fLocNumEQ; i++)
	{
		double pivot = fKD[i];
		
		if (pivot < kSmall)
		{
			if (firstline)
			{
				fOut << "\nZero or negative pivots:\n\n";
				firstline = 0;
			}
		
			fOut << setw(kIntWidth) << i + 1;
			fOut << setw(d_width)   << pivot << '\n';
		}
	}
	
	if (!firstline)
		fOut << '\n';
}

void CCNSMatrixT::PrintAllPivots(void) const
{
	if (fCheckCode != GlobalMatrixT::kAllPivots) return;
	int d_width = OutputWidth(fOut, fKD);

	fOut << "\nAll pivots:\n\n";
	
	for (int i = 0; i < fLocNumEQ; i++)
	{
		fOut << setw(kIntWidth) << i + 1;
		fOut << setw(d_width) << fKD[i] << '\n';
	}
	fOut << '\n';
}

void CCNSMatrixT::PrintLHS(bool force) const
{
	if (!force && fCheckCode != GlobalMatrixT::kPrintLHS)
		return;
		
	fOut << "\nLHS matrix:\n\n";
	fOut << (*this) << "\n\n";
}

/* test if {row,col} is within the skyline */
int CCNSMatrixT::InSkyline(int row, int col) const
{
	/* range checks */
	if (row < 0 || row >= fLocNumEQ) throw ExceptionT::kGeneralFail;
	if (col < 0 || col >= fLocNumEQ) throw ExceptionT::kGeneralFail;

	if (row == col)      /* element on diagonal */
		return 1;
	else if (row > col)  /* element in lower triangle */
	{
		int bandwidth = BandWidth(row);
		int offset    = row - col;
		
		if (offset > bandwidth)
			return 0;
		else
			return 1;
	}
	else                 /* element in upper triangle */
	{
		int bandwidth = BandWidth(col);
		int offset    = col - row;
		
		if (offset > bandwidth)
			return 0;
		else
			return 1;
	}
}

/* element accessor - for assembly */
double& CCNSMatrixT::operator()(int row, int col) const
{
	const char caller[] = "CCNSMatrixT::operator()";

	/* range checks */
	if (row < 0 || row >= fLocNumEQ) ExceptionT::OutOfRange(caller);
	if (col < 0 || col >= fLocNumEQ) ExceptionT::OutOfRange(caller);

	if (row == col)      /* element on diagonal */
		return fKD[row];
	else if (row > col)  /* element in lower triangle */
	{
		int bandwidth = BandWidth(row);
		int offset    = row - col;
		
		/* over skyline */
		if (offset > bandwidth) ExceptionT::OutOfRange(caller);

		/* row major storage below the diagonal */
		return fKL[famax[row] + bandwidth - offset];
	}
	else                 /* element in upper triangle */
	{
		int bandwidth = BandWidth(col);
		int offset    = col - row;
		
		/* over skyline */
		if (offset > bandwidth) ExceptionT::OutOfRange(caller);
		
		/* col major storage above the diagonal */
		return fKU[famax[col] + bandwidth - offset];
	}
}

/**************************************************************************
* Private
**************************************************************************/

/* (re-) compute the matrix structure and return the bandwidth
* and mean bandwidth */
void CCNSMatrixT::ComputeSize(int& num_nonzero, int& mean_bandwidth, int& bandwidth)
{
	/* clear diags/columns heights */
	fNumberOfTerms = 0;
	for (int i = 0; i < fLocNumEQ; i++)
		famax[i] = 0;
		
	/* compute column heights */
	const iArray2DT* peq;
	fEqnos.Top();
	while (fEqnos.Next(peq))
		SetSkylineHeights(*peq);		

	const RaggedArray2DT<int>* prageq;
	fRaggedEqnos.Top();
	while (fRaggedEqnos.Next(prageq))
		SetSkylineHeights(*prageq);		

	/* skyline indices */
	mean_bandwidth = 0;
	bandwidth = 0;
	famax[0]  = 0; //first equation has no elements in fKU or fKS
	if (fLocNumEQ > 1)
	{
		int last_h = 0;
		for (int eq = 1; eq < fLocNumEQ; eq++)
		{
			/* find max column height */
			bandwidth = (famax[eq] > bandwidth) ? famax[eq] : bandwidth;
		
			int temp_h = famax[eq];
			famax[eq]  = famax[eq-1] + last_h;
			last_h = temp_h;
		}

		/* size of tri blocks stored at end */
		famax[fLocNumEQ] = famax[fLocNumEQ-1] + last_h;
		num_nonzero = 2*famax[fLocNumEQ] + fLocNumEQ;

		/* final dimensions */
		mean_bandwidth = int(ceil(double(num_nonzero)/fLocNumEQ));
		bandwidth += 1;
	}
}

/* computes the column heights for the given equation list */
void CCNSMatrixT::SetSkylineHeights(const iArray2DT& eqnos)
{
	int nel = eqnos.MajorDim(); /* number of elements */
	int nee = eqnos.MinorDim(); /* number of element equations */

	for (int j = 0; j < nel; j++)
	{
		const int* eleqnos = eqnos(j);
		int  min = fLocNumEQ;
	
		/* find the smallest eqno > 0 */
		for (int k = 0; k < nee; k++)
		{
			int eq = eleqnos[k];
		
			if (eq > 0 && eq < min)
				min = eq;
		}
	
		/* set skyline height */
		for (int i = 0; i < nee; i++)
		{
			int eq = eleqnos[i];
		
			if (eq > 0)
			{
				int height = eq - min; /* this is the number of elements */
									   /* ABOVE the diagonal */
			
				int& currheight = famax[--eq];
				
				if (height > currheight)
					currheight = height;
			}
		}
	}
}

/* computes the column heights for the given equation list */
void CCNSMatrixT::SetSkylineHeights(const RaggedArray2DT<int>& eqnos)
{
	/* loop over elements */
	int nel = eqnos.MajorDim();
	for (int j = 0; j < nel; j++)
	{
		int nee = eqnos.MinorDim(j);
		const int* eleqnos = eqnos(j);
	
		/* find the smallest eqno > 0 */
		int min = fLocNumEQ;
		for (int k = 0; k < nee; k++)
		{
			int eq = eleqnos[k];
		
			if (eq > 0 && eq < min)
				min = eq;
		}
	
		/* set column height */
		for (int i = 0; i < nee; i++)
		{
			int eq = eleqnos[i];
		
			if (eq > 0)
			{
				int height = eq - min; /* this is the number of elements */
									   /* ABOVE the diagonal */
			
				int& currheight = famax[--eq];
				
				if (height > currheight)
					currheight = height;
			}
		}
	}
}

/* returns number of non-zero elements */
int CCNSMatrixT::NumberOfFilled(void)
{
	/* make all 0.0 */
	Clear();

	const iArray2DT* peq;
	fEqnos.Top();
	while (fEqnos.Next(peq))
		FillWithOnes(*peq);		

	const RaggedArray2DT<int>* prageq;
	fRaggedEqnos.Top();
	while (fRaggedEqnos.Next(prageq))
		FillWithOnes(*prageq);		

	/* count number non-zero */
	int count = 0;
	double* p = fMatrix;
	for (int i = 0; i < fNumberOfTerms; i++)
		if (*p++ > 0.0)
			count++;

	return count;
}

/* place 1.0 in elements that will be filled */
void CCNSMatrixT::FillWithOnes(const iArray2DT& eqnos)
{
	int nel = eqnos.MajorDim(); /* number of elements */
	int nee = eqnos.MinorDim(); /* number of element equations */

	/* make element stiffness filled with ones */
	ElementMatrixT elstiff(nee, ElementMatrixT::kSymmetric);
	elstiff = 1.0;

	iArrayT localeqnos;
	for (int i = 0; i < nel; i++)
	{
		eqnos.RowAlias(i,localeqnos);
		Assemble(elstiff, localeqnos);
	}
}

/* place 1.0 in elements that will be filled */
void CCNSMatrixT::FillWithOnes(const RaggedArray2DT<int>& eqnos)
{
	/* number of elements */
	int nel = eqnos.MajorDim();

	/* make element stiffness filled with ones */
	int maxlength = eqnos.MaxMinorDim();
	dArrayT space(maxlength*maxlength);
	space = 1.0;

	ElementMatrixT elstiff(ElementMatrixT::kSymmetric);

	iArrayT localeqnos;
	for (int i = 0; i < nel; i++)
	{
		int nee = eqnos.MinorDim(i);
		elstiff.Set(nee, nee, space.Pointer());
	
		eqnos.RowAlias(i,localeqnos);
		Assemble(elstiff, localeqnos);
	}
}

/* Decompose K in LU form */
void CCNSMatrixT::SolNonSymSysSkyLine( double* KS, double* KI, double* KD,
	int* maxa, int neq)
{
	int it, firstp, fpvec, firstpscalp, len, i, posd1, posd2, pos1, pos2, pos, posd;
	double *vecLI, *vecCI, *vecLS, *vecCS, *vecLD, *vecCD;
	double sumI, sumS, sumD;

	int fail = 0;
	for( it = 1; it < neq; it++)
	{
		firstp = it - (maxa[it+1] - maxa[it]);
		posd1 = maxa[it + 1];

		for( i = firstp; i < it; i++)
		{
			fpvec = i - (maxa[i+1] - maxa[i]);
			firstpscalp = (fpvec > firstp) ? fpvec:firstp;
			len = (i - 1) - firstpscalp + 1;

			posd2 = maxa[ i + 1];
			pos1  = it - firstpscalp;
			pos2  =  i - firstpscalp;

			vecLI = &KI[posd1 - pos1];
			vecCI = &KS[posd2 - pos2];
			vecLS = &KI[posd2 - pos2];
			vecCS = &KS[posd1 - pos1];

			sumI = Dot( vecLI, vecCI, len);
			sumS = Dot( vecLS, vecCS, len);

			pos = posd1 - (it - i);

			/* skipped with zero pivot */
			if (KD[i] != 0.0)
				KI[pos] = (KI[pos] - sumI) / KD[i];
			else
				fail = 1;
				
			KS[pos] -= sumS;
		}

		len = (it - 1) - firstp + 1;
		posd = posd1 - (it - firstp);

		vecLD = &KI[posd];
		vecCD = &KS[posd];

		sumD = Dot( vecLD, vecCD, len);
		KD[it] -= sumD;
	}

	if (fail)
	{
		cout << "\n CCNSMatrixT::SolNonSymSysSkyLine: factorization is approximate due to zero";
		cout << " values on the diagonal" << endl;
	}
}

void CCNSMatrixT::solvLT(double* KI, double* F, int* maxa,int neq)
{
	int i, firstp, len;
	double *vecA, *vecB;
	double sum;

	for( i = 0; i < neq; i++)
	{
		firstp = i - (maxa[i+1] - maxa[i]);
		vecA = &KI[maxa[i]];
		vecB =  &F[firstp];
		len = (i - 1) - firstp + 1;
		sum = Dot(vecA, vecB, len);
		F[i] -= sum;
	}
}

void CCNSMatrixT::solvUT(double* KS, double* KD, double* u,
double* F, int* maxa, int neq)
{
	int i, j, it, firstp;

	int fail = 0;
	if (KD[neq - 1] != 0.0)
		u[neq - 1] = F[neq - 1] / KD[neq - 1];
	else
		fail = 1;
		
	for( i = neq - 2; i > -1; i--)
	{
		it = i + 1;
		firstp = it - (maxa[it + 1] - maxa[it]);
		for( j = firstp; j < it; j++)
			F[j] -= u[it] * KS[maxa[it + 1] - it + j];
			
		if (KD[i] != 0.0)	
			u[i] = F[i] / KD[i];
		else
			fail = 1;
	}

	if (fail)
	{
		cout << "\n CCNSMatrixT::BackSubstitute: solution is approximate due to zero";
		cout << " values on the diagonal" << endl;
	}
}
