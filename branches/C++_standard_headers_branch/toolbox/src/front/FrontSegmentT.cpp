/* $Id: FrontSegmentT.cpp,v 1.7 2003-11-04 01:21:14 paklein Exp $ */
/* created: paklein (03/19/1999) */
#include "FrontSegmentT.h"
#include "ArrayT.h"
#include <math.h>

using namespace Tahoe;

/* vector functions */
inline static void CrossProduct(const double* A, const double* B, double* AxB)
{   AxB[0] = A[1]*B[2] - A[2]*B[1];
	AxB[1] = A[2]*B[0] - A[0]*B[2];
	AxB[2] = A[0]*B[1] - A[1]*B[0]; };

inline static double Dot(const double* A, const double* B)
{ return A[0]*B[0] + A[1]*B[1] + A[2]*B[2]; };

inline static void Normalize(double* A)
{
	double length =  sqrt(Dot(A,A));
	A[0] /= length;	
	A[1] /= length;	
	A[2] /= length;	
};

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<FrontSegmentT*>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
FrontSegmentT::FrontSegmentT(const double* A, const double* B, const double* C)
{
	Reset(A, B, C);
}
FrontSegmentT::FrontSegmentT(const double* A, const double* B,
	const FrontSegmentT& source)
{
	Reset(A, B, source);
}	

/* reset segment coordinates */
void FrontSegmentT::Reset(const double* A, const double* B, const double* C)
{
	/* copy point data */
	fA[0] = A[0];
	fA[1] = A[1];
	fA[2] = A[2];

	fB[0] = B[0];
	fB[1] = B[1];
	fB[2] = B[2];

	/* facet edges */
	fN_t[0] = fB[0] - fA[0];
	fN_t[1] = fB[1] - fA[1];
	fN_t[2] = fB[2] - fA[2];
	Normalize(fN_t);

	double AC[3];	
	AC[0] = C[0] - fA[0];
	AC[1] = C[1] - fA[1];
	AC[2] = C[2] - fA[2];
	
	/* facet normal */
	double N_z[3];
	CrossProduct(AC, fN_t, N_z);
	
	/* "outward" */
	CrossProduct(N_z, fN_t, fN_n);
	Normalize(fN_n);
}

void FrontSegmentT::Reset(const double* A, const double* B, const FrontSegmentT& source)
{
	/* endpoints */
	fA[0] = A[0];
	fA[1] = A[1];
	fA[2] = A[2];
	
	fB[0] = B[0];
	fB[1] = B[1];
	fB[2] = B[2];

	/* directions */
	fN_n[0] = source.fN_n[0];
	fN_n[1] = source.fN_n[1];
	fN_n[2] = source.fN_n[2];

	fN_t[0] = source.fN_t[0];
	fN_t[1] = source.fN_t[1];
	fN_t[2] = source.fN_t[2];
}

/* segment length */
double FrontSegmentT::Length(void) const
{
	double dx = fB[0] - fA[0];
	double dy = fB[1] - fA[1];
	double dz = fB[2] - fA[2];

	return sqrt(dx*dx + dy*dy + dz*dz);
}

/* compute midpoint */
void FrontSegmentT::MidPoint(double* x) const // x must be length 3
{
	x[0] = 0.5*(fB[0] + fA[0]);
	x[1] = 0.5*(fB[1] + fA[1]);
	x[2] = 0.5*(fB[2] + fA[2]);
}
