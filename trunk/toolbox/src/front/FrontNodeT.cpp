/* $Id: FrontNodeT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (03/19/1999)                                          */

#include "FrontNodeT.h"

#include <math.h>

#include "dMatrixT.h"

/* constants */
const double Pi = acos(-1.0);

/* vector functions */
inline static void CrossProduct(const double* A, const double* B, double* AxB)
{   AxB[0] = A[1]*B[2] - A[2]*B[1];
	AxB[1] = A[2]*B[0] - A[0]*B[2];
	AxB[2] = A[0]*B[1] - A[1]*B[0];
};

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
const bool ArrayT<FrontNodeT*>::fByteCopy = true;

/* constructor */
FrontNodeT::FrontNodeT(int dim, const double* x, const double* v_n,
	const double* v_t, double cone, double da, int num_pts)
{
	Reset(dim, x, v_n, v_t, cone, da, num_pts);
}
	
/* write point data to output */
void FrontNodeT::Write(ostream& out) const
{
	/* coordinates */
	out << setw(kDoubleWidth) << fx[0];
	out << setw(kDoubleWidth) << fx[1];
	if (fdim == 3)
		out << setw(kDoubleWidth) << fx[2];
	out << '\n';

	/* dimension */
	out << ftip_pts.MajorDim() << '\n';

	/* sampling point coordinates */
	ftip_pts.WriteNumbered(out);
	
	/* transformation tensors */
	dMatrixT Q;
	for (int i = 0; i < fQ.MajorDim(); i++)
	{
		Q.Set(fdim, fdim, fQ(i));
		out << Q << '\n';
	}
}

/* reset functions */
void FrontNodeT::Reset2D(const double* x, const double* v_n, const double* v_t,
	double cone, double da, int num_pts)
{
#pragma unused(v_t) // only needed for 3D

	/* copy coordinates */
	fx[0] = x[0];
	fx[1] = x[1];
	fx[2] = 0.0;

	/* allocate space - points symmetric about straight ahead */
	ftip_pts.Allocate(2*num_pts + 1, 2);
	fQ.Allocate(2*num_pts + 1, 2*2);

	/* relative angle */
	double t0 = atan2(v_n[1], v_n[0]);

	cone *= (Pi/180.0); // to radians
	double  t = (num_pts == 0) ? 0.0 : t0 - cone;
	double dt = (num_pts == 0) ? 0.0 : cone/num_pts;
	dMatrixT Q;
	for (int i = 0; i < ftip_pts.MajorDim(); i++)
	{
		/* set alias */
		Q.Set(2, 2, fQ(i));
	
		/* transform in local frame */
		Q(0,0) = cos(t); Q(0,1) =-sin(t);
		Q(1,0) = sin(t); Q(1,1) = cos(t);
	
		/* sampling point */
		double * x_s = ftip_pts(i);
		x_s[0] = fx[0] + da*Q(0,0);
		x_s[1] = fx[1] + da*Q(1,0);
		
		/* next */
		t += dt;
	}
}

void FrontNodeT::Reset3D(const double* x, const double* v_n, const double* v_t,
	double cone, double da, int num_pts)
{
	/* copy coordinates */
	fx[0] = x[0];
	fx[1] = x[1];
	fx[2] = x[2];
	
	/* directions */
	double N_n[3];
	N_n[0] = v_n[0];
	N_n[1] = v_n[1];
	N_n[2] = v_n[2];
	Normalize(N_n);

	double N_t[3];
	N_t[0] = v_t[0];
	N_t[1] = v_t[1];
	N_t[2] = v_t[2];
	Normalize(N_t);

	/* local coordinate frame */
	CrossProduct(N_n, N_t, fN_nt);
	dMatrixT Q_0(3);
	Q_0(0,0) = N_n[0];
	Q_0(1,0) = N_n[1];
	Q_0(2,0) = N_n[2];

	Q_0(0,1) = N_t[0];
	Q_0(1,1) = N_t[1];
	Q_0(2,1) = N_t[2];

	Q_0(0,2) = fN_nt[0];
	Q_0(1,2) = fN_nt[1];
	Q_0(2,2) = fN_nt[2];

#if __option(extended_errorcheck)
	double det = Q_0.Det();
	if ((det - 1.0) > 1.0e-10)
	{
		cout << "\n FrontNodeT::Reset: local coordinate basis error: det = " << det;
		throw eGeneralFail;
	}
#endif	

	/* allocate space - points symmetric about straight ahead */
	ftip_pts.Allocate(2*num_pts + 1, 3);
	fQ.Allocate(2*num_pts + 1, 3*3);

	/* temp space for coordinate transforms */
	dMatrixT Q_x(3); // in local frame
	dMatrixT Q_X;    // in global frame
	Q_x = 0.0;
	Q_x(1,1) = 1.0;  // rotation axis
	
	cone *= (Pi/180.0); // to radians
	double  t = (num_pts == 0) ? 0.0 : -cone;
	double dt = (num_pts == 0) ? 0.0 : cone/num_pts;
	for (int i = 0; i < ftip_pts.MajorDim(); i++)
	{
		/* transform in local frame */
		Q_x(0,0) = cos(t); Q_x(0,2) =-sin(t);
		Q_x(2,0) = sin(t); Q_x(2,2) = cos(t);
	
		/* to global frame */
		Q_X.Set(3, 3, fQ(i));
		Q_X.MultAB(Q_0, Q_x);
	
		/* sampling point */
		double * x_s = ftip_pts(i);
		x_s[0] = fx[0] + da*Q_X(0,0);
		x_s[1] = fx[1] + da*Q_X(1,0);
		x_s[2] = fx[2] + da*Q_X(2,0);
		
		/* next */
		t += dt;
	}
}
