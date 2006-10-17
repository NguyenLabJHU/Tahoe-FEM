/* $Id: HexahedronT.cpp,v 1.13 2006-10-17 21:17:09 regueiro Exp $ */
/* created: paklein (10/22/1997) */
#include "HexahedronT.h"
#include <math.h>
#include "ExceptionT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"
#include "LocalArrayT.h"

using namespace Tahoe;

/* parameters */
const int kHexnsd         = 3;
const int kNumVertexNodes = 8;
const int kNumFacets      = 6;
const double sqrt3 = sqrt(3.0);
const double sqrt15 = sqrt(15.0);

/* vertex node coordinates */
const double ra[8] = {-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0};
const double sa[8] = {-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0};
const double ta[8] = {-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0};

/* constructor */
HexahedronT::HexahedronT(int numnodes): GeometryBaseT(numnodes, kNumFacets) 
{
	const char caller[] = "HexahedronT::HexahedronT";
	const double ra20[20] = {-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, -1.0, 1.0, 1.0, -1.0};
	const double sa20[20] = {-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, -1.0, 1.0, 1.0};
	const double ta20[20] = {-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};
// I have copied above lines to add 7 nodes/* Davoud */
	const double ra27[27] = {-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0};
	const double sa27[27] = {-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, -1.0, 1.0, 1.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0};
	const double ta27[27] = {-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	fCoords.Dimension(3, numnodes);
	double* x = fCoords(0);
	double* y = fCoords(1);
	double* z = fCoords(2);

	for (int i = 0; i< numnodes; i++)
	{
		/*
		x[i] = ra20[i];
		y[i] = sa20[i];
		z[i] = ta20[i];
		*/
// I have copied above lines to add 7 nodes/* Davoud */
		x[i] = ra27[i];
		y[i] = sa27[i];
		z[i] = ta27[i];
	}
}

/* evaluate the shape functions and gradients. */
void HexahedronT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na) const
{
	const char caller[] = "HexahedronT::EvaluateShapeFunctions";

#if __option(extended_errorcheck)
	if (coords.Length() != 3 ||
	        Na.Length() != fNumNodes) ExceptionT::SizeMismatch(caller);
	if (fNumNodes != kNumVertexNodes && 
	    fNumNodes != 20 && fNumNodes != 27) ExceptionT::GeneralFail(caller);
#endif

	/* coordinates */	
	double r = coords[0];
	double s = coords[1];
	double t = coords[2];
#if 0
	if (r < -1.0 || r > 1.0) ExceptionT::OutOfRange(caller);
	if (s < -1.0 || s > 1.0) ExceptionT::OutOfRange(caller);
	if (t < -1.0 || t > 1.0) ExceptionT::OutOfRange(caller);
#endif

	/* vertex nodes */
	double* na  = Na.Pointer();
	for (int lnd = 0; lnd < kNumVertexNodes; lnd++)
	{
		double tempr1 = 1.0 + ra[lnd]*r;
		double temps1 = 1.0 + sa[lnd]*s;
		double tempt1 = 1.0 + ta[lnd]*t;			
		*na++  = 0.125*tempr1*temps1*tempt1;
	}
	
	/* higher order shape functions */
	if (fNumNodes == 20)
	{
		na  = Na.Pointer();
	
		/* linear factors */
		double r_min = 1.0 - r, r_max = 1.0 + r;
		double s_min = 1.0 - s, s_max = 1.0 + s;
		double t_min = 1.0 - t, t_max = 1.0 + t;

		/* bubbles */
		double r2 = 1.0 - r*r;
		double s2 = 1.0 - s*s;
		double t2 = 1.0 - t*t;
		
		/* local node number */
		int lnd = 8;

		/* node 9 */
		double N = 0.25*r2*s_min*t_min;
		na[lnd] = N;  
		lnd++;
		N *= 0.5;
		/* correct nodes {1,2} */
		na[0] -= N;   
		na[1] -= N; 

		/* node 10 */
		N = 0.25*r_max*s2*t_min;		
		na[lnd] = N; 
		lnd++;
		N *= 0.5;
		/* correct nodes {2,3} */
		na[1] -= N;   
		na[2] -= N; 

		/* node 11 */
		N = 0.25*r2*s_max*t_min;
		na[lnd] = N; 
		lnd++;
		N *= 0.5;
		/* correct nodes {3,4} */
		na[2] -= N;   
		na[3] -= N; 

		/* node 12 */
		N = 0.25*r_min*s2*t_min;
		na[lnd] = N; 
		lnd++;
		N *= 0.5;
		/* correct nodes {4,1} */
		na[3] -= N;   
		na[0] -= N; 

		/* node 13 */
		N = 0.25*r2*s_min*t_max;
		na[lnd] = N; 
		lnd++;
		N *= 0.5;
		/* correct nodes {5,6} */
		na[4] -= N;   
		na[5] -= N; 

		/* node 14 */
		N = 0.25*r_max*s2*t_max;		
		na[lnd] = N; 
		lnd++;
		N *= 0.5;
		/* correct nodes {6,7} */
		na[5] -= N;   
		na[6] -= N; 

		/* node 15 */
		N = 0.25*r2*s_max*t_max;
		na[lnd] = N; 
		lnd++;
		N *= 0.5;
		/* correct nodes {7,8} */
		na[6] -= N;   
		na[7] -= N; 

		/* node 16 */
		N = 0.25*r_min*s2*t_max;
		na[lnd] = N; 
		lnd++;
		N *= 0.5;
		/* correct nodes {8,5} */
		na[7] -= N;   
		na[4] -= N; 

		/* node 17 */
		N = 0.25*r_min*s_min*t2;
		na[lnd] = N;
		lnd++;
		N *= 0.5;
		/* correct nodes {1,5} */
		na[0] -= N;   
		na[4] -= N; 

		/* node 18 */
		N = 0.25*r_max*s_min*t2;
		na[lnd] = N; 
		lnd++;
		N *= 0.5;
		/* correct nodes {2,6} */
		na[1] -= N;   
		na[5] -= N; 

		/* node 19 */
		N = 0.25*r_max*s_max*t2;
		na[lnd] = N; 
		lnd++;
		N *= 0.5; 
		/* correct nodes {3,7} */
		na[2] -= N;   
		na[6] -= N; 

		/* node 20 */
		N = 0.25*r_min*s_max*t2;
		na[lnd] = N; 
		lnd++;
		N *= 0.5;
		/* correct nodes {4,8} */
		na[3] -= N;   
		na[7] -= N; 
	}
	else if (fNumNodes == 27)
	{
// instead of creating shape functions magically(by defining shape functions for the new 7 nodes and correcting prviously defined shape functions)
// they will be defined directly. /* Davoud */
		na  = Na.Pointer();
	
		/* base quadratic factors */
		double r1=0.5*r*(r-1); double s1=0.5*s*(s-1); double t1=0.5*t*(t-1);
		double r2=0.5*r*(r+1); double s2=0.5*s*(s+1); double t2=0.5*t*(t+1);
		double r3=1-r*r;       double s3=1-s*s;       double t3=1-t*t;

		/* local node number */
		int lnd = 0;

		/* node 1 */
		na[lnd] = r1*s1*t1;  
		lnd++;

		/* node 2 */
		na[lnd] = r2*s1*t1; 
		lnd++;

		/* node 3 */
		na[lnd] = r2*s2*t1; 
		lnd++;

		/* node 4 */
		na[lnd] = r1*s2*t1; 
		lnd++;

		/* node 5 */
		na[lnd] = r1*s1*t2; 
		lnd++;

		/* node 6 */
		na[lnd] = r2*s1*t2; 
		lnd++;

		/* node 7 */
		na[lnd] = r2*s2*t2; 
		lnd++;

		/* node 8 */
		na[lnd] = r1*s2*t2; 
		lnd++;

		/* node 9 */
		na[lnd] = r3*s1*t1;
		lnd++;

		/* node 10 */
		na[lnd] = r2*s3*t1; 
		lnd++;

		/* node 11 */
		na[lnd] = r3*s2*t1; 
		lnd++;

		/* node 12 */
		na[lnd] = r1*s3*t1; 
		lnd++;

		/* node 13 */
		na[lnd] = r3*s1*t2; 
		lnd++;

		/* node 14 */
		na[lnd] = r2*s3*t2; 
		lnd++;

		/* node 15 */
		na[lnd] = r3*s2*t2; 
		lnd++;

		/* node 16 */
		na[lnd] = r1*s3*t2; 
		lnd++;

		/* node 17 */
		na[lnd] = r1*s1*t3; 
		lnd++;

		/* node 18 */
		na[lnd] = r2*s1*t3; 
		lnd++;

		/* node 19 */
		na[lnd] = r2*s2*t3; 
		lnd++;

		/* node 20 */
		na[lnd] = r1*s2*t3; 
		lnd++;

		/* node 21 */
		na[lnd] = r3*s3*t1; 
		lnd++;

		/* node 22 */
		na[lnd] = r3*s3*t2; 
		lnd++;

		/* node 23 */
		na[lnd] = r3*s1*t3; 
		lnd++;

		/* node 24 */
		na[lnd] = r3*s2*t3; 
		lnd++;

		/* node 25 */
		na[lnd] = r1*s3*t3; 
		lnd++;

		/* node 26 */
		na[lnd] = r2*s3*t3; 
		lnd++;

		/* node 27 */
		na[lnd] = r3*s3*t3; 
		lnd++;
	}
}

/* evaluate the shape functions and gradients. */
void HexahedronT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na, 
	dArray2DT& DNa) const
{
	const char caller[] = "HexahedronT::EvaluateShapeFunctions";

#if __option(extended_errorcheck)
	if (coords.Length() != 3 ||
	        Na.Length() != fNumNodes ||
	     DNa.MajorDim() != 3 ||
	     DNa.MinorDim() != fNumNodes) ExceptionT::SizeMismatch(caller);
	if (fNumNodes != kNumVertexNodes && 
	    fNumNodes != 20 && fNumNodes != 27) ExceptionT::GeneralFail(caller);
#endif

	/* coordinates */	
	double r = coords[0];
	double s = coords[1];
	double t = coords[2];
#if 0
	if (r < -1.0 || r > 1.0) ExceptionT::OutOfRange(caller, "r = %15.12e", r);
	if (s < -1.0 || s > 1.0) ExceptionT::OutOfRange(caller, "s = %15.12e", s);
	if (t < -1.0 || t > 1.0) ExceptionT::OutOfRange(caller, "t = %15.12e", t);
#endif

	/* vertex nodes */
	double* na  = Na.Pointer();
	double* nax = DNa(0);
	double* nay = DNa(1);
	double* naz = DNa(2);
	for (int lnd = 0; lnd < kNumVertexNodes; lnd++)
	{
		double tempr1 = 1.0 + ra[lnd]*r;
		double temps1 = 1.0 + sa[lnd]*s;
		double tempt1 = 1.0 + ta[lnd]*t;
			
		*na++  = 0.125*tempr1*temps1*tempt1;
		*nax++ = 0.125*ra[lnd]*temps1*tempt1;
		*nay++ = 0.125*tempr1*sa[lnd]*tempt1;
		*naz++ = 0.125*tempr1*temps1*ta[lnd];
	}
	
	/* higher order shape functions */
	if (fNumNodes == 20)
	{
		na  = Na.Pointer();
		nax = DNa(0);
		nay = DNa(1);
		naz = DNa(2);
	
		/* linear factors */
		double r_min = 1.0 - r, r_max = 1.0 + r;
		double s_min = 1.0 - s, s_max = 1.0 + s;
		double t_min = 1.0 - t, t_max = 1.0 + t;

		/* bubbles */
		double r2 = 1.0 - r*r;
		double s2 = 1.0 - s*s;
		double t2 = 1.0 - t*t;
		
		/* local node number */
		int lnd = 8;

		/* node 9 */
		double N  = 0.25*r2*s_min*t_min;
		double Nx =-0.50*r*s_min*t_min;
		double Ny =-0.25*r2*t_min;
		double Nz =-0.25*r2*s_min;
		na[lnd] = N; nax[lnd] = Nx; nay[lnd] = Ny; naz[lnd] = Nz; lnd++;
		N *= 0.5; Nx *= 0.5; Ny *= 0.5; Nz *= 0.5;
		/* correct nodes {1,2} */
		 na[0] -= N;   na[1] -= N; 
		nax[0] -= Nx; nax[1] -= Nx; 
		nay[0] -= Ny; nay[1] -= Ny; 
		naz[0] -= Nz; naz[1] -= Nz;

		/* node 10 */
		N  = 0.25*r_max*s2*t_min;		
		Nx = 0.25*s2*t_min;
		Ny =-0.50*r_max*s*t_min;
		Nz =-0.25*r_max*s2;
		na[lnd] = N; nax[lnd] = Nx; nay[lnd] = Ny; naz[lnd] = Nz; lnd++;
		N *= 0.5; Nx *= 0.5; Ny *= 0.5; Nz *= 0.5;
		/* correct nodes {2,3} */
		 na[1] -= N;   na[2] -= N; 
		nax[1] -= Nx; nax[2] -= Nx; 
		nay[1] -= Ny; nay[2] -= Ny; 
		naz[1] -= Nz; naz[2] -= Nz;

		/* node 11 */
		N  = 0.25*r2*s_max*t_min;
		Nx =-0.50*r*s_max*t_min;
		Ny = 0.25*r2*t_min;
		Nz =-0.25*r2*s_max;
		na[lnd] = N; nax[lnd] = Nx; nay[lnd] = Ny; naz[lnd] = Nz; lnd++;
		N *= 0.5; Nx *= 0.5; Ny *= 0.5; Nz *= 0.5;
		/* correct nodes {3,4} */
		 na[2] -= N;   na[3] -= N; 
		nax[2] -= Nx; nax[3] -= Nx; 
		nay[2] -= Ny; nay[3] -= Ny; 
		naz[2] -= Nz; naz[3] -= Nz;

		/* node 12 */
		N  = 0.25*r_min*s2*t_min;
		Nx =-0.25*s2*t_min;
		Ny =-0.50*r_min*s*t_min;
		Nz =-0.25*r_min*s2;
		na[lnd] = N; nax[lnd] = Nx; nay[lnd] = Ny; naz[lnd] = Nz; lnd++;
		N *= 0.5; Nx *= 0.5; Ny *= 0.5; Nz *= 0.5;
		/* correct nodes {4,1} */
		 na[3] -= N;   na[0] -= N; 
		nax[3] -= Nx; nax[0] -= Nx; 
		nay[3] -= Ny; nay[0] -= Ny; 
		naz[3] -= Nz; naz[0] -= Nz;

		/* node 13 */
		N  = 0.25*r2*s_min*t_max;
		Nx =-0.50*r*s_min*t_max;
		Ny =-0.25*r2*t_max;
		Nz = 0.25*r2*s_min;
		na[lnd] = N; nax[lnd] = Nx; nay[lnd] = Ny; naz[lnd] = Nz; lnd++;
		N *= 0.5; Nx *= 0.5; Ny *= 0.5; Nz *= 0.5;
		/* correct nodes {5,6} */
		 na[4] -= N;   na[5] -= N; 
		nax[4] -= Nx; nax[5] -= Nx; 
		nay[4] -= Ny; nay[5] -= Ny; 
		naz[4] -= Nz; naz[5] -= Nz;

		/* node 14 */
		N  = 0.25*r_max*s2*t_max;		
		Nx = 0.25*s2*t_max;
		Ny =-0.50*r_max*s*t_max;
		Nz = 0.25*r_max*s2;
		na[lnd] = N; nax[lnd] = Nx; nay[lnd] = Ny; naz[lnd] = Nz; lnd++;
		N *= 0.5; Nx *= 0.5; Ny *= 0.5; Nz *= 0.5;		
		/* correct nodes {6,7} */
		 na[5] -= N;   na[6] -= N; 
		nax[5] -= Nx; nax[6] -= Nx; 
		nay[5] -= Ny; nay[6] -= Ny; 
		naz[5] -= Nz; naz[6] -= Nz;

		/* node 15 */
		N  = 0.25*r2*s_max*t_max;
		Nx =-0.50*r*s_max*t_max;
		Ny = 0.25*r2*t_max;
		Nz = 0.25*r2*s_max;
		na[lnd] = N; nax[lnd] = Nx; nay[lnd] = Ny; naz[lnd] = Nz; lnd++;
		N *= 0.5; Nx *= 0.5; Ny *= 0.5; Nz *= 0.5;
		/* correct nodes {7,8} */
		 na[6] -= N;   na[7] -= N; 
		nax[6] -= Nx; nax[7] -= Nx; 
		nay[6] -= Ny; nay[7] -= Ny; 
		naz[6] -= Nz; naz[7] -= Nz;

		/* node 16 */
		N  = 0.25*r_min*s2*t_max;
		Nx =-0.25*s2*t_max;
		Ny =-0.50*r_min*s*t_max;
		Nz = 0.25*r_min*s2;
		na[lnd] = N; nax[lnd] = Nx; nay[lnd] = Ny; naz[lnd] = Nz; lnd++;
		N *= 0.5; Nx *= 0.5; Ny *= 0.5; Nz *= 0.5;
		/* correct nodes {8,5} */
		 na[7] -= N;   na[4] -= N; 
		nax[7] -= Nx; nax[4] -= Nx; 
		nay[7] -= Ny; nay[4] -= Ny; 
		naz[7] -= Nz; naz[4] -= Nz;

		/* node 17 */
		N  = 0.25*r_min*s_min*t2;
		Nx =-0.25*s_min*t2;
		Ny =-0.25*r_min*t2;
		Nz =-0.50*r_min*s_min*t;
		na[lnd] = N; nax[lnd] = Nx; nay[lnd] = Ny; naz[lnd] = Nz; lnd++;
		N *= 0.5; Nx *= 0.5; Ny *= 0.5; Nz *= 0.5;
		/* correct nodes {1,5} */
		 na[0] -= N;   na[4] -= N; 
		nax[0] -= Nx; nax[4] -= Nx; 
		nay[0] -= Ny; nay[4] -= Ny; 
		naz[0] -= Nz; naz[4] -= Nz;

		/* node 18 */
		N  = 0.25*r_max*s_min*t2;
		Nx = 0.25*s_min*t2;
		Ny =-0.25*r_max*t2;
		Nz =-0.50*r_max*s_min*t;
		na[lnd] = N; nax[lnd] = Nx; nay[lnd] = Ny; naz[lnd] = Nz; lnd++;
		N *= 0.5; Nx *= 0.5; Ny *= 0.5; Nz *= 0.5;
		/* correct nodes {2,6} */
		 na[1] -= N;   na[5] -= N; 
		nax[1] -= Nx; nax[5] -= Nx; 
		nay[1] -= Ny; nay[5] -= Ny; 
		naz[1] -= Nz; naz[5] -= Nz;

		/* node 19 */
		N  = 0.25*r_max*s_max*t2;
		Nx = 0.25*s_max*t2;
		Ny = 0.25*r_max*t2;
		Nz =-0.50*r_max*s_max*t;
		na[lnd] = N; nax[lnd] = Nx; nay[lnd] = Ny; naz[lnd] = Nz; lnd++;
		N *= 0.5; Nx *= 0.5; Ny *= 0.5; Nz *= 0.5;
		/* correct nodes {3,7} */
		 na[2] -= N;   na[6] -= N; 
		nax[2] -= Nx; nax[6] -= Nx; 
		nay[2] -= Ny; nay[6] -= Ny; 
		naz[2] -= Nz; naz[6] -= Nz;

		/* node 20 */
		N  = 0.25*r_min*s_max*t2;
		Nx =-0.25*s_max*t2;
		Ny = 0.25*r_min*t2;
		Nz =-0.50*r_min*s_max*t;
		na[lnd] = N; nax[lnd] = Nx; nay[lnd] = Ny; naz[lnd] = Nz; lnd++;
		N *= 0.5; Nx *= 0.5; Ny *= 0.5; Nz *= 0.5;
		/* correct nodes {4,8} */
		 na[3] -= N;   na[7] -= N; 
		nax[3] -= Nx; nax[7] -= Nx; 
		nay[3] -= Ny; nay[7] -= Ny; 
		naz[3] -= Nz; naz[7] -= Nz;
	}
	if (fNumNodes == 27)
	{
// instead of creating derivatives magically(by defining derivatives for the new 7 nodes and correcting prviously defined derivatives)
// they will be defined directly. /* Davoud */
		na  = Na.Pointer();
		nax = DNa(0);
		nay = DNa(1);
		naz = DNa(2);
	
		/* base quadratic factors */
		double r1=0.5*r*(r-1); double s1=0.5*s*(s-1); double t1=0.5*t*(t-1);
		double r2=0.5*r*(r+1); double s2=0.5*s*(s+1); double t2=0.5*t*(t+1);
		double r3=1-r*r;       double s3=1-s*s;       double t3=1-t*t;

		/* base derivatives */
		double dr1=0.5*(2*r-1);double ds1=0.5*(2*s-1);double dt1=0.5*(2*t-1);
		double dr2=0.5*(2*r+1);double ds2=0.5*(2*s+1);double dt2=0.5*(2*t+1);
		double dr3=-2*r;       double ds3=-2*s;       double dt3=-2*t;
		
		/* local node number */
		int lnd = 0;

		/* node 1 */
		na[lnd] = r1*s1*t1;  
		nax[lnd] = dr1*s1*t1;  
		nay[lnd] = r1*ds1*t1;
		naz[lnd] = r1*s1*dt1;
		lnd++;
		
		/* node 2 */
		na[lnd] = r2*s1*t1; 
		nax[lnd] = dr2*s1*t1;  
		nay[lnd] = r2*ds1*t1;
		naz[lnd] = r2*s1*dt1;
		lnd++;

		/* node 3 */
		na[lnd] = r2*s2*t1; 
		nax[lnd] = dr2*s2*t1;  
		nay[lnd] = r2*ds2*t1;
		naz[lnd] = r2*s2*dt1;
		lnd++;

		/* node 4 */
		na[lnd] = r1*s2*t1; 
		nax[lnd] = dr1*s2*t1;  
		nay[lnd] = r1*ds2*t1;
		naz[lnd] = r1*s2*dt1;
		lnd++;

		/* node 5 */
		na[lnd] = r1*s1*t2; 
		nax[lnd] = dr1*s1*t2;  
		nay[lnd] = r1*ds1*t2;
		naz[lnd] = r1*s1*dt2;
		lnd++;

		/* node 6 */
		na[lnd] = r2*s1*t2; 
		nax[lnd] = dr2*s1*t2;  
		nay[lnd] = r2*ds1*t2;
		naz[lnd] = r2*s1*dt2;
		lnd++;

		/* node 7 */
		na[lnd] = r2*s2*t2;
		nax[lnd] = dr2*s2*t2;  
		nay[lnd] = r2*ds2*t2;
		naz[lnd] = r2*s2*dt2; 
		lnd++;

		/* node 8 */
		na[lnd] = r1*s2*t2; 
		nax[lnd] = dr1*s2*t2;  
		nay[lnd] = r1*ds2*t2;
		naz[lnd] = r1*s2*dt2;
		lnd++;

		/* node 9 */
		na[lnd] = r3*s1*t1;
		nax[lnd] = dr3*s1*t1;  
		nay[lnd] = r3*ds1*t1;
		naz[lnd] = r3*s1*dt1;
		lnd++;

		/* node 10 */
		na[lnd] = r2*s3*t1;
		nax[lnd] = dr2*s3*t1;  
		nay[lnd] = r2*ds3*t1;
		naz[lnd] = r2*s3*dt1; 
		lnd++;

		/* node 11 */
		na[lnd] = r3*s2*t1; 
		nax[lnd] = dr3*s2*t1;  
		nay[lnd] = r3*ds2*t1;
		naz[lnd] = r3*s2*dt1;
		lnd++;

		/* node 12 */
		na[lnd] = r1*s3*t1; 
		nax[lnd] = dr1*s3*t1;  
		nay[lnd] = r1*ds3*t1;
		naz[lnd] = r1*s3*dt1;
		lnd++;

		/* node 13 */
		na[lnd] = r3*s1*t2; 
		nax[lnd] = dr3*s1*t2;  
		nay[lnd] = r3*ds1*t2;
		naz[lnd] = r3*s1*dt2;
		lnd++;

		/* node 14 */
		na[lnd] = r2*s3*t2; 
		nax[lnd] = dr2*s3*t2;  
		nay[lnd] = r2*ds3*t2;
		naz[lnd] = r2*s3*dt2;
		lnd++;

		/* node 15 */
		na[lnd] = r3*s2*t2; 
		nax[lnd] = dr3*s2*t2;  
		nay[lnd] = r3*ds2*t2;
		naz[lnd] = r3*s2*dt2;
		lnd++;

		/* node 16 */
		na[lnd] = r1*s3*t2; 
		nax[lnd] = dr1*s3*t2;  
		nay[lnd] = r1*ds3*t2;
		naz[lnd] = r1*s3*dt2;
		lnd++;

		/* node 17 */
		na[lnd] = r1*s1*t3; 
		nax[lnd] = dr1*s1*t3;  
		nay[lnd] = r1*ds1*t3;
		naz[lnd] = r1*s1*dt3;
		lnd++;

		/* node 18 */
		na[lnd] = r2*s1*t3; 
		nax[lnd] = dr2*s1*t3;  
		nay[lnd] = r2*ds1*t3;
		naz[lnd] = r2*s1*dt3;
		lnd++;

		/* node 19 */
		na[lnd] = r2*s2*t3;
		nax[lnd] = dr2*s2*t3;  
		nay[lnd] = r2*ds2*t3;
		naz[lnd] = r2*s2*dt3; 
		lnd++;

		/* node 20 */
		na[lnd] = r1*s2*t3;
		nax[lnd] = dr1*s2*t3;  
		nay[lnd] = r1*ds2*t3;
		naz[lnd] = r1*s2*dt3; 
		lnd++;

		/* node 21 */
		na[lnd] = r3*s3*t1; 
		nax[lnd] = dr3*s3*t1;  
		nay[lnd] = r3*ds3*t1;
		naz[lnd] = r3*s3*dt1;
		lnd++;

		/* node 22 */
		na[lnd] = r3*s3*t2;
		nax[lnd] = dr3*s3*t2;  
		nay[lnd] = r3*ds3*t2;
		naz[lnd] = r3*s3*dt2; 
		lnd++;

		/* node 23 */
		na[lnd] = r3*s1*t3; 
		nax[lnd] = dr3*s1*t3;  
		nay[lnd] = r3*ds1*t3;
		naz[lnd] = r3*s1*dt3;
		lnd++;

		/* node 24 */
		na[lnd] = r3*s2*t3;
		nax[lnd] = dr3*s2*t3;  
		nay[lnd] = r3*ds2*t3;
		naz[lnd] = r3*s2*dt3; 
		lnd++;

		/* node 25 */
		na[lnd] = r1*s3*t3; 
		nax[lnd] = dr1*s3*t3;  
		nay[lnd] = r1*ds3*t3;
		naz[lnd] = r1*s3*dt3;
		lnd++;

		/* node 26 */
		na[lnd] = r2*s3*t3; 
		nax[lnd] = dr2*s3*t3;  
		nay[lnd] = r2*ds3*t3;
		naz[lnd] = r2*s3*dt3;
		lnd++;

		/* node 27 */
		na[lnd] = r3*s3*t3; 
		nax[lnd] = dr3*s3*t3;  
		nay[lnd] = r3*ds3*t3;
		naz[lnd] = r3*s3*dt3;
		lnd++;
	}
}

/* compute local shape functions and derivatives */
void HexahedronT::SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x,
	dArrayT& weights) const
{
	const char caller[] = "HexahedronT::SetLocalShape";

	/* dimensions */
	int numnodes  = Na.MinorDim();
	int numint    = weights.Length();
	int nsd       = Na_x[0].MajorDim();

	/* dimension checks */
	if (numnodes != 8 && numnodes != 20 && numnodes != 27)
		ExceptionT::GeneralFail(caller, "unsupported number of element nodes: %d", numnodes);

	if (numint != 1 &&
	    numint != 8 &&
	    numint != 9 &&
	    numint != 27 &&
	    numint != 64)
	    ExceptionT::GeneralFail(caller, "unsupported number of integration points: %d", numint);
	
	if (nsd != kHexnsd) ExceptionT::GeneralFail(caller);

	/* initialize */
	Na = 0.0;
	for (int i = 0; i < Na_x.Length(); i++)
		Na_x[i] = 0.0;

	/* integration point coordinates */
	double xa_64[64], ya_64[64], za_64[64];
	const double *xa, *ya, *za;
	double 	g;
	
	/* integration weights */
	switch (numint)
	{
		case 1:	
			
		g = 0.0;
		xa = ra;
		ya = sa;
		za = ta;
		weights[0] = 8.0;
		break;

		case 8:
	
		g = 1.0/sqrt3;
		xa = ra;
		ya = sa;
		za = ta;
		weights = 1.0;
		break;

		case 9:
		{
			/* first 8 quadrature points */
			for (int i = 0; i < 8; i++)
			{
				xa_64[i] = ra[i];
				ya_64[i] = sa[i];
				za_64[i] = ta[i];
				weights[i] = 5.0/9.0;
			}
			
			/* the center point */
			xa_64[8] = 0.0;
			ya_64[8] = 0.0;
			za_64[8] = 0.0;
			weights[8] = 32.0/9.0;

			/* set pointers */
			xa = xa_64;
			ya = ya_64;
			za = za_64;
			g  = sqrt(3.0/5.0);
			break;
		}
		case 27:
		{
			/* coordinates */
			double b1 = sqrt(3.0/5.0);
			double b_1D[3] = {-b1, 0.0, b1};
			
			/* weights */
			double w1 = 5.0/9.0;
			double w2 = 8.0/9.0;
			double w_1D[3] = {w1, w2, w1};
			int x_i = 0;
			int y_i = 0;
			int z_i = 0;
			for (int i = 0; i < 27; i++)
			{
				xa_64[i]   = b_1D[x_i];
				ya_64[i]   = b_1D[y_i];
				za_64[i]   = b_1D[z_i];
				weights[i] = w_1D[x_i]*w_1D[y_i]*w_1D[z_i];
							
				if (++x_i == 3)
				{
					x_i = 0;
					if (++y_i == 3)
					{
						y_i = 0;
						z_i++;
					}
				}
			}						
		
			xa = xa_64;
			ya = ya_64;
			za = za_64;
			g  = 1.0;		
			break;
			
		}	
		case 64:
		{
			/* coordinates */
			double b1 = sqrt((3.0 - 2.0*sqrt(6.0/5.0))/7.0);
			double b2 = sqrt((3.0 + 2.0*sqrt(6.0/5.0))/7.0);
			double b_1D[4] = {-b2,-b1, b1, b2};

			/* weights */					
			double w1 = (18.0 + sqrt(30.0))/36.0;
			double w2 = (18.0 - sqrt(30.0))/36.0;
			double w_1D[4] = {w2, w1, w1, w2};

			int x_i = 0;
			int y_i = 0;
			int z_i = 0;
			for (int i = 0; i < 64; i++)
			{
				xa_64[i]   = b_1D[x_i];
				ya_64[i]   = b_1D[y_i];
				za_64[i]   = b_1D[z_i];
				weights[i] = w_1D[x_i]*w_1D[y_i]*w_1D[z_i];
							
				if (++x_i == 4)
				{
					x_i = 0;
					if (++y_i == 4)
					{
						y_i = 0;
						z_i++;
					}
				}
			}						
			
			xa = xa_64;
			ya = ya_64;
			za = za_64;
			g  = 1.0;		
			break;
		}	
		default:
			ExceptionT::GeneralFail(caller);
	}
	
	/* evaluate shape functions at integration points */
	dArrayT coords(3), na;
	for (int i = 0; i < numint; i++)	
	{
		/* ip coordinates */
		coords[0] = g*xa[i];
		coords[1] = g*ya[i];
		coords[2] = g*za[i];
		
		/* shape function values */
		Na.RowAlias(i, na);
	
		/* evaluate (static binding) */
		HexahedronT::EvaluateShapeFunctions(coords, na, Na_x[i]);
	}
}

/* compute gradients of the "bubble" modes */
void HexahedronT::BubbleModeGradients(ArrayT<dArray2DT>& Na_x) const
{
	const char caller[] = "HexahedronT::BubbleModeGradients";

	/* limit integration rules */
	int nip = Na_x.Length();
	if (nip != 8 && nip != 9)
		ExceptionT::GeneralFail(caller, "only 8 or 9 point rules defined: %d", nip);

	/* integration rules */
	double	ra[9] = {-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0, 0.0};
	double  sa[9] = {-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0, 0.0};
	double  ta[9] = {-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0, 0.0};
	double g = (nip == 8) ? 1.0/sqrt3 : sqrt(3.0/5.0);

	/* set gradients */
	for (int i = 0; i < nip; i++)
	{
		dArray2DT& na_x = Na_x[i];
		if (na_x.MajorDim() != 3 || na_x.MinorDim() != 3)
			ExceptionT::GeneralFail(caller, "gradients array must be 3x3: %dx%d", 
				na_x.MajorDim(), na_x.MinorDim());

		/* integration point coordinates */
		double r = g*ra[i];
		double s = g*sa[i];
		double t = g*ta[i];
	
		/* set gradients */
		double* nax = na_x(0);
		double* nay = na_x(1);
		double* naz = na_x(2);
		nax[0] = r;
		nax[1] = 0.0;
		nax[2] = 0.0;
		nay[0] = 0.0;
		nay[1] = s;
		nay[2] = 0.0;
		naz[0] = 0.0;
		naz[1] = 0.0;
		naz[2] = t;
	}
}

/* return the local node numbers for each facet of the element
* numbered to produce at outward normal in the order: vertex
* nodes, mid-edge nodes, mid-face nodes */
void HexahedronT::NodesOnFacet(int facet, iArrayT& facetnodes) const
{
	const char caller[] = "HexahedronT::NodesOnFacet";
	if (fNumNodes != 8 && fNumNodes != 20 && fNumNodes != 27)
		ExceptionT::GeneralFail(caller, "only implemented 8 and 20 and 27 element nodes: %d", fNumNodes);

#if __option(extended_errorcheck)
	if (facet < 0 || facet > 5) ExceptionT::OutOfRange(caller);
#endif

	/* nodes-facet data */
	int dat8[] = {0,3,2,1,
		      4,5,6,7,
		      0,1,5,4,
		      1,2,6,5,
		      2,3,7,6,
		      3,0,4,7};

	int dat20[] = {0,3,2,1,11,10, 9, 8,
		       4,5,6,7,12,13,14,15,
		       0,1,5,4, 8,17,12,16,
		       1,2,6,5, 9,18,13,17,
		       2,3,7,6,10,19,14,18,
		       3,0,4,7,11,16,15,19};

	int dat27[] = {0,3,2,1,11,10, 9, 8,20,
		       4,5,6,7,12,13,14,15,21,
		       0,1,5,4, 8,17,12,16,22,
		       1,2,6,5, 9,18,13,17,25,
		       2,3,7,6,10,19,14,18,23,
		       3,0,4,7,11,16,15,19,24};

	/* collect facet data */		
	iArrayT tmp;
	if (fNumNodes == 8)
		tmp.Set(4, dat8 + facet*4);
	else if (fNumNodes == 20)
		tmp.Set(8, dat20 + facet*8);
	else
	    tmp.Set(9, dat27 + facet*9); // I'm not sure about this line/*Davoud*/
	/* (allocate and) copy in */
	facetnodes = tmp;
}

/* return the local node numbers for each edge of element */
void HexahedronT::NodesOnEdges(iArray2DT& nodes_on_edges) const
{
	/* edges in 8-node hex */
	int dat8[12*2] = {
		0,1,
		1,2,
		2,3,
		3,0,
		4,5,
		5,6,
		6,7,
		7,4,
		0,4,
		1,5,
		2,6,
		3,7
	};

	/* edges in 20-node hex */
	int dat20[12*3] = {
		0,8,1,
		1,9,2,
		2,10,3,
		3,11,0,
		4,12,5,
		5,13,6,
		6,14,7,
		7,15,4,
		0,16,4,
		1,17,5,
		2,18,6,
		3,19,7
	};

	/* edges in 27-node hex */
	int dat27[12*3] = {
		0,8,1,
		1,9,2,
		2,10,3,
		3,11,0,
		4,12,5,
		5,13,6,
		6,14,7,
		7,15,4,
		0,16,4,
		1,17,5,
		2,18,6,
		3,19,7
	};

	iArray2DT tmp;
	if (fNumNodes == 8)
		tmp.Alias(12, 2, dat8);
	else if (fNumNodes == 20)
		tmp.Alias(12, 3, dat20);
	else if (fNumNodes == 27)
		tmp.Alias(12, 3, dat27);
	else
		ExceptionT::OutOfRange("HexahedronT::NodesOnEdges");

	/* copy in */
	nodes_on_edges = tmp;
}

void HexahedronT::NumNodesOnFacets(iArrayT& num_nodes) const
{
//TEMP
	if (fNumNodes != 8 && fNumNodes != 20 && fNumNodes != 27)
		ExceptionT::GeneralFail("HexahedronT::NumNodesOnFacets", "only implemented 8 and 20 and 27 element nodes: %d", fNumNodes);

	num_nodes.Dimension(6);
	if (fNumNodes == 8)
		num_nodes = 4;
	else if(fNumNodes == 20) 
		num_nodes = 8;
	else
	    num_nodes = 9; //I'm not sure about this line/*Davoud*/
}

/* returns the nodes on each facet needed to determine neighbors
* across facets */
void HexahedronT::NeighborNodeMap(iArray2DT& facetnodes) const
{
	/* nodes-facet data */
	int dat8[] = {0,3,2,1,
		          4,5,6,7,
		          0,1,5,4,
		          1,2,6,5,
		          2,3,7,6,
		          3,0,4,7};
	iArray2DT temp(6, 4, dat8);

	facetnodes = temp;
}

/* return geometry and number of nodes on each facet */
void HexahedronT::FacetGeometry(ArrayT<CodeT>& facet_geom, iArrayT& facet_nodes) const
{
	if (fNumNodes != 8 && fNumNodes != 20 && fNumNodes != 27)
		ExceptionT::GeneralFail("HexahedronT::FacetGeometry", "only implemented for 8 and 20 and 27 nodes: %d", fNumNodes);

	facet_geom.Dimension(fNumFacets);
	facet_geom = kQuadrilateral;
	
	facet_nodes.Dimension(fNumFacets);

	if (fNumNodes == 8)	
		facet_nodes = 4;
	else if(fNumNodes == 20) 
		facet_nodes = 8;
	else
	    facet_nodes = 9; //I'm not sure about this line/*Davoud*/
}

/* set the values of the nodal extrapolation matrix */
void HexahedronT::SetExtrapolation(dMatrixT& extrap) const
{
	const char caller[] = "HexahedronT::SetExtrapolation";

	/* dimensions */
	int numnodes = extrap.Rows();
	int numint   = extrap.Cols();

	/* dimension checks */
	if (numnodes < 8 || numnodes > 27) ExceptionT::GeneralFail(caller);

	/* initialize */
	extrap = 0.0;
	
	switch (numint)
	{
		case 1:	
			
		extrap = 1.0;
		break;

		case 8:
		{
			/* smoothin matrix data: [max nen] x [nip = 8] */
			double data_160[27*8] = {
  0.125*(1. + 3.*sqrt3)   , 0.125*(1. + sqrt3)      , 0.125*(1. - 1.*sqrt3)   , 0.125*(1. + sqrt3)
, 0.125*(1. + sqrt3)      , 0.125*(1. - 1.*sqrt3)   , 0.125*(1. - 3.*sqrt3)   , 0.125*(1. - 1.*sqrt3)
, 0.125*(1. + 2.*sqrt3)   , 0.125                   , 0.125                   , 0.125*(1. + 2.*sqrt3)
, 0.125                   , 0.125*(1. - 2.*sqrt3)   , 0.125*(1. - 2.*sqrt3)   , 0.125
, 0.125*(1. + 2.*sqrt3)   , 0.125                   , 0.125*(1. - 2.*sqrt3)   , 0.125
,0.0,0.0,0.0,0.0,0.0,0.0,0.0
, 0.125*(1. + sqrt3)      , 0.125*(1. + 3.*sqrt3)   , 0.125*(1. + sqrt3)      , 0.125*(1. - 1.*sqrt3)
, 0.125*(1. - 1.*sqrt3)   , 0.125*(1. + sqrt3)      , 0.125*(1. - 1.*sqrt3)   , 0.125*(1. - 3.*sqrt3)
, 0.125*(1. + 2.*sqrt3)   , 0.125*(1. + 2.*sqrt3)   , 0.125                   , 0.125
, 0.125                   , 0.125                   , 0.125*(1. - 2.*sqrt3)   , 0.125*(1. - 2.*sqrt3)
, 0.125                   , 0.125*(1. + 2.*sqrt3)   , 0.125                   , 0.125*(1. - 2.*sqrt3)
,0.0,0.0,0.0,0.0,0.0,0.0,0.0
, 0.125*(1. - 1.*sqrt3)   , 0.125*(1. + sqrt3)      , 0.125*(1. + 3.*sqrt3)   , 0.125*(1. + sqrt3)
, 0.125*(1. - 3.*sqrt3)   , 0.125*(1. - 1.*sqrt3)   , 0.125*(1. + sqrt3)      , 0.125*(1. - 1.*sqrt3)
, 0.125                   , 0.125*(1. + 2.*sqrt3)   , 0.125*(1. + 2.*sqrt3)   , 0.125
, 0.125*(1. - 2.*sqrt3)   , 0.125                   , 0.125                   , 0.125*(1. - 2.*sqrt3)
, 0.125*(1. - 2.*sqrt3)   , 0.125                   , 0.125*(1. + 2.*sqrt3)   , 0.125
,0.0,0.0,0.0,0.0,0.0,0.0,0.0
, 0.125*(1. + sqrt3)      , 0.125*(1. - 1.*sqrt3)   , 0.125*(1. + sqrt3)      , 0.125*(1. + 3.*sqrt3)
, 0.125*(1. - 1.*sqrt3)   , 0.125*(1. - 3.*sqrt3)   , 0.125*(1. - 1.*sqrt3)   , 0.125*(1. + sqrt3)
, 0.125                   , 0.125                   , 0.125*(1. + 2.*sqrt3)   , 0.125*(1. + 2.*sqrt3)
, 0.125*(1. - 2.*sqrt3)   , 0.125*(1. - 2.*sqrt3)   , 0.125                   , 0.125
, 0.125                   , 0.125*(1. - 2.*sqrt3)   , 0.125                   , 0.125*(1. + 2.*sqrt3)
,0.0,0.0,0.0,0.0,0.0,0.0,0.0
, 0.125*(1. + sqrt3)      , 0.125*(1. - 1.*sqrt3)   , 0.125*(1. - 3.*sqrt3)   , 0.125*(1. - 1.*sqrt3)
, 0.125*(1. + 3.*sqrt3)   , 0.125*(1. + sqrt3)      , 0.125*(1. - 1.*sqrt3)   , 0.125*(1. + sqrt3)
, 0.125                   , 0.125*(1. - 2.*sqrt3)   , 0.125*(1. - 2.*sqrt3)   , 0.125
, 0.125*(1. + 2.*sqrt3)   , 0.125                   , 0.125                   , 0.125*(1. + 2.*sqrt3)
, 0.125*(1. + 2.*sqrt3)   , 0.125                   , 0.125*(1. - 2.*sqrt3)   , 0.125
,0.0,0.0,0.0,0.0,0.0,0.0,0.0
, 0.125*(1. - 1.*sqrt3)   , 0.125*(1. + sqrt3)      , 0.125*(1. - 1.*sqrt3)   , 0.125*(1. - 3.*sqrt3)
, 0.125*(1. + sqrt3)      , 0.125*(1. + 3.*sqrt3)   , 0.125*(1. + sqrt3)      , 0.125*(1. - 1.*sqrt3)
, 0.125                   , 0.125                   , 0.125*(1. - 2.*sqrt3)   , 0.125*(1. - 2.*sqrt3)
, 0.125*(1. + 2.*sqrt3)   , 0.125*(1. + 2.*sqrt3)   , 0.125                   , 0.125
, 0.125                   , 0.125*(1. + 2.*sqrt3)   , 0.125                   , 0.125*(1. - 2.*sqrt3)
,0.0,0.0,0.0,0.0,0.0,0.0,0.0
, 0.125*(1. - 3.*sqrt3)   , 0.125*(1. - 1.*sqrt3)   , 0.125*(1. + sqrt3)      , 0.125*(1. - 1.*sqrt3)
, 0.125*(1. - 1.*sqrt3)   , 0.125*(1. + sqrt3)      , 0.125*(1. + 3.*sqrt3)   , 0.125*(1. + sqrt3)
, 0.125*(1. - 2.*sqrt3)   , 0.125                   , 0.125                   , 0.125*(1. - 2.*sqrt3)
, 0.125                   , 0.125*(1. + 2.*sqrt3)   , 0.125*(1. + 2.*sqrt3)   , 0.125
, 0.125*(1. - 2.*sqrt3)   , 0.125                   , 0.125*(1. + 2.*sqrt3)   , 0.125
,0.0,0.0,0.0,0.0,0.0,0.0,0.0
, 0.125*(1. - 1.*sqrt3)   , 0.125*(1. - 3.*sqrt3)   , 0.125*(1. - 1.*sqrt3)   , 0.125*(1. + sqrt3)
, 0.125*(1. + sqrt3)      , 0.125*(1. - 1.*sqrt3)   , 0.125*(1. + sqrt3)      , 0.125*(1. + 3.*sqrt3)
, 0.125*(1. - 2.*sqrt3)   , 0.125*(1. - 2.*sqrt3)   , 0.125                   , 0.125
, 0.125                   , 0.125                   , 0.125*(1. + 2.*sqrt3)   , 0.125*(1. + 2.*sqrt3)
, 0.125                   , 0.125*(1. - 2.*sqrt3)   , 0.125                   , 0.125*(1. + 2.*sqrt3)
,0.0,0.0,0.0,0.0,0.0,0.0,0.0
	};

			/* copy block out */	
			dMatrixT smooth_160(27, 8, data_160);
			smooth_160.CopyBlock(0, 0, extrap);
			break;
		}
		case 9:
		{
			/* smoothin matrix data: [max nen] x [nip = 9] */
			double data_180[27*9] = {
5.799038105676658, 2.566987298107781, 3.433012701892219, 2.566987298107781, 
2.566987298107781, 3.433012701892219, 3.200961894323342, 3.433012701892219, 
0.808012701892219, -0.375, -0.375, 0.808012701892219, 
-0.375, -0.05801270189221924, -0.05801270189221924, -0.375, 
0.808012701892219, -0.375, -0.05801270189221924, -0.375, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
2.566987298107781, 5.799038105676658, 2.566987298107781, 3.433012701892219, 
3.433012701892219, 2.566987298107781, 3.433012701892219, 3.200961894323342, 
0.808012701892219, 0.808012701892219, -0.375, -0.375, 
-0.375, -0.375, -0.05801270189221924, -0.05801270189221924, 
-0.375, 0.808012701892219, -0.375, -0.05801270189221924, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
3.433012701892219, 2.566987298107781, 5.799038105676658, 2.566987298107781, 
3.200961894323342, 3.433012701892219, 2.566987298107781, 3.433012701892219, 
-0.375, 0.808012701892219, 0.808012701892219, -0.375, 
-0.05801270189221924, -0.375, -0.375, -0.05801270189221924, 
-0.05801270189221924, -0.375, 0.808012701892219, -0.375,
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
2.566987298107781, 3.433012701892219, 2.566987298107781, 5.799038105676658, 
3.433012701892219, 3.200961894323342, 3.433012701892219, 2.566987298107781, 
-0.375, -0.375, 0.808012701892219, 0.808012701892219, 
-0.05801270189221924, -0.05801270189221924, -0.375, -0.375, 
-0.375, -0.05801270189221924, -0.375, 0.808012701892219,
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
2.566987298107781, 3.433012701892219, 3.200961894323342, 3.433012701892219, 
5.799038105676658, 2.566987298107781, 3.433012701892219, 2.566987298107781, 
-0.375, -0.05801270189221924, -0.05801270189221924, -0.375, 
0.808012701892219, -0.375, -0.375, 0.808012701892219, 
0.808012701892219, -0.375, -0.05801270189221924, -0.375, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
3.433012701892219, 2.566987298107781, 3.433012701892219, 3.200961894323342, 
2.566987298107781, 5.799038105676658, 2.566987298107781, 3.433012701892219, 
-0.375, -0.375, -0.05801270189221924, -0.05801270189221924, 
0.808012701892219, 0.808012701892219, -0.375, -0.375, 
-0.375, 0.808012701892219, -0.375, -0.05801270189221924, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
3.200961894323342, 3.433012701892219, 2.566987298107781, 3.433012701892219, 
3.433012701892219, 2.566987298107781, 5.799038105676658, 2.566987298107781, 
-0.05801270189221924, -0.375, -0.375, -0.05801270189221924, 
-0.375, 0.808012701892219, 0.808012701892219, -0.375, 
-0.05801270189221924, -0.375, 0.808012701892219, -0.375, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
3.433012701892219, 3.200961894323342, 3.433012701892219, 2.566987298107781, 
2.566987298107781, 3.433012701892219, 2.566987298107781, 5.799038105676658, 
-0.05801270189221924, -0.05801270189221924, -0.375, -0.375, 
-0.375, -0.375, 0.808012701892219, 0.808012701892219, 
-0.375, -0.05801270189221924, -0.375, 0.808012701892219, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
-26., -26., -26., -26., 
-26., -26., -26., -26., 
1., 1., 1., 1., 
1., 1., 1., 1., 
1., 1., 1., 1.,
0.0,0.0,0.0,0.0,0.0,0.0,0.0
};

			/* copy block out */	
			dMatrixT smooth_180(27, 9, data_180);
			smooth_180.CopyBlock(0, 0, extrap);
			break;
		}
		case 27:
		{
			/* smoothin matrix data: [max nen] x [nip = 27] */
			double data_540[27*27] = {
0.00925925925925926*(175. + 45.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 
0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(175. - 45.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 
0, 0, 0, 0, 
0, 0, 0, 0, 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(-80. - 20.*sqrt15), 0.00925925925925926*(-80. - 20.*sqrt15), -0.1851851851851852, -0.1851851851851852, 
-0.1851851851851852, -0.1851851851851852, 0.00925925925925926*(-80. + 20.*sqrt15), 0.00925925925925926*(-80. + 20.*sqrt15), 
0.05555555555555556*(20. + 5.*sqrt15), 0, 0.2777777777777778, 0, 
0.2777777777777778, 0, 0.05555555555555556*(20. - 5.*sqrt15), 0, 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(175. + 45.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 
0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(175. - 45.*sqrt15), 
0, 0, 0, 0, 
0, 0, 0, 0, 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(-80. - 20.*sqrt15), -0.1851851851851852, -0.1851851851851852, 0.00925925925925926*(-80. - 20.*sqrt15), 
-0.1851851851851852, 0.00925925925925926*(-80. + 20.*sqrt15), 0.00925925925925926*(-80. + 20.*sqrt15), -0.1851851851851852, 
0, 0.2777777777777778, 0, 0.05555555555555556*(20. + 5.*sqrt15), 
0, 0.05555555555555556*(20. - 5.*sqrt15), 0, 0.2777777777777778, 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 
0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 
0.05555555555555556*(-10. - 2.*sqrt15), 0.05555555555555556*(-10. - 2.*sqrt15), 0.05555555555555556*(-10. - 2.*sqrt15), 0.05555555555555556*(-10. - 2.*sqrt15), 
0.05555555555555556*(-10. + 2.*sqrt15), 0.05555555555555556*(-10. + 2.*sqrt15), 0.05555555555555556*(-10. + 2.*sqrt15), 0.05555555555555556*(-10. + 2.*sqrt15), 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
-0.1851851851851852, 0.00925925925925926*(-80. - 20.*sqrt15), 0.00925925925925926*(-80. - 20.*sqrt15), -0.1851851851851852, 
0.00925925925925926*(-80. + 20.*sqrt15), -0.1851851851851852, -0.1851851851851852, 0.00925925925925926*(-80. + 20.*sqrt15), 
0, 0.05555555555555556*(20. + 5.*sqrt15), 0, 0.2777777777777778, 
0, 0.2777777777777778, 0, 0.05555555555555556*(20. - 5.*sqrt15), 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(175. + 45.*sqrt15), 
0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(175. - 45.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 
0, 0, 0, 0, 
0, 0, 0, 0, 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
-0.1851851851851852, -0.1851851851851852, 0.00925925925925926*(-80. - 20.*sqrt15), 0.00925925925925926*(-80. - 20.*sqrt15), 
0.00925925925925926*(-80. + 20.*sqrt15), 0.00925925925925926*(-80. + 20.*sqrt15), -0.1851851851851852, -0.1851851851851852, 
0.2777777777777778, 0, 0.05555555555555556*(20. + 5.*sqrt15), 0, 
0.05555555555555556*(20. - 5.*sqrt15), 0, 0.2777777777777778, 0, 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(175. + 45.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 
0.00925925925925926*(175. - 45.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 
0, 0, 0, 0, 
0, 0, 0, 0, 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(-80. - 20.*sqrt15), -0.1851851851851852, 0.00925925925925926*(-80. + 20.*sqrt15), -0.1851851851851852, 
0.00925925925925926*(-80. - 20.*sqrt15), -0.1851851851851852, 0.00925925925925926*(-80. + 20.*sqrt15), -0.1851851851851852, 
0, 0, 0, 0, 
0, 0, 0, 0, 
0.05555555555555556*(20. + 5.*sqrt15), 0.2777777777777778, 0.05555555555555556*(20. - 5.*sqrt15), 0.2777777777777778, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 
0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 
0.05555555555555556*(-10. - 2.*sqrt15), 0, 0.05555555555555556*(-10. + 2.*sqrt15), 0, 
0.05555555555555556*(-10. - 2.*sqrt15), 0, 0.05555555555555556*(-10. + 2.*sqrt15), 0, 
0.05555555555555556*(-10. - 2.*sqrt15), 0.05555555555555556*(-10. - 2.*sqrt15), 0.05555555555555556*(-10. + 2.*sqrt15), 0.05555555555555556*(-10. + 2.*sqrt15), 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
-0.1851851851851852, 0.00925925925925926*(-80. - 20.*sqrt15), -0.1851851851851852, 0.00925925925925926*(-80. + 20.*sqrt15), 
-0.1851851851851852, 0.00925925925925926*(-80. - 20.*sqrt15), -0.1851851851851852, 0.00925925925925926*(-80. + 20.*sqrt15), 
0, 0, 0, 0, 
0, 0, 0, 0, 
0.2777777777777778, 0.05555555555555556*(20. + 5.*sqrt15), 0.2777777777777778, 0.05555555555555556*(20. - 5.*sqrt15), 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 
0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 
0, 0.05555555555555556*(-10. + 2.*sqrt15), 0, 0.05555555555555556*(-10. - 2.*sqrt15), 
0, 0.05555555555555556*(-10. + 2.*sqrt15), 0, 0.05555555555555556*(-10. - 2.*sqrt15), 
0.05555555555555556*(-10. - 2.*sqrt15), 0.05555555555555556*(-10. + 2.*sqrt15), 0.05555555555555556*(-10. + 2.*sqrt15), 0.05555555555555556*(-10. - 2.*sqrt15), 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
-0.2962962962962963, -0.2962962962962963, -0.2962962962962963, -0.2962962962962963, 
-0.2962962962962963, -0.2962962962962963, -0.2962962962962963, -0.2962962962962963, 
0.4444444444444445, 0.4444444444444445, 0.4444444444444445, 0.4444444444444445, 
0.4444444444444445, 0.4444444444444445, 0.4444444444444445, 0.4444444444444445, 
0.4444444444444445, 0.4444444444444445, 0.4444444444444445, 0.4444444444444445, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 
0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 
0, 0.05555555555555556*(-10. - 2.*sqrt15), 0, 0.05555555555555556*(-10. + 2.*sqrt15), 
0, 0.05555555555555556*(-10. - 2.*sqrt15), 0, 0.05555555555555556*(-10. + 2.*sqrt15), 
0.05555555555555556*(-10. + 2.*sqrt15), 0.05555555555555556*(-10. - 2.*sqrt15), 0.05555555555555556*(-10. - 2.*sqrt15), 0.05555555555555556*(-10. + 2.*sqrt15), 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
-0.1851851851851852, 0.00925925925925926*(-80. + 20.*sqrt15), -0.1851851851851852, 0.00925925925925926*(-80. - 20.*sqrt15), 
-0.1851851851851852, 0.00925925925925926*(-80. + 20.*sqrt15), -0.1851851851851852, 0.00925925925925926*(-80. - 20.*sqrt15), 
0, 0, 0, 0, 
0, 0, 0, 0, 
0.2777777777777778, 0.05555555555555556*(20. - 5.*sqrt15), 0.2777777777777778, 0.05555555555555556*(20. + 5.*sqrt15), 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 
0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 
0.05555555555555556*(-10. + 2.*sqrt15), 0, 0.05555555555555556*(-10. - 2.*sqrt15), 0, 
0.05555555555555556*(-10. + 2.*sqrt15), 0, 0.05555555555555556*(-10. - 2.*sqrt15), 0, 
0.05555555555555556*(-10. + 2.*sqrt15), 0.05555555555555556*(-10. + 2.*sqrt15), 0.05555555555555556*(-10. - 2.*sqrt15), 0.05555555555555556*(-10. - 2.*sqrt15), 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(-80. + 20.*sqrt15), -0.1851851851851852, 0.00925925925925926*(-80. - 20.*sqrt15), -0.1851851851851852, 
0.00925925925925926*(-80. + 20.*sqrt15), -0.1851851851851852, 0.00925925925925926*(-80. - 20.*sqrt15), -0.1851851851851852, 
0, 0, 0, 0, 
0, 0, 0, 0, 
0.05555555555555556*(20. - 5.*sqrt15), 0.2777777777777778, 0.05555555555555556*(20. + 5.*sqrt15), 0.2777777777777778, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(175. - 45.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 
0.00925925925925926*(175. + 45.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 
0, 0, 0, 0, 
0, 0, 0, 0, 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
-0.1851851851851852, -0.1851851851851852, 0.00925925925925926*(-80. + 20.*sqrt15), 0.00925925925925926*(-80. + 20.*sqrt15), 
0.00925925925925926*(-80. - 20.*sqrt15), 0.00925925925925926*(-80. - 20.*sqrt15), -0.1851851851851852, -0.1851851851851852, 
0.2777777777777778, 0, 0.05555555555555556*(20. - 5.*sqrt15), 0, 
0.05555555555555556*(20. + 5.*sqrt15), 0, 0.2777777777777778, 0, 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(175. - 45.*sqrt15), 
0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(175. + 45.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 
0, 0, 0, 0, 
0, 0, 0, 0, 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
-0.1851851851851852, 0.00925925925925926*(-80. + 20.*sqrt15), 0.00925925925925926*(-80. + 20.*sqrt15), -0.1851851851851852, 
0.00925925925925926*(-80. - 20.*sqrt15), -0.1851851851851852, -0.1851851851851852, 0.00925925925925926*(-80. - 20.*sqrt15), 
0, 0.05555555555555556*(20. - 5.*sqrt15), 0, 0.2777777777777778, 
0, 0.2777777777777778, 0, 0.05555555555555556*(20. + 5.*sqrt15), 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 0.00925925925925926*(40. - 8.*sqrt15), 
0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 0.00925925925925926*(40. + 8.*sqrt15), 
0.05555555555555556*(-10. + 2.*sqrt15), 0.05555555555555556*(-10. + 2.*sqrt15), 0.05555555555555556*(-10. + 2.*sqrt15), 0.05555555555555556*(-10. + 2.*sqrt15), 
0.05555555555555556*(-10. - 2.*sqrt15), 0.05555555555555556*(-10. - 2.*sqrt15), 0.05555555555555556*(-10. - 2.*sqrt15), 0.05555555555555556*(-10. - 2.*sqrt15), 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(-80. + 20.*sqrt15), -0.1851851851851852, -0.1851851851851852, 0.00925925925925926*(-80. + 20.*sqrt15), 
-0.1851851851851852, 0.00925925925925926*(-80. - 20.*sqrt15), 0.00925925925925926*(-80. - 20.*sqrt15), -0.1851851851851852, 
0, 0.2777777777777778, 0, 0.05555555555555556*(20. - 5.*sqrt15), 
0, 0.05555555555555556*(20. + 5.*sqrt15), 0, 0.2777777777777778, 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(175. - 45.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 
0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(175. + 45.*sqrt15), 
0, 0, 0, 0, 
0, 0, 0, 0, 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(-80. + 20.*sqrt15), 0.00925925925925926*(-80. + 20.*sqrt15), -0.1851851851851852, -0.1851851851851852, 
-0.1851851851851852, -0.1851851851851852, 0.00925925925925926*(-80. - 20.*sqrt15), 0.00925925925925926*(-80. - 20.*sqrt15), 
0.05555555555555556*(20. - 5.*sqrt15), 0, 0.2777777777777778, 0, 
0.2777777777777778, 0, 0.05555555555555556*(20. + 5.*sqrt15), 0, 
0, 0, 0, 0, 
0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.00925925925925926*(175. - 45.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(25. - 5.*sqrt15), 
0.00925925925925926*(25. - 5.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 0.00925925925925926*(175. + 45.*sqrt15), 0.00925925925925926*(25. + 5.*sqrt15), 
0, 0, 0, 0, 
0, 0, 0, 0, 
0, 0, 0, 0,
0.0,0.0,0.0,0.0,0.0,0.0,0.0	
			};

			/* copy block out */
			dMatrixT smooth_540(27, 27, data_540);
			smooth_540.CopyBlock(0, 0, extrap);
			break;
		}
		default:		
			ExceptionT::GeneralFail(caller, "no nodal extrapolation with Gauss rule: %d", numint);
	}
}

/* integration point gradient matrix */
void HexahedronT::IPGradientTransform(int ip, dMatrixT& transform) const
{
	const char caller[] = "HexahedronT::IPGradientTransform";

	/* dimensions */
	int nsd = transform.Rows();
	int nip = transform.Cols();
	if (nsd != 3) ExceptionT::SizeMismatch(caller);
	
	if (nip == 1)
		transform = 0.0;
	else if (nip == 8) {
		double a = sqrt(3.0)/2.0;
		double m0[3*8] = {-a, -a, -a, a, 0, 0, 0, 0, 0, 0, a, 0, 0, 0, a, 0, 0, 0, 0, 0, 0, 0, 0, 0};
		double m1[3*8] = {-a, 0, 0, a, -a, -a, 0, a, 0, 0, 0, 0, 0, 0, 0, 0, 0, a, 0, 0, 0, 0, 0, 0};
		double m2[3*8] = {0, 0, 0, 0, -a, 0, a, a, -a, -a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, a, 0, 0, 0};
		double m3[3*8] = {0, -a, 0, 0, 0, 0, a, 0, 0, -a, a, -a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, a};
		double m4[3*8] = {0, 0, -a, 0, 0, 0, 0, 0, 0, 0, 0, 0, -a, -a, a, a, 0, 0, 0, 0, 0, 0, a, 0};
		double m5[3*8] = {0, 0, 0, 0, 0, -a, 0, 0, 0, 0, 0, 0, -a, 0, 0, a, -a, a, 0, a, 0, 0, 0, 0};
		double m6[3*8] = {0, 0, 0, 0, 0, 0, 0, 0, -a, 0, 0, 0, 0, 0, 0, 0, -a, 0, a, a, a, -a, 0, 0};
		double m7[3*8] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -a, 0, -a, 0, 0, 0, 0, a, 0, 0, -a, a, a};
		double* m[8] = {m0, m1, m2, m3, m4, m5, m6, m7};
		ArrayT<double*> m_array(8, m);
		dMatrixT trans(3, 8, m_array[ip]);
		transform = trans;		
	}
	else
		ExceptionT::GeneralFail(caller, "unsupported number of integration points %d", nip);	
}

/* return true if the given point is within the domain */
bool HexahedronT::PointInDomain(const LocalArrayT& coords, const dArrayT& point) const
{
#if __option(extended_errorcheck)
		if (coords.NumberOfNodes() != 8) 
			ExceptionT::GeneralFail("HexahedronT::PointInDomain", "expecting 8 element nodes: %d", coords.NumberOfNodes());
#endif

	/* nodes-facet data - ordered for outward normals */
	int dat8[] = {
		0, 3, 2, 1,
		4, 5, 6, 7,
		0, 1, 5, 4,
		1, 2, 6, 5,
		2, 3, 7, 6,
		3, 0, 4, 7};

	/* method: check all faces and see of point lies inside, breaking
	*          each face into 2 triangular facets */
	bool in_domain = true;
	int* facet_nodes = dat8;
	for (int i = 0; in_domain && i < 6; i++)
	{
		/* facet 1 */
		double ab_0 = coords(facet_nodes[1], 0) - coords(facet_nodes[0], 0);
		double ab_1 = coords(facet_nodes[1], 1) - coords(facet_nodes[0], 1);
		double ab_2 = coords(facet_nodes[1], 2) - coords(facet_nodes[0], 2);
		double ab_max = (fabs(ab_0) > fabs(ab_1)) ? fabs(ab_0) : fabs(ab_1);
		ab_max = (fabs(ab_2) > ab_max) ? fabs(ab_2) : ab_max;

		double ac_0 = coords(facet_nodes[3], 0) - coords(facet_nodes[0], 0);
		double ac_1 = coords(facet_nodes[3], 1) - coords(facet_nodes[0], 1);
		double ac_2 = coords(facet_nodes[3], 2) - coords(facet_nodes[0], 2);
		double ac_max = (fabs(ac_0) > fabs(ac_1)) ? fabs(ac_0) : fabs(ac_1);
		ac_max = (fabs(ac_2) > ac_max) ? fabs(ac_2) : ac_max;

		double L_ref = (ab_max > ac_max) ? ab_max : ac_max;

		double ap_0 = point[0] - coords(facet_nodes[0], 0);
		double ap_1 = point[1] - coords(facet_nodes[0], 1);
		double ap_2 = point[2] - coords(facet_nodes[0], 2);
			
		/* vector triple product */
		double ac_ab_0 = ac_1*ab_2 - ac_2*ab_1;
		double ac_ab_1 = ac_2*ab_0 - ac_0*ab_2;
		double ac_ab_2 = ac_0*ab_1 - ac_1*ab_0;			
		double triple_product = ac_ab_0*ap_0 + ac_ab_1*ap_1 + ac_ab_2*ap_2;
		in_domain = (triple_product/(L_ref*L_ref*L_ref)) > -kSmall;

		/* facet 2 */
		if (in_domain) {

			ab_0 = coords(facet_nodes[3], 0) - coords(facet_nodes[2], 0);
			ab_1 = coords(facet_nodes[3], 1) - coords(facet_nodes[2], 1);
			ab_2 = coords(facet_nodes[3], 2) - coords(facet_nodes[2], 2);
			ab_max = (fabs(ab_0) > fabs(ab_1)) ? fabs(ab_0) : fabs(ab_1);
			ab_max = (fabs(ab_2) > ab_max) ? fabs(ab_2) : ab_max;
		
			ac_0 = coords(facet_nodes[1], 0) - coords(facet_nodes[2], 0);
			ac_1 = coords(facet_nodes[1], 1) - coords(facet_nodes[2], 1);
			ac_2 = coords(facet_nodes[1], 2) - coords(facet_nodes[2], 2);
			ac_max = (fabs(ac_0) > fabs(ac_1)) ? fabs(ac_0) : fabs(ac_1);
			ac_max = (fabs(ac_2) > ac_max) ? fabs(ac_2) : ac_max;
			
			L_ref = (ab_max > ac_max) ? ab_max : ac_max;
			
			ap_0 = point[0] - coords(facet_nodes[2], 0);
			ap_1 = point[1] - coords(facet_nodes[2], 1);
			ap_2 = point[2] - coords(facet_nodes[2], 2);

			/* vector triple product */
			ac_ab_0 = ac_1*ab_2 - ac_2*ab_1;
			ac_ab_1 = ac_2*ab_0 - ac_0*ab_2;
			ac_ab_2 = ac_0*ab_1 - ac_1*ab_0;			
			triple_product = ac_ab_0*ap_0 + ac_ab_1*ap_1 + ac_ab_2*ap_2;
			in_domain = (triple_product/(L_ref*L_ref*L_ref)) > -kSmall;
		}
		
		facet_nodes += 4;
	}
	
	return in_domain;
}

/* subdomain geometry */
GeometryT::CodeT HexahedronT::NodalSubDomainGeometry(void) const
{
	/* limited support */
	if (fNumNodes != 8)
		ExceptionT::GeneralFail("HexahedronT::NodalSubDomainGeometry",
			"unsupported number of nodes %d", fNumNodes);

	return GeometryT::kHexahedron;
}

/* number of nodes defining the nodal subdomain */
int HexahedronT::NodalSubDomainNumPoints(void) const
{
	/* limited support */
	if (fNumNodes != 8)
		ExceptionT::GeneralFail("HexahedronT::NodalSubDomainGeometry",
			"unsupported number of nodes %d", fNumNodes);

	return 8;
}
	
/* compute the coordinates of the points defining the nodal subdomain */
void HexahedronT::NodalSubDomainCoordinates(const LocalArrayT& coords, int node,
	LocalArrayT& subdomain_coords) const
{
	const char caller[] = "HexahedronT::NodalSubDomainCoordinates";

#if __option(extended_errorcheck)
	/* limited support */
	if (fNumNodes != 8)
		ExceptionT::GeneralFail(caller, "unsupported number of nodes %d", fNumNodes);

	/* checks */
	if (coords.NumberOfNodes() != fNumNodes || 
		coords.MinorDim() != 3 ||
		node < 0 || node >= fNumNodes ||
		subdomain_coords.MinorDim() != 3 ||
		subdomain_coords.NumberOfNodes() != HexahedronT::NodalSubDomainNumPoints())
		ExceptionT::SizeMismatch(caller);
#endif

	int me[8] = {0,1,2,3,4,5,6,7};
	int x[8] = {1,0,3,2,5,4,7,6};
	int y[8] = {3,2,1,0,7,6,5,4};
	int z[8] = {4,5,6,7,0,1,2,3};
	int xy[8] = {2,3,0,1,6,7,4,5};
	int zx[8] = {5,4,7,6,1,0,3,2};
	int zy[8] = {7,6,5,4,3,2,1,0};
	int xyz[8] = {6,7,4,5,2,3,0,1};
	
	int* seq[8][8] = {
		{me, x, xy, y, z, zx, xyz, zy},
		{me, y, xy, x, z, zy, xyz, zx},
		{me, x, xy, y, z, zx, xyz, zy},
		{me, y, xy, x, z, zy, xyz, zx},
		{me, y, xy, x, z, zy, xyz, zx},
		{me, x, xy, y, z, zx, xyz, zy},
		{me, y, xy, x, z, zy, xyz, zx},
		{me, x, xy, y, z, zx, xyz, zy}
	};

	const double* px = coords(0);
	const double* py = coords(1);
	const double* pz = coords(2);

	int** seq_node = seq[node];
	int lnd;

	/* self */
	lnd = seq_node[0][node];
	subdomain_coords(lnd,0) = px[node];
	subdomain_coords(lnd,1) = py[node];
	subdomain_coords(lnd,2) = pz[node];

	/* edge-center */
	lnd = seq_node[1][node];
	subdomain_coords(lnd,0) = 0.5*(px[node] + px[lnd]);
	subdomain_coords(lnd,1) = 0.5*(py[node] + py[lnd]);
	subdomain_coords(lnd,2) = 0.5*(pz[node] + pz[lnd]);

	/* face-center */
	lnd = seq_node[2][node];
	int n1 = seq_node[1][node];
	int n3 = seq_node[3][node];
	subdomain_coords(lnd,0) = 0.25*(px[node] + px[lnd] + px[n1] + px[n3]);
	subdomain_coords(lnd,1) = 0.25*(py[node] + py[lnd] + py[n1] + py[n3]);
	subdomain_coords(lnd,2) = 0.25*(pz[node] + pz[lnd] + pz[n1] + pz[n3]);

	/* edge-center */
	lnd = seq_node[3][node];
	subdomain_coords(lnd,0) = 0.5*(px[node] + px[lnd]);
	subdomain_coords(lnd,1) = 0.5*(py[node] + py[lnd]);
	subdomain_coords(lnd,2) = 0.5*(pz[node] + pz[lnd]);

	/* edge-center */
	lnd = seq_node[4][node];
	subdomain_coords(lnd,0) = 0.5*(px[node] + px[lnd]);
	subdomain_coords(lnd,1) = 0.5*(py[node] + py[lnd]);
	subdomain_coords(lnd,2) = 0.5*(pz[node] + pz[lnd]);

	/* face-center */
	lnd = seq_node[5][node];
	int n4 = seq_node[4][node];
	subdomain_coords(lnd,0) = 0.25*(px[node] + px[lnd] + px[n1] + px[n4]);
	subdomain_coords(lnd,1) = 0.25*(py[node] + py[lnd] + py[n1] + py[n4]);
	subdomain_coords(lnd,2) = 0.25*(pz[node] + pz[lnd] + pz[n1] + pz[n4]);

	/* body-center */
	lnd = seq_node[6][node];
	subdomain_coords(lnd,0) = 0.125*(px[0] + px[1] + px[2] + px[3] + px[4] + px[5] + px[6] + px[7]);
	subdomain_coords(lnd,1) = 0.125*(py[0] + py[1] + py[2] + py[3] + py[4] + py[5] + py[6] + py[7]);
	subdomain_coords(lnd,2) = 0.125*(pz[0] + pz[1] + pz[2] + pz[3] + pz[4] + pz[5] + pz[6] + pz[7]);

	/* face-center */
	lnd = seq_node[7][node];
	subdomain_coords(lnd,0) = 0.25*(px[node] + px[lnd] + px[n3] + px[n4]);
	subdomain_coords(lnd,1) = 0.25*(py[node] + py[lnd] + py[n3] + py[n4]);
	subdomain_coords(lnd,2) = 0.25*(pz[node] + pz[lnd] + pz[n3] + pz[n4]);
}
