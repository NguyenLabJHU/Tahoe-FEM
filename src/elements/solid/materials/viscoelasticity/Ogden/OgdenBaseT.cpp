/* $Id: OgdenBaseT.cpp,v 1.1 2006-10-30 23:32:05 thao Exp $ */
/* created: paklein (10/01/2000) */
#include "OgdenBaseT.h"

#include <iostream.h>
#include <math.h>

using namespace Tahoe;

/* constructor */
OgdenBaseT::OgdenBaseT(ifstreamT& in, const FSMatSupportT& support):
	ParameterInterfaceT("Ogden_base"),
//	FSSolidMatT(in, support),
	fSpectralDecomp(NumSD()),
	fSpectralDecomp3D(3),
	fb(NumSD()),
	fEigs(NumSD()),
	fdWdE(NumSD()),
	fddWddE(NumSD()),
	fModMat(dSymMatrixT::NumValues(NumSD())),
	fModulus(dSymMatrixT::NumValues(NumSD())),
	fStress(NumSD())
{
ExceptionT::GeneralFail("OgdenBaseT::OgdenBaseT", "out of date");
}

/* class specific initializations */
void OgdenBaseT::Initialize(void)
{
ExceptionT::GeneralFail("OgdenBaseT::Initialize", "out of date");
#if 0
	/* initial modulus */
	fEigs = 1.0;
	ddWddE(fEigs, fdWdE, fddWddE);
	double lambda = fddWddE(0,1);
	double mu = 0.5*(fddWddE(0,0) - fddWddE(0,1));

	if (NumSD() == 2 && PurePlaneStress())
		IsotropicT::Set_PurePlaneStress_mu_lambda(mu, lambda);
	else
	{
		double kappa = lambda + 2.0/3.0*mu;
		IsotropicT::Set_mu_kappa(mu, kappa);
	}
#endif
}

/* modulus */
const dMatrixT& OgdenBaseT::C_IJKL(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fModulus.SetToScaled(Fmat.Det(), PullBack(Fmat, OgdenBaseT::c_ijkl()));
	return fModulus;
} 
/**< \todo compute directly in spatial representation rather than transforming */
	
/* stresses */
const dSymMatrixT& OgdenBaseT::S_IJ(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fStress.SetToScaled(Fmat.Det(), PullBack(Fmat, OgdenBaseT::s_ij()));
	return fStress;
}
/**< \todo compute directly in spatial representation rather than transforming */

/* material description */
const dMatrixT& OgdenBaseT::c_ijkl(void)
{
        s_ij();
	/* stretch */
	Compute_b(fb);

	/* spectral decomposition */
	fSpectralDecomp.SpectralDecomp_Jacobi(fb, false);
	/* principal values */
	
       	const dArrayT& eigenstretch = fSpectralDecomp.Eigenvalues();
	ddWddE(eigenstretch, fdWdE, fddWddE);

	/* axial */
       	fModulus = fSpectralDecomp.EigsToRank4(fddWddE);

	/* shear terms */
      	const ArrayT<dArrayT>& eigenvectors = fSpectralDecomp.Eigenvectors();
	double diff_stretch, dtde;
	double small = 10e-12;

       	if (NumSD() == 2)
	{
	  const double& l0 = eigenstretch[0];
	  const double& l1 = eigenstretch[1];
	  
	  diff_stretch = l0-l1;
	  
	  if (fabs(diff_stretch) > small)
	    dtde = (fdWdE[0]*l1 - fdWdE[1]*l0)/diff_stretch;
	  else
	    dtde = 0.5*(fddWddE(0,0) - fddWddE(0,1));
	  MixedRank4_2D(eigenvectors[0], eigenvectors[1], fModMat);
	  fModulus.AddScaled(2.0*dtde, fModMat);
	}
	else
	{
	  const double& l0 = eigenstretch[0];
	  const double& l1 = eigenstretch[1];
	  const double& l2 = eigenstretch[2];
	  
	  diff_stretch =  l0 - l1;
	  if (fabs(diff_stretch) > small)
	    dtde = (fdWdE[0]*l1 - fdWdE[1]*l0)/diff_stretch;
	  else
	    dtde = 0.5*(fddWddE(0,0) - fddWddE(0,1))-fdWdE[0];
	  MixedRank4_3D(eigenvectors[0], eigenvectors[1], fModMat);
	  fModulus.AddScaled(2.0*dtde, fModMat);
	  
	  diff_stretch = l0 - l2;
	  if (fabs(diff_stretch) > small)
	    dtde = (fdWdE[0]*l2 - fdWdE[2]*l0)/diff_stretch;
	  else
	    dtde = 0.5*(fddWddE(0,0) - fddWddE(0,2))-fdWdE[2];
	  MixedRank4_3D(eigenvectors[0], eigenvectors[2], fModMat);
	  fModulus.AddScaled(2.0*dtde, fModMat);
	  
	  diff_stretch = l1 - l2;
	  if (fabs(diff_stretch) > small)
	    dtde = (fdWdE[1]*l2 - fdWdE[2]*l1)/diff_stretch;
	  else
	    dtde = 0.5*(fddWddE(1,1) - fddWddE(1,2))-fdWdE[1];;
	  MixedRank4_3D(eigenvectors[1], eigenvectors[2], fModMat);
	  fModulus.AddScaled(2.0*dtde, fModMat);
	}
	/*	cout <<"\n b: "<<fb;
       	cout <<"\n n0: "<<eigenvectors[0];
	cout <<"\n n1: "<<eigenvectors[1];
	if (NumSD() == 3)
	cout <<"\n n2: "<<eigenvectors[2];

	cout <<"\n stretch: "<<eigenstretch; */
	/*	cout <<"\n stress: "<<fStress;
		cout <<"\n modulus: "<<fModulus;*/
	
	return fModulus;
}

const dSymMatrixT& OgdenBaseT::s_ij(void)
{
	/* stretch */
	Compute_b(fb);

	/* spectral decomposition */
       	fSpectralDecomp.SpectralDecomp_Jacobi(fb, false);
       	dWdE(fSpectralDecomp.Eigenvalues(), fdWdE);

       	fStress = fSpectralDecomp.EigsToRank2(fdWdE);

	return(fStress);
}

/* return the pressure associated with the last call to stress */
double OgdenBaseT::Pressure(void) const
{
	return dArrayT::Dot(fSpectralDecomp.Eigenvalues(), fdWdE)/
			sqrt(fSpectralDecomp.Eigenvalues().Product());
}

/*************************************************************************
* Private
*************************************************************************/

/* construct symmetric rank-4 mixed-direction tensor (6.1.44) */
void OgdenBaseT::MixedRank4_2D(const dArrayT& a, const dArrayT& b,
	dMatrixT& rank4_ab) const
{
#if __option(extended_errorcheck)
	if (a.Length() != 2 ||
	    b.Length() != 2 ||
	    rank4_ab.Rows() != 3 ||
	    rank4_ab.Cols() != 3) throw ExceptionT::kSizeMismatch;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11;

//	z1 = A(1.);
//	z2 = A(2.);
//	z3 = B(1.);
//	z4 = B(2.);

	
	z1 = a[0];
	z2 = a[1];
	z3 = b[0];
	z4 = b[1];

       	z5 = z1*z1;
	z6 = z2*z2;
	z7 = z3*z3;
	z8 = 2.*z1*z2*z3*z4;
	z9 = z4*z4;
	z3 = 2.*z3*z4;
	z4 = z3*z5;
	z3 = z3*z6;
	z10 = 2.*z1*z2*z7;
	z11 = 2.*z5*z7;
	z7 = z6*z7;
	z1 = 2.*z1*z2*z9;
	z2 = z5*z9;
	z5 = 2.*z6*z9;
	z4 = z10 + z4;
	z1 = z1 + z3;
	z2 = z2 + z7 + z8;
	z3 = 0.5*z4;
	z1 = 0.5*z1;
	z2 = 0.5*z2;
       

	//{{z11, z8, z3},
	// {z8, z5, z1},
	// {z3, z1, z2}}

	/*	double a1 = a[0];
	double a2 = a[1];
	double b1 = b[0];
	double b2 = b[1];

       	double d11 = 0.5*(a1*b1*a1*b1 + a1*b1*b1*a1 + b1*a1*b1*a1 + b1*a1*a1*b1);
	double d12 = 0.5*(a1*b1*a2*b2 + a1*b1*b2*a2 + b1*a1*b2*a2 + b1*a1*a2*b2); 
	double d13 = 0.5*(a1*b1*a1*b2 + a1*b1*b1*a2 + b1*a1*b1*a2 + b1*a1*a1*b2);
	double d21 = 0.5*(a2*b2*a1*b1 + a2*b2*b1*a1 + b2*a2*b1*a1 + b2*a2*a1*b1);
	double d22 = 0.5*(a2*b2*a2*b2 + a2*b2*b2*a2 + b2*a2*b2*a2 + b2*a2*a2*b2);
	double d23 = 0.5*(a2*b2*a1*b2 + a2*b2*b1*a2 + b2*a2*b1*a2 + b2*a2*a1*b2);
	double d31 = 0.5*(a1*b2*a1*b1 + a1*b2*b1*a1 + b1*a2*b1*a1 + b1*a2*a1*b1);
	double d32 = 0.5*(a1*b2*a2*b2 + a1*b2*b2*a2 + b1*a2*b2*a2 + b1*a2*a2*b2);
	double d33 = 0.5*(a1*b2*a1*b2 + a1*b2*b1*a2 + b1*a2*b1*a2 + b1*a2*a1*b2);

	p[0]=d11;
	p[1]=d21;
	p[2]=d31;
	p[3]=d12;
	p[4]=d22;
	p[5]=d32;
	p[6]=d13;
	p[7]=d23;
	p[3]=d33; */

	double* p = rank4_ab.Pointer();
      	*p++ = z11;
	*p++ = z8;
	*p++ = z3;
	*p++ = z8;
	*p++ = z5;
	*p++ = z1;
	*p++ = z3;
	*p++ = z1;
	*p   = z2;
	
}

void OgdenBaseT::MixedRank4_3D(const dArrayT& a, const dArrayT& b,
	dMatrixT& rank4_ab) const
{
#if __option(extended_errorcheck)
	if (a.Length() != 3 ||
	    b.Length() != 3 ||
	    rank4_ab.Rows() != 6 ||
	    rank4_ab.Cols() != 6) throw ExceptionT::kSizeMismatch;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24;
	double z25, z26, z27, z28, z29, z30, z31, z32, z33, z34, z35, z36;
	double z37, z38, z39, z40, z41;

//	z1 = A(1.);
//	z2 = A(2.);
//	z3 = A(3.);
//	z4 = B(1.);
//	z5 = B(2.);
//	z6 = B(3.);

	z1 = a[0];	
	z2 = a[1];
	z3 = a[2];
	z4 = b[0];
	z5 = b[1];
	z6 = b[2];
	z7 = z1*z1;
	z8 = z2*z2;
	z9 = z3*z3;
	z10 = 2.*z1*z4;
	z11 = z10*z2;
	z12 = z2*z4;
	z13 = z1*z12;
	z14 = z4*z4;
	z15 = z10*z5;
	z15 = z15*z3;
	z16 = z11*z5;
	z17 = z12*z5;
	z18 = 2.*z17;
	z17 = z17*z3;
	z18 = z18*z3;
	z19 = z1*z4*z5;
	z19 = z19*z3;
	z20 = z5*z5;
	z11 = z11*z6;
	z13 = z13*z6;
	z10 = z10*z3*z6;
	z12 = z12*z3*z6;
	z21 = 2.*z12;
	z22 = z1*z2*z5*z6;
	z23 = 2.*z22;
	z24 = z1*z3*z5*z6;
	z25 = 2.*z24;
	z26 = 2.*z2*z3*z5*z6;
	z27 = z6*z6;
	z28 = 2.*z1*z14;
	z29 = 2.*z1*z2;
	z30 = z14*z2;
	z11 = z11 + z15;
	z15 = z1*z20*z3;
	z31 = 2.*z2*z20*z3;
	z18 = z18 + z23;
	z21 = z21 + z25;
	z23 = z1*z2*z27;
	z1 = 2.*z1*z27*z3;
	z25 = 2.*z2*z27*z3;
	z2 = z2*z28;
	z28 = z28*z3;
	z29 = z20*z29;
	z3 = z3*z30;
	z30 = 2.*z14*z7;
	z32 = z20*z7;
	z33 = z27*z7;
	z34 = 2.*z4*z5*z7;
	z35 = 2.*z4*z6*z7;
	z7 = z5*z6*z7;
	z36 = z14*z8;
	z37 = 2.*z20*z8;
	z38 = z27*z8;
	z39 = 2.*z4*z5*z8;
	z40 = z4*z6*z8;
	z8 = 2.*z5*z6*z8;
	z14 = z14*z9;
	z20 = z20*z9;
	z27 = 2.*z27*z9;
	z41 = z4*z5*z9;
	z4 = 2.*z4*z6*z9;
	z5 = 2.*z5*z6*z9;
	z6 = 0.5*z11;
	z9 = 0.5*z18;
	z11 = 0.5*z21;
	z2 = z2 + z34;
	z18 = z28 + z35;
	z3 = z13 + z19 + z3 + z7;
	z7 = z16 + z32 + z36;
	z13 = z29 + z39;
	z15 = z15 + z17 + z22 + z40;
	z8 = z31 + z8;
	z14 = z10 + z14 + z33;
	z17 = z20 + z26 + z38;
	z12 = z12 + z23 + z24 + z41;
	z1 = z1 + z4;
	z4 = z25 + z5;
	z2 = 0.5*z2;
	z5 = 0.5*z18;
	z3 = 0.5*z3;
	z7 = 0.5*z7;
	z13 = 0.5*z13;
	z15 = 0.5*z15;
	z8 = 0.5*z8;
	z14 = 0.5*z14;
	z17 = 0.5*z17;
	z12 = 0.5*z12;
	z1 = 0.5*z1;
	z4 = 0.5*z4;
	
	//{{z30, z16, z10,  z6,  z5,  z2},
	// {z16, z37, z26,  z8,  z9, z13},
	// {z10, z26, z27,  z4,  z1, z11},
	// { z6,  z8,  z4, z17, z12, z15},
	// { z5,  z9,  z1, z12, z14,  z3},
	// { z2, z13, z11, z15,  z3,  z7}}
	
	double* p = rank4_ab.Pointer();
	*p++ = z30;
	*p++ = z16;
	*p++ = z10;
	*p++ = z6;
	*p++ = z5;
	*p++ = z2;
	*p++ = z16;
	*p++ = z37;
	*p++ = z26;
	*p++ = z8;
	*p++ = z9;
	*p++ = z13;
	*p++ = z10;
	*p++ = z26;
	*p++ = z27;
	*p++ = z4;
	*p++ = z1;
	*p++ = z11;
	*p++ = z6;
	*p++ = z8;
	*p++ = z4;
	*p++ = z17;
	*p++ = z12;
	*p++ = z15;
	*p++ = z5;
	*p++ = z9;
	*p++ = z1;
	*p++ = z12;
	*p++ = z14;
	*p++ = z3;
	*p++ = z2;
	*p++ = z13;
	*p++ = z11;
	*p++ = z15;
	*p++ = z3;
	*p  = z7;
}
