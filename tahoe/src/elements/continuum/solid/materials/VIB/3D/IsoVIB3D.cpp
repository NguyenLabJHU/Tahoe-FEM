/* $Id: IsoVIB3D.cpp,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (03/15/1998)                                          */
/* 3D Isotropic VIB solver using spectral decomposition formulation       */

#include "IsoVIB3D.h"

#include <math.h>
#include <iostream.h>
#include "Constants.h"

#include "ElasticT.h"
#include "fstreamT.h"
#include "C1FunctionT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"

/* point generators */
#include "LatLongPtsT.h"
#include "IcosahedralPtsT.h"

/* constructors */
IsoVIB3D::IsoVIB3D(ifstreamT& in, const ElasticT& element):
	FDStructMatT(in, element),
	VIB(in, 3, 3, 6),
	fEigs(3),
	fEigmods(3),
	fSpectral(3),
	fModulus(dSymMatrixT::NumValues(3))
{	
	/* construct point generator */
	int gencode;
	in >> gencode;
	switch (gencode)
	{
		case SpherePointsT::kLatLong:
			fSphere = new LatLongPtsT(in);
			break;
	
		case SpherePointsT::kIcosahedral:
			fSphere = new IcosahedralPtsT(in);
			break;
			
		default:
			throw eBadInputValue;
	}
	if (!fSphere) throw eOutOfMemory;

	/* set tables */
	Construct();
}

/* destructor */
IsoVIB3D::~IsoVIB3D(void) { delete fSphere; }

/* print parameters */
void IsoVIB3D::Print(ostream& out) const
{
	/* inherited */
	FDStructMatT::Print(out);
	VIB::Print(out);

	fSphere->Print(out);
}

/* print name */
void IsoVIB3D::PrintName(ostream& out) const
{
	/* inherited */
	FDStructMatT::PrintName(out);
	VIB::PrintName(out);

	out << "    Isotropic/Principal Stretch Formulation\n";

	/* integration rule */
	fSphere->PrintName(out);
}

/* modulus */
const dMatrixT& IsoVIB3D::c_ijkl(void)
{
	/* principal stretches */
	C().PrincipalValues(fEigs);

	/* stretched bonds */
	ComputeLengths(fEigs);

	/* derivatives of the potential */
	fPotential->MapDFunction(fLengths, fdU);
	fPotential->MapDDFunction(fLengths, fddU);	

	/* initialize kernel pointers */
	double* pdU  = fdU.Pointer();
	double* pddU = fddU.Pointer();
	double* pl   = fLengths.Pointer();
	double* pj   = fjacobian.Pointer();

	/* stress */
	double *ps1  = fStressTable(0);
	double *ps2  = fStressTable(1);
	double *ps3  = fStressTable(2);

	/* modulus */
	double *pc11  = fModuliTable(0);
	double *pc22  = fModuliTable(1);
	double *pc33  = fModuliTable(2);

	double *pc23  = fModuliTable(3);
	double *pc13  = fModuliTable(4);
	double *pc12  = fModuliTable(5);
	
	/* PK2 principal values */	
	double s1 = 0.0;
	double s2 = 0.0;
	double s3 = 0.0;

	fEigmods = 0.0;
	double& c11 = fEigmods(0,0);
	double& c22 = fEigmods(1,1);
	double& c33 = fEigmods(2,2);

	double& c23 = fEigmods(1,2);
	double& c13 = fEigmods(0,2);
	double& c12 = fEigmods(0,1);

	for (int i = 0; i < fLengths.Length(); i++)
	{
		double sfactor = (*pj)*(*pdU)/(*pl);
		double cfactor = (*pj++)*((*pddU++)/((*pl)*(*pl)) -
		                          (*pdU++)/((*pl)*(*pl)*(*pl)));	
		pl++;
	
		s1 += sfactor*(*ps1++);
		s2 += sfactor*(*ps2++);
		s3 += sfactor*(*ps3++);

		c11 += cfactor*(*pc11++);
		c22 += cfactor*(*pc22++);
		c33 += cfactor*(*pc33++);

		c23 += cfactor*(*pc23++);
		c13 += cfactor*(*pc13++);
		c12 += cfactor*(*pc12++);
	}

	/* (material) -> (spatial) */
	double J = sqrt(fEigs.Product());
	if (J < kSmall) throw eBadJacobianDet;
	
	c11 *= (fEigs[0]*fEigs[0]/J);
	c22 *= (fEigs[1]*fEigs[1]/J);
	c33 *= (fEigs[2]*fEigs[2]/J);

	c23 *= (fEigs[1]*fEigs[2]/J);
	c13 *= (fEigs[0]*fEigs[2]/J);
	c12 *= (fEigs[0]*fEigs[1]/J);

	//TEMP - need special cases for equibiaxial cases
	// trap nearly repeated roots
	double d01 = fEigs[0] - fEigs[1];
	double d02 = fEigs[0] - fEigs[2];
	double d12 = fEigs[1] - fEigs[2];

	/* special treatment of equitriaxial */
//	if (fabs(fEigs[0]-fEigs[1]) < kSmall &&
//	    fabs(fEigs[1]-fEigs[2]) < kSmall &&
//	    fabs(fEigs[2]-fEigs[0]) < kSmall)
	if (fabs(d01*d02) < kSmall ||
	    fabs(d01*d12) < kSmall ||
	    fabs(d02*d12) < kSmall) // TEMP catch all degenerate cases
	{
		fModulus = 0.0;
		
		fModulus(0,0) = c11;
		fModulus(1,1) = c22;
		fModulus(2,2) = c33;

		/* using Cauchy symmetry */
		fModulus(3,3) = fModulus(1,2) =
		               fModulus(2,1) = c23;
		fModulus(4,4) = fModulus(0,2) =
		               fModulus(2,0) = c13;
		fModulus(5,5) = fModulus(0,1) =
		               fModulus(1,0) = c12;
	}
	/* 3 distinct roots */
	else
	{
		/* Cauchy stress principal values */
		fEigs[0] *= (s1/J);
		fEigs[1] *= (s2/J);
		fEigs[2] *= (s3/J);

		/* additional diagonal contribution */
		c11 += 2.0*fEigs[0];
		c22 += 2.0*fEigs[1];
		c33 += 2.0*fEigs[2];

		/* set spectral decomp of b */
		const dSymMatrixT& b_3D = b();
		fSpectral.DecompAndModPrep(b_3D, true);

		/* construct moduli */		
		fModulus = fSpectral.EigsToRank4(fEigmods);
		fModulus.AddScaled(2.0*fEigs[0],fSpectral.SpatialTensor(b_3D, 0));
		fModulus.AddScaled(2.0*fEigs[1],fSpectral.SpatialTensor(b_3D, 1));
		fModulus.AddScaled(2.0*fEigs[2],fSpectral.SpatialTensor(b_3D, 2));
	}
	
	return fModulus;
}
	
/* stress */
const dSymMatrixT& IsoVIB3D::s_ij(void)
{
	/* principal stretches */
	C().PrincipalValues(fEigs);

	/* stretched bonds */
	ComputeLengths(fEigs);

	/* derivatives of the potential */
	fPotential->MapDFunction(fLengths, fdU);

	/* initialize kernel pointers */
	double* pdU = fdU.Pointer();
	double* pl  = fLengths.Pointer();
	double* pj  = fjacobian.Pointer();

	double* p1  = fStressTable(0);
	double* p2  = fStressTable(1);
	double* p3  = fStressTable(2);
	
	/* PK2 principal values */	
	double s1 = 0.0;
	double s2 = 0.0;
	double s3 = 0.0;
	for (int i = 0; i < fLengths.Length(); i++)
	{
		double factor = (*pj++)*(*pdU++)/(*pl++);
		
		s1 += factor*(*p1++);
		s2 += factor*(*p2++);
		s3 += factor*(*p3++);
	}

	/* PK2 -> Cauchy (with thickness) */
	double J = sqrt(fEigs.Product());
	if (J < kSmall) throw eBadJacobianDet;
	
	fEigs[0] *= (s1/J);
	fEigs[1] *= (s2/J);
	fEigs[2] *= (s3/J);

	/* set spectral decomp of b */
	fSpectral.SpectralDecomp(b(), true);

	/* build stress */
	return fSpectral.EigsToRank2(fEigs);
}

/* material description */
const dMatrixT& IsoVIB3D::C_IJKL(void)
{
	/* could construct directly in material description */
	return c_to_C(IsoVIB3D::c_ijkl());
}

const dSymMatrixT& IsoVIB3D::S_IJ(void)
{
	/* could construct directly in material description */
	return  s_to_S(IsoVIB3D::s_ij());
}

/* strain energy density */
double IsoVIB3D::StrainEnergyDensity(void)
{
	/* principal stretches */
	C().PrincipalValues(fEigs);

	/* stretched bonds */
	ComputeLengths(fEigs);

	/* update potential table */
	fPotential->MapFunction(fLengths,fU);

	/* sum contributions */
	double  energy = 0.0;
	double* pU = fU.Pointer();
	double* pj = fjacobian.Pointer();	
	for (int i = 0; i < fLengths.Length(); i++)
		energy += (*pU++)*(*pj++);
	
	return energy;
}

/***********************************************************************
* Protected
***********************************************************************/

/* strained lengths in terms of the Lagrangian stretch eigenvalues */
void IsoVIB3D::ComputeLengths(const dArrayT& eigs)
{
	double C1 = eigs[0];
	double C2 = eigs[1];
	double C3 = eigs[2];

	/* initialize kernel pointers */
	double* pl = fLengths.Pointer();
	double* s1 = fStressTable(0);
	double* s2 = fStressTable(1);
	double* s3 = fStressTable(2);
		
	for (int i = 0; i < fLengths.Length(); i++)
		*pl++ = sqrt(C1*(*s1++) + C2*(*s2++) + C3*(*s3++));
}

/***********************************************************************
* Private
***********************************************************************/

/* initialize angle tables */
void IsoVIB3D::Construct(void)
{
	/* fetch points */
	const dArray2DT& points = fSphere->SpherePoints(0.0, 0.0);
	int numpoints = points.MajorDim();
	
	/* allocate memory */
	Allocate(numpoints);
	
	/* fetch jacobians */
	fjacobian = fSphere->Jacobians();
	
	/* set pointers */
	double *s1 = fStressTable(0);
	double *s2 = fStressTable(1);
	double *s3 = fStressTable(2);

	double *c11 = fModuliTable(0);
	double *c22 = fModuliTable(1);
	double *c33 = fModuliTable(2);

	double *c23 = fModuliTable(3);
	double *c13 = fModuliTable(4);
	double *c12 = fModuliTable(5);

	for (int i = 0; i < numpoints; i++)
	{
		/* direction cosines */
		double *xsi = points(i);

		double xsi1 = xsi[0];
		double xsi2 = xsi[1];
		double xsi3 = xsi[2];
		
		/* stress angle tables */
		s1[i] = xsi1*xsi1;
		s2[i] = xsi2*xsi2;
		s3[i] = xsi3*xsi3;
	
		/* moduli angle tables */
		c11[i] = s1[i]*s1[i];
		c22[i] = s2[i]*s2[i];
		c33[i] = s3[i]*s3[i];

		c23[i] = s2[i]*s3[i];
		c13[i] = s1[i]*s3[i];
		c12[i] = s1[i]*s2[i];
	}
}
