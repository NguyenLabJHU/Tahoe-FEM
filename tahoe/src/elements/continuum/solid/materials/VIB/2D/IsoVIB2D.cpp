/* $Id: IsoVIB2D.cpp,v 1.1.1.1 2001-01-29 08:20:24 paklein Exp $ */
/* created: paklein (11/08/1997)                                          */
/* 2D Isotropic VIB solver using spectral decomposition formulation       */

#include "IsoVIB2D.h"

#include <math.h>
#include <iostream.h>
#include "Constants.h"

#include "ElasticT.h"
#include "C1FunctionT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"

/* point generator */
#include "EvenSpacePtsT.h"

/* constructors */
IsoVIB2D::IsoVIB2D(ifstreamT& in, const ElasticT& element):
	FDStructMatT(in, element),
	Material2DT(in, kPlaneStress),
	VIB(in, 2, 2, 3),
	fEigs(2),
	fEigmods(2),
	fSpectral(2),
	fModulus(dSymMatrixT::NumValues(2))
{
	/* point generator */
	fCircle = new EvenSpacePtsT(in);

	/* set tables */
	Construct();
}

/* destructor */
IsoVIB2D::~IsoVIB2D(void)
{
	delete fCircle;
}

/* print parameters */
void IsoVIB2D::Print(ostream& out) const
{
	/* inherited */
	FDStructMatT::Print(out);
	Material2DT::Print(out);
	VIB::Print(out);

	fCircle->Print(out);
}

/* modulus */
const dMatrixT& IsoVIB2D::c_ijkl(void)
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
	double* ps1  = fStressTable(0);
	double* ps2  = fStressTable(1);

	/* modulus */
	double* pc11  = fModuliTable(0);
	double* pc22  = fModuliTable(1);
	double* pc12  = fModuliTable(2);
	
	/* PK2 principal values */	
	double s1 = 0.0;
	double s2 = 0.0;

	fEigmods = 0.0;
	double& c11 = fEigmods(0,0);
	double& c22 = fEigmods(1,1);
	double& c12 = fEigmods(0,1);

	for (int i = 0; i < fLengths.Length(); i++)
	{
		double sfactor = (*pj)*(*pdU)/(*pl);
		double cfactor = (*pj++)*((*pddU++)/((*pl)*(*pl)) -
		                          (*pdU++)/((*pl)*(*pl)*(*pl)));
		pl++;

		s1 += sfactor*(*ps1++);
		s2 += sfactor*(*ps2++);

		c11 += cfactor*(*pc11++);
		c22 += cfactor*(*pc22++);
		c12 += cfactor*(*pc12++);
	}

	/* (material) -> (spatial) (with thickness) */
	double J = sqrt(fEigs[0]*fEigs[1])/fThickness;

	if (fabs(fEigs[0]-fEigs[1]) < kSmall)
	{
		fModulus = 0.0;

		double k = fEigs[0]*fEigs[0]/J;
		fModulus(0,0) = fModulus(1,1) = c11*k;
fModulus(0,1) = fModulus(1,0) = c12*k;
fModulus(2,2) = 0.5*(c11 - c12)*k;

//		fModulus(2,2) = fModulus(0,1) =
//		               fModulus(1,0) = c12*k; //Cauchy symmetry
	}
	else
	{
		fModulus = 0.0;
	
		/* (mat mod) -> (spat mod) */
		c11 *= (fEigs[0]*fEigs[0]/J);
		c22 *= (fEigs[1]*fEigs[1]/J);
		c12 *= (fEigs[0]*fEigs[1]/J);

		/* Cauchy stress principal values in fEigs */
		fEigs[0] *= (s1/J);
		fEigs[1] *= (s2/J);

		/* additional diagonal contribution */
		c11 += 2.0*fEigs[0];
		c22 += 2.0*fEigs[1];

		/* set spectral decomp of b */
		const dSymMatrixT& b_2D = b();
		fSpectral.DecompAndModPrep(b_2D, false);

		/* construct moduli */
		fModulus = fSpectral.EigsToRank4(fEigmods);		
		fModulus.AddScaled(2.0*fEigs[0],fSpectral.SpatialTensor(b_2D, 0));
		fModulus.AddScaled(2.0*fEigs[1],fSpectral.SpatialTensor(b_2D, 1));
	}
	
	return fModulus;
}
	
/* stress */
const dSymMatrixT& IsoVIB2D::s_ij(void)
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

	double *p1  = fStressTable(0);
	double *p2  = fStressTable(1);
	
	/* PK2 principal values */	
	double s1 = 0.0;
	double s2 = 0.0;
	for (int i = 0; i < fLengths.Length(); i++)
	{
		double factor = (*pj++)*(*pdU++)/(*pl++);
		s1 += factor*(*p1++);
		s2 += factor*(*p2++);
	}

	/* PK2 -> Cauchy (with thickness) */
	double J = sqrt(fEigs[0]*fEigs[1])/fThickness;
	fEigs[0] *= (s1/J);
	fEigs[1] *= (s2/J);

	/* set spectral decomp of b */
	fSpectral.SpectralDecomp(b(), false);

	/* build stress */
	return fSpectral.EigsToRank2(fEigs);
}

/* material description */
const dMatrixT& IsoVIB2D::C_IJKL(void)
{
	/* spatial tangent modulus */
	const dMatrixT& modulus = c_ijkl();
	
	/* tranform to material */
	return c_to_C(modulus);  	
}

const dSymMatrixT& IsoVIB2D::S_IJ(void)
{
	/* Cauchy stress */
	const dSymMatrixT& cauchy = s_ij();

	/* convert to PK2 */
	return s_to_S(cauchy);
}

//TEMP
const dSymMatrixT& IsoVIB2D::CurvatureTensor(void)
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
	double* ps1  = fStressTable(0);
	double* ps2  = fStressTable(1);

	/* modulus */
	double* pc11  = fModuliTable(0);
	double* pc22  = fModuliTable(1);
	double* pc12  = fModuliTable(2);
	
	/* PK2 principal values */	
	double s1 = 0.0;
	double s2 = 0.0;

	fEigmods = 0.0;
	double& c11 = fEigmods(0,0);
	double& c22 = fEigmods(1,1);
	double& c12 = fEigmods(0,1);

	for (int i = 0; i < fLengths.Length(); i++)
	{
		double sfactor = (*pj)*(*pdU)/(*pl);
		double cfactor = (*pj++)*((*pddU++)/((*pl)*(*pl)) -
		                          (*pdU++)/((*pl)*(*pl)*(*pl)));
		pl++;

		s1 += sfactor*(*ps1++);
		s2 += sfactor*(*ps2++);

		c11 += cfactor*(*pc11++);
		c22 += cfactor*(*pc22++);
		c12 += cfactor*(*pc12++);
	}
	
	double L_1 = sqrt(fEigs[0]);
	double L_2 = sqrt(fEigs[1]);

	fEigmods[0] = L_1*L_1*c11 + s1;
	fEigmods[1] = L_2*L_2*c22 + s2;
	fEigmods[2] = L_1*L_2*c12;
	
	return fEigmods;
}

/* strain energy density */
double IsoVIB2D::StrainEnergyDensity(void)
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
	
	return energy*fThickness;
}

/***********************************************************************
* Protected
***********************************************************************/

/* print name */
void IsoVIB2D::PrintName(ostream& out) const
{
	/* inherited */
	FDStructMatT::PrintName(out);
	VIB::PrintName(out);

	out << "    Isotropic/Principal Stretch Formulation\n";

	/* integration rule */
	fCircle->PrintName(out);
}

/* strained lengths in terms of the Lagrangian stretch eigenvalues */
void IsoVIB2D::ComputeLengths(const dArrayT& eigs)
{
	double C1 = eigs[0];
	double C2 = eigs[1];

	/* initialize kernel pointers */
	double* pl = fLengths.Pointer();
	double* s1 = fStressTable(0);
	double* s2 = fStressTable(1);
		
	for (int i = 0; i < fLengths.Length(); i++)
		*pl++ = sqrt(C1*(*s1++) + C2*(*s2++));
}

/***********************************************************************
* Private
***********************************************************************/

/* Initialize angle tables */
void IsoVIB2D::Construct(void)
{
	/* fetch points */
	const dArray2DT& points = fCircle->CirclePoints(0.0);
	int numpoints = points.MajorDim();
	
	/* allocate memory */
	Allocate(numpoints);
	
	/* fetch jacobians */
	fjacobian = fCircle->Jacobians();
	
	/* set pointers */
	double *s1 = fStressTable(0);
	double *s2 = fStressTable(1);

	double *c11 = fModuliTable(0);
	double *c22 = fModuliTable(1);
	double *c12 = fModuliTable(2);

	for (int i = 0; i < numpoints; i++)
	{
		/* direction cosines */
		double *xsi = points(i);
		double cosi = xsi[0];
		double sini = xsi[1];
		
		/* stress angle tables */
		s1[i] = cosi*cosi;
		s2[i] = sini*sini;
	
		/* moduli angle tables */
		c11[i] = s1[i]*s1[i];
		c22[i] = s2[i]*s2[i];
		c12[i] = s2[i]*s1[i];
	}
}
