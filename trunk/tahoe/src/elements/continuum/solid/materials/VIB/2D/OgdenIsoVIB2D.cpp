/* $Id: OgdenIsoVIB2D.cpp,v 1.1.1.1 2001-01-29 08:20:24 paklein Exp $ */
/* created: paklein (11/08/1997)                                          */
/* 2D Isotropic VIB using Ogden's spectral formulation                    */

#include "OgdenIsoVIB2D.h"

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
OgdenIsoVIB2D::OgdenIsoVIB2D(ifstreamT& in, const ElasticT& element):
	OgdenIsotropicT(in, element),
	Material2DT(in, kPlaneStress),
	VIB(in, 2, 2, 3),
	fCircle(NULL)
{
	/* point generator */
	fCircle = new EvenSpacePtsT(in);

	/* set tables */
	Construct();
}

/* destructor */
OgdenIsoVIB2D::~OgdenIsoVIB2D(void)
{
	delete fCircle;
}

/* print parameters */
void OgdenIsoVIB2D::Print(ostream& out) const
{
	/* inherited */
	OgdenIsotropicT::Print(out);
	Material2DT::Print(out);
	VIB::Print(out);

	fCircle->Print(out);
}

/* print name */
void OgdenIsoVIB2D::PrintName(ostream& out) const
{
	/* inherited */
	OgdenIsotropicT::PrintName(out);
	VIB::PrintName(out);
	out << "    Odgen principal stretch formulation\n";

	/* integration rule */
	fCircle->PrintName(out);
}

/* strain energy density */
double OgdenIsoVIB2D::StrainEnergyDensity(void)
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

/* principal values given principal stretches */
void OgdenIsoVIB2D::dWdE(const dArrayT& eigenstretch, dArrayT& eigenstress)
{
	/* stretched bonds */
	ComputeLengths(eigenstretch);

	/* derivatives of the potential */
	fPotential->MapDFunction(fLengths, fdU);

	/* initialize kernel pointers */
	double* pdU = fdU.Pointer();
	double* pl  = fLengths.Pointer();
	double* pj  = fjacobian.Pointer();

	double *p0  = fStressTable(0);
	double *p1  = fStressTable(1);
	
	/* PK2 principal values */	
	double& s0 = eigenstress[0] = 0.0;
	double& s1 = eigenstress[1] = 0.0;
	for (int i = 0; i < fLengths.Length(); i++)
	{
		double factor = (*pj++)*(*pdU++)/(*pl++);
		s0 += factor*(*p0++);
		s1 += factor*(*p1++);
	}

	/* thickness */
	s0 *= fThickness;
	s1 *= fThickness;
}

void OgdenIsoVIB2D::ddWddE(const dArrayT& eigenstretch, dArrayT& eigenstress,
	dSymMatrixT& eigenmod)
{
	/* stretched bonds */
	ComputeLengths(eigenstretch);

	/* derivatives of the potential */
	fPotential->MapDFunction(fLengths, fdU);
	fPotential->MapDDFunction(fLengths, fddU);	

	/* initialize kernel pointers */
	double* pdU  = fdU.Pointer();
	double* pddU = fddU.Pointer();
	double* pl   = fLengths.Pointer();
	double* pj   = fjacobian.Pointer();

	/* stress */
	double* ps0  = fStressTable(0);
	double* ps1  = fStressTable(1);

	/* modulus */
	double* pc00  = fModuliTable(0);
	double* pc11  = fModuliTable(1);
	double* pc01  = fModuliTable(2);
	
	/* PK2 principal values */	
	double& s0 = eigenstress[0] = 0.0;
	double& s1 = eigenstress[1] = 0.0;

	double& c00 = eigenmod(0,0) = 0.0;
	double& c11 = eigenmod(1,1) = 0.0;
	double& c01 = eigenmod(0,1) = 0.0;
	for (int i = 0; i < fLengths.Length(); i++)
	{
		double sfactor = (*pj)*(*pdU)/(*pl);
		double cfactor = (*pj++)*((*pddU++)/((*pl)*(*pl)) -
		                          (*pdU++)/((*pl)*(*pl)*(*pl)));
		pl++;

		s0 += sfactor*(*ps0++);
		s1 += sfactor*(*ps1++);

		c00 += cfactor*(*pc00++);
		c11 += cfactor*(*pc11++);
		c01 += cfactor*(*pc01++);
	}

	/* thickness */
	s0  *= fThickness;
	s1  *= fThickness;
	c00 *= fThickness;
	c11 *= fThickness;
	c01 *= fThickness;
}

/* strained lengths in terms of the Lagrangian stretch eigenvalues */
void OgdenIsoVIB2D::ComputeLengths(const dArrayT& eigs)
{
	double C0 = eigs[0];
	double C1 = eigs[1];

	/* initialize kernel pointers */
	double* pl = fLengths.Pointer();
	double* s0 = fStressTable(0);
	double* s1 = fStressTable(1);
		
	for (int i = 0; i < fLengths.Length(); i++)
		*pl++ = sqrt(C0*(*s0++) + C1*(*s1++));
}

/***********************************************************************
* Private
***********************************************************************/

/* Initialize angle tables */
void OgdenIsoVIB2D::Construct(void)
{
	/* fetch points */
	const dArray2DT& points = fCircle->CirclePoints(0.0);
	int numpoints = points.MajorDim();
	
	/* allocate memory */
	Allocate(numpoints);
	
	/* fetch jacobians */
	fjacobian = fCircle->Jacobians();
	
	/* set pointers */
	double *s0 = fStressTable(0);
	double *s1 = fStressTable(1);

	double *c00 = fModuliTable(0);
	double *c11 = fModuliTable(1);
	double *c01 = fModuliTable(2);

	for (int i = 0; i < numpoints; i++)
	{
		/* direction cosines */
		double *xsi = points(i);
		double cosi = xsi[0];
		double sini = xsi[1];
		
		/* stress angle tables */
		s0[i] = cosi*cosi;
		s1[i] = sini*sini;
	
		/* moduli angle tables */
		c00[i] = s0[i]*s0[i];
		c11[i] = s1[i]*s1[i];
		c01[i] = s0[i]*s1[i];
	}
}
