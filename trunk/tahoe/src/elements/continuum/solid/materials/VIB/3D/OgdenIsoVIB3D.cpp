/* $Id: OgdenIsoVIB3D.cpp,v 1.9 2003-11-21 22:46:38 paklein Exp $ */
/* created: paklein (11/08/1997) */
#include "OgdenIsoVIB3D.h"

#include <math.h>
#include <iostream.h>
#include "toolboxConstants.h"
#include "C1FunctionT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "ifstreamT.h"

/* point generators */
#include "LatLongPtsT.h"
#include "IcosahedralPtsT.h"

using namespace Tahoe;

/* constructors */
OgdenIsoVIB3D::OgdenIsoVIB3D(ifstreamT& in, const FSMatSupportT& support):
	OgdenIsotropicT(in, support),
	VIB(in, 3, 3, 6),
	fSphere(NULL)
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
			throw ExceptionT::kBadInputValue;
	}
	if (!fSphere) throw ExceptionT::kOutOfMemory;

	/* set tables */
	Construct();
}

/* destructor */
OgdenIsoVIB3D::~OgdenIsoVIB3D(void) { delete fSphere; }

/* print parameters */
void OgdenIsoVIB3D::Print(ostream& out) const
{
	/* inherited */
	OgdenIsotropicT::Print(out);
	VIB::Print(out);

	fSphere->Print(out);
}

/* print name */
void OgdenIsoVIB3D::PrintName(ostream& out) const
{
	/* inherited */
	OgdenIsotropicT::PrintName(out);
	VIB::PrintName(out);
	out << "    Odgen principal stretch formulation\n";

	/* integration rule */
	fSphere->PrintName(out);
}

/* strain energy density */
double OgdenIsoVIB3D::StrainEnergyDensity(void)
{
	/* stretch */
	Compute_C(fC);

	/* principal stretches */
	fC.PrincipalValues(fEigs);

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

/* principal values given principal values of the stretch tensors,
 * i.e., the principal stretches squared */
void OgdenIsoVIB3D::dWdE(const dArrayT& eigenstretch2, dArrayT& eigenstress)
{
	/* stretched bonds */
	ComputeLengths(eigenstretch2);

	/* derivatives of the potential */
	fPotential->MapDFunction(fLengths, fdU);

	/* initialize kernel pointers */
	double* pdU = fdU.Pointer();
	double* pl  = fLengths.Pointer();
	double* pj  = fjacobian.Pointer();

	double *p0  = fStressTable(0);
	double *p1  = fStressTable(1);
	double *p2  = fStressTable(2);
	
	/* PK2 principal values */	
	double& s0 = eigenstress[0] = 0.0;
	double& s1 = eigenstress[1] = 0.0;
	double& s2 = eigenstress[2] = 0.0;
	for (int i = 0; i < fLengths.Length(); i++)
	{
		double factor = (*pj++)*(*pdU++)/(*pl++);
		s0 += factor*(*p0++);
		s1 += factor*(*p1++);
		s2 += factor*(*p2++);
	}
}

void OgdenIsoVIB3D::ddWddE(const dArrayT& eigenstretch2, dArrayT& eigenstress,
	dSymMatrixT& eigenmod)
{
	/* stretched bonds */
	ComputeLengths(eigenstretch2);

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
	double* ps2  = fStressTable(2);

	/* modulus */
	double* pc00  = fModuliTable(0);
	double* pc11  = fModuliTable(1);
	double* pc22  = fModuliTable(2);

	double* pc12  = fModuliTable(3);
	double* pc02  = fModuliTable(4);
	double* pc01  = fModuliTable(5);
	
	/* PK2 principal values */	
	double& s0 = eigenstress[0] = 0.0;
	double& s1 = eigenstress[1] = 0.0;
	double& s2 = eigenstress[2] = 0.0;

	double& c00 = eigenmod(0,0) = 0.0;
	double& c11 = eigenmod(1,1) = 0.0;
	double& c22 = eigenmod(2,2) = 0.0;

	double& c12 = eigenmod(1,2) = 0.0;
	double& c02 = eigenmod(0,2) = 0.0;
	double& c01 = eigenmod(0,1) = 0.0;
	for (int i = 0; i < fLengths.Length(); i++)
	{
		double sfactor = (*pj)*(*pdU)/(*pl);
		double cfactor = (*pj++)*((*pddU++)/((*pl)*(*pl)) -
		                          (*pdU++)/((*pl)*(*pl)*(*pl)));
		pl++;

		s0 += sfactor*(*ps0++);
		s1 += sfactor*(*ps1++);
		s2 += sfactor*(*ps2++);

		c00 += cfactor*(*pc00++);
		c11 += cfactor*(*pc11++);
		c22 += cfactor*(*pc22++);

		c12 += cfactor*(*pc12++);
		c02 += cfactor*(*pc02++);
		c01 += cfactor*(*pc01++);
	}
}

/* strained lengths in terms of the Lagrangian stretch eigenvalues */
void OgdenIsoVIB3D::ComputeLengths(const dArrayT& eigs)
{
	double C0 = eigs[0];
	double C1 = eigs[1];
	double C2 = eigs[2];

	/* initialize kernel pointers */
	double* pl = fLengths.Pointer();
	double* s0 = fStressTable(0);
	double* s1 = fStressTable(1);
	double* s2 = fStressTable(2);
	for (int i = 0; i < fLengths.Length(); i++)
		*pl++ = sqrt(C0*(*s0++) + C1*(*s1++) + C2*(*s2++));
}

/***********************************************************************
* Private
***********************************************************************/

/* Initialize angle tables */
void OgdenIsoVIB3D::Construct(void)
{
	/* fetch points */
	const dArray2DT& points = fSphere->SpherePoints(0.0, 0.0);
	int numpoints = points.MajorDim();
	
	/* allocate memory */
	Dimension(numpoints);
	
	/* fetch jacobians */
	fjacobian = fSphere->Jacobians();
	
	/* set pointers */
	double *s0 = fStressTable(0);
	double *s1 = fStressTable(1);
	double *s2 = fStressTable(2);

	double *c00 = fModuliTable(0);
	double *c11 = fModuliTable(1);
	double *c22 = fModuliTable(2);

	double *c12 = fModuliTable(3);
	double *c02 = fModuliTable(4);
	double *c01 = fModuliTable(5);

	for (int i = 0; i < numpoints; i++)
	{
		/* direction cosines */
		const double *xsi = points(i);

		double xsi0 = xsi[0];
		double xsi1 = xsi[1];
		double xsi2 = xsi[2];
		
		/* stress angle tables */
		s0[i] = xsi0*xsi0;
		s1[i] = xsi1*xsi1;
		s2[i] = xsi2*xsi2;
	
		/* moduli angle tables */
		c00[i] = s0[i]*s0[i];
		c11[i] = s1[i]*s1[i];
		c22[i] = s2[i]*s2[i];

		c12[i] = s1[i]*s2[i];
		c02[i] = s0[i]*s2[i];
		c01[i] = s0[i]*s1[i];
	}
}
