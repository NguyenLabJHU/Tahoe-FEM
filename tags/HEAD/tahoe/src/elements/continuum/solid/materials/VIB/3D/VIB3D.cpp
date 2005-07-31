/* $Id: VIB3D.cpp,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (04/20/1997)                                          */
/* Base class for general 3D probabolistic Cauchy-Born materials.         */

#include "VIB3D.h"

#include <math.h>
#include <iostream.h>

#include "Constants.h"
#include "ExceptionCodes.h"

#include "fstreamT.h"
#include "C1FunctionT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"

/* point generators */
#include "LatLongPtsT.h"
#include "IcosahedralPtsT.h"
#include "FCCPtsT.h"

/* constructors */
VIB3D::VIB3D(ifstreamT& in, const ElasticT& element):
	NL_E_MatT(in, element),
	VIB_E_MatT(in, 3)
{
	/* construct point generator */
	int pointcode;
	in >> pointcode;
	switch(pointcode)
	{
		case SpherePointsT::kLatLong:
			fSphere = new LatLongPtsT(in);
			break;
			
		case SpherePointsT::kIcosahedral:
			fSphere = new IcosahedralPtsT(in);
			break;

		case SpherePointsT::kFCC:
		{
			int num_shells;
			double bond_length;
			in >> num_shells >> bond_length;
			fSphere = new FCCPtsT(num_shells, bond_length);
			break;		
		}	
		default:
		
			throw eBadInputValue;
	}
	
	if (!fSphere) throw(eOutOfMemory);
	
	/* default construction */
	SetAngles(0.0, 0.0);
}

/* destructor */
VIB3D::~VIB3D(void) { delete fSphere; }

/* print parameters */
void VIB3D::Print(ostream& out) const
{
	/* inherited */
	NL_E_MatT::Print(out);
	VIB_E_MatT::Print(out);
	
	fSphere->Print(out);
}

/* set angle offset - for testing onset of amorphous behavior */
void VIB3D::SetAngles(double phi, double theta)
{
	/* fetch points */
	const dArray2DT& points = fSphere->SpherePoints(phi,theta);
	int numpoints = points.MajorDim();
	
	/* allocate memory */
	Allocate(numpoints);
	
	/* fetch jacobians */
	fjacobian = fSphere->Jacobians();
	
	/* set STRESS angle table pointers */
	double *S11, *S22, *S33, *S23, *S13, *S12;
	SetStressPointers3D(S11,S22,S33,S23,S13,S12);
	
	/* set MODULI angle table pointers */
	double *C11, *C12, *C13, *C14, *C15;
	double *C16, *C22, *C23, *C24, *C25;
	double *C26, *C33, *C34, *C35, *C36;
	SetModuliPointers3D(C11, C12, C13, C14, C15,
C16, C22, C23, C24, C25,
C26, C33, C34, C35, C36);

	/* set tables */
	for (int i = 0; i < numpoints; i++)
	{
		/* direction cosines */
		double* xsi = points(i);
		double xsi1 = xsi[0];
		double xsi2 = xsi[1];
		double xsi3 = xsi[2];
			
		/* STRESS angle tables */
		S11[i] = xsi1*xsi1;
		S22[i] = xsi2*xsi2;
		S33[i] = xsi3*xsi3;
		S23[i] = xsi2*xsi3;
		S13[i] = xsi1*xsi3;
		S12[i] = xsi1*xsi2;
	
		/* MODULI angle tables */
		C11[i] = S11[i]*S11[i];
		C12[i] = S11[i]*S22[i];
		C13[i] = S11[i]*S33[i];
		C14[i] = S11[i]*S23[i];
		C15[i] = S11[i]*S13[i];
		C16[i] = S11[i]*S12[i];

		C22[i] = S22[i]*S22[i];
		C23[i] = S22[i]*S33[i];
		C24[i] = S22[i]*S23[i];
		C25[i] = S22[i]*S13[i];
		C26[i] = S22[i]*S12[i];

		C33[i] = S33[i]*S33[i];
		C34[i] = S33[i]*S23[i];
		C35[i] = S33[i]*S13[i];
		C36[i] = S33[i]*S12[i];
	}
}

/***********************************************************************
* Protected
***********************************************************************/

/* print name */
void VIB3D::PrintName(ostream& out) const
{
	/* inherited */
	NL_E_MatT::PrintName(out);
	VIB_E_MatT::PrintName(out);

	fSphere->PrintName(out);
}

/* compute the symetric Cij reduced index matrix */
void VIB3D::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{
	/* fill length table */
	ComputeLengths(E);

	/* derivatives of the potential */
	fPotential->MapDFunction(fLengths, fdU);
	fPotential->MapDDFunction(fLengths, fddU);

	/* initialize */
	moduli = 0.0;

	/* references to the stress components */
	double& C11 = moduli(0,0);
	double& C12 = moduli(0,1);
	double& C13 = moduli(0,2);
	double& C14 = moduli(0,3);
	double& C15 = moduli(0,4);
	double& C16 = moduli(0,5);

	double& C22 = moduli(1,1);
	double& C23 = moduli(1,2);
	double& C24 = moduli(1,3);
	double& C25 = moduli(1,4);
	double& C26 = moduli(1,5);

	double& C33 = moduli(2,2);
	double& C34 = moduli(2,3);
	double& C35 = moduli(2,4);
	double& C36 = moduli(2,5);

	/* initialize kernel pointers */
	double* pl   = fLengths.Pointer();
	double* pdU  = fdU.Pointer();
	double* pddU = fddU.Pointer();
	double* pj   = fjacobian.Pointer();
	
	/* set pointers */
	double *p11, *p12, *p13, *p14, *p15;
	double *p16, *p22, *p23, *p24, *p25;
	double *p26, *p33, *p34, *p35, *p36;
	SetModuliPointers3D(p11, p12, p13, p14, p15,
p16, p22, p23, p24, p25,
p26, p33, p34, p35, p36);
	
	for (int i = 0; i < fLengths.Length(); i++)
	{
		double factor = (*pj++)*((*pddU++)/((*pl)*(*pl)) -
		                         (*pdU++)/((*pl)*(*pl)*(*pl)));
		pl++;

		C11 += factor*(*p11++);
		C12 += factor*(*p12++);
		C13 += factor*(*p13++);
		C14 += factor*(*p14++);
		C15 += factor*(*p15++);
		C16 += factor*(*p16++);

		C22 += factor*(*p22++);
		C23 += factor*(*p23++);
		C24 += factor*(*p24++);
		C25 += factor*(*p25++);
		C26 += factor*(*p26++);

		C33 += factor*(*p33++);
		C34 += factor*(*p34++);
		C35 += factor*(*p35++);
		C36 += factor*(*p36++);
	}

	/* Cauchy symmetry */
	moduli(3,3) = C23;
	moduli(4,4) = C13;
	moduli(5,5) = C12;

	moduli(4,5) = C14;
	moduli(3,5) = C25;
	moduli(3,4) = C36;
	
	/* make symmetric */
	moduli.CopySymmetric();
}

/* compute the symetric 2nd Piola-Kirchhoff reduced index vector */
void VIB3D::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	/* fill length table */
	ComputeLengths(E);

	/* derivatives of the potential */
	fPotential->MapDFunction(fLengths, fdU);

	/* initialize */
	PK2 = 0.0;

	/* references to the PK2 components */
	double& s11 = PK2[0];
	double& s22 = PK2[1];
	double& s33 = PK2[2];
	double& s23 = PK2[3];
	double& s13 = PK2[4];
	double& s12 = PK2[5];

	/* initialize kernel pointers */
	double* pl  = fLengths.Pointer();
	double* pdU = fdU.Pointer();
	double* pj  = fjacobian.Pointer();
	
	/* set pointers */
	double *p11, *p22, *p33, *p23, *p13, *p12;
	SetStressPointers3D(p11,p22,p33,p23,p13,p12);

	for (int i = 0; i < fLengths.Length(); i++)
	{
		double factor = (*pj++)*(*pdU++)/(*pl++);

		s11 += factor*(*p11++);
		s22 += factor*(*p22++);
		s33 += factor*(*p33++);
		s23 += factor*(*p23++);
		s13 += factor*(*p13++);
		s12 += factor*(*p12++);
	}
}

/* returns the strain energy density for the specified strain */
double VIB3D::ComputeEnergyDensity(const dSymMatrixT& E)
{
	return( VIBEnergyDensity(E) );
}
