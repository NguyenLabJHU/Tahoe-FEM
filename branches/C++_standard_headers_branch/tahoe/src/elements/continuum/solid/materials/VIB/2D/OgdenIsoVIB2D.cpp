/* $Id: OgdenIsoVIB2D.cpp,v 1.12 2004-07-15 08:27:45 paklein Exp $ */
/* created: paklein (11/08/1997) */
#include "OgdenIsoVIB2D.h"

#include <math.h>
#include <iostream.h>
#include "toolboxConstants.h"
#include "C1FunctionT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"

/* point generator */
#include "EvenSpacePtsT.h"

using namespace Tahoe;

/* constructor */
OgdenIsoVIB2D::OgdenIsoVIB2D(void):
	ParameterInterfaceT("Ogden_isotropic_VIB_2D"),
	VIB(2, 2, 3),
	fCircle(NULL)
{

}

/* destructor */
OgdenIsoVIB2D::~OgdenIsoVIB2D(void) { delete fCircle; }

/* strain energy density */
double OgdenIsoVIB2D::StrainEnergyDensity(void)
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

/* describe the parameters needed by the interface */
void OgdenIsoVIB2D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	OgdenIsotropicT::DefineParameters(list);
	VIB::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStress);

	/* integration points */
	ParameterT points(ParameterT::Integer, "n_points");
	points.AddLimit(1, LimitT::LowerInclusive);
	list.AddParameter(points);
}

/* information about subordinate parameter lists */
void OgdenIsoVIB2D::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	OgdenIsotropicT::DefineSubs(sub_list);
	VIB::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* OgdenIsoVIB2D::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* sub = OgdenIsotropicT::NewSub(name);
	if (sub) return sub;
	else return VIB::NewSub(name);
}

/* describe the parameters needed by the interface */
void OgdenIsoVIB2D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	OgdenIsotropicT::TakeParameterList(list);
	VIB::TakeParameterList(list);

	/* point generator */
	int points = list.GetParameter("n_points");
	fCircle = new EvenSpacePtsT(points);
	Construct();
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* principal values given principal values of the stretch tensors,
 * i.e., the principal stretches squared */
void OgdenIsoVIB2D::dWdE(const dArrayT& eigenstretch2, dArrayT& eigenstress)
{
	/* stretched bonds */
	ComputeLengths(eigenstretch2);

	/* derivatives of the potential */
	fPotential->MapDFunction(fLengths, fdU);
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
}

void OgdenIsoVIB2D::ddWddE(const dArrayT& eigenstretch2, dArrayT& eigenstress,
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
		double cfactor = (*pj)*((*pddU)/((*pl)*(*pl)) -
		                        (*pdU)/((*pl)*(*pl)*(*pl)));
		pl++; pj++; pddU++; pdU++;

		s0 += sfactor*(*ps0++);
		s1 += sfactor*(*ps1++);

		c00 += cfactor*(*pc00++);
		c11 += cfactor*(*pc11++);
		c01 += cfactor*(*pc01++);
	}
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
	VIB::Dimension(numpoints);
	
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
		const double *xsi = points(i);
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
