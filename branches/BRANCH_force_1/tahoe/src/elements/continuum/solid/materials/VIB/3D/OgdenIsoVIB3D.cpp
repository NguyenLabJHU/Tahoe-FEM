/* $Id: OgdenIsoVIB3D.cpp,v 1.10 2004-07-15 08:27:51 paklein Exp $ */
/* created: paklein (11/08/1997) */
#include "OgdenIsoVIB3D.h"

#include <math.h>
#include <iostream.h>
#include "toolboxConstants.h"
#include "C1FunctionT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"


/* point generators */
#include "VIB3D.h"
#include "LatLongPtsT.h"
#include "IcosahedralPtsT.h"
#include "FCCPtsT.h"

using namespace Tahoe;

/* constructor */
OgdenIsoVIB3D::OgdenIsoVIB3D(void):
	ParameterInterfaceT("Ogden_isotropic_VIB"),
	VIB(3, 3, 6),
	fSphere(NULL)
{

}

/* destructor */
OgdenIsoVIB3D::~OgdenIsoVIB3D(void) { delete fSphere; }

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

/* information about subordinate parameter lists */
void OgdenIsoVIB3D::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	OgdenIsotropicT::DefineSubs(sub_list);
	VIB::DefineSubs(sub_list);

	/* choice of integration schemes */
	sub_list.AddSub("sphere_integration_choice", ParameterListT::Once, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* OgdenIsoVIB3D::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* sub = OgdenIsotropicT::NewSub(name);
	if (sub) 
		return sub;
	else if (name == "sphere_integration_choice")
	{
		/* use other VIB material to construct point generator */
		VIB3D vib;
		return vib.NewSub(name);
	}
	else 
		return VIB::NewSub(name);
}

/* accept parameter list */
void OgdenIsoVIB3D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	OgdenIsotropicT::TakeParameterList(list);
	VIB::TakeParameterList(list);

	/* use other VIB material to construct integration rule */
	VIB3D vib;
	const ParameterListT& points = list.GetListChoice(vib, "sphere_integration_choice");
	if (points.Name() == "latitude_longitude")
	{
		int n_phi = points.GetParameter("n_phi");
		int n_theta = points.GetParameter("n_theta");
		fSphere = new LatLongPtsT(n_phi, n_theta);
	}
	else if (points.Name() == "icosahedral")
	{
		int np = points.GetParameter("points");
		fSphere = new IcosahedralPtsT(np);
	}
	else if (points.Name() == "fcc_points")
	{
		int num_shells = points.GetParameter("shells");
		double bond_length = points.GetParameter("nearest_neighbor_distance");
		fSphere = new FCCPtsT(num_shells, bond_length);
	}
	else
		ExceptionT::GeneralFail("OgdenIsoVIB3D::TakeParameterList", "unrecognized point scheme \"%s\"", points.Name().Pointer());	

	/* set tables */
	Construct();	
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
	VIB::Dimension(numpoints);
	
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
