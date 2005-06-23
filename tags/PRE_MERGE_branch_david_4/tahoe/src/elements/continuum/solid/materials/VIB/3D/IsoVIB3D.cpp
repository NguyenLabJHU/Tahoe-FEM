/* $Id: IsoVIB3D.cpp,v 1.11 2004-07-15 08:27:51 paklein Exp $ */
/* created: paklein (03/15/1998) */
#include "IsoVIB3D.h"

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
IsoVIB3D::IsoVIB3D(void):
	ParameterInterfaceT("isotropic_VIB"),
	VIB(3, 3, 6),
	fSpectral(3),
	fSphere(NULL)
{	

}

/* destructor */
IsoVIB3D::~IsoVIB3D(void) { delete fSphere; }

/* modulus */
const dMatrixT& IsoVIB3D::c_ijkl(void)
{
	/* stretch */
	Compute_b(fb);
	
	/* compute spectral decomposition */
	fSpectral.SpectralDecomp_Jacobi(fb, false);
	fEigs = fSpectral.Eigenvalues();

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
	if (J <= 0.0) throw ExceptionT::kBadJacobianDet;
	
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

		/* set constribution due to b */
		fSpectral.PerturbRoots();
		fSpectral.ModulusPrep(fb);

		/* construct moduli */		
		fModulus = fSpectral.EigsToRank4(fEigmods);
		fModulus.AddScaled(2.0*fEigs[0],fSpectral.SpatialTensor(fb, 0));
		fModulus.AddScaled(2.0*fEigs[1],fSpectral.SpatialTensor(fb, 1));
		fModulus.AddScaled(2.0*fEigs[2],fSpectral.SpatialTensor(fb, 2));
	}
	
	return fModulus;
}
	
/* stress */
const dSymMatrixT& IsoVIB3D::s_ij(void)
{
	/* stretch */
	Compute_b(fb);

	/* compute spectral decomposition */
	fSpectral.SpectralDecomp_Jacobi(fb, false);
	fEigs = fSpectral.Eigenvalues();

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
	if (J <= kSmall) throw ExceptionT::kBadJacobianDet;
	
	fEigs[0] *= (s1/J);
	fEigs[1] *= (s2/J);
	fEigs[2] *= (s3/J);

	/* build stress */
	return fSpectral.EigsToRank2(fEigs);
}

/* material description */
const dMatrixT& IsoVIB3D::C_IJKL(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fModulus.SetToScaled(Fmat.Det(), PullBack(Fmat, IsoVIB3D::c_ijkl()));
	return fModulus;
}
/**< \todo construct directly in material description */

const dSymMatrixT& IsoVIB3D::S_IJ(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fStress.SetToScaled(Fmat.Det(), PullBack(Fmat, IsoVIB3D::s_ij()));
	return fStress;
}
/**< \todo construct directly in material description */

/* strain energy density */
double IsoVIB3D::StrainEnergyDensity(void)
{
	/* stretch */
	Compute_b(fb);

	/* principal stretches */
	fb.PrincipalValues(fEigs);

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
void IsoVIB3D::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);
	VIB::DefineSubs(sub_list);

	/* choice of integration schemes */
	sub_list.AddSub("sphere_integration_choice", ParameterListT::Once, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* IsoVIB3D::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* sub = FSSolidMatT::NewSub(name);
	if (sub) 
		return sub;
	else if (name == "sphere_integration_choice")
	{
		/* use other VIB material to construct point generator */
		VIB3D vib;
		return vib.NewSub(name);
	}	
	else /* inherited */
		return VIB::NewSub(name);
}

/* accept parameter list */
void IsoVIB3D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSSolidMatT::TakeParameterList(list);
	VIB::TakeParameterList(list);

	/* dimension work space */
	fEigs.Dimension(3);
	fEigmods.Dimension(3);
	fb.Dimension(3);
	fModulus.Dimension(dSymMatrixT::NumValues(3));
	fStress.Dimension(3);

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
		ExceptionT::GeneralFail("IsoVIB3D::TakeParameterList", "unrecognized point scheme \"%s\"", points.Name().Pointer());	

	/* set tables */
	Construct();	
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
	VIB::Dimension(numpoints);
	
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
		const double *xsi = points(i);

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
