/* $Id: IsoVIB2D.cpp,v 1.10 2004-07-15 08:27:45 paklein Exp $ */
/* created: paklein (11/08/1997) */
#include "IsoVIB2D.h"

#include <math.h>
#include "toolboxConstants.h"
#include "C1FunctionT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"

/* point generator */
#include "EvenSpacePtsT.h"

using namespace Tahoe;

/* constructors */
IsoVIB2D::IsoVIB2D(void):
	ParameterInterfaceT("isotropic_VIB_2D"),
	VIB(2, 2, 3),
	fCircle(NULL),
	fSpectral(2)
{

}

/* destructor */
IsoVIB2D::~IsoVIB2D(void) { delete fCircle; }

/* modulus */
const dMatrixT& IsoVIB2D::c_ijkl(void)
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
	double J = sqrt(fEigs[0]*fEigs[1]);

	if (fabs(fEigs[0]-fEigs[1]) < kSmall)
	{
		fModulus = 0.0;

		double k = fEigs[0]*fEigs[0]/J;
		fModulus(0,0) = fModulus(1,1) = c11*k;
		fModulus(0,1) = fModulus(1,0) = c12*k;
		fModulus(2,2) = 0.5*(c11 - c12)*k;
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

		/* set constribution due to b */
		fSpectral.PerturbRoots();
		fSpectral.ModulusPrep(fb);

		/* construct moduli */
		fModulus = fSpectral.EigsToRank4(fEigmods);		
		fModulus.AddScaled(2.0*fEigs[0], fSpectral.SpatialTensor(fb, 0));
		fModulus.AddScaled(2.0*fEigs[1], fSpectral.SpatialTensor(fb, 1));
	}
	
	return fModulus;
}
	
/* stress */
const dSymMatrixT& IsoVIB2D::s_ij(void)
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
	double J = sqrt(fEigs[0]*fEigs[1]);
	fEigs[0] *= (s1/J);
	fEigs[1] *= (s2/J);

	/* build stress */
	return fSpectral.EigsToRank2(fEigs);
}

/* material description */
const dMatrixT& IsoVIB2D::C_IJKL(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fModulus.SetToScaled(Fmat.Det(), PullBack(Fmat, c_ijkl()));
	return fModulus;
}
/**< \todo construct directly in material description */

const dSymMatrixT& IsoVIB2D::S_IJ(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fStress.SetToScaled(Fmat.Det(), PullBack(Fmat, s_ij()));
	return fStress;
}
/**< \todo construct directly in material description */

//TEMP
const dSymMatrixT& IsoVIB2D::CurvatureTensor(void)
{
	/* stretch */
	Compute_b(fb);

	/* principal stretches */
	fb.PrincipalValues(fEigs);

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

/* describe the parameters needed by the interface */
void IsoVIB2D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FSSolidMatT::DefineParameters(list);
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
void IsoVIB2D::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);
	VIB::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* IsoVIB2D::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* sub = FSSolidMatT::NewSub(name);
	if (sub) return sub;
	else return VIB::NewSub(name);
}

/* accept parameter list */
void IsoVIB2D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSSolidMatT::TakeParameterList(list);
	VIB::TakeParameterList(list);

	/* dimension work space */
	fEigs.Dimension(2);
	fEigmods.Dimension(2);
	fb.Dimension(2);
	fModulus.Dimension(dSymMatrixT::NumValues(2));
	fStress.Dimension(2);

	/* point generator */
	int points = list.GetParameter("n_points");
	fCircle = new EvenSpacePtsT(points);
	Construct();
}

/***********************************************************************
 * Protected
 ***********************************************************************/

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
	VIB::Dimension(numpoints);
	
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
		const double *xsi = points(i);
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
