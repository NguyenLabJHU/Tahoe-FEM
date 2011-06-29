/* $Id: VIB2D.cpp,v 1.12 2005-03-16 10:20:48 paklein Exp $ */
/* created: paklein (04/09/1997) */
#include "VIB2D.h"

#include <math.h>
#include <iostream.h>

#include "toolboxConstants.h"
#include "C1FunctionT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "ParameterContainerT.h"

/* point generators */
#include "CirclePointsT.h"
#include "GaussPtsT.h"
#include "EvenSpacePtsT.h"

using namespace Tahoe;

/* constants */
const double Pi = acos(-1.0);

/* constructor */
VIB2D::VIB2D(void):
	ParameterInterfaceT("VIB_2D"),
	VIB_E_MatT(2),
	fCircle(NULL)
{
	/* set default value */
	fConstraint = kPlaneStress;
}

/* destructor */
VIB2D::~VIB2D(void) { delete fCircle; }

/* set angle offset - for testing onset of amorphous behavior */
void VIB2D::SetAngle(double angleoffset)
{
	/* fetch points */
	const dArray2DT& points = fCircle->CirclePoints(angleoffset);
	int numpoints = points.MajorDim();
	
	/* allocate memory */
	VIB::Dimension(numpoints);
	
	/* fetch jacobians */
	fjacobian = fCircle->Jacobians();
	
	/* set pointers */
	double *s11, *s22, *s12;	
	SetStressPointers2D(s11, s22, s12);

	double *c11, *c22, *c26, *c16, *c12;	
	SetModuliPointers2D(c11, c22, c26, c16, c12);
	
	for (int i = 0; i < numpoints; i++)
	{
		/* direction cosines */
		const double *xsi = points(i);
		double cosi = xsi[0];
		double sini = xsi[1];
		
		/* stress angle tables */
		s11[i] = cosi*cosi;
		s22[i] = sini*sini;
		s12[i] = sini*cosi; 	
	
		/* moduli angle tables */
		c11[i] = s11[i]*s11[i];
		c22[i] = s22[i]*s22[i];
		c26[i] = s22[i]*s12[i];
		c16[i] = s11[i]*s12[i];
		c12[i] = s22[i]*s11[i];
	}
}

//TEMP - microscopic test of stability
void VIB2D::Perturb(dArrayT& dU, double eps)
{
	/* sampling points (even number) */
	int  npts = (2*dU.Length()+1)/2;
	double da = 2.0*Pi/npts;

// if the number of lengths is ODD need to make bond vectors
// for twice as many bonds since deformation will be different
// for each centrosymmetric pair so sampling is actually different
	
	/* work space */
	dArrayT lengths(npts);
	dArrayT lx(npts), ly(npts);
	dArrayT du(npts);
	dArrayT fx(npts), fy(npts);

	/* get deformation gradient */
	const dMatrixT& Fcurr = F();

	fx = 0.0;
	fy = 0.0;
	double a = 0.0;
	for (int i = 0; i < npts; i++)
	{
		/* perturbation */
		double dx = eps*cos(a);
		double dy = eps*sin(a);
		
		/* compute deformed "bonds" */
		double aj = 0.0;
		double L0 = 1.0;
		for (int j = 0; j < npts; j++)
		{
			double Xj = L0*cos(aj);
			double Yj = L0*sin(aj);
		
			double lxj = Fcurr(0,0)*Xj + Fcurr(0,1)*Yj - dx;
			double lyj = Fcurr(1,0)*Xj + Fcurr(1,1)*Yj - dy;
			
			/* store lengths */
			lx[j] = lxj;
			ly[j] = lyj;
			lengths[j] = sqrt(lxj*lxj + lyj*lyj);
			
			/* next bond */
			aj += da;
		}
		
		/* derivatives of the potential */
		fPotential->MapDFunction(lengths, du);

		/* integrate (sum) force */
		double& fxi = fx[i];
		double& fyi = fy[i];
		for (int k = 0; k < npts; k++)
		{
			double Dubyl = du[k]/lengths[k];
		
			fxi += Dubyl*lx[k];
			fyi += Dubyl*ly[k];
		}
		
		/* energy release */
		dU[i] = fxi*dx + fyi*dy;
	
		/* next perturbation direction */
		a += da;
	}
}

/* describe the parameters needed by the interface */
void VIB2D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	NL_E_MatT::DefineParameters(list);
	VIB_E_MatT::DefineParameters(list);
}

/* information about subordinate parameter lists */
void VIB2D::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	NL_E_MatT::DefineSubs(sub_list);
	VIB_E_MatT::DefineSubs(sub_list);

	/* choice of integration schemes */
	sub_list.AddSub("circle_integration_choice", ParameterListT::Once, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* VIB2D::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* sub = NL_E_MatT::NewSub(name);
	if (sub) 
		return sub;
	if (name == "circle_integration_choice")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		/* evenly spaced */
		ParameterContainerT even("even_spaced");
		ParameterT n_even(ParameterT::Integer, "n_even");
		n_even.AddLimit(1, LimitT::LowerInclusive);
		even.AddParameter(n_even);
		choice->AddSub(even);

		/* Gauss integration */
		ParameterContainerT gauss("Gauss_points");
		ParameterT n_gauss(ParameterT::Integer, "n_gauss");
		n_gauss.AddLimit(9, LimitT::Only);
		n_gauss.AddLimit(10, LimitT::Only);
		gauss.AddParameter(n_gauss);
		choice->AddSub(gauss);
	
		return choice;
	}
	else
		return VIB_E_MatT::NewSub(name);
}

/* describe the parameters needed by the interface */
void VIB2D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	NL_E_MatT::TakeParameterList(list);
	VIB_E_MatT::TakeParameterList(list);

	/* construct integration scheme */
	const ParameterListT& points = list.GetListChoice(*this, "circle_integration_choice");
	if (points.Name() == "even_spaced")
	{
		int n = points.GetParameter("n_even");
		fCircle = new EvenSpacePtsT(n);
	}
	else if (points.Name() == "Gauss_points")
	{
		int n = points.GetParameter("n_gauss");
		fCircle = new GaussPtsT(n);
	}
	else
		ExceptionT::GeneralFail("VIB2D::TakeParameterList", "unrecognized point scheme \"%s\"", points.Name().Pointer());	

	/* default construction */
	SetAngle(0.0);

	/* set the unstressed energy density */
	SetReferenceEnergy();
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* compute the symetric Cij reduced index matrix */
void VIB2D::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
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
	double& C22 = moduli(1,1);
	double& C66 = moduli(2,2);

	double& C26 = moduli(1,2);
	double& C16 = moduli(0,2);
	double& C12 = moduli(0,1);

	/* initialize kernel pointers */
	double* pdU  = fdU.Pointer();
	double* pddU = fddU.Pointer();
	double* pl   = fLengths.Pointer();
	double* pj   = fjacobian.Pointer();

	/* set pointers */
	double *p11, *p22, *p26, *p16, *p12;	
	SetModuliPointers2D(p11, p22, p26, p16, p12);
	
	int N = fLengths.Length();
	for (int i = 0; i < N; i++)
	{
		double factor = (*pj++)*((*pddU++)/((*pl)*(*pl)) -
		                         (*pdU++)/((*pl)*(*pl)*(*pl)));
		pl++;

		C11 += factor*(*p11++);
		C22 += factor*(*p22++);

		C26 += factor*(*p26++);
		C16 += factor*(*p16++);
		C12 += factor*(*p12++);
	}

	/* Cauchy symmetry */
	C66 = C12;
	
	/* make symmetric */
	moduli.CopySymmetric();
}

/*
* Compute the symetric 2nd Piola-Kirchhoff reduced index vector.
*/
void VIB2D::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
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
	double& s12 = PK2[2];

	/* initialize kernel pointers */
	double* pdU = fdU.Pointer();
	double* pl  = fLengths.Pointer();
	double* pj  = fjacobian.Pointer();
	
	/* set pointers */
	double *p11, *p22, *p12;	
	SetStressPointers2D(p11, p22, p12);
	
	int N = fLengths.Length();
	for (int i = 0; i < N; i++)
	{
		double factor = (*pj++)*(*pdU++)/(*pl++);

		s11 += factor*(*p11++);
		s22 += factor*(*p22++);
		s12 += factor*(*p12++);
	}
}

/*
* Returns the strain energy density for the specified strain.
*/
double VIB2D::ComputeEnergyDensity(const dSymMatrixT& E)
{
	return( VIBEnergyDensity(E) );
}
