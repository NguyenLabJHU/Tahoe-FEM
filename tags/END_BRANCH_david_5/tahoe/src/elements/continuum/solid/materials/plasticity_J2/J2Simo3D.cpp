/* $Id: J2Simo3D.cpp,v 1.15 2004-07-15 08:28:54 paklein Exp $ */
/* created: paklein (06/22/1997) */
#include "J2Simo3D.h"
#include "ElementCardT.h"
#include "StringT.h"
using namespace Tahoe;

/* constants */
const double sqrt23 = sqrt(2.0/3.0);

/* constructor */
J2Simo3D::J2Simo3D(void):
	ParameterInterfaceT("Simo_J2")
{

}

/* form of tangent matrix (symmetric by default) */
GlobalT::SystemTypeT J2Simo3D::TangentType(void) const
{
	return GlobalT::kNonSymmetric;
}

/* update internal variables */
void J2Simo3D::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Update(element, Mu(), NumIP());
}

/* reset internal variables to last converged solution */
void J2Simo3D::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Reset(element, NumIP());
}

/* modulus */
const dMatrixT& J2Simo3D::c_ijkl(void)
{
	/* Compute F_mechanical and f_relative 3D */
	ComputeGradients();

	int ip = CurrIP();
	int nip = NumIP();
	ElementCardT& element = CurrentElement();

	/* compute isochoric elastic stretch */
	double J = fFmech.Det();
	const dSymMatrixT& b_els = TrialElasticState(fFmech, ffrel, element, nip, ip);
	
	/* elastic tangent modulus */
	ComputeModuli(J, b_els, fModulus);
	
	/* elastoplastic correction */
	fModulus += ModuliCorrection(element, Mu(), nip, ip);
	
	return fModulus;
}

/* stress */
const dSymMatrixT& J2Simo3D::s_ij(void)
{
	/* Compute F_total and f_relative 3D */
	ComputeGradients();

	int ip = CurrIP();
	int nip = NumIP();
	ElementCardT& element = CurrentElement();

	/* compute isochoric elastic stretch */
	double J = fFmech.Det();
	const dSymMatrixT& b_els = TrialElasticState(fFmech, ffrel, element, nip, ip);
	
	/* elastic stress */
	ComputeCauchy(J, b_els, fStress);

	/* modify Cauchy stress (return mapping) */
	double mu = Mu();
#pragma message("make elastic its a parameter")
	int iteration = fFSMatSupport->GroupIterationNumber();
	if (iteration > -1 && PlasticLoading(element, mu, ip)) /* 1st iteration is elastic */
//	if (PlasticLoading(element, ip)) /* no iteration is elastic */
	{
		/* element not yet plastic */
		if (!element.IsAllocated())
		{
			/* allocate element storage */
			AllocateElement(element, nip);
		
			/* set trial state and load data */
			TrialElasticState(fFmech, ffrel, element, nip, ip);
			
			/* set the loading state */
			PlasticLoading(element, mu, ip);
		}
	
		/* apply correction due to the return mapping */
		fStress += StressCorrection(element, mu, ip);
	}

	return fStress;
}

/* returns the strain energy density for the specified strain */
double J2Simo3D::StrainEnergyDensity(void)
{
	/* Compute F_total and f_relative */
	ComputeGradients();

	/* compute isochoric elastic stretch */
	double J = fFmech.Det();
	const dSymMatrixT& b_els = TrialElasticState(fFmech, ffrel,
		CurrentElement(), NumIP(), CurrIP());

	return ComputeEnergy(J, b_els);
}

/* incremental heat generation */
double J2Simo3D::IncrementalHeat(void)
{
	/* trust the "current" element is already loaded */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated())
		return fInternal[kHeatIncr];
	else
		return 0.0;
}

/* returns the number of output variables */
int J2Simo3D::NumOutputVariables(void) const { return 4; }

/* returns labels for output variables */
void J2Simo3D::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set labels */
	labels.Dimension(4);
	labels[0] = "alpha";
	labels[1] = "norm_beta";
	labels[2] = "VM_Kirch";
	labels[3] = "press";
}

/* compute output variables */
void J2Simo3D::ComputeOutput(dArrayT& output)
{
	/* check */
	if (output.Length() < 4)
		ExceptionT::GeneralFail("J2Simo3D::ComputeOutput", "expecting 4 output variables: %d", output.Length());

	/* compute Cauchy stress (load state variables) */
	dSymMatrixT stress = s_ij();
	
	/* pressure */
	output[3] = stress.Trace()/3.0;

	/* Cauchy -> relative stress = dev[Kirchhoff] - beta */
	stress *= fFmech.Det();
	stress.Deviatoric();

	/* apply update to state variable data */
	bool has_internal = CurrentElement().IsAllocated();
	if (has_internal)
	{
		/* get flags */
		iArrayT& flags = CurrentElement().IntegerData();
		if (flags[CurrIP()] == kIsPlastic)
		{
			/* factors */
			double alpha = fInternal[kalpha];
			double dgamma = fInternal[kdgamma];
			double mu_bar_bar = fInternal[kmu_bar_bar];
			double k = 2.0*mu_bar_bar*dgamma/Mu();
		
			/* update variables */
			alpha += sqrt23*dgamma;
			fbeta_bar_trial_.SetToCombination(1.0, fbeta_bar, k*dH(alpha)/3.0, fUnitNorm);

			/* write output */
			output[0] = alpha;
			output[1] = sqrt(fbeta_bar_trial_.ScalarProduct());
			stress -= fbeta_bar_trial_;
		}
		else
		{
			/* write output */
			output[0] = fInternal[kalpha];
			output[1] = sqrt(fbeta_bar.ScalarProduct());
			stress -= fbeta_bar;
		}
	}
	else
	{
		output[0] = 0.0;
		output[1] = 0.0;
	}

	/* ||dev[t]|| */
	output[2] = sqrt(stress.ScalarProduct())/sqrt23;
}

/* information about subordinate parameter lists */
void J2Simo3D::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SimoIso3D::DefineSubs(sub_list);
	J2SimoC0HardeningT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* J2Simo3D::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* sub = SimoIso3D::NewSub(name);
	if (sub)
		return sub;
	else
		return J2SimoC0HardeningT::NewSub(name);
}

/* accept parameter list */
void J2Simo3D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SimoIso3D::TakeParameterList(list);
	J2SimoC0HardeningT::TakeParameterList(list);

	/* dimension work space */
	fFmech.Dimension(3);
	ffrel.Dimension(3);
	fF_temp.Dimension(3);
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* compute F_mechanical and f_relative */
void J2Simo3D::ComputeGradients(void)
{
	/* mechanical part of the gradient */
	fFmech = F_mechanical();

	/* relative deformation gradient */
	fF_temp.Inverse(F_mechanical_last());
	ffrel.MultAB(fFmech, fF_temp);
}
