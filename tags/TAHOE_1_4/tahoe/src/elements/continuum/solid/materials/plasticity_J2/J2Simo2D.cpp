/* $Id: J2Simo2D.cpp,v 1.13 2003-11-21 22:46:48 paklein Exp $ */
/* created: paklein (06/22/1997) */
#include "J2Simo2D.h"
#include "StringT.h"
#include "ElementCardT.h"

using namespace Tahoe;

/* constants */
const double sqrt23 = sqrt(2.0/3.0);

/* constructor */
J2Simo2D::J2Simo2D(ifstreamT& in, const FSMatSupportT& support):
	SimoIso2D(in, support),
	J2SimoC0HardeningT(in, NumIP(), Mu()),
	fFmech(3),
	ffrel(3),
	fF_temp(2),
	fFmech_2D(2),
	ffrel_2D(2)	
{

}

/* form of tangent matrix (symmetric by default) */
GlobalT::SystemTypeT J2Simo2D::TangentType(void) const
{
	return GlobalT::kNonSymmetric;
}

/* update internal variables */
void J2Simo2D::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Update(element);
}

/* reset internal variables to last converged solution */
void J2Simo2D::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Reset(element);
}

/* print parameters */
void J2Simo2D::Print(ostream& out) const
{
	/* inherited */
	SimoIso2D::Print(out);
	J2SimoC0HardeningT::Print(out);
}

/* modulus */
const dMatrixT& J2Simo2D::c_ijkl(void)
{
	/* Compute F_mechanical and f_relative 3D */
	ComputeGradients();

	int ip = CurrIP();
	ElementCardT& element = CurrentElement();

	/* compute isochoric elastic stretch */
	double J = fFmech.Det();
	const dSymMatrixT& b_els = TrialElasticState(fFmech, ffrel, element, ip);
	
	/* 3D elastic stress */
	ComputeModuli(J, b_els, fModulus);
	
	/* elastoplastic correction */
	fModulus += ModuliCorrection(element, ip);
	
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(fModulus);
	fModulus2D *= fThickness;

	return fModulus2D;
}

/* stress */
const dSymMatrixT& J2Simo2D::s_ij(void)
{
	/* Compute F_total and f_relative 3D */
	ComputeGradients();

	int ip = CurrIP();
	ElementCardT& element = CurrentElement();

	/* compute isochoric elastic stretch */
	double J = fFmech.Det();
	const dSymMatrixT& b_els = TrialElasticState(fFmech, ffrel, element, ip);
	
	/* 3D elastic stress */
	ComputeCauchy(J, b_els, fStress);

	/* modify Cauchy stress (return mapping) */
	if (PlasticLoading(element, ip))
	{
		/* element not yet plastic */
		if (!element.IsAllocated())
		{
			/* allocate element storage */
			AllocateElement(element);
		
			/* set trial state and load data */
			TrialElasticState(fFmech, ffrel, element, ip);
			
			/* set the loading state */
			PlasticLoading(element, ip);
		}
	
		/* apply correction due to the return mapping */
		fStress += StressCorrection(element, ip);
	}
	
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(fStress);
	fStress2D *= fThickness;

	return fStress2D;
}

/* returns the strain energy density for the specified strain */
double J2Simo2D::StrainEnergyDensity(void)
{
	/* Compute F_total and f_relative 3D */
	ComputeGradients();

	/* compute isochoric elastic stretch */
	double J = fFmech.Det();
	const dSymMatrixT& b_els = TrialElasticState(fFmech, ffrel,
		CurrentElement(), CurrIP());

	return fThickness*ComputeEnergy(J, b_els);
}

/* incremental heat generation */
double J2Simo2D::IncrementalHeat(void)
{
	/* trust the "current" element is already loaded */
	const ElementCardT& element = CurrentElement();
	if (element.IsAllocated())
		return fInternal[kHeatIncr];
	else
		return 0.0;
}

/** returns the number of output variables */
int J2Simo2D::NumOutputVariables(void) const { return 4; }

/** returns labels for output variables */
void J2Simo2D::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set labels */
	labels.Dimension(4);
	labels[0] = "alpha";
	labels[1] = "norm_beta";
	labels[2] = "VM_Kirch";
	labels[3] = "press";
}

/** compute output variables */
void J2Simo2D::ComputeOutput(dArrayT& output)
{
	/* check */
	if (output.Length() < 4) {
		cout << "\n J2Simo2D::ComputeOutput: expecting 4 output variables: " 
		     << output.Length() << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* compute Cauchy stress (load state variables) */
	s_ij();
	dSymMatrixT stress = fStress;

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
			double k = 2.0*mu_bar_bar*dgamma/fmu;
		
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

/***********************************************************************
* Protected
***********************************************************************/

/* print name */
void J2Simo2D::PrintName(ostream& out) const
{
	/* inherited */
	SimoIso2D::PrintName(out);
	J2SimoC0HardeningT::PrintName(out);
}

/***********************************************************************
* Private
***********************************************************************/

/* compute F_mechanical and f_relative */
void J2Simo2D::ComputeGradients(void)
{
	/* compute relative deformation gradient */
	fFmech_2D = F_mechanical();
	fF_temp.Inverse(F_mechanical_last());
	ffrel_2D.MultAB(fFmech_2D, fF_temp);

	/* 2D -> 3D */
	fFmech.Rank2ExpandFrom2D(fFmech_2D);
	ffrel.Rank2ExpandFrom2D(ffrel_2D);
	
	/* approximate effect of plane constraint */
	if (HasThermalStrain())
	{
		double F_inv_last = (F_thermal_inverse_last())(0,0);
		double F_inv = (F_thermal_inverse())(0,0);
		
		/* assuming isotropic thermal strain */
		fFmech(2,2) = F_inv;
		ffrel(2,2) = F_inv/F_inv_last;
	}
	else
	{
		/* plane strain */
		fFmech(2,2) = 1.0;
		ffrel(2,2) = 1.0;
	}
}
