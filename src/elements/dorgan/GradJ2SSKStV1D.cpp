/* $Id: GradJ2SSKStV1D.cpp,v 1.4 2004-07-27 21:13:57 rdorgan Exp $ */
#include "GradJ2SSKStV1D.h"
#include "GradSSMatSupportT.h"

#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* parameters */
const int kNSD = 1;

/* element output data */
const int kNumOutput = 8;
static const char* Labels[kNumOutput] = {
	"epsilon",     // total strain
	"alpha",       // equivalent plastic strain
	"r",           // isotropic hardening
	"r,x",         // gradient isotropic hardening
	"r,xx",        // Laplacian isotropic hardening
	"R",           // isotropic hardening conjugate force
	"R,x",         // gradient isotropic hardening conjugate force
	"R,xx"};       // Laplacian isotropic hardening conjugate force

/* constructor */
GradJ2SSKStV1D::GradJ2SSKStV1D(void):
	ParameterInterfaceT ("grad_small_strain_StVenant_J2_1D"),
	HookeanMatT(kNSD),
	fStress(kNSD),
	fModulus(dSymMatrixT::NumValues(kNSD)),

	fStress_3D(3),
	
	fOffDiagonalModulus_bh (dSymMatrixT::NumValues(kNSD), 1),
	fOffDiagonalModulus_hb (1, dSymMatrixT::NumValues(kNSD)),
	fGradientModulus_hh    (1),
	fGradientModulus_hp    (1,kNSD),
	fGradientModulus_hq    (1)
{

}

/* update internal variables */
void GradJ2SSKStV1D::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) 
		Update(element, NumIP());
}

/* reset internal variables to last converged solution */
void GradJ2SSKStV1D::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) 
		Reset(element, NumIP());
}

/* modulus */
const dMatrixT& GradJ2SSKStV1D::c_ijkl(void)
{
	/* elastoplastic correction */
	fModulus.SumOf(HookeanMatT::Modulus(), ModuliCorrection(CurrentElement(), Young(), fk_r, fc_r, NumIP(), CurrIP()));	
	return fModulus;
}

/* off diagonal moduli for Kar */
const dMatrixT& GradJ2SSKStV1D::odm_bh_ij()
{
	/* off-diagonal modulus */
	fOffDiagonalModulus_bh = 0.;
	fOffDiagonalModulus_bh.AddScaled( 1.0, ODModuliCorrection_bh(CurrentElement(), Young(), fk_r, fc_r, NumIP(), CurrIP()));
	return fOffDiagonalModulus_bh;
}

/* off diagonal moduli for K_ra */
const dMatrixT& GradJ2SSKStV1D::odm_hb_ij()
{
	/* off-diagonal modulus */
	fOffDiagonalModulus_hb = 0.;
	fOffDiagonalModulus_hb.AddScaled( 1.0, ODModuliCorrection_bh(CurrentElement(), Young(), fk_r, fc_r, NumIP(), CurrIP()));
	return fOffDiagonalModulus_hb;
}

/* moduli for local term in K_hh */
const dMatrixT& GradJ2SSKStV1D::gm_hh()
{
	/* off-diagonal modulus */
	fGradientModulus_hh = 0.;
	fGradientModulus_hh.AddScaled( 1.0, GModuliCorrection_hh(CurrentElement(), Young(), fk_r, fc_r, NumIP(), CurrIP()));
	return fGradientModulus_hh;
}

/* moduli for local term in K_hp */
const dMatrixT& GradJ2SSKStV1D::gm_hp()
{
	/* off-diagonal modulus */
	fGradientModulus_hp = 0.;
	fGradientModulus_hp.AddScaled( 1.0, GModuliCorrection_hp(CurrentElement(), Young(), fk_r, fc_r, NumIP(), CurrIP()));
	return fGradientModulus_hp;			
}

/* moduli for local term in K_hq */
const dMatrixT& GradJ2SSKStV1D::gm_hq()
{
	/* off-diagonal modulus */
	fGradientModulus_hq = 0.;
	fGradientModulus_hq.AddScaled( 1.0, GModuliCorrection_hq(CurrentElement(), Young(), fk_r, fc_r, NumIP(), CurrIP()));
	return fGradientModulus_hq;			
}

/* stress */
const dSymMatrixT& GradJ2SSKStV1D::s_ij(void)
{
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_els = ElasticStrain(e_tot, element, NumIP(), ip);

	/* elastic stress */
	HookeanStress(e_els, fStress);

	/* modify Cauchy stress (return mapping) */
	int iteration = fGradSSMatSupport->GroupIterationNumber();
	if (iteration > -1) /* elastic iteration */
		StressCorrection(fStress, element, Young(), fk_r, fc_r, NumIP(), ip);

	return fStress;	
}

/* yield criteria moduli */
double GradJ2SSKStV1D::yc()
{
	/* yield condition */
	return YCModuli(fStress, CurrentElement(), Young(), fk_r, fc_r, NumIP(), CurrIP());
}

/* returns the strain energy density for the specified strain */
double GradJ2SSKStV1D::StrainEnergyDensity(void)
{
	return HookeanEnergy(ElasticStrain(e(), CurrentElement(), NumIP(), CurrIP()));		
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default. */
int GradJ2SSKStV1D::NumOutputVariables(void) const  { return kNumOutput; }
void GradJ2SSKStV1D::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void GradJ2SSKStV1D::ComputeOutput(dArrayT& output)
{
	output = 0.0;

	/* stress tensor (loads element data and sets fStress) */
	s_ij();

	/* total strain */
	output[0] = e()[0];

	const ElementCardT& element = CurrentElement();
	if (element.IsAllocated())
	{
		/* plastic strain */
		output[1] = fInternal[kalpha];
		
		/* status flags */
		const iArrayT& flags = element.IntegerData();
		if (flags[CurrIP()] == kIsPlastic) // output with update
			output[1] += fInternal[kdgamma];
	}
	else
		output[1] = 0.0;

	/* isotropic hardening */
	output[2] = GradSSSolidMatT::Lambda();

	/* gradient isotropic hardening */
	output[3] = GradLambda();

	/* Laplacian isotropic hardening */
	output[4] = LapLambda();

	/* isotropic hardening conjugate force */
	output[5] = K(output[2]) - K(0.0);

	/* gradient isotropic hardening conjugate force */
	output[6] = Grad1R(output[2], output[3], output[4]);

	/* Laplacian isotropic hardening conjugate force */
	output[7] = Grad2R(output[2], output[3], output[4]);
}

/* implementation of the ParameterInterfaceT interface */
void GradJ2SSKStV1D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	IsotropicT::DefineParameters(list);
	HookeanMatT::DefineParameters(list);
	GradJ2SSC0Hardening1DT::DefineParameters(list);

	ParameterT isotropic_hardening_length_scale(fc_r, "isotropic_hardening_length_scale");
	isotropic_hardening_length_scale.SetDefault(0.0);
	isotropic_hardening_length_scale.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(isotropic_hardening_length_scale);

	ParameterT isotropic_hardening_nonassociative(fk_r, "isotropic_hardening_nonassociative");
	isotropic_hardening_nonassociative.SetDefault(0.0);
	isotropic_hardening_nonassociative.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(isotropic_hardening_nonassociative);
}

/* information about subordinate parameter lists */
void GradJ2SSKStV1D::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	IsotropicT::DefineSubs(sub_list);
	HookeanMatT::DefineSubs(sub_list);
	GradJ2SSC0Hardening1DT::DefineSubs(sub_list);
}

/* return the description of the given inline subordinate parameter list */
void GradJ2SSKStV1D::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	/* inherited */
	if (sub_lists.Length() == 0)
		IsotropicT::DefineInlineSub(name, order, sub_lists);
	if (sub_lists.Length() == 0)
		HookeanMatT::DefineInlineSub(name, order, sub_lists);
	if (sub_lists.Length() == 0)
		GradJ2SSC0Hardening1DT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* GradJ2SSKStV1D::NewSub(const StringT& name) const
{
	ParameterInterfaceT* sub = NULL;

	sub = IsotropicT::NewSub(name);
	if (sub) return sub;

	sub = HookeanMatT::NewSub(name);
	if (sub) return sub;
	
	return GradJ2SSC0Hardening1DT::NewSub(name);
}

/* accept parameter list */
void GradJ2SSKStV1D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	IsotropicT::TakeParameterList(list);
	HookeanMatT::TakeParameterList(list);
	GradJ2SSC0Hardening1DT::TakeParameterList(list);
	
	/* length scale in nonlocal measure of isotropic hardening */
	fc_r = list.GetParameter("isotropic_hardening_length_scale");
	
	/* coefficient in nonassociative plasticity for isotropic hardening  */
	fk_r = list.GetParameter("isotropic_hardening_nonassociative");
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* set modulus */
void GradJ2SSKStV1D::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli1D(modulus);
}
