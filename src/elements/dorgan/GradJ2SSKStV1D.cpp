/* $Id: GradJ2SSKStV1D.cpp,v 1.5 2004-08-04 22:02:13 rdorgan Exp $ */
#include "GradJ2SSKStV1D.h"
#include "GradSSMatSupportT.h"
#include "iArrayT.h"
#include "ElementCardT.h"
#include "StringT.h"

#include <math.h>

using namespace Tahoe;

/* parameters */
const int    kNumInternal = 4;
const int    kNSD         = 1;
const double kYieldTol    = 1.0e-10;

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
	GradSSSolidMatT(),  // remove

	fStress(kNSD),
	fElasticStrain(kNSD),

	fK(NULL),

	fNumIP(-1),
	
	fModulus               (dSymMatrixT::NumValues(kNSD)),
	fOffDiagonalModulus_bh (dSymMatrixT::NumValues(kNSD), 1),
	fOffDiagonalModulus_hb (1, dSymMatrixT::NumValues(kNSD)),
	fGradientModulus_hh    (1),
	fGradientModulus_hp    (1,kNSD),
	fGradientModulus_hq    (1)
{

}

/** destructor */
GradJ2SSKStV1D::~GradJ2SSKStV1D(void) { delete fK; }

/* update internal variables */
void GradJ2SSKStV1D::UpdateHistory(void)
{
	ElementCardT& element = CurrentElement();

	/* update plastic variables */
	for (int ip = 0; ip < fNumIP; ip++)
	{
		LoadData(element, ip);

		/* update state */
		fPlasticStrain_0 = fPlasticStrain_j;
		fNorm_0          = fNorm_j;
		fInternal_0      = fInternal_j;
	}
}

/* reset internal variables to last converged solution */
void GradJ2SSKStV1D::ResetHistory(void)
{
	ElementCardT& element = CurrentElement();

	for (int ip = 0; ip < fNumIP; ip++)
	{
		LoadData(element, ip);

		/* reset state */
		fPlasticStrain_j = fPlasticStrain_0;
		fNorm_j          = fNorm_0;
		fInternal_j      = fInternal_0;
	}
}

/* modulus */
const dMatrixT& GradJ2SSKStV1D::c_ijkl(void)
{
	/* elastic modulus */
	fModulus = Young();
	return fModulus;
}

/* off diagonal moduli for Kar */
const dMatrixT& GradJ2SSKStV1D::odm_bh_ij()
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	/* off-diagonal modulus */
	fOffDiagonalModulus_bh = 0.;
	fOffDiagonalModulus_bh.AddScaled(-Young(), fNorm_j);
	return fOffDiagonalModulus_bh;
}

/* off diagonal moduli for K_ra */
const dMatrixT& GradJ2SSKStV1D::odm_hb_ij()
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	/* off-diagonal modulus */
	fOffDiagonalModulus_hb = 0.;
	fOffDiagonalModulus_hb.AddScaled(-Young(), fNorm_j);
	return fOffDiagonalModulus_hb;
}

/* moduli for local term in K_hh */
const dMatrixT& GradJ2SSKStV1D::gm_hh()
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	double fLambda     = GradSSSolidMatT::Lambda();
	double fGradLambda = GradLambda();
	double fLapLambda  = LapLambda();
	
	double fR     = K(fLambda) - K(0.0);
	double fdR    = dK(fLambda);
	double fddR   = ddK(fLambda);
	double fdddR  = ddK(fLambda);
	
	double fGrad1R = Grad1R(fLambda, fGradLambda, fLapLambda);
	double fGrad2R = Grad2R(fLambda, fGradLambda, fLapLambda);
	double fGrad3R = Grad3R(fLambda, fGradLambda, fLapLambda);
	double fGrad4R = Grad4R(fLambda, fGradLambda, fLapLambda);

	/* gradient modulus */
	if (!fInternal_j[kWeakened])
	{
		//		fGradientModulus_hh[0] = Young() + fdR*(1-fk_r*fR-2*fc_r*fk_r*fGrad2R-fc_r*fc_r*fk_r*fGrad4R)-2*fc_r*fk_r*fddR*fGradLambda*(fGrad1R+fc_r*fGrad3R);
		fGradientModulus_hh[0] = Young()+(1-fk_r*fR+fc_r*fk_r*fGrad2R)*(fdR+fc_r*fdddR*pow(fGradLambda,2)+fc_r*fddR*fLapLambda);
		fGradientModulus_hh[0] -= fc_r*fk_r*(2*(fGrad1R+fc_r*fGrad3R)*fddR*fGradLambda+(fGrad2R+fc_r*fGrad4R)*fdR);
	}
	else
	{
		//		fGradientModulus_hh[0] = Young() + fdR*(-2*fc_r*fk_r*fGrad2R-fc_r*fc_r*fk_r*fGrad4R)-2*fc_r*fk_r*fddR*fGradLambda*(fGrad1R+fc_r*fGrad3R);
		fGradientModulus_hh[0] = Young()+(1-fk_r*fR+fc_r*fk_r*fGrad2R)*(fc_r*fdddR*pow(fGradLambda,2)+fc_r*fddR*fLapLambda);
		fGradientModulus_hh[0] += (fc_r*fk_r*fGrad2R)*(fdR);
		fGradientModulus_hh[0] -= fc_r*fk_r*(2*(fGrad1R+fc_r*fGrad3R)*fddR*fGradLambda+(fGrad2R+fc_r*fGrad4R)*fdR);
	}

	return fGradientModulus_hh;
}

/* moduli for gradient term in K_hp */
const dMatrixT& GradJ2SSKStV1D::gm_hp()
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	double fLambda     = GradSSSolidMatT::Lambda();
	double fGradLambda = GradLambda();
	double fLapLambda  = LapLambda();
	
	double fR     = K(fLambda) - K(0.0);
	double fdR    = dK(fLambda);
	double fddR   = ddK(fLambda);
	
	double fGrad1R = Grad1R(fLambda, fGradLambda, fLapLambda);
	double fGrad2R = Grad2R(fLambda, fGradLambda, fLapLambda);
	double fGrad3R = Grad3R(fLambda, fGradLambda, fLapLambda);

	/* gradient modulus */
	fGradientModulus_hp[0] = 2*fc_r*fddR*fGradLambda*(1-fk_r*fR-fk_r*fc_r*fGrad2R)-2*fk_r*fc_r*fdR*(fGrad1R+fc_r*fGrad3R);
	return fGradientModulus_hp;
}

/* moduli for gradient term in K_hq */
const dMatrixT& GradJ2SSKStV1D::gm_hq()
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	double fLambda     = GradSSSolidMatT::Lambda();
	double fGradLambda = GradLambda();
	double fLapLambda  = LapLambda();
	
	double fR     = K(fLambda) - K(0.0);
	double fdR    = dK(fLambda);
	
	double fGrad2R = Grad2R(fLambda, fGradLambda, fLapLambda);

	/* gradient modulus */
	fGradientModulus_hq = fc_r*fdR*(1-fk_r*fR-fk_r*fc_r*fGrad2R);
	return fGradientModulus_hq;
}

/* stress */
const dSymMatrixT& GradJ2SSKStV1D::s_ij(void)
{
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();

	/* load and compute internal variables */
	LoadData(element, ip);

	double& ftrial = fInternal_j[kYieldCrt];
	double& fYieldStep = fInternal_j[kYieldStep];

	/* compute trial elastic strain */
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_trl = ElasticStrain(e_tot, fPlasticStrain_0);

	/* compute trial elastic stress */
	HookeanStress(e_trl, fStress);

	/* compute trial yield condition */
	ftrial = YieldCondition(Lambda_last(), LapLambda_last());

	int iteration = fGradSSMatSupport->GroupIterationNumber();

	if (iteration <= -1)     // elastic iteration
	{
		fNorm_j = 0.;        // elastic
		ftrial = 0.;
		fYieldStep = fInternal_0[kYieldStep];
	}
	else if (ftrial > kYieldTol || fYieldStep)   // check elastic predictor
	{
		if (fStress[0] >= 0.)
			fNorm_j = 1.0;   // plastic
		else
			fNorm_j = -1.0;  // plastic
		fYieldStep = 1.;
	}
	else
	{
		fNorm_j = 0.;    // elastic
		ftrial = 0.;
	}

	/* update plastic strain */
	fPlasticStrain_j = fPlasticStrain_0;
	fPlasticStrain_j.AddScaled(del_Lambda()*fYieldStep, fNorm_j);

	/* plastic corrector */
	fStress.AddScaled(-Young()*del_Lambda()*fYieldStep, fNorm_j);

	/* update yield condition */
	if (fYieldStep && iteration > -1) ftrial = YieldCondition(GradSSSolidMatT::Lambda(), LapLambda()); // plastic

	fInternal_j[kStrain] =  e_tot[0];

	return fStress;
}

/* yield criteria moduli */
double GradJ2SSKStV1D::yc()
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	/* yield condition */
	return fInternal_j[kYieldCrt];
}

/* returns the strain energy density for the specified strain */
double GradJ2SSKStV1D::StrainEnergyDensity(void)
{
	/* load internal variables */
	LoadData(CurrentElement(), CurrIP());

	/* compute elastic strain */
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_els = ElasticStrain(e_tot, fPlasticStrain_j);

	return HookeanEnergy(e_els);		
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
	/* stress tensor (loads element data and sets fStress) */
	s_ij();

	/* total strain */
	output[0] = fInternal_j[kStrain];

	/* equivalent plastic strain */
	output[1] = fPlasticStrain_j[0];

	/* isotropic hardening */
	output[2] = GradSSSolidMatT::Lambda();

	/* gradient isotropic hardening */
	output[3] = GradLambda();

	/* Laplacian isotropic hardening */
	output[4] = LapLambda();

	/* isotropic hardening conjugate force */
	output[5] = K(output[2]) - K(0.0);

	/* gradient isotropic hardening conjugate force */
	output[6] = Grad1R(GradSSSolidMatT::Lambda(),GradLambda(),LapLambda());

	/* Laplacian isotropic hardening conjugate force */
	output[7] = Grad2R(GradSSSolidMatT::Lambda(),GradLambda(),LapLambda());

	/* elastic strain */
	//	output[8] = fInternal_j[kStrain] - fPlasticStrain_j[0];

	/* yield stress */
	//	output[2] = K(output[2]);

	/* gradient enhanced yield stress */
	//	output[3] = K(output[2]) + fc_r*dK(output[2])*LapLambda();
}

/* implementation of the ParameterInterfaceT interface */
void GradJ2SSKStV1D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	IsotropicT::DefineParameters(list);
	GradSSSolidMatT::DefineParameters(list);
	ParameterInterfaceT::DefineParameters(list);

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
	GradSSSolidMatT::DefineSubs(sub_list);
	ParameterInterfaceT::DefineSubs(sub_list);

	/* hardening function */
	sub_list.AddSub("hardening_function_choice", ParameterListT::Once, true);
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

	if (name == "hardening_function_choice")
	{
		order = ParameterListT::Choice;
	
		/* function types */
		sub_lists.AddSub("linear_function");
		sub_lists.AddSub("cubic_spline");
		sub_lists.AddSub("linear_exponential");
		sub_lists.AddSub("power_law");
		sub_lists.AddSub("piecewise_linear");
	}
	else /* inherited */
		ParameterInterfaceT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* GradJ2SSKStV1D::NewSub(const StringT& name) const
{
	ParameterInterfaceT* sub = NULL;

	/* try each base class */
	sub = IsotropicT::NewSub(name);
	if (sub) return sub;
	
	sub = HookeanMatT::NewSub(name);
	if (sub) return sub;
	
	sub = GradSSSolidMatT::NewSub(name);
	if (sub) return sub;

	/* try to construct C1 function */
	C1FunctionT* function = C1FunctionT::New(name);
	if (function)
		return function;
	else /* inherited */
		return ParameterInterfaceT::NewSub(name);
}

/* accept parameter list */
void GradJ2SSKStV1D::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "GradJ2SSKStV1D::TakeParameterList";

	fNumIP = NumIP();

	/* inherited */
	IsotropicT::TakeParameterList(list);
	HookeanMatT::TakeParameterList(list);
	GradSSSolidMatT::TakeParameterList(list);
	
	/* length scale in nonlocal measure of isotropic hardening */
	fc_r = list.GetParameter("isotropic_hardening_length_scale");
	
	/* coefficient in nonassociative plasticity for isotropic hardening  */
	fk_r = list.GetParameter("isotropic_hardening_nonassociative");

	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* construct hardening function */
	const ParameterListT& hardening = list.GetListChoice(*this, "hardening_function_choice");
	fK = C1FunctionT::New(hardening.Name());
	if (!fK) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", hardening.Name().Pointer());
	fK->TakeParameterList(hardening);

	/* allocate space for all elements */
	AllocateAllElements();
}

/*************************************************************************
* Private
*************************************************************************/

/* incremental change in field */
double GradJ2SSKStV1D::del_Lambda(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	/* increment in field */
	return GradSSSolidMatT::Lambda() - Lambda_last();
}

/* incremental change in laplacian of field */
double GradJ2SSKStV1D::del_GradLambda(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	/* increment in field */
	return GradLambda() - GradLambda_last();
}

/* incremental change in laplacian of field */
double GradJ2SSKStV1D::del_LapLambda(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	/* increment in field */
	return LapLambda() - LapLambda_last();
}

/* set modulus */
void GradJ2SSKStV1D::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli1D(modulus);
}

/* return a pointer to a new element object constructed with
* the data from element */
void GradJ2SSKStV1D::AllocateAllElements(void)
{
	/* determine storage */
	int i_size = 0;
	// remove fFlags!
	i_size += fNumIP; //fFlags

	int d_size = 0;
	int dim = dSymMatrixT::NumValues(kNSD);
	d_size += dim;	  //fPlasticStrain_0
	d_size += dim;	  //fPlasticStrain_j
	d_size += dim;	  //fNorm_0
	d_size += dim;	  //fNorm_n
	d_size += kNumInternal; //fInternal_0
	d_size += kNumInternal; //fInternal_j

	d_size *= fNumIP;

	/* allocate space for all elements */
	for (int el = 0; el < NumElements(); el++)
	{
		/* get pointer to element el */
		ElementCardT& element = ElementCard(el);

		/* construct new element */
		element.Dimension(fNumIP, d_size);
	
		/* initialize values */
		element.IntegerData() = 1;  // REMOVE FLAG INTEGER DATA
		element.DoubleData()  = 0.0;
	}
}

/* load element data for the specified integration point */
void GradJ2SSKStV1D::LoadData(const ElementCardT& element, int ip)
{
	/* fetch arrays */
	const dArrayT& d_array = element.DoubleData();
	
	/* decode */
	dSymMatrixT::DimensionT sdim = dSymMatrixT::int2DimensionT(kNSD);
	int dim   = dSymMatrixT::NumValues(kNSD);
	int block = 4*dim + 2*kNumInternal;
	int dex   = ip*block;

	fPlasticStrain_0.Alias(sdim,         &d_array[dex]                );
	fPlasticStrain_j.Alias(sdim,         &d_array[dex += dim]         );
	fNorm_0.Alias         (sdim,         &d_array[dex += dim]         );
	fNorm_j.Alias         (sdim,         &d_array[dex += dim]         );
	fInternal_0.Alias     (kNumInternal, &d_array[dex += dim]         );
	fInternal_j.Alias     (kNumInternal, &d_array[dex += kNumInternal]);
}

/* returns elastic strain */
const dSymMatrixT& GradJ2SSKStV1D::ElasticStrain(const dSymMatrixT& totalstrain,
	const dSymMatrixT& plasticstrain)
{	
	/* compute elastic strain */
	fElasticStrain.DiffOf(totalstrain, plasticstrain);

	return fElasticStrain;
}	

/** 1st gradient of Isotropic Hardening conjugate force */
double GradJ2SSKStV1D::Grad1R(double lambda, double gradlambda, double laplambda)
{
#pragma unused(laplambda)

	double fdR    = dK(lambda);

	/* increment in field */
	return fdR*gradlambda;
}
	
/** 2nd gradient of Isotropic Hardening conjugate force */
double GradJ2SSKStV1D::Grad2R(double lambda, double gradlambda, double laplambda)
{
	double fdR    = dK(lambda);
	double fddR   = ddK(lambda);

	/* increment in field */
	return fddR*pow(gradlambda,2) + fdR*laplambda;
}
	
/** 3rd gradient of Isotropic Hardening conjugate force */
double GradJ2SSKStV1D::Grad3R(double lambda, double gradlambda, double laplambda)
{
	double fddR   = ddK(lambda);
	double fdddR  = dddK(lambda);

	/* increment in field */
	return fdddR*pow(gradlambda,3) + 3*fddR*gradlambda*laplambda;
}
	
/** 4th gradient of Isotropic Hardening conjugate force */
double GradJ2SSKStV1D::Grad4R(double lambda, double gradlambda, double laplambda)
{
	double fddR     = ddK(lambda);
	double fdddR    = dddK(lambda);
	double fddddR   = ddddK(lambda);

	/* increment in field */
	return fddddR*pow(gradlambda,4) + 6*fdddR*pow(gradlambda,2)*laplambda + 3*fddR*pow(laplambda,2);
}

/* yield criteria moduli */
double GradJ2SSKStV1D::YieldCondition(double lambda, double laplambda)
{
	//	double flambda = GradSSSolidMatT::Lambda();

	if (( K(lambda) < 1e-03 || fInternal_j[kWeakened] ) && 0.)
	{			
		if (!fInternal_0[kWeakened])
		{				
			cout << "Element [" << CurrElementNumber() << "], ip [" << CurrIP() << "]:  ";
			cout << "Isotropic Hardening: " << lambda;
			cout << "       Yield Stress: " << K(lambda) << endl;
		}

		fInternal_j[kWeakened] = 1;
		
		/* check yield strength */
		return  fabs(fStress[0]) - ( 1e-03 + fc_r*dK(lambda)*laplambda);
	}
	else
	{
		fInternal_j[kWeakened] = 0;

		/* check yield strength */
		return  fabs(fStress[0]) - ( K(lambda) + fc_r*dK(lambda)*laplambda);
	}
}
