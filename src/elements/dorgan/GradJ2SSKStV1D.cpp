/* $Id: GradJ2SSKStV1D.cpp,v 1.9 2004-11-30 23:06:24 rdorgan Exp $ */
#include "GradJ2SSKStV1D.h"
#include "GradSSMatSupportT.h"
#include "ElementCardT.h"
#include "StringT.h"

#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "iArrayT.h"

#include <math.h>

using namespace Tahoe;

/* parameters */
const int    kNumInternal = 7;
const int    kNSD         = 1;
const double kYieldTol    = 1.0e-10;
const double kYieldStressSmall = 1.0e-03;

/* element output data */
const int kNumOutput = 7;
static const char* Labels[kNumOutput] = {
	"epsilon",     // equivalent total strain
	"alpha",       // equivalent plastic strain
	"r",           // isotropic hardening
	"r,xx",        // Laplacian isotropic hardening
	"R",           // isotropic hardening conjugate force
	"R,xx",        // Laplacian isotropic hardening conjugate force
	"del_r"};      // increment isotropic hardening

/* constructor */
GradJ2SSKStV1D::GradJ2SSKStV1D(void):
	ParameterInterfaceT ("grad_small_strain_StVenant_J2_1D"),
	HookeanMatT(kNSD),

	fStress(kNSD),
	fElasticStrain(kNSD),

	fK(NULL),

	fNumIP(-1),
	
	fModulus               (dSymMatrixT::NumValues(kNSD)),
	fOffDiagonalModulus_bh (dSymMatrixT::NumValues(kNSD),1),
	fOffDiagonalModulus_hb (1,dSymMatrixT::NumValues(kNSD)),
	fGradientModulus_hh    (1),
	fGradientModulus_hp    (1,kNSD),
	fGradientModulus_hq    (1),

	fTensorTemp1(dSymMatrixT::NumValues(kNSD),1),
	fTensorTemp2(kNSD,1)
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
		fUnitNorm_0      = fUnitNorm_j;
		fInternal_0      = fInternal_j;
		fGradIsoHard_0   = fGradIsoHard_j;
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
		fUnitNorm_j      = fUnitNorm_0;
		fInternal_j      = fInternal_0;
		fGradIsoHard_j   = fGradIsoHard_0;
	}
}

/* update flag describing a weakened ip */
void GradJ2SSKStV1D::UpdateWeakened(const ElementCardT& element, int ip)
{
	LoadData(element, ip);

	/* update state */
	fInternal_j[kWeakened] = 1;
}

/* reset flag describing a weakened ip */
void GradJ2SSKStV1D::ResetWeakened(const ElementCardT& element, int ip)
{
	LoadData(element, ip);

	/* update state */
	fInternal_j[kWeakened] = 0;
}

/* modulus */
const dMatrixT& GradJ2SSKStV1D::c_ijkl(void)
{
	/* elastic modulus */
	fModulus = HookeanMatT::Modulus();
	return fModulus;
}

/* off diagonal moduli for Kar */
const dMatrixT& GradJ2SSKStV1D::odm_bh_ij(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	/* off-diagonal modulus */
	fOffDiagonalModulus_bh = 0.;
	fTensorTemp1.Outer(HookeanMatT::Modulus(), fUnitNorm_j);
	fOffDiagonalModulus_bh.AddScaled(-1.0, fTensorTemp1);
	return fOffDiagonalModulus_bh;
}

/* off diagonal moduli for K_ra */
const dMatrixT& GradJ2SSKStV1D::odm_hb_ij(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	/* off-diagonal modulus */
	fOffDiagonalModulus_hb = 0.;
	fTensorTemp1.Outer(HookeanMatT::Modulus(), fUnitNorm_j);
	fOffDiagonalModulus_hb.AddScaled(-1.0, fTensorTemp1);
	fOffDiagonalModulus_hb.Transpose(fOffDiagonalModulus_hb);
	return fOffDiagonalModulus_hb;
}

/* moduli for local term in K_hh */
const dMatrixT& GradJ2SSKStV1D::gm_hh(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	double& fIsoHard     = fInternal_j[kIsoHard];
	double& fLapIsoHard  = fInternal_j[kLapIsoHard];
	double& fYieldStep   = fInternal_j[kYieldStep];

	double fR     = K(fIsoHard) - K(0.0);
	double fdR    = dK(fIsoHard);
	double fddR   = ddK(fIsoHard);
	double fdddR  = dddK(fIsoHard);
	
	dMatrixT fGrad1R = Grad1R(fIsoHard, fGradIsoHard_j, fLapIsoHard);
	double   fGrad2R = Grad2R(fIsoHard, fGradIsoHard_j, fLapIsoHard);
	dMatrixT fGrad3R = Grad3R(fIsoHard, fGradIsoHard_j, fLapIsoHard);
	double   fGrad4R = Grad4R(fIsoHard, fGradIsoHard_j, fLapIsoHard);

	fGradientModulus_hh = 0.;

	if (fYieldStep > 0.5)
	{
		dMatrixT temp(1);
		
		/* gradient modulus */
 		fTensorTemp1.Outer(HookeanMatT::Modulus(), fUnitNorm_j);
		fTensorTemp1.Multx(fUnitNorm_j, temp);
		fGradientModulus_hh[0] = temp[0] + (1-fk_r*fR+fc_r*fk_r*fGrad2R) * (fdR+fc_r*fdddR*fGradIsoHard_j.ScalarProduct()+fc_r*fddR*fLapIsoHard) - fc_r*fk_r*(2*(dMatrixT::Dot(fGrad1R,fGradIsoHard_j)+fc_r*dMatrixT::Dot(fGrad3R,fGradIsoHard_j))*fddR + (fGrad2R+fc_r*fGrad4R)*fdR);

		if (fInternal_0[kWeakened] > 0.5 || fInternal_j[kWeakened] > 0.5)
			fGradientModulus_hh[0] -= (1-fk_r*fR)*fdR;
	}
	else
		fGradientModulus_hh[0] = Young();
//		fGradientModulus_hh[0] = pow(Young(),4);
		
	return fGradientModulus_hh;
}

/* moduli for gradient term in K_hp */
const dMatrixT& GradJ2SSKStV1D::gm_hp(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	double& fIsoHard     = fInternal_j[kIsoHard];
	double& fLapIsoHard  = fInternal_j[kLapIsoHard];

	double fR     = K(fIsoHard) - K(0.0);
	double fdR    = dK(fIsoHard);
	double fddR   = ddK(fIsoHard);
	
	dMatrixT fGrad1R = Grad1R(fIsoHard, fGradIsoHard_j, fLapIsoHard);
	double   fGrad2R = Grad2R(fIsoHard, fGradIsoHard_j, fLapIsoHard);
	dMatrixT fGrad3R = Grad3R(fIsoHard, fGradIsoHard_j, fLapIsoHard);

	/* gradient modulus */
	fGradientModulus_hp.SetToScaled(2*fc_r*fddR*(1-fk_r*fR-fk_r*fc_r*fGrad2R), fGradIsoHard_j);
	fGradientModulus_hp.AddScaled(-2*fk_r*fc_r*fdR,fGrad1R);
	fGradientModulus_hp.AddScaled(-2*fk_r*fc_r*fc_r*fdR,fGrad3R);
	return fGradientModulus_hp;
}

/* moduli for gradient term in K_hq */
const dMatrixT& GradJ2SSKStV1D::gm_hq(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	double& fIsoHard     = fInternal_j[kIsoHard];
	double& fLapIsoHard  = fInternal_j[kLapIsoHard];

	double fR     = K(fIsoHard) - K(0.0);
	double fdR    = dK(fIsoHard);
	
	double fGrad2R = Grad2R(fIsoHard, fGradIsoHard_j, fLapIsoHard);

	/* gradient modulus */
	fGradientModulus_hq = fc_r*fdR*(1-fk_r*fR-fk_r*fc_r*fGrad2R);
	return fGradientModulus_hq;
}

/* stress */
const dSymMatrixT& GradJ2SSKStV1D::s_ij(void)
{
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();

	/* load internal variables */
	LoadData(element, ip);

	double& ftrial       = fInternal_j[kYieldCrt];
	double& fYieldStep   = fInternal_j[kYieldStep];
	double& fIsoHard     = fInternal_j[kIsoHard];
	double& fLapIsoHard  = fInternal_j[kLapIsoHard];
	
	/* compute trial state */
	fIsoHard         = fInternal_0[kIsoHard];
	fGradIsoHard_j   = fGradIsoHard_0;
	fLapIsoHard      = fInternal_0[kLapIsoHard];
	fPlasticStrain_j = fPlasticStrain_0;

	/* increment isotropic hardening */
	fIsoHard       += del_Lambda();
	fGradIsoHard_j += del_GradLambda();
	fLapIsoHard    += del_LapLambda();

	/* compute trial elastic strain */
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_trl = ElasticStrain(e_tot, element, ip);

	/* compute trial elastic stress */
	HookeanStress(e_trl, fStress);

	/* compute trial yield condition */
	ftrial = YieldCondition(fIsoHard, fGradIsoHard_j, fLapIsoHard);

	int iteration = fGradSSMatSupport->GroupIterationNumber();

	/* elastic iteration */
	if (iteration <= -1)
	{
		fUnitNorm_j = 0.;
		ftrial = 0.;
		fYieldStep = 0;
	}
	/* plastic iteration */
	else if (ftrial > kYieldTol || fYieldStep > 0.5)   // check elastic predictor
	{
		fYieldStep = 1.;

		/* compute unit normal */
		fUnitNorm_j.SetToScaled(1.0/sqrt(fStress.ScalarProduct()), fStress);
		
		/* update plastic strain */
		fPlasticStrain_j.AddScaled(del_Lambda(), fUnitNorm_j);

		/* plastic corrector */
		fTensorTemp1.Outer(HookeanMatT::Modulus(), fUnitNorm_j);
		fStress.AddScaled(-del_Lambda(), fTensorTemp1);

 		/* update yield condition */
 		ftrial = YieldCondition(fIsoHard, fGradIsoHard_j, fLapIsoHard);
	}
	/* elastic iteration */
	else
	{
		fUnitNorm_j = 0.;
		ftrial = 0.;
	}

	return fStress;
}

/* unit normal */
const dSymMatrixT& GradJ2SSKStV1D::n_ij(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	/* yield condition */
	return fUnitNorm_j;
}

/* yield criteria moduli */
double GradJ2SSKStV1D::yc(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	/* yield condition */
	return fInternal_j[kYieldCrt];
}

/* yield strenth */
double GradJ2SSKStV1D::ys(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	/* yield condition */
	return fInternal_j[kYieldStrength];
}

/** returns 1 if the ip has weakened during the current time step, 0 otherwise */
int GradJ2SSKStV1D::weakened(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	if (fInternal_j[kYieldStrength] < kYieldStressSmall && fInternal_0[kWeakened] < 0.5)
		return 1;
	else
		return 0;
}

/* returns the strain energy density for the specified strain */
double GradJ2SSKStV1D::StrainEnergyDensity(void)
{
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();

	/* load internal variables */
	LoadData(element, ip);

	/* compute elastic strain */
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_els = ElasticStrain(e_tot, element, ip);

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
	output[0] = sqrt(e().ScalarProduct());

	/* equivalent plastic strain */
	output[1] = sqrt(fPlasticStrain_j.ScalarProduct());

	/* isotropic hardening */
	output[2] = fInternal_j[kIsoHard];

	/* Laplacian isotropic hardening */
	output[3] = fInternal_j[kLapIsoHard];

	/* isotropic hardening conjugate force */
	output[4] = K(output[2]) - K(0.0);

	/* Laplacian isotropic hardening conjugate force */
	output[5] = Grad2R(output[2],fGradIsoHard_j,output[3]);

	/* increment isotropic hardening */
	output[6] = fInternal_j[kIsoHard] - fInternal_0[kIsoHard];
}

/* implementation of the ParameterInterfaceT interface */
void GradJ2SSKStV1D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	GradSSSolidMatT::DefineParameters(list);
	IsotropicT::DefineParameters(list);
	HookeanMatT::DefineParameters(list);

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
	GradSSSolidMatT::DefineSubs(sub_list);
	IsotropicT::DefineSubs(sub_list);
	HookeanMatT::DefineSubs(sub_list);

	/* hardening function */
	sub_list.AddSub("hardening_function_choice", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void GradJ2SSKStV1D::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	/* inherited */
	GradSSSolidMatT::DefineInlineSub(name, order, sub_lists);
	IsotropicT::DefineInlineSub(name, order, sub_lists);
	HookeanMatT::DefineInlineSub(name, order, sub_lists);

	if (name == "hardening_function_choice")
	{
		order = ParameterListT::Choice;
	
		/* function types */
		sub_lists.AddSub("linear_function");
		sub_lists.AddSub("cubic_spline");
		sub_lists.AddSub("linear_exponential");
		sub_lists.AddSub("power_law_2");
		sub_lists.AddSub("piecewise_linear");
	}
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* GradJ2SSKStV1D::NewSub(const StringT& name) const
{
	ParameterInterfaceT* sub = NULL;

	/* try each base class */
	sub = GradSSSolidMatT::NewSub(name);
	if (sub) return sub;

	sub = IsotropicT::NewSub(name);
	if (sub) return sub;
	
	sub = HookeanMatT::NewSub(name);
	if (sub) return sub;
	
	/* try to construct C1 function */
	return C1FunctionT::New(name);
}

/* accept parameter list */
void GradJ2SSKStV1D::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "GradJ2SSKStV1D::TakeParameterList";

	fNumIP = NumIP();

	/* inherited */
	GradSSSolidMatT::TakeParameterList(list);
	IsotropicT::TakeParameterList(list);
	HookeanMatT::TakeParameterList(list);
	
	/* construct hardening function */
	const ParameterListT& hardening = list.GetListChoice(*this, "hardening_function_choice");
	fK = C1FunctionT::New(hardening.Name());
	if (!fK) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", hardening.Name().Pointer());
	fK->TakeParameterList(hardening);

	/* length scale in nonlocal measure of isotropic hardening */
	fc_r = list.GetParameter("isotropic_hardening_length_scale");
	
	/* coefficient in nonassociative plasticity for isotropic hardening  */
	fk_r = list.GetParameter("isotropic_hardening_nonassociative");

	/* allocate space for all elements */
	AllocateAllElements();
}

/*************************************************************************
* Protected
*************************************************************************/

/* returns elastic strain */
const dSymMatrixT& GradJ2SSKStV1D::ElasticStrain(const dSymMatrixT& totalstrain,
	const ElementCardT& element, int ip)
{	
	/* load internal variables */
	LoadData(element, ip);

	/* compute elastic strain */
	fElasticStrain.DiffOf(totalstrain, fPlasticStrain_j);
	
	return fElasticStrain;
}

/*************************************************************************
* Private
*************************************************************************/

/* incremental change in field */
double GradJ2SSKStV1D::del_Lambda(void)
{
	/* increment in field */
	return GradSSSolidMatT::Lambda()[0] - Lambda_last()[0];
}

/* incremental change in gradient of field */
dMatrixT GradJ2SSKStV1D::del_GradLambda(void)
{
	/* increment in field */
	fTensorTemp2 = GradLambda();
	fTensorTemp2.AddScaled(-1.0, GradLambda_last());
	return fTensorTemp2;
}

/* incremental change in laplacian of field */
double GradJ2SSKStV1D::del_LapLambda(void)
{
	/* increment in field */
	return LapLambda()[0] - LapLambda_last()[0];
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
	d_size += dim;	        //fPlasticStrain_0
	d_size += dim;	        //fPlasticStrain_j
	d_size += dim;	        //fUnitNorm_0
	d_size += dim;          //fUnitNorm_n
	d_size += kNumInternal; //fInternal_0
	d_size += kNumInternal; //fInternal_j
	d_size += kNSD*1;       //fGradIsoHard_0
	d_size += kNSD*1;       //fGradIsoHard_j

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
	int block = 4*dim + 2*kNumInternal + 2*kNSD*1;
	int dex   = ip*block;

	fPlasticStrain_0.Alias(sdim,         &d_array[dex]                );
	fPlasticStrain_j.Alias(sdim,         &d_array[dex += dim]         );
	fUnitNorm_0.Alias     (sdim,         &d_array[dex += dim]         );
	fUnitNorm_j.Alias     (sdim,         &d_array[dex += dim]         );
	fInternal_0.Alias     (kNumInternal, &d_array[dex += dim]         );
	fInternal_j.Alias     (kNumInternal, &d_array[dex += kNumInternal]);
	fGradIsoHard_0.Alias  (kNSD,1,       &d_array[dex += kNumInternal]);
	fGradIsoHard_j.Alias  (kNSD,1,       &d_array[dex += kNSD        ]);
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
dMatrixT GradJ2SSKStV1D::Grad1R(double lambda, dMatrixT gradlambda, double laplambda)
{
#pragma unused(laplambda)

	double fdR    = dK(lambda);

	fTensorTemp2.SetToScaled(fdR, gradlambda);
	return fTensorTemp2;
}
	
/** 2nd gradient of Isotropic Hardening conjugate force */
double GradJ2SSKStV1D::Grad2R(double lambda, dMatrixT gradlambda, double laplambda)
{
	double fdR    = dK(lambda);
	double fddR   = ddK(lambda);

	return fddR*gradlambda.ScalarProduct() + fdR*laplambda;
}
	
/** 3rd gradient of Isotropic Hardening conjugate force */
dMatrixT GradJ2SSKStV1D::Grad3R(double lambda, dMatrixT gradlambda, double laplambda)
{
	double fddR   = ddK(lambda);
	double fdddR  = dddK(lambda);

	fTensorTemp2.SetToScaled(fdddR*gradlambda.ScalarProduct() + 3*fddR*laplambda, gradlambda);
	return fTensorTemp2;
}
	
/** 4th gradient of Isotropic Hardening conjugate force */
double GradJ2SSKStV1D::Grad4R(double lambda, dMatrixT gradlambda, double laplambda)
{
	double fddR     = ddK(lambda);
	double fdddR    = dddK(lambda);
	double fddddR   = ddddK(lambda);

	/* increment in field */
	return fddddR*pow(gradlambda.ScalarProduct(),2) + 6*fdddR*gradlambda.ScalarProduct()*laplambda + 3*fddR*pow(laplambda,2);
}

/* yield criteria moduli */
double GradJ2SSKStV1D::YieldCondition(double isohard, dMatrixT gradisohard, double lapisohard)
{
	double& fYieldStrength = fInternal_j[kYieldStrength];
	fYieldStrength = K(isohard);
	
	if ( (fYieldStrength < kYieldStressSmall || fInternal_0[kWeakened] > 0.5 || fInternal_j[kWeakened] > 0.5) && false)
	{			
		if (fInternal_0[kWeakened] < 0.5)
		{				
			cout << "Element [" << CurrElementNumber() << "], ip [" << CurrIP() << "]:  ";
			cout << "Isotropic Hardening: " << isohard;
			cout << "       Yield Stress: " << fYieldStrength << endl;

			fInternal_j[kWeakened] = 1;
		}

		/* check yield strength */
		return  sqrt(fStress.ScalarProduct()) - ( kYieldStressSmall + fc_r*Grad2R(isohard, gradisohard, lapisohard));
	}
	else
		/* check yield strength */
		return  sqrt(fStress.ScalarProduct()) - ( fYieldStrength + fc_r*Grad2R(isohard, gradisohard, lapisohard));
}
