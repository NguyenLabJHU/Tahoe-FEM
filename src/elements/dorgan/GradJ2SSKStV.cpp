/* $Id: GradJ2SSKStV.cpp,v 1.1 2004-09-02 18:25:04 rdorgan Exp $ */
#include "GradJ2SSKStV.h"
#include "GradSSMatSupportT.h"
#include "ElementCardT.h"
#include "StringT.h"

#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "iArrayT.h"

#include <math.h>

using namespace Tahoe;

/* parameters */
const double sqrt23 = sqrt(2.0/3.0);
const int    kNumInternal = 7;
const int    kNSD         = 3;
const double kYieldTol    = 1.0e-10;
const double kYieldStressSmall = 1.0e-03;

/* element output data */
const int kNumOutput = 5;
static const char* Labels[kNumOutput] = {
	"alpha",       // equivalent plastic strain
	"VM",          // Von Mises stress
	"press",       // pressure
	"r",           // isotropic hardening
	"r,xx"};       // Laplacian isotropic hardening

/* constructor */
GradJ2SSKStV::GradJ2SSKStV(void):
	ParameterInterfaceT ("grad_small_strain_StVenant_J2"),
	HookeanMatT(kNSD),

	fStress(kNSD),
	fRelStress(kNSD),
	fElasticStrain(kNSD),

	fK(NULL),

	fNumIP(-1),
	
	//	fModulus               (dSymMatrixT::NumValues(kNSD)),
	fOffDiagonalModulus_bh (dSymMatrixT::NumValues(kNSD),1),
	fOffDiagonalModulus_hb (1,dSymMatrixT::NumValues(kNSD)),
	fGradientModulus_hh    (1),
	fGradientModulus_hp    (1,kNSD),
	fGradientModulus_hq    (1),

	fTensorTemp1(dSymMatrixT::NumValues(kNSD),1),
	fTensorTemp2(kNSD,1),
	fTensorTemp3(dSymMatrixT::NumValues(kNSD)),
	fTensorTemp4(kNSD)
{
	const char caller[] = "GradSmallStrainT::GradJ2SSKStV";
}

/** destructor */
GradJ2SSKStV::~GradJ2SSKStV(void) { delete fK; }

/* update internal variables */
void GradJ2SSKStV::UpdateHistory(void)
{
	ElementCardT& element = CurrentElement();

	/* update plastic variables */
	for (int ip = 0; ip < fNumIP; ip++)
	{
		LoadData(element, ip);

		/* update state */
		fPlasticStrain_0 = fPlasticStrain_j;
		fInternal_0      = fInternal_j;
		fGradIsoHard_0   = fGradIsoHard_j;
	}
}

/* reset internal variables to last converged solution */
void GradJ2SSKStV::ResetHistory(void)
{
	ElementCardT& element = CurrentElement();

	for (int ip = 0; ip < fNumIP; ip++)
	{
		LoadData(element, ip);

		/* reset state */
		fPlasticStrain_j = fPlasticStrain_0;
		fInternal_j      = fInternal_0;
		fGradIsoHard_j   = fGradIsoHard_0;
	}
}

/* update flag describing a weakened ip */
void GradJ2SSKStV::UpdateWeakened(const ElementCardT& element, int ip)
{
	LoadData(element, ip);

	/* update state */
	fInternal_j[kWeakened] = 1;
}

/* reset flag describing a weakened ip */
void GradJ2SSKStV::ResetWeakened(const ElementCardT& element, int ip)
{
	LoadData(element, ip);

	/* update state */
	fInternal_j[kWeakened] = 0;
}

/* modulus */
const dMatrixT& GradJ2SSKStV::c_ijkl(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	return fEModulus;
}

/* off diagonal moduli for Kar */
const dMatrixT& GradJ2SSKStV::odm_bh_ij(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	/* off-diagonal modulus */
	fOffDiagonalModulus_bh = 0.;
	fEModulus.Multx(fUnitNorm,fTensorTemp1);
	fOffDiagonalModulus_bh.AddScaled(-1.0, fTensorTemp1);
	return fOffDiagonalModulus_bh;
}

/* off diagonal moduli for K_ra */
const dMatrixT& GradJ2SSKStV::odm_hb_ij(void)
{
	/* off-diagonal modulus */
	fOffDiagonalModulus_hb = 0.;
	fOffDiagonalModulus_hb.Transpose(odm_bh_ij());
	return fOffDiagonalModulus_hb;
}

/* moduli for local term in K_hh */
const dMatrixT& GradJ2SSKStV::gm_hh(void)
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
		double temp;
 
		/* gradient modulus */
		fEModulus.Multx(fUnitNorm,fTensorTemp1);
		temp = dMatrixT::Dot(fUnitNorm, fTensorTemp1);
		fGradientModulus_hh[0] = temp + sqrt23*(sqrt23-fk_r*fR+fc_r*fk_r*fGrad2R) * (fdR+fc_r*fdddR*fGradIsoHard_j.ScalarProduct()+fc_r*fddR*fLapIsoHard) - sqrt23*fc_r*fk_r*(2*fddR*(dMatrixT::Dot(fGrad1R,fGradIsoHard_j)+fc_r*dMatrixT::Dot(fGrad3R,fGradIsoHard_j)) + (fGrad2R+fc_r*fGrad4R)*fdR);

		if (fInternal_0[kWeakened] > 0.5 || fInternal_j[kWeakened] > 0.5)
			fGradientModulus_hh[0] -= sqrt23*(sqrt23-fk_r*fR)*fdR;
	}
	else
		fGradientModulus_hh[0] = Young();
		
	return fGradientModulus_hh;
}

/* moduli for gradient term in K_hp */
const dMatrixT& GradJ2SSKStV::gm_hp(void)
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
	fGradientModulus_hp.SetToScaled(2*fc_r*sqrt23*fddR*(sqrt23-fk_r*fR-fk_r*fc_r*fGrad2R), fGradIsoHard_j);
	fGradientModulus_hp.AddScaled(-2*fk_r*fc_r*sqrt23*fdR,fGrad1R);
	fGradientModulus_hp.AddScaled(-2*fk_r*fc_r*sqrt23*fc_r*fdR,fGrad3R);
	return fGradientModulus_hp;
}

/* moduli for gradient term in K_hq */
const dMatrixT& GradJ2SSKStV::gm_hq(void)
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
	fGradientModulus_hq = fc_r*fdR*sqrt23*(sqrt23-fk_r*fR-fk_r*fc_r*fGrad2R);
	return fGradientModulus_hq;
}

/* stress */
const dSymMatrixT& GradJ2SSKStV::s_ij(void)
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
	fIsoHard += sqrt23*del_Lambda();
	fGradIsoHard_j.AddScaled(sqrt23, del_GradLambda());
	fLapIsoHard += sqrt23*del_LapLambda();

	/* compute trial elastic strain */
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_trl = ElasticStrain(e_tot, element, ip);

	/* compute trial elastic stress */
	HookeanStress(e_trl, fStress);

	/* compute trial relative stress */
	fRelStress.Deviatoric(fStress);

	/* compute trial yield condition */
	ftrial = YieldCondition(fIsoHard, fGradIsoHard_j, fLapIsoHard);

	int iteration = fGradSSMatSupport->GroupIterationNumber();

	if (iteration <= -1)
	{
		/* elastic iteration */
		fYieldStep = 0;

		fUnitNorm = 0.;
		ftrial = 0.;
	}
	else if (ftrial > kYieldTol || fYieldStep > 0.5)   // check elastic predictor
	{
		/* plastic iteration */
		fYieldStep = 1.;

		/* compute unit normal */
		fUnitNorm.SetToScaled(1.0/sqrt(fRelStress.ScalarProduct()), fRelStress);
		
 		/* update plastic strain */
  		fPlasticStrain_j.AddScaled(del_Lambda(), fUnitNorm);

 		/* plastic corrector */
		HookeanMatT::Modulus().Multx(fUnitNorm,fTensorTemp1);
 		fStress.AddScaled(-del_Lambda(), fTensorTemp1);
		fRelStress.Deviatoric(fStress);

 		/* update yield condition */
 		ftrial = YieldCondition(fIsoHard, fGradIsoHard_j, fLapIsoHard);
	}
	else
	{
		/* elastic iteration */
		fUnitNorm = 0.;
		ftrial = 0.;
	}

	/* compute consistent elastic tangent moduli */
	TangentModulus();
	
	return fStress;
}

/* unit normal */
const dSymMatrixT& GradJ2SSKStV::n_ij(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	/* yield condition */
	return fUnitNorm;
}

/* yield criteria moduli */
double GradJ2SSKStV::yc(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	/* yield condition */
	return fInternal_j[kYieldCrt];
}

/* yield strenth */
double GradJ2SSKStV::ys(void)
{
	/* load internal variables */
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	LoadData(element, ip);

	/* yield condition */
	return fInternal_j[kYieldStrength];
}

/** returns 1 if the ip has weakened during the current time step, 0 otherwise */
int GradJ2SSKStV::weakened(void)
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
double GradJ2SSKStV::StrainEnergyDensity(void)
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
int GradJ2SSKStV::NumOutputVariables(void) const  { return kNumOutput; }
void GradJ2SSKStV::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Dimension(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void GradJ2SSKStV::ComputeOutput(dArrayT& output)
{
	/* stress tensor (loads element data and sets fStress) */
	s_ij();

	/* pressure */
	output[2] = fStress.Trace()/3.0;
	
	/* deviatoric Von Mises stress */
	fStress.Deviatoric();
	double J2 = fStress.Invariant2();
	J2 = (J2 < 0.0) ? 0.0 : J2;
	output[1] = sqrt(3.0*J2);

	/* equivalent plastic strain */
	output[0] = sqrt23*sqrt(fPlasticStrain_j.ScalarProduct());

	/* isotropic hardening */
	output[3] = fInternal_j[kIsoHard];

	/* Laplacian isotropic hardening */
	output[4] = fInternal_j[kLapIsoHard];
}

/* implementation of the ParameterInterfaceT interface */
void GradJ2SSKStV::DefineParameters(ParameterListT& list) const
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
void GradJ2SSKStV::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	GradSSSolidMatT::DefineSubs(sub_list);
	IsotropicT::DefineSubs(sub_list);
	HookeanMatT::DefineSubs(sub_list);

	/* hardening function */
	sub_list.AddSub("hardening_function_choice", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void GradJ2SSKStV::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
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
		sub_lists.AddSub("power_law");
		sub_lists.AddSub("piecewise_linear");
	}
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* GradJ2SSKStV::NewSub(const StringT& name) const
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
void GradJ2SSKStV::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "GradJ2SSKStV::TakeParameterList";

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
const dSymMatrixT& GradJ2SSKStV::ElasticStrain(const dSymMatrixT& totalstrain,
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
double GradJ2SSKStV::del_Lambda(void)
{
	/* increment in field */
	return GradSSSolidMatT::Lambda()[0] - Lambda_last()[0];
}

/* incremental change in gradient of field */
dMatrixT GradJ2SSKStV::del_GradLambda(void)
{
	fTensorTemp2 = 0.;
	
	/* increment in field */
	/* third component set to zero for 2D */
	for (int sd = 0; sd < NumSD(); sd++)
		fTensorTemp2[sd] = GradLambda()[sd] - GradLambda_last()[sd];
	
	return fTensorTemp2;
}

/* incremental change in laplacian of field */
double GradJ2SSKStV::del_LapLambda(void)
{
	/* increment in field */
	return LapLambda()[0] - LapLambda_last()[0];
}

/* set modulus */
void GradJ2SSKStV::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}

/* return a pointer to a new element object constructed with
* the data from element */
void GradJ2SSKStV::AllocateAllElements(void)
{
	/* determine storage */
	int i_size = 0;
	// remove fFlags!
	i_size += fNumIP; //fFlags

	int d_size = 0;
	int dim = dSymMatrixT::NumValues(kNSD);
	d_size += dim;          //fUnitNorm
	d_size += dim;	        //fPlasticStrain_0
	d_size += dim;	        //fPlasticStrain_j
	d_size += kNumInternal; //fInternal_0
	d_size += kNumInternal; //fInternal_j
	d_size += kNSD*1;       //fGradIsoHard_0
	d_size += kNSD*1;       //fGradIsoHard_j
	d_size += dim*dim;      //fEModulus

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
void GradJ2SSKStV::LoadData(const ElementCardT& element, int ip)
{
	/* fetch arrays */
	const dArrayT& d_array = element.DoubleData();
	
	/* decode */
	dSymMatrixT::DimensionT sdim = dSymMatrixT::int2DimensionT(kNSD);
	int dim   = dSymMatrixT::NumValues(kNSD);
	int block = 3*dim + 2*kNumInternal + 2*kNSD*1 + dim*dim;
	int dex   = ip*block;

	fUnitNorm.Alias       (sdim,         &d_array[dex                ]);
	fPlasticStrain_0.Alias(sdim,         &d_array[dex += dim         ]);
	fPlasticStrain_j.Alias(sdim,         &d_array[dex += dim         ]);
	fInternal_0.Alias     (kNumInternal, &d_array[dex += dim         ]);
	fInternal_j.Alias     (kNumInternal, &d_array[dex += kNumInternal]);
	fGradIsoHard_0.Alias  (kNSD,1,       &d_array[dex += kNumInternal]);
	fGradIsoHard_j.Alias  (kNSD,1,       &d_array[dex += kNSD        ]);
	fEModulus.Alias       (dim,dim,      &d_array[dex += kNSD        ]);
}

/** 1st gradient of Isotropic Hardening conjugate force */
dMatrixT GradJ2SSKStV::Grad1R(double lambda, dMatrixT gradlambda, double laplambda)
{
#pragma unused(laplambda)

	double fdR    = dK(lambda);

	fTensorTemp2.SetToScaled(fdR, gradlambda);
	return fTensorTemp2;
}
	
/** 2nd gradient of Isotropic Hardening conjugate force */
double GradJ2SSKStV::Grad2R(double lambda, dMatrixT gradlambda, double laplambda)
{
	double fdR    = dK(lambda);
	double fddR   = ddK(lambda);


	return fddR*gradlambda.ScalarProduct() + fdR*laplambda;
}
	
/** 3rd gradient of Isotropic Hardening conjugate force */
dMatrixT GradJ2SSKStV::Grad3R(double lambda, dMatrixT gradlambda, double laplambda)
{
	double fddR   = ddK(lambda);
	double fdddR  = dddK(lambda);

	fTensorTemp2.SetToScaled(fdddR*gradlambda.ScalarProduct() + fddR*laplambda, gradlambda);
	return fTensorTemp2;
}
	
/** 4th gradient of Isotropic Hardening conjugate force */
double GradJ2SSKStV::Grad4R(double lambda, dMatrixT gradlambda, double laplambda)
{
	double fddR     = ddK(lambda);
	double fdddR    = dddK(lambda);
	double fddddR   = ddddK(lambda);

	/* increment in field */
	return fddddR*pow(gradlambda.ScalarProduct(),2) + 2*fdddR*gradlambda.ScalarProduct()*laplambda + fddR*pow(laplambda,2);
}

/* computes the consistent elastic tangent moduli */
void GradJ2SSKStV::TangentModulus()
{
	fEModulus = HookeanMatT::Modulus();

	/* moduli corrections */
	/* H = [Modulus^-1 + dUnitNormal/dsigma * del_Lambda]^-1 */
	if (fInternal_j[kYieldCrt] > 0.5 && true)
	{
		fEModulus.Inverse();

		/* compute constants */
		double thetahat = del_Lambda()/sqrt(fRelStress.ScalarProduct());
		
		//		fTensorTemp3.ReducedIndexDeviatoric();
		//		fEModulus.AddScaled(thetahat, fTensorTemp3);

		fTensorTemp3.ReducedIndexI();
		fEModulus.AddScaled(thetahat, fTensorTemp3);

		fTensorTemp4.Identity();
		fTensorTemp3.Outer(fTensorTemp4, fTensorTemp4);
		fEModulus.AddScaled(-1.0/3.0*thetahat, fTensorTemp3);

		fTensorTemp3.Outer(fUnitNorm,fUnitNorm);
		fEModulus.AddScaled(-1.0*thetahat, fTensorTemp3);

		fEModulus.Inverse();
	}
}

/* yield criteria moduli */
double GradJ2SSKStV::YieldCondition(double isohard, dMatrixT gradisohard, double lapisohard)
{
	double& fYieldStrength = fInternal_j[kYieldStrength];
	fYieldStrength = K(isohard);
	
	if ( (fYieldStrength < kYieldStressSmall || fInternal_0[kWeakened] > 0.5 || fInternal_j[kWeakened] > 0.5) && false )
	{			
		if (fInternal_0[kWeakened] < 0.5)
		{				
			cout << "Element [" << CurrElementNumber() << "], ip [" << CurrIP() << "]:  ";
			cout << "Isotropic Hardening: " << isohard;
			cout << "       Yield Stress: " << fYieldStrength << endl;

			fInternal_j[kWeakened] = 1;
		}

		/* check yield strength */
		return  sqrt(fRelStress.ScalarProduct()) - sqrt23 * ( kYieldStressSmall + fc_r*Grad2R(isohard, gradisohard, lapisohard));
	}
	else
		/* check yield strength */
		return  sqrt(fRelStress.ScalarProduct()) - sqrt23 * ( fYieldStrength + fc_r*Grad2R(isohard, gradisohard, lapisohard));
}
