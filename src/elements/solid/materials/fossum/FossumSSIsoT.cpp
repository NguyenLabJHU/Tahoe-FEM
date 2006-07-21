/* 3-invariant, single-surface dilation/compaction plasticity model
 * with isotropic and kinematic hardening
 * Implemented 8/02 Craig Foster
 */

#include "FossumSSIsoT.h"
#include "SSEnhLocMatSupportT.h"

using namespace Tahoe;

const int kNumInternal = 3; // number of internal state variables
const double sqrt3 = sqrt(3.0);
const double kYieldTol = 1.0e-10;
const int kNSD = 3;

/* element output data */
const int kNumOutput = 10;
static const char* Labels[kNumOutput] = {
	"alpha11",  // back stress
	"alpha22",  
	"alpha33",
	"alpha23",
	"alpha13",
	"alpha12",
	"kappa",
	"meanstress",
	"J2",
	"J3"
};

/*constructor*/
FossumSSIsoT::FossumSSIsoT(void):
	ParameterInterfaceT("Fossum_small_strain"),
	fA(-1.0),         
	fB(-1.0),
	fC(-1.0),
	fTheta(-1.0),
	fR(-1.0),        
	fKappa0(0.0),
	fW(0.0),
	fD1(0.0), 
	fD2(0.0),
	fCalpha(-1.0),   
	fPsi(-1.0),
	fN(-1.0),

    HookeanMatT(3),
      
    fElasticStrain(kNSD),
    fModuliCorr(dSymMatrixT::NumValues(kNSD)),
    fModuliCorrPerfPlas(dSymMatrixT::NumValues(kNSD)),
    One(kNSD),
    fBackStress(kNSD),
    fDeltaAlpha(kNSD),
    fKappaDummy(0),
    fKappaCapped(fKappaDummy),

    /*spectral decomp parameters */
    spectre(kNSD),
    m(kNSD),
    principalEqStress(kNSD),

    /*viscous paramters*/
	fTimeFactor(1.0),
        
	fSSEnhLocMatSupport(NULL)
{

    /* allocate space for principal dirs m */
    for (int i = 0; i < 3; i++)
    {
        m[i].Dimension(kNSD);
    }

    /* initialize constant tensor */
    One.Identity();

}


/* destructor */
FossumSSIsoT::~FossumSSIsoT(void) { }

/*describe parameters needed by the interface*/
void FossumSSIsoT::DefineParameters(ParameterListT& list) const  
{      
	/* inherited */
	SSIsotropicMatT::DefineParameters(list);
	HookeanMatT::DefineParameters(list);
 
  ParameterT A(fA, "shear_surface_parameter__A__stress");    
  A.AddLimit(0.0, LimitT::LowerInclusive); 
  list.AddParameter(A);

  list.AddParameter(fB, "shear_surface_parameter__B__1_by_stress"); 
  list.AddParameter(fC, "shear_surface_parameter__C__stress");
  
  ParameterT theta(fTheta, "shear_surface_parameter__theta__radians");    
  theta.AddLimit(0.0, LimitT::LowerInclusive); 
  list.AddParameter(theta);

  ParameterT R(fR, "cap_ratio__R__dimensionless");    
  R.AddLimit(0.0, LimitT::LowerInclusive); 
  list.AddParameter(R);

  list.AddParameter(fKappa0, "initial_cap_position__kappa_0__stress");

  list.AddParameter(fW, "cap_growth_parameter__W__dimensionless"); 
  list.AddParameter(fD1, "cap_growth_parameter__D1__1_by_stress");
  list.AddParameter(fD2, "cap_growth_parameter__D2__1_by_stress_squared");

  ParameterT calpha(fCalpha, "shear_surface_growth_factor__c_alpha__stress");
  calpha.AddLimit(0.0, LimitT::LowerInclusive);
  list.AddParameter(calpha);     

  ParameterT N(fN, "shear_surface_offset__N__stress");
  N.AddLimit(0.0, LimitT::LowerInclusive);
  list.AddParameter(N);

  ParameterT psi(fPsi, "triaxial_compression_to_extension_strength_ratio__psi__dimensionless");
  psi.AddLimit(0.5, LimitT::LowerInclusive);//.69 for Gudheus, but future version may be as low as .5
  psi.AddLimit(2.0, LimitT::UpperInclusive);
  list.AddParameter(psi);

  ParameterT fluidity(fFluidity, "fluidity_parameter__eta__dimensionless");
  fluidity.AddLimit(0.0, LimitT::LowerInclusive);
  list.AddParameter(fluidity);

  ParameterT debug(fFossumDebug, "local_debug_parameter");
  debug.SetDefault(false);
  list.AddParameter(debug);
}

/* accept parameter list */
void FossumSSIsoT::TakeParameterList(const ParameterListT& list)
{
	fNumIP = NumIP();

    /* inherited */
    SSIsotropicMatT::TakeParameterList(list);
    HookeanMatT::TakeParameterList(list);
    
    // dimension
    fStress.Dimension(3);
    fSigma.Dimension(3);
    fStrain.Dimension(3);
    fModulus.Dimension(dSymMatrixT::NumValues(3));
    fModulusPerfPlas.Dimension(dSymMatrixT::NumValues(3));
    fModulusContinuum.Dimension(dSymMatrixT::NumValues(3));
    fModulusContinuumPerfPlas.Dimension(dSymMatrixT::NumValues(3));

    fA = list.GetParameter("shear_surface_parameter__A__stress");
    fB = list.GetParameter("shear_surface_parameter__B__1_by_stress");
    fC = list.GetParameter("shear_surface_parameter__C__stress");
    fTheta = list.GetParameter("shear_surface_parameter__theta__radians");
    fR = list.GetParameter("cap_ratio__R__dimensionless");
    fKappa0 = list.GetParameter("initial_cap_position__kappa_0__stress");
    fW = list.GetParameter("cap_growth_parameter__W__dimensionless");
    fD1 = list.GetParameter("cap_growth_parameter__D1__1_by_stress");
    fD2 = list.GetParameter("cap_growth_parameter__D2__1_by_stress_squared");
    fCalpha = list.GetParameter("shear_surface_growth_factor__c_alpha__stress");
    fN = list.GetParameter("shear_surface_offset__N__stress");
    fPsi = list.GetParameter("triaxial_compression_to_extension_strength_ratio__psi__dimensionless");
    fFluidity = list.GetParameter("fluidity_parameter__eta__dimensionless");
    fFossumDebug = list.GetParameter("local_debug_parameter");

    fmu = Mu();
    flambda = Lambda();
    
    /* cast to small strain embedded discontinuity material pointer */
	fSSEnhLocMatSupport = TB_DYNAMIC_CAST(const SSEnhLocMatSupportT*, fSSMatSupport);
}

void FossumSSIsoT::DefineSubs(SubListT& sub_list) const
{
    /* inherited */
    SSIsotropicMatT::DefineSubs(sub_list);
    HookeanMatT::DefineSubs(sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FossumSSIsoT::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* params = SSIsotropicMatT::NewSub(name);
	if (params) 
		return params;
	else
		return HookeanMatT::NewSub(name);
}


        
/*  protected: */

/*
 * Returns the value of the yield function given the
 * Cauchy stress vector and state variables, where alpha is
 * the back stress and kappa an internal state variable
 */
double FossumSSIsoT::YieldCondition(const dSymMatrixT& stress, 
									const double kappa, dSymMatrixT& backStress)
{
	dSymMatrixT equivalentStress(kNSD); 
	double I1, J2, J3;

	equivalentStress = 0.0;
	equivalentStress.DiffOf(stress,backStress);

	/* calculate stress invariants */
	I1 = equivalentStress.Trace();
	/* make deviatoric */
	equivalentStress.PlusIdentity(-I1/3.0);
	J2 = .5*equivalentStress.ScalarProduct();
	J3 = equivalentStress.Det();

	return YieldFn(I1, J2, J3, kappa);
}

double FossumSSIsoT::YieldFn(double I1, double J2, double J3, double kappa)
{
	return YieldFnGamma(J2, J3)*YieldFnGamma(J2, J3)*J2 - YieldFnFfMinusN(I1)*YieldFnFfMinusN(I1)*YieldFnFc(I1, kappa);
}

/*-----------------------------------------------------------------*/
//DPSSKStV-like fns

/* form of tangent matrix (symmetric by default) */
GlobalT::SystemTypeT FossumSSIsoT::TangentType(void) const { return GlobalT::kNonSymmetric; }

/* update internal variables */
void FossumSSIsoT::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Update(element);
}

/* reset internal variables to last converged solution */
void FossumSSIsoT::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Reset(element);
}

/* elastic modulus */
const dMatrixT& FossumSSIsoT::ce_ijkl(void)
{
	fModulusCe = HookeanMatT::Modulus();
	return fModulusCe;
}

/* continuum modulus */
const dMatrixT& FossumSSIsoT::con_ijkl(void)
{
	//elastoplastic correction
	fModulusContinuum.SumOf( HookeanMatT::Modulus(), ModuliCorrection(CurrentElement(),CurrIP()) );
	return fModulusContinuum;
}

/* perfectly plastic continuum modulus */
const dMatrixT& FossumSSIsoT::con_perfplas_ijkl(void)
{
	//elastoplastic correction
	fModulusContinuumPerfPlas.SumOf( HookeanMatT::Modulus(), ModuliCorrPerfPlas(CurrentElement(),CurrIP()) );
	return fModulusContinuumPerfPlas;
}

/* returns the strain energy density for the specified strain */
double FossumSSIsoT::StrainEnergyDensity(void)
{
	return HookeanEnergy(ElasticStrain(e(), CurrentElement(), CurrIP()));
}

/* returns the number of variables computed for nodal extrapolation
 * during for element output, ie. internal variables. Returns 0
 * by default. */
int FossumSSIsoT::NumOutputVariables(void) const  { return kNumOutput; } 

void FossumSSIsoT::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Allocate(kNumOutput);
        
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++) labels[i] = Labels[i];
}


void FossumSSIsoT::ComputeOutput(dArrayT& output)
{
	const ElementCardT& element = CurrentElement();
	int ip = CurrIP();
	dMatrixT Ce = HookeanMatT::Modulus();

	/*OUTPUT FOR ALPHA, KAPPA */ 
	if (element.IsAllocated())
	{
	  	LoadData(element, ip);
		for (int i = 0; i < 6 ; i++) output [i] = fBackStress [i] + fDeltaAlpha[i];
		output [6] = fInternal[kkappa] + fInternal[kdeltakappa];
	}
	else
	{       
		output [0] = 0.0;
		output [1] = 0.0;
		output [2] = 0.0;
		output [3] = 0.0;
		output [4] = 0.0;
		output [5] = 0.0;
		output [6] = fKappa0;
	}
	
	/* stress tensor (load state) */
	const dSymMatrixT& stress = s_ij();

	/* pressure */
	output[7] = fStress.Trace()/3.0;

	/* J2 invariant */
	fStress.Deviatoric();
	double J2 = fStress.Invariant2();
	J2 = (J2 < 0.0) ? 0.0 : J2;
	output[8] = J2;
        
	/* J3 - 3rd invariant of deviatoric stress */
	double J3 = fStress.Det();
	output[9] = J3;
	
	/* Clean up - restore full stress state */ 
	fStress.AddScaled(output[7], One);
}

/*************************************************************************
* Protected
*************************************************************************/

/* set modulus */
void FossumSSIsoT::SetModulus(dMatrixT& modulus)
{
	IsotropicT::ComputeModuli(modulus);
}

/*-----------------------------------------------------------------------*/
//LinHardT

/* returns elastic strain */
const dSymMatrixT& FossumSSIsoT::ElasticStrain(const dSymMatrixT& totalstrain, 
					const ElementCardT& element, int ip)
{
	/* remove plastic strain */
	if (element.IsAllocated()) 
	{
		/* load internal variables */
		LoadData(element, ip);

		/* compute elastic strain */
		fElasticStrain.DiffOf(totalstrain, fPlasticStrain);
        
		return fElasticStrain;
	}        
	/* no plastic strain */
	else        
		return totalstrain;
}


/* correction for continuum elasto-plastic tangent */
const dMatrixT& FossumSSIsoT::ModuliCorrection(const ElementCardT& element, int ip)
{
	/* initialize */
	fModuliCorr = 0.0;
	
	dArrayT hardeningFns(7), dfdq(7);
	dSymMatrixT dfdSigma(3), dfdAlpha(3);
	dMatrixT Ce = HookeanMatT::Modulus();

	if (element.IsAllocated() && (element.IntegerData())[ip] == kIsPlastic)
	{
	
		/* load internal state variables */
		LoadData(element,ip);
		
		double kappa = fInternal[kkappa] + fInternal[kdeltakappa];
		dSymMatrixT alpha(3); 
		alpha.SumOf(fBackStress, fDeltaAlpha);

		/*Find Invariants */
		double I1 = 0.0, J2 =0.0, J3 = 1.0;
		for (int i = 0; i < kNSD; i++)
			I1 += principalEqStress[i];
		for (int i = 0; i < kNSD; i++)
		{
			J2 += 0.5*(principalEqStress[i] - I1/3.0) * (principalEqStress[i] - I1/3.0);
			J3 *= (principalEqStress[i] - I1/3.0);
		}

		dfdSigma = DfdSigma(I1,J2,J3, kappa, principalEqStress, m);
		dfdAlpha = DfdAlpha(I1,J2,J3, kappa, principalEqStress, m);
		hardeningFns = Hardening(I1,J2,J3, kappa, principalEqStress, m, alpha);                             
		for (int i = 0; i < 6; i++) dfdq[i] = dfdAlpha[i];
		dfdq[6] = dfdKappa(I1, kappa);

		double chi = dfdSigma.B_ij_A_ijkl_B_kl(Ce) - dfdq.Dot(dfdq,hardeningFns);
		dSymMatrixT CeTimesdfdSigma(3);
		CeTimesdfdSigma.A_ijkl_B_kl(Ce, dfdSigma);
		dMatrixT corr(6);
		corr.Outer(CeTimesdfdSigma, CeTimesdfdSigma);

		fModuliCorr.AddScaled(-1.0/chi, corr);
		
	}

	return fModuliCorr;
}       

/* correction for continuum elasto-perfectly-plastic tangent */
const dMatrixT& FossumSSIsoT::ModuliCorrPerfPlas(const ElementCardT& element, int ip)
{
	/* initialize */
	fModuliCorrPerfPlas = 0.0;

	dSymMatrixT dfdSigma(3);
	dMatrixT Ce = HookeanMatT::Modulus();

	if (element.IsAllocated() && (element.IntegerData())[ip] == kIsPlastic)
	{
	
		/* load internal state variables */
		LoadData(element,ip);
		
		double kappa = fInternal[kkappa] + fInternal[kdeltakappa];

		/*Find Invariants */
		double I1 = 0.0, J2 =0.0, J3 = 1.0;
		for (int i = 0; i < kNSD; i++)
			I1 += principalEqStress[i];
		for (int i = 0; i < kNSD; i++)
		{
			J2 += 0.5*(principalEqStress[i] - I1/3.0) * (principalEqStress[i] - I1/3.0);
			J3 *= (principalEqStress[i] - I1/3.0);
		}

		dfdSigma = DfdSigma(I1,J2,J3, kappa, principalEqStress, m);

		double chi = dfdSigma.B_ij_A_ijkl_B_kl(Ce);
		dSymMatrixT CeTimesdfdSigma(3);
		CeTimesdfdSigma.A_ijkl_B_kl(Ce, dfdSigma);
		dMatrixT corr(6);
		corr.Outer(CeTimesdfdSigma, CeTimesdfdSigma);

		fModuliCorrPerfPlas.AddScaled(-1.0/chi, corr);		
	}

	return fModuliCorrPerfPlas;
}       

                 
/* return a pointer to a new plastic element object constructed with
 * the data from element */
void FossumSSIsoT::AllocateElement(ElementCardT& element)
{
	/* determine storage */
	int i_size = 0;
	i_size += fNumIP; //fFlags
	i_size += fNumIP; //kappaCapped

	int d_size = 0;
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fPlasticStrain
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fStress
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; // fBackStress
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; // fDeltaAlpha
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; // fStrain
	d_size += kNSD * dSymMatrixT::NumValues(kNSD)*fNumIP; // m

	d_size += kNSD * fNumIP;              // principal eq stresses
	d_size += kNumInternal*fNumIP;        //fInternal

	/* construct new plastic element */
	element.Dimension(i_size, d_size);
        
	/* initialize values */
	element.IntegerData() = kIsElastic;
	element.DoubleData()  = 0.0;  // initialize all double types to 0.0

	/* initialize kappa to kappa0 */
	for (int i=0; i<fNumIP; i++)
	  {
		element.DoubleData() [((5 + kNSD) * dSymMatrixT::NumValues(kNSD) + kNSD)* fNumIP + kNumInternal * i + kkappa] = fKappa0;
	  } 
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* element level data */
void FossumSSIsoT::Update(ElementCardT& element)
{
	/* get flags */
	iArrayT& Flags = element.IntegerData();

	/* check if reset state */
	if (Flags[0] == kReset)
	{
		Flags = kIsElastic; //don't update again
		return; 
	}

	/* update plastic variables */
	for (int ip = 0; ip < fNumIP; ip++)
		if (Flags[ip] == kIsPlastic) /* plastic update */
		{
			/* do not repeat if called again. */
			Flags[ip] = kIsElastic;
			/* NOTE: ComputeOutput writes the updated internal variables
			*       for output even during iteration output, which is
			*       called before UpdateHistory */

			/* fetch element data */
			LoadData(element, ip);
        
			/* plastic increment */
			double& dgamma = fInternal[kdgamma];

			double I1, J2, J3;
			I1 = 0.0; J2 = 0.0; J3 = 1.0;

			for (int i = 0; i < kNSD; i++) I1 += principalEqStress[i];
			
			for (int i = 0; i < kNSD; i++)
			{
				J2 += 0.5 * (principalEqStress[i] - I1/3.0) * (principalEqStress[i] - I1/3.0);
				J3 *= (principalEqStress[i] - I1/3.0);
			}

			fInternal[kkappa] += fInternal[kdeltakappa];
			fBackStress += fDeltaAlpha;
			
			//alt for visco
			fPlasticStrain = 0.0;

			/* Elastic Strain */
			double meanStrain = fStress.Trace()/(9.0*Kappa());
			fStress.Deviatoric();
			fPlasticStrain = fStress;
			fPlasticStrain /= (2.0*fmu);
			
			fPlasticStrain [0] += meanStrain;
			fPlasticStrain [1] += meanStrain;
			fPlasticStrain [2] += meanStrain;

			/* subtract elastic strain from total strain */
			fPlasticStrain *= -1.0;
			fPlasticStrain += fStrain;

			/*keep cap from retracting ????*/
			//fKappa0 = fInternal[kkappa]; 
		}
}



/* resets to the last converged solution */
void FossumSSIsoT::Reset(ElementCardT& element)
{
	/* flag not to update again */
	(element.IntegerData()) = kReset;
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* load element data for the specified integration point */
void FossumSSIsoT::LoadData(const ElementCardT& element, int ip)
{
	/* check */
	if (!element.IsAllocated()) throw ExceptionT::kGeneralFail;

	/* fetch arrays */
	const dArrayT& d_array = element.DoubleData();
        
	/* decode */
	dSymMatrixT::DimensionT dim = dSymMatrixT::int2DimensionT(kNSD);
	int stressdim = dSymMatrixT::NumValues(kNSD);
	int offset    = stressdim*fNumIP;
	int dex       = ip*stressdim;
        
	fPlasticStrain.Alias(         dim, &d_array[           dex]);
	fStress.Alias       (         dim, &d_array[  offset + dex]);   
	fBackStress.Alias   (         dim, &d_array[2*offset + dex]);
	fDeltaAlpha.Alias   (         dim, &d_array[3*offset + dex]);
	fStrain.Alias       (         dim, &d_array[4*offset + dex]);
	for (int i = 0; i < kNSD; i++)
	  m[i].Alias        (         dim, &d_array[(5+i)*offset + dex]);

	principalEqStress.Alias(     kNSD, &d_array[(5 + kNSD)*offset + ip * kNSD ]); 
	fInternal.Alias     (kNumInternal, &d_array[(5 + kNSD)*offset + fNumIP * kNSD + ip*kNumInternal]);

	const iArrayT& i_array = element.IntegerData();
	fKappaCapped = i_array [fNumIP + ip];

}


/*---------------------------------------------------------------*/
/* The following functions are  auxiliary functions for the yield function */

double FossumSSIsoT::YieldFnFfMinusN(double I1)
{
	return fA - fC*exp(fB*I1) - fTheta * I1 - fN;
}

double FossumSSIsoT::YieldFnFc(double I1, const double kappa)
{
	return 1 - HeavisideFn(kappa - I1)*(I1 - kappa)*(I1 - kappa)/((Xfn(kappa) - kappa)*(Xfn(kappa) - kappa));
}

int FossumSSIsoT::HeavisideFn(double arg)
{
	if (arg < 0)
		return 0;
	else
		return 1;
}

double FossumSSIsoT::Xfn(const double kappa)
{
	return kappa - fR * YieldFnFf(kappa);
}

double FossumSSIsoT::YieldFnFf(double I1)
{
	return fA - fC*exp(fB*I1) - fTheta * I1;
}

/* --------------------------------------------------------------------*/

/* stress */
const dSymMatrixT& FossumSSIsoT::s_ij(void)
{
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();

#ifdef ENHANCED_STRAIN_LOC_DEV	
	int element_locflag = 0;
	if (element.IsAllocated()) 
	{
		element_locflag = fSSEnhLocMatSupport->ElementLocflag();
	}
	if ( element_locflag == 2 )
	{
		fStress = fSSEnhLocMatSupport->ElementStress(ip);
	}
	else
	{
		fStress = sigma_ij();
	}
#else
	fStress = sigma_ij();
#endif

	return fStress;
}


/*-------------------------------------------------------------*/
/*Return Mapping algorithm */

/* stress for non-localized elements */
const dSymMatrixT& FossumSSIsoT::sigma_ij(void)
{

  double yieldFnTol = 1.0e-8; 
  int ip = CurrIP();
  ElementCardT& element = CurrentElement();
 
  /* strains and elastic stress*/
  /* Note ElasticStrain loads ISV's if element is allocated */
  const dSymMatrixT& e_els = ElasticStrain(e(), element, ip);

  /* if stress has been solved for, do not resolve */
  /*
  if (fSSMatSupport->RunState() != GlobalT::kFormRHS && element.IsAllocated())
    return fStress;
  */

  fStrain = e_els;

  if (element.IsAllocated())
    fStrain += fPlasticStrain; //really hackish but nec for 2D
  
  HookeanStress(e_els, fSigma);
  /* working ISV's for iteration */
  dSymMatrixT workingBackStress(kNSD); 
  double workingKappa;
  if (element.IsAllocated())
    {   
      workingBackStress = fBackStress;
      workingKappa = fInternal[kkappa];
    }
  else
    {
      workingBackStress = 0.0;
      workingKappa = fKappa0;
    }

  /*check for yielding */
  double initialYieldCheck;
  initialYieldCheck = YieldCondition(fSigma, workingKappa, workingBackStress);
  if ( initialYieldCheck < yieldFnTol)
    {
      if (element.IsAllocated())
	  {
	    int &flag = (element.IntegerData())[ip]; 
	    flag = kIsElastic;
	  }
      return fSigma;
    }
  else
    {
      /* ELSE plastic loading */

      /* allocate element for plastic variables if not already done */	
      if (!element.IsAllocated()) 
		{
			AllocateElement(element);
			LoadData(element, ip); //resets fSigma and fStrain to 0 on 1st allocation
			HookeanStress(e_els, fSigma);
			fStrain = e_els; //no plastic strain yet
		}

      int &flag = (element.IntegerData())[ip]; 
      flag = kIsPlastic; 

      /*initialize increment variables */
      fInternal[kdgamma] = 0.0;
      fDeltaAlpha = 0.0;
      fInternal[kdeltakappa] = 0.0;

      /*spectral decomposition of equivalent stress*/  
      dSymMatrixT eqTrialStress(kNSD);
      eqTrialStress.DiffOf(fSigma, workingBackStress);
      
      spectre.SpectralDecomp_Jacobi(eqTrialStress, true);
      
      for (int i = 0; i < 3; i++)
		{
			m[i].Outer(spectre.Eigenvectors() [i]);
		}
      
      principalEqStress = spectre.Eigenvalues();
      dArrayT iterationVars(7);

      if(StressPointIteration(initialYieldCheck, iterationVars, workingBackStress, workingKappa))
		{
			/*update stress and ISV's */
			fSigma.AddScaled(iterationVars[0],m[0]);
			fSigma.AddScaled(iterationVars[1],m[1]);
			fSigma.AddScaled(iterationVars[2],m[2]);
	  
			fDeltaAlpha.AddScaled(iterationVars[3],m[0]);
			fDeltaAlpha.AddScaled(iterationVars[4],m[1]);
			fDeltaAlpha.AddScaled(-iterationVars[3]-iterationVars[4],m[2]);
	  
			fInternal[kdeltakappa] += iterationVars[5];
	  
			fInternal[kdgamma] = iterationVars[6];	  
		}
      else //i.e Newton iteration did not converge to valid solution
		{
			cout << "FossumSSIsoT::s_ij, failed to find valid solution\n" << flush;
			throw ExceptionT::kGeneralFail;
		}
    }

  //Rate-Dependence effects. Duvaut-Lions formulation. See Simo and Hughes, p/217
  if (fFluidity != 0.0) //fluidity param = 0 => inviscid case
    {
      double dt = fSSMatSupport->TimeStep(); 
      
      //strains from previous time step
      const dSymMatrixT& e_els_last = ElasticStrain(e_last(), element, ip);
      dSymMatrixT e_tot_last = e_els_last;
      e_tot_last += fPlasticStrain;      

      // stress from previous time step
      dSymMatrixT fSigma_last(3);
      HookeanStress(e_els_last, fSigma_last);
	
      // strain and elastic stress increment
      dSymMatrixT delta_e(3);
      delta_e.DiffOf(fStrain, e_tot_last);
      dSymMatrixT elastic_stress_increment(3);
      HookeanStress(delta_e, elastic_stress_increment);

      /*Set time factor previously set to 1.0*/
      fTimeFactor = 1 - exp(-1*dt/fFluidity);

      fSigma *= fTimeFactor;
      fSigma.AddScaled(1- fTimeFactor, fSigma_last);
      fSigma.AddScaled(fTimeFactor*fFluidity/dt, elastic_stress_increment);
      
      // update isv's
      fDeltaAlpha *= fTimeFactor;
      fInternal[kdeltakappa] *= fTimeFactor;
    }
    
  return fSigma;
}

/* Solver for inviscid stress, via N-R iteration */ 
bool FossumSSIsoT::StressPointIteration(double initialYieldCheck, dArrayT& iterationVars, dSymMatrixT workingBackStress, double workingKappa)
{
	int ip = CurrIP();

	/* initialize iteration variables */
	int newtonCounter = 0, maxIter = 20;
	dArrayT iterationVarsIncr(7), residual(7), residual0(7);
	dSymMatrixT workingStress(kNSD);
	double I1 = 0.0, J2 = 0.0, J3 = 1.0;
	LAdMatrixT dRdX (7); 

	iterationVars = 0.0;
	residual = 0.0;
	residual [6] = initialYieldCheck;
	residual0 = residual;
	
	for (int i = 0; i < kNSD; i++) I1 += principalEqStress[i];

	for (int i = 0; i < kNSD; i++)
	{
		J2 += 0.5*(principalEqStress[i] - I1/3.0) * (principalEqStress[i] - I1/3.0);
		J3 *= (principalEqStress[i] - I1/3.0);
	}

	/*local Newton iteration on variables*/
	while(!ResidualIsConverged(residual, residual0))
	  {
	    /* check to see if not converging */
	    if (newtonCounter++ > maxIter && fFossumDebug)
	      {
			cout << "FossumSSIsoT::StressPointIteration, Newton Iteration failed to converge\n" << flush;
			return false;
	      }

	    /* form dR/dx */
	    dRdX = 0.0; 
	    dRdX = FormdRdX(I1, J2, J3, principalEqStress, workingKappa, fSigma, workingBackStress, iterationVars [6], m);
	    		
		/*solve for dx*/
		iterationVarsIncr.Copy(residual.Pointer());
		iterationVarsIncr*=-1.0;
		iterationVarsIncr = CondenseAndSolve(dRdX, residual);
		  	
		//Do not allow kappa to increase
		if (workingKappa + iterationVarsIncr[5] > fInternal[kkappa])
		  { 
		    iterationVarsIncr = CapKappa(residual, dRdX, workingKappa); 
		    fKappaCapped = true;
		  }
		else
		  fKappaCapped = false;

		/*incr x = x + dx */
		iterationVars += iterationVarsIncr;

		/*update working stress and backstress*/
		principalEqStress [0] += iterationVarsIncr [0] - iterationVarsIncr [3];
		principalEqStress [1] += iterationVarsIncr [1] - iterationVarsIncr [4];
		principalEqStress [2] += iterationVarsIncr [2] + iterationVarsIncr [3] +iterationVarsIncr [4];

		/*reset*/
		I1 = 0.0; 
		J2 = 0.0; 
		J3 = 1.0;

	    for (int i = 0; i < kNSD; i++) I1 += principalEqStress[i];
	    
	    for (int i = 0; i < kNSD; i++)
		  {
			J2 += 0.5 * (principalEqStress[i] - I1/3.0) * (principalEqStress[i] - I1/3.0);
			J3 *= (principalEqStress[i] - I1/3.0);
		  }
	    
	    workingStress.AddScaled(iterationVarsIncr[0],m[0]);
		workingStress.AddScaled(iterationVarsIncr[1],m[1]);
		workingStress.AddScaled(iterationVarsIncr[2],m[2]);
	
		workingBackStress.AddScaled(iterationVarsIncr[3],m[0]);
		workingBackStress.AddScaled(iterationVarsIncr[4],m[1]);
		workingBackStress.AddScaled(-iterationVarsIncr[3]-iterationVarsIncr[4],m[2]);
		workingKappa = fInternal[kkappa] + iterationVars [5];
	    
		/*form new residual */
		residual = 0.0;

	    for (int i = 0; i < kNSD; i++) 
	      {
			for (int j = 0; j < kNSD; j++)
				residual [i] += iterationVars[6] *ElasticConstant(i,j) 
		               *dfdSigmaA(I1, J2, J3, principalEqStress [j], workingKappa);
			residual [i] += iterationVars [i];
	      }

	    for (int i = kNSD; i < 2*kNSD - 1; i++)
	      {
			residual [i] = iterationVars [6] *fCalpha * Galpha(workingBackStress) * dfdDevStressA (I1, J2, J3, principalEqStress [i - kNSD]);
			residual [i] -= iterationVars [i];
	      } 

	    if (fKappaCapped)
	      residual [5] = 0.0;
	    else   
	      residual [5] = iterationVars [6] * KappaHardening(I1, workingKappa) - iterationVars [5];

	    residual [6] = YieldFn(I1, J2, J3, workingKappa);
	} //end while loop


	/* Check that failure surface has not been violated, 
	i.e. J_2^alpha < N, i.e. G^/alpha >= 0 */
	double kin_tol = -1.0e-12;

	if ( Galpha(workingBackStress) < kin_tol)
	  {
	    if (fFossumDebug)
			cout << "FossumSSIsoT::StressPointIteration, Back stress growth limit violated.  Spurious solution obtained.\n" << flush;
	    return false;
	  }

	return true;
}

bool FossumSSIsoT::ResidualIsConverged(dArrayT& residual, dArrayT& residual0)
{
	double yieldFnTol = 1.0e-12, stressTol = 1.0e-10;
	int i;

	for (i = 0; i < 6; i++)
		if (residual0 [i] == 0.0) residual0 [i] = fabs(residual [i]);

	for (i = 0; i < 6; i++)	
		if (fabs(residual [i]) > stressTol && fabs(residual [i]) > stressTol*residual0 [i])
			return false;

	if (fabs(residual [6]) > stressTol && fabs(residual [6]) > stressTol*residual0 [6])
		return false;

	return true;
}

dArrayT FossumSSIsoT::CapKappa(const dArrayT &residual, const LAdMatrixT &dRdX, const double kappa)
{

	dArrayT iterationVarsIncr(7);
	iterationVarsIncr[5] = fInternal[kkappa] - kappa;

	LAdMatrixT modified_dRdX(6); 
	modified_dRdX = 0.0;

	dArrayT modified_residual(6), B(6);

	for (int i=0; i<5; i++)
	{
		modified_residual[i] = residual [i];
		B[i] = dRdX(i,5);
		for (int j=0; j<5; j++) modified_dRdX(i,j) = dRdX(i,j);
		modified_dRdX(i,5) = dRdX(i,6);
		modified_dRdX(5,i) = dRdX(6,i);
	}

	modified_residual[5] = residual[6];
	B[5] = dRdX(6,5);

	modified_residual.AddScaled(-iterationVarsIncr[5],B);

	modified_residual = CondenseAndSolve(modified_dRdX, modified_residual);
  
	for (int i=0; i<5; i++) iterationVarsIncr[i] = modified_residual[i];

	iterationVarsIncr[6] = modified_residual[5];

	return iterationVarsIncr;
}

dArrayT FossumSSIsoT::CondenseAndSolve(const LAdMatrixT& dRdX, const dArrayT& residual)
{
	int condensedMatrixSize = dRdX.Rows() - 1;

	/* break down for static condensation */
	dArrayT U(condensedMatrixSize), V(condensedMatrixSize), Rtilde(condensedMatrixSize);
	LAdMatrixT A(condensedMatrixSize);
  
	for (int i = 0; i < condensedMatrixSize; i++)
	{
		for (int j = 0; j < condensedMatrixSize; j++) A(i,j) = -dRdX(i,j);
		U[i] = -dRdX(i,condensedMatrixSize);	
		V[i] = dRdX(condensedMatrixSize,i);
		Rtilde [i] = residual [i];
	}	

	dArrayT AinvU(condensedMatrixSize), AinvRtilde(condensedMatrixSize),
				iterationVarsIncr(condensedMatrixSize + 1);
	LAdMatrixT Acopy(condensedMatrixSize);
  
	Acopy.Copy(A.Pointer());
	AinvU.Copy(U.Pointer());
	Acopy.LinearSolve(AinvU);
  
	Acopy.Copy(A.Pointer());
	AinvRtilde.Copy(Rtilde.Pointer()); 
	Acopy.LinearSolve(AinvRtilde);

	iterationVarsIncr[condensedMatrixSize] = (residual[condensedMatrixSize] + V.Dot(V, AinvRtilde))/V.Dot(V, AinvU);

	Rtilde.AddScaled(-iterationVarsIncr[condensedMatrixSize], U);
	A.LinearSolve(Rtilde);

	for (int i = 0; i < condensedMatrixSize; i++) 
		iterationVarsIncr [i] = Rtilde[i];

	return iterationVarsIncr;
}


double FossumSSIsoT::ElasticConstant(int i, int j)
{
	if (i==j) 
		return flambda + 2*fmu;
	else 
		return flambda;
}


/* hardening functions and related derivatives*/
double FossumSSIsoT::Galpha(dSymMatrixT alpha)
{
	// if offset fN = 0, already at failure surface, no growth of back stress 
	if (fN == 0.0) return 0.0;
	
	//alpha already deviatoric
	double J2alpha = .5 * alpha.ScalarProduct();

	return 1 - (sqrt(J2alpha))/fN;
}

double FossumSSIsoT::KappaHardening(double I1, double kappa)
{
	return 3 * dfdI1 (I1, kappa) / (dPlasticVolStraindX(kappa) * dXdKappa(kappa));
}

double FossumSSIsoT::dfdDevStressA (double I1, double J2, double J3, double sigmaA)
{
	return dfdJ2(J2, J3) * (sigmaA - I1/3.0) + dfdJ3(J2, J3) * ((sigmaA - I1/3.0) * (sigmaA - I1/3.0) - 2.0 * J2 / 3.0);
}

double FossumSSIsoT::dfdSigmaA(double I1, double J2, double J3, double sigmaA, double kappa)
{
	return dfdI1(I1, kappa) + dfdDevStressA(I1, J2, J3, sigmaA);
}

double FossumSSIsoT::dfdI1(double I1, double kappa)
{
	return -2*YieldFnFfMinusN(I1)*YieldFnFc(I1, kappa)*dFfdI1(I1) - YieldFnFfMinusN(I1) * YieldFnFfMinusN(I1) * dFcdI1(I1, kappa);
}

double FossumSSIsoT::dFfdI1(double I1)
{
	return -fB*fC*exp(fB*I1) - fTheta;
}

double FossumSSIsoT::dFcdI1(double I1, double kappa)
{
	return HeavisideFn(kappa - I1) * -2 * (I1 - kappa) / ((Xfn(kappa) - kappa)*(Xfn(kappa) - kappa));
}

double FossumSSIsoT::dfdJ3(double J2, double J3)
{
	return - YieldFnGamma (J2, J3) * (1 - 1/ fPsi) * 3 * sqrt3 / (2 * sqrt (J2)); 
}

double FossumSSIsoT::dPlasticVolStraindX(double kappa) 
{
	double XLessX0 = Xfn(kappa) - Xfn(fKappa0);
	double returnValue = fW;
	returnValue *= fD1 - 2*fD2*(XLessX0);
	returnValue *= exp (fD1* XLessX0);
	returnValue *= exp ( -fD2 * XLessX0 * XLessX0);

	return returnValue;
}

double FossumSSIsoT::dXdKappa(double kappa)
{
	return  1 - fR * dFfdI1(kappa); 
}

/* Matrix for Stress point NR iteration*/
LAdMatrixT FossumSSIsoT::FormdRdX(double I1, double J2, double J3, dArrayT principalEqStress, double workingKappa, dSymMatrixT workingStress, dSymMatrixT workingBackStress, double dGamma, ArrayT<dSymMatrixT> m)
{
	int A, B, C;
	LAdMatrixT dRdX (7);
	dRdX = 0.0;

	for (A = 0; A < kNSD; A++)
		for (B = 0; B < kNSD; B++)
		{
			for (C = 0; C < kNSD; C++)
				dRdX (A,B) += dGamma * ElasticConstant(A,C) * d2fdSigmaBdSigmaC (I1, J2, J3, principalEqStress[B], principalEqStress[C], B, C, workingKappa); 
			dRdX(A,B) += KroneckerDelta (A, B);

		}

	for (A = 0; A < kNSD; A++)
		for (B = 0; B < kNSD - 1; B++)
			for (C = 0; C < kNSD; C++)
				dRdX (A, B + kNSD) += dGamma * ElasticConstant(A,C) * (- d2fdSigmaBdSigmaC (I1, J2, J3, principalEqStress[B], principalEqStress[C], B, C, workingKappa) 
						+ d2fdSigmaBdSigmaC (I1, J2, J3, principalEqStress[C], principalEqStress[2], C, 2, workingKappa));
 
	for (A = 0; A < kNSD; A++)
		for (C = 0; C < kNSD; C++)
		{
			dRdX (A,5) += dGamma * ElasticConstant(A,C) * d2fdSigmaCdKappa (I1, workingKappa); 
			dRdX (A, 6) += ElasticConstant(A,C)*dfdSigmaA(I1, J2, J3, principalEqStress [C], workingKappa);
		}
      
	/* ------ */
	//dR(alpha A)/d(dSigma B)
 
	for (A = 0; A < kNSD - 1; A++)
		for (B = 0; B < kNSD; B++)
		{ 
			dRdX (A + kNSD, B) = Galpha(workingBackStress) * d2fdDevStressdSigmaB (I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B); 
			dRdX (A + kNSD, B) *= fCalpha * dGamma; 
		}

	//dR(alpha A)/d(dAlpha B)
	for (A = 0; A < kNSD - 1; A++)
		for (B = 0; B < kNSD -1 ; B++)
		{ 

			dRdX (A + kNSD, B + kNSD) = - Galpha(workingBackStress) * d2fdDevStressdSigmaB (I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B); 
			dRdX (A + kNSD, B + kNSD) +=  Galpha(workingBackStress) * d2fdDevStressdSigmaB (I1, J2, J3, principalEqStress[A], principalEqStress[2], A, 2);
			dRdX (A + kNSD, B + kNSD) +=  dGalphadAlphaB (workingBackStress, principalEqStress, B, m) * dfdDevStressA(I1, J2, J3, principalEqStress[A]);
			dRdX (A + kNSD, B + kNSD) *= fCalpha * dGamma;
			dRdX (A + kNSD, B + kNSD) -= KroneckerDelta (A, B);
		}

	//dR(Alpha A)/d(Kappa) = 0

	//dR(alpha A)/d(dGamma B)
	for (A = 0; A < kNSD - 1; A++)
	{
		dRdX (A + kNSD, 6) = fCalpha * Galpha(workingBackStress) *
				dfdDevStressA(I1, J2, J3, principalEqStress[A]);
	}
	// ------------ 

	for ( A = 0; A < kNSD; A++)
		dRdX (5, A) = 3 * dGamma * d2fdI1dI1(I1, workingKappa) 
				/ (dPlasticVolStraindX( workingKappa) * dXdKappa (workingKappa));
                 
	double c1 = dPlasticVolStraindX(workingKappa);
	double c2 = dXdKappa(workingKappa); 
	double c3 = d2fdI1dKappa (I1, workingKappa);
	double c4 =  dfdI1(I1, workingKappa);
	double c5 = d2PlasticVolStraindXdX(workingKappa);
	double c6 = d2XdKappadKappa(workingKappa);

	dRdX (5, 5) = dPlasticVolStraindX(workingKappa) * dXdKappa(workingKappa) * d2fdI1dKappa (I1, workingKappa);
	dRdX (5, 5) -= dfdI1(I1, workingKappa) * dPlasticVolStraindX(workingKappa) * d2XdKappadKappa(workingKappa);
	dRdX (5, 5) -= dfdI1(I1, workingKappa) * dXdKappa(workingKappa) * dXdKappa(workingKappa) * d2PlasticVolStraindXdX(workingKappa);
	dRdX (5, 5) *= 3 * dGamma / (dPlasticVolStraindX(workingKappa) * dPlasticVolStraindX(workingKappa) * dXdKappa (workingKappa) * dXdKappa (workingKappa));
 	dRdX (5, 5) -= 1;
 
	dRdX (5, 6) = KappaHardening(I1, workingKappa);

	for ( A = 0; A < kNSD; A++)
		dRdX (6, A) = dfdSigmaA(I1, J2, J3, principalEqStress[A], workingKappa);

 
	for ( A = 0; A < kNSD - 1; A++)
		dRdX (6, A + kNSD) = -dfdSigmaA (I1, J2, J3, principalEqStress[A], workingKappa)
				+ dfdSigmaA (I1, J2, J3, principalEqStress[2], workingKappa);
 
	dRdX (6, 5) = dfdKappa(I1, workingKappa);

	return dRdX;
}

/* derivative for FromdRdX */
int FossumSSIsoT::KroneckerDelta (int A, int B)
{
	if ( A == B ) 
		return 1;
	else 
		return 0;
} 

double FossumSSIsoT::d2fdSigmaBdSigmaC (double I1, double J2, double J3, double principalEqStressA, double principalEqStressB, int A, int B, double kappa)
{
	return d2fdI1dI1(I1, kappa) + d2fdDevStressdSigmaB(I1, J2, J3, principalEqStressA, principalEqStressB, A, B);
}

double FossumSSIsoT::d2fdDevStressdSigmaB (double I1, double J2, double J3, double principalEqStressA, double principalEqStressB, int A, int B)
{
	double returnValue = 0.0, xiA = principalEqStressA - I1/3.0, xiB = principalEqStressB - I1/3.0;
  
	returnValue += dfdJ2(J2, J3) * (KroneckerDelta(A,B) - 1.0/3.0);
	returnValue += dfdJ3(J2, J3) * (2 * xiA * KroneckerDelta (A, B)  - 2.0/3.0 *(xiA + xiB));
	returnValue += d2fdJ2dJ2 (J2, J3) * xiA * xiB;
	returnValue += d2fdJ2dJ3 (J2, J3) * (xiA * (xiB*xiB - 2.0/3.0 * J2) + xiB * (xiA*xiA - 2.0/3.0 * J2)); 
	returnValue += d2fdJ3dJ3 (J2, J3) * (xiA * xiA  - 2.0/3.0 * J2) * (xiB * xiB  - 2.0/3.0 * J2);

	return returnValue;
}


double FossumSSIsoT::d2fdI1dI1(double I1, double kappa)
{
	double Fc = YieldFnFc(I1, kappa);
	double FfMinusN = YieldFnFfMinusN(I1);
	double dFfbydI1 = dFfdI1(I1);

	return -1*(2 * Fc * dFfbydI1 * dFfbydI1 + 4 * FfMinusN * dFfbydI1 * dFcdI1(I1, kappa)
          + 2 * FfMinusN * Fc * d2FfdI1dI1(I1) + FfMinusN*FfMinusN*d2FcdI1dI1(I1, kappa));
}

double FossumSSIsoT::d2FfdI1dI1(double I1)
{
	return - fB * fB * fC * exp (fB * I1);
}

double FossumSSIsoT::d2FcdI1dI1(double I1, double kappa)
{
	double XMinusKappa = Xfn (kappa) - kappa;

	return HeavisideFn( kappa - I1) * -2 / (XMinusKappa * XMinusKappa);
}

double FossumSSIsoT::d2fdJ2dJ2 (double J2, double J3)
{
	double gamma = YieldFnGamma(J2, J3);
	double dGdJ2 = dGammadJ2(J2, J3);

	return 2*(2 * gamma * dGdJ2 + J2 * dGdJ2 * dGdJ2 + J2 * gamma * d2GammadJ2dJ2(J2, J3));
}

double FossumSSIsoT::d2GammadJ2dJ2(double J2, double J3)
{
	return (1 - 1/fPsi) * -45*sqrt3 * J3/ (16 * J2*J2*J2*sqrt(J2));
}

double FossumSSIsoT::d2fdJ2dJ3 (double J2, double J3)
{
	double gamma = YieldFnGamma(J2, J3);
	double dGdJ3 = dGammadJ3 (J2);

	return 2 * (gamma * dGdJ3 + J2 * dGammadJ2(J2, J3) * dGdJ3 + J2 * gamma * d2GammadJ2dJ3(J2) );
}

double FossumSSIsoT::d2GammadJ2dJ3 (double J2)
{
	return ( 1 - 1/fPsi) * 9 * sqrt3 / (8 * J2 * J2 * sqrt(J2));
}

double FossumSSIsoT::dGammadJ3(double J2)
{
	return ( 1 - 1/fPsi) * -3 * sqrt3 / (4 * J2 * sqrt(J2));
}

double FossumSSIsoT::d2fdJ3dJ3 (double J2, double J3)
{
	double dGdJ3 = dGammadJ3(J2);

	return 2 * J2 * dGdJ3 * dGdJ3; 
}

double FossumSSIsoT::d2fdSigmaCdKappa (double I1, double kappa)
{
	double FfMinusN = YieldFnFfMinusN(I1);
  
	double c1 = dFfdI1(I1);
	double c2 = dFcdKappa(I1, kappa);
	double c3 = d2FcdI1dKappa(I1, kappa);

	return -2 * FfMinusN * dFfdI1(I1) * dFcdKappa(I1, kappa)
			- FfMinusN * FfMinusN * d2FcdI1dKappa(I1, kappa); 
}

double FossumSSIsoT::d2FcdI1dKappa(double I1, double kappa)
{
	double XminusKappa = Xfn(kappa) - kappa;

	return 2 * HeavisideFn(kappa - I1) * kappa 
		* (XminusKappa + 2 * fR * (I1 - kappa) * (fB * fC * exp (fB * kappa) + fTheta))
		/ (XminusKappa * XminusKappa * XminusKappa);
}

double FossumSSIsoT::dGalphadAlphaB (dSymMatrixT alpha, dArrayT principalEqStress, int B, ArrayT<dSymMatrixT> m)
{
	/* N = 0 => Galpha is identically 0, does not change with stress or alpha */
	if (fN == 0) return 0;

	dSymMatrixT nB(3);
	nB.DiffOf(m[B], m[2]);

	double J2alpha = .5 * alpha.ScalarProduct();

	if (J2alpha == 0.0) return -1.0/(sqrt(2.0) * fN);

	return -1.0/(2*fN*sqrt(J2alpha))*InnerProduct(alpha, nB);
}

double FossumSSIsoT::d2fdI1dKappa (double I1, double kappa)
{
	return d2fdSigmaCdKappa(I1, kappa);
}

double FossumSSIsoT::dFcdKappa(double I1, double kappa)
{
	double XminusKappa = Xfn(kappa) - kappa;

	return 2 * HeavisideFn(kappa - I1) * kappa * ( I1 - kappa) * ((XminusKappa) + fR * (I1 - kappa) * ( fTheta + fB * fC * exp (fB * kappa))) / (XminusKappa * XminusKappa * XminusKappa);
}

double FossumSSIsoT::d2XdKappadKappa( double kappa)
{
	return fR * fB * fB * fC * exp(fB * kappa);
}

double FossumSSIsoT::d2PlasticVolStraindXdX(double kappa)
{
	double XminusX0 = Xfn (kappa) - Xfn (fKappa0);
	double work = fD1 - 2*fD2*XminusX0;

	return fW * ((-2*fD2 + work * work) * exp ( fD1 * XminusX0 - fD2 * XminusX0 * XminusX0));  
}

double FossumSSIsoT::dfdKappa(double I1, double kappa)
{
	double FfMinusN =YieldFnFf(I1) - fN;

	return -FfMinusN*FfMinusN*dFcdKappa(I1, kappa);
}

double FossumSSIsoT::InnerProduct(dSymMatrixT A, dSymMatrixT B)
{
	double returnValue = 0.0;
	int i, j;
 
	for (i = 0; i < kNSD; i++) returnValue += A[i]*B[i];
 
	for (i = kNSD; i < 2 * kNSD; i++) returnValue += 2*A[i]*B[i];

	return returnValue;
}

/*---------------------------------------------------------*/

/* consistent tangent modulus */
const dMatrixT& FossumSSIsoT::c_ijkl(void)
{
	return SSSolidMatT::c_ijkl();
}


/* consistent perfectly plastic modulus */
const dMatrixT& FossumSSIsoT::c_perfplas_ijkl(void)
{
	dMatrixT elasticModulus(6), elasticCompliance(6);
	dMatrixT generalizedModulus(14), generalizedCompliance(14);
	dMatrixT d2fdSigmadSigma(6), d2fdqdq(7), dhdq(7);
	dMatrixT d2fdSigmadq(6,7), dhdSigma(7,6);
	dMatrixT Ce = HookeanMatT::Modulus();

	ElementCardT& element = CurrentElement();
	int ip = CurrIP();


	if (!element.IsAllocated() || (element.IntegerData())[ip] != kIsPlastic)
	{
		fModulusPerfPlas = HookeanMatT::Modulus();
		return fModulusPerfPlas;
	}
	else
	{         
		 /* load internal state variables */
		LoadData(element,ip);
		
		double kappa = fInternal[kkappa] + fInternal[kdeltakappa];
		dSymMatrixT alpha(3); 
		alpha.SumOf(fBackStress, fDeltaAlpha);

		double I1 = 0.0, J2 =0.0, J3 = 1.0;

		for (int i = 0; i < kNSD; i++)
			I1 += principalEqStress[i];

		for (int i = 0; i < kNSD; i++)
		{
			J2 += 0.5*(principalEqStress[i] - I1/3.0) * (principalEqStress[i] - I1/3.0);
			J3 *= (principalEqStress[i] - I1/3.0);
		}

		/*form generalized modulus */
		elasticCompliance.Inverse(Ce);
		generalizedCompliance = 0.0;
       
		d2fdSigmadSigma = D2fdSigmadSigma(I1, J2, J3, kappa, principalEqStress, m);

		for (int i=0; i<6; i++)
			for (int j=0; j<6; j++)
				generalizedCompliance(i,j) = elasticCompliance(i,j) + fInternal[kdgamma]*d2fdSigmadSigma(i,j);

		for (int i=6; i<13; i++)
			for (int j=6; j<13; j++)
				generalizedCompliance(i,j) = KroneckerDelta(i,j);

		dArrayT df(13), dr(13), hardeningFns(7), dfdq(7);
		dSymMatrixT dfdSigma(3), dfdAlpha(3);

		dfdSigma = DfdSigma(I1,J2,J3, kappa, principalEqStress, m);

		for (int i=0; i<6; i++)
		{
			generalizedCompliance(13, i) = dfdSigma[i];
			generalizedCompliance(i, 13) = dfdSigma[i];
		}

		generalizedCompliance(13,12) = 0.0;

		//double shear 
		for (int i=3; i<6; i++)
		{
			generalizedCompliance(i + 6, i + 6) -= 1.0;
			generalizedCompliance(13, i) *= 2.0;
			generalizedCompliance(13,i+6) *= 2.0;
			generalizedCompliance(i, 13) *= 2.0;
			generalizedCompliance(i+6, 13) *= 2.0;
			for (int j = 0; j < 6; j++)
			{
				generalizedCompliance(i, j) -= elasticCompliance(i,j);
				generalizedCompliance(i, j) *= 2.0;
				generalizedCompliance(i, j) += elasticCompliance(i,j);

				generalizedCompliance(j, i) -= elasticCompliance(j,i);
				generalizedCompliance(j, i) *= 2.0;
				generalizedCompliance(j, i) += elasticCompliance(j,i);

				generalizedCompliance(i, j + 6) *= 2.0;
				generalizedCompliance(j, i + 6) *= 2.0;

				generalizedCompliance(i + 6, j) *= 2.0;
				generalizedCompliance(j + 6, i) *= 2.0;

				generalizedCompliance(i + 6, j + 6) *= 2.0;
				generalizedCompliance(j + 6, i + 6) *= 2.0;
			}
			generalizedCompliance(i + 6, i + 6) += 2.0;
		}

		generalizedModulus = generalizedCompliance.Inverse();

		/* pull out modulus */
		for (int i=0; i<6; i++)
			for (int j=0; j<6; j++)
				fModulusPerfPlas(i,j) = generalizedModulus(i,j);
      
		dMatrixT fModulus2PerfPlas(6);   
		fModulus2PerfPlas = Ce;

		double chi = dfdSigma.B_ij_A_ijkl_B_kl(Ce);
		dSymMatrixT CeTimesdfdSigma(3);

		CeTimesdfdSigma.A_ijkl_B_kl(Ce, dfdSigma);
		dMatrixT corr(6);

		corr.Outer(CeTimesdfdSigma, CeTimesdfdSigma);
		fModulus2PerfPlas.AddScaled(-1.0/chi, corr);
	} //end else
 
	return fModulusPerfPlas;
}

/* tensor-valued derivatives for consistent tangent */
dMatrixT FossumSSIsoT::D2fdSigmadSigma(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m)
{
	int A,B,i,j;
	double d2;
	dMatrixT d2fdSigmadSigma(6);

	/*initialize */
	d2fdSigmadSigma = 0.0;

	for (A=0; A<3; A++)
		for (B=0; B<3; B++)//for (B=A; B<3; B++)
		{
			d2 = d2fdSigmaBdSigmaC(I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B, kappa);
			for (i=0; i<6; i++)
				for (j=0; j<6; j++)
				{
					d2fdSigmadSigma(i,j) += d2 * (m[A]) [i] * (m[B]) [j]; 
				}
		}
	return d2fdSigmadSigma;
}

dSymMatrixT FossumSSIsoT::DfdSigma(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m)
{
	int A, i;
	double dfdA;
	dSymMatrixT dfdSigma(3);

	dfdSigma = 0.0;

	for (A=0; A<3; A++)
	{
		dfdA = dfdSigmaA(I1, J2, J3, principalEqStress[A], kappa);
		for (i=0; i<6; i++)//for (i=0; i<3; i++)
		//for (j=i; j<3; j++)
			dfdSigma[i] += dfdA *(m[A])[i]; 
	}

	return dfdSigma;
}

dSymMatrixT FossumSSIsoT::DfdAlpha(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m)
{
	int A, i;
	double dfdA;
	dSymMatrixT dfdAlpha(3);
	ArrayT<dSymMatrixT> n(2);  

	for (int i=0; i<6; i++)
		dfdAlpha[i] = - DfdSigma(I1,J2,J3, kappa, principalEqStress, m)[i];

	return dfdAlpha;
}

/*------*/

dArrayT FossumSSIsoT::Hardening(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m, dSymMatrixT alpha) 
{
	int i;
	dArrayT hardening(7);
	dSymMatrixT dfdDevStress(3);
  
	dfdDevStress = 0.0;
  
	for (int A = 0; A < kNSD; A++)
		dfdDevStress.AddScaled(dfdDevStressA(I1, J2, J3, principalEqStress[A]), m[A]) ;
  
	for (i=0; i<6; i++)
		hardening [i] = fCalpha * Galpha(alpha) * dfdDevStress [i];
   
	if (!fKappaCapped)
		hardening [6] = KappaHardening(I1, kappa);

	return hardening;
}
