/* 3-invariant, single-surface dilation/compaction plasticity model
 * with isotropic and kinematic hardeneing
 * Implemented 8/02 Craig Foster
 */

#include "FossumSSIsoT.h"
#include "SSMatSupportT.h"

//#include <iostream.h>
//#include "ifstreamT.h"
#include "dSymMatrixT.h"

#include "ElementCardT.h"
#include "StringT.h"
#include "DetCheckT.h"

#include <math.h>

#include "iArrayT.h"
#include "ArrayT.h"

using namespace Tahoe;

const int kNumInternal = 4; // number of internal state variables
const double sqrt3 = sqrt(3.0);
//const double sqrt32 = sqrt(3.0/2.0);
const double kYieldTol = 1.0e-10;
const int kNSD = 3;

/* element output data */
const int kNumOutput = 11;
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
	"J3",
	"loccheck"

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
	fPsi(1.0),
	fN(-1.0),
	fNumIP(NumIP()),
	fmu(Mu()),
	flambda(Lambda()),
//	SSSolidMatT(in, support),
//	IsotropicT(in),


	HookeanMatT(3),
	fStress(3),
	fModulus(dSymMatrixT::NumValues(3)),
	fModulusPerfPlas(dSymMatrixT::NumValues(3)),   
	fModulusContinuum(dSymMatrixT::NumValues(3)),
	fModulusContinuumPerfPlas(dSymMatrixT::NumValues(3)),     
	fElasticStrain(kNSD),
	fStressCorr(kNSD),
	fModuliCorr(dSymMatrixT::NumValues(kNSD)),
	fModuliCorrPerfPlas(dSymMatrixT::NumValues(kNSD)),
	fDevStress(kNSD),
	fMeanStress(0.0),
	fDevStrain(kNSD), 
	fTensorTemp(dSymMatrixT::NumValues(kNSD)),
	IdentityTensor2(kNSD),
	One(kNSD),
	fBackStress(kNSD),
	fDeltaAlpha(kNSD),

	/*spectral decomp parameters */
	spectre(kNSD),
	m(kNSD),
	principalEqStress(kNSD),

	/*viscous paramters*/
	//fStressInviscid(kNSD),
        fTimeFactor(1.0)

{
	/* allocate space for principal dirs m */
	for (int i = 0; i < 3; i++)
	{
		m[i].Dimension(kNSD);
	}

	/* read parameters -old tahoe 1_x style*/
	/*
	in >> fA;      
	if (fA < 0.0 ) {cout << "Bad value for A\n" << flush; throw ExceptionT::kBadInputValue;}   
	in >> fB;
	in >> fC;
	in >> fTheta; 
	if (fTheta < 0.0 ) {cout << "Bad value for theta\n" << flush; throw ExceptionT::kBadInputValue;}
	in >> fR;     
	if (fR < 0.0 ) {cout << "Bad value for R\n" << flush; throw ExceptionT::kBadInputValue;}  
	in >> fKappa0;
	in >> fW;
	in >> fD1; 
	in >> fD2;
	in >> fCalpha; 
	if (fCalpha < 0.0 ) {cout<< "Bad value for Calpha, Calpha = "<< fCalpha << "\n"<< flush; throw ExceptionT::kBadInputValue;}
	in >> fPsi;     
	if (fPsi < 0.0 ) {cout << "Bad value for Psi\n" << flush; throw ExceptionT::kBadInputValue;}
	in >> fN;
	if (fN < 0.0) {cout << "Bad value for N\n" << flush; throw ExceptionT::kBadInputValue;}
	in >> fFluidity;
        if (fFluidity < 0.0) {cout << "Bad value for fluidity\n" << flush; throw ExceptionT::kBadInputValue;}
	in >> fFossumDebug;
        if (fFossumDebug != 0 && fFossumDebug != 1) {cout << "Bad value for fluidity\n" << flush; throw ExceptionT::kBadInputValue;}
	*/


	 

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
 
  ParameterT A(fA, "Shear_surface_parameter_A");    
  A.AddLimit(0.0, LimitT::LowerInclusive); 
  list.AddParameter(A);

  list.AddParameter(fB, "Shear_surface_parameter_B"); 
  list.AddParameter(fC, "Shear_surface_parameter_C");
  
  ParameterT theta(fTheta, "Shear_surface_parameter_theta");    
  theta.AddLimit(0.0, LimitT::LowerInclusive); 
  list.AddParameter(theta);

  ParameterT R(fR, "Cap_Ratio_R");    
  R.AddLimit(0.0, LimitT::LowerInclusive); 
  list.AddParameter(R);

  list.AddParameter(fW, "Cap_growth_Parameter_W"); 
  list.AddParameter(fD1, "Cap_growth_Parameter_D1");
  list.AddParameter(fD2, "Cap_growth_Parameter_D2");

  ParameterT calpha(fCalpha, "shear_surface_growth_factor_c_alpha");
  calpha.AddLimit(0.0, LimitT::LowerInclusive);
  list.AddParameter(calpha);     

  ParameterT N(fN, "shear_surface_offset_N");
  N.AddLimit(0.0, LimitT::LowerInclusive);
  list.AddParameter(N);

  ParameterT fluidity(fFluidity, "fluidity_parameter_eta");
  fluidity.AddLimit(0.0, LimitT::LowerInclusive);
  list.AddParameter(fluidity);

  ParameterT debug(fFossumDebug, "local_debug_parameter");
  debug.SetDefault(false);
  list.AddParameter(debug);
}

/* accept parameter list */
void FossumSSIsoT::TakeParameterList(const ParameterListT& list)
{
    /* inherited */
    SSIsotropicMatT::TakeParameterList(list);
    HookeanMatT::TakeParameterList(list);

    fA = list.GetParameter("Shear_surface_parameter_A");
    fB = list.GetParameter("Shear_surface_parameter_B");
    fC = list.GetParameter("Shear_surface_parameter_C");
    fTheta = list.GetParameter("Shear_surface_parameter_theta");
    fR = list.GetParameter("Cap_Ratio_R");
    fW = list.GetParameter("Cap_growth_Parameter_W");
    fD1 = list.GetParameter("Cap_growth_Parameter_D1");
    fD2 = list.GetParameter("Cap_growth_Parameter_D2");
    fCalpha = list.GetParameter("shear_surface_growth_factor_c_alpha");
    fN = list.GetParameter("shear_surface_offset_N");
    fFluidity = list.GetParameter("fluidity_parameter_eta");
    fFossumDebug = list.GetParameter("local_debug_parameter");
}

void FossumSSIsoT::DefineSubs(SubListT& sub_list) const
{
    /* inherited */
    SSIsotropicMatT::DefineSubs(sub_list);
    HookeanMatT::DefineSubs(sub_list);

    /* parameters for Drucker-Prager plasticity */
    //sub_list.AddSub("Fossum_small_strain");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */


ParameterInterfaceT* FossumSSIsoT::NewSub(const StringT& name) const
{
  //if (name == "Fossum_small_strain_2D")
  //  return new FossumSSIso2DT;
  //else
  //  {
      /* inherited */
      ParameterInterfaceT* params = SSIsotropicMatT::NewSub(name);
      if (params) 
	return params;
      else
	return HookeanMatT::NewSub(name);
      //      }
}


#if 0
/* write parameters to stream */
void FossumSSIsoT::Print(ostream& out) const
{
	SSSolidMatT::Print(out);
	IsotropicT::Print(out);

	out << "Material parameter for F_f A.............. = " << fA << endl; 
	out << "Material parameter for F_f B.............. = " << fB << endl;
	out << "Material parameter for F_f C.............. = " << fC << endl;
	out << "Material parameter for F_f theta.......... = " << fTheta << endl;
	out << "Ratio of radii of cap function R.......... = " << fR << endl;        
	out << "Intial cap locator kappa0................. = " << fKappa0 << endl;
	out << "Cap hardening parameter W................. = " << fW << endl;
	out << "Cap hardening parameter D1................ = " << fD1 << endl;
	out << "Cap hardening parameter D2................ = " << fD2 << endl;
	out << "Back stress growth rate factor C_alpha.... = " << fCalpha << endl;
	out << "Ratio of Tensile to Compressive Str psi... = " << fPsi << endl;
	out << "Offset from yield sfce to failure sfce N.. = " << fN << endl; 
        out << "Fluidity parameter eta.................... = " << fFluidity <<endl;
	if (fFossumDebug == 0)
	  out << "Fossum model debugging is off" << endl;
	else
	  out << "Fossum model debugging is on" << endl;
 
}

void FossumSSIsoT::PrintName(ostream& out) const
{
	SSSolidMatT::PrintName(out);
	out << "Geomaterial plasticity model: 3-invariant, single-surface dilation\\\n";
	out << "compaction plasticity model w/ isotropic and kinematic hardening\n"; 
	out << "    Kirchhoff-St.Venant\n";
}
#endif
        
/*  protected: */

/*
 * Returns the value of the yield function given the
 * Cauchy stress vector and state variables, where alpha is
 * the back stress and kappa an internal state variable
 */
double FossumSSIsoT::YieldCondition(const dSymMatrixT& stress, 
									const double kappa, dSymMatrixT& backStress)
{
	// not ready yet
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

/* initialization */
void FossumSSIsoT::Initialize(void)
{
ExceptionT::GeneralFail("FossumSSIsoT::Initialize", "out of date");
#if 0
	/* inherited */
	HookeanMatT::Initialize();
#endif
}

/* form of tangent matrix (symmetric by default) */
GlobalT::SystemTypeT FossumSSIsoT::TangentType(void) const { return GlobalT::kSymmetric; }

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

/*
* Test for localization using "current" values for Cauchy
* stress and the spatial tangent moduli. Returns 1 if the
* determinant of the acoustic tensor is negative and returns
* the normal for which the determinant is minimum. Returns 0
* of the determinant is positive.
*/

int FossumSSIsoT::IsLocalized(dArrayT& normal)
{
	DetCheckT checker(fStress, fModulus, fModulusCe);
	//checker.SetElementGroup(ContinuumElement());
	checker.SetfStructuralMatSupport(*fSSMatSupport);

	int loccheck = checker.IsLocalized(normal);
	return loccheck;
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
	int i, ip = CurrIP();
	LoadData(element, ip);
	dMatrixT Ce = HookeanMatT::Modulus();

	/*OUTPUT FOR ALPHA, KAPPA */ 
	if (element.IsAllocated())
	{
		for (i = 0; i < 6 ; i++) output [i] = fBackStress [i];
		output [6] = fInternal[kkappa];
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

	if (element.IsAllocated())
	//if (0)  //to disable localization check
	{
		const iArrayT& flags = element.IntegerData();
		if (flags[CurrIP()] == kIsPlastic)
		{
			// check for localization
			// compute modulus 

			//const dMatrixT& modulus = c_ijkl();

			//const dMatrixT& modulus = c_perfplas_ijkl();
			//const dMatrixT& modulus = con_ijkl();
			const dMatrixT& modulus = con_perfplas_ijkl();

			/* localization condition checker */
			DetCheckT checker(stress, modulus, Ce);

			AutoArrayT <dArrayT> normals;
			AutoArrayT <dArrayT> slipdirs;
			normals.Dimension(3);
			slipdirs.Dimension(3);
			output[10] = checker.IsLocalized_SS(normals,slipdirs);

		}
		else
		{
			output[10] = 0.0;
		}
	}
	else
	{
		output[10] = 0.0;
	}
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

/* return correction to stress vector computed by mapping the
 * stress back to the yield surface, if needed */
const dSymMatrixT& FossumSSIsoT::StressCorrection(
		const dSymMatrixT& trialstrain, ElementCardT& element, int ip)
{
	/*not ready yet*/

	/*initialize*/      
	fStressCorr = 0.0;
	return fStressCorr;
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
		
		// is this called after Update???
		//double kappa = fInternal[kkappa];
		double kappa = fInternal[kkappa] + fDeltaKappa;
		dSymMatrixT alpha(3); 
		//alpha = fBackStress;
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
		//cout << "dfdSigma = \n" << dfdSigma << endl;
		//cout << "dfdq = \n" << dfdq << endl;
		//cout << "hardeningFns = \n" << hardeningFns << endl;

		double chi = dfdSigma.B_ij_A_ijkl_B_kl(Ce) - dfdq.Dot(dfdq,hardeningFns);
		dSymMatrixT CeTimesdfdSigma(3);
		CeTimesdfdSigma.A_ijkl_B_kl(Ce, dfdSigma);
		dMatrixT corr(6);
		corr.Outer(CeTimesdfdSigma, CeTimesdfdSigma);
  
		//cout << "chi = " << chi << endl;
		//cout << "corr = \n" << corr << endl;

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
		
		// is this called after Update???
		//double kappa = fInternal[kkappa];
		double kappa = fInternal[kkappa] + fDeltaKappa;

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
		//cout << "dfdSigma = \n" << dfdSigma << endl;

		double chi = dfdSigma.B_ij_A_ijkl_B_kl(Ce);
		dSymMatrixT CeTimesdfdSigma(3);
		CeTimesdfdSigma.A_ijkl_B_kl(Ce, dfdSigma);
		dMatrixT corr(6);
		corr.Outer(CeTimesdfdSigma, CeTimesdfdSigma);
  
		//cout << "chi = " << chi << endl;
		//cout << "corr = \n" << corr << endl;

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

	int d_size = 0;
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fPlasticStrain
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; //fUnitNorm
	d_size += dSymMatrixT::NumValues(kNSD)*fNumIP; // fBackStress
	d_size += kNumInternal*fNumIP;        //fInternal

	/* construct new plastic element */
	element.Dimension(i_size, d_size);
        
	/* initialize values */
	element.IntegerData() = kIsElastic;
	element.DoubleData()  = 0.0;  // initialize all double types to 0.0

	/* initialize kappa to kappa0 */
	int i;

	for (i=0; i<fNumIP; i++)
		element.DoubleData() [3*dSymMatrixT::NumValues(kNSD)*fNumIP + i*kNumInternal] = fKappa0;
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
			//cout << "kdgamma = " << fInternal[kdgamma] << endl;

			double I1, J2, J3;
			I1 = 0.0; J2 = 0.0; J3 = 1.0;

			for (int i = 0; i < kNSD; i++) I1 += principalEqStress[i];
			
			for (int i = 0; i < kNSD; i++)
			{
				J2 += 0.5 * (principalEqStress[i] - I1/3.0) * (principalEqStress[i] - I1/3.0);
				J3 *= (principalEqStress[i] - I1/3.0);
			}
			/*
			fInternal[kkappa] += dgamma*KappaHardening(I1, fInternal[kkappa]);
			fBackStress.AddScaled(dgamma*fCalpha*Galpha(fStress, J2), 
						DfdDevStress(I1, J2, J3, principalEqStress, m));
			*/

			fInternal[kkappa] += fDeltaKappa;
			fBackStress += fDeltaAlpha;

			/* plastic strain increment */
			fPlasticStrain.AddScaled(dgamma, DfdSigma(I1, J2, J3, fInternal[kkappa], principalEqStress, m));

			fKappa0 = fInternal[kkappa];
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
	fUnitNorm.Alias     (         dim, &d_array[  offset + dex]);   
	fBackStress.Alias   (         dim, &d_array[2*offset + dex]);  
	fInternal.Alias     (kNumInternal, &d_array[3*offset + ip*kNumInternal]);
}

/* returns 1 if the trial elastic strain state lies outside of the 
 * yield surface */
int FossumSSIsoT::PlasticLoading(const dSymMatrixT& trialstrain, 
	ElementCardT& element, int ip)
{
	/* not yet plastic */
	//if (!element.IsAllocated()) 
	//return( YieldCondition(DeviatoricStress(trialstrain,element),
	//		   MeanStress(trialstrain,element), 0.0) > kYieldTol );
	/* already plastic */
	//else
	if (element.IsAllocated()) //added 
	{
		/* get flags */
		iArrayT& Flags = element.IntegerData();
                
		/* load internal variables */
		LoadData(element, ip);
                
		//fInternal[kftrial] = YieldCondition(DeviatoricStress(trialstrain,element), MeanStress(trialstrain,element),fInternal[kalpha]);

		/* plastic */
		if (fInternal[kftrial] > kYieldTol)
		{                
			/* compute unit normal */
			double& norm = fInternal[kstressnorm];

			norm = sqrt(fDevStress.ScalarProduct());
			fUnitNorm.SetToScaled(1.0/norm, fDevStress);
                
			/* set flag */
			Flags[ip] = kIsPlastic;
        
			return 1;
		}
		/* elastic */
		else
		{
			/* set flag */
			Flags[ip] = kIsElastic; //removed to avoid restting 7/01            
			return 0;
		}
	}
	
	//if element is not allocated, assume elastic
	//assumed elastic by default
	return 0;
}       

/* Computes the stress corresponding to the given element
 * and elastic strain.  The function returns a reference to the
 * stress in fDevStress */
dSymMatrixT& FossumSSIsoT::DeviatoricStress(const dSymMatrixT& trialstrain,
											const ElementCardT& element)
{
	#pragma unused(element)

	/* deviatoric strain */
	fDevStrain.Deviatoric(trialstrain);

	/* compute deviatoric elastic stress */
	fDevStress.SetToScaled(2.0*fmu,fDevStrain);

	return fDevStress;
}

/* computes the hydrostatic (mean) stress */
double FossumSSIsoT::MeanStress(const dSymMatrixT& trialstrain,
								const ElementCardT& element)
{
	#pragma unused(element)

	//fMeanStress = fkappa*trialstrain.Trace();
	fMeanStress = 0.0;
	return fMeanStress;
}

/*---------------------------------------------------------------*/

double FossumSSIsoT::YieldFnFfMinusN(double I1)
{
	return fA - fC*exp(fB*I1) - fTheta * I1 - fN;
}

double FossumSSIsoT::YieldFnFc(double I1, const double kappa)
{
	return 1 - HeavisideFn(Lfn(kappa) - I1)*(I1 - Lfn(kappa))*(I1 - Lfn(kappa))/((Xfn(kappa) - Lfn(kappa))*(Xfn(kappa) - Lfn(kappa)));
}

int FossumSSIsoT::HeavisideFn(double arg)
{
	if (arg < 0)
		return 0;
	else
		return 1;
}

double FossumSSIsoT::Lfn(const double kappa)
{
	//if (kappa < fKappa0)
	return kappa;
	//else
	//return fKappa0;
	//return HeavisideFn(fKappa0 - kappa) * kappa;
}

double FossumSSIsoT::Xfn(const double kappa)
{
	double L = Lfn(kappa);

	return L - fR * YieldFnFf(L);
}

double FossumSSIsoT::YieldFnFf(double I1)
{
	return fA - fC*exp(fB*I1) - fTheta * I1;
}

/*-------------------------------------------------------------*/
/*Return Mapping algorithm */

/* stress */
const dSymMatrixT& FossumSSIsoT::s_ij(void)
{
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	const dSymMatrixT& e_tot = e();
	const dSymMatrixT& e_els = ElasticStrain(e_tot, element, ip);
	double yieldFnTol = 1.0e-8;

	/* elastic stress */
	HookeanStress(e_els, fStress);

	/* modify Cauchy stress (return mapping) */
	//fStress += StressCorrection(e_els, element, ip);
	//return fStress; 

	if (!element.IsAllocated()) AllocateElement(element);

	LoadData(element, ip);

	dSymMatrixT backStress_n; //load alpha_n
	double kappa_n; //load kappa  
	backStress_n = fBackStress;
	kappa_n = fInternal[kkappa];

	/*
	cout << "fStress =\n" << fStress << endl;
	cout << "fBackStress =\n" << fBackStress << endl;
	cout << "fInternal[kkappa] = " << fInternal[kkappa] << endl << endl; 
	*/

	double initialYieldCheck;
	int &flag = (element.IntegerData())[ip]; 
	/* if elastic, return elastic stress */
	initialYieldCheck = YieldCondition(fStress, kappa_n, backStress_n);
	if ( initialYieldCheck < yieldFnTol)
	{
		flag = kIsElastic; 
		return fStress; 
	}
	/* ELSE plastic loading */	
	
	flag = kIsPlastic; 
  
	dSymMatrixT eqTrialStress(kNSD);
	eqTrialStress.DiffOf(fStress, backStress_n);

	/*spectral decomposition of equivalent stress*/
	//SpectralDecompT spectre(kNSD);
	spectre.SpectralDecomp_Jacobi(eqTrialStress, true);

	//ArrayT<dSymMatrixT> m(3);
	int i, j;
		
	for (i = 0; i < 3; i++)
	{
		//m[i].Dimension(kNSD);
		m[i].Outer(spectre.Eigenvectors() [i]);
	}

	/* initialize iteration variables */
	int newtonCounter = 0, maxIter = 20;
	dArrayT iterationVars(7), iterationVarsIncr(7), residual(7), residual0(7);
	dSymMatrixT workingStress(kNSD), workingBackStress(kNSD);
	double I1 = 0.0, J2 = 0.0, J3 = 1.0;
	double workingKappa;
	//dArrayT principalEqStress;
	LAdMatrixT dRdX (7); //LAdMatrixT ??????

	iterationVars = 0.0;
	residual = 0.0;
	residual [6] = initialYieldCheck;// 4.720200e+02;//
	residual0 = residual;

	/*
	if (ip == 0 && fFossumDebug)
	{
	       cout << "\n\n\n Initial residual = \n";
		cout << residual << endl;
	}
	*/

	workingStress = fStress;
	workingBackStress = fBackStress;
	workingKappa = fInternal[kkappa];

	principalEqStress = spectre.Eigenvalues();
	//cout << " principalEqStress = \n" << principalEqStress << endl;
	//cout << " m0 = \n" << m[0] << endl;
	//cout << " m1 = \n" << m[1] << endl;
	//cout << " m2 = \n" << m[2] << endl;

	for (i = 0; i < kNSD; i++) I1 += principalEqStress[i];

	for (i = 0; i < kNSD; i++)
	{
		J2 += 0.5*(principalEqStress[i] - I1/3.0) * (principalEqStress[i] - I1/3.0);
		J3 *= (principalEqStress[i] - I1/3.0);
	}

	/*local Newton iteration on variables*/
	while(!ResidualIsConverged(residual, residual0))
	  {
	    /* check to see if not converging */
	    if (newtonCounter++ > maxIter)
	      {
		cout << "FossumSSIsoT::s_ij, Newton Iteration failed to converge\n" << flush;
		throw ExceptionT::kGeneralFail;
	      }

	    /* form dR/dx */
	    dRdX = 0.0; 
	    dRdX = FormdRdX(I1, J2, J3, principalEqStress, workingKappa, workingStress, workingBackStress, iterationVars [6], m);
      
       
	        if (ip == 0 && fFossumDebug)
	      {
		//cout << "residual =\n" << residual << endl;
		//cout << "iterationVars = " << iterationVars << endl;
		//cout << "dRdX = \n" << dRdX << endl << flush;
	      }
	    

	    
		/* break down for static condensation */
		/*	
		dArrayT U(6), V(6), Rtilde(6);
		LAdMatrixT A(6);
		
		for (i = 0; i < 6; i++)
		{
			for (j = 0; j < 6; j++)	A(i,j) = -dRdX(i,j);
			U[i] = -dRdX(i,6);	
			V[i] = dRdX(6,i);
			Rtilde [i] = residual [i];
		}	
		*/
		
		/*solve for dx*/
		iterationVarsIncr.Copy(residual.Pointer());
		iterationVarsIncr*=-1.0;
		//LAdMatrixT dRdXCopy(dRdX);
		//dRdXCopy.Copy(dRdX.Pointer());
		//dRdXCopy.LinearSolve(iterationVarsIncr);

		//if (fFossumDebug)
		//  {
		//cout << "iterationVarsIncr = \n" << iterationVarsIncr << endl;

		iterationVarsIncr = CondenseAndSolve(dRdX, residual);
		  
		/*
		if (fFossumDebug)
		{
			cout << "dRdX = \n" << dRdX << endl << flush;
			cout << "iterationVarsIncr = \n" << iterationVarsIncr << endl;
		}
		*/


		//iterationVarsIncr = residual;

		/*	
		dArrayT AinvU(6), AinvRtilde(6);
		AinvU.Copy(U.Pointer());
		AinvRtilde.Copy(Rtilde.Pointer());
		LAdMatrixT Acopy(6), Acopy2(6);
		Acopy.Copy(A.Pointer());
		Acopy2.Copy(A.Pointer());

		//cout << "U = \n" << AinvU << endl;
		//cout << "Rtilde = \n" << AinvRtilde << endl;

		//cout << "A = \n" << A << endl;
		Acopy.LinearSolve(AinvU);
		//cout << "A = \n" << A << endl;
		A.LinearSolve(AinvRtilde);
		//cout << "A = \n" << A << endl;

		//cout << "AinvU = \n" << AinvU << endl;
		//cout << "AinvRtilde = \n" << AinvRtilde << endl;

		cout << "V = \n" << V << endl;
		iterationVarsIncr[6] = (residual[6] + V.Dot(V, AinvRtilde))/V.Dot(V, AinvU);
		cout << "V = \n" << V << endl;
		Rtilde.AddScaled(-iterationVarsIncr[6], U);
		//Rtilde*=-1;
		Acopy2.LinearSolve(Rtilde); 
		*/

		/* If kappa = kappa0 and kappa wants to increase, hardening is 0 and resolve*/
		/* !!!!!!!!!!!!!!!!!!!!!!!!!!!! undo this to work on capping kappa */
		/*	
		if ( workingKappa - fKappa0 > -1.0e-10*fKappa0 && Rtilde[5] >= 0.0 )
		{ 			
			for (i = 0; i < 7; i++)
			{
				dRdX(i, 5) = 0.0;
				dRdX(5, i) = 0.0;
			}
			dRdX(5,5) = -1.0;
				
			for (i = 0; i < 6; i++)
			{
				for (j = 0; j < 6; j++)	A(i,j) = -dRdX(i,j);
				U[i] = -dRdX(i,6);	
				V[i] = dRdX(6,i);
				Rtilde [i] = residual [i];
			}	
		
			AinvU.Copy(U.Pointer());
			AinvRtilde.Copy(Rtilde.Pointer());
		
			A.LinearSolve(AinvU);
			A.LinearSolve(AinvRtilde);

			iterationVarsIncr[6] = (residual[6] + V.Dot(V, AinvRtilde))/V.Dot(V, AinvU); 
			Rtilde.AddScaled(-iterationVarsIncr[6], U);	
			A.LinearSolve(Rtilde);

		}
		*/
		
        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
		/* and this too */
		/*		
		if (workingKappa + Rtilde[5] > fKappa0)
		{
			double kappaIncr = fKappa0 - workingKappa;
        	
			dArrayT U(5), V(5), Rtilde2(5), AinvU(5), AinvRtilde(5);
			LAdMatrixT A(5);
			for (i = 0; i < 5; i++)
			{
				for (j = 0; j < 5; j++)	A(i,j) = -dRdX(i,j);
				U[i] = -dRdX(i,6);	
				V[i] = dRdX(6,i);
				Rtilde [i] = residual [i] + dRdX(i,5)*kappaIncr;
			}	
			cout << "A=\n" << A << endl;
			
			AinvU.Copy(U.Pointer());
			AinvRtilde.Copy(Rtilde.Pointer());	

			A.LinearSolve(AinvU);
			A.LinearSolve(AinvRtilde);

			iterationVarsIncr[6] = (residual[6] + V.Dot(V, AinvRtilde))/V.Dot(V, AinvU);		
			Rtilde2.AddScaled(-iterationVarsIncr[6], U);	
			A.LinearSolve(Rtilde2);
			
			for (i = 0; i < 5; i++) Rtilde[i] = Rtilde2 [i];
			Rtilde[5] = kappaIncr;			
		}
		*/
        
        
		//for (i = 0; i < 6; i++)
		//iterationVarsIncr[i] = Rtilde[i];
		
		/*
		if (newtonCounter > 2)
			for(i = 0; i < 7; i++)
				for(j = 0; j < 7; j++)
					if ( i != 5 || j != 5)
						dRdX(i,j) = -KroneckerDelta(i,j);
						
		dRdX.LinearSolve(residual);
		iterationVarsIncr = residual;
		*/	
	
		//Do not allow kappa to exceed kappa0
		if (workingKappa + iterationVarsIncr[5] > fKappa0) iterationVarsIncr = CapKappa(residual, dRdX, workingKappa); 
		

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

	    for (i = 0; i < kNSD; i++) I1 += principalEqStress[i];
	    
	    for (i = 0; i < kNSD; i++)
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


	    for (i = 0; i < kNSD; i++) 
	      {
		for (j = 0; j < kNSD; j++)
		  residual [i] += iterationVars[6] *ElasticConstant(i,j) 
		               *dfdSigmaA(I1, J2, J3, principalEqStress [j], workingKappa);
		residual [i] += iterationVars [i];
	      }

	    for (i = kNSD; i < 2*kNSD - 1; i++)
	      {
		residual [i] = iterationVars [6] *fCalpha * Galpha(workingBackStress) * dfdDevStressA (I1, J2, J3, principalEqStress [i - kNSD]);
		residual [i] -= iterationVars [i];
	      } 

	    if ( workingKappa - fKappa0 > -1.0e-12*fabs(fKappa0))
	      residual [5] = 0.0;
	    else   
	      residual [5] = iterationVars [6] * KappaHardening(I1, workingKappa) - iterationVars [5];

	    residual [6] = YieldFn(I1, J2, J3, workingKappa);

	    /*
	    if (ip == 0 && fFossumDebug)
	      {
		//cout << "dRdX = \n" << dRdX << endl << flush;
		//cout << "iterationVars = \n" << iterationVars << endl;
		//cout << "\nresidual =\n" << residual << endl;
		//cout << "workingKappa = " << workingKappa << endl;
		//cout << "fKappa0 = " << fKappa0 << endl;



		
		double r = 0.0;
		for (int i = 0; i < 7; i++)
		  r += residual[i]*residual[i];
		r = sqrt(r);
		cout << "r = " << r << endl;
		
	      }
	    */
	
	} //end while loop


	/*
	cout << "workingStress =\n" << workingStress << endl; 
	cout << "workingBackStress =\n" << workingBackStress << endl;
	cout << "workingKappa = " << workingKappa << endl;
	*/

	/*update variables and exit */
	fStress.AddScaled(iterationVars[0],m[0]);
	fStress.AddScaled(iterationVars[1],m[1]);
	fStress.AddScaled(iterationVars[2],m[2]);

	fDeltaAlpha = 0.0;
	fDeltaAlpha.AddScaled(iterationVars[3],m[0]);
	fDeltaAlpha.AddScaled(iterationVars[4],m[1]);
	fDeltaAlpha.AddScaled(-iterationVars[3]-iterationVars[4],m[2]);

	fDeltaKappa = iterationVars[5];
	
	fInternal[kdgamma] = iterationVars[6];
	
	//fKappa0 = fInternal[kkappa];

	//cout << "fStress = \n" << fStress << endl;
	//cout << "fBackStress = \n" << fBackStress << endl;
	//cout << "fInternal[kkappa] = " << fInternal[kkappa] << endl;
	
	/*
	cout << "fStress [0] = " << fStress [0] << endl;
	cout << "fStress [1] = " << fStress [1] << endl;
	cout << "fStress [2] = " << fStress [2] << endl;
	cout << "fStress [3] = " << fStress [3] << endl;
	cout << "fStress [4] = " << fStress [4] << endl;
	cout << "fStress [5] = " << fStress [5] << endl;
	*/
	
	//fStressInviscid = fStress;

	if (fFossumDebug && ip == 0)
	    cout << " fStress Inviscid = \n" << fStress << endl;

// Rate-Dependence effects. Duvaut-Lions formulation. See Simo and Hughes, p/217
	
	if (fFluidity != 0.0) //fluidity param = 0 => inviscid case
	  {
	    double dt = fSSMatSupport -> TimeStep(); 
            //double time_factor = exp(-1*dt/fFluidity);

	    const dSymMatrixT& e_tot_last = e_last();

	    //cout << "e_tot_last = \n" << e_tot_last << flush;

	    const dSymMatrixT& e_els_last = ElasticStrain(e_tot_last, element, ip);
	    
	    // stress from previous time step
	    dSymMatrixT fStress_last(3);
	    HookeanStress(e_els_last, fStress_last);

	    dSymMatrixT delta_e(3);
	    delta_e.DiffOf(e_tot, e_tot_last);
	    dSymMatrixT elastic_stress_increment(3);
	    HookeanStress(delta_e, elastic_stress_increment);

	    //fStress = time_factor * fStress_last + (1 - time_factor) *fStress
	    //  + (1 - time_factor)/(dt/fFluidity) * elastic_stress_increment;
            fStress *= fTimeFactor;
	    fStress.AddScaled(1- fTimeFactor, fStress_last);
	    fStress.AddScaled((fTimeFactor)*fFluidity/dt, elastic_stress_increment);

	    // update isv's
	    fDeltaAlpha *= fTimeFactor;
	    fDeltaKappa *= fTimeFactor;
	  }

	if (fFossumDebug && ip == 0)
	  {
	    cout << " fStress = \n" << fStress << endl;
	  }


	return fStress;
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
	iterationVarsIncr[5] = fKappa0 - kappa;

	LAdMatrixT modified_dRdX(6); 
	modified_dRdX = 0.0;

	dArrayT modified_residual(6), B(6);
	//dArrayT modified_dX(6);

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

	//cout << "modified_residual = \n" << modified_residual << endl;
	//cout << "B = \n" << B << endl;
	//cout << "modified_dRdX = \n" << modified_dRdX << endl;

	modified_residual.AddScaled(-iterationVarsIncr[5],B);

	//dArrayT modified_residual_copy(6);
	//modified_residual_copy.Copy(modified_residual.Pointer());

	//if (fFossumDebug) modified_residual *= -1;

	//cout << "modified_residual = \n" << modified_residual << endl;
	//cout << "modified_dRdX = \n" << modified_dRdX << endl;

	modified_residual = CondenseAndSolve(modified_dRdX, modified_residual);

	//cout << "modified_residual = \n" << modified_residual << endl;

	//else
	
	//cout << "modified_residual_copy = \n" << modified_residual_copy << endl;
	//cout << "modified_dRdX = \n" << modified_dRdX << endl;

	//modified_residual *= -1;
	//modified_dRdX.LinearSolve(modified_residual);
 
	//cout << "modified_residual_copy = \n" << modified_residual_copy << endl;
  
	for (int i=0; i<5; i++) iterationVarsIncr[i] = modified_residual[i];

	iterationVarsIncr[6] = modified_residual[5];

	//cout << "iterationVarsIncr = \n" << iterationVarsIncr << endl;

	return iterationVarsIncr;
}

dArrayT FossumSSIsoT::CondenseAndSolve(const LAdMatrixT& dRdX, const dArrayT& residual)
{
	int condensedMatrixSize = dRdX.Rows() - 1;

	//cout << "dRdX =\n" << dRdX << endl;
	//cout << "residual =\n" << residual << endl; 

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

	//cout << "A = \n" << A << endl;
	//cout << "U = \n" << U << endl;
	//cout << "V = \n" << V << endl;
	//cout << "Rtilde = \n" << Rtilde << endl;

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
	//Rtilde*=-1;
	A.LinearSolve(Rtilde);

	for (int i = 0; i < condensedMatrixSize; i++) iterationVarsIncr [i] = Rtilde[i];

	//cout << "iterationVarsIncr = \n" << iterationVarsIncr << endl;

	return iterationVarsIncr;
}


double FossumSSIsoT::ElasticConstant(int i, int j)
{
	double lambda = flambda, mu = fmu;
	if (i==j) return flambda + 2*fmu;
	else return flambda;
}

/*
double FossumSSIsoT::Galpha(dSymMatrixT workingStress, double J2)
{
	// if offset fN = 0, already at failure surface, no growth of back stress 
	if (fN == 0.0) return 0.0;
	
	dSymMatrixT devStress;
	double meanStress = workingStress.Trace()/3.0;

	devStress = workingStress;
	devStress.PlusIdentity(-meanStress);

	//cout << " devStress = \n" << devStress << endl;
	//cout << " J2 = \n" << J2 << endl; 

	return 1- (sqrt(.5 * devStress.ScalarProduct()) - sqrt(J2))/fN;
}
*/

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

	//cout << "dPlasticVolStraindX = " << dPlasticVolStraindX(kappa) << endl;
	//cout << "dXdKappa = " << dXdKappa(kappa) << endl;

	return 3 * dfdI1 (I1, kappa) / (dPlasticVolStraindX(kappa) * dXdKappa(kappa));
}

double FossumSSIsoT::dfdDevStressA (double I1, double J2, double J3, double sigmaA)
{
	// cout << dfdJ2(J2, J3) << endl;
	//cout << (sigmaA - I1/3.0) << endl;
	//cout << dfdJ3(J2, J3) << endl;
	//cout <<  ((sigmaA - I1/3.0) * (sigmaA - I1/3.0) - 2.0 * J2 / 3.0) << endl << endl;

	return dfdJ2(J2, J3) * (sigmaA - I1/3.0) + dfdJ3(J2, J3) * ((sigmaA - I1/3.0) * (sigmaA - I1/3.0) - 2.0 * J2 / 3.0);
}

double FossumSSIsoT::dfdSigmaA(double I1, double J2, double J3, double sigmaA, double kappa)
{
	//cout << dfdI1(I1, kappa) << endl;
	//cout << dfdDevStressA(I1, J2, J3, sigmaA) << endl <<  endl;

	return dfdI1(I1, kappa) + dfdDevStressA(I1, J2, J3, sigmaA);
}

double FossumSSIsoT::dfdI1(double I1, double kappa)
{
	//cout << dFfdI1(I1) << endl;
	//cout << dFcdI1(I1, kappa) << endl << endl;

	return -2*YieldFnFfMinusN(I1)*YieldFnFc(I1, kappa)*dFfdI1(I1) - YieldFnFfMinusN(I1) * YieldFnFfMinusN(I1) * dFcdI1(I1, kappa);
}

double FossumSSIsoT::dFfdI1(double I1)
{
	return -fB*fC*exp(fB*I1) - fTheta;
}

double FossumSSIsoT::dFcdI1(double I1, double kappa)
{
	//cout << I1 << endl;
	//cout << Lfn(kappa) << endl;
	//cout << Xfn(kappa) << endl << endl;

	return HeavisideFn(Lfn(kappa) - I1) * -2 * (I1 - Lfn(kappa)) / ((Xfn(kappa) - Lfn(kappa))*(Xfn(kappa) - Lfn(kappa)));
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
	return dLdKappa(kappa) * ( 1 - fR * dFfdI1(Lfn(kappa))); 
}

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

	//cout << " Galpha(workingStress, J2) = " << Galpha(workingStress, J2) << endl;
 
	for (A = 0; A < kNSD - 1; A++)
		for (B = 0; B < kNSD; B++)
		{ 
			//cout << "d2fdDevStressdSigmaB = " << d2fdDevStressdSigmaB (I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B) << endl;
			//cout << "dGalphadSigmaB = " <<  dGalphadSigmaB (workingStress, m[B], principalEqStress[B], I1, J2) << endl;

			dRdX (A + kNSD, B) = Galpha(workingBackStress) * d2fdDevStressdSigmaB (I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B) + dGalphadSigmaB (workingStress, m[B], principalEqStress[B], I1, J2) * dfdDevStressA(I1, J2, J3, principalEqStress[A]);
			dRdX (A + kNSD, B) *= fCalpha * dGamma; 
		}

	//dR(alpha A)/d(dAlpha B)
	for (A = 0; A < kNSD - 1; A++)
		for (B = 0; B < kNSD -1 ; B++)
		{ 
			//cout << "dGalphadAlphaB = " << dGalphadAlphaB (J2, principalEqStress[B], principalEqStress[2]) << endl;

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
		//cout << "Galpha =\n" << Galpha(workingStress, J2) << endl;
		//cout << "dfdDevStressA = \n" << dfdDevStressA(I1, J2, J3, principalEqStress[A]);
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

	//cout << "principalEqStress =\n" << principalEqStress << endl;

	for ( A = 0; A < kNSD; A++)
		dRdX (6, A) = dfdSigmaA(I1, J2, J3, principalEqStress[A], workingKappa);

 
	for ( A = 0; A < kNSD - 1; A++)
		dRdX (6, A + kNSD) = -dfdSigmaA (I1, J2, J3, principalEqStress[A], workingKappa)
				+ dfdSigmaA (I1, J2, J3, principalEqStress[2], workingKappa);
 
	dRdX (6, 5) = dfdKappa(I1, workingKappa);
 
	//for elastic-perfectly plastic only
	//dRdX(3,3) = -1;
	//dRdX(4,4) = -1;
	//dRdX(5,5) = -1;

	return dRdX;
}

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
	double XMinusL = Xfn (kappa) - Lfn (kappa);

	return HeavisideFn( Lfn(kappa) - I1) * -2 / (XMinusL * XMinusL);
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
	/*
	cout << FfMinusN << " " << dFfdI1(I1) << " " << dFcdKappa(I1, kappa) << " "
			<< d2FcdI1dKappa(I1, kappa) << " " <<
			-2 * FfMinusN * dFfdI1(I1) * dFcdKappa(I1, kappa) - FfMinusN * FfMinusN * d2FcdI1dKappa(I1, kappa) << endl;
	*/
	return -2 * FfMinusN * dFfdI1(I1) * dFcdKappa(I1, kappa)
			- FfMinusN * FfMinusN * d2FcdI1dKappa(I1, kappa); //was - - 
}

double FossumSSIsoT::d2FcdI1dKappa(double I1, double kappa)
{
	double L = Lfn(kappa);
	double X = Xfn(kappa);
	/*
	cout << X << " " << L << " " << kappa << endl << X - L << " " << I1 - L << " " << exp(fB * L) << " " <<  2 * HeavisideFn(L - I1) * dLdKappa(kappa) 
		* ((X-L) + 2 * fR * (I1 - L) * (fB * fC * exp (fB * L) + fTheta  ))
		/ ((X-L) * (X-L) * (X-L)) <<  endl;
	*/

	return 2 * HeavisideFn(L - I1) * dLdKappa(kappa) 
		* ((X-L) + 2 * fR * (I1 - L) * (fB * fC * exp (fB * L) + fTheta  ))
		/ ((X-L) * (X-L) * (X-L));
}

double FossumSSIsoT::dGalphadSigmaB (dSymMatrixT workingStress, dSymMatrixT principalDirectionB ,double principalEqStressB, double I1, double J2)
{
	// N = 0 => Galpha is identically 0, does not change with stress or alpha 

	return 0.0;
	/*
	if (fN == 0) return 0;

	dSymMatrixT workingDevStress(3);
	workingDevStress.Deviatoric(workingStress);

	return ( InnerProduct(workingDevStress, principalDirectionB)/sqrt(.5*workingDevStress.ScalarProduct()) - (principalEqStressB - I1/3.0)/sqrt(J2))/(-2*fN);
	*/
}

double FossumSSIsoT::dGalphadAlphaB (dSymMatrixT alpha, dArrayT principalEqStress, int B, ArrayT<dSymMatrixT> m)
{
	/* N = 0 => Galpha is identically 0, does not change with stress or alpha */
	if (fN == 0) return 0;

	dSymMatrixT nB(3);
	nB.DiffOf(m[B], m[2]);

	double J2alpha = .5 * alpha.ScalarProduct();

	/*
	if (fFossumDebug)
	{
		cout << "J2alpha = " << J2alpha << endl;
		cout << "nB = \n" << nB << endl;
		cout << " innerproduct = " << InnerProduct(alpha, nB);
		cout << "alpha = \n" << alpha << endl;
	}
	*/

	if (J2alpha == 0.0) return -1.0/(sqrt(2.0) * fN);

	return -1.0/(2*fN*sqrt(J2alpha))*InnerProduct(alpha, nB);

	//return (principalEqStress3 - principalEqStressB)/(fN*2*sqrt(J2));
}

double FossumSSIsoT::d2fdI1dKappa (double I1, double kappa)
{
	return d2fdSigmaCdKappa(I1, kappa);
}

double FossumSSIsoT::dFcdKappa(double I1, double kappa)
{
	double X = Xfn(kappa), L = Lfn(kappa);

	return 2 * HeavisideFn(L - I1) * dLdKappa(kappa) * ( I1 - L) * ((X - L) + fR * (I1 - L) * ( fTheta + fB * fC * exp (fB * L))) / ((X - L) * (X - L) * (X - L));
}

double FossumSSIsoT::dLdKappa (double kappa)
{
	//return HeavisideFn(fKappa0 - kappa);
	return 1.0;
}

double FossumSSIsoT::d2XdKappadKappa( double kappa)
{
	return fR * fB * fB * fC * exp(fB * Lfn(kappa));
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
 
	/* 
	cout << returnValue << " ";
	returnValue = 0.0;

	for (i = 0; i < kNSD; i++)
		for (j = 0; j < kNSD; j++)
			returnValue += A(i,j)*B(i,j);

	cout << returnValue << endl;
	*/

	return returnValue;
}

/*---------------------------------------------------------*/

/* modulus */
const dMatrixT& FossumSSIsoT::c_ijkl(void)
{
	//fModulus.SumOf(HookeanMatT::Modulus(),
	// ModuliCorrection(CurrentElement(), CurrIP()));
	//int i,j;
	dMatrixT elasticModulus(6), elasticCompliance(6);
	dMatrixT generalizedModulus(14), generalizedCompliance(14);
	dMatrixT d2fdSigmadSigma(6), d2fdqdq(7), dhdq(7);
	dMatrixT d2fdSigmadq(6,7), dhdSigma(7,6);
	dMatrixT Ce = HookeanMatT::Modulus();

	ElementCardT& element = CurrentElement();
	int ip = CurrIP();

	/*
	int flag;
	if (!element.IsAllocated()) flag = -1;
	else flag = (element.IntegerData())[ip];
	*/

	if (!element.IsAllocated() || (element.IntegerData())[ip] != kIsPlastic)
	{
		fModulus = HookeanMatT::Modulus();
		return fModulus;
	}
	else
	{         
		 /* load internal state variables */
		LoadData(element,ip);
		
		double kappaInviscid = fInternal[kkappa] + fDeltaKappa/fTimeFactor;
		dSymMatrixT alphaInviscid = fBackStress; 
		alphaInviscid.AddScaled(1.0/fTimeFactor, fDeltaAlpha);

		/* load stress state, spectral dirs., invariants */
		//dSymMatrixT eqStress(3);
		//eqStress.DiffOf(fStress, fBackStress);

		/*spectral decomposition of equivalent stress*/
		// change so don't have to recalculate these 
	
		/*
		SpectralDecompT spectre(kNSD);
		spectre.SpectralDecomp_Jacobi(eqStress, true);

		ArrayT<dSymMatrixT> m(3);
	
		for (i = 0; i < 3; i++)
		{
			m[i].Allocate(kNSD);
			m[i].Outer(spectre.Eigenvectors() [i]);
		}
		*/
       
		/*Find Invariants */
		//dArrayT principalEqStress (3);
		double I1 = 0.0, J2 =0.0, J3 = 1.0;

		//principalEqStress = spectre.Eigenvalues();

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
       
		d2fdSigmadSigma = D2fdSigmadSigma(I1, J2, J3, kappaInviscid, principalEqStress, m);
		d2fdSigmadq = D2fdSigmadq(I1, J2, J3, kappaInviscid, principalEqStress, m);
		d2fdqdq = D2fdqdq(I1, J2, J3, kappaInviscid, principalEqStress, m);
		dhdSigma = DhdSigma(I1,J2,J3, kappaInviscid, principalEqStress, m, alphaInviscid);
		dhdq = Dhdq(I1,J2,J3, kappaInviscid, principalEqStress, m, alphaInviscid);
 
		/*
		cout << "elasticCompliance = \n" << elasticCompliance << endl;
		cout << "d2fdSigmadSigma = \n" << d2fdSigmadSigma << endl;
		cout << "d2fdqdq = \n" << d2fdqdq << endl;
		//cout << "fModulus = \n" << fModulus << endl;
		*/
		//d2fdSigmadSigma = 0.0;
		//d2fdSigmadq = 0.0;
		//d2fdqdq = 0.0;

		/*
		if (fFossumDebug)
		{
			dMatrixT dummy(7);
			dummy.SetToScaled(-fInternal[kdgamma], dhdq);
			cout << "-fInternal[kdgamma]*dhdq = \n" << dummy << endl;
		}
		*/

		for (int i=0; i<6; i++)
			for (int j=0; j<6; j++)
				generalizedCompliance(i,j) = elasticCompliance(i,j) + fInternal[kdgamma]*d2fdSigmadSigma(i,j);

		for (int i=0; i<6; i++)
			for (int j=6; j<13; j++)
			{
				generalizedCompliance(i,j) = fInternal[kdgamma]*d2fdSigmadq(i, j-6);
				generalizedCompliance(j,i) = -fInternal[kdgamma]*dhdSigma(j-6,i);
			}

		for (int i=6; i<13; i++)
			for (int j=6; j<13; j++)
				generalizedCompliance(i,j) = KroneckerDelta(i,j) - fInternal[kdgamma]*dhdq(i - 6, j - 6);

		//cout << "Generalized Elastic Modulus =\n" << generalizedModulus << endl;

		dArrayT df(13), dr(13), hardeningFns(7), dfdq(7);
		dSymMatrixT dfdSigma(3), dfdAlpha(3);

		dfdSigma = DfdSigma(I1,J2,J3, kappaInviscid, principalEqStress, m);
		dfdAlpha = DfdAlpha(I1,J2,J3, kappaInviscid, principalEqStress, m);
		hardeningFns = Hardening(I1,J2,J3, kappaInviscid, principalEqStress, m, alphaInviscid);                             
		for (int i = 0; i < 6; i++) dfdq[i] = dfdAlpha[i];

		//if (!fFossumDebug)
		//if ( kappaInviscid - fKappa0 <= -1.0e-12*fabs(fKappa0))
		dfdq[6] = dfdKappa(I1, kappaInviscid);

		//double shear components   
		/*
		for (int i=3; i<6; i++)
		{
			dfdSigma[i] += dfdSigma[i];
			//dfdAlpha[i] += dfdAlpha [i];
			hardeningFns[i] += hardeningFns[i];
		}
		*/

		for (int i=0; i<6; i++)
		{
			generalizedCompliance(13, i) = dfdSigma[i];
			generalizedCompliance(13,i+6) = dfdAlpha [i];
			generalizedCompliance(i, 13) = dfdSigma[i];
		}

		generalizedCompliance(13,12) = 0.0;

		//if (!fFossumDebug)
		//if ( kappaInviscid - fKappa0 <= -1.0e-12*fabs(fKappa0))
		generalizedCompliance(13,12) = dfdKappa(I1, kappaInviscid);

		for (int i=0; i<7; i++)
			generalizedCompliance(i+6, 13) = -1.0 *  hardeningFns [i];

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
				//generalizedCompliance(i, j) *= 2.0; 
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

		/*
		if (fFossumDebug)
		{
			//cout << "elasticCompliance = \n" << elasticCompliance << endl;
			cout << "fInternal[kdgamma] = " << fInternal[kdgamma] << endl;
			//cout << "d2fdSigmadSigma = \n" << d2fdSigmadSigma << endl;

			//cout << "fmu = " << fmu << endl;
			//cout << "flambda = " << flambda << endl;
			cout << "generalizedCompliance = \n" << generalizedCompliance << endl;
		}
		*/

		//generalized elastic consistent tangent
		generalizedModulus = generalizedCompliance.Inverse();

		/* pull out modulus */
		for (int i=0; i<6; i++)
			for (int j=0; j<6; j++)
				fModulus(i,j) = generalizedModulus(i,j);

		/*
		for (int i=3; i<6; i++)
			for (int j=0; j<6; j++)
				fModulus(i,j) *= 0.5;
		*/ 

     
      
		//cout << "fModulus = \n" << fModulus << endl;
  
		/*
		dMatrixT comp (6);  
		comp.Inverse(fModulus);
		cout << "compl = \n" << comp << endl << flush;    
		*/
 
		//cout << "fModulus = \n" << fModulus << endl;
		//continuum modulus
		/*
		dMatrixT fModulus2(6);   
		fModulus2 = Ce;

		double chi = dfdSigma.B_ij_A_ijkl_B_kl(Ce) - dfdq.Dot(dfdq,hardeningFns);
		dSymMatrixT CeTimesdfdSigma(3);

		CeTimesdfdSigma.A_ijkl_B_kl(Ce, dfdSigma);
		dMatrixT corr(6);

		//cout << "dFdSigma = \n" << dfdSigma << endl;
		//cout << "dfdq = \n" << dfdq << endl;
		//cout << "hardeningFns = \n" << hardeningFns << endl;

		corr.Outer(CeTimesdfdSigma, CeTimesdfdSigma);
  
		//cout << "chi = " << chi << endl;
		//cout << "corr = \n" << corr << endl;

		fModulus2.AddScaled(-1.0/chi, corr);
		
		if (fFossumDebug) cout << "fModulus2 = \n" << fModulus2 << endl;
  
		*/

	} //end else

	if (fFossumDebug) cout << "fModulus = \n" << fModulus << endl; 
	
	// Rate-Dependence effects. Duvaut-Lions formulation. See Simo and Hughes, p/217
	//double dt = fSSMatSupport -> TimeStep();
	fModulus *= fTimeFactor;
	fModulus.AddScaled(fTimeFactor*fFluidity/(fSSMatSupport -> TimeStep()), Ce);

	if (fFossumDebug) cout << "fModulus = \n" << fModulus << endl; 

 
	return fModulus;
}


/* perfectly plastic modulus */
const dMatrixT& FossumSSIsoT::c_perfplas_ijkl(void)
{
	//fModulusPerfPlas.SumOf(HookeanMatT::Modulus(),
	//ModCorrPerfPlas(CurrentElement(), CurrIP()));
	//int i,j;
	dMatrixT elasticModulus(6), elasticCompliance(6);
	dMatrixT generalizedModulus(14), generalizedCompliance(14);
	dMatrixT d2fdSigmadSigma(6), d2fdqdq(7), dhdq(7);
	dMatrixT d2fdSigmadq(6,7), dhdSigma(7,6);
	dMatrixT Ce = HookeanMatT::Modulus();

	ElementCardT& element = CurrentElement();
	int ip = CurrIP();

	/*
	int flag;
	if (!element.IsAllocated()) flag = -1;
	else flag = (element.IntegerData())[ip];
	*/

	if (!element.IsAllocated() || (element.IntegerData())[ip] != kIsPlastic)
	{
		fModulusPerfPlas = HookeanMatT::Modulus();
		return fModulusPerfPlas;
	}
	else
	{         
		 /* load internal state variables */
		LoadData(element,ip);
		
		double kappa = fInternal[kkappa] + fDeltaKappa;
		dSymMatrixT alpha(3); 
		alpha.SumOf(fBackStress, fDeltaAlpha);

		/* load stress state, spectral dirs., invariants */
		//dSymMatrixT eqStress(3);
		//eqStress.DiffOf(fStress, fBackStress);

		/*spectral decomposition of equivalent stress*/
		// change so don't have to recalculate these 

		/*
		SpectralDecompT spectre(kNSD);
		spectre.SpectralDecomp_Jacobi(eqStress, true);

		ArrayT<dSymMatrixT> m(3);

		for (i = 0; i < 3; i++)
		{
			m[i].Allocate(kNSD);
			m[i].Outer(spectre.Eigenvectors() [i]);
		}
		*/

		/*Find Invariants */
		//dArrayT principalEqStress (3);
		double I1 = 0.0, J2 =0.0, J3 = 1.0;

		//principalEqStress = spectre.Eigenvalues();

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
		/*
		d2fdSigmadq = D2fdSigmadq(I1, J2, J3, kappa, principalEqStress, m);
		d2fdqdq = D2fdqdq(I1, J2, J3, kappa, principalEqStress, m);
		dhdSigma = DhdSigma(I1,J2,J3, kappa, principalEqStress, m, alpha);
		dhdq = Dhdq(I1,J2,J3, kappa, principalEqStress, m, alpha);
		*/
 
		/*
		cout << "elasticCompliance = \n" << elasticCompliance << endl;
		cout << "d2fdSigmadSigma = \n" << d2fdSigmadSigma << endl;
		cout << "d2fdqdq = \n" << d2fdqdq << endl;
		//cout << "fModulusPerfPlas = \n" << fModulusPerfPlas << endl;
		*/
		//d2fdSigmadSigma = 0.0;
		//d2fdSigmadq = 0.0;
		//d2fdqdq = 0.0;

		/*
		if (fFossumDebug)
		{
			dMatrixT dummy(7);
			dummy.SetToScaled(-fInternal[kdgamma], dhdq);
			cout << "-fInternal[kdgamma]*dhdq = \n" << dummy << endl;
		}
		*/

		for (int i=0; i<6; i++)
			for (int j=0; j<6; j++)
				generalizedCompliance(i,j) = elasticCompliance(i,j) + fInternal[kdgamma]*d2fdSigmadSigma(i,j);

		/*
		for (int i=0; i<6; i++)
			for (int j=6; j<13; j++)
			{
				generalizedCompliance(i,j) = fInternal[kdgamma]*d2fdSigmadq(i, j-6);
				generalizedCompliance(j,i) = -fInternal[kdgamma]*dhdSigma(j-6,i);
			}
		*/

		for (int i=6; i<13; i++)
			for (int j=6; j<13; j++)
				//generalizedCompliance(i,j) = KroneckerDelta(i,j) - fInternal[kdgamma]*dhdq(i - 6, j - 6);
				generalizedCompliance(i,j) = KroneckerDelta(i,j);

		//cout << "Generalized Elastic Modulus =\n" << generalizedModulus << endl;

		dArrayT df(13), dr(13), hardeningFns(7), dfdq(7);
		dSymMatrixT dfdSigma(3), dfdAlpha(3);

		dfdSigma = DfdSigma(I1,J2,J3, kappa, principalEqStress, m);
		/*
		dfdAlpha = DfdAlpha(I1,J2,J3, kappa, principalEqStress, m);
		hardeningFns = Hardening(I1,J2,J3, kappa, principalEqStress, m, alpha);                             
		for (int i = 0; i < 6; i++) dfdq[i] = dfdAlpha[i];
		*/

		//if (!fFossumDebug)
		//if ( kappa - fKappa0 <= -1.0e-12*fabs(fKappa0))
		//dfdq[6] = dfdKappa(I1, kappa);

		//double shear components   
		/*
		for (int i=3; i<6; i++)
		{
			dfdSigma[i] += dfdSigma[i];
			//dfdAlpha[i] += dfdAlpha [i];
			hardeningFns[i] += hardeningFns[i];
		}
		*/

		for (int i=0; i<6; i++)
		{
			generalizedCompliance(13, i) = dfdSigma[i];
			//generalizedCompliance(13,i+6) = dfdAlpha [i];
			generalizedCompliance(i, 13) = dfdSigma[i];
		}

		generalizedCompliance(13,12) = 0.0;

		//if (!fFossumDebug)
		//if ( kappa - fKappa0 <= -1.0e-12*fabs(fKappa0))
		//generalizedCompliance(13,12) = dfdKappa(I1, kappa);

		//for (int i=0; i<7; i++)
		//	generalizedCompliance(i+6, 13) = -1.0 *  hardeningFns [i];

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
				//generalizedCompliance(i, j) *= 2.0; 
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

		if (fFossumDebug)
		{
			//cout << "elasticCompliance = \n" << elasticCompliance << endl;
			cout << "fInternal[kdgamma] = " << fInternal[kdgamma] << endl;
			//cout << "d2fdSigmadSigma = \n" << d2fdSigmadSigma << endl;

			//cout << "fmu = " << fmu << endl;
			//cout << "flambda = " << flambda << endl;
			cout << "generalizedCompliance = \n" << generalizedCompliance << endl;
		}

		//generalized elastic consistent tangent
		generalizedModulus = generalizedCompliance.Inverse();

		/* pull out modulus */
		for (int i=0; i<6; i++)
			for (int j=0; j<6; j++)
				fModulusPerfPlas(i,j) = generalizedModulus(i,j);

		/*
		for (int i=3; i<6; i++)
			for (int j=0; j<6; j++)
				fModulusPerfPlas(i,j) *= 0.5;
		*/ 

		if (fFossumDebug) cout << "fModulusPerfPlas = \n" << fModulusPerfPlas << endl;      
      
		//cout << "fModulusPerfPlas = \n" << fModulusPerfPlas << endl;
  
		/*
		dMatrixT comp (6);  
		comp.Inverse(fModulusPerfPlas);
		cout << "compl = \n" << comp << endl << flush;    
		*/
 
		//cout << "fModulusPerfPlas = \n" << fModulusPerfPlas << endl;
		//continuum modulus
     
		dMatrixT fModulus2PerfPlas(6);   
		fModulus2PerfPlas = Ce;

		//double chi = dfdSigma.B_ij_A_ijkl_B_kl(Ce) - dfdq.Dot(dfdq,hardeningFns);
		double chi = dfdSigma.B_ij_A_ijkl_B_kl(Ce);
		dSymMatrixT CeTimesdfdSigma(3);

		CeTimesdfdSigma.A_ijkl_B_kl(Ce, dfdSigma);
		dMatrixT corr(6);

		//cout << "dFdSigma = \n" << dfdSigma << endl;
		//cout << "dfdq = \n" << dfdq << endl;
		//cout << "hardeningFns = \n" << hardeningFns << endl;

		corr.Outer(CeTimesdfdSigma, CeTimesdfdSigma);
  
		//cout << "chi = " << chi << endl;
		//cout << "corr = \n" << corr << endl;

		fModulus2PerfPlas.AddScaled(-1.0/chi, corr);
	 
		if (fFossumDebug) cout << "fModulus2PerfPlas = \n" << fModulus2PerfPlas << endl;
  
	} //end else

	//return Ce;   
	return fModulusPerfPlas;
}



dMatrixT FossumSSIsoT::D2fdSigmadSigma(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m)
{
	int A,B,i,j;
	double d2;
	dMatrixT d2fdSigmadSigma(6);

	/*initialize */
	d2fdSigmadSigma = 0.0;

	//cout << "principalEqStress = \n" << principalEqStress << endl;
	//cout << "I1 " << I1 << endl;
	//cout << "J2 " << J2 << endl;
	//cout << "J3 " << J3 << endl;
	//cout << "kappa " << kappa << endl;

	for (A=0; A<3; A++)
		for (B=0; B<3; B++)//for (B=A; B<3; B++)
		{
			d2 = d2fdSigmaBdSigmaC(I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B, kappa);
			/*
			for (i=0; i<6; i++)
				for (j=i; j<6; j++)
				{
					d2fdSigmadSigma(i,j) += d2 * (m[A]) [i] * (m[B]) [j];
					if ( i != j)
					{
						d2fdSigmadSigma(i,j) += d2 * (m[B]) [i] * (m[A]) [j];
						d2fdSigmadSigma(j,i) = d2fdSigmadSigma(i,j);
					}
				}			
			*/
			for (i=0; i<6; i++)
				for (j=0; j<6; j++)
				{
					d2fdSigmadSigma(i,j) += d2 * (m[A]) [i] * (m[B]) [j]; 
				}
		}
	return d2fdSigmadSigma;
}

dMatrixT FossumSSIsoT::D2fdSigmadq(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m)
{
	int A,B, i,j;
	double d2;
	dMatrixT d2fdSigmadq(6,7);
	ArrayT<dSymMatrixT> n(2);

	d2fdSigmadq = 0.0;

	/*
	for (i = 0; i < 2; i++)
	{
		n[i].Allocate(kNSD);
		n[i].DiffOf( m[i], m[2]);
	}

	for (A=0; A<3; A++)
		for (B=0; B<2; B++)
		{
			d2 = -d2fdDevStressdSigmaB(I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B);
			d2 += d2fdDevStressdSigmaB(I1, J2, J3, principalEqStress[A], principalEqStress[2], A, 2);
			for (i=0; i<6; i++)
				for (j=0; j<6; j++)
					d2fdSigmadq(i,j) += d2 * (m[A]) [i] * (n[B]) [j];	     
		}
	*/

	dMatrixT d2fdSigmadSigma(6);

	d2fdSigmadSigma =  D2fdSigmadSigma(I1, J2, J3, kappa, principalEqStress, m);

	for (i=0; i<6; i++)
		for (j=0; j<6; j++)
			d2fdSigmadq(i,j) = - d2fdSigmadSigma(i,j);	  

	//if (!fFossumDebug)
	//if ( kappa - fKappa0 <= -1.0e-12*fabs(fKappa0))
	//{
	for (A=0; A<3; A++)
		for (i=0; i<6; i++)
			d2fdSigmadq(i,6) += d2fdSigmaCdKappa(I1, kappa) * (m[A]) [i];
	//}

	return d2fdSigmadq;
}

dMatrixT FossumSSIsoT::D2fdqdq(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m)
{
	int A,B, i,j;
	double d2;
	dMatrixT d2fdqdq (7);
	ArrayT<dSymMatrixT> n(2);  

	d2fdqdq = 0.0;
  
	/*
	for (i = 0; i < 2; i++)
	{
		n[i].Allocate(kNSD);
		n[i].DiffOf( m[i], m[2]);
	}

	// d2fdAlphadAlpha 
	for (A=0; A<2; A++)
		for (B=A; B<2; B++)
		{
			d2 = d2fdAlphaBdAlphaC(I1, J2, J3, principalEqStress, A, B, kappa);
			for (i=0; i<6; i++)
				for (j=i; j<6; j++)
				{
					d2fdqdq(i,j) += d2 * (n[A]) [i] * (n[B]) [j];
					if ( i != j)
					{
						d2fdqdq(i,j) += d2 * (n[B]) [i] * (n[A]) [j];
						d2fdqdq(j,i) = d2fdqdq(i,j);
					}
				}
		}
	*/

	dMatrixT d2fdSigmadSigma(6);

	d2fdSigmadSigma =  D2fdSigmadSigma(I1, J2, J3, kappa, principalEqStress, m);

	for (i=0; i<6; i++)
		for (j=0; j<6; j++)
			d2fdqdq(i,j) =  d2fdSigmadSigma(i,j);

	/* d2fdAlphadKappa = 0 */
	/* d2fdKappadKappa */
	//comment out for perfectplasticity
	//if (!fFossumDebug)
	//if ( kappa - fKappa0 <= -1.0e-12*fabs(fKappa0))
	d2fdqdq(6,6) = D2fdKappadKappa(I1, kappa);

	return d2fdqdq;
}

/*
double FossumSSIsoT::d2fdAlphaBdAlphaC(double I1, double J2, double J3, dArrayT  principalEqStress, int B, int C, double kappa)
{
	return d2fdSigmaBdSigmaC(I1, J2, J3, principalEqStress[B], principalEqStress[C], B, C, kappa)
		- d2fdSigmaBdSigmaC(I1, J2, J3, principalEqStress[B], principalEqStress[2], B, 2, kappa)
		- d2fdSigmaBdSigmaC(I1, J2, J3, principalEqStress[2], principalEqStress[C], 2, C, kappa)
		+ d2fdSigmaBdSigmaC(I1, J2, J3, principalEqStress[2], principalEqStress[2], 2, 2, kappa);
}
*/

double FossumSSIsoT::D2fdKappadKappa(double I1, double kappa)
{
	double Ff =YieldFnFf(I1);
	
	//if (kappa - fKappa0 > -10e-12 * fKappa0)
	//  return 0.0;
	//else
	return -Ff*Ff*D2FcdKappadKappa(I1, kappa);
}

double FossumSSIsoT::D2FcdKappadKappa(double I1, double kappa)
{
	double returnValue = 0.0;
	double L = Lfn(kappa);
	double XminusL = Xfn(kappa) - L;
	double I1minusL = I1 - L;
	double dFfdL = dFfdI1(L);
	double d2FfdL2 = d2FfdI1dI1(L);

	returnValue = fR * dFfdL*(-I1minusL*XminusL + 3*I1minusL*I1minusL + XminusL); 
	returnValue += XminusL * (- XminusL -2*I1minusL - fR * d2FfdL2);  
	returnValue *= -2*HeavisideFn(-I1minusL)*dLdKappa(kappa)/(XminusL*XminusL*XminusL*XminusL);

	return returnValue;
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
	
	/*  
	dfdAlpha = 0.0;

	for (i = 0; i < 2; i++)
	{
		n[i].Allocate(kNSD);
		n[i].DiffOf( m[i], m[2]);
	}

	for (A=0; A<2; A++)
	{
		dfdA = - dfdSigmaA(I1, J2, J3, principalEqStress[A], kappa)
				+ dfdSigmaA(I1, J2, J3, principalEqStress[2], kappa);
		for (i=0; i<6; i++)
			dfdAlpha [i] += dfdA *(n[A])[i]; 
	}
	*/

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
   
	if ( kappa - fKappa0 <= -1.0e-12*fabs(fKappa0))
		hardening [6] = KappaHardening(I1, kappa);

	return hardening;
}

dMatrixT FossumSSIsoT::DhdSigma(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m, dSymMatrixT alpha)
{
	int A, B, i, j;
	dMatrixT dhdSigma(7, 6);
	ArrayT<dSymMatrixT> n(2);
	double d2;	
	
	dhdSigma = 0.0;	

	/*	
	for (i = 0; i < 2; i++)
	{
		n[i].Allocate(kNSD);
		n[i].DiffOf( m[i], m[2]);
	}
	*/  

	for (A=0; A<3; A++)
		for (B=0; B<3; B++)
		{
			d2 = fCalpha * Galpha(alpha) * d2fdDevStressdSigmaB(I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B);
			//d2 += fCalpha * dfdDevStressA(I1, J2, J3, principalEqStress[A])
			//  * dGalphadSigmaB (stress, m[B], principalEqStress[B], I1, J2);
			for (i=0; i<6; i++)
				for (j=0; j<6; j++)
					dhdSigma(i,j) += d2 * (m[A]) [i] * (m[B]) [j];	     
		}

	if ( kappa - fKappa0 <= -1.0e-12*fabs(fKappa0))
		for (i=0; i<3; i++)  //only volumetric stresses affect
			dhdSigma(6,i) =  3.0 * d2fdI1dI1(I1, kappa) / (dPlasticVolStraindX( kappa) * dXdKappa ( kappa));

	return dhdSigma;
}

dMatrixT FossumSSIsoT::Dhdq(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m, dSymMatrixT alpha)
{
	int i,j;
	dMatrixT dhdq (7);
	dMatrixT d2fdSdAlpha = D2fdDevStressdAlpha(I1, J2, J3, principalEqStress);
	ArrayT<dSymMatrixT> n(2); 

	dhdq = 0.0;

	/*
	for (i = 0; i < 2; i++)
	{
		n[i].Allocate(kNSD);
		n[i].DiffOf( m[i], m[2]);
	}
	*/ 

	dSymMatrixT dfdDevStress(3);
	dfdDevStress = 0.0;
  
	for (int A = 0; A < kNSD; A++)
		dfdDevStress.AddScaled(dfdDevStressA(I1, J2, J3, principalEqStress[A]), m[A]);

	dSymMatrixT dGalphadAlpha(3);
	dGalphadAlpha = 0.0;

	//for (int B = 0; B < 2; B++)
	//  dGalphadAlpha.AddScaled(dGalphadAlphaB (J2, principalEqStress[B], principalEqStress[2]),n[B]);

	/*
	if (fN != 0)
	{
		for (int A = 0; A < kNSD; A++)
			dGalphadAlpha.AddScaled((principalEqStress[A] - I1/3.0)/(-2*fN*sqrt(J2)),m[A]);
	}
	*/

	if (fN != 0)
	{
		double J2alpha = .5 * alpha.ScalarProduct();

		if (J2alpha == 0.0)
		{
			dSymMatrixT ones(3);
			for (int i = 0; i < 3; i++)
			{
				ones [i] = 1.0;
				ones [i+3] = 2.0;
			}
			dGalphadAlpha.AddScaled(-1.0/(fN*sqrt(2.0)),ones);
		}
		else
			dGalphadAlpha.AddScaled(-1.0/(2.0*fN*sqrt(J2alpha)),alpha);
	}
	/*
	cout << "Galpha = " << Galpha(alpha) << endl;
	cout << "d2fdSdAlpha = \n" << d2fdSdAlpha << endl;
	cout << "dfdDevStress = \n" << dfdDevStress << endl;
	cout << "dGalphadAlpha = \n" << dGalphadAlpha << endl;
	*/

	for (i=0; i<6; i++)
		for (j=0; j<6; j++) 
			dhdq(i,j)  = fCalpha * (Galpha(alpha) * d2fdSdAlpha(i,j)
					+ dfdDevStress[i] * dGalphadAlpha [j]);
 
	if ( kappa - fKappa0 <= -1.0e-12*fabs(fKappa0))
	{
		dhdq(6,6) = dPlasticVolStraindX(kappa) * dXdKappa(kappa) * d2fdI1dKappa (I1, kappa);
		dhdq(6,6) -= dfdI1(I1, kappa) * dPlasticVolStraindX(kappa) * d2XdKappadKappa(kappa);
		dhdq(6,6) -= dfdI1(I1, kappa) * dXdKappa(kappa) * dXdKappa(kappa) * d2PlasticVolStraindXdX(kappa);
		dhdq(6,6) *= 3 / (dPlasticVolStraindX(kappa) * dPlasticVolStraindX(kappa) * dXdKappa (kappa) * dXdKappa (kappa));
	}
  		
	return dhdq;
}

dSymMatrixT FossumSSIsoT::DfdDevStress(double I1, double J2, double J3, dArrayT principalEqStress, ArrayT<dSymMatrixT> m)
{
	dSymMatrixT dfdDevStress(3);
	dfdDevStress = 0.0;
  
	for (int A = 0; A < kNSD; A++)
		dfdDevStress.AddScaled(dfdDevStressA(I1, J2, J3, principalEqStress[A]), m[A]) ;

	return dfdDevStress;
}

dMatrixT FossumSSIsoT::D2fdDevStressdAlpha(double I1, double J2, double J3, dArrayT principalEqStress)
{
	double d2;
	dMatrixT d2fdDevStressdAlpha(6);
	d2fdDevStressdAlpha = 0.0;

	for (int A=0; A<3; A++)
		for (int B=0; B<3; B++)
		{
			d2 =  d2fdDevStressdSigmaB (I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B); 
			for (int i=0; i<6; i++)
				for (int j=0; j<6; j++)
					d2fdDevStressdAlpha(i,j) -= d2 * (m[A]) [i] * (m[B]) [j];
		}

	return d2fdDevStressdAlpha;
}
