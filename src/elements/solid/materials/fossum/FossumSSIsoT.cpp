/* 3-invariant, single-surface dilation/compaction plasticity model
 * with isotropic and kinematic hardeneing
 * Implemented 8/02 Craig Foster
 */

//random comment

#include "FossumSSIsoT.h"
#include <iostream.h>
#include "fstreamT.h"
#include "dSymMatrixT.h"

#include "ElementCardT.h"
#include "StringT.h"
#include "DetCheckT.h"
#include <iostream.h>


#include "DPSSLinHardT.h"
#include <iostream.h>
#include <math.h>

#include "iArrayT.h"
#include "ElementCardT.h"
#include "StringT.h"

#include "SpectralDecompT.h"
#include "LAdMatrixT.h"
#include "ArrayT.h"

using namespace Tahoe;

const int    kNumInternal = 4; // number of internal state variables
const double sqrt3       = sqrt(3.0);
//const double sqrt32       = sqrt(3.0/2.0);
const double kYieldTol    = 1.0e-10;
const int    kNSD         = 3;

/*constructor*/

FossumSSIsoT::FossumSSIsoT(ifstreamT& in, const SmallStrainT& element):
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

  SSStructMatT(in, element),
  IsotropicT(in),
  HookeanMatT(3),
  fStress(3),
  fModulus(dSymMatrixT::NumValues(3)),
  fModulusdisc(dSymMatrixT::NumValues(3)),     

  fElasticStrain(kNSD),
  fStressCorr(kNSD),
  fModuliCorr(dSymMatrixT::NumValues(kNSD)),
  fModuliCorrDisc(dSymMatrixT::NumValues(kNSD)), //for disc check
  fDevStress(kNSD),
  fMeanStress(0.0),
  fDevStrain(kNSD), 
  fTensorTemp(dSymMatrixT::NumValues(kNSD)),
  IdentityTensor2(kNSD),
  One(kNSD),
  fBackStress(kNSD)

{
  /* read parameters */
  in >> fA;      if (fA < 0.0 ) {cout << "Bad value for A\n" << flush; 
                                    throw eBadInputValue;}   
  in >> fB;
  in >> fC;
  in >> fTheta; if (fTheta < 0.0 ) {cout << "Bad value for theta\n" << flush; 
                                    throw eBadInputValue;}
  in >> fR;     if (fR < 0.0 ) {cout << "Bad value for R\n" << flush; 
                                    throw eBadInputValue;}  
  in >> fKappa0;
  in >> fW;
  in >> fD1; 
  in >> fD2;
  in >> fCalpha; if (fCalpha < 0.0 ) {cout<< "Bad value for Calpha\n"<< flush; 
                                    throw eBadInputValue;}
  in >> fPsi;     if (fPsi < 0.0 ) {cout << "Bad value for Psi\n" << flush; 
                                    throw eBadInputValue;}
  in >> fN;       if (fN < 0.0) {cout << "Bad value for N\n" << flush;
                                     throw eBadInputValue;}

  /* initialize constant tensor */
  One.Identity();
}


/* destructor */
FossumSSIsoT::~FossumSSIsoT(void) { }

/* write parameters to stream */
void FossumSSIsoT::Print(ostream& out) const
{
  SSStructMatT::Print(out);
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
}

void FossumSSIsoT::PrintName(ostream& out) const
{
  SSStructMatT::PrintName(out);

  out << "Fossum plasticity model: a 3-invariant, single-surface dilation\\\n";
  out << "compaction plasticity model w/ isotropic and kinematic hardening\n"; 

  out << "    Kirchhoff-St.Venant\n";
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
  // not ready yet
  dSymMatrixT equivalentStress; 
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

/* element output data */
const int kNumOutput = 12;
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
            "loccheck",
            "loccheckd", // localization check
};

/* initialization */
void FossumSSIsoT::Initialize(void)
{
        /* inherited */
        HookeanMatT::Initialize();
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



/*discontinuous modulus */
const dMatrixT& FossumSSIsoT::cdisc_ijkl(void)
{
        /* elastoplastic correction */
        fModulusdisc.SumOf(HookeanMatT::Modulus(),
        ModuliCorrDisc(CurrentElement(), CurrIP()));
        return fModulusdisc;
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
        DetCheckT checker(fStress, fModulus);
        //checker.SetElementGroup(ContinuumElement());

        int loccheck= checker.IsLocalized(normal);
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
        for (int i = 0; i < kNumOutput; i++)
                labels[i] = Labels[i];
}

void FossumSSIsoT::ComputeOutput(dArrayT& output)
{

  /* ADD OUTPUT FOR ALPHA, KAPPA */        
output [0] = 0.0;
output [1] = 0.0;
output [2] = 0.0;
output [3] = 0.0;
output [4] = 0.0;
output [5] = 0.0;
output [6] = 0.0;

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

        /* stress-like internal variable alpha */
        const ElementCardT& element = CurrentElement();
        if (element.IsAllocated())
        {
	  //output[0] = fInternal[kalpha];
                const iArrayT& flags = element.IntegerData();
                if (flags[CurrIP()] == kIsPlastic)
                  {
		    //output[0] -= fH_prime*fInternal[kdgamma];
                        
                        // check for localization
                        // compute modulus 
                        const dMatrixT& modulus = c_ijkl();

                        // continuous localization condition checker
                        DetCheckT checker(stress, modulus);
                        dArrayT normal(stress.Rows());
                        output[10] = checker.IsLocalized_SS(normal);
        
                        /* compute discontinuous bifurcation modulus */
                        const dMatrixT& modulusdisc = cdisc_ijkl();

                        /* discontinuous localization condition checker */
                        DetCheckT checkerdisc(stress, modulusdisc);
                        dArrayT normaldisc(stress.Rows());
                        output[11] = checkerdisc.IsLocalized_SS(normaldisc);
                       
		  }
		else
		  {
		    output[10] = 0.0;
		    output[11] = 0.0;
		  }
	}
	else
	  {
	    output[10] = 0.0;
	    output[11] = 0.0;
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
        const dSymMatrixT& trialstrain, 
        ElementCardT& element, int ip)
{
  /*not ready yet*/

  /*initialize*/      
  fStressCorr = 0.0;
  
  return fStressCorr;
}

/* return the correction to moduli due to plasticity (if any)
 *
 * Note: Return mapping occurs during the call to StressCorrection.
 *       The element passed in is already assumed to carry current
 *       internal variable values */
const dMatrixT& FossumSSIsoT::ModuliCorrection(const ElementCardT& element, 
        int ip)
{
        /* initialize */
fModuliCorr = 0.0;

//NOT READY YET

        return fModuliCorr;
}       


/* return the correction to modulus Cep~, checking for discontinuous
 *   bifurcation */


const dMatrixT& FossumSSIsoT::ModuliCorrDisc(const ElementCardT& element, 
        int ip)
{
        /* initialize */

fModuliCorrDisc = 0.0;

//NOT READY YET


        return fModuliCorrDisc;
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
        element.Allocate(i_size, d_size);
        
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
                        //double& dgamma = fInternal[kdgamma];
                        //cout << "kdgamma = " << fInternal[kdgamma] << endl;
                
                        /* internal state variable */
                        //fInternal[kalpha] -= fH_prime*dgamma;

                        /* dev plastic strain increment        */
                        //fPlasticStrain.AddScaled( sqrt32*dgamma, fUnitNorm);
                
                /* vol plastic strain increment        */
                        //fPlasticStrain.AddScaled( fdilation*dgamma/sqrt(3.0), One );
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
        if (!element.IsAllocated()) throw eGeneralFail;

        /* fetch arrays */
        dArrayT& d_array = element.DoubleData();
        
        /* decode */
        int stressdim = dSymMatrixT::NumValues(kNSD);
        int offset    = stressdim*fNumIP;
        int dex       = ip*stressdim;
        
        fPlasticStrain.Set(        kNSD, &d_array[           dex]);
             fUnitNorm.Set(        kNSD, &d_array[  offset + dex]);   
             fBackStress.Set(        kNSD, &d_array[2*offset + dex]);  
             fInternal.Set(kNumInternal, &d_array[3*offset + ip*kNumInternal]);
}

/* returns 1 if the trial elastic strain state lies outside of the 
 * yield surface */
int FossumSSIsoT::PlasticLoading(const dSymMatrixT& trialstrain, 
        const ElementCardT& element, int ip)
{
  /* not yet plastic */
  //if (!element.IsAllocated()) 
  //return( YieldCondition(DeviatoricStress(trialstrain,element),
  //		   MeanStress(trialstrain,element), 0.0) > kYieldTol );
        /* already plastic */
  //      else
  if (element.IsAllocated()) //added 
        {
        /* get flags */
        iArrayT& Flags = element.IntegerData();
                
        /* load internal variables */
        LoadData(element, ip);
                
	//      fInternal[kftrial] = YieldCondition(DeviatoricStress(trialstrain,element), MeanStress(trialstrain,element),fInternal[kalpha]);

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

  fMeanStress = fkappa*trialstrain.Trace();
  return fMeanStress;
}

/*---------------------------------------------------------------*/
/* Auxiliary Functions for yield function */
double FossumSSIsoT::YieldFnGamma(double J2, double J3)
{
  double sin3Beta = -3*sqrt3*J3/(2*J2);

  return .5*(1 + sin3Beta + 1/fPsi*(1 - sin3Beta)); 
}

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
  return HeavisideFn(fKappa0 - kappa) * kappa;
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

	LoadData(element, ip);

	dSymMatrixT backStress_n; //load alpha_n
	double kappa_n; //load kappa  
	backStress_n = fBackStress;
	kappa_n = fInternal[kkappa];

	double initialYieldCheck;     
  

	/* if elastic, return elastic stress */
	initialYieldCheck = YieldCondition(fStress, kappa_n, backStress_n);
	if ( initialYieldCheck < yieldFnTol)
	  return fStress;

	dSymMatrixT eqTrialStress;
	eqTrialStress.DiffOf(fStress, backStress_n);

	/*spectral decomposition of equivalent stress*/
	SpectralDecompT spectre(kNSD);
	spectre.SpectralDecomp_Jacobi(eqTrialStress, true);

	ArrayT<dSymMatrixT> m(3);

	m[0].Outer(spectre.Eigenvectors() [0]);
	m[1].Outer(spectre.Eigenvectors() [1]);
	m[2].Outer(spectre.Eigenvectors() [2]);

	/* initialize iteration variables */
	int i, j, newtonCounter = 0, maxIter = 30;
	dArrayT iterationVars(7), iterationVarsIncr(7), residual(7), residual0(7);
        dSymMatrixT workingStress(kNSD), workingBackStress(kNSD);
	double I1 = 0.0, J2 = 0.0, J3 = 1.0;
	double workingKappa;
	dArrayT principalEqStress;
	LAdMatrixT dRdX (7); //LAdMatrixT ??????

	iterationVars = 0.0;
	residual = 0.0;
	residual [6] = initialYieldCheck;
	residual0 = residual;

	workingStress = fStress;
	workingBackStress = fBackStress;
	workingKappa = fKappa0;

	principalEqStress = spectre.Eigenvalues();

	for (i = 0; i < kNSD; i++)
	  I1 += principalEqStress[i];

	for (i = 0; i < kNSD; i++)
	  {
	    J2 += 0.5*(principalEqStress[i] - I1/3) * (principalEqStress[i] - I1/3);
	    J3 *= (principalEqStress[i] - I1/3);
	  }

	/*local Newton iteration on variables*/
	while(!ResidualIsConverged(residual, residual0))
	  {
	    /* check to see if not converging */
	    if (newtonCounter++ > maxIter)
	      {
		cout << "FossumSSIsoT::s_ij, Newton Iteration failed to converge\n" << flush;
		throw eGeneralFail;
	      }

	    /* form dR/dx */
	    dRdX = 0.0; //implement
	    dRdX = FormdRdX(I1, J2, J3, principalEqStress, workingKappa, workingStress, iterationVars [6], m);

	    /*solve for dx*/
	    dRdX.LinearSolve(residual);

	    /*incr x = x + dx */
	    iterationVars += iterationVarsIncr;

	    /*update working stress and backstress*/
	    principalEqStress [0] += iterationVarsIncr [0] - iterationVarsIncr [3];
	    principalEqStress [1] += iterationVarsIncr [1] - iterationVarsIncr [4];
	    principalEqStress [2] += iterationVarsIncr [2] 
	                          + iterationVarsIncr [3] +iterationVarsIncr [4];

	    /*reset*/
            I1 = 0.0; J2 = 0.0; J3 = 1.0;

	    for (i = 0; i < kNSD; i++)
	      I1 += principalEqStress[i];
	    
	    for (i = 0; i < kNSD; i++)
	      {
		J2 += 0.5 * (principalEqStress[i] - I1/3) * (principalEqStress[i] - I1/3);
		J3 *= (principalEqStress[i] - I1/3);
	      }

	    workingStress.AddScaled(iterationVarsIncr[0],m[0]);
	    workingStress.AddScaled(iterationVarsIncr[1],m[1]);
	    workingStress.AddScaled(iterationVarsIncr[2],m[2]);
	
	    workingBackStress.AddScaled(iterationVarsIncr[3],m[0]);
	    workingBackStress.AddScaled(iterationVarsIncr[4],m[1]);
	    workingBackStress.AddScaled(-iterationVarsIncr[3]-iterationVars[4],m[2]);
	    workingKappa = fKappa0 + iterationVars [5];

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
		residual [i] = iterationVars [6] *fCalpha * Galpha(workingStress, J2) * dfdDevStressA (I1, J2, J3, principalEqStress [i - kNSD]);
		residual [i] -= iterationVars [i];
	      } 
	    residual [5] = iterationVars [6] * KappaHardening(I1, workingKappa) - iterationVars [5];
	    residual [6] = YieldFn(I1, J2, J3, workingKappa);

	  }
	/*update variables and exit */
	fStress.AddScaled(iterationVars[0],m[0]);
	fStress.AddScaled(iterationVars[1],m[1]);
	fStress.AddScaled(iterationVars[2],m[2]);

	fBackStress.AddScaled(iterationVars[3],m[0]);
	fBackStress.AddScaled(iterationVars[4],m[1]);
	fBackStress.AddScaled(-iterationVars[3]-iterationVars[4],m[2]);

	fInternal[kkappa] += iterationVars[5];
	fInternal[kdgamma] = iterationVars[6];

	return fStress;
}

bool FossumSSIsoT::ResidualIsConverged(dArrayT& residual, dArrayT& residual0)
{
	double yieldFnTol = 1.0e-8, stressTol = 1.0e-8;
	int i;

	for (i = 0; i < 6; i++)
	  if (residual0 [i] == 0.0)
	    residual0 [i] = residual [i];

	for (i = 0; i < 6; i++)	
	  if (residual [i] > stressTol && residual [i] > stressTol*residual0 [i])
	    return false;

	if (residual [6] > stressTol && residual [6] > stressTol*residual0 [6])
	  return false;

	return true;
}

double FossumSSIsoT::ElasticConstant(int i, int j)
  {
    if (i==j)
      return flambda + 2*fmu;
    else
      return flambda;
  }


double FossumSSIsoT::Galpha(dSymMatrixT workingStress, double J2)
{
  dSymMatrixT devStress;
  double pressure = (workingStress [1] + workingStress [2] + workingStress [3])/3;

  devStress = workingStress;
  devStress.PlusIdentity(-pressure);

  return 1- sqrt(.5 * devStress.ScalarProduct() - sqrt(J2))/fN;
}

double FossumSSIsoT::KappaHardening(double I1, double kappa)
{
return 3 * dfdI1 (I1, kappa) / (dPlasticVolStraindX(kappa) * dXdKappa(kappa));
}

double FossumSSIsoT::dfdDevStressA (double I1, double J2, double J3, double sigmaA)
{
 return dfdJ2(J2, J3) * (sigmaA - I1/3) + dfdJ3(J2, J3) * ((sigmaA - I1/3) * (sigmaA - I1/3) - 2* J3 / 3);
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
  return HeavisideFn(Lfn(kappa) - I1) * -2 * (I1 - Lfn(kappa)) / ((Xfn(kappa) - Lfn(kappa))*(Xfn(kappa) - Lfn(kappa)));
}

double FossumSSIsoT::dfdJ2(double J2, double J3)
{
  return YieldFnGamma (J2, J3) * YieldFnGamma (J2, J3) + 2 * J2 * YieldFnGamma (J2, J3)* dGammadJ2 (J2, J3);
}

double FossumSSIsoT::dGammadJ2 (double J2, double J3)
{
  return 9 * sqrt3 * J3 * ( 1 - 1/fPsi) / (8 * pow(J2,2.5));
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

LAdMatrixT FossumSSIsoT::FormdRdX(double I1, double J2, double J3, dArrayT principalEqStress, double workingKappa, dSymMatrixT workingStress, double dGamma, ArrayT<dSymMatrixT> m)
{
  int A, B, C;
  LAdMatrixT dRdX (7);

  dRdX = 0.0;

  for (A = 0; A < kNSD; A++)
    for (B = 0; B < kNSD; B++)
      {
      for (C = 0; C < kNSD; C++)
	dRdX (A,B) += dGamma * ElasticConstant(A,C) * d2fdSigmaBdSigmaC (I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B, workingKappa); 
      dRdX += KroneckerDelta (A, B);
      }

  for (A = 0; A < kNSD; A++)
    for (B = 0; B < kNSD - 1; B++)
      for (C = 0; C < kNSD; C++)
	dRdX (A, B + kNSD) += dGamma * ElasticConstant(A,C) * (- d2fdSigmaBdSigmaC (I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B, workingKappa) + d2fdSigmaBdSigmaC (I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B, workingKappa));
  
 for (A = 0; A < kNSD; A++)
      for (C = 0; C < kNSD; C++)
	{
	  dRdX (A,5) += dGamma * ElasticConstant(A,C) * d2fdSigmaCdKappa (I1, workingKappa); 
	  dRdX (A, 6) += ElasticConstant(A,C)*dfdSigmaA(I1, J2, J3, principalEqStress [A], workingKappa);
	}
      
 /* ------ */

 for (A = 0; A < kNSD - 1; A++)
   for (B = 0; B < kNSD; B++)
     { 
       dRdX (A + kNSD, B) = Galpha(workingStress, J2) * d2fdDevStressdSigmaB (I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B) + dGalphadSigmaB (workingStress, m[B], principalEqStress[B], I1, J2) * dfdDevStressA(I1, J2, J3, principalEqStress[B]);
       dRdX (A + kNSD, B) *= fCalpha * dGamma; 
     }

 
 for (A = 0; A < kNSD - 1; A++)
   for (B = 0; B < kNSD -1 ; B++)
     { 
       dRdX (A + kNSD, B + kNSD) = - Galpha(workingStress, J2) * d2fdDevStressdSigmaB (I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B) +  Galpha(workingStress, J2) * d2fdDevStressdSigmaB (I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B) + dGalphadAlphaB (J2, principalEqStress[B], principalEqStress[2]) * dfdDevStressA(I1, J2, J3, principalEqStress[B]);
       dRdX (A + kNSD, B + kNSD) *= fCalpha * dGamma;
       dRdX (A + kNSD, B + kNSD) -= KroneckerDelta (A, B);
     }

 //dR(Alpha A)/d(Kappa) = 0

 for (A = 0; A < kNSD - 1; A++)
   dRdX (A + kNSD, 6) = fCalpha * Galpha(workingStress, J2) *
                        dfdDevStressA(I1, J2, J3, principalEqStress[A]);

 /* ------------ */

 for ( A = 0; A < kNSD; A++)
   dRdX (5, A) = 3 * dGamma * d2fdI1dI1(I1, workingKappa) 
                 / (dPlasticVolStraindX( workingKappa) * dXdKappa (workingKappa));

 dRdX (5, 5) = dPlasticVolStraindX(workingKappa) * dXdKappa(workingKappa) * d2fdI1dKappa (I1, workingKappa);
 dRdX (5, 5) -= dfdI1(I1, workingKappa) * dPlasticVolStraindX(workingKappa) * d2XdKappadKappa(workingKappa);
 dRdX (5, 5) -= dfdI1(I1, workingKappa) * dXdKappa(workingKappa) * d2PlasticVolStraindXdX(workingKappa);
 dRdX (5, 5) *= 3 * dGamma / (dPlasticVolStraindX(workingKappa) * dPlasticVolStraindX(workingKappa) * dXdKappa (workingKappa) * dXdKappa (workingKappa));
 dRdX (5, 5) -= 1;
 
 dRdX (5, 6) = KappaHardening(I1, workingKappa);

  /*------------- */

 for ( A = 0; A < kNSD; A++)
   dRdX (6, A) = dfdSigmaA(I1, J2, J3, principalEqStress[A], workingKappa);

 for ( A = 0; A < kNSD - 1; A++)
   dRdX (6, A + kNSD) = -dfdSigmaA (I1, J2, J3, principalEqStress[A], workingKappa)
                        + dfdSigmaA (I1, J2, J3, principalEqStress[2], workingKappa);

 dRdX (6, 5) = dfdKappa(I1, workingKappa);

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
  double returnValue = 0.0;
  
  returnValue += dfdJ2(J2, J3) * (KroneckerDelta(A,B));
  returnValue += dfdJ3(J2, J3) * (2 * (principalEqStressA - I1/3) * KroneckerDelta (A, B) - 1/3) - 2/3 *(principalEqStressB - I1/3);
  returnValue += d2fdJ2dJ2 (J2, J3) * (principalEqStressA - I1/3) * (principalEqStressB - I1/3);
returnValue += d2fdJ2dJ3 (J2, J3) * ((principalEqStressA - I1/3) * (principalEqStressA + principalEqStressB - 2/3 * I1) - 2/3 * J2); 
returnValue += d2fdJ3dJ3 (J2, J3) * ((principalEqStressA - I1/3) *(principalEqStressA - I1/3)  - 2/3 * J2);

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
  return (1 - 1/fPsi) * -45*sqrt3 * J3/ (16 * pow(J2, 3.5));
}

double FossumSSIsoT::d2fdJ2dJ3 (double J2, double J3)
{
  double gamma = YieldFnGamma(J2, J3);
  double dGdJ3 = dGammadJ3 (J2);

  return 2 * (gamma * dGdJ3 + J2 * dGammadJ2(J2, J3) * dGdJ3 + J2 * gamma * d2GammadJ2dJ3(J2) );
}

double FossumSSIsoT::d2GammadJ2dJ3 (double J2)
{
  return ( 1 - 1/fPsi) * 9 * sqrt3 / (8 * pow(J2, 2.5));
}

double FossumSSIsoT::dGammadJ3(double J2)
{
  return ( 1 - 1/fPsi) * 3 * sqrt3 / (4 * pow(J2, 2.5));
}

double FossumSSIsoT::d2fdJ3dJ3 (double J2, double J3)
{
  double dGdJ3 = dGammadJ3(J2);

  return 2 * J2 * dGdJ3 * dGdJ3; 
}

double FossumSSIsoT::d2fdSigmaCdKappa (double I1, double kappa)
{
  double FfMinusN = YieldFnFfMinusN(I1);

  return -2 * FfMinusN * dFfdI1(I1) * dFcdKappa(I1, kappa)
         - FfMinusN * FfMinusN * d2FcdI1dKappa(I1, kappa);
}

double FossumSSIsoT::d2FcdI1dKappa(double I1, double kappa)
{
  double L = Lfn(kappa);
  double X = Xfn(kappa);

  return 2 * HeavisideFn(L - I1) * dLdKappa(kappa) 
         * ((X-L) + 2 * fR * (I1 - L) * (fB * fC * exp (fB * L) + fTheta  )
	    / ((X-L) * (X-L) * (X-L)));
}

double FossumSSIsoT::dGalphadSigmaB (dSymMatrixT workingStress, dSymMatrixT principalDirectionB ,double principalEqStressB, double I1, double J2)
{
  dSymMatrixT workingDevStress;
  workingDevStress.Deviatoric(workingStress);

  return ( InnerProduct(workingDevStress, principalDirectionB)/sqrt(workingStress.ScalarProduct()) - (principalEqStressB - I1/3)/sqrt(2*J2))/fN;
}

double FossumSSIsoT::dGalphadAlphaB (double J2, double principalEqStressB, double principalEqStress3)
{
  return (principalEqStressB - principalEqStress3)/(fN*sqrt(2*J2));
}

double FossumSSIsoT::d2fdI1dKappa (double I1, double kappa)
{
  return -2 * YieldFnFfMinusN (I1) * dFcdKappa (I1, kappa) * dFfdI1 (I1)
         - pow(YieldFnFfMinusN (I1),2) * d2FcdI1dKappa (I1, kappa);
}

double FossumSSIsoT::dFcdKappa(double I1, double kappa)
{
  double X = Xfn(kappa), L = Lfn(kappa);

  return -2 * HeavisideFn(L - I1) * dLdKappa(kappa) * ( I1 - L) * ((X - L) + fR * (I1 - L) * ( fTheta + fB * fC * exp (fB * L))) / ((X - L) * (X - L) * (X - L));
}

double FossumSSIsoT::dLdKappa (double kappa)
{
  return HeavisideFn(fKappa0 - kappa);
}

double FossumSSIsoT::d2XdKappadKappa( double kappa)
{
  return fR * fB * fB * fC * exp(fB * Lfn(kappa));
}

double FossumSSIsoT::d2PlasticVolStraindXdX(double kappa)
{
double XminusX0 = Xfn (kappa) - Xfn (fKappa0);

 return fW * (fD1 - 2 * fD2 * ( XminusX0 - 1) * exp ( fD1 * XminusX0 - fD2 * XminusX0 * XminusX0));  
}

double FossumSSIsoT::dfdKappa(double I1, double kappa)
{
  double Ff =YieldFnFf(I1);

  return -Ff*Ff*dFcdKappa(I1, kappa);
}

double FossumSSIsoT::InnerProduct(dSymMatrixT A, dSymMatrixT B)
{
 double returnValue = 0.0;
 int i;

 for (i = 0; i < kNSD; i++)
   returnValue += A[i]*B[i];
 
 for (i = kNSD; i < 2 * kNSD; i++)
   returnValue += 2*A[i]*B[i];

 return returnValue;
}

/*---------------------------------------------------------*/


/* modulus */
const dMatrixT& FossumSSIsoT::c_ijkl(void)
{
  //fModulus.SumOf(HookeanMatT::Modulus(),
  // ModuliCorrection(CurrentElement(), CurrIP()));
  int i,j;
  dSymMatrixT elasticModulus, elasticCompliance;
  dSymMatrixT generalizedModulus(12), generalizedCompliance(12);
  dSymMatrixT d2fdSigmadSigma(6), d2fdqdq(7);
  dMatrixT  d2fdSigmadq(6,7);

  ElementCardT& element = CurrentElement();
  int ip = CurrIP();

  if (!element.IsAllocated() || (element.IntegerData())[ip] != kIsPlastic)
    return HookeanMatT::Modulus();
  else

     {         

       /* load internal state variables */
       LoadData(element,ip);

       /* load stress state, spectral dirs., invariants */
	dSymMatrixT eqStress;
	eqStress.DiffOf(fStress, fBackStress);

	/*spectral decomposition of equivalent stress*/
	// change so don't have to recalculate these 
	SpectralDecompT spectre(kNSD);
	spectre.SpectralDecomp_Jacobi(eqStress, true);

	ArrayT<dSymMatrixT> m(3);

	m[0].Outer(spectre.Eigenvectors() [0]);
	m[1].Outer(spectre.Eigenvectors() [1]);
	m[2].Outer(spectre.Eigenvectors() [2]);
       
	/*Find Invariants */
	dArrayT principalEqStress (3);
	double I1 = 0.0, J2 =0.0, J3 = 1.0;

	principalEqStress = spectre.Eigenvalues();

	for (i = 0; i < kNSD; i++)
	  I1 += principalEqStress[i];

	for (i = 0; i < kNSD; i++)
	  {
	    J2 += 0.5*(principalEqStress[i] - I1/3) * (principalEqStress[i] - I1/3);
	    J3 *= (principalEqStress[i] - I1/3);
	  }

	/*form generalized modulus */
	dSymMatrixT Ce = Ce.FromMatrix(HookeanMatT::Modulus());

       elasticCompliance.Inverse(Ce);
       generalizedCompliance = 0.0;
       
       d2fdSigmadSigma = D2fdSigmadSigma(I1, J2, J3, fInternal[kkappa], principalEqStress, m);
       d2fdSigmadq = D2fdSigmadq(I1, J2, J3, fInternal[kkappa], principalEqStress, m);
       d2fdqdq = D2fdqdq(I1, J2, J3, fInternal[kkappa], principalEqStress, m);

       for (i=0; i<6; i++)
	 for (j=i; j<6; j++)
	   generalizedCompliance(i,j) = 
	     elasticCompliance(i,j) + fInternal[kdgamma]*d2fdSigmadSigma(i,j);

       for (i=0; i<6; i++)
	 for (j=6; j<13; j++)
	   generalizedCompliance(i,j) = 
	     fInternal[kdgamma]*d2fdSigmadq(i, j - 6);

       for (i=6; i<13; i++)
	 for (j=i; j<13; j++)
	   generalizedCompliance(i,j) = 
	     -KroneckerDelta(i,j) + fInternal[kdgamma]*d2fdqdq(i,j);

       //generalized elastic consistent tangent
       generalizedModulus = generalizedCompliance.Inverse();

       // form vector df
       dArrayT df(13);
       dSymMatrixT dfdSigma(6), dfdAlpha(6);

       dfdSigma = DfdSigma(I1,J2,J3, fInternal[kkappa], principalEqStress, m);
       dfdAlpha = DfdAlpha(I1,J2,J3, fInternal[kkappa], principalEqStress, m);

       for (i=0; i<6; i++)
	 {
	   df[i] = dfdSigma[i];
	   df[i+6] = dfdAlpha [i];
	 }

       /* double shear components */
       for (i=3; i<6; i++)
	 {
	   df[i] += dfdSigma[i];
	   df[i+6] += dfdAlpha [i];
	 }

       df [12] = dfdKappa(I1, fInternal[kkappa]);

       /* A = A - A:df @ A:df / df:A:df */
       //rename work, 3 places
       dSymMatrixT work;
       dArrayT Atimesdf;

       generalizedModulus.Multx(df, Atimesdf);
       work.Outer(Atimesdf);

       generalizedModulus.AddScaled(-1/(generalizedModulus.MultmBn(df, df)), work);

	 /* pull out modulus */
       for (i=0; i<6; i++)
	 for (j=i; j<6; j++)
	   fModulus(i,j) = generalizedModulus(i,j);

     }
  return fModulus;
}

dSymMatrixT FossumSSIsoT::D2fdSigmadSigma(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m)
{
  int A,B,i,j;
  double d2;
  dSymMatrixT d2fdSigmadSigma(6);

  /*initialize */
  d2fdSigmadSigma = 0.0;

  for (A=0; A<3; A++)
    for (B=A; B<3; B++)
      {
	d2 = d2fdSigmaBdSigmaC(I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B, kappa);
	for (i=0; i<6; i++)
	  for (j=i; j<6; j++)
	    {
	      d2fdSigmadSigma(i,j) += d2 * (m[A]) [i] * (m[B]) [j];
	      if ( i != j)
		d2fdSigmadSigma(i,j) += d2 * (m[B]) [i] * (m[A]) [j];
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

  n[0].DiffOf( m[0], m[2]);
  n[1].DiffOf( m[1], m[2]);

 for (A=0; A<3; A++)
    for (B=0; B<2; B++)
      {
	d2 = d2fdDevStressdSigmaB(I1, J2, J3, principalEqStress[A], principalEqStress[B], A, B);
	for (i=0; i<6; i++)
	  for (j=0; j<6; j++)
	      d2fdSigmadq(i,j) += d2 * (m[A]) [i] * (n[B]) [j];	     
      }
 
 for (A=0; A<3; A++)
   for (i=0; i<6; i++)
     d2fdSigmadq(i,6) += d2fdSigmaCdKappa(I1, kappa) * (m[A]) [i];

  return d2fdSigmadq;
}

dSymMatrixT FossumSSIsoT::D2fdqdq(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m)
{
  int A,B, i,j;
  double d2;
  dSymMatrixT d2fdqdq;
  ArrayT<dSymMatrixT> n(2);  

  d2fdqdq = 0.0;

  n[0].DiffOf( m[0], m[2]);
  n[1].DiffOf( m[1], m[2]);

  /* d2fdAlphadAlpha */
  for (A=0; A<2; A++)
    for (B=A; B<2; B++)
      {
	d2 = d2fdAlphaBdAlphaC(I1, J2, J3, principalEqStress, A, B, kappa);
	for (i=0; i<6; i++)
	  for (j=i; j<6; j++)
	    {
	      d2fdqdq(i,j) += d2 * (n[A]) [i] * (n[B]) [j];
	      if ( i != j)
		d2fdqdq(i,j) += d2 * (n[B]) [i] * (n[A]) [j];
	    }
      }

  /* d2fdAlphadKappa = 0 */

  /* d2fdKappadKappa */
  d2fdqdq(6,6) = D2fdKappadKappa(I1, kappa);

 return d2fdqdq;
}


double FossumSSIsoT::d2fdAlphaBdAlphaC(double I1, double J2, double J3, dArrayT  principalEqStress, int B, int C, double kappa)
{
  return d2fdSigmaBdSigmaC(I1, J2, J3, principalEqStress[B], principalEqStress[C], B, C, kappa)
    - d2fdSigmaBdSigmaC(I1, J2, J3, principalEqStress[B], principalEqStress[2], B, 2, kappa)
    - d2fdSigmaBdSigmaC(I1, J2, J3, principalEqStress[2], principalEqStress[C], 2, C, kappa)
    + d2fdSigmaBdSigmaC(I1, J2, J3, principalEqStress[2], principalEqStress[2], 2, 2, kappa);
}

double FossumSSIsoT::D2fdKappadKappa(double I1, double kappa)
{
  double Ff =YieldFnFf(I1);

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
  int A, i,j;
  double dfdA;
  dSymMatrixT dfdSigma(3);

  dfdSigma = 0.0;

  for (A=0; A<3; A++)
    {
      dfdA = dfdSigmaA(I1, J2, J3, principalEqStress[A], kappa);
	for (i=0; i<3; i++)
	  for (j=i; j<3; j++)
	    dfdSigma(i,j) += dfdA *(m[A])(i,j); 
    }

  return dfdSigma;
}

dSymMatrixT FossumSSIsoT::DfdAlpha(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m)
{
  int A, i,j;
  double dfdA;
  dSymMatrixT dfdAlpha(3);
  ArrayT<dSymMatrixT> n(2);  

  n[0].DiffOf( m[0], m[2]);
  n[1].DiffOf( m[1], m[2]);

  for (A=0; A<2; A++)
    {
      dfdA = - dfdSigmaA(I1, J2, J3, principalEqStress[A], kappa)
	     + dfdSigmaA(I1, J2, J3, principalEqStress[2], kappa);
	for (i=0; i<3; i++)
	  for (j=i; j<3; j++)
	    dfdAlpha(i,j) += dfdA *(n[A])(i,j); 
    }
  return dfdAlpha;
}



