/* $Id: FossumSSIsoT.h,v 1.7 2003-02-03 04:40:29 paklein Exp $ */
#ifndef _FOSSUM_SS_ISOT_H_
#define _FOSSUM_SS_ISOT_H_

#include "SolidMaterialT.h"
#include "SSSolidMatT.h"
#include "HookeanMatT.h"
#include "SSHookeanMatT.h"

#include "IsotropicT.h"

/* direct members */
#include "dMatrixT.h"
#include "dArrayT.h"
#include "LAdMatrixT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class dSymMatrixT;

/* forward declarations */
class ElementCardT;

/** 3-invariant, single-surface dilation/compaction plasticity model
 * with isotropic and kinematic hardeneing
 * Implemented 8/02 Craig Foster
 */
class FossumSSIsoT: public SSSolidMatT,
                    public IsotropicT,
                    public HookeanMatT
{
  public:

        /* constructor */
        FossumSSIsoT(ifstreamT& in, const SSMatSupportT& support);
//      FossumSSIsoT(ifstreamT& in, const SSMatSupportT& support, int num_ip, double mu, double lambda);

        /* destructor */
        virtual ~FossumSSIsoT(void);

        /* write parameters to stream */
        virtual void Print(ostream& out) const;
        virtual void PrintName(ostream& out) const;
        
  protected:

        /*
         * Returns the value of the yield function given the
         * Cauchy stress vector and state variables, where alpha is
         * the deviatoric stress-like internal state variable
         */
        double YieldCondition(const dSymMatrixT& stress, 
                              const double kappa, dSymMatrixT& alpha);
	double YieldFn(double I1, double J2, double J3, double kappa);

  protected:
  
  double fA;         // Material Parameter for F_f part of yield fn
  double fB;
  double fC;
  double fTheta;
      
  double fR;         // Ratio of principal radii of ellipse F_c
  double fKappa0;    // Determines starting position of hardening cap

  double fW, fD1, fD2; // Isotropic cap hardening function parameters

  double fCalpha;    // determines rate of growth of back stress tensor

  double fPsi;       // Ratio of failure strength in tension to f. s. in compr.
  double fN;         // offset from initial yield to failure


/*------------------------------------------------------------------*/
//DPSSKStV-like fns

/* initialization */
virtual void Initialize(void);

/* form of tangent matrix (symmetric by default) */
virtual GlobalT::SystemTypeT TangentType(void) const;

/* update internal variables */
virtual void UpdateHistory(void);

/* reset internal variables to last converged solution */
virtual void ResetHistory(void);


        /* modulus */
        virtual const dMatrixT& c_ijkl(void);
        virtual const dMatrixT& cdisc_ijkl(void);
        
        /* stress */
        virtual const dSymMatrixT& s_ij(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. */
	virtual double Pressure(void) const { return fStress.Trace()/3.0; };

        /* returns the strain energy density for the specified strain */
        virtual double StrainEnergyDensity(void);

        /* returns the number of variables computed for nodal extrapolation
         * during for element output, ie. internal variables */
        virtual int  NumOutputVariables(void) const;
        virtual void OutputLabels(ArrayT<StringT>& labels) const;
        virtual void ComputeOutput(dArrayT& output);

        /*
         * Test for localization using "current" values for Cauchy
         * stress and the spatial tangent moduli. Returns 1 if the
         * determinant of the acoustic tensor is negative and returns
         * the normal for which the determinant is minimum. Returns 0
         * of the determinant is positive.
         */
         int IsLocalized(dArrayT& normal);

protected:

        /* set modulus */

        virtual void SetModulus(dMatrixT& modulus); 
         int loccheck;
 
  private:
  
        /* return values */
        dSymMatrixT        fStress;
        dMatrixT        fModulus;
        dMatrixT        fModulusdisc;

/*--------------------------------------------------------------------*/
// hardening stuff


  protected:

        /* status flags */
        enum LoadingStatusT {kIsPlastic = 0,
                             kIsElastic = 1,
                             kReset = 3}; // indicate not to repeat update

        /* returns elastic strain (3D) */
        virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain, 
                const ElementCardT& element, int ip);
                        
        /* return correction to stress vector computed by mapping the
         * stress back to the yield surface, if needed */
        const dSymMatrixT& StressCorrection(const dSymMatrixT& trialstrain, 
                ElementCardT& element, int ip); 

        /* return the correction to moduli due to plasticity (if any)
         *
         * Note: Return mapping occurs during the call to StressCorrection.
         *       The element passed in is already assumed to carry current
         *       internal variable values */
        const dMatrixT& ModuliCorrection(const ElementCardT& element, int ip); 

        /* Modulus for checking discontinuous bifurcation */

        const dMatrixT& ModuliCorrDisc(const ElementCardT& element, int ip);

        /* return a pointer to a new plastic element object constructed with
         * the data from element */
        void AllocateElement(ElementCardT& element);

        enum InternalVariablesT {kkappa = 0,  // stress-like internal state variable
				 kstressnorm = 1,  // norm of stress
				 kdgamma = 2,  // consistency parameter
				 kftrial = 3}; // yield function value
        /* element level data */
        void Update(ElementCardT& element);
        void Reset(ElementCardT& element);

        /* returns 1 if the trial elastic strain state lies outside of the 
         * yield surface */
        int PlasticLoading(const dSymMatrixT& trialstrain, 
                           const ElementCardT& element, 
                           int ip);

        /* computes the deviatoric stress corresponding to the given element
         * and elastic strain.  The function returns a reference to the
         * stress in fDevStress */
        dSymMatrixT& DeviatoricStress(const dSymMatrixT& trialstrain, 
                const ElementCardT& element);

        /* computes the hydrostatic (mean) stress. */
        double MeanStress(const dSymMatrixT& trialstrain,
                const ElementCardT& element);

  private:

        /* load element data for the specified integration point */
        void LoadData(const ElementCardT& element, int ip);

  protected:

        /* element level internal state variables */
        dSymMatrixT fPlasticStrain; //total plastic strain (deviatoric and volumetric)
        dSymMatrixT fUnitNorm;      //unit normal to the yield surface
        dArrayT     fInternal;      //internal variables

  private:

        /* number of integration points */
        int fNumIP;

        /* material parameters **/
        double fmu;
        double flambda;
        double fkappa;
        double fX_H;
        double fX;
        double fMeanStress;
	dSymMatrixT fBackStress;  

        /* return values */
        dSymMatrixT        fElasticStrain;
        dSymMatrixT        fStressCorr;
        dMatrixT        fModuliCorr;
        dMatrixT        fModuliCorrDisc;
                
        /* work space */
        dSymMatrixT fDevStress;
        dSymMatrixT fDevStrain; /* deviatoric part of the strain tensor */
        dSymMatrixT IdentityTensor2;  

        dMatrixT      fTensorTemp;
        dSymMatrixT   One;  
        
 private:

	double YieldFnGamma(double J2, double J3);
	double YieldFnFfMinusN(double I1);
	double YieldFnFc(double I1, const double kappa);
	int HeavisideFn(double arg);
	double Lfn(const double kappa);
	double Xfn(const double kappa);
	double YieldFnFf(double I1);

	bool ResidualIsConverged(dArrayT& residual, dArrayT& residual0);
	double ElasticConstant(int i, int j);	
	double Galpha(dSymMatrixT workingStress, double J2);
	double KappaHardening(double I1, double kappa);
	double dfdDevStressA (double I1, double J2, double J3, double sigmaA);
	double dfdSigmaA(double I1, double J2, double J3, double sigmaA, double kappa);
	double dfdI1(double I1, double kappa);
	double dFfdI1(double I1);
	double dFcdI1(double I1, double kappa);
	double dfdJ2(double J2, double J3);
	double dGammadJ2 (double J2, double J3);
	double dfdJ3(double J2, double J3);
	double dPlasticVolStraindX(double kappa);
	double dXdKappa(double kappa);

LAdMatrixT FormdRdX(double I1, double J2, double J3, dArrayT principalEqStress, double workingKappa, dSymMatrixT workingStress, double dGamma, ArrayT<dSymMatrixT> m);

int KroneckerDelta (int A, int B);
 double d2fdSigmaBdSigmaC (double I1, double J2, double J3, double principalEqStressA, double principalEqStressB, int A, int B, double kappa);
double d2fdDevStressdSigmaB (double I1, double J2, double J3, double principalEqStressA, double principalEqStressB, int A, int B);
double d2fdI1dI1(double I1, double kappa);
double d2FfdI1dI1(double I1);
double d2FcdI1dI1(double I1, double kappa);
double d2fdJ2dJ2 (double J2, double J3);
double d2GammadJ2dJ2(double J2, double J3);
double d2fdJ2dJ3 (double J2, double J3);
double d2GammadJ2dJ3 (double J2);
double dGammadJ3(double J2);
double d2fdJ3dJ3 (double J2, double J3);
double d2fdSigmaCdKappa (double I1, double kappa);
double d2FcdI1dKappa(double I1, double kappa);
double dGalphadSigmaB (dSymMatrixT workingStress, dSymMatrixT principalDirectionB ,double principalEqStressB, double I1, double J2);
double dGalphadAlphaB (double J2, double principalEqStressB, double principalEqStress3);
double d2fdI1dKappa (double I1, double kappa);
double dFcdKappa (double I1, double kappa);
double dLdKappa (double kappa);
double d2XdKappadKappa( double kappa);
double d2PlasticVolStraindXdX(double kappa);
double dfdKappa(double I1, double kappa);
double InnerProduct(dSymMatrixT A, dSymMatrixT B);


dMatrixT D2fdSigmadSigma(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m);
dMatrixT D2fdSigmadq(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m);
dMatrixT D2fdqdq(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m);
double d2fdAlphaBdAlphaC(double I1, double J2, double J3, dArrayT
principalEqStress, int B, int C, double kappa);
double D2fdKappadKappa(double I1, double kappa);
double D2FcdKappadKappa(double I1, double kappa);
dSymMatrixT DfdSigma(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m);
dSymMatrixT DfdAlpha(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m);
 dArrayT Dfdq(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m);

dArrayT Hardening(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m, dSymMatrixT stress);
dMatrixT DhdSigma(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m, dSymMatrixT stress);
dMatrixT Dhdq(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m, dSymMatrixT stress);

};

} // namespace Tahoe 

#endif/*_FOSSUM_SS_ISOT_H_*/
