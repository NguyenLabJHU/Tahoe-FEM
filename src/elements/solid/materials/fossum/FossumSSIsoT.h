/* 3-invariant, single-surface dilation/compaction plasticity model
 * with isotropic and kinematic hardening
 * Implemented 8/02 Craig Foster
 */

#ifndef _FOSSUM_SS_ISOT_H_
#define _FOSSUM_SS_ISOT_H_

#include "SolidMaterialT.h"
#include "SSIsotropicMatT.h"
#include "HookeanMatT.h"
#include "SSHookeanMatT.h"
#include "SpectralDecompT.h"

/* direct members */
#include "dMatrixT.h"
#include "dArrayT.h"
#include "LAdMatrixT.h"

// Primitive-like fns

namespace Tahoe {

/* forward declarations */
class dSymMatrixT;

/* forward declarations */
class ElementCardT;

 class FossumSSIsoT: public SSIsotropicMatT, public HookeanMatT //, public ParameterInterfaceT
{

public:

	/* constructor */
	// FossumSSIsoT(ifstreamT& in, const SmallStrainT& element, int num_ip, double mu, double lambda);
	FossumSSIsoT(void);

	/* destructor */
	virtual ~FossumSSIsoT(void);

	/* required parameter flags */
	virtual bool Need_Strain_last(void) const {return true;};

	/** material has history variables */
	virtual bool HasHistory(void) const { return true; };

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineParameters(ParameterListT& list) const;

	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/*
	* Returns the value of the yield function given the
	* Cauchy stress vector and state variables, where alpha is
	* the deviatoric stress-like internal state variable
	*/
	double YieldCondition(const dSymMatrixT& stress, const double kappa, dSymMatrixT& alpha);
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
        double fFluidity;   //fluidity parameter, relation time = fFluidity/(2*fmu)
	bool fFossumDebug;

	/* initialization */
	virtual void Initialize(void);

	/* form of tangent matrix (symmetric by default) */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* update internal variables */
	virtual void UpdateHistory(void);

	/* reset internal variables to last converged solution */
	virtual void ResetHistory(void);

	/* moduli */
	virtual const dMatrixT& c_ijkl(void);
	virtual const dMatrixT& ce_ijkl(void);
	virtual const dMatrixT& con_ijkl(void);
	virtual const dMatrixT& c_perfplas_ijkl(void);
	virtual const dMatrixT& con_perfplas_ijkl(void);

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
	* stress and the spatial tangent moduli. Returns true if the
	* determinant of the acoustic tensor is negative and returns
	* the normals and slipdirs. Returns false if the determinant is positive.
	*/
	bool IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs);

protected:

	/* set modulus */

	virtual void SetModulus(dMatrixT& modulus); 
	int loccheck;
 
private:
  
	/* return values */
	dSymMatrixT fStress;
	dSymMatrixT fStrain;
	dMatrixT fModulus;
	dMatrixT fModulusCe;
	dMatrixT fModulusPerfPlas;
	dMatrixT fModulusContinuum;
	dMatrixT fModulusContinuumPerfPlas;


/*--------------------------------------------------------------------*/
// hardening stuff


protected:

	/* status flags */
	enum LoadingStatusT {kIsPlastic = 0,
						kIsElastic = 1,
						kReset = 3}; // indicate not to repeat update

	/* returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain, const ElementCardT& element, int ip);
                        
	/* return correction to stress vector computed by mapping the
	* stress back to the yield surface, if needed */
	const dSymMatrixT& StressCorrection(const dSymMatrixT& trialstrain, ElementCardT& element, int ip); 

	/* return the correction to moduli due to plasticity (if any)
	*
	* Note: Return mapping occurs during the call to StressCorrection.
	*       The element passed in is already assumed to carry current
	*       internal variable values */
	const dMatrixT& ModuliCorrection(const ElementCardT& element, int ip); 

	/* Modulus for checking discontinuous bifurcation */
	const dMatrixT& ModuliCorrPerfPlas(const ElementCardT& element, int ip);

	/* return a pointer to a new plastic element object constructed with
	* the data from element */
	void AllocateElement(ElementCardT& element);

	enum InternalVariablesT {kkappa = 0,  // stress-like internal state variable
				 kdeltakappa = 1,  //increment of kappa 
				 kdgamma = 2}; // consistency parameter
							
	/* element level data */
	void Update(ElementCardT& element);
	void Reset(ElementCardT& element);

	/* returns 1 if the trial elastic strain state lies outside of the 
	* yield surface */

	//int PlasticLoading(const dSymMatrixT& trialstrain, ElementCardT& element, int ip);

	/* computes the deviatoric stress corresponding to the given element
	* and elastic strain.  The function returns a reference to the
	* stress in fDevStress */
	dSymMatrixT& DeviatoricStress(const dSymMatrixT& trialstrain, const ElementCardT& element);

	/* computes the hydrostatic (mean) stress. */
	double MeanStress(const dSymMatrixT& trialstrain, const ElementCardT& element);

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
	//double fDeltaKappa;
	double fX_H;
	double fX;
	double fMeanStress;
	dSymMatrixT fBackStress;
	dSymMatrixT fDeltaAlpha;  

	/* spectral decomp parameters*/
	SpectralDecompT spectre;
	ArrayT<dSymMatrixT> m;
	dArrayT principalEqStress;

	/* return values */
	dSymMatrixT fElasticStrain;
	dSymMatrixT fStressCorr;
	dMatrixT fModuliCorr;
	dMatrixT fModuliCorrPerfPlas;
                
	/* work space */
	dSymMatrixT fDevStress;
	dSymMatrixT fDevStrain; /* deviatoric part of the strain tensor */
	dSymMatrixT IdentityTensor2;  

	dMatrixT fTensorTemp;
	dSymMatrixT One;  
	double fTimeFactor;
	//dSymMatrixT fStressInviscid;

	int &fKappaCapped;
	int fKappaDummy;
        
private:

	double YieldFnGamma(double J2, double J3);
	double YieldFnFfMinusN(double I1);
	double YieldFnFc(double I1, const double kappa);
	int HeavisideFn(double arg);
	double Lfn(const double kappa);
	double Xfn(const double kappa);
	double YieldFnFf(double I1);

	bool StressPointIteration(double initialYieldCheck, dArrayT& iterationVars, dSymMatrixT workingBackStress, double workingKappa);
	bool ResidualIsConverged(dArrayT& residual, dArrayT& residual0);
	dArrayT CapKappa(const dArrayT &residual, const LAdMatrixT &dRdX, const double kappa);
	dArrayT CondenseAndSolve(const LAdMatrixT& dRdX, const dArrayT& residual);

	double ElasticConstant(int i, int j);	
	//double Galpha(dSymMatrixT workingStress, double J2);
	double Galpha(dSymMatrixT alpha);
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

	LAdMatrixT FormdRdX(double I1, double J2, double J3, dArrayT principalEqStress, double workingKappa, dSymMatrixT workingStress, dSymMatrixT workingBackStress, double dGamma, ArrayT<dSymMatrixT> m);

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
	//double dGalphadAlphaB (double J2, double principalEqStressB, double principalEqStress3);
	double dGalphadAlphaB (dSymMatrixT alpha, dArrayT principalEqStress, int B, ArrayT<dSymMatrixT> m);
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
	//double d2fdAlphaBdAlphaC(double I1, double J2, double J3, dArrayT principalEqStress, int B, int C, double kappa);
	double D2fdKappadKappa(double I1, double kappa);
	double D2FcdKappadKappa(double I1, double kappa);
	dSymMatrixT DfdSigma(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m);
	dSymMatrixT DfdAlpha(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m);

	dArrayT Hardening(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m, dSymMatrixT alpha);
	dMatrixT DhdSigma(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m, dSymMatrixT alpha);
	dMatrixT Dhdq(double I1, double J2, double J3, double kappa, dArrayT principalEqStress, ArrayT<dSymMatrixT> m, dSymMatrixT alpha);
	dSymMatrixT DfdDevStress(double I1, double J2, double J3, dArrayT principalEqStress, ArrayT<dSymMatrixT> m);
	dMatrixT D2fdDevStressdAlpha(double I1, double J2, double J3, dArrayT principalEqStress);

};

const double sqrt3__ = sqrt(3.0);

/* Auxiliary Functions for yield function */
inline double FossumSSIsoT::YieldFnGamma(double J2, double J3)
{
	if (J2 <= 0.0)
		return 1.0;   //limit as s_ij -> 0

	double sin3Beta = -3.0*sqrt3__*J3/(2.0*J2*sqrt(J2));

	return .5*(1 + sin3Beta + 1/fPsi*(1.0 - sin3Beta)); 
}

inline double FossumSSIsoT::dGammadJ2 (double J2, double J3)
{
	return 9 * sqrt3__ * J3 * ( 1 - 1/fPsi) / (8 * J2*J2*sqrt(J2));
}

inline double FossumSSIsoT::dfdJ2(double J2, double J3)
{
	return YieldFnGamma (J2, J3) * YieldFnGamma (J2, J3) + 2 * J2 * YieldFnGamma (J2, J3)* dGammadJ2 (J2, J3);
}


} // namespace Tahoe 

#endif/*_FOSSUM_SS_ISOT_H_*/
