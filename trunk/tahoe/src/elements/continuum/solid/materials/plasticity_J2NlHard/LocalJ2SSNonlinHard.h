/* $Id: LocalJ2SSNonlinHard.h,v 1.8 2004-07-15 08:29:01 paklein Exp $ */
#ifndef _LOCAL_J2_SS_NONLIN_HARD_H_
#define _LOCAL_J2_SS_NONLIN_HARD_H_

/* base classes */
#include "SSIsotropicMatT.h"
#include "HookeanMatT.h"

/* direct members */
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"

#include "GlobalT.h"

namespace Tahoe {

class ElementCardT;
class ifstreamT;

class LocalJ2SSNonlinHard: public SSIsotropicMatT, public HookeanMatT
{
public:

	/* constructor */
	LocalJ2SSNonlinHard(ifstreamT& in, const SSMatSupportT& support);

	/* initialization */
	virtual void Initialize(void);

	/* form of tangent matrix (symmetric by default) */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* update internal variables */
	virtual void UpdateHistory(void);

	/* reset internal variables to last converged solution */
	virtual void ResetHistory(void);
	
	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. */
	virtual double Pressure(void) const { return fStress.Trace()/3.0; };
	/*@}*/

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);

	/* returns the number of variables computed for nodal extrapolation
	 * during for element output, ie. internal variables */
	virtual int NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);

	/* returns elastic strain */
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain,
		const ElementCardT& element, int ip);
			
	

	/* solve for the state at selected ip */
	void SolveState(ElementCardT& element, int ip);

	/* return pointers to new element objects constructed with the data from elements */
	void AllocateAllElements();

	/* indexes to access internal variable (scalars) array */
	enum InternalVariablesT {kIsotHard = 0,  // isotropic hardening
                                 kdelLmbda = 1,  // consistency parameter
			         kYieldCrt = 2}; // yield criteria

	/* returns the value of the yield function given the relative stress and isotropic hardening */
	double YieldCondition(const dSymMatrixT& relstress, double isotropic) const;

private:
	/* status flags */
	enum LoadingStatusT {kIsPlastic = 0,
                             kIsElastic = 1,
                                 kReset = 3}; // indicator not to repeat update

	/* load element data for the specified integration point */
	void LoadData(const ElementCardT& element, int ip);

	/* computes the increment in the plasticity parameter */
	void IncrementPlasticParameter(double& varLambda);

	/* computes the increments in the stress and internal variables */
	void IncrementState(const double& varLambda);

	/* computes the unit normal and the yield condition */
	void UpdateState();

	/* computes the consistent tangent moduli */
	void TangentModuli();

protected:

	/* element level internal variables at current time step*/
	dSymMatrixT fStress;        //stress
	dSymMatrixT fPlstStrn;      //plastic strain
	dSymMatrixT fUnitNorm;      //unit normal to the stress surface
	dSymMatrixT fKineHard;      //stress surface "center", kinematic hardening
	dArrayT     fInternal;      //internal variables

	/* element level internal variables at previous time step*/
	dSymMatrixT fStress_n;      //stress
	dSymMatrixT fPlstStrn_n;    //plastic strain
	dSymMatrixT fUnitNorm_n;    //unit normal to the stress surface
	dSymMatrixT fKineHard_n;    //stress surface "center", kinematic hardening
	dArrayT     fInternal_n;    //internal variables

private:

	/* number of integration points */
	int fNumIP;

	/* material parameters **/
	double fmu;

    	/* material constants */
        double yield;
        double k1, k2, k3, k4;

	/* return values */
	dSymMatrixT	fElasticStrain;
	dMatrixT	fModulus;
	dMatrixT	fModuliCorr;

        /* general workspaces */
	dSymMatrixT fRelStress;
        dSymMatrixT fsymmatx1;
        dMatrixT    fmatx1;
        dMatrixT    fmatx2;
        dMatrixT    fmatx3;
	dMatrixT    ftnsr1;

};

} // namespace Tahoe
#endif /* _LOCAL_J2_SS_NONLIN_HARD_H_ */
