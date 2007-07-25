/* $Id: J2SSKStV2DPlaneStress.h,v 1.1 2006-07-21 20:03:19 tdnguye Exp $ */
/* created: paklein (06/18/1997) */
#ifndef _J2_SS_KSTV_2DPlaneStress_H_
#define _J2_SS_KSTV_2DPlaneStress_H_

#include "SSSolidMatT.h"
#include "IsotropicT.h"
#include "HookeanMatT.h"
#include "J2SSC0HardeningT.h"

namespace Tahoe {

class J2SSKStV2DPlaneStress: public SSSolidMatT,
							 public IsotropicT,
							 public HookeanMatT,
							 public J2SSC0HardeningT
{
public:

	/* constructor */
	J2SSKStV2DPlaneStress(void);

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
	virtual double Pressure(void) const { return fStress2D.Trace()/3.0; };
	/*@}*/

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);

	/* returns the number of variables computed for nodal extrapolation
	 * during for element output, ie. internal variables */
	virtual int NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);
	
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	virtual double YieldCondition(const dSymMatrixT& relstress, double alpha) const;
	
protected:

	/* returns elastic strain */
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain,
		const ElementCardT& element, int ip);
			
	/* return the correction to stress vector computed by the mapping the
	 * stress back to the yield surface, if needed */
	virtual const dSymMatrixT& StressCorrection(const dSymMatrixT& trialstrain, ElementCardT& element, int ip);

	/* return the correction to moduli due to plasticity (if any)
	 *
	 * Note: Return mapping occurs during the call to StressCorrection.
	 *       The element passed in is already assumed to carry current
	 *       internal variable values */
	virtual const dMatrixT& ModuliCorrection(const ElementCardT& element, int ip);

	/* return a pointer to a new plastic element object constructed with
	 * the data from element */
	virtual void AllocateElement(ElementCardT& element);

		/* load element data for the specified integration point */
	virtual void LoadData(const ElementCardT& element, int ip);

	virtual int PlasticLoading(const dSymMatrixT& trialstrain,
	ElementCardT& element, int ip);

	virtual void Update(ElementCardT& element);
	/* computes the relative stress corresponding for the given element
	 * and elastic strain.  The functions returns a reference to the
	 * relative stress in fRelStress */
	virtual dSymMatrixT& RelativeStress(const dSymMatrixT& totalstrain, const ElementCardT& element);

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);
	
			/*internal variable accessors*/
	virtual const iArrayT& InternalDOF(void) const;
	virtual const dArrayT& InternalStrainVars(void);
	virtual const dArrayT& InternalStressVars(void);


private:

	/* return values */
	dSymMatrixT	fStress2D;
	dMatrixT	fModulus2D;

	/* element level internal variables */
	dSymMatrixT fPlasticStrain2D; //deviatoric part of the plastic strain
	dSymMatrixT fRelStress2D;      //unit normal to the stress surface
	dSymMatrixT fTrialRelStress2D;      //unit normal to the stress surface
	dSymMatrixT fBeta2D;          //stress surface "center", kinematic hardening

	/* return values */
	dSymMatrixT	fElasticStrain2D;
	dMatrixT	fModuliCorr2D;
	
	/*work space*/
	dArrayT fVec1;
	dArrayT fVec2;
	dMatrixT    fP;
	dMatrixT	fModTemp;
	dMatrixT    fModTheta;
	
	/** number of elastic iterations */
	int fmax_iteration;
};

} // namespace Tahoe 
#endif /* _J2_SS_KSTV_2D_H_ */
