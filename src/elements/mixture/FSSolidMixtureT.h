/* $Id: FSSolidMixtureT.h,v 1.1 2004-10-21 18:48:44 paklein Exp $ */
#ifndef _FS_SOLID_MIX_T_H_
#define _FS_SOLID_MIX_T_H_

/* base class */
#include "FSSolidMatT.h"

namespace Tahoe {

/* forward declarations */
class FSSolidMixtureSupportT;

/** base class for finite deformation solid composed of a mixture */
class FSSolidMixtureT: public FSSolidMatT
{
public:

	/** constructor */
	FSSolidMixtureT(void);

	/** destructor */
	virtual ~FSSolidMixtureT(void);

	/** set the material support or pass NULL to clear */
//	virtual void SetFSSolidMixtureSupport(const FSSolidMixtureSupportT* support);

	/** finite strain mixture materials support */
//	const FSSolidMixtureSupportT& FSSolidMixtureSupport(void) const;

	/** get nodal concentrations over the current element */
	void UpdateConcentrations(void);

	/** \return mass density */
	virtual double Density(void);

	/** \name spatial representation */
	/*@{*/
	/** strain energy density */
	virtual double StrainEnergyDensity(void);
	
	/** total material tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** partial material tangent modulus */
	const dMatrixT& c_ijkl(int i);

	/** total Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/** partial Cauchy stress */
	const dSymMatrixT& s_ij(int i);

	/** return the pressure associated with the last call to 
	 * FSSolidMixtureT::s_ij. See SolidMaterialT::Pressure
	 * for more information. */
	virtual double Pressure(void) const { return fStress.Trace()/3.0; };
	/*@}*/

	/** \name history variables, see ContinuumMaterialT for more information */
	/*@{*/
	/** has history variables */
	virtual bool HasHistory(void) const { return true; };
	
	/** all integration points need to be initialized */
	virtual bool NeedsPointInitialization(void) const { return true; };
	
	/** initialization. Called per integration point for every
	 * -# reference concentration */
	virtual void PointInitialize(void);
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** return the specified stress function or NULL */
	FSSolidMatT* New(const StringT& name) const;

protected:

	/** support for finite strain mixture materials */
//	const FSSolidMixtureSupportT* fFSSolidMixtureSupport;

	/** element concentration */
	LocalArrayT fConc;
	
	/** integration point concentrations */
	dArrayT fIPConc;

	/** array of stored energy functions */
	ArrayT<FSSolidMatT*> fStressFunctions;
	
	/** support for stress functions */
	FSMatSupportT* fStressSupport;

	ArrayT<dMatrixT> fF_species;
	dMatrixT fF_growth_inv;
};

#if 0
/* finite strain materials support */
inline const FSSolidMixtureSupportT& FSSolidMixtureT::FSSolidMixtureSupport(void) const
{ 
#if __option(extended_errorcheck)
	if (!fFSSolidMixtureSupport) 
		ExceptionT::GeneralFail("FSSolidMixtureT::FSSolidMixtureSupport", "pointer not set");
#endif
	return *fFSSolidMixtureSupport; 
}
#endif

} /* namespace Tahoe */

#endif /* _FS_SOLID_MIX_T_H_ */
