/* $Id: FSSolidMixtureT.h,v 1.8 2005-05-10 17:57:07 paklein Exp $ */
#ifndef _FS_SOLID_MIX_T_H_
#define _FS_SOLID_MIX_T_H_

/* base class */
#include "FSSolidMatT.h"

namespace Tahoe {

/* forward declarations */
//class FSSolidMixtureSupportT;

/** base class for finite deformation solid composed of a mixture */
class FSSolidMixtureT: public FSSolidMatT
{
public:

	/** constructor */
	FSSolidMixtureT(void);

	/** destructor */
	virtual ~FSSolidMixtureT(void);

	/** form of tangent matrix. \return symmetric by default */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/** set the material support or pass NULL to clear */
//	virtual void SetFSSolidMixtureSupport(const FSSolidMixtureSupportT* support);

	/** finite strain mixture materials support */
//	const FSSolidMixtureSupportT& FSSolidMixtureSupport(void) const;

	/** \name concentration */
	/*@{*/
	/** concentration enum */
	enum ConcentrationT {
		kReference,
		kCurrent
	};
	void SetConcentration(int i, ConcentrationT conc);
	ConcentrationT Concentration(int i) const { return fConcentration[i]; };
	/*@}*/

	/** need to compute objective velocity gradient in mass balance
	 * \note this is really only needed when FSSolidMixtureT::fConcentration ==
	 *       FSSolidMixtureT::kCurrent; however, the element class configures
	 *       itself based on the state of the material after construction,
	 *       while the concentration type of each associated MixtureSpeciesT
	 *       isn't determined until later. */
	virtual bool Need_F_last(void) const { return true; };

	/** get all nodal concentrations over the current element */
	void UpdateConcentrations(void);

	/** update the specific nodal concentrations over the current element */
	void UpdateConcentrations(int i);

	/** return the index of the species associated with the given field name, or
	 * -1 if the field is not found */
	int SpeciesIndex(const StringT& field_name) const;

	/** density varies with position */
	virtual bool HasChangingDensity(void) const { return true; };

	/** mass density. Method does retrieve current values of the nodal concentrations. */
	virtual double Density(void);
	
	/** mass density of the given species */
	double Density(int i);

	/** variation of the Cauchy for the given species with concentration */
	const dSymMatrixT& ds_ij_dc(int i);

	/** \name spatial representation */
	/*@{*/
	/** strain energy density. Method does retrieve current values of the nodal concentrations. */
	virtual double StrainEnergyDensity(void);
	
	/** total material tangent modulus. Method does retrieve current values of the nodal concentrations. */
	virtual const dMatrixT& c_ijkl(void);

	/** partial material tangent modulus. Method does not retrieve current values of the nodal 
	 * concentrations. These can be updated with FSSolidMixtureT::UpdateConcentrations. */
	const dMatrixT& c_ijkl(int i);

	/** total Cauchy stress. Method does retrieve current values of the nodal concentrations. */
	virtual const dSymMatrixT& s_ij(void);

	/** partial Cauchy stress. Method does not retrieve current values of the nodal 
	 * concentrations. These can be updated with FSSolidMixtureT::UpdateConcentrations. */
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

	/** \name material output variables.
	 * See ContinuumMaterialT for more documentation */
	/*@{*/
	/** return the number of constitutive model output values. The material
	 * returns the partial stresses. */
	virtual int NumOutputVariables(void) const;

	/** return the labels for model output parameters */
	virtual void OutputLabels(Tahoe::ArrayT<StringT>& labels) const;

	/** return material output variables */
	virtual void ComputeOutput(Tahoe::dArrayT& output);
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
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

	/** compute integration point concentrations. Interpolate the nodal values
	 * to the current integration point and convert all concentrations to
	 * reference concentration if FSSolidMixtureT::Concentration is
	 * FSSolidMixtureT::kCurrent.
	 * \param c_nodal nodal concentrations
	 * \param c_ip returns with reference concentrations at the current
	 *        integration point. */
	void IPConcentration(const LocalArrayT& c_nodal, dArrayT& c_ip) const;

protected:

	/** support for finite strain mixture materials */
//	const FSSolidMixtureSupportT* fFSSolidMixtureSupport;

	/** element concentration. Nodal values of the species concentrations. The array can
	 * contain either reference or current concentrations. The type is indicated in
	 * FSSolidMixtureT::fConcentration. */
	LocalArrayT fConc;
	
	/** integration point concentrations. If updated using FSSolidMixtureT::IPConcentration,
	 * all concentrations will be with respect to the reference configuration. */
	dArrayT fIPConc;

	/** concentration type for each species */
	ArrayT<ConcentrationT> fConcentration;
	bool fHasCurrent;

	/** array of stored energy functions */
	ArrayT<FSSolidMatT*> fStressFunctions;

	/** concentration field for each stress function */
	ArrayT<const FieldT*> fFields;
	
	/** support for stress functions */
	FSMatSupportT* fStressSupport;

	ArrayT<dMatrixT> fF_species;
	dMatrixT fF_growth_inv;
	dSymMatrixT fs_ij_tmp;
	dSymMatrixT fI;
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
