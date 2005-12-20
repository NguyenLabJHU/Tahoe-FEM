/* $Id: CurrMixtureSpeciesT.h,v 1.1 2005-12-20 17:24:04 thao Exp $ */
#ifndef _CURR_MIXTURE_SPECIES_T_H_
#define _CURR_MIXTURE_SPECIES_T_H_

/* base class */
#include "NLDiffusionElementT.h"

namespace Tahoe {

class UpdatedLagMixtureT;
class Q1P0MixtureT;

/** class to handle transport with component of mixture */
class CurrMixtureSpeciesT: public NLDiffusionElementT
{
public:
	
	/** constructor */
	CurrMixtureSpeciesT(const ElementSupportT& support);

	/** write element output */
	virtual void WriteOutput(void);

	/** the flux velocity */
	const dArray2DT& FluxVelocity(void) const { return fFluxVelocity; };

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);

	/*@}*/

protected:

	/** method used to compute stress gradient */
	enum GradientOptionT {
		kGlobalProjection,
		kElementProjection
	};

	/** concentration enum */
	enum ConcentrationT {
		kReference,
		kCurrent
	};

	enum SpeciesT {
		kSolid,
		kFluid,
		kSolute
	};
	
	/** returns species type*/
	const SpeciesT ReturnSpecies(void) const {return fSpecies;}

	/** allocate and initialize shape function objects */
	virtual void SetShape(void);

	/** allocate and initialize local arrays */
	virtual void SetLocalArrays(void);

	/** compute shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/** \name element loop operations */
	/*@{*/
	/** reset loop */
	virtual void Top(void);
	
	/** advance to next element. \return true if there is another element, 
	 * false otherwise */ 
	virtual bool NextElement(void);
	/*@}*/

	/** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
	/*@{*/
	/** form group contribution to the stiffness matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

	/** form group contribution to the residual */
	virtual void RHSDriver(void);
	/*@}*/

	/** calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

	/** form the element stiffness matrix */
	virtual void FormStiffness(double constK);

	/** compute the mass flux and flux velocities */
	void ComputeMassFlux(void);

	/** compute the flux velocities and their variation with concentration */
//	void ComputeDMassFlux(void);

	/** compute the divergence tensor field given the values at the integration points */
	void ComputeDivergence(const dMatrixT& ip_grad_transform, const ArrayT<dMatrixT>& tensor_ip,
		dArrayT& div) const;
	
	/*project background velocities from ip to nodes*/
	void ProjectV(void);

protected:

	/** method used to compute stress gradient */
	GradientOptionT fGradientOption;

	/** concentration type */
	ConcentrationT fConcentration;

	/*phase type*/
	SpeciesT fSpecies;
	
	/** write total species mass to output */
	bool fOutputMass;

	/** background solid */
	UpdatedLagMixtureT* fUpdatedLagMixture;
	Q1P0MixtureT*       fQ1P0Mixture;

	/** background species */
	CurrMixtureSpeciesT* fBackgroundSpecies;
	
	/** index of the species within the mixture */
	int fIndex;

	/** flux velocities */
	dArray2DT fFluxVelocity;

	/** mass flux */
	dArray2DT fMassFlux;

	/** concentration specific driving force*/
	dArray2DT fDrivingForce;

	/** concentration specific driving force*/
	dArray2DT fDivBackgroundVel;

	/*ip specific stresses*/
	ArrayT<dMatrixT> ftau_ip;

	/** nodal specific stresses */
	dArray2DT ftau_avg;
	
	dArray2DT fv_bg_avg;
	
	/** \name element-by-element stress projection */
	/*@{*/
	/*grad_x to calculate div_x tau from integration point values*/
	ArrayT<dMatrixT> fip_gradient;
	/*@}*/
	
	/** \name work space */
	/*@{*/
	dMatrixT fNEEmat;
	dMatrixT fNSDmat1, fNSDmat2, fNSDmat3;
	/*@}*/

	/** current coords with local ordering */
	LocalArrayT fLocCurrCoords;	
};

} /* namespace Tahoe */

#endif /* _CURR_MIXTURE_SPECIES_T_H_ */
