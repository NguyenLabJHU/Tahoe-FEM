/* $Id: CurrMixtureSpeciesT.h,v 1.6 2006-01-12 01:51:14 thao Exp $ */
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

	/** compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);

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
	void ComputeMassFlux(bool compute_dmass_flux);
	
	/*project background velocities from ip to nodes*/
	void ProjectV(void);

	/** \name variation of divergence with respect to the nodal values */
	/*@{*/
	/** using element-by-element projection */
	void ComputeDDivergence(const dMatrixT& ip_grad_transform, const ArrayT<dMatrixT>& tensor_ip,
		dMatrixT& d_div) const;

	void ComputeDDivergence(const ArrayT<dMatrixT>& tensor_ip, dMatrixT& d_div) const;
	/*@}*/

	/** driver for calculating output values */
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	                           const iArrayT& e_codes, dArray2DT& e_values);

private:

	/** \name construct output labels array */
	/*@{*/
	virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void GenerateOutputLabels(const iArrayT& n_counts,
		ArrayT<StringT>& n_labels, const iArrayT& e_counts, ArrayT<StringT>& e_labels) const;
	/*@}*/

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

	/** variation in mass flux with concentration */
	dArray2DT fDMassFlux;

	/** concentration specific driving force*/
	dArray2DT fDrivingForce;

	/** concentration specific driving force*/
	dArray2DT fDivBackgroundVel;

	/*ip specific stresses*/
	ArrayT<dMatrixT> fcauchy_ip;
    ArrayT<dMatrixT> fdcauchy_ip;

	/** nodal specific stresses */
	dArray2DT fcauchy_avg;
	
	dArray2DT fv_bg_avg;
	
	/** \name element-by-element stress projection */
	/*@{*/
	/*grad_x to calculate div_x cauchy from integration point values*/
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
