/* $Id: BridgingScaleT.h,v 1.3 2002-07-18 17:45:23 paklein Exp $ */
#ifndef _BRIDGING_SCALE_T_H_
#define _BRIDGING_SCALE_T_H_

/* direct members */
#include "ElasticT.h"
#include "RaggedArray2DT.h"

namespace Tahoe {

/* forward declarations */
class RodT;

/** base class for elements using shape functions */
class BridgingScaleT: public ElementBaseT
{
public:

	/** constructor */
	BridgingScaleT(const ElementSupportT& support, const FieldT& field,
		const RodT& particle,
		const ElasticT& solid);

	/** destructor */
	virtual ~BridgingScaleT(void);

	/* accessors */
	const ShapeFunctionT& ShapeFunction(void) const;

	/** element coordinates.
	 * \return initial nodal coordinates of current element: [nen] x [nsd] */
	const LocalArrayT& InitialCoordinates() const;
	
	/** element displacements.
	 * \return nodal displacements of current element: [nen] x [ndof] */
	const LocalArrayT& Displacements() const;

	/** initialization. called immediately after constructor */
	virtual void Initialize(void);

	/** collecting element group equation numbers. This call from the FEManagerT
	 * is a signal to the element group that the equation system is up to date
	 * for the current time increment. See ElementBaseT::Equations for more
	 * information. */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/* initialize/finalize time increment */
	virtual void InitStep(void);
	virtual void CloseStep(void);
	virtual void ResetStep(void); // restore last converged state

	/** register self for output */
	virtual void RegisterOutput(void);

	/** send output */
	virtual void WriteOutput(IOBaseT::OutputModeT mode);

protected:

	/** allocate and initialize local arrays */
	virtual void SetLocalArrays(void);

	/** form the residual force vector. computes contribution from natural
	 * boundary conditions */
	virtual void RHSDriver(void);

	/** write element group parameters to out */
	virtual void PrintControlData(ostream& out) const;
	
	/* element data */
	virtual void EchoOutputCodes(ifstreamT& in, ostream& out) = 0;

	/** write all current element information to the stream. used to generate
	 * debugging information after runtime errors */
	virtual void CurrElementInfo(ostream& out) const;

	/* driver for calculating output values */
	virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const = 0;
	virtual void SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const = 0;
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	                           const iArrayT& e_codes, dArray2DT& e_values) = 0;


private:

	/* construct output labels array */
	virtual void GenerateOutputLabels(
		const iArrayT& n_codes, ArrayT<StringT>& n_labels, 
		const iArrayT& e_codes, ArrayT<StringT>& e_labels) const = 0;

	/* bridging scale-related computational functions */
	
	/* computes error caused by projecting solution onto FEM basis space */
	void ComputeError(void);

	/* computes "fine scale" displacement, ie MD - overlap due to FEM */
	void ComputeFineScaleU(void);

protected:

	/** "element" group calculating particle solution */
	const RodT& fParticle;
	
	/** continuum group solving displacements */
	const ElasticT& fSolid; 

	/** list of particles per element: [n_cell] x [n_part_i] */
	RaggedArray2DT<int> fParticlesInCell;

	/* output control */
	iArrayT	fNodalOutputCodes;
	iArrayT	fElementOutputCodes;
	
	/* arrays with local ordering */
	LocalArrayT	fLocInitCoords;	/**< initial coords with local ordering */
	LocalArrayT fLocDisp;	    /**< displacements with local ordering  */ 
	
	/* work space */
	dArrayT fNEEvec; /**< work space vector: [element DOF] */
	dArrayT fDOFvec; /**< work space vector: [nodal DOF]   */
//	dArrayT fNSDvec; /**< work space vector: [nodal dim]   */

private:

	/** output ID */
	int fOutputID;
};

/* inlines */

/* accessors */
inline const ShapeFunctionT& BridgingScaleT::ShapeFunction(void) const
{
	return fSolid.ShapeFunction();
}

inline const LocalArrayT& BridgingScaleT::InitialCoordinates() const
{	
	return fLocInitCoords;
}

inline const LocalArrayT& BridgingScaleT::Displacements() const
{
	return fLocDisp;
}

} // namespace Tahoe 
#endif /* _BRIDGING_SCALE_T_H_ */
