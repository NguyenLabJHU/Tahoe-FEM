/* $Id: BridgingScaleT.h,v 1.22 2003-01-29 07:34:28 paklein Exp $ */
#ifndef _BRIDGING_SCALE_T_H_
#define _BRIDGING_SCALE_T_H_

/* direct members */
#include "SolidElementT.h"
#include "RaggedArray2DT.h"
#include "ElementMatrixT.h"
#include "CCSMatrixT.h"

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
		const SolidElementT& solid);

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
	virtual void WriteOutput(void);

	/** form of tangent matrix, symmetric by default */
	virtual GlobalT::SystemTypeT TangentType(void) const { return GlobalT::kSymmetric; };

	/** accumulate the residual force on the specified node
	 * \param node test node
	 * \param force array into which to assemble to the residual force */
	virtual void AddNodalForce(const FieldT&, int, dArrayT&) {};

	/** returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void) { return 0; };

	/** compute specified output parameter and send for smoothing */
	virtual void SendOutput(int) {};

protected:

	/** echo element connectivity data. No connectivities need to be read */
	virtual void EchoConnectivityData(ifstreamT&, ostream&) {};

	/* called by FormRHS and FormLHS */
	virtual void LHSDriver(GlobalT::SystemTypeT) {};
	virtual void RHSDriver(void) {};

	/** allocate and initialize local arrays */
	virtual void SetLocalArrays(void);

	/** write element group parameters to out */
	virtual void PrintControlData(ostream& out) const;

	/** write all current element information to the stream. used to generate
	 * debugging information after runtime errors */
	virtual void CurrElementInfo(ostream& out) const;

private:

	/** \name bridging scale-related computational functions */
	/*@{*/
	/** compute coarse and fine scale displacement fields */
	void CoarseFineFields(void);
	/*@}*/

protected:

	/** "element" group calculating particle solution */
	const RodT& fParticle;
	iArray2DT fParticlesUsed;
	
	/** continuum group solving displacements */
	const SolidElementT& fSolid;
	iArray2DT fSolidNodesUsed;

	/** list of particles per element: [n_cell] x [n_part_i] */
	RaggedArray2DT<int> fParticlesInCell;
	
	/** take fParticlesInCell, now have list of inverse mappings per element:
	 *  [n_cell] x [n_inversemap_i] */
	RaggedArray2DT<double> fInverseMapInCell;

	/** Nodal degrees of freedom */
	dArrayT fUx, fUy;
	dMatrixT fWtempU;
	dArray2DT fFineScaleU;

	int fTotalNodes;
	iArray2DT fConnect, fAtomConnect;
	ElementMatrixT fElMatU;
	CCSMatrixT fGlobalMass;

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

	/** \name output ID's */
	/*@{*/
	int fParticleOutputID;
	int fSolidOutputID;
	/*@}*/
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
