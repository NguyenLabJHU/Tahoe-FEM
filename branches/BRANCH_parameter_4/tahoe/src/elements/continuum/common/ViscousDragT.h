/* $Id: ViscousDragT.h,v 1.2.26.2 2004-07-13 16:42:30 paklein Exp $ */
#ifndef _VISCOUS_DRAG_T_H_
#define _VISCOUS_DRAG_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "InverseMapT.h"

namespace Tahoe {

/** viscous drag. Apply drag to the nodes in a given element block
 * that's proportional to the nodal velocity and weighted by the
 * mass associated with the node. */
class ViscousDragT: public ElementBaseT
{
public:

	/** constructor */
	ViscousDragT(const ElementSupportT& support);

	/** form of tangent matrix, symmetric by default */
	virtual GlobalT::SystemTypeT TangentType(void) const { return GlobalT::kDiagonal; };

	/** collecting element group equation numbers */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** accumulate the residual force on the specified node */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);

	/** returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void) { return 0.0; };

	/** \name writing output */
	/*@{*/	
	/** register element for output. An interface to indicate the element group
	 * must create an OutputSetT and register it with FEManagerT::RegisterOutput
	 * to obtain an output ID that is used to write data to the current
	 * output destination. */
	virtual void RegisterOutput(void);

	/** write element output. An interface to indicate the element group
	 * gather nodal and element data and send it for output with
	 * FEManagerT::WriteOutput */
	virtual void WriteOutput(void);

	/** compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);
	/*@}*/
	
protected:

	/** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
	/*@{*/
	/** form group contribution to the stiffness matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

	/** form group contribution to the residual */
	virtual void RHSDriver(void);
	/*@}*/

	/** override to disable */
	virtual void EchoConnectivityData(ifstreamT&, ostream&) { };

private:

	/** viscosity per unit volume */
	double fViscosity;
	
	/** element block ID */
	StringT fID;
	
	/** nodes used */
	iArray2DT fNodesUsed;

	/** nodal mass */
	dArrayT fNodalMass;

	/** drag force */
	dArray2DT fDragForce;
	
	/** equations */
	iArray2DT fEqnos;
	
	/** map for global number to local index */
	InverseMapT fGlobalToLocal;
	
	/** output ID */
	int fOutputID;
	
	/** incremental viscous dissipation for the entire group */
	double fIncrementalDissipation;
};

} /* namespace Tahoe */

#endif /* _VISCOUS_DRAG_T_H_ */
