/* $Id: ViscousDragT.h,v 1.1.2.1 2003-09-10 13:35:27 paklein Exp $ */
#ifndef _VISCOUS_DRAG_T_H_
#define _VISCOUS_DRAG_T_H_

#include "ElementBaseT.h"

namespace Tahoe {

/** viscous drag. Apply drag to the nodes in a given element block
 * that's proportional to the nodal velocity and weighted by the
 * mass associated with the node. */
class ViscousDragT: public ElementBaseT
{
public:

	/** constructor */
	ViscousDragT(const ElementSupportT& support, const FieldT& field);

	/** class initialization */
	virtual void Initialize(void);

	/** collecting element group equation numbers */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
	
protected:

	/** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
	/*@{*/
	/** form group contribution to the stiffness matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

	/** form group contribution to the residual */
	virtual void RHSDriver(void);
	/*@}*/

	/** override to disable */
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out) { };

private:

	/** viscosity per unit volume */
	double fViscosity;
	
	/** element block ID */
	StringT fID;
	
	/** nodes used */
	iArrayT fNodesUsed;

	/** nodal mass */
	dArrayT fNodalMass;

	/** drag force */
	dArray2DT fDragForce;
	
	/** equations */
	iArray2DT fEqnos;
};

} /* namespace Tahoe */

#endif /* _VISCOUS_DRAG_T_H_ */
