/* $Id: MFPenaltyContact2DT.h,v 1.1 2003-11-04 01:31:05 paklein Exp $ */
#ifndef _MF_PENALTY_CONTACT2D_T_H_
#define _MF_PENALTY_CONTACT2D_T_H_

/* base class */
#include "PenaltyContact2DT.h"

/* direct members */
#include "nMatrixGroupT.h"
#include "VariArrayT.h"

namespace Tahoe {

/* forward declarations */
class MeshFreeSupportT;

/** penalty-based striker-on-facet formulation for meshfree striker nodes */
class MFPenaltyContact2DT: public PenaltyContact2DT
{
public:

	/** constructor */
	MFPenaltyContact2DT(const ElementSupportT& support, const FieldT& field);

protected:

	/** print element group data */
	virtual void PrintControlData(ostream& out) const;
		 	
	/** construct the effective mass matrix. Not implemeneted. */
	virtual void LHSDriver(GlobalT::SystemTypeT);

	/* construct the residual force vector */
	virtual void RHSDriver(void);

	/** \name steps in setting contact configuration */
	/*@{*/
	/** set contact surfaces and strikers. This implementation is based on
	 * ContactT::EchoConnectivityData, but taking the strikers from the
	 * list of meshfree nodes, so that all meshfree nodes can be potential
	 * strikers. Also, since the connectivities will be oddly-shaped,
	 * they will not be registered with the ModelManagerT. */
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);

	/** set "internal" data. This implementation is bsed on Contact2DT::SetActiveInteractions,
	 * but modified to account for compute the current coordinates using the
	 * meshfree, kernel representation for the nodal displacements. */
	virtual bool SetActiveInteractions(void);

	/** set "external" data to send to FEManager */
	virtual void SetConnectivities(void);
	/*@}*/
		
protected:

	/** \name meshfree element group */
	/*@{*/
	int fGroupNumber;
	const ElementBaseT* fElementGroup;

	/** meshfree support from MFPenaltyContact2DT::fElementGroup */
	const MeshFreeSupportT* fMeshFreeSupport;
	/*@}*/

	/** \name dynamic memory managers */
	/*@{*/
	/** manager for the Contact2DT derivative arrays */
	nMatrixGroupT<double> fdvT_man;

	/** residual vector */
	VariArrayT<double> fRHS_man;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _MF_PENALTY_CONTACT2D_T_H_ */
