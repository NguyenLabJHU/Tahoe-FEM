/* $Id: MFPenaltyContact2DT.h,v 1.2.32.4 2004-07-13 16:42:29 paklein Exp $ */
#ifndef _MF_PENALTY_CONTACT2D_T_H_
#define _MF_PENALTY_CONTACT2D_T_H_

/* base class */
#include "PenaltyContact2DT.h"

/* direct members */
#include "nMatrixGroupT.h"
#include "VariArrayT.h"
#include "InverseMapT.h"

namespace Tahoe {

/* forward declarations */
class MeshFreeSupportT;

/** penalty-based striker-on-facet formulation for meshfree striker nodes */
class MFPenaltyContact2DT: public PenaltyContact2DT
{
public:

	/** constructor */
	MFPenaltyContact2DT(const ElementSupportT& support);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** construct the effective mass matrix. Not implemeneted. */
	virtual void LHSDriver(GlobalT::SystemTypeT);

	/* construct the residual force vector */
	virtual void RHSDriver(void);

	/** \name steps in setting contact configuration */
	/*@{*/
	/** Echo contact bodies and striker nodes. After the read section, should 
	 * have valid nodes/facet connectivities for the local database. */
	virtual void ExtractContactGeometry(const ParameterListT& list);

	/** set "internal" data. This implementation is bsed on Contact2DT::SetActiveInteractions,
	 * but modified to account for compute the current coordinates using the
	 * meshfree, kernel representation for the nodal displacements. */
	virtual bool SetActiveInteractions(void);
	/*@}*/

	/** compute the current coordinates of the given meshless striker nodes.
	 * Current striker coordinates are written into ContactT::fStrikerCoords. */
	void ComputeStrikerCoordinates(const ArrayT<int>& strikers);

	/** set derivative arrays given the array of shape functions for the
	 * nodes in the neighborhood of the meshfree striker. */
	void SetDerivativeArrays(const dArrayT& mf_shape);
	
protected:

	/** \name meshfree element group */
	/*@{*/
	const ElementBaseT* fElementGroup;

	/** meshfree support from MFPenaltyContact2DT::fElementGroup */
	MeshFreeSupportT* fMeshFreeSupport;
	/*@}*/

	/** map from global ID to meshfree node index */
	InverseMapT fNodeToMeshFreePoint;

	/** map from global ID to active striker index */
	InverseMapT fNodeToActiveStriker;

	/** \name dynamic memory managers */
	/*@{*/
	/** striker coordinate work space */
	nVariArray2DT<double> fStrikerCoords_man;
	
	/** manager for the Contact2DT derivative arrays */
	nMatrixGroupT<double> fdvT_man;

	/** residual vector */
	VariArrayT<double> fRHS_man;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _MF_PENALTY_CONTACT2D_T_H_ */
