/* $Id: MFPenaltyContact2DT.h,v 1.6 2005-01-28 02:45:12 paklein Exp $ */
#ifndef _MF_PENALTY_CONTACT2D_T_H_
#define _MF_PENALTY_CONTACT2D_T_H_

/* base class */
#include "PenaltyContact2DT.h"

/* direct members */
#include "nMatrixGroupT.h"
#include "VariArrayT.h"
#include "InverseMapT.h"
#include "RaggedArray2DT.h"

namespace Tahoe {

/* forward declarations */
class MeshFreeSupportT;
class SCNIMFT;

/** penalty-based striker-on-facet formulation for meshfree striker nodes */
class MFPenaltyContact2DT: public PenaltyContact2DT
{
public:

	/** constructor */
	MFPenaltyContact2DT(const ElementSupportT& support);

	/** \name writing output */
	/*@{*/	
	/** register element for output */
	virtual void RegisterOutput(void);

	/** write output */
	virtual void WriteOutput(void);
	/*@}*/

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

	/** \name meshfree element group
	 * Depending on the type of meshless formulation, one of the pointers 
	 * MFPenaltyContact2DT::fMeshFreeSupport or MFPenaltyContact2DT::fSCNI
	 * will be non-NULL and will provide access to the meshless shape functions
	 * need to calculate the penetration and distribution of contact forces. */
	/*@{*/
	const ElementBaseT* fElementGroup;

	/** meshfree support from MFPenaltyContact2DT::fElementGroup if the group is
	 * derived from meshless classes using MeshFreeSupportT */
	MeshFreeSupportT* fMeshFreeSupport;
	
	/** meshfree support from MFPenaltyContact2DT::fElementGroup if the group is
	 * based on the smoothed, conforming, nodal integration method */
	const SCNIMFT* fSCNI;
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

	/** \name SCNI data 
	 * one row for each SCNI striker node */
	/*@{*/
	AutoArrayT<int> fSCNI_tmp;
	iArrayT fSCNI_LocalID;
	RaggedArray2DT<int> fSCNI_Support;
	RaggedArray2DT<double> fSCNI_Phi;
	/*@}*/

	/** \name output contact forces */
	/*@{*/
	int fOutputID;
	bool fOutputForce;
	iArray2DT fNodesUsed2D;
	dArray2DT fForce;
	InverseMapT fNodesUsed_inv;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _MF_PENALTY_CONTACT2D_T_H_ */
