/* $Id: MFLagMultT.h,v 1.1 2004-05-06 18:54:47 cjkimme Exp $ */
#ifndef _MF_AUG_LAG_MULT_T_H_
#define _MF_AUG_LAG_MULT_T_H_

/* base classes */
#include "FBC_ControllerT.h"
#include "DOFElementT.h"

/* direct members */
#include "ElementMatrixT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "ScheduleT.h"
#include "KBC_CardT.h"


namespace Tahoe {

/* forward declarations */
class XDOF_ManagerT;
class FieldT;
class StringT;
class SCNIMFT;

/** rigid barrier with augmented Lagrangian enforcement of
 * non-interpenetration */
class MFAugLagMultT: public FBC_ControllerT, public DOFElementT
{
public:

	/* constructor */
	MFAugLagMultT(FEManagerT& fe_manager, XDOF_ManagerT* XDOF_nodes, const FieldT& field,
		const dArray2DT& coords, const dArray2DT& disp);

	/* input processing */
	virtual void EchoData(ifstreamT& in, ostream& out);

	/* initialize data */
	virtual void Initialize(void);
	virtual void SetEquationNumbers(void);

	/* append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
	virtual void Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
		AutoArrayT<const iArray2DT*>& equivalent_nodes) const;

	/* restarts */
	virtual void ReadRestart(istream& in);
	virtual void WriteRestart(ostream& out) const;

	/* finalize step */
	virtual void CloseStep(void);
	
	/* apply force */
	virtual void ApplyRHS(void);

	/** tangent
	 * \param sys_type "maximum" tangent type needed by the solver. The GlobalT::SystemTypeT
	 *        enum is ordered by generality. The solver should indicate the most general
	 *        system type that is actually needed. */
	virtual void ApplyLHS(GlobalT::SystemTypeT sys_type);

	/* returns the array for the DOF tags needed for the current config */
	virtual void SetDOFTags(void);
	virtual iArrayT& DOFTags(int tag_set);

	/* generate nodal connectivities - does nothing here */
	virtual void GenerateElementData(void);

	/* return the contact elements */
	virtual const iArray2DT& DOFConnects(int tag_set) const;

	/* restore the DOF values to the last converged solution */
	virtual void ResetDOF(dArray2DT& DOF, int tag_set) const;

	/* returns 1 if group needs to reconfigure DOF's, else 0 */
	virtual int Reconfigure(void);

	/** restore any state data to the previous converged state */
	virtual void ResetState(void) { };

	/** return the equation group to which the generate degrees of
	 * freedom belong. */
	virtual int Group(void) const;	

private:

	/* accumulate the constraint force vector fConstraintForce */
	virtual void ComputeConstraintValues(double kforce);

private:

	/** nodemanager */
	XDOF_ManagerT* fXDOF_Nodes;
	
	/** the field */
	const FieldT& fField;
	
	int fNumConstrainedDOFs;
	
	/** \name references to NodeManagerT data */
	/*@{*/
	const iArray2DT& rEqnos;  /**< nodal equation numbers */
	const dArray2DT& rCoords; /**< nodal coordinates */
	const dArray2DT& rDisp;   /**< nodal displacement */
	/*@}*/
	
	/** \name contact force node and equation numbers */
	/*@{*/
	iArrayT fConstraintNodes;
	iArrayT fConstraintEqnos;

	/** shallow version of PenaltyRegionT::fContactForce2D */
	dArrayT fConstraintForce;
	dArray2DT fConstraintForce2D;
	
	/** contact equation sets (shallow copy of contact node equs) */
	iArray2DT fConstraintEqnos2D;
	iArray2DT fConstraintTags;
	/*@}*/
	
	/** \name Augmented multiplier info */
	/*@{*/
	iArrayT fConstraintDOFtags; /**< constrained DOF tags and DOF's */
	iArrayT fFloatingDOF;    /**< 1 if multiplier is attacted to node that has KBC's */
	dArrayT fLastDOF;        /**< multiplier history */ 
	/*@}*/
	
	/** \name communication data with meshfree elements */
	/*@{*/
	int fBlockID;
	/*@}*/
	
	/** penalty stiffness */
	double fk;
	
	/* workspace */
	ElementMatrixT fLHS;  //tangent matrix
	
	/** \name storage for Kinematic boundary condition constraints */
	/*@{*/
	iArrayT fConstrainedDOFs, fScheduleNums, fNumConstraints;
	ArrayT<StringT> fNodeSetIDs;
	dArrayT fScales;
	ArrayT<KBC_CardT::CodeT> fCodes;
	/*@}*/
	
	/** values of constrained displacements */
	dArrayT fConstraintValues;
	
	SCNIMFT* mfElemGroup;
};

} // namespace Tahoe 
#endif /* _MF_AUG_LAG_MULT_T_H_ */
