/* $Id: MFAugLagMultT.h,v 1.6 2005-04-12 15:34:07 paklein Exp $ */
#ifndef _MF_AUG_LAG_MULT_T_H_
#define _MF_AUG_LAG_MULT_T_H_

#include "ElementsConfig.h"
#ifdef CONTINUUM_ELEMENT

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
#include "RaggedArray2DT.h"
#include "nVariMatrixT.h"
#include "VariArrayT.h"

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

	/** constructor */
	MFAugLagMultT(void);

	/** initialize data */
	virtual void SetEquationNumbers(void);

	/* append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
	virtual void Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
		AutoArrayT<const iArray2DT*>& equivalent_nodes) const;

	/* initial condition/restart functions
	 *
	 * Set to initial conditions.  The restart functions
	 * should read/write any data that overrides the default
	 * values */
	virtual void InitialCondition(void);
	virtual void ReadRestart(istream& in);
	virtual void WriteRestart(ostream& out) const;

	/* initialize/finalize step */
	virtual void InitStep(void);
	virtual void CloseStep(void);
	
	/* apply force */
	virtual void ApplyRHS(void);
	
	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/** tangent
	 * \param sys_type "maximum" tangent type needed by the solver. The GlobalT::SystemTypeT
	 *        enum is ordered by generality. The solver should indicate the most general
	 *        system type that is actually needed. */
	virtual void ApplyLHS(GlobalT::SystemTypeT sys_type);

	/* reset to the last known solution */
	virtual void Reset(void);

	/** \name writing results */
	/*@{*/
	/** register data for output */
	virtual void RegisterOutput(void);

	/** write results */
	virtual void WriteOutput(ostream& out) const;
	/*@}*/

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
	virtual void ResetState(void) {};

	/** return the equation group to which the generate degrees of
	 * freedom belong. */
	virtual int Group(void) const;	

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/* accumulate the constraint force vector fConstraintForce */
	virtual void ComputeConstraintValues(double kforce);
	
	/** initialize data structures and communicate with meshfree elements */
	void ChatWithElementGroup(void);

protected:

	/** nodemanager */
//	XDOF_ManagerT* fXDOF_Nodes;
	
	int fNumConstrainedDOFs; // This the number of LaGrange multipliers
	
	/** \name references to NodeManagerT data */
	/*@{*/
//	const iArray2DT& rEqnos;  /**< nodal equation numbers */
//	const dArray2DT& rDisp;   /**< nodal displacement */
	/*@}*/
	
	/** \name contact force node and equation numbers */
	/*@{*/
	//iArrayT fConstraintNodes, fLocallyNumberedConstraintNodes;
	iArrayT fConstraintEqnos;
	int fNumEqs;

	/** shallow version of PenaltyRegionT::fContactForce2D */
	dArrayT fConstraintForce;
	
	/** constraint equation sets (shallow copy of contact node equs) */
	iArray2DT fConstraintTags;
	iArray2DT fConstraintEqnos2D;
	iArray2DT fFlatEqNos;
	/*@}*/
	
	/** \name Augmented multiplier info */
	/*@{*/
	iArrayT fConstraintDOFtags; /**< constrained DOF tags and DOF's */
	dArrayT fLastDOF;        /**< multiplier history */ 
	/*@}*/
	
	/** \name communication data with meshfree elements */
	/*@{*/
	int fBlockID;
	/*@}*/
	
	/** penalty stiffness */
	double fk;
	
	/* RHS workspace */
	dArrayT fRHS;
	VariArrayT<double> fRHS_wrapper;
	
	/* LHS workspaces */
	nVariMatrixT<double> fLHS_wrapper;
	ElementMatrixT fLHS; 
	VariArrayT<int> fRowEqs_wrapper;
	iArrayT fRowEqs;
	nVariMatrixT<double> fOtherLHS_wrapper;
	ElementMatrixT fOtherLHS;
	
	/** \name storage for Kinematic boundary condition constraints */
	/*@{*/
	iArrayT fConstrainedDOFs, fScheduleNums;
	ArrayT<StringT> fNodeSetIDs;
	RaggedArray2DT<int> fNodeSets, fEqNos;
	RaggedArray2DT<int> fLocallyNumberedNodeSets;
	iArray2DT fFlattenedNodeSets;
	iArrayT fLocalFlatNodes;
	iArrayT fUnionOfNodes;
	iArrayT fKey;
	dArrayT fScales;
	ArrayT<KBC_CardT::CodeT> fCodes;
	/*@}*/
	
	/** a flattened version of fNodeSets containing the size of each node's support */
	iArrayT fSupportSizes;
	
	/** values of constrained displacements */
	dArrayT fConstraintValues;
	
	/** the group furnishing MF shape functions for the contrained nodes */
	SCNIMFT* mfElemGroup;
	
	RaggedArray2DT<int> fsupport;
	RaggedArray2DT<double> fphi;
	
	/** \name writing results */
	/*@{*/	
	/** output ID */
	int fOutputID;
	
	/** "connectivities" for output, just alias of PenaltyRegionT::fContactNodes */
	
	/*@}*/	
};

} // namespace Tahoe 

#endif /* CONTINUUM_ELEMENT */

#endif /* _MF_AUG_LAG_MULT_T_H_ */
