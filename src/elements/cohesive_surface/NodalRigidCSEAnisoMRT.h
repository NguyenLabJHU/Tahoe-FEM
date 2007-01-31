/* $Id: NodalRigidCSEAnisoMRT.h,v 1.2 2007-01-31 03:07:06 skyu Exp $ */
#ifndef _NODAL_RIGID_CSE_ANISO_MR_T_H_
#define _NODAL_RIGID_CSE_ANISO_MR_T_H_

/* base class */
#include "CSEAnisoT.h"
#include "DOFElementT.h"

/* direct members */
#include "VariArrayT.h"
#include "nVariArray2DT.h"

namespace Tahoe {

/** forward declarations */
class MR_NodalRP2DT;

/** Cohesive surface elements with rigid constraints */
class NodalRigidCSEAnisoMRT: public CSEAnisoT, public DOFElementT
{
public:

	/** constructor */
	NodalRigidCSEAnisoMRT(const ElementSupportT& support);

	/** destructor */
	virtual ~NodalRigidCSEAnisoMRT(void);

	/** collecting element connectivities for the field */
	virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;

	/** append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** close current time increment */
	virtual void CloseStep(void);

	/** \name implementation of the DOFElementT interface */
	/*@{*/	
	/* returns the array for the DOF tags needed for the current config */
	virtual void SetDOFTags(void);
	virtual iArrayT& DOFTags(int tag_set);

	/** generate nodal connectivities */
	virtual void GenerateElementData(void);

	/** return the contact elements */
	virtual const iArray2DT& DOFConnects(int tag_set) const;

	/** restore the DOF values to the last converged solution */
	virtual void ResetDOF(dArray2DT& DOF, int tag_set) const;

	/** returns 1 if group needs to reconfigure DOF's, else 0 */
	virtual int Reconfigure(void);

	/** restore any state data to the previous converged state. */
	virtual void ResetState(void);
	
	/** the group */
	virtual int Group(void) const { return CSEAnisoT::Group(); }
	/*@}*/

	/** \name writing output */
	/*@{*/
	virtual void RegisterOutput(void);
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values, const iArrayT& e_codes, dArray2DT& e_values);
	/*@}*/

	/** \name restart functions */
	/*@{*/
	/** write restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::ReadRestart implementation. */
	virtual void WriteRestart(ostream& out) const;

	/** read restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::WriteRestart implementation. */
	virtual void ReadRestart(istream& in);
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** constraint status */
	enum ConstraintStatusT {
		  kFree = 0, /**< traction computed from cohesive relations */
		kActive = 1, /**< traction computed from constraint force */
		kFailed = 2, /**< skip calculations because cohesive relations reached failure */
      kUnneeded = 3  /**< not constrained because of applied boundary conditions */
	};

	/* tangent matrix and force vector */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);
	virtual void RHSDriver(void);

	/** write all current element information to the stream */
	virtual void CurrElementInfo(ostream& out) const;

protected:

	/** regularization parameter */
	double fr; 

	/** constraints */
	iArrayT fConstraintXDOFTags;

	/** connectivities for each constraint */
	iArray2DT fXDOFConnectivities;

	/** equations for each constraint */
	iArray2DT fXDOFEqnos;

	/** flags indicating whether the given constraint is active. Dimensions
	 * of the array are: [nel] x [nip*ndof] */
	Array2DT<char> fConstraintStatus;

	/** \name constraint history. Two snapshots of the constraint values
	 * are needed because DOF values need to be recovered differently for:
	 * -# restoring last converged solution following NodalRigidCSEAnisoMRT::ResetState,
	 *    after the solution for the current step has failed. 
	 * -# forwarding the current solution after NodalRigidCSEAnisoMRT::Reconfigure
	 *    changes the number of active constraints */
	/*@{*/
	Array2DT<char> fConstraintStatus_last;
	AutoArrayT<double> fConstraints_last;

	Array2DT<char> fConstraintStatus_n;
	AutoArrayT<double> fConstraints_n;
	/*@}*/

	/** \name node pair information */
	/*@{*/
	dArrayT   fCZNodeAreas; /**< the tributary area associated with every pair */
	iArray2DT fCZNodePairs;
	iArray2DT fCZNodePairPoints; /**< "point" connectivities needed to write output */
	iArray2DT fCZNodePairDispEqnos;
	dArrayT   fCZDirection;
	int fCurrPair; /**< index used for looping through pairs */
	/*@}*/

	/** \name state variable storage arrays. 
	 * arrays have dimensions: [npairs] x [nvar] */
	/*@{*/
	dArray2DT fStateVariables;
	
	/** previous converged solution */
	dArray2DT fStateVariables_n;
	/*@}*/

	/** the cohesive relation */
	MR_NodalRP2DT* fCZRelation;

	/** coordinate transformation */
	dMatrixT Q;
	ArrayT<dMatrixT> dQ;
	
private:

	/** \name dynamic work space managers */
	/*@{*/
	VariArrayT<int> fConstraintXDOFTags_man;
	nVariArray2DT<int> fXDOFConnectivities_man;		
	nVariArray2DT<int> fXDOFEqnos_man;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _NODAL_RIGID_CSE_ANISO_MR_T_H_ */
