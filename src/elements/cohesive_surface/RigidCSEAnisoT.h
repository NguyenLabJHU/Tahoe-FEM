/* $Id: RigidCSEAnisoT.h,v 1.3 2006-05-21 17:47:59 paklein Exp $ */
#ifndef _RIGID_CSE_ANISO_T_H_
#define _RIGID_CSE_ANISO_T_H_

/* base class */
#include "CSEAnisoT.h"
#include "DOFElementT.h"

/* direct members */
#include "VariArrayT.h"
#include "nVariArray2DT.h"

namespace Tahoe {

/** Cohesive surface elements with rigid constraints */
class RigidCSEAnisoT: public CSEAnisoT, public DOFElementT
{
public:

	/** constructor */
	RigidCSEAnisoT(const ElementSupportT& support);

	/** destructor */
	virtual ~RigidCSEAnisoT(void);

	/** collecting element connectivities for the field */
	virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;

	/** append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** close current time increment */
	virtual void CloseStep(void);

	/** restore the element group to its state at the beginning of the
	 * current time step. */
	virtual GlobalT::RelaxCodeT ResetStep(void); 

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

	/** restore any state data to the previous converged state */
	virtual void ResetState(void) { };

	/** the group */
	virtual int Group(void) const { return CSEAnisoT::Group(); }
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
		kFree = 0,
		kActive = 1
	};

	/* tangent matrix and force vector */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);
	virtual void RHSDriver(void);

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

	/** \name history */
	/*@{*/
	Array2DT<char> fConstraintStatus_last;
	AutoArrayT<double> fConstraints_last;
	/*@}*/

	/** \name local arrays */
	/*@{*/
	LocalArrayT fDisp;
	LocalArrayT fLastDisp;
	/*@}*/

	/** shape functions for the constraints */
	SurfaceShapeT* fConstraintShapes;

private:

	/** \name dynamic work space managers */
	/*@{*/
	VariArrayT<int> fConstraintXDOFTags_man;
	nVariArray2DT<int> fXDOFConnectivities_man;		
	nVariArray2DT<int> fXDOFEqnos_man;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _RIGID_CSE_ANISO_T_H_ */
