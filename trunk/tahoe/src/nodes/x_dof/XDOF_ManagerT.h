/* $Id: XDOF_ManagerT.h,v 1.2 2001-08-15 18:36:21 paklein Exp $ */
/* created: paklein (06/01/1998) */

#ifndef _XDOF_MANAGER_T_H_
#define _XDOF_MANAGER_T_H_

/* direct members */
#include "dArrayT.h"
#include "AutoArrayT.h"
#include "VariArrayT.h"

/* forward declarations */
class DOFElementT;
class iArray2DT;
class dArray2DT;
class iArrayT;

/** base class for manager of degrees of freedom requested
 * by DOFElementT's. Element groups must be derived from the
 * DOFElementT class and must register themselves using the
 * XDOF_ManagerT::Register funciton. */
class XDOF_ManagerT
{
public:
	
	/** constructor */
	XDOF_ManagerT(void);	

	/** destructor */
	virtual ~XDOF_ManagerT(void);

	/** add element group to list. Each group can request an arbitrary
	 * number of tag sets. Each set can have a different number of DOF's
	 * per tag requested. Requests for tags are collected by the XDOF_ManagerT
	 * using the DOFElementT interface.
	 * \param group pointer to the DOFElementT class. Each group should
	 *        only register once.
	 * \param numDOF array of the number of degrees of freedom per tag 
	 *        in each set of tags. The length of the array is the number
	 *        of tag sets the group requires. */
	virtual void Register(DOFElementT* group, const iArrayT& numDOF);
	
	/** get equations of the element DOF's.
	 * \param group pointer to the DOFElementT requesting equation numbers
	 * \param tag_set set number with the DOFElementT */
	virtual const iArray2DT& XDOF_Eqnos(const DOFElementT* group, int tag_set) const;

	/** get values of the element DOF's.
	 * \param group pointer to the DOFElementT requesting DOF values
	 * \param tag_set set number with the DOFElementT */
	virtual const dArray2DT& XDOF(const DOFElementT* group, int tag_set) const;

protected:

	/** call groups to reset external DOF's */
	void Reset(void);

	/** update DOF's using the global update vector */
	void Update(const dArrayT& update);

	/** (self-)configure element group */
	void ConfigureElementGroup(int group_number, int& tag_num);

	/** assign equation numbers */
	void SetEquations(int& num_eq);

	/** remove external DOF's from first slot of each row */
	void CheckEquationNumbers(ostream& out, iArray2DT& eqnos);

	/** resolve index of the tag set.
	 * \param group pointer to the DOFElementT
	 * \param tag_set set number for the element */
	int TagSetIndex(const DOFElementT* group, int tag_set) const;

protected:

	/** registered element groups */
	AutoArrayT<DOFElementT*> fDOFElements;
	
	/** number of tag sets by group */
	AutoArrayT<int> fNumTagSets;

	/** global equations numbers for each tag set */
	AutoArrayT<iArray2DT*> fXDOF_Eqnos;

	/** values of the degrees of freedom by set */
	AutoArrayT<dArray2DT*> fXDOFs;
};

#endif /* _XDOF_MANAGER_T_H_ */
