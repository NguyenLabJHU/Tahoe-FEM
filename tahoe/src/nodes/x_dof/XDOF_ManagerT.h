/* $Id: XDOF_ManagerT.h,v 1.1.1.1 2001-01-29 08:20:40 paklein Exp $ */
/* created: paklein (06/01/1998)                                          */
/* base class which defines the interface for a manager                   */
/* of DOF's comprised of FE DOF's plus extenal (element)                  */
/* DOF's                                                                  */

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

class XDOF_ManagerT
{
public:
	
	/* constructor */
	XDOF_ManagerT(void);	

	/* destructor */
	virtual ~XDOF_ManagerT(void);

	/* add element group to list */
	virtual void Register(DOFElementT* group, int numDOF);
	
	/* get equation numbers for the specified group */
	virtual const iArray2DT& XDOF_Eqnos(const DOFElementT* group) const;

	/* returns reference to current constrain values */
	virtual const dArray2DT& XDOF(const DOFElementT* group) const;

protected:

	/* call groups to reset external DOF's */
	void Reset(void);

	/* update DOF's using the global update vector */
	void Update(const dArrayT& update);

	/* (self-)configure element group */
	void ConfigureElementGroup(int group_number, int& tag_num);

	/* assign equation numbers */
	void SetEquations(int& num_eq);

	/* remove external DOF's from first slot of each row */
	void CheckEquationNumbers(ostream& out, iArray2DT& eqnos);

protected:

	/* augmented Lagrangian element groups */
	AutoArrayT<DOFElementT*> fDOFElements;

	/* DOF and global equations */
	AutoArrayT<iArray2DT*> fXDOF_Eqnos;
	AutoArrayT<dArray2DT*> fXDOFs;
};

#endif /* _XDOF_MANAGER_T_H_ */
