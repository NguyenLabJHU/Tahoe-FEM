/* $Id: NodeManagerT.h,v 1.1.1.1 2001-01-29 08:20:39 paklein Exp $ */
/* created: paklein (05/23/1996)                                          */
/* Field variables plus averging                                          */

#ifndef _NODEMANAGER_T_H_
#define _NODEMANAGER_T_H_

/* base classes */
#include "GroupAverageT.h"
#include "NodeManagerPrimitive.h"

/* forward declarations */
template <class TYPE> class AutoArrayT;
class dArray2DT;
class StringT;

class NodeManagerT: public NodeManagerPrimitive, public GroupAverageT
{
public:
	
	/* constructor */
	NodeManagerT(FEManagerT& fe_manager);	

protected:

	/* initialization */
	virtual void AllocateGlobal(void);

	/* BC Controllers */
	virtual KBC_ControllerT* NewKBC_Controller(int code);

private:

	/* return the number of DOF's */
	virtual int DegreesOfFreedom(int nsd) const;
	
};

#endif /* _NODEMANAGER_T_H_ */
