/* $Id: TorsionKBCT.h,v 1.4 2004-07-15 08:31:21 paklein Exp $ */
#ifndef _TORSION_KBC_T_H_
#define _TORSION_KBC_T_H_

/* base class */
#include "KBC_ControllerT.h"

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"
#include "ScheduleT.h"

namespace Tahoe {

/** torsion boundary condition controller. */
class TorsionKBCT: public KBC_ControllerT
{
public:

	/** constructor */
	TorsionKBCT(const BasicSupportT& support);

	/** set initial conditions */
	void InitialCondition(void);

	/** compute updated prescibed displacements */
	virtual void InitStep(void);

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

protected:

	/** time increment */
	double fStartTime;

	/** rotation rate (rad/s) */
	double fw;
	
	/** \name axis of rotation */
	/*@{*/
	/** direction */
	int fAxis;
	
	/** point on the axis */
	dArrayT fPoint;
	/*@}*/

	/** \name nodes */
	/*@{*/
	ArrayT<StringT> fID_List;
	iArrayT fNodes;
	/*@}*/

	/** runtime data */
	ScheduleT fDummySchedule;
};

} /* namespace Tahoe */

#endif /* _TORSION_KBC_T_H_ */
