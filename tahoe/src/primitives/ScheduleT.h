/* $Id: ScheduleT.h,v 1.6 2004-07-15 08:31:36 paklein Exp $ */
/* created: paklein (05/24/1996) */
#ifndef _SCHEDULE_T_H_
#define _SCHEDULE_T_H_

/* direct members */
#include "ParameterInterfaceT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class C1FunctionT;

/** schedule functions comprised of a C1FunctionT and a cached time and function value */
class ScheduleT: public ParameterInterfaceT
{
public:

	/** constructors */
	/*@{*/
	ScheduleT(void); 

	/** construct schedule with constant value */
	ScheduleT(double value); 
	/*@}*/

	/** destructor */
	~ScheduleT(void); 

	/** set schedule to the given time */
	void SetTime(double time);

	/** \name current values */
	/*@{*/
	/** schedule value at the given time */
	double Value(void) const { return fCurrentValue; };
	
	/** get value at the given time. Call does not change the internal time,
	 * which must be set with ScheduleT::SetTime */
	double Value(double time) const;
	
	/** the internal time */
	double Time(void) const { return fCurrentTime; };
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/* information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/** function data */
	C1FunctionT* fFunction;
	
	/** \name current time and value */
	/*@{*/
	double fCurrentTime;
	double fCurrentValue;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _SCHEDULE_T_H_ */
