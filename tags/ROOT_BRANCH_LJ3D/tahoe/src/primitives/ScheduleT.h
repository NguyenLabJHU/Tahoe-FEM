/* $Id: ScheduleT.h,v 1.5 2003-10-28 07:12:14 paklein Exp $ */
/* created: paklein (05/24/1996) */

#ifndef _SCHEDULE_T_H_
#define _SCHEDULE_T_H_

/* direct members */
#include "dArrayT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** the class formerly known as LoadTime. Piecewise linear function */
class ScheduleT
{
public:

	/** \name constructors */
	/*@{*/
	/** constant value for all time */
	ScheduleT(double value); 
	ScheduleT(int numpts);
	ScheduleT(const dArrayT& times, const dArrayT& values);
	/*@}*/

	/* I/O */
	void Read(ifstreamT& in);
	void Write(ostream& out) const;

	/** set schedule to the given time */
	void SetTime(double time);

	/** \name current values */
	/*@{*/
	/** schedule value at the given time */
	double Value(void) const;
	
	/** get value at the given time. Call does not change the internal time,
	 * which must be set with ScheduleT::SetTime */
	double Value(double time) const;
	
	/** the internal time */
	double Time(void) const { return fCurrentTime; };
	/*@}*/

private:

	/* check that times are sequential */
	void CheckSequential(void) const;

private:

	/** \name function data */
	/*@{*/
	dArrayT fTime;
	dArrayT fValue;
	/*@}*/
	
	/** \name current time and value */
	/*@{*/
	double fCurrentTime;
	double fCurrentValue;
	/*@}*/
};

/* inlines */
inline double ScheduleT::Value(void) const { return fCurrentValue; }
inline double ScheduleT::Value(double time) const
{
	/* non-const temporary */
	ScheduleT* non_const_this = (ScheduleT*) this;
	double curr_time = non_const_this->Time();
	non_const_this->SetTime(time);
	double curr_value = non_const_this->Value();
	non_const_this->SetTime(curr_time);
	
	return curr_value;
}

} // namespace Tahoe 
#endif /* _SCHEDULE_T_H_ */
