/* $Id: ScheduleT.h,v 1.4 2002-07-05 22:28:33 paklein Exp $ */
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

	/* interpolate load factor to time tim */
	void SetTime(double time);

	/* return current value of the load factor */
	double Value(void) const;
	double Value(double time);

private:

	/* check that times are sequential */
	void CheckSequential(void) const;

private:

	/** \name function data */
	/*@{*/
	dArrayT fTime;
	dArrayT fValue;
	/*@}*/
	
	/* current value */
	double fCurrentValue;
};

/* inlines */
inline double ScheduleT::Value(void) const { return fCurrentValue; }
inline double ScheduleT::Value(double time)
{
	SetTime(time);
	return fCurrentValue;
}

} // namespace Tahoe 
#endif /* _SCHEDULE_T_H_ */
