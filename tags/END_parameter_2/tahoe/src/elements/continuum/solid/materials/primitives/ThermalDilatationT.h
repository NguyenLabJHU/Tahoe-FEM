/* $Id: ThermalDilatationT.h,v 1.5.40.1 2004-02-18 16:33:51 paklein Exp $ */
/* created: paklein (08/25/1996) */
#ifndef _THERMALDILAT_H_
#define _THERMALDILAT_H_

#include "Environment.h"

namespace Tahoe {

/* forward declarations */
class ScheduleT;

class ThermalDilatationT
{
public:

	/* constructor */
	ThermalDilatationT(void);
	
	/* to set LTf pointer */
	int ScheduleNum(void) const;
	void SetSchedule(const ScheduleT* LTf);

	/** set the schedule number */
	void SetScheduleNum(int schedule) { LTfnum = schedule; };

	/** set the dilation */
	void SetPercentElongation(double perccent_elongation) { fPercentElongation = perccent_elongation; };

	/* returns true if active */
	bool IsActive(void) const;

	/* returns the current elongation factor */
	double PercentElongation(void) const;
							
private:
	
	double fPercentElongation;
	int LTfnum;
	const ScheduleT* LTfPtr;	
};

/* inline functions */

/* returns true if active */
inline bool ThermalDilatationT::IsActive(void) const { return fPercentElongation != 0.0; }

/* set LTf pointer */
inline int ThermalDilatationT::ScheduleNum(void) const { return LTfnum; }
inline void ThermalDilatationT::SetSchedule(const ScheduleT* LTf)
{ 
	LTfPtr = LTf; 
	if (!LTfPtr) 
		fPercentElongation = 0.0;
}

} /* namespace Tahoe */

#endif /* _THERMALDILAT_H_ */
