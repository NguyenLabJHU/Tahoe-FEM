/* $Id: ThermalDilatationT.h,v 1.5 2002-07-05 22:28:27 paklein Exp $ */
/* created: paklein (08/25/1996) */

#ifndef _THERMALDILAT_H_
#define _THERMALDILAT_H_

#include "Environment.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class ScheduleT;

class ThermalDilatationT
{
public:

	/* constructor */
	ThermalDilatationT(ifstreamT& in);
	
	/* to set LTf pointer */
	int ScheduleNum(void) const;
	void SetSchedule(const ScheduleT* LTf);

	/* returns true if active */
	bool IsActive(void) const;
	
	/* I/O functions */
	void Print(ostream& out) const;

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

} // namespace Tahoe 
#endif /* _THERMALDILAT_H_ */
