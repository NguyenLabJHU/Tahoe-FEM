/* $Id: TimeSequence.h,v 1.3 2002-07-05 22:28:33 paklein Exp $ */
/* created: paklein (05/22/1996)                                          */

#ifndef _TIMESEQ_H_
#define _TIMESEQ_H_

#include "Environment.h"

#include "ios_fwd_decl.h"


namespace Tahoe {

/* forward declarations */
class ifstreamT;

class TimeSequence
{
friend class TimeManagerT;

public:

	/* constructor */
	TimeSequence(void);
	
	/* I/O */
	void Read(ifstreamT& in);
	void Write(ostream& out) const;

private:
	
	int	   fNumSteps;
	int	   fOutputInc;
	int	   fMaxCuts;
	double fTimeStep;
};

} // namespace Tahoe 
#endif /* _TIMESEQ_H_ */
