/* $Id: TimeSequence.h,v 1.1.1.1.10.1 2002-06-27 18:04:02 cjkimme Exp $ */
/* created: paklein (05/22/1996)                                          */

#ifndef _TIMESEQ_H_
#define _TIMESEQ_H_

#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"

namespace Tahoe {

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
