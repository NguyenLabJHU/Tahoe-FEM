/* $Id: TimeSequence.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (05/22/1996)                                          */

#ifndef _TIMESEQ_H_
#define _TIMESEQ_H_

#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"
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

#endif /* _TIMESEQ_H_ */
