/* $Id: TimeSequence.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (05/22/1996)                                          */

#include "TimeSequence.h"

#include "fstreamT.h"

#include "Constants.h"
#include "ExceptionCodes.h"

/* constructor */
TimeSequence::TimeSequence(void) { }

/* I/O operators */
void TimeSequence::Read(ifstreamT& in)
{
	in >> fNumSteps;  if (fNumSteps     < 0) throw eBadInputValue;
	in >> fOutputInc;
	in >> fMaxCuts;	  if (fMaxCuts  <  0  ) throw eBadInputValue;
	in >> fTimeStep;  if (fTimeStep <= 0.0) throw eBadInputValue;
}

void TimeSequence::Write(ostream& out) const
{
	out << " Number of time steps. . . . . . . . . . . . . . = " << fNumSteps        << '\n';
	out << " Output print increment (< 0: current step size) = " << fOutputInc       << '\n';
	out << " Maximum number of load step cuts. . . . . . . . = " << fMaxCuts         << '\n';
	out << " Time step . . . . . . . . . . . . . . . . . . . = " << fTimeStep        << '\n';
}
