/* $Id: TimeSequence.cpp,v 1.3 2002-09-12 17:50:07 paklein Exp $ */
/* created: paklein (05/22/1996)                                          */

#include "TimeSequence.h"

#include "fstreamT.h"

#include "toolboxConstants.h"
#include "ExceptionCodes.h"

/* constructor */

using namespace Tahoe;

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
