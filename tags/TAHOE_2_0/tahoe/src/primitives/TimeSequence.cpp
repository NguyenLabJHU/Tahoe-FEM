/* $Id: TimeSequence.cpp,v 1.6 2004-07-15 08:31:36 paklein Exp $ */
/* created: paklein (05/22/1996) */
#include "TimeSequence.h"

#include "ifstreamT.h"

using namespace Tahoe;

/* constructor */
TimeSequence::TimeSequence(void) { }

/* I/O operators */
void TimeSequence::Read(ifstreamT& in)
{
	in >> fNumSteps;  if (fNumSteps     < 0) throw ExceptionT::kBadInputValue;
	in >> fOutputInc;
	in >> fMaxCuts;	  if (fMaxCuts  <  0  ) throw ExceptionT::kBadInputValue;
	in >> fTimeStep;  if (fTimeStep <= 0.0) throw ExceptionT::kBadInputValue;
}

void TimeSequence::Write(ostream& out) const
{
	out << " Number of time steps. . . . . . . . . . . . . . = " << fNumSteps        << '\n';
	out << " Output print increment (< 0: current step size) = " << fOutputInc       << '\n';
	out << " Maximum number of load step cuts. . . . . . . . = " << fMaxCuts         << '\n';
	out << " Time step . . . . . . . . . . . . . . . . . . . = " << fTimeStep        << '\n';
}
