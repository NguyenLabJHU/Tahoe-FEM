/* $Id: TimeManagerT.cpp,v 1.1.1.1 2001-01-29 08:20:21 paklein Exp $ */
/* created: paklein (05/23/1996)                                          */

#include "TimeManagerT.h"

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "fstreamT.h"

#include "Constants.h"
#include "FEManagerT.h"
#include "TimeSequence.h"
#include "dArrayT.h"
#include "NodeManagerT.h"

/* step cut status flags */
const int kDecreaseStep =-1;
const int kSameStep     = 0;
const int kIncreaseStep = 1;

/* controllers */
#include "HHTalpha.h"

/* constructor */
TimeManagerT::TimeManagerT(FEManagerT& FEM):
	theBoss(FEM),

	/* runtime data for the current sequence */
	fCurrentSequence(-1),
	fStepNum(0),
	fTime(0.0),
	fNumStepCuts(0),
	fStepCutStatus(kSameStep),
	fOutputCount(0),
	
	fNumSteps(0), fOutputInc(1), fMaxCuts(0), fTimeStep(1.0),
	fIsTimeShifted(0), fTimeShift(0.0)
{
	ifstreamT& in  = theBoss.Input();
	ostream&   out = theBoss.Output();

	/* Time sequences - allocate memory and echo */
	int num_sequences;
	in >> num_sequences;	
	if (num_sequences < 1) throw eBadInputValue;
	fSequences.Allocate(num_sequences);
	EchoTimeSequences(in, out);
	
	/* Loadtime functions - allocate memory and echo */	
	int num_LTf;
	in >> num_LTf;
	if (num_LTf < 1) throw eBadInputValue;
	fLTf.Allocate(num_LTf); // add: f(t) = 1.0

	EchoLoadTime(in, out);
	
	/* console variables */
	iSetName("time");
	iAddVariable("num_steps", fNumSteps);
	iAddVariable("output_inc", fOutputInc);
	iAddVariable("max_step_cuts", fMaxCuts);
	iAddVariable("time_step", fTimeStep);
}

/* run through the time sequences.  NextSequence returns 0
* if there are no more time sequences */
void TimeManagerT::Top(void)
{
	fCurrentSequence = -1;
}

bool TimeManagerT::NextSequence(void)
{
	fCurrentSequence++;
	
	/* initialize next sequence */
	if (fCurrentSequence < fSequences.Length())
	{
		TimeSequence& curr_sequence = fSequences[fCurrentSequence];
		
		/* runtime data for the current sequence */
		fStepNum       = 0;
		fTime          = 0.0;
		fNumStepCuts   = 0;
		fStepCutStatus = kSameStep;
		fOutputCount   = 0;
		
		/* copy from current sequence */
		fNumSteps  = curr_sequence.fNumSteps;
		fOutputInc = curr_sequence.fOutputInc;
		fMaxCuts   = curr_sequence.fMaxCuts;
		fTimeStep  = curr_sequence.fTimeStep;
		
		/* broadcast time step */
		theBoss.SetTimeStep(fTimeStep);

		/* print headers */
		theBoss.Output() << "\n T i m e   S e q u e n c e : ";
		theBoss.Output() << fCurrentSequence + 1 << "\n\n";
		cout << "\n T i m e   S e q u e n c e : ";
		cout << fCurrentSequence + 1 << endl;
		return true;
	}
	else
		return false;
}

/* advance one time step, return 0 when if the sequence
* has ended */
bool TimeManagerT::Step(void)
{
	/* check that time has not been shifted */
	if (fIsTimeShifted) throw eGeneralFail;

	if (fStepNum < fNumSteps)
	{
		/* adjust time step */
		if (fStepCutStatus == kDecreaseStep)
			DoDecreaseStep();
		else if (fStepCutStatus == kIncreaseStep)
			DoIncreaseStep();
			
		/* advance time */
		fStepNum++;
		IncrementTime(fTimeStep);
		
		/* print less often for explicit */
		GlobalT::AnalysisCodeT analysiscode = theBoss.Analysis();
		bool is_explicit = (analysiscode == GlobalT::kLinExpDynamic ||
		                    analysiscode == GlobalT::kNLExpDynamic    ||
		                    analysiscode == GlobalT::kVarNodeNLExpDyn ||
		                    analysiscode == GlobalT::kNLExpDynKfield);
		
		/* verbose flag */
		bool write_header = !is_explicit     ||
		                    fOutputInc == -1 ||
		                    fOutputInc == 0  ||
		                    fmod(fStepNum, fOutputInc) == 0;

		
		/* step header */
		int start_count = 5;
		if (write_header || (is_explicit && fStepNum <= start_count))
	{
			cout << "\n     Time: " << fTime << '\n';
			cout <<   "Step size: " << fTimeStep << '\n';
			cout <<   "     Step: " << fStepNum << " of " << fNumSteps << '\n';
			cout << endl;
		}
		if (is_explicit && fStepNum == start_count && fOutputInc > 1 && fStepNum < fNumSteps)
			cout << " (marker at output increments only)" << endl;
				
		return true;
	}
	else
		/* reached end of time sequence */
		return false;
}

void TimeManagerT::ResetStep(void)
{
	/* check that time has not been shifted */
	if (fIsTimeShifted) throw eGeneralFail;

	/* too far */
	if (fStepNum == 0)
	{
		cout << "\n TimeManagerT::ResetStep: already at the start time"<< endl;
		throw eGeneralFail;
	}
	/* return to previous time */
	else
	{		
		IncrementTime(-fTimeStep);
		fStepNum--;
	}
}

/* load control functions */
bool TimeManagerT::DecreaseLoadStep(void)
{
	if (fNumStepCuts >= fMaxCuts)
	{
		cout << "\n TimeManagerT::DecreaseLoadStep: load increment cut ";
		cout << fMaxCuts << " times.\n";
		cout << "\n Exiting time sequence\n" << endl;

		/* end current time sequence */
		fStepNum = fNumSteps;		
		return false;
	}
	else
	{
		/* set flag */
		fStepCutStatus = kDecreaseStep;		
		return true;
	}
}

bool TimeManagerT::IncreaseLoadStep(void)
{
	/* set flag */
	if (fNumStepCuts > 0  &&		/* only as big as original */
		IsEven(fStepNum)  && 		/* step number divisible by 2 */		
		IsEven(fNumSteps) && 		/* total steps divisible by 2 */
	    fStepNum < fNumSteps )		/* not the last step */
	{
		fStepCutStatus = kIncreaseStep;
		return true;
	}
	else
		return false;
}

/* return a pointer to the specified LoadTime function */
LoadTime* TimeManagerT::GetLTf(int num) const
{
	/* range check */
	if (num < 0 || num >= fLTf.Length())
	{
		cout << "\n TimeManagerT::GetLTf: function number " << num << " is out of\n"
		     <<   "     range {" << 0 << "," << fLTf.Length() - 1 << "}" << endl;
		throw eOutOfRange;
	}

	return fLTf[num];
}

double TimeManagerT::LoadFactor(int nLTf) const
{
	return fLTf[nLTf]->LoadFactor();
}

/* initialize/restart functions
*
* Initialize functions reset all kinematic data to the
* default initial system state.  The restart functions
* should read/write any data that overrides the default
* values */
void TimeManagerT::ReadRestart(istream& restart_in)
{	
	int sequencenumber;
	restart_in >> sequencenumber;
	if (sequencenumber != fCurrentSequence) throw eBadInputValue;

	/* total desired simulation time */
	double tot_time = fTimeStep*fNumSteps;
	
	restart_in >> fTime >> fNumStepCuts >> fTimeStep;
	cout << " Restarting from time: " << fTime << '\n';

	/* reset number of steps to preserve the total time */
	double float_steps = tot_time/fTimeStep;
	fNumSteps = int((2.0*float_steps + 1.0)/2.0);
	
	/* set load factors */
	for (int i = 0; i < fLTf.Length(); i++)
		fLTf[i]->SetLoadFactor(fTime);
}

void TimeManagerT::WriteRestart(ostream& restart_out) const
{
	restart_out << fCurrentSequence << '\n';
	restart_out << fTime            << '\n';
	restart_out << fNumStepCuts     << '\n';
	restart_out << fTimeStep        << '\n';
}

/* finalize step (trigger output) */
void TimeManagerT::CloseStep(void) //TEMP? - let FEManager control/monitor output?
{
	if (fOutputInc != 0 &&
	   (++fOutputCount >= ((fOutputInc > 0) ? fOutputInc :-fOutputInc) || // at increment
	    fStepNum == fNumSteps)) // at end
	{
		/* run time output */
		theBoss.WriteOutput(IOBaseT::kAtInc);
			
		/* reset count */
		fOutputCount = 0;
	}
}

/************************************************************************
* Private
************************************************************************/

void TimeManagerT::EchoTimeSequences(ifstreamT& in, ostream& out)
{
	int num_seq = fSequences.Length();
	out << "\n T i m e   S e q u e n c e   D a t a :\n\n";
	out << " Number of time sequences  . . . . . . . . . . . = " << num_seq;
	out << "\n\n";
	
	for (int i = 0; i < num_seq; i++)
	{
		int seqnum;
		in >> seqnum;	
		if (seqnum < 1 ||
		    seqnum > num_seq) throw eBadInputValue;

		out << " Sequence number . . . . . . . . . . . . . . . . = ";
		out << seqnum << '\n';
		
		/* echo data */
		seqnum--;
		fSequences[seqnum].Read(in);
		fSequences[seqnum].Write(out);
		out << '\n';
	}
}

void TimeManagerT::EchoLoadTime(ifstreamT& in, ostream& out)
{
	int num_LTf = fLTf.Length();
	out << "\n L o a d - T i m e   F u n c t i o n   D a t a :\n\n";
	out << " Number of load-time functions . . . . . . . . . = " << num_LTf << '\n';

	for (int i = 0; i < num_LTf; i++)
	{
		int LTfnum, numpts;
		in >> LTfnum >> numpts;

		/* checks */
		if (LTfnum < 1 || LTfnum > num_LTf) throw eBadInputValue;
		if (numpts < 1) throw eBadInputValue;

		out << " Loadtime function number. . . . . . . . . . . . = ";
		out << LTfnum << "\n\n";
		
		/* echo data */
		LTfnum--;
		fLTf[LTfnum] = new LoadTime(numpts);
		if (!fLTf[LTfnum]) throw(eOutOfMemory);

		fLTf[LTfnum]->Read(in);
		fLTf[LTfnum]->Write(out);
		out << '\n';
	}
}

/* increment the time and reset the load factors */
void TimeManagerT::IncrementTime(double dt)
{
	/* Increment the time */
	fTime += dt;

	/* set load factors */
	for (int i = 0; i < fLTf.Length(); i++)
		fLTf[i]->SetLoadFactor(fTime);
}

/* returns 1 if the number is even, otherwise returns 0	*/
int TimeManagerT::IsEven(int number) const
{
	return (fmod(number,2) < 0.5) ? 1 : 0;
}

/* adjust time stepping parameters */
void TimeManagerT::DoDecreaseStep(void)
{
	fNumStepCuts++;
	fStepCutStatus = kSameStep;

	/* set new time step parameters */
	fTimeStep    /= 2;
	fOutputInc   *= (fOutputInc < 0) ? 1 : 2;
	fNumSteps    *= 2;
	fStepNum     *= 2;
	fOutputCount *= 2;
	
	/* broadcast time step change */
	theBoss.SetTimeStep(fTimeStep);
	cout << "\n TimeManagerT: step reduced to: " << fTimeStep << endl;
}

void TimeManagerT::DoIncreaseStep(void)
{
	fNumStepCuts--;
	fStepCutStatus = kSameStep;
	
	/* set new time step parameters */
	fTimeStep    *= 2;
	fOutputInc   /= (fOutputInc < 0) ? 1 : 2;
	fNumSteps    /= 2;	
	fStepNum     /= 2;
	fOutputCount /= 2;

	/* broadcast time step change */
	theBoss.SetTimeStep(fTimeStep);

	cout << "\n TimeManagerT: step increased to: " << fTimeStep << endl;
}

/* to allow LinearHHTalpha to adjust the time.  LinearHHTalpha must
* also restore then time when finished */
void TimeManagerT::ShiftTime(double dt)
{
	/* set flag */
	fIsTimeShifted = 1;
	fTimeShift   += dt;

	IncrementTime(dt);
}

/* reset the time back to what it was before the calls to IncrementTime */
void TimeManagerT::ResetTime(void)
{
	/* unset flag */
	fIsTimeShifted = 0;

	IncrementTime(-fTimeShift);
	fTimeShift = 0.0;
}
