/* $Id: TimeManagerT.cpp,v 1.26 2009-04-23 15:01:09 tdnguye Exp $ */
/* created: paklein (05/23/1996) */
#include "TimeManagerT.h"

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "toolboxConstants.h"
#include "FEManagerT.h"
#include "dArrayT.h"
#include "NodeManagerT.h"

using namespace Tahoe;

/* constructor */
TimeManagerT::TimeManagerT(FEManagerT& FEM):
	ParameterInterfaceT("time"),
	fFEManager(FEM),

	/* runtime data for the current sequence */
	fStepNum(0),
	fTime(0.0),
	fNumStepCuts(0),
	fStepCutStatus(kSameStep),
	
	fNumSteps(0), fOutputInc(1), fMaxCuts(0), fTimeStep(1.0),
	fIsTimeShifted(0), fTimeShift(0.0),
	fImpExp(IntegratorT::kImplicit)
{

}

/* set to initial conditions */
void TimeManagerT::InitialCondition(void)
{
	/* broadcast time step */
	fFEManager.SetTimeStep(fTimeStep);

	/* nodes */
	const NodeManagerT* nodes = fFEManager.NodeManager();
	if (!nodes) ExceptionT::GeneralFail("TimeManagerT::InitialCondition");

	/* see if all integrators are explicit */
	fImpExp = IntegratorT::kExplicit;
	int num_groups = fFEManager.NumGroups();
	for (int i = 0; i < num_groups && fImpExp == IntegratorT::kExplicit; i++)
		fImpExp = nodes->ImplicitExplicit(i);
}

/* advance one time step, return 0 when if the sequence
* has ended */
bool TimeManagerT::Step(void)
{
	/* check that time has not been shifted */
	if (fIsTimeShifted) throw ExceptionT::kGeneralFail;

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
		bool is_explicit = fImpExp == IntegratorT::kExplicit;

		/* verbose flag */
		bool write_header = !is_explicit     ||
		                    fOutputInc == -1 ||
		                    fOutputInc == 0  ||
		                    fStepNum % fOutputInc == 0;

		
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
	const char caller[] = "TimeManagerT::ResetStep";

	/* check that time has not been shifted */
	if (fIsTimeShifted) ExceptionT::GeneralFail(caller);

	/* too far */
	if (fStepNum == 0)
		 ExceptionT::GeneralFail(caller, "already at the start time");
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

/* return a pointer to the specified ScheduleT function */
ScheduleT* TimeManagerT::Schedule(int num) const
{
	/* range check */
	if (num < 0 || num >= fSchedule.Length())
		ExceptionT::OutOfRange("TimeManagerT::Schedule", "schedule %d is out of range {0,%d}",
			num, fSchedule.Length() - 1);

	return fSchedule[num];
}

double TimeManagerT::ScheduleValue(int num) const
{
	return fSchedule[num]->Value();
}

/* initialize/restart functions
*
* Initialize functions reset all kinematic data to the
* default initial system state.  The restart functions
* should read/write any data that overrides the default
* values */
void TimeManagerT::ReadRestart(istream& restart_in)
{	
#if 0 
	/* total desired simulation time */ 
	double tot_time = fTimeStep*fNumSteps; 
         
	restart_in >> fTime >> fNumStepCuts >> fTimeStep; 
	cout << " Restarting from time: " << fTime << '\n';

	/* reset number of steps to preserve the total time */ 
	double float_steps = tot_time/fTimeStep; 
	fNumSteps = int((2.0*float_steps + 1.0)/2.0); 
#endif 
 
	//TEMP - do no use the time step written in the restart file
	//       will become restart option
	double time_step; 
	restart_in >> fTime >> fNumStepCuts >> time_step; 
	cout << " Restarting from time: " << fTime << '\n';  
	cout << "   Previous step cuts: " << fNumStepCuts << '\n';  
	cout << "   Previous step size: " << time_step << '\n';  
	cout << "    Current step size: " << fTimeStep << '\n';
        
	/* set load factors */
	for (int i = 0; i < fSchedule.Length(); i++)
		fSchedule[i]->SetTime(fTime);
}

void TimeManagerT::WriteRestart(ostream& restart_out) const
{
	restart_out << fTime            << '\n';
	restart_out << fNumStepCuts     << '\n';
	restart_out << fTimeStep        << '\n';
}

/* return true if output should be written for the current step */
bool TimeManagerT::WriteOutput(void) const
{
	if (fabs(fTimeStep) > kSmall && /* no output if clock is not running */
		fOutputInc != 0 &&
	   		(fStepNum % fOutputInc == 0 || /* at increment */
	    	 fStepNum == fNumSteps) /* at end */
	)
		return true;
	else
		return false;
}

/* describe the parameters needed by the interface */
void TimeManagerT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	ParameterT scaling(ParameterT::Enumeration, "scaling");
	scaling.AddEnumeration("linear", kLinear);
	scaling.AddEnumeration("logarithmic", kLog);
	scaling.SetDefault(kLinear);
	list.AddParameter(scaling);
	
	ParameterT init_time(ParameterT::Double, "init_time");
	init_time.AddLimit(0, LimitT::LowerInclusive);
	init_time.SetDefault(0.0);
	list.AddParameter(init_time);
	
	/* number of steps */
	ParameterT num_steps(ParameterT::Integer, "num_steps");
	num_steps.AddLimit(0, LimitT::LowerInclusive);
	list.AddParameter(num_steps);

	/* results output increment */
	ParameterT output_inc(ParameterT::Integer, "output_inc");
	output_inc.SetDefault(0);
	list.AddParameter(output_inc);

	/* maximum number of time increment cuts */
	ParameterT max_step_cuts(ParameterT::Integer, "max_step_cuts");
	max_step_cuts.SetDefault(0);
	max_step_cuts.AddLimit(0, LimitT::LowerInclusive);
	list.AddParameter(max_step_cuts);

	/* (initial) time increment */
	ParameterT time_step(ParameterT::Double, "time_step");
	time_step.AddLimit(0, LimitT::LowerInclusive);
	list.AddParameter(time_step);
}

/* accept parameter list */
void TimeManagerT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	fTimeScaling = list.GetParameter("scaling");
	fInitTime = list.GetParameter("init_time");
	
	/* set time parameters */
	fNumSteps  = list.GetParameter("num_steps");
	fOutputInc = list.GetParameter("output_inc");
	fMaxCuts   = list.GetParameter("max_step_cuts");
	fTimeStep  = list.GetParameter("time_step");	

	/* runtime data for the current sequence */
	fStepNum       = 0;
	fTime          = 0.0;
	fNumStepCuts   = 0;
	fStepCutStatus = kSameStep;

	/* console variables */
	iSetName("time");
	iAddVariable("num_steps", fNumSteps);
	iAddVariable("output_inc", fOutputInc);
	iAddVariable("max_step_cuts", fMaxCuts);
	iAddVariable("time_step", fTimeStep);

	/* construct schedule functions */
	fSchedule.Dimension(list.NumLists("schedule_function"));
	for (int i = 0; i < fSchedule.Length(); i++)
	{
		fSchedule[i] = new ScheduleT;
		fSchedule[i]->TakeParameterList(*(list.List("schedule_function", i)));
		fSchedule[i]->SetTime(0.0);
	}
}

/* information about subordinate parameter lists */
void TimeManagerT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* schedule function */
	sub_list.AddSub("schedule_function", ParameterListT::Any);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* TimeManagerT::NewSub(const StringT& name) const
{
	if (name == "schedule_function")
		return new ScheduleT;
	else /* inherited */
		return ParameterInterfaceT::NewSub(name);
}

/************************************************************************
 * Private
 ************************************************************************/

/* increment the time and reset the load factors */
void TimeManagerT::IncrementTime(double dt)
{
	/* Increment the time */
	if (fTimeScaling == kLinear)
		fTime += dt;
	else 
	{
		if (StepNumber() == 1)
			fTime = fInitTime;
		else
			fTime *= pow(10, dt);
	}
		
	/* set load factors */
	for (int i = 0; i < fSchedule.Length(); i++)
		fSchedule[i]->SetTime(fTime);
}

/* returns 1 if the number is even, otherwise returns 0	*/
int TimeManagerT::IsEven(int number) const
{
	return int(number % 2 == 0);
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
	
	/* broadcast time step change */
	fFEManager.SetTimeStep(fTimeStep);
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

	/* broadcast time step change */
	fFEManager.SetTimeStep(fTimeStep);
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
