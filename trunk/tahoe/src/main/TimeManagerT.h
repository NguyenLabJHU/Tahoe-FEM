/* $Id: TimeManagerT.h,v 1.11 2003-10-28 07:36:49 paklein Exp $ */
/* created: paklein (05/23/1996) */

#ifndef _TIMEMANAGER_T_H_
#define _TIMEMANAGER_T_H_

/* base class */
#include "iConsoleObjectT.h"
#include "ParameterInterfaceT.h"

/* direct members */
#include "StringT.h"
#include "pArrayT.h"
#include "ScheduleT.h"
#include "iAutoArrayT.h"
#include "IOBaseT.h"
#include "TimeSequence.h"
#include "IntegratorT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class FEManagerT;
class CoordinatorT;
class nLinearHHTalpha;
class NodeManagerT;

class TimeManagerT: public iConsoleObjectT, public ParameterInterfaceT
{
	/* time shifters */
	friend class LinearHHTalpha;
	friend class NLHHTalpha;

public:

	/** enum of integrator types */
	enum CodeT {
		kLinearStatic = 0,
		      kStatic = 1,
           kTrapezoid = 2,
           kLinearHHT = 3,
		kNonlinearHHT = 4,
		  kExplicitCD = 5,
		kVerlet = 6,
		kGear6 = 7
	};

	/** stream extraction operator */
	friend istream& operator>>(istream& in, TimeManagerT::CodeT& code);

	/** constructor */
	TimeManagerT(FEManagerT& FEM);

	/** initialization */
	void Initialize(void);

	/* run through the time sequences */
	void Top(void);
	bool NextSequence(void);

	/* time sequence */
	bool Step(void);
	void ResetStep(void);
	const int& StepNumber(void) const;
	const int& NumberOfSteps(void) const;
	const double& Time(void) const;
	const double& TimeStep(void) const;

	/** set the time step. Use at your own risk. TimeManagerT manages its
	 * own time step. This method can be used to override this time step
	 * control. However, the TimeManagerT's state will be inconsistent until
	 * the step is restored. */
	void SetTimeStep(double dt) { fTimeStep = dt; };

	/* load control functions (returns true if successful) */
	bool DecreaseLoadStep(void);
	bool IncreaseLoadStep(void);
	
	/* return a pointer to the specified ScheduleT function */

	/** \name schedule information */
	/*@{*/
	int NumSchedule(void) const;
	ScheduleT* Schedule(int num) const;
	double ScheduleValue(int num) const;
	/*@}*/

	/* accessors */
	int SequenceNumber(void) const;
	int NumSequences(void) const;
			
	/* initialize/restart functions
	 *
	 * Initialize functions reset all kinematic data to the
	 * default initial system state.  The restart functions
	 * should read/write any data that overrides the default
	 * values */
	void ReadRestart(istream& in);
	void WriteRestart(ostream& out) const;

	/** return true if output should be written for the current step */
	bool WriteOutput(void) const;

	/** return a pointer to a integrator of the specified type */
	IntegratorT* New_Integrator(CodeT type) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:	

	/** step cut status flags */
	enum StatusT { kDecreaseStep =-1,
                       kSameStep = 0,
                   kIncreaseStep = 1};
	
	/* output functions */
	void EchoTimeSequences(ifstreamT& in, ostream& out);
	void EchoSchedule(ifstreamT& in, ostream& out);

	/* increment the time and reset the load factors */
	void IncrementTime(double dt);

	/* returns 1 if the number is even, otherwise returns 0 */
	int IsEven(int number) const;
	
	/* adjust time stepping parameters */
	void DoDecreaseStep(void);
	void DoIncreaseStep(void);

private:

	FEManagerT& theBoss;
	
	ArrayT<TimeSequence> fSequences;
	pArrayT<ScheduleT*>  fSchedule;

	/** \name copied from current sequence */
	/*@{*/
	int	   fNumSteps;
	int	   fOutputInc;
	int	   fMaxCuts;
	double fTimeStep;
	/*@}*/
	
	/** \name runtime data for the current sequence */
	/*@{*/
	int	   fCurrentSequence;
	int	   fStepNum;
	double fTime;
	int    fNumStepCuts;
	int	   fStepCutStatus;
	/*@}*/

	/* time stepper */
	int	   fIsTimeShifted;
	double fTimeShift;
	
	/** will be IntegratorT::kExplicit if all integrators are explicit
	 * otherwise will be IntegratorT::kImplicit */
	IntegratorT::ImpExpFlagT fImpExp;

/* functions for time shifters */
private:

	/* to allow LinearHHTalpha to adjust the time.  LinearHHTalpha must
	 * call ResetTime when finished.  MUST call ResetTime before the next call
	 * to Step */
	void ShiftTime(double dt);
	
	/* reset the time back to what it was before the calls to IncrementTime */
	void ResetTime(void);
};

/* inlines */
inline const double& TimeManagerT::Time(void) const { return fTime; }
inline const double& TimeManagerT::TimeStep(void) const { return fTimeStep; }
inline const int& TimeManagerT::StepNumber(void) const { return fStepNum; }
inline const int& TimeManagerT::NumberOfSteps(void) const { return fNumSteps; }

inline int TimeManagerT::NumSchedule(void) const { return fSchedule.Length() ; }
inline int TimeManagerT::SequenceNumber(void) const { return fCurrentSequence; }
inline int TimeManagerT::NumSequences(void) const { return fSequences.Length(); }
} // namespace Tahoe 
#endif /* _TIMEMANAGER_T_H_ */
