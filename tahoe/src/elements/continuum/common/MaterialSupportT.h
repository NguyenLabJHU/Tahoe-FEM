/* $Id: MaterialSupportT.h,v 1.2.8.3 2002-10-30 09:18:11 paklein Exp $ */
#ifndef _MATERIAL_SUPPORT_T_H_
#define _MATERIAL_SUPPORT_T_H_

/* direct members */
#include "GlobalT.h"
#include "LocalArrayT.h"

namespace Tahoe {

/* forward declarations */
class ContinuumElementT;
class ElementCardT;
class ScheduleT;

/** support for the Tahoe materials classes. */
class MaterialSupportT
{
public:

	/** constructor */
	MaterialSupportT(int nsd, int ndof, int nip);

	/** destructor */
	virtual ~MaterialSupportT(void);

	/** \name dimensions */
	/*@{*/
	/** number of spatial dimensions */
	int NumSD(void) const { return fNumSD; };
	
	/** number of degrees of freedom (per node) */
	int NumDOF(void) const { return fNumDOF; };

	/** stress evaluation points per element */
	int NumIP(void) const { return fNumIP; };
	/*@}*/

	/** \name run time status */
	/*@{*/
	/** return a const reference to the run state flag */
	const GlobalT::StateT RunState(void) const;

	/** current stress evaluation point within the element. If
	 * no source for the current point is set using 
	 * MaterialSupportT::SetCurrIP, will return 0. */
	int CurrIP(void) const;

	/** the iteration number for the current time increment. If
	 * no source for the iteration number is set using 
	 * MaterialSupportT::SetIterationNumber, will return -1. */
	int IterationNumber(void) const;
	
	/** the current simulation time */
	double Time(void) const;

	/** the simulation time increment */
	double TimeStep(void) const;

	/** the simulation time increment number */
	int StepNumber(void) const;

	/** set the source for the run state flag */
	void SetRunState(const GlobalT::StateT& run_state);

	/** set the source for the current stress evaluation point */
	void SetCurrIP(const int& curr_ip);

	/** set the source for the iteration number */
	void SetIterationNumber(const int& iter);

	/** set source for the current simulation time */
	void SetTime(double& time);

	/** set source for the simulation time increment */
	void SetTimeStep(double& time_step);

	/** set source for the simulation time increment number */
	void SetStepNumber(int& step_number);
	/*@}*/
	
	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	const ContinuumElementT* ContinuumElement(void) const;

	/** set the element group pointer */
	virtual void SetContinuumElement(const ContinuumElementT* p);

	int NumElements(void) const;
	int CurrElementNumber(void) const;
	ElementCardT& ElementCard(int card) const;
	ElementCardT& CurrentElement(void) const;

	/** return a pointer the specified local array, or NULL if the array is not
	 * available. During calls the materials routines these will contain the
	 * values for the current element. */
	virtual const LocalArrayT* LocalArray(LocalArrayT::TypeT t) const;

	/** return a pointer to the specified LoadTime function */
	const ScheduleT* Schedule(int num) const;

	/** interpolate the given field to the current integration point. Returns true if the
	 * field is available, false otherwise. */
	bool Interpolate(const LocalArrayT& u, dArrayT& u_ip) const;

	/** interpolate the given field to the given integration point. Returns true if the
	 * field is available, false otherwise. */
	bool Interpolate(const LocalArrayT& u, dArrayT& u_ip, int ip) const;
	/*@}*/

  private:
  
  	/** \name dimensions */
  	/*@{*/
	/** number of degrees of spatial dimensions */
	int fNumSD;

	/** number of degrees of freedom */
	int fNumDOF;
	
	/** number of integration points */
	int fNumIP;
  	/*@}*/
  	
  	/** sources for run time information */
  	/*@{*/
	const GlobalT::StateT* fRunState;
	const int* fCurrIP;
	const int* fIterationNumber;
	const double* fTime;
	const double* fTimeStep;
	const int* fStepNumber;
  	/*@}*/
  
  	/** pointer to the continuum element */
  	const ContinuumElementT* fContinuumElement;
};

/* inlines functions */
inline const ContinuumElementT* MaterialSupportT::ContinuumElement(void) const
{
	return fContinuumElement;
}

inline void MaterialSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	fContinuumElement = p;
}

/* run time status */
inline const GlobalT::StateT MaterialSupportT::RunState(void) const
{
	if (fRunState) return *fRunState;
	else return GlobalT::kNone;
}

inline int MaterialSupportT::CurrIP(void) const
{
	if (fCurrIP) return *fCurrIP;
	else return 0;
}

inline int MaterialSupportT::IterationNumber(void) const
{
	if (fIterationNumber) return *fIterationNumber;
	else return -1;
}

inline double MaterialSupportT::Time(void) const
{
	if (fTime) return *fTime;
	else return 0.0;
}

inline double MaterialSupportT::TimeStep(void) const
{
	if (fTimeStep) return *fTimeStep;
	else return 0.0;
}

inline int MaterialSupportT::StepNumber(void) const
{
	if (fStepNumber) return *fStepNumber;
	else return -1;
}

inline void MaterialSupportT::SetRunState(const GlobalT::StateT& run_state)
{
	fRunState = &run_state;
}

inline void MaterialSupportT::SetCurrIP(const int& curr_ip)
{
	fCurrIP = &curr_ip;
}

inline void MaterialSupportT::SetIterationNumber(const int& iter)
{
	fIterationNumber = &iter;
}

inline void MaterialSupportT::SetTime(double& time)
{
	fTime = &time;
}

inline void MaterialSupportT::SetTimeStep(double& time_step)
{
	fTimeStep = &time_step;
}

inline void MaterialSupportT::SetStepNumber(int& step_number)
{
	fStepNumber = &step_number;
}

} /* namespace Tahoe */
#endif /* _SS_HOOKEAN_MAT_H_ */
