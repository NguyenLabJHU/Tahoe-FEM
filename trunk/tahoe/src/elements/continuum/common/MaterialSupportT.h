/* $Id: MaterialSupportT.h,v 1.11 2004-06-26 05:54:58 paklein Exp $ */
#ifndef _MATERIAL_SUPPORT_T_H_
#define _MATERIAL_SUPPORT_T_H_

/* direct members */
#include "GlobalT.h"
#include "LocalArrayT.h"
#include "AutoArrayT.h"
#include "ElementCardT.h"

namespace Tahoe {

/* forward declarations */
class ContinuumElementT;
class ElementCardT;
class ScheduleT;
class ifstreamT;
class ofstreamT;
class CommunicatorT;

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

	/** \name multiprocessor support */
	/*@{*/
	/** the number of processes */
	int Size(void) const { return fSize; };

	/** the rank of this process */
	int Rank(void) const { return fRank; };

	/** the low-level global communicator, or NULL if it doesn't exist */
	const CommunicatorT* Communicator(void) const { return fCommunicator; };

	/** the low-level communicator only including processes with non-zero numbers
	 * of elements, or NULL if it doesn't exist */
	const CommunicatorT* GroupCommunicator(void) const { return fGroupCommunicator; };
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

	/** number of steps in the simulation for the current step size */
	int NumberOfSteps(void) const;

	/** set the source for the run state flag */
	void SetRunState(const GlobalT::StateT& run_state);

	/** set the source for the current stress evaluation point */
	void SetCurrIP(const int& curr_ip);

	/** set the source for the iteration number */
	void SetIterationNumber(const int& iter);

	/** set source for the current simulation time */
	void SetTime(const double& time);

	/** set source for the simulation time increment */
	void SetTimeStep(const double& time_step);

	/** set source for the simulation time increment number */
	void SetStepNumber(const int& step_number);

	/** set the source for the number of steps */
	void SetNumberOfSteps(const int& number_of_steps);
	/*@}*/
	
	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	const ContinuumElementT* ContinuumElement(void) const;
	
	/** set the source for element cards */
	void SetElementCards(AutoArrayT<ElementCardT>* element_cards);

	/** return the number of elements. If the element cards pointer
	 * is not set with MaterialSupportT::SetElementCards, this will return 0 */
	int NumElements(void) const;

	/** return the current element.  If the element cards pointer
	 * is not set with MaterialSupportT::SetElementCards, this will return -1 */
	int CurrElementNumber(void) const;

	/** return the specified card.  If the element cards pointer
	 * is not set with MaterialSupportT::SetElementCards, this will return NULL */
	//ElementCardT* ElementCard(int card);
	ElementCardT* ElementCard(int card) const;

	/** return the current card.  If the element cards pointer
	 * is not set with MaterialSupportT::SetElementCards, this will return NULL */
	ElementCardT* CurrentElement(void) const;

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

	/** \name set host code information */
	/*@{*/
	/** set the element group pointer */
	virtual void SetContinuumElement(const ContinuumElementT* p);

	/** set pointer local array */
	virtual void SetLocalArray(const LocalArrayT& array);
	/*@}*/

	/** \name input/output streams */
	/*@{*/
	/** the parameters stream */
	ifstreamT& Input(void) const;

	/** the echo file */
	ofstreamT& Output(void) const;
	/*@}*/

  private:
  
  	/** \name dimensions */
  	/*@{*/
	/** number of spatial dimensions */
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
	const int* fNumberOfSteps;
  	/*@}*/

	/** \name multiprocessor information */
	/*@{*/
	int fSize;
	int fRank;
	
	/** global communicator */
	const CommunicatorT* fCommunicator;

	/** communicator including only processes with non-zero numbers of elements */
	const CommunicatorT* fGroupCommunicator;
	/*@}*/

	/** pointer to element card information */
	AutoArrayT<ElementCardT>* fElementCards;	
  
  	/** pointer to the continuum element */
  	const ContinuumElementT* fContinuumElement;

	/** \name pointers to local arrays */
	/*@{*/
	const LocalArrayT* fInitCoords;
	const LocalArrayT* fDisp;
	/*@}*/
};

/* inlines functions */
inline const ContinuumElementT* MaterialSupportT::ContinuumElement(void) const
{
	return fContinuumElement;
}

/* set the source for element cards */
inline void MaterialSupportT::SetElementCards(AutoArrayT<ElementCardT>* element_cards)
{
	fElementCards = element_cards;
}

/* return the number of elements */
inline int MaterialSupportT::NumElements(void) const
{
	if (fElementCards) 
		return fElementCards->Length();
	else
		return 0;
}

/* return the current element */
inline int MaterialSupportT::CurrElementNumber(void) const
{
	if (fElementCards) 
		return fElementCards->Position();
	else
		return -1;
}

/* return the specified card */
inline ElementCardT* MaterialSupportT::ElementCard(int card) const
{
	if (fElementCards) 
		return fElementCards->Pointer(card);
	else
		return NULL;
}

/* return the current */
inline ElementCardT* MaterialSupportT::CurrentElement(void) const
{
	if (fElementCards && fElementCards->InRange()) 
		return &(fElementCards->Current());
	else
		return NULL;
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

inline int MaterialSupportT::NumberOfSteps(void) const
{
	if (fNumberOfSteps) return *fNumberOfSteps;
	else return 0;
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

inline void MaterialSupportT::SetTime(const double& time)
{
	fTime = &time;
}

inline void MaterialSupportT::SetTimeStep(const double& time_step)
{
	fTimeStep = &time_step;
}

inline void MaterialSupportT::SetStepNumber(const int& step_number)
{
	fStepNumber = &step_number;
}

inline void MaterialSupportT::SetNumberOfSteps(const int& number_of_steps)
{
	fNumberOfSteps = &number_of_steps;
}

} /* namespace Tahoe */
#endif /* _SS_HOOKEAN_MAT_H_ */
