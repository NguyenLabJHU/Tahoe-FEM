/* $Id: FBC_CardT.h,v 1.7 2008-05-12 22:29:16 regueiro Exp $ */
/* created: paklein (06/15/1996) */
#ifndef _FBC_CARD_T_H_
#define _FBC_CARD_T_H_

#ifdef __DEVELOPMENT__
#include "DevelopmentElementsConfig.h"
#endif

namespace Tahoe {

/* forward declarations */
class ScheduleT;

/** nodal force boundary condition information */
class FBC_CardT
{
public:
	
	/* constructor */
	FBC_CardT(void);

	/* modifiers */
	void SetValues(int node, int dof, const ScheduleT* schedule, double value);
	void SplitForce(void);

#ifdef DEM_COUPLING_DEV
	void AddValues(int node, int dof, const ScheduleT* schedule, double value);
	void ClearValues();
#endif

	/* return the node and DOF number specified for the force */
	void Destination(int& node, int& dof) const;
	int Node(void) const;
	int DOF(void) const;
	double Value(void) const;
	const ScheduleT* Schedule(void) const { return fSchedule; };
	
	/* return the current value */
	double CurrentValue(void) const;

private:
	
	int fNode; /* need node number and dof number b/c */
	int fDOF;  /* fDestination not set at input time  */
	double fValue;				

	const ScheduleT* fSchedule;		
};

/* inlines */

/* return the node and DOF number specified for the force */
inline void FBC_CardT::Destination(int& node, int& dof) const
{
	node = fNode;
	dof  = fDOF;
}
inline int FBC_CardT::Node(void) const   { return fNode; }
inline int FBC_CardT::DOF(void) const    { return fDOF;  }

#ifdef DEM_COUPLING_DEV
inline double FBC_CardT::Value(void) const    { return fValue;  }
#endif

} // namespace Tahoe 
#endif /* _FBC_CARD_T_H_ */
