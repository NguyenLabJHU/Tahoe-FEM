/* $Id: FBC_CardT.h,v 1.6.20.1 2005-06-08 17:21:42 paklein Exp $ */
/* created: paklein (06/15/1996) */
#ifndef _FBC_CARD_T_H_
#define _FBC_CARD_T_H_

/* direct members */
#include "StringT.h"

namespace Tahoe {

/* forward declarations */
class ScheduleT;

/** nodal force boundary condition information */
class FBC_CardT
{
public:

	/** usage mode */
	enum ModeT {
		kUndefined,
		kNode,
		kSet
	};
	
	/* constructor */
	FBC_CardT(void);

	/* modifiers */
	void SetValues(int node, int dof, const ScheduleT* schedule, double value);
	void SetValues(const StringT& ID, int dof, const ScheduleT* schedule, double value);
	void SplitForce(void);
	
	/* return the node and DOF number specified for the force */
	void Destination(int& node, int& dof) const;
	int Node(void) const;
	const StringT& ID(void) const { return fID; };
	int DOF(void) const;
	ModeT Mode(void) const { return fMode; };
	const ScheduleT* Schedule(void) const { return fSchedule; };
	
	/* return the current value */
	double CurrentValue(void) const;

private:
	
	/** \name node or node set */
	/*@{*/
	int fNode;
	StringT fID;
	/*@}*/

	int fDOF;  /* fDestination not set at input time  */
	double fValue;
	ModeT  fMode;

	const ScheduleT* fSchedule;		
};

/* inlines */

/* return the node and DOF number specified for the force */
inline void FBC_CardT::Destination(int& node, int& dof) const {
	node = fNode;
	dof  = fDOF;
}
inline int FBC_CardT::Node(void) const   { return fNode; }
inline int FBC_CardT::DOF(void) const    { return fDOF;  }
} // namespace Tahoe 
#endif /* _FBC_CARD_T_H_ */
