/* $Id: FBC_CardT.h,v 1.5 2002-07-05 22:28:33 paklein Exp $ */
/* created: paklein (06/15/1996) */

#ifndef _FBC_CARD_T_H_
#define _FBC_CARD_T_H_

#include "Environment.h"
#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class NodeManagerT;
class ifstreamT;
class ScheduleT;

/** nodal force boundary condition information */
class FBC_CardT
{
public:
	
	/* constructor */
	FBC_CardT(void);

	/* modifiers */
	void SetValues(const NodeManagerT& theBoss, ifstreamT& in);
	void SetValues(const NodeManagerT& theBoss, int node, int dof, int schedule,
		double value);
	void SplitForce(void);
	
	/* return the node and DOF number specified for the force */
	void Destination(int& node, int& dof) const;
	int Node(void) const;
	int DOF(void) const;
	
	/* return the current value */
	double CurrentValue(void) const;

	/* output */
	static void WriteHeader(ostream& out);
	void WriteValues(ostream& out) const;

private:
	
	int fNode; /* need node number and dof number b/c */
	int fDOF;  /* fDestination not set at input time  */
	int fSchedNum;
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
} // namespace Tahoe 
#endif /* _FBC_CARD_T_H_ */
