/* $Id: FBC_CardT.h,v 1.2 2001-07-16 22:41:02 rrsettg Exp $ */
/* created: paklein (06/15/1996)                                          */
/* Adds direct link to RHS vector to speed calculation of                 */
/* nodal contribution to the residual force vector.                       */

#ifndef _FBC_CARD_T_H_
#define _FBC_CARD_T_H_

#include "Environment.h"
#include "ios_fwd_decl.h"

/* forward declarations */
class NodeManagerPrimitive;
class ifstreamT;
class LoadTime;

class FBC_CardT
{
public:
	
	/* constructor */
	FBC_CardT(void);

	/* modifiers */
	void SetValues(const NodeManagerPrimitive& theBoss, ifstreamT& in);
	void SetValues(const NodeManagerPrimitive& theBoss, int node, int dof, int nLTf,
		double value);
	void SplitForce(void);
	
	/* return the node and DOF number specified for the force */
	void Destination(int& node, int& dof) const;
	int Node(void) const;
	int DOF(void) const;
	
	/* return the current value */
	double CurrentValue(void) const;

	/* output */
	void WriteHeader(ostream& out) const;
	void WriteValues(ostream& out) const;

private:
	
	int fNode; /* need node number and dof number b/c */
	int fDOF;  /* fDestination not set at input time  */
	int fLTf;
	double fValue;				

	const LoadTime* fLTfPtr;		
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
#endif /* _FBC_CARD_T_H_ */
