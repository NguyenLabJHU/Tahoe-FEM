/* $Id: FBC_CardT.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
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
	
	/* return the node and DOF number specified for the force */
	void Destination(int& node, int& dof) const;

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

#endif /* _FBC_CARD_T_H_ */
