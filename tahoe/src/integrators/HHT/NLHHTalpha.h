/* $Id: NLHHTalpha.h,v 1.2 2002-04-02 23:19:20 paklein Exp $ */
/* created: paklein (10/11/1996) */

#ifndef _NL_HHT_ALPHA_H_
#define _NL_HHT_ALPHA_H_

/* base classes */
#include "nNLHHTalpha.h"
#include "eNLHHTalpha.h"

/* forward declarations */
class NodeManagerT;
class TimeManagerT;

class NLHHTalpha: public nNLHHTalpha, public eNLHHTalpha
{
public:

	/** constructor */
	NLHHTalpha(TimeManagerT& TM, ifstreamT& in, ostream& out, bool auto2ndorder);

	/* take responsibility for forming the nodal contribution
	 * to the RHS vector:
	 *
	 *                     F(t_n+1+alpha)
	 */
	virtual void FormNodalForce(NodeManagerT* nodeboss) const;
	  	
protected:  	
	
	/* recalculate time stepping constants */
	virtual void ComputeParameters(void);

private:

	TimeManagerT& 	TimeBoss;
	double			fTimeShift; /* t_n+1+alpha */
	
};

#endif /* _NL_HHT_ALPHA_H_ */
