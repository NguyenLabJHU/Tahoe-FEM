/* $Id: NLHHTalpha.h,v 1.5 2004-07-15 08:30:28 paklein Exp $ */
/* created: paklein (10/11/1996) */
#ifndef _NL_HHT_ALPHA_H_
#define _NL_HHT_ALPHA_H_

/* base classes */
#include "nNLHHTalpha.h"
#include "eNLHHTalpha.h"

namespace Tahoe {

/* forward declarations */
class NodeManagerT;
class TimeManagerT;

class NLHHTalpha: public nNLHHTalpha, public eNLHHTalpha
{
public:

	/** constructor */
	NLHHTalpha(double alpha);

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

	TimeManagerT* 	fTimeBoss;
	double			fTimeShift; /* t_n+1+alpha */
	
};

} // namespace Tahoe 
#endif /* _NL_HHT_ALPHA_H_ */
