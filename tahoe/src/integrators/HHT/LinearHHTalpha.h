/* $Id: LinearHHTalpha.h,v 1.4 2002-07-05 22:27:52 paklein Exp $ */
/* created: paklein (10/11/1996) */

#ifndef _LINEAR_HHT_ALPHA_H_
#define _LINEAR_HHT_ALPHA_H_

/* base classes */
#include "nLinearHHTalpha.h"
#include "eLinearHHTalpha.h"

namespace Tahoe {

/* forward declarations */
class NodeManagerT;
class TimeManagerT;

class LinearHHTalpha: public nLinearHHTalpha, public eLinearHHTalpha
{
public:

	/** constructor */
	LinearHHTalpha(TimeManagerT& TM, ifstreamT& in, ostream& out, bool auto2ndorder);

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

} // namespace Tahoe 
#endif /* _LINEAR_HHT_ALPHA_H_ */
