/* $Id: LinearHHTalpha.h,v 1.4.56.1 2004-04-08 07:33:37 paklein Exp $ */
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
	LinearHHTalpha(double alpha);

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

	TimeManagerT* fTimeBoss;
	double        fTimeShift; /* t_n+1+alpha */
	
};

} // namespace Tahoe 
#endif /* _LINEAR_HHT_ALPHA_H_ */
