/* $Id: KBC_PrescribedT.h,v 1.4 2004-07-22 08:32:01 paklein Exp $ */
#ifndef _KBC_PRESCRIBED_T_H_
#define _KBC_PRESCRIBED_T_H_

/* base class */
#include "KBC_ControllerT.h"

namespace Tahoe {

/** controller to insert arbitrary kinematic boundary conditions.
 * Motion is prescribed by creating  */
class KBC_PrescribedT: public KBC_ControllerT
{
public:

	/** constructor */
	KBC_PrescribedT(const BasicSupportT& support);

	/** initialize data */
	virtual void Initialize(ifstreamT& in);

	/** set to initial conditions */
	virtual void InitialCondition(void);
};

} /* namespace Tahoe */

#endif /* _KBC_PRESCRIBED_T_H_ */
