/* $Id: KBC_PrescribedT.h,v 1.1.2.1 2003-02-10 02:12:23 paklein Exp $ */
#ifndef _KBC_PRESCRIBED_T_H_
#define _KBC_PRESCRIBED_T_H_

/* base class */
#include "KBC_ControllerT.h"

namespace Tahoe {

/** controller to insert arbitrary kinematic boundary conditions */
class KBC_PrescribedT: public KBC_ControllerT
{
public:

	/** constructor */
	KBC_PrescribedT(NodeManagerT& node_manager);

	/** initialize data */
	virtual void Initialize(ifstreamT& in);

	/** set to initial conditions */
	virtual void InitialCondition(void);
};

} /* namespace Tahoe */

#endif /* _KBC_PRESCRIBED_T_H_ */
