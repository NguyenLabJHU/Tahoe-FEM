/* $Id: FieldSupportT.h,v 1.5.16.4 2004-03-31 16:14:37 paklein Exp $ */
#ifndef _FIELD_SUPPORT_T_H_
#define _FIELD_SUPPORT_T_H_

/* base class */
#include "BasicSupportT.h"

namespace Tahoe {

/* forward declarations */
class FBC_ControllerT;
class KBC_ControllerT;

/** support for FieldT. Provides a limited interface to get 
 * information in and out of FieldT's. */
class FieldSupportT: public BasicSupportT
{
public:

	/** constructor */
	FieldSupportT(const FEManagerT& fe);

	/** \name construct BC controllers
	 * Construct new kinematic or force boundary condition controllers. Responsibility 
	 * for deleteting instantiations resides with the client who requested them.
	 */
	/*@{*/
	KBC_ControllerT* NewKBC_Controller(FieldT& field, int code) const;
	FBC_ControllerT* NewFBC_Controller(FieldT& field, int code) const;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _FIELD_SUPPORT_T_H_ */
