/* $Id: ParameterListT.h,v 1.1 2002-09-03 07:04:33 paklein Exp $ */
#ifndef _PARAMETER_LIST_T_H_
#define _PARAMETER_LIST_T_H_

/* direct members */
#include "ParameterT.h"
#include "AutoArrayT.h"

namespace Tahoe {

/** list of parameters */
class ParameterListT
{
public:

	/** constructor */
	ParameterListT(void);
	
protected:

	/** parameters name space */
	StringT fNamespace;

	/** list of parameters */
	AutoArrayT<ParameterT> fParameters;
};

} // namespace Tahoe 
#endif /* _PARAMETER_LIST_T_H_ */
