/* $Id: ParameterInterfaceT.h,v 1.1 2003-04-22 18:32:16 paklein Exp $ */
#ifndef _PARAMETER_INTERFACE_T_H_
#define _PARAMETER_INTERFACE_T_H_

namespace Tahoe {

/* forward declarations */
class ParameterListT;
class StringT;

/** abstract interface for classes which define and use parameters */
class ParameterInterfaceT
{
public:

	/** constructor */
	ParameterInterfaceT(void) {};

	/** identifier */
	virtual const StringT& Name(void) const = 0;

	/** accept completed parameter list */
	virtual void SetParameters(const ParameterListT& list);
	
	/** build complete parameter list description */
	virtual void DefineParameters(ParameterListT& list);
};

} /* namespace Tahoe */

#endif /* _PARAMETER_SOURCE_T_H_ */
