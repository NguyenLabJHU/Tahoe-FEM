/* $Id: lawn.h,v 1.1.2.3 2003-05-04 22:13:39 paklein Exp $ */
#ifndef _LAWN_H_
#define _LAWN_H_

/* base class */
#include "ParameterInterfaceT.h"

using namespace Tahoe;

class lawn: public ParameterInterfaceT
{
public:

	/** constructor */
	lawn(void);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

};

#endif /* _LAWN_H_ */
