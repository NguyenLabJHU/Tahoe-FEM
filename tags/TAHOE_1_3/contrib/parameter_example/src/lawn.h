/* $Id: lawn.h,v 1.2 2003-05-04 22:49:50 paklein Exp $ */
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
