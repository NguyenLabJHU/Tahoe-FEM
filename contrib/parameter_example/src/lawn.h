/* $Id: lawn.h,v 1.1.2.2 2003-05-03 09:08:27 paklein Exp $ */
#ifndef _LAWN_H_
#define _LAWN_H_

/* base class */
#include "ParameterInterfaceT.h"

using namespace Tahoe;

class lawn: public ParameterInterfaceT
{
public:

	lawn(void);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void SetParameters(const ParameterListT& list);
	/*@}*/

private:

};

#endif /* _LAWN_H_ */
