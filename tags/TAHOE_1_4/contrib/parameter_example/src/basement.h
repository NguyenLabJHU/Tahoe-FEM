/* $Id: basement.h,v 1.2 2003-05-04 22:49:50 paklein Exp $ */
#ifndef _BASEMENT_H_
#define _BASEMENT_H_

/* base class */
#include "ParameterInterfaceT.h"

using namespace Tahoe;

class basement: public ParameterInterfaceT
{
public:

	/** constructor */
	basement(const StringT& name);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	double height_;
	double length_;
	double width_;
};

#endif /* _BASEMENT_H_ */
