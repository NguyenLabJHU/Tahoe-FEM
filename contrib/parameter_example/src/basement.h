/* $Id: basement.h,v 1.1.2.1 2003-05-03 09:08:27 paklein Exp $ */
#ifndef _BASEMENT_H_
#define _BASEMENT_H_

/* base class */
#include "ParameterInterfaceT.h"

using namespace Tahoe;

class basement: public ParameterInterfaceT
{
public:

	basement(const StringT& name);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void SetParameters(const ParameterListT& list);
	/*@}*/

private:

	double height_;
	double length_;
	double width_;
};

#endif /* _BASEMENT_H_ */
