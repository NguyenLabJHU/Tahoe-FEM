/* $Id: roof.h,v 1.2 2003-05-04 22:49:50 paklein Exp $ */
#ifndef _ROOF_H_
#define _ROOT_H_

/* base class */
#include "ParameterInterfaceT.h"

using namespace Tahoe;

class roof: public ParameterInterfaceT
{
public:

	enum style {
		undefined,
		shingle,
		slate
	};

	roof(void);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/** roof style as enumerated type */
	style style_;
};

#endif /* _HOUSE_H_ */
