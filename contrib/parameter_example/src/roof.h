/* $Id: roof.h,v 1.1.2.2 2003-05-03 09:08:27 paklein Exp $ */
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
	virtual void SetParameters(const ParameterListT& list);
	/*@}*/

private:

	/** roof style as enumerated type */
	style style_;
};

#endif /* _HOUSE_H_ */
