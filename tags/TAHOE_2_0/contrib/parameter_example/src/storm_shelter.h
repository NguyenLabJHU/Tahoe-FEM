/* $Id: storm_shelter.h,v 1.2 2003-05-04 22:49:50 paklein Exp $ */
#ifndef _STORM_SHELTER_H_
#define _STORM_SHELTER_H_

/* base class */
#include "basement.h"

using namespace Tahoe;

class storm_shelter: public basement
{
public:

	/** constructor */
	storm_shelter(void);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	bool auxiliary_power_;
	bool first_aid_kit_;
	double stored_water_;
};

#endif /* _STORM_SHELTER_H_ */
