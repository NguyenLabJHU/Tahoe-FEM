/* $Id: storm_shelter.h,v 1.1.2.1 2003-05-03 09:08:27 paklein Exp $ */
#ifndef _STORM_SHELTER_H_
#define _STORM_SHELTER_H_

/* base class */
#include "basement.h"

using namespace Tahoe;

class storm_shelter: public basement
{
public:

	storm_shelter(void);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void SetParameters(const ParameterListT& list);
	/*@}*/

private:

	bool auxiliary_power_;
	bool first_aid_kit_;
	double stored_water_;
};

#endif /* _STORM_SHELTER_H_ */
