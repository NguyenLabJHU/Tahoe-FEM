/* $Id: driveway.h,v 1.1.2.2 2003-05-03 09:08:27 paklein Exp $ */
#ifndef _DRIVEWAY_H_
#define _DRIVEWAY_H_

/* base class */
#include "ParameterInterfaceT.h"

using namespace Tahoe;

class driveway: public ParameterInterfaceT
{
public:

	enum surface {
		undefined,
		dirt,
		gravel,
		cobblestone,
		asphalt,
		concrete
	};

	driveway(void);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void SetParameters(const ParameterListT& list);
	/*@}*/

private:

	/** length of the driveway */
	double length_;

	/** driveway surface as enumerated type */
	surface surface_;
};

#endif /* _DRIVEWAY_H_ */
