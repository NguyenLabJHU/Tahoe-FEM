/* $Id: room.h,v 1.1.2.2 2003-05-03 09:08:27 paklein Exp $ */
#ifndef _ROOM_H_
#define _ROOM_H_

/* base class */
#include "ParameterInterfaceT.h"

using namespace Tahoe;

class room: public ParameterInterfaceT
{
public:

	room(const StringT& name);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void SetParameters(const ParameterListT& list);
	/*@}*/

private:

	/** \name dimensions */
	/*@{*/
	int length;
	int width;
	/*@}*/
};

#endif /* _ROOM_H_ */
