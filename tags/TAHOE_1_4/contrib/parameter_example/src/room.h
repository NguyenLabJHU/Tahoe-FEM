/* $Id: room.h,v 1.2 2003-05-04 22:49:50 paklein Exp $ */
#ifndef _ROOM_H_
#define _ROOM_H_

/* base class */
#include "ParameterInterfaceT.h"

using namespace Tahoe;

class room: public ParameterInterfaceT
{
public:

	/** constructor */
	room(const StringT& name);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/** \name dimensions */
	/*@{*/
	int length;
	int width;
	/*@}*/
};

#endif /* _ROOM_H_ */
