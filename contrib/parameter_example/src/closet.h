/* $Id: closet.h,v 1.1.2.2 2003-05-03 09:08:27 paklein Exp $ */
#ifndef _CLOSET_H_
#define _CLOSET_H_

/* base class */
#include "room.h"

using namespace Tahoe;

class closet: public room
{
public:

	closet(void);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void SetParameters(const ParameterListT& list);
	/*@}*/

private:

	bool has_shelf_;
	bool has_bar_;
};

#endif /* _CLOSET_H_ */
