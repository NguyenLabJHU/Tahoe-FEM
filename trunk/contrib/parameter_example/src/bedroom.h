/* $Id: bedroom.h,v 1.3 2003-08-14 01:22:43 paklein Exp $ */
#ifndef _BEDROOM_H_
#define _BEDROOM_H_

/* base class */
#include "room.h"

/* forward declarations */
class window;

using namespace Tahoe;

class bedroom: public room
{
public:

	/** constructor */
	bedroom(void);

	/** destructor */
	~bedroom(void);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);

	virtual void DefineSubs(SubListT& sub_list) const;
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;
	/*@}*/

private:

	int floor;
	
	/** windows */
	ArrayT<window*> windows_;
};

#endif /* _BEDROOM_H_ */
