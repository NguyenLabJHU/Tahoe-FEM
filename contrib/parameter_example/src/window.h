/* $Id: window.h,v 1.2 2003-05-04 22:49:50 paklein Exp $ */
#ifndef _WINDOW_H_
#define _WINDOW_H_

/* base class */
#include "ParameterInterfaceT.h"

using namespace Tahoe;

class window: public ParameterInterfaceT
{
public:

	/** constructor */
	window(void);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	double width_;
	double height_;

};

#endif /* _WINDOW_H_ */
