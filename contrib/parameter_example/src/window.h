/* $Id: window.h,v 1.1.2.1 2003-05-03 09:08:27 paklein Exp $ */
#ifndef _WINDOW_H_
#define _WINDOW_H_

/* base class */
#include "ParameterInterfaceT.h"

using namespace Tahoe;

class window: public ParameterInterfaceT
{
public:

	window(void);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void SetParameters(const ParameterListT& list);
	/*@}*/

private:

	double width_;
	double height_;

};

#endif /* _WINDOW_H_ */
