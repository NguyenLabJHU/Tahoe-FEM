/* $Id: crawl_space.h,v 1.2 2003-05-04 22:49:50 paklein Exp $ */
#ifndef _CRAWLSPACE_H_
#define _CRAWLSPACE_H_

/* base class */
#include "basement.h"

using namespace Tahoe;

class crawl_space: public basement
{
public:

	/** constructor */
	crawl_space(void);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	bool sump_pump_;
};

#endif /* _CRAWLSPACE_H_ */
