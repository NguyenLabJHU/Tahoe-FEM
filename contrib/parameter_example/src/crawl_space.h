/* $Id: crawl_space.h,v 1.1.2.1 2003-05-03 09:08:27 paklein Exp $ */
#ifndef _CRAWLSPACE_H_
#define _CRAWLSPACE_H_

/* base class */
#include "basement.h"

using namespace Tahoe;

class crawl_space: public basement
{
public:

	crawl_space(void);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void SetParameters(const ParameterListT& list);
	/*@}*/

private:

	bool sump_pump_;
};

#endif /* _CRAWLSPACE_H_ */
