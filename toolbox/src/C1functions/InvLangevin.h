/* $Id: InvLangevin.h,v 1.1 2007-03-08 18:43:48 tdnguye Exp $ */

#ifndef _InvLangevin_H_
#define _InvLangevin_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class InvLangevin: public C1FunctionT
{
public:

	/* constructor */
	InvLangevin(void);
	
	/* returning values */
	virtual double Function(double x) const;
	virtual double DFunction(double x) const;
	virtual double DDFunction(double x) const;

	/* returning values in groups */
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;

};


} // namespace Tahoe 
#endif /* _InvLangevin_H_ */
