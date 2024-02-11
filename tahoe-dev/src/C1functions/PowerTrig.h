/* $Id: PowerTrig.h,v 1.1 2006-08-18 18:43:15 thao Exp $ */

#ifndef _PowerTrig_H_
#define _PowerTrig_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class PowerTrig: public C1FunctionT
{
public:

	/* constructor */
	PowerTrig(double A, double B, double C, double n, double phi);
	PowerTrig(void);
	
	/* returning values */
	virtual double Function(double x) const;
	virtual double DFunction(double x) const;
	virtual double DDFunction(double x) const;

	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);

private:
	/* parameters */
	double fA;
	double fB;
	double fC;
	double fphi;
	double fn;
};


} // namespace Tahoe 
#endif /* _PowerTrig_H_ */
