/* $Id: LanirFiber.h,v 1.1 2009-04-23 14:51:52 thao Exp $ */

#ifndef _LanirFiber_H_
#define _LanirFiber_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe {

class LanirFiber: public C1FunctionT
{
public:

	/* constructor */
	LanirFiber(double K, double alpha, double beta);
	LanirFiber(void);

	/* returning values */
	virtual double Function(double x) const;
	virtual double DFunction(double x) const;
	virtual double DDFunction(double x) const;

	/* returning values in groups */
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;

	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);

	/*derivatives of DFunction wrt to parameters*/
	virtual void ParameterDerivatives(double x, dArrayT& data) const;
	
	private:
	double PDF(double x) const;
	
private:
	/* parameters */
	double fK;  /*fiber stiffness*/
	/*gamma distribution*/
	double falpha; /*shape parameter*/
	double fbeta; /*scale parameter*/
	
	int fN;
};


} // namespace Tahoe 
#endif /* _LanirFiber_H_ */
