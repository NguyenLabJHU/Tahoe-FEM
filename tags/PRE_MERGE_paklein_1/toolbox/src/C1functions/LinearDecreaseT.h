/* $Id: LinearDecreaseT.h,v 1.2 2002-07-02 19:56:31 cjkimme Exp $ */
/* created: paklein (03/25/1999)                                          */

#ifndef _LINDEC_T_H_
#define _LINDEC_T_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class LinearDecreaseT: public C1FunctionT
{
public:

	/* constructor */
	LinearDecreaseT(double A, double L);

	/* I/O */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* returning values */
	virtual double Function(double x) const;
	virtual double DFunction(double x) const;
	virtual double DDFunction(double x) const;

	/* returning values in groups */
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;

private:

	/* potential parameters */
	double fA;
	double fL;
};

/* inlines */

/* returning values */
inline double LinearDecreaseT::Function(double x) const 
{ return ((1-fA*x/fL)); }
inline double LinearDecreaseT::DFunction(double x) const 
{ 
  #pragma unused(x)
  return (-fA/fL); 
}
inline double LinearDecreaseT::DDFunction(double x) const
{
#pragma unused(x)
	return (0);
}

} // namespace Tahoe 
#endif /* _LINDEC_T_H_ */



