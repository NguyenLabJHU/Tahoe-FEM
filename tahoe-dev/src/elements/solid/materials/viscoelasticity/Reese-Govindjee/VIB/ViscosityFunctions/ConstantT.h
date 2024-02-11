/* $Id: ConstantT.h,v 1.1 2006-10-30 23:32:05 thao Exp $ */
/* created: paklein (03/25/1999)                                          */

#ifndef _CONST_T_H_
#define _CONST_T_H_

/* base class */
#include "ViscFuncT.h"

namespace Tahoe {

class ConstantT: public ViscFuncT
{
public:

	/* constructor */
	ConstantT(double A);

	/* I/O */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* returning values */
	virtual double Function(double Jv, double Je) const;
	virtual double DFuncDJv(double Jv, double Je) const;
	virtual double DFuncDJe(double Jv, double Je) const;

private:

	/* potential parameters */
	double fA;};

/* inlines */

/* returning values */
inline double ConstantT::Function(double Jv, double Je) const 
{ 
  #pragma unused(Jv, Je)
  return fA; 
}
inline double ConstantT::DFuncDJv(double Jv, double Je) const 
{ 
  #pragma unused(Jv, Je)
  return (0.0); 
}
inline double ConstantT::DFuncDJe(double Jv, double Je) const
{
  #pragma unused(Jv, Je)
  return (0.0);
}

}
#endif /* _CONST_T_H_ */



