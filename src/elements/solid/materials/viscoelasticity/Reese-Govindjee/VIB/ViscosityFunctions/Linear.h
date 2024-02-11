/* $Id: Linear.h,v 1.1 2006-10-30 23:32:05 thao Exp $ */

#ifndef _LINEAR_H_
#define _LINEAR_H_

/* base class */
#include "ViscFuncT.h"

namespace Tahoe {

class Linear: public ViscFuncT
{
public:

	/* constructor */
	Linear(double A, double B);

	/* I/O */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* returning values */
	virtual double Function(double Jv, double Je) const;
	virtual double DFuncDJv(double Jv, double Je) const;
	virtual double DFuncDJe(double Jv, double Je) const;

private:

	/* potential parameters */
	double fA;
	double fB;
};

/* inlines */

/* returning values */
inline double Linear::Function(double Jv, double Je) const 
{ 
  double J = Je*Jv;
  return (fA*J+fB); 
}

inline double Linear::DFuncDJv(double Jv, double Je) const 
{ 
#pragma unused(Jv)
  return (fA*Je); 
}
inline double Linear::DFuncDJe(double Jv, double Je) const
{ 
#pragma unused(Je)
  return (fA*Jv); 
}
}
#endif /* _LINEAR_T_H_ */



