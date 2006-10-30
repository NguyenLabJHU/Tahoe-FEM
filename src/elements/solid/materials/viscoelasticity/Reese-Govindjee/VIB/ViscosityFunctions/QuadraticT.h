/* $Id: QuadraticT.h,v 1.1 2006-10-30 23:32:05 thao Exp $ */
/* created: paklein (03/25/1999)                                          */

#ifndef _QUAD_T_H_
#define _QUAD_T_H_

/* base class */
#include "ViscFuncT.h"

namespace Tahoe {

class QuadraticT: public ViscFuncT
{
public:

	/* constructor */
	QuadraticT(double A, double B, double C);

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
	double fC;
};

/* inlines */

/* returning values */
inline double QuadraticT::Function(double Jv, double Je) const 
{
  double J = Je*Jv;
  return (fA*(J-fB)*(J-fB)+fC); 
}
inline double QuadraticT::DFuncDJv(double Jv, double Je) const 
{ 
  return (2*fA*(Je*Jv-fB)*Je); 
}
inline double QuadraticT::DFuncDJe(double Jv, double Je) const
{
  return (2*fA*(Je*Jv-fB)*Jv); 
}
}
#endif /* _QUAD_T_H_ */
