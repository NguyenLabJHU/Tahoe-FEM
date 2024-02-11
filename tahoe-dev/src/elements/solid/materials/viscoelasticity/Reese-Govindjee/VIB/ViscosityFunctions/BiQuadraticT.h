/* $Id: BiQuadraticT.h,v 1.1 2006-10-30 23:32:05 thao Exp $ */
/* created: paklein (03/25/1999)                                          */

#ifndef _BIQUAD_T_H_
#define _BIQUAD_T_H_

/* base class */
#include "ViscFuncT.h"

namespace Tahoe {

class BiQuadraticT: public ViscFuncT
{
public:

	/* constructor */
	BiQuadraticT(double A1, double A2, double B, double C);

	/* I/O */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* returning values */
	virtual double Function(double Jv, double Je) const;
	virtual double DFuncDJv(double Jv, double Je) const;
	virtual double DFuncDJe(double Jv, double Je) const;

private:

	/* potential parameters */
	double fA1;
	double fA2;
	double fB;
	double fC;
};

/* inlines */

/* returning values */
inline double BiQuadraticT::Function(double Jv, double Je) const 
{
	double J=Je*Jv;
	if (J <= fB)
		return(fA1*(J-fB)*(J-fB)+fC);
	else
		return(fA2*(J-fB)*(J-fB)+fC);
}
inline double BiQuadraticT::DFuncDJv(double Jv, double Je) const 
{ 
        double J=Je*Jv;
	if (J <= fB)
		return(2*fA1*(J-fB)*Je);
	else
		return (2*fA2*(J-fB)*Je);

}
inline double BiQuadraticT::DFuncDJe(double Jv, double Je) const 
{ 
        double J=Je*Jv;
	if (J <= fB)
		return(2*fA1*(J-fB)*Jv);
	else
		return (2*fA2*(J-fB)*Jv);
}

}
#endif /* BIQUAD_T_H_ */
