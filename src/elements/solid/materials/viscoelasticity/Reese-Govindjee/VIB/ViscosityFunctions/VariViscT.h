//* $Id: VariViscT.h,v 1.1 2006-10-30 23:32:05 thao Exp $ */
/* created: paklein (03/25/1999)                                          */

#ifndef _VariVisc_T_H_
#define _VariVisc_T_H_

/* base class */
#include <math.h>
#include <iostream.h>

#include "ViscFuncT.h"

namespace Tahoe {

class VariViscT: public ViscFuncT
{
public:

	/* constructor */
	VariViscT(double n,double d,double z, double Jcrit, double Jo=1.00);

	/* I/O */
	void Print(ostream& out) const;
	void PrintName(ostream& out) const;
	
	/* returning values */
	double Function(double Jv, double Je) const;
	double DFuncDJv(double Jv, double Je) const;
	double DFuncDJe(double Jv, double Je) const;
private:

	/* potential parameters */
	double fno;
	double fdo;
	double fJo;
	double fz;
	double fJcrit;
};

/* inlines */

/* returning values */
inline double VariViscT::Function(double Jv, double Je) const 
{
  //        double error = 1e-5*fno;
        double error = 1e-3*fno;
	double d = fdo*(1+exp((Je-fJcrit)/fz));
	double fun = fno/exp((Jv-fJo)/d)+error;
	return fun; 
}

inline double VariViscT::DFuncDJv(double Jv, double Je) const
{
        double d = fdo*(1+exp((Je-fJcrit)/fz));
	double dfun = -fno/d/exp((Jv-fJo)/d);
	return (dfun);	

}

inline double VariViscT::DFuncDJe(double Jv, double Je) const
{
        double d = fdo*(1+exp((Je-fJcrit)/fz));
	double dfun = fno/d/exp((Jv-fJo)/d);
	return (dfun*(Jv-fJo)*(1-fdo/d)/fz);
}    
}

#endif /* VariVisc_T_H_ */
