/* $Id: ParabolaT.h,v 1.6 2007-04-09 21:06:15 tdnguye Exp $ */
/* created: paklein (03/25/1999)                                          */

#ifndef _PARABOLA_T_H_
#define _PARABOLA_T_H_

/* base class */
#include "C1FunctionT.h"


namespace Tahoe {

class ParabolaT: public C1FunctionT
{
public:

        /* constructor */
        ParabolaT(double k, double B=0.0, double l0=0.0);
		ParabolaT(void);

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

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@{*/

private:

        /* potential parameters */
        double fk;
	double fl0;
	double fB;
};

/* inlines */

/* returning values */
inline double ParabolaT::Function(double x) const { 
  return 0.5*fk*(x-fl0)*(x-fl0)-0.5*fk*fB; }
inline double ParabolaT::DFunction(double x) const { return fk*(x-fl0); }
inline double ParabolaT::DDFunction(double x) const
{
#pragma unused(x)
        return fk;
}

} // namespace Tahoe 
#endif /* _PARABOLA_T_H_ */
