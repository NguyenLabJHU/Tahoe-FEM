/* $Id: ModSmithFerrante.h,v 1.4 2004-06-07 21:37:33 paklein Exp $ */
#ifndef _MOD_SMITH_FERRANTE_H_
#define _MOD_SMITH_FERRANTE_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe {

/** SmithFerrante potential modified with a quadratic response for x < 0 */
class ModSmithFerrante: public C1FunctionT
{
public:

	/** \name constructors */
	/*@{*/
	ModSmithFerrante(void);
	ModSmithFerrante(double A, double B);
	/*@}*/

	/** \name I/O */
	/*@{*/
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	/*@}*/
	
	/** \name returning values */
	/*@{*/
	virtual double Function(double x) const;
	virtual double DFunction(double x) const;
	virtual double DDFunction(double x) const;
	/*@}*/

	/** \name returning values in groups */
	/*@{*/
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/** \name potential parameters */
	/*@{*/
	double fA;
	double fB;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _MOD_SMITH_FERRANTE_H_ */
