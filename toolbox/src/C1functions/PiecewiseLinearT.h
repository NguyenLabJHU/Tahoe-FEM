/* $Id: PiecewiseLinearT.h,v 1.2 2004-07-12 21:48:02 paklein Exp $ */
#ifndef _PIECEWISE_LINEAR_T_H_
#define _PIECEWISE_LINEAR_T_H_

/* base class */
#include "C1FunctionT.h"

/* direct members */
#include "dRangeArrayT.h"

namespace Tahoe {

/** interface for a piecewise linear function */
class PiecewiseLinearT: public C1FunctionT
{
public:

	/** \name constructors */
	/*@{*/
	PiecewiseLinearT(void);	
	PiecewiseLinearT(const dArray2DT& points);
	/*@}*/

	/* I/O */
	virtual void Print(ostream& out) const;     	    	   	
	virtual void PrintName(ostream& out) const;     	    	   	
	    	   	    	
	/** \name returning values */
	/*@{*/
	virtual double Function(double x) const;
	virtual double DFunction(double x) const;
	virtual double DDFunction(double x) const;
	/*@}*/

	/* returning values in groups - returns refence to out to allow:
	 *
	 *	dArrayT& goodname = pfunc->MapFunction(in, tempspace);
	 */
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;

	/** Return 0th, 1st, and 2nd derivative in the respective fields of the dArrayT */  	
	virtual void SetAll(double x, dArrayT& data) const;   	

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists
	 * \param sub_lists description of subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** accept parameter list.
	 * \param list input parameter list, which should be validated using ParameterInterfaceT::ValidateParameterList
	 *        to ensure the list conforms to the description defined by the interface. */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@{*/

protected:

	/* non-virtual function calls */
	double function(double x) const;
	double Dfunction(double x) const;
	double DDfunction(double x) const;
	void all_functions(double x, double& f, double& Df, double& DDf) const;

private:

	/** make sure knots are valid */
	void CheckKnots(const dRangeArrayT& x, const dArrayT& y) const;

protected:   	

	dRangeArrayT fXPoints;	
	dArrayT      fYPoints;	
};

} /* namespace Tahoe */

#endif /* _PIECEWISE_LINEAR_T_H_ */
