/* $Id: Triantafyllidis.h,v 1.1.6.1 2002-06-27 18:00:41 cjkimme Exp $ */

#ifndef _TRIANTAFYLLIDIS_H_
#define _TRIANTAFYLLIDIS_H_

/* base class */
#include "C1FunctionT.h"

/** the potential used by Triantafyllidis and Bardenhagen in \a Journal \a of
 * \a Elasticity (1993) */

namespace Tahoe {

class Triantafyllidis: public C1FunctionT
{
public:

	/** constructor */
	Triantafyllidis(double A);

	/** I/O */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* function evaluation */
	virtual double Function(double x) const;   /**< the function */
	virtual double DFunction(double x) const;  /**< the first derivative */
	virtual double DDFunction(double x) const; /**< the second derivative */

	/* returning values in groups. derived classes should define
	 * their own non-virtual function called within this functon
	 * which maps in to out w/o requiring a virtual function call
	 * everytime. Default behavior is just to map the virtual functions
	 * above */
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;

private:

	/* scaling parameter */
	double fA;
};

} // namespace Tahoe 
#endif /* _TRIANTAFYLLIDIS_H_ */
