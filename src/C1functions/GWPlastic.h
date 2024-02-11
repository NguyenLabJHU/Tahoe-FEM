/* $Id: GWPlastic.h,v 1.6 2003-11-21 22:54:23 paklein Exp $ */
#ifndef _GW_PLASTIC_H_
#define _GW_PLASTIC_H_

/* base class */
#include "C1FunctionT.h"

namespace Tahoe {

class GWPlastic: public C1FunctionT
{
public:

	/*
	 * Constructor
	 */
	GWPlastic( double MU, double SIGMA,
	double MODULUS, double YIELD, double LENGTHSCALE, double ASPERITYAREA,
	double ADHESION_ENERGY, double ADHESION_MODULUS); 
	
	/*
	 * Destructor
	 */
	~GWPlastic();

	/*
	 * Reset parameters
	 */
	void ResetParameters(double DMIN);

	/*
	 * I/O
	 */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/*
	 * Returning values
	 */
	virtual double Function(double x) const;
	virtual double DFunction(double x) const;
	virtual double DDFunction(double x) const;

	/*
	 * Returning values in groups - derived classes should define
	 * their own non-virtual function called within this functon
	 * which maps in to out w/o requiring a virtual function call
	 * everytime. Default behavior is just to map the virtual functions
	 * above.
	 */
	virtual dArrayT& MapFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDFunction(const dArrayT& in, dArrayT& out) const;
	virtual dArrayT& MapDDFunction(const dArrayT& in, dArrayT& out) const;

	double PlasticArea(double dmin) const;
	double DPlasticArea(double dmin) const;
private:

	/* moment function */
	C1FunctionT *fmoment0, *fmoment1;

	/* parameters */
	double fE; // elastic modulus
	double fY; // yield value
	double fL; // lengthscale
	double fa0; // asperity area
	double fM; // mean
	double fS; // standard deviation
	double fW; // adhesion energy
	double fK; // adhesion modulus

	double fdmin; // mininum approach
	double fdc; // critical approach plastic
	double fda; // critical approach adhesion
	double fAe; // elastic area
	double fAp; // plastic area
};

} // namespace Tahoe 
#endif /* _GW_PLASTIC_H_ */
