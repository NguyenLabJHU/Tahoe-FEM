/*
 * File: dArrayT.cpp
 */

/*
 * created      : PAK (08/11/96)
 * last modified: PAK (08/20/97)
 */

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "Constants.h"
#include "dArrayT.h"

/* constructor */ 
dArrayT::dArrayT(void) { }
dArrayT::dArrayT(int length): nArrayT<double>(length) { }
dArrayT::dArrayT(int length, double* p): nArrayT<double>(length,p) { }
dArrayT::dArrayT(const dArrayT& source): nArrayT<double>(source) { }

/* L2 norm of the vector */
double dArrayT::Magnitude(void) const
{
	register double magsqr = 0.0;
	double* p = Pointer();

	for (int i = 0; i < Length(); i++)
	{
		magsqr += (*p)*(*p); 
		p++;
	}
	
	return sqrt(magsqr);
}

/*
 * Returns the univector in the direction of *this.  If no argument
 * is passed in, *this is scaled and a referenc to *this is returned.
 */
dArrayT& dArrayT::UnitVector(const dArrayT& vector)
{
	SetToScaled(1.0/vector.Magnitude(), vector);
	return *this;
}
