/* $Id: J2PrimitiveT.cpp,v 1.3.50.1 2004-06-09 23:17:56 paklein Exp $ */
/* created: paklein (02/17/1997)                                          */
/* Base class for a J2 plastic material with linear kinematic/            */
/* isotropic hardening laws defined by:                                   */
/* 		H(a) = (1 - ftheta) fH_bar a                                         */
/* K(a) = fYield + ftheta fH_bar a                                        */
/* 		where a is the internal hardening variable                           */

#include "J2PrimitiveT.h"

#include <iostream.h>
#include <math.h>

#include "fstreamT.h"
#include "dSymMatrixT.h"


using namespace Tahoe;

const double sqrt23 = sqrt(2.0/3.0);

/* constructor */
J2PrimitiveT::J2PrimitiveT(ifstreamT& in):
	fYield(0.0),
	ftheta(-1.0),
	fH_bar(-1.0)
{
	/* read parameters */
	in >> fYield;	if (fYield <= 0.0) throw ExceptionT::kBadInputValue;
	in >> fH_bar;	if (fH_bar < 0.0) throw ExceptionT::kBadInputValue;
	in >> ftheta;	if (ftheta < 0.0 || ftheta > 1.0) throw ExceptionT::kBadInputValue;
}

/* destructor */
J2PrimitiveT::~J2PrimitiveT(void) { }

/***********************************************************************
 * Protected
 ***********************************************************************/

/* returns the value value of the yield function given the
* relative stress vector and state variables, where  alpha
* represents isotropic hardening.  NOTE: the relative stress
* should already contain the correction for any kinematic
* hardening. */
double J2PrimitiveT::YieldCondition(const dSymMatrixT& relstress,
	double alpha) const
{
	return sqrt(relstress.ScalarProduct()) - sqrt23*K(alpha);
}
