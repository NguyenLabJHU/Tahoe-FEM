/* $Id: PairPropertyT.cpp,v 1.3.24.1 2003-11-04 19:47:18 bsun Exp $ */
#include "PairPropertyT.h"
#include <stddef.h>

using namespace Tahoe;

/* constructor */
PairPropertyT::PairPropertyT(void)
{
	SetName("pair_property");
}

/* return Paradyn-style coefficients table */
bool PairPropertyT::getParadynTable(const double** coeff, double& dr, int& row_size, int& num_rows) const
{
	*coeff = NULL;
	dr = 1.0;
	row_size = 0;
	num_rows = 0;
	return false;
}
