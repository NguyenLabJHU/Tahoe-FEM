/* $Id: PairPropertyT.cpp,v 1.2 2002-12-04 18:55:30 paklein Exp $ */
#include "PairPropertyT.h"

using namespace Tahoe;

/* constructor */
PairPropertyT::PairPropertyT(void) {}

/* return Paradyn-style coefficients table */
bool PairPropertyT::getParadynTable(const double** coeff, double& dr, int& row_size, int& num_rows) const
{
	*coeff = NULL;
	dr = 1.0;
	row_size = 0;
	num_rows = 0;
	return false;
}
