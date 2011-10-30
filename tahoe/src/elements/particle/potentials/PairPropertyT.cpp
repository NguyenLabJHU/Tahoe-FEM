/* $Id: PairPropertyT.cpp,v 1.4.70.1 2011-10-30 06:26:11 bcyansfn Exp $ */
#include "PairPropertyT.h"
#include <cstddef>

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
