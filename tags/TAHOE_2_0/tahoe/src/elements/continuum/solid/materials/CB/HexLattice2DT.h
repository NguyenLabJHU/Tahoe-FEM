/* $Id: HexLattice2DT.h,v 1.3 2004-07-15 08:26:42 paklein Exp $ */
#ifndef _HEX_LATTICE_2D_T_H_
#define _HEX_LATTICE_2D_T_H_

/* base class */
#include "CBLatticeT.h"

namespace Tahoe {

/** a 2D hexagonal lattice. The number of shells of neighbors can be
 * selected up to 5th nearest. The bonds are scaled such that the
 * nearest neighbor bond distance is 1. */
class HexLattice2DT: public CBLatticeT
{
public:

	/** constructor */
	HexLattice2DT(int nshells);

	/** number of shells */
	int NumShells(void) const { return fNumShells; };

protected:

	/** initialize bond table values */
	virtual void LoadBondTable(void);

private:

	/** number of shells of neighbors */
	int fNumShells;
};

} /* namespace Tahoe */

#endif /* _HEX_LATTICE_2D_T_H_ */
