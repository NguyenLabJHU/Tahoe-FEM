/* $Id: Lattice1DT.h,v 1.1.4.1 2004-04-24 19:57:32 paklein Exp $ */
#ifndef _LATTICE_1D_T_H_
#define _LATTICE_1D_T_H_

/* base class */
#include "CBLatticeT.h"

namespace Tahoe {

/** a 1D lattice */
class Lattice1DT: public CBLatticeT
{
public:

	/** constructor */
	Lattice1DT(int nshells);

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
