/* $Id: FCCLatticeT.h,v 1.2 2003-03-31 23:14:38 paklein Exp $ */
#ifndef _FCC_LATTICE_T_H_
#define _FCC_LATTICE_T_H_

/* base class */
#include "CBLatticeT.h"

namespace Tahoe {

/** a 3D hexagonal lattice. The number of shells of neighbors can be
 * selected up to 5th nearest. The bonds are scaled such that the
 * nearest neighbor bond distance is 1. */
class FCCLatticeT: public CBLatticeT
{
public:

	/** the Q matrix passed into this constructor is used to rotate the
	 * bond vectors into the orientation prescribed by Q.  No check is
	 * performed on the orthogonality of Q, only its dimensions.  Q is
	 * deep copied.  Q is defined as:
	 *
	 *			Q = d x_natural / d x_global
	 *
	 * So that the vectors are transformed by:
	 *
	 *			r_global = Transpose[Q].r_natural
	 */
	FCCLatticeT(const dMatrixT& Q, int nshells);

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

#endif /* _FCC_LATTICE_T_H_ */
