/* $Id: HexLattice2DT.h,v 1.5 2005-08-30 07:16:49 jzimmer Exp $ */
#ifndef _HEX_LATTICE_2D_T_H_
#define _HEX_LATTICE_2D_T_H_

/* base class */
#include "CBLatticeT.h"
#include "ParameterInterfaceT.h"

namespace Tahoe {

/** a 2D hexagonal lattice. The number of shells of neighbors can be
 * selected up to 5th nearest. The bonds are scaled such that the
 * nearest neighbor bond distance is 1.
 * \note the lattice can be initialized either by
 * -# HEXLattice2DT::Initialize
 * -# HEXLattice2DT::TakeParameterList */
class HexLattice2DT: public CBLatticeT, public ParameterInterfaceT
{
public:

        /** orientation codes - for crystal axes rotated wrt global axes*/
        enum OrientationCodeT {
        kHEX2Dnatural = 0,
        kHEX2D90 = 1,
        kRotationAngle = 2};

	/** constructor */
	HexLattice2DT(int nshells);

	/** number of shells */
	int NumShells(void) const { return fNumShells; };

        /** \name implementation of the ParameterInterfaceT interface */
        /*@{*/
        /** information about subordinate parameter lists */
        virtual void DefineSubs(SubListT& sub_list) const;

        /** a pointer to the ParameterInterfaceT of the given subordinate */
        virtual ParameterInterfaceT* NewSub(const StringT& name) const;

        /** accept parameter list */
        virtual void TakeParameterList(const ParameterListT& list);
        /*@}*/

        /** set the transformation matrix for the given orientation, where Q is
         * defined as:
         *
         *                      Q = d x_natural / d x_global
         *
         * So that the vectors are transformed by:
         *
         *                      r_global = Transpose[Q].r_natural
         */
        static void SetQ(const ParameterListT& list, dMatrixT& Q);

protected:

	/** initialize bond table values */
	virtual void LoadBondTable(void);

private:

	/** number of shells of neighbors */
	int fNumShells;
};

} /* namespace Tahoe */

#endif /* _HEX_LATTICE_2D_T_H_ */
