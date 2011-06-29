/* $Id: FCCLatticeT.h,v 1.3 2004-07-15 08:26:42 paklein Exp $ */
#ifndef _FCC_LATTICE_T_H_
#define _FCC_LATTICE_T_H_

/* base classes */
#include "CBLatticeT.h"
#include "ParameterInterfaceT.h"

namespace Tahoe {

/** a 3D hexagonal lattice. The number of shells of neighbors can be
 * selected up to 5th nearest. The bonds are scaled such that the
 * nearest neighbor bond distance is 1. 
 * \note the lattice can be initialized either by
 * -# FCCLatticeT::Initialize
 * -# FCCLatticeT::TakeParameterList */
class FCCLatticeT: public CBLatticeT, public ParameterInterfaceT
{
public:

	/** orientation codes - for crystal axes rotated wrt global axes*/
	enum OrientationCodeT {
	kFCC3Dnatural = 0,
	    kFCC3D110 = 1,
	  kFCC3D111_a = 2, /**< x,y,z = <-1 1 0> <-1-1 2> < 1 1 1> */
	  kFCC3D111_b = 3, /**< x,y,z = < 1-1 0> < 1 1-2> < 1 1 1> */
	  kEulerAngles = 4};

	/** constructor */
	FCCLatticeT(int nshells);

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
	 *			Q = d x_natural / d x_global
	 *
	 * So that the vectors are transformed by:
	 *
	 *			r_global = Transpose[Q].r_natural
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

#endif /* _FCC_LATTICE_T_H_ */
