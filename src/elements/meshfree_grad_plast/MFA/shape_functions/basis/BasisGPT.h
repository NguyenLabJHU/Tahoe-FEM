/* $Id: BasisGPT.h,v 1.2 2004-07-14 19:51:08 kyonten Exp $ */
/* created: paklein (12/10/1999)                                          */
/* base class for basis functions                                         */

#ifndef _BASIS_GP_T_H_
#define _BASIS_GP_T_H_

/* direct members */
#include "dMatrixT.h"
#include "dArray2DT.h"
#include "nArray2DGroupT.h"

namespace Tahoe {

/* forward declarations */
class dArrayT;

class BasisGPT
{
public:

	/** constructor */
	BasisGPT(int complete, int nsd);

	/** destructor */
	virtual ~BasisGPT(void) { };

	/* return the number of basis functions */
	virtual int BasisDimension(void) const = 0;

	/* evaluate basis functions at coords */
	virtual void SetBasis(const dArray2DT& coords, int order) = 0;

	/* return values */
	const dArray2DT& P(void) const;
	const dArray2DT& DP(int component) const;
	const dArray2DT& DDP(int component) const;
	const dArray2DT& DDDP(int component) const; // kyonten

protected:

	/* dimension work space */
	virtual void Dimension(int num_nodes);

protected:

	int fComplete;
	int fNumSD;

	dArray2DT fP;           // [nbasis] x [nnd]
	ArrayT<dArray2DT> fDP;  // [nsd] x [nbasis] x [nnd]
	ArrayT<dArray2DT> fDDP; // [nstr] x [nbasis] x [nnd]
	ArrayT<dArray2DT> fDDDP; // [nstr] x [nbasis] x [nnd] // kyonten 
							// difference b/ween nstr and nsd?  

	/* dynamic workspace manager */
	nArray2DGroupT<double> fArray2DGroup1; // [nbasis] x [nnd]
};

/* return values */
inline const dArray2DT& BasisGPT::P(void) const
{
	return fP;
}

inline const dArray2DT& BasisGPT::DP(int component) const
{
	return fDP[component];
}

inline const dArray2DT& BasisGPT::DDP(int component) const
{
	return fDDP[component];
}

inline const dArray2DT& BasisGPT::DDDP(int component) const // kyonten
{
	return fDDDP[component];
}
} // namespace Tahoe 
#endif /* _BASIS_GP_T_H_ */
