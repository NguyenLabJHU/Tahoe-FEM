/* $Id: BasisT.h,v 1.5 2004-10-30 20:54:28 raregue Exp $ */
/* created: paklein (12/10/1999)                                          */
/* base class for basis functions                                         */

#ifndef _BASIS_T_H_
#define _BASIS_T_H_

/* direct members */
#include "dMatrixT.h"
#include "dArray2DT.h"
#include "nArray2DGroupT.h"

namespace Tahoe {

/* forward declarations */
class dArrayT;

class BasisT
{
public:

	/** constructor */
	BasisT(int complete, int nsd);

	/** destructor */
	virtual ~BasisT(void) { };

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
							// nstr = number of stress/strain dimensions (in vector form)
							// nsd = number of spatial dimensions 

	/* dynamic workspace manager */
	nArray2DGroupT<double> fArray2DGroup1; // [nbasis] x [nnd]
};

/* return values */
inline const dArray2DT& BasisT::P(void) const
{
	return fP;
}

inline const dArray2DT& BasisT::DP(int component) const
{
	return fDP[component];
}

inline const dArray2DT& BasisT::DDP(int component) const
{
	return fDDP[component];
}

inline const dArray2DT& BasisT::DDDP(int component) const // kyonten
{
	return fDDDP[component];
}

} // namespace Tahoe 
#endif /* _BASIS_T_H_ */
