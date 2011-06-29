/* $Id: TensorTransformT.h,v 1.6 2004-07-15 08:29:20 paklein Exp $ */
#ifndef _TENSOR_TRANSFORM_T_H_
#define _TENSOR_TRANSFORM_T_H_

/* direct members */
#include "dMatrixT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/** class to handle tensor transformations */
class TensorTransformT
{
public:

	/** constructors */
	TensorTransformT(int dim);
	TensorTransformT(void) { };
	
	/** set the dimension */
	void Dimension(int dim);

	/** rank 2 tensor transformations */
	const dSymMatrixT& PushForward(const dMatrixT& fwd, const dSymMatrixT& a);
	const dSymMatrixT& PullBack(const dMatrixT& fwd, const dSymMatrixT& a);

	/** rank 4 tensor transformations */
	const dMatrixT& PushForward(const dMatrixT& fwd, const dMatrixT& a);
	const dMatrixT& PullBack(const dMatrixT& fwd, const dMatrixT& a);
	
private:

	/** rank 4 tensor transformations */
	void FFFFC_2D(const dMatrixT& F, dMatrixT& C);
	void FFFFC_3D(const dMatrixT& F, dMatrixT& C);
	void FFFFC_2D_Z(const dMatrixT& F, dMatrixT& C) const;
			
private:

	/* return values */
	dSymMatrixT fRank2;
	dMatrixT    fRank4;
	
	/* work space */
	dMatrixT fPull;
	dMatrixT fTransform;
	dMatrixT fOuter;
	dMatrixT fRedMat;
};

/* inlines */
inline const dSymMatrixT& TensorTransformT::PushForward(const dMatrixT& fwd, 
	const dSymMatrixT& a)
{
	/* transformation */
	fRank2.MultQBQT(fwd, a);
	return fRank2;
}

inline const dSymMatrixT& TensorTransformT::PullBack(const dMatrixT& fwd, 
	const dSymMatrixT& a)
{
	/* inverse transform */
	fPull.Inverse(fwd);

	/* transformation */
	fRank2.MultQBQT(fPull, a);
	return fRank2;
}

} // namespace Tahoe 
#endif /* _TENSOR_TRANSFORM_T_H_ */
