/* $Id: Anisotropic2DT.h,v 1.4 2004-07-15 08:29:19 paklein Exp $ */
/* created: paklein (06/11/1997) */
#ifndef _ANISOTROPIC2D_T_H_
#define _ANISOTROPIC2D_T_H_

#include "Environment.h"

/* direct members */
#include "dSymMatrixT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class Rotate2DT;
class dMatrixT;

/** base class for 2D anisotropic materials */
class Anisotropic2DT
{
public:

	/** constructor */
	Anisotropic2DT(ifstreamT& in);
	Anisotropic2DT(void);

	/** set the rotation angle */
	void SetRotation(double angle);

	/* destructor */
	virtual ~Anisotropic2DT(void);
	
	/* I/O functions */
	void Print(ostream& out) const;

protected:

	/* transformation tensor */
	const dMatrixT& Q(void) const;

	/* return a reference to the transformed vector.  Note, returns a
	 * references to the argument if !fIsRotated */
	const dArrayT& TransformIn(const dArrayT& vector);	
	const dArrayT& TransformOut(const dArrayT& vector);	
		
	/* return a reference to the transformed reduced index matrix (stress or
	 * strain).  Note, returns a references to the argument if !fIsRotated */
	const dSymMatrixT& TransformIn(const dSymMatrixT& redmat);	
	const dSymMatrixT& TransformOut(const dSymMatrixT& redmat);	
	
	/* 4th rank tensor tranformation - use for cases where moduli are constant
	 * and could therefore be stored in their transformed state */
	void TransformIn(dMatrixT& redtensor);	
	void TransformOut(dMatrixT& redtensor);	
	
private:

	/** called by constructor */
	void Construct(istream& in);
		
protected:
	
	/** coordinate transformer */
	Rotate2DT* fRotator;	
};

} // namespace Tahoe 
#endif /* _ANISOTROPIC2D_T_H_ */
