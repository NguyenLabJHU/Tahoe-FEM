/* $Id: Anisotropic2DT.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (06/11/1997)                                          */
/* Base class for 2D anisotropic materials                                */

#ifndef _ANISOTROPIC2D_T_H_
#define _ANISOTROPIC2D_T_H_

#include "Environment.h"

/* direct members */
#include "dSymMatrixT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class Rotate2DT;
class dMatrixT;

class Anisotropic2DT
{
public:

	/* constructor */
	Anisotropic2DT(ifstreamT& in);

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

	/* called by constructor */
	void Construct(istream& in);
		
protected:
			
	/* coordinate transformer */
	int		   fIsRotated;
	Rotate2DT* fRotator;
	
};

#endif /* _ANISOTROPIC2D_T_H_ */
