/* $Id: OrthoMLS3DT.h,v 1.2 2002-07-02 19:56:56 cjkimme Exp $ */
/* created: paklein (07/03/1998)                                          */

#ifndef _ORTHO_MLS_3D_T_H_
#define _ORTHO_MLS_3D_T_H_

/* base class */
#include "OrthoMLSSolverT.h"


namespace Tahoe {

class OrthoMLS3DT: public OrthoMLSSolverT
{
public:

	/* constructor */
	OrthoMLS3DT(int complete);
	
protected:

	/* return the number of monomial terms for the given completeness */
	virtual int NumberOfMonomials(int completeness) const;

	/* evaluate monomials and derivatives at coords */
	virtual void SetMonomials(const dArrayT& coords, dArrayT& p, dArray2DT& Dp);

};

} // namespace Tahoe 
#endif /* _ORTHO_MLS_3D_T_H_ */
