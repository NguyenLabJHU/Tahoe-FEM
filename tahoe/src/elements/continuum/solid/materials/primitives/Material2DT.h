/* $Id: Material2DT.h,v 1.1.1.1.10.1 2002-06-27 18:03:51 cjkimme Exp $ */
/* created: paklein (02/15/1997)                                          */
/* Mix-in for 2D materials.                                               */

#ifndef _MATERIAL_2D_T_H_
#define _MATERIAL_2D_T_H_

#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"

namespace Tahoe {

class ifstreamT;

class Material2DT
{
public:

	/* 2D constrain options */
	enum ConstraintOptionT {kPlaneStress = 1,
                            kPlaneStrain = 2};
	friend istream& operator>>(istream& in, Material2DT::ConstraintOptionT& option);

	/* constructor */
	Material2DT(ifstreamT& in);
	Material2DT(ifstreamT& in, ConstraintOptionT constraintopt);
	
	/* I/O functions */
	void Print(ostream& out) const;

	/* accessors */
	double Thickness(void) const;
	int ConstraintOption(void) const;
	
protected:

	double fThickness;
	ConstraintOptionT fConstraintOption;	
};

} // namespace Tahoe 
#endif /* _MATERIAL_2D_T_H_ */
