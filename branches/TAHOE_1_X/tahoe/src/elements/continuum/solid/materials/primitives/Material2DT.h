/* $Id: Material2DT.h,v 1.5 2003-04-12 22:10:48 thao Exp $ */
/* created: paklein (02/15/1997) */

#ifndef _MATERIAL_2D_T_H_
#define _MATERIAL_2D_T_H_

#include "Environment.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** mix-in for 2D materials */
class Material2DT
{
public:

	/** \name 2D constrain options */
	enum ConstraintOptionT {kPlaneStress = 1, /**< plane stress */
                            kPlaneStrain = 2  /**< plane strain */};

	/** stream extraction operators */
	friend istream& operator>>(istream& in, Material2DT::ConstraintOptionT& option);

	/** \name constructor */
	/*@{*/
	Material2DT(ifstreamT& in);
	Material2DT(ifstreamT& in, ConstraintOptionT constraintopt);
	Material2DT(void);
        /*@}*/
	
	/** I/O functions */
	void Print(ostream& out) const;

	/** \name accessors */
	/*@{*/
	double Thickness(void) const { return fThickness; };
	ConstraintOptionT ConstraintOption(void) const { return fConstraintOption; };
	/*@}*/
	
protected:

	double fThickness;
	ConstraintOptionT fConstraintOption;	
};

} // namespace Tahoe 
#endif /* _MATERIAL_2D_T_H_ */
