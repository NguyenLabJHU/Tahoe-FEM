/* $Id: MeshFreeT.h,v 1.6 2002-07-05 22:28:37 paklein Exp $ */
/* created: paklein (12/08/1999)                                          */

#ifndef _MESHFREE_T_H_
#define _MESHFREE_T_H_

#include "ios_fwd_decl.h"


namespace Tahoe {

/** class to define enumerations for meshfree methods */
class MeshFreeT
{
public:

	/** mesh free formulation code */
	enum FormulationT {kEFG = 0, /**< Element Free Galerkin method */
	                  kRKPM = 1  /**< Reproducing Kernel Particle Method */};

	/** input extraction operator */
	friend istream& operator>>(istream& in, MeshFreeT::FormulationT& code);

	/** window types */
	enum WindowTypeT {kGaussian = 0, /**< Guassian with spherical support */
	               kCubicSpline = 1, /**< cubic spline with spherical support */
	                     kBrick = 2  /**< product of Gaussians with quad/hex support */ };	                 

	/** input extraction operator */
	friend istream& operator>>(istream& in, MeshFreeT::WindowTypeT& code);
};

} // namespace Tahoe 
#endif /* _MESHFREE_T_H_ */
