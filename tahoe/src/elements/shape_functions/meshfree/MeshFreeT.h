/* $Id: MeshFreeT.h,v 1.3 2001-07-03 01:35:50 paklein Exp $ */
/* created: paklein (12/08/1999)                                          */

#ifndef _MESHFREE_T_H_
#define _MESHFREE_T_H_

#include "ios_fwd_decl.h"

class MeshFreeT
{
public:

	/* mesh free formulation code */
	enum FormulationT {kEFG = 0,
	                  kRKPM = 1};
	friend istream& operator>>(istream& in, MeshFreeT::FormulationT& code);

	/* window types */
	enum WindowTypeT {kGaussian = 0,
	               kCubicSpline = 1,
	                     kBrick = 2};	                 
	friend istream& operator>>(istream& in, MeshFreeT::WindowTypeT& code);
};

#endif // _MESHFREE_T_H_
