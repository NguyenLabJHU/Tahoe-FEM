/* File: GeometryT.h */

#ifndef _GEOMETRY_T_H_
#define _GEOMETRY_T_H_

#include "Environment.h"
#include "ios_fwd_decl.h"

class GeometryT
{
  public:

	/* geometries supported by the derived classes */
	enum GeometryCode {kNone          =-2,
                       kPoint         =-1,
	                   kLine          = 0,
	                   kQuadrilateral = 1,
	                   kTriangle      = 2,
	                   kHexahedron    = 3,
	                   kTetrahedron   = 4,
	                   kPentahedron   = 5}; //not implemented, only CSE output
	friend istream& operator>>(istream& in, GeometryT::GeometryCode& code);
};

#endif // _GEOMETRY_T_H_
