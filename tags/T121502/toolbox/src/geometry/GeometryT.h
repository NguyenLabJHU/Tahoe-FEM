/* $Id: GeometryT.h,v 1.3 2002-10-20 22:32:08 paklein Exp $ */
/* created: paklein (10/10/1999) */

#ifndef _GEOMETRY_T_H_
#define _GEOMETRY_T_H_

#include "Environment.h"
#include "ios_fwd_decl.h"
#include "ExceptionT.h"

namespace Tahoe {

/* forward declarations */
class GeometryBaseT;

/** class to define enumerations for element geometries and
 * associated operations */
class GeometryT
{
public:

	/** geometry types */
	enum CodeT {kNone          =-2,
	            kPoint         =-1,
	            kLine          = 0,
	            kQuadrilateral = 1,
	            kTriangle      = 2,
	            kHexahedron    = 3,
	            kTetrahedron   = 4,
	            kPentahedron   = 5  /**< not implemented, only cohesive surface output */};

	/** geometry_code -> nsd macro: of the parent domain */
	static int GeometryToNumSD(GeometryT::CodeT code);
	
	/** geometry names */
	static const char* fNames[8];
	
	/** convert GeometryT::CodeT to a string */
	static const char* ToString(GeometryT::CodeT code);
	
	/** return a pointer to a new GeometryBaseT. User is responsible for deleting class. */
	static GeometryBaseT* NewGeometry(GeometryT::CodeT geometry, int nen);
};

/** stream extraction operator for GeometryT::CodeT */
istream& operator>>(istream& in, GeometryT::CodeT& code);

/* convert GeometryT::CodeT to a string */
inline const char* GeometryT::ToString(GeometryT::CodeT code)
{
#if __option(extended_errorcheck)
	/* range check */
	if (code <  kNone || code > kPentahedron) throw ExceptionT::kOutOfRange;
#endif
	return fNames[code + 2];
}

} // namespace Tahoe 
#endif // _GEOMETRY_T_H_