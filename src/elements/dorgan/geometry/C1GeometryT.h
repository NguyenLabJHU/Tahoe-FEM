/* $Id: C1GeometryT.h,v 1.1 2003-09-29 19:58:59 rdorgan Exp $ */
#ifndef _C1_GEOMETRY_T_H_
#define _C1_GEOMETRY_T_H_

#include "Environment.h"
#include "ios_fwd_decl.h"
#include "ExceptionT.h"

namespace Tahoe {

/* forward declarations */
class C1GeometryBaseT;

/** class to define enumerations for element geometries and
 * associated operations */
class C1GeometryT
{
public:

	/** geometry types */
	enum CodeT {kNone          =-2,
	            kPoint         =-1,
	            kC1Line        = 0};

	/** geometry_code -> nsd macro: of the parent domain */
	static int GeometryToNumSD(C1GeometryT::CodeT code);
	
	/** geometry names */
	static const char* fNames[3];
	
	/** convert C1GeometryT::CodeT to a string */
	static const char* ToString(C1GeometryT::CodeT code);
	
	/** return a pointer to a new C1GeometryBaseT. User is responsible for deleting class. */
	static C1GeometryBaseT* NewGeometry(C1GeometryT::CodeT geometry, int nen);
};

/** stream extraction operator for C1GeometryT::CodeT */
istream& operator>>(istream& in, C1GeometryT::CodeT& code);

/* convert C1GeometryT::CodeT to a string */
inline const char* C1GeometryT::ToString(C1GeometryT::CodeT code)
{
#if __option(extended_errorcheck)
	/* range check */
	if (code <  kNone || code > kC1Line) throw ExceptionT::kOutOfRange;
#endif
	return fNames[code + 2];
}

} // namespace Tahoe 
#endif // _C1_GEOMETRY_T_H_
