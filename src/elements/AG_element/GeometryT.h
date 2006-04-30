/* $Id: GeometryT.h,v 1.1 2006-04-30 04:09:15 regueiro Exp $ */
 /* created: paklein (10/10/1999) */

 #ifndef _GEOMETRY_T_H_
 #define _GEOMETRY_T_H_

 #include "Environment.h"
 #include "ios_fwd_decl.h"
 #include "ExceptionT.h"
 
 #include "ParameterListT.h"

 namespace Tahoe {

 /* forward declarations */
 class GeometryBaseT;
 class ParameterInterfaceT;
 class ParameterListT;
 class StringT;

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
                 kPentahedron   = 5, /*< not implemented, only cohesive surface output */
                 kQuadrilateralAG	= 6 };   
	
     /** convert int to GeometryT::CodeT */
     static CodeT int2CodeT(int i);

     /** convert string to GeometryT::CodeT */
     static CodeT string2CodeT(const char* name);

     /** geometry_code -> nsd macro: of the parent domain */
     static int GeometryToNumSD(GeometryT::CodeT code);

     /** geometry names */
     static const char* fNames[9];

     /** convert GeometryT::CodeT to a string */
     static const char* ToString(GeometryT::CodeT code);

     /** return a pointer to a new GeometryBaseT. User is responsible for deleting class. */
     static GeometryBaseT* New(GeometryT::CodeT geometry, int nen);

     /** return a description of the given geometry name or NULL if the name
      * does not match any of the geometry names. The host code is responsible
      * for deleting the returned object. */
     static ParameterInterfaceT* New(const StringT& name);
      
    /** pointer to AG parameter list */
	//static ParameterListT* fQuadAG;
	//static ParameterListT& fQuadAG;
	     
 };

 /** stream extraction operator for GeometryT::CodeT */
 istream& operator>>(istream& in, GeometryT::CodeT& code);

 /* convert GeometryT::CodeT to a string */
 inline const char* GeometryT::ToString(GeometryT::CodeT code)
 {
 #if __option(extended_errorcheck)
     /* range check */
     if (code <  kNone || code > kQuadrilateralAG) throw ExceptionT::kOutOfRange;
 #endif
     return fNames[code + 2];
 }

 } // namespace Tahoe
 #endif // _GEOMETRY_T_H_

