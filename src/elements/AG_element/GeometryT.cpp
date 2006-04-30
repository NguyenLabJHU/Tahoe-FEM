/* $Id: GeometryT.cpp,v 1.1 2006-04-30 04:09:15 regueiro Exp $ */
 /* created: paklein (10/10/1999) */
 #include "GeometryT.h"

 #include <iostream.h>
 #include "ExceptionT.h"
 #include "ArrayT.h"

 /* geometries */
 #include "LineT.h"
 #include "QuadT.h"
 #include "QuadAGT.h"
 #include "TriT.h"
 #include "HexahedronT.h"
 #include "TetrahedronT.h"
 #include "PentahedronT.h"
 #include "ParameterContainerT.h"

 #include <string.h>

 using namespace Tahoe;
 
 namespace Tahoe {
 DEFINE_TEMPLATE_STATIC const bool ArrayT<GeometryT::CodeT>::fByteCopy = true;
 } /* namespace Tahoe */

 namespace Tahoe {
 /* initialize static geometry names array */
 const char* GeometryT::fNames[9] =
           {"none",
           "point",
            "line",
   "quadrilateral",
        "triangle",
      "hexahedron",
     "tetrahedron",
     "pentahedron",
 "quadrilateralAG"};

 /* convert string to GeometryT::CodeT */
 GeometryT::CodeT GeometryT::string2CodeT(const char* name)
 {
     CodeT codes[9] = {kNone, kPoint, kLine, kQuadrilateral, kTriangle,
         kHexahedron, kTetrahedron, kPentahedron, kQuadrilateralAG};
     for (int i = 0; i < 9; i++)
         if (strcmp(name, fNames[i]) == 0)
             return codes[i];

     /* fail */
     ExceptionT::GeneralFail("GeometryT::string2CodeT", "known geometry \"%s\"", name);
     return kNone;
 }

 } /* namespace Tahoe */

 namespace Tahoe {

 /* convert int to GeometryT::CodeT */
 GeometryT::CodeT GeometryT::int2CodeT(int i) {
     switch (i)
     {
         case GeometryT::kNone:
             return GeometryT::kNone;
         case GeometryT::kPoint:
             return GeometryT::kPoint;
         case GeometryT::kLine:
             return GeometryT::kLine;
         case GeometryT::kQuadrilateral:
             return GeometryT::kQuadrilateral;
         case GeometryT::kTriangle:
             return GeometryT::kTriangle;
         case GeometryT::kHexahedron:
             return GeometryT::kHexahedron;
         case GeometryT::kTetrahedron:
             return GeometryT::kTetrahedron;
         case GeometryT::kPentahedron:
             return GeometryT::kPentahedron;
         case GeometryT::kQuadrilateralAG:
             return GeometryT::kQuadrilateralAG;
         default:
             ExceptionT::GeneralFail("GeometryT::int2CodeT", "unknown code %d", i);
     }
     return GeometryT::kNone;
 }

 istream& operator>>(istream& in, GeometryT::CodeT& code)
 {
     int i_code;
     in >> i_code;
     code = GeometryT::int2CodeT(i_code);

     return in;
 }
 } /* namespace Tahoe */

 /* geometry_code -> nsd macro */
 namespace Tahoe {
 int GeometryT::GeometryToNumSD(GeometryT::CodeT code)
 {
     /* set spatial dimension */
     if (code == GeometryT::kPoint)
         return 0;
     else if (code == GeometryT::kLine)
         return 1;
     else if (code == GeometryT::kQuadrilateral ||
              code == GeometryT::kTriangle)
         return 2 ;
	else if (code == GeometryT::kQuadrilateralAG)
         return 2 ;
     else if (code == GeometryT::kHexahedron ||
              code == GeometryT::kTetrahedron)
         return 3 ;
     else
       throw ExceptionT::kBadInputValue;

     return 0;
 }
 } /* namespace Tahoe */

 /* return a pointer to a new GeometryBaseT. User is responsible for deleting class. */
 namespace Tahoe {
 GeometryBaseT* GeometryT::New(GeometryT::CodeT geometry, int nen)
 {
     const char caller[] = "GeometryT::New";
     GeometryBaseT* geom = NULL;
     try {
     switch (geometry)
     {
         case GeometryT::kLine:
             geom = new LineT(nen);
             break;

         case GeometryT::kQuadrilateral:
             geom = new QuadT(nen);
             break;

         case GeometryT::kTriangle:
             geom = new TriT(nen);
             break;

         case GeometryT::kHexahedron:
             geom = new HexahedronT(nen);
             break;

         case GeometryT::kTetrahedron:
             geom = new TetrahedronT(nen);
             break;

         case GeometryT::kPentahedron:
             geom = new PentahedronT(nen);
             break;

         case GeometryT::kQuadrilateralAG:
             geom = new QuadAGT(nen);
             break;
             
         default:
             ExceptionT::GeneralFail(caller, "unknown geometry code %d", geometry);
     }
     }
     catch (ExceptionT::CodeT exc) {
         ExceptionT::Throw(exc, caller);
     }
     return geom;
 }
 } /* namespace Tahoe */

 /* return a description of the given geometry name or NULL */
 namespace Tahoe {
 ParameterInterfaceT* GeometryT::New(const StringT& name)
 {
     if (name == GeometryT::ToString(GeometryT::kLine))
     {
         ParameterContainerT* line = new ParameterContainerT(name);

         /* integration rules */
         ParameterT num_ip(ParameterT::Integer, "num_ip");
         num_ip.AddLimit(1, LimitT::Only);
         num_ip.AddLimit(2, LimitT::Only);
         num_ip.AddLimit(3, LimitT::Only);
         num_ip.AddLimit(4, LimitT::Only);
         num_ip.SetDefault(2);
         line->AddParameter(num_ip);

         return line;
     }
     else if (name == GeometryT::ToString(GeometryT::kQuadrilateral))
     {
         ParameterContainerT* quad = new ParameterContainerT(name);

         /* integration rules */
         ParameterT num_ip(ParameterT::Integer, "num_ip");
         num_ip.AddLimit(1, LimitT::Only);
         num_ip.AddLimit(4, LimitT::Only);
         num_ip.AddLimit(5, LimitT::Only);
         num_ip.AddLimit(9, LimitT::Only);
         num_ip.AddLimit(16, LimitT::Only);
         num_ip.SetDefault(4);
         quad->AddParameter(num_ip);

         return quad;
     }
     else if (name == GeometryT::ToString(GeometryT::kQuadrilateralAG))
     {
         ParameterContainerT* quadAG = new ParameterContainerT(name);

         /* integration rules */
         ParameterT num_ip(ParameterT::Integer, "num_ip");
         num_ip.AddLimit(1, LimitT::Only);
         num_ip.AddLimit(4, LimitT::Only);
         num_ip.AddLimit(5, LimitT::Only);
         num_ip.AddLimit(9, LimitT::Only);
         num_ip.AddLimit(16, LimitT::Only);
         num_ip.SetDefault(9);
         quadAG->AddParameter(num_ip);
				
		// select m and n
		/*
		ParameterT mparam(ParameterT::Double, "m");
		mparam.SetDefault(1.0);
		fQuadAG.AddParameter(mparam);
		ParameterT nparam(ParameterT::Double, "n");
		nparam.SetDefault(1.0);
		fQuadAG.AddParameter(nparam);
		*/
		
		/* AG Type Parameter: AG-Edge or AG-Corner
		If AGType = 1 means AG-Edge and AGType = 0 means AG-Corner */
		/*
		ParameterT agtype(ParameterT::Integer, "AGType");
		agtype.SetDefault(1);
		fQuadAG.AddParameter(agtype);
		*/
		
		return quadAG;
     }
     else if (name == GeometryT::ToString(GeometryT::kTriangle))
     {
         ParameterContainerT* tri = new ParameterContainerT(name);

         /* integration rules */
         ParameterT num_ip(ParameterT::Integer, "num_ip");
         num_ip.AddLimit(1, LimitT::Only);
         num_ip.AddLimit(4, LimitT::Only);
         num_ip.AddLimit(6, LimitT::Only);
         num_ip.SetDefault(1);
         tri->AddParameter(num_ip);

         return tri;
     }
     else if (name == GeometryT::ToString(GeometryT::kHexahedron))
     {
         ParameterContainerT* hex = new ParameterContainerT(name);

         /* integration rules */
         ParameterT num_ip(ParameterT::Integer, "num_ip");
         num_ip.AddLimit(1, LimitT::Only);
         num_ip.AddLimit(8, LimitT::Only);
         num_ip.AddLimit(9, LimitT::Only);
         num_ip.AddLimit(27, LimitT::Only);
         num_ip.AddLimit(64, LimitT::Only);
         num_ip.SetDefault(8);
         hex->AddParameter(num_ip);

         return hex;
     }
     else if (name == GeometryT::ToString(GeometryT::kTetrahedron))
     {
         ParameterContainerT* tet = new ParameterContainerT(name);

         /* integration rules */
         ParameterT num_ip(ParameterT::Integer, "num_ip");
         num_ip.AddLimit(1, LimitT::Only);
         num_ip.AddLimit(4, LimitT::Only);
         num_ip.SetDefault(1);
         tet->AddParameter(num_ip);

         return tet;
     }
     else
         return NULL;
 }
 } /* namespace Tahoe */

