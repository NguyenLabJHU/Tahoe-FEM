/* $Id: GeometryT.cpp,v 1.3 2002-10-20 22:32:08 paklein Exp $ */
/* created: paklein (10/10/1999) */
#include "GeometryT.h"

#include <iostream.h>
#include "ExceptionT.h"
#include "ArrayT.h"

/* geometries */
#include "LineT.h"
#include "QuadT.h"
#include "TriT.h"
#include "HexahedronT.h"
#include "TetrahedronT.h"
#include "PentahedronT.h"

using namespace Tahoe;

namespace Tahoe {
const bool ArrayT<GeometryT::CodeT>::fByteCopy = true;
} /* namespace Tahoe */

namespace Tahoe { 
/* initialize static geometry names array */
const char* GeometryT::fNames[8] = 
          {"None",
          "Point",
           "Line",
  "Quadrilateral",
       "Triangle",
     "Hexahedron",
    "Tetrahedron",
    "Pentahedron"};
} /* namespace Tahoe */ 

namespace Tahoe {
istream& operator>>(istream& in, GeometryT::CodeT& code)
{
	int i_code;
	in >> i_code;
	switch (i_code)
	{
		case GeometryT::kNone:
			code = GeometryT::kNone;
			break;
		case GeometryT::kPoint:
			code = GeometryT::kPoint;
			break;
		case GeometryT::kLine:
			code = GeometryT::kLine;
			break;
		case GeometryT::kQuadrilateral:
			code = GeometryT::kQuadrilateral;
			break;
		case GeometryT::kTriangle:
			code = GeometryT::kTriangle;
			break;
		case GeometryT::kHexahedron:
			code = GeometryT::kHexahedron;
			break;
		case GeometryT::kTetrahedron:
			code = GeometryT::kTetrahedron;
			break;
		case GeometryT::kPentahedron:
			code = GeometryT::kPentahedron;
			break;
		default:
			cout << "\n operator>>GeometryT::CodeT: unknown code: "
			<< i_code<< endl;
			throw ExceptionT::kBadInputValue;	
	}
	return in;
}
} // namespace Tahoe 

/* geometry_code -> nsd macro */
namespace Tahoe {
int GeometryT::GeometryToNumSD(GeometryT::CodeT code)
{
	/* set spatial dimension */
	if (code == GeometryT::kLine)
		return 1;
	else if (code == GeometryT::kQuadrilateral ||
	         code == GeometryT::kTriangle)
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
GeometryBaseT* GeometryT::NewGeometry(GeometryT::CodeT geometry, int nen)
{
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

		default:
			cout << "\n GeometryT::NewGeometry: unknown geometry code: " << geometry << endl;
			throw ExceptionT::kGeneralFail;			
	}
	}
	catch (ExceptionT::CodeT exception) {
		cout << "\n GeometryT::NewGeometry: caught exception: " <<  ExceptionT::ToString(exception) << endl;
		throw exception;
	}
	return geom;
}
} /* namespace Tahoe */
