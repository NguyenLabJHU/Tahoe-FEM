/* $Id: C1GeometryT.cpp,v 1.1 2003-09-29 19:58:59 rdorgan Exp $ */
#include "C1GeometryT.h"

#include <iostream.h>
#include "ExceptionT.h"
#include "ArrayT.h"

/* geometries */
#include "C1LineT.h"

using namespace Tahoe;

namespace Tahoe {
const bool ArrayT<C1GeometryT::CodeT>::fByteCopy = true;
} /* namespace Tahoe */

namespace Tahoe { 
/* initialize static geometry names array */
const char* C1GeometryT::fNames[3] = 
          {"None",
          "Point",
           "C1Line"};
} /* namespace Tahoe */ 

namespace Tahoe {
istream& operator>>(istream& in, C1GeometryT::CodeT& code)
{
	int i_code;
	in >> i_code;
	switch (i_code)
	{
		case C1GeometryT::kNone:
			code = C1GeometryT::kNone;
			break;
		case C1GeometryT::kPoint:
			code = C1GeometryT::kPoint;
			break;
		case C1GeometryT::kC1Line:
			code = C1GeometryT::kC1Line;
			break;
		default:
			cout << "\n operator>>C1GeometryT::CodeT: unknown code: "
			<< i_code<< endl;
			throw ExceptionT::kBadInputValue;	
	}
	return in;
}
} // namespace Tahoe 

/* geometry_code -> nsd macro */
namespace Tahoe {
int C1GeometryT::GeometryToNumSD(C1GeometryT::CodeT code)
{
	/* set spatial dimension */
	if (code == C1GeometryT::kC1Line)
		return 1;
	else
	  throw ExceptionT::kBadInputValue;

	return 0;
}
} /* namespace Tahoe */

/* return a pointer to a new C1GeometryBaseT. User is responsible for deleting class. */
namespace Tahoe {
C1GeometryBaseT* C1GeometryT::NewGeometry(C1GeometryT::CodeT geometry, int nen)
{
	C1GeometryBaseT* geom = NULL;
	try {
	switch (geometry)
	{
		case C1GeometryT::kC1Line:		
			geom = new C1LineT(nen);
			break;
	
		default:
			cout << "\n C1GeometryT::NewGeometry: unknown geometry code: " << geometry << endl;
			throw ExceptionT::kGeneralFail;			
	}
	}
	catch (ExceptionT::CodeT exception) {
		cout << "\n C1GeometryT::NewGeometry: caught exception: " <<  ExceptionT::ToString(exception) << endl;
		throw exception;
	}
	return geom;
}
} /* namespace Tahoe */
