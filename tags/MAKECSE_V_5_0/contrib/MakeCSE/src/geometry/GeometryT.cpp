/* File: GeometryT.cpp */

#include "GeometryT.h"

#include <iostream.h>
#include "ExceptionCodes.h"

istream& operator>>(istream& in, GeometryT::GeometryCode& code)
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
			cout << "\n operator>>GeometryT::GeometryCode: unknown code: " 
			<< i_code<< endl;
			throw eBadInputValue;	
	}
	return in;
}
