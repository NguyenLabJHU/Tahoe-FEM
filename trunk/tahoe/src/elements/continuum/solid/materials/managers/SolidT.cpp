/* $Id: SolidT.cpp,v 1.1 2001-04-27 10:53:30 paklein Exp $ */
/* created: paklein (03/10/2001)                                          */

#include "SolidT.h"

#include <iostream.h>
#include "ExceptionCodes.h"

/* stream extraction operator */ 
istream& operator>>(istream& in, SolidT::SolidT& code)
{
	int i_code;
	in >> i_code;
	switch (i_code)
	{
		case SolidT::kSSKStV:
			code = SolidT::kSSKStV;
			break;
		case SolidT::kFDKStV:
			code = SolidT::kFDKStV;
			break;
		case SolidT::kSSCubic:
			code = SolidT::kSSCubic;
			break;
		case SolidT::kFDCubic:
			code = SolidT::kFDCubic;
			break;
		case SolidT::kSimoIso:
			code = SolidT::kSimoIso;
			break;
		case SolidT::kQuadLog:
			code = SolidT::kQuadLog;
			break;
		case SolidT::kQuadLogOgden:
			code = SolidT::kQuadLogOgden;
			break;
		case SolidT::kJ2SSKStV:
			code = SolidT::kJ2SSKStV;
			break;
		case SolidT::kJ2Simo:
			code = SolidT::kJ2Simo;
			break;
		case SolidT::kJ2QL:
			code = SolidT::kJ2QL;
			break;
		case SolidT::kDPSSKStV:
			code = SolidT::kDPSSKStV;
			break;
		case SolidT::kLJTr2D:
			code = SolidT::kLJTr2D;
			break;
		case SolidT::kLJFCC111:
			code = SolidT::kLJFCC111;
			break;
		case SolidT::kmodCauchyBornDC:
			code = SolidT::kmodCauchyBornDC;
			break;
		case SolidT::kVIB:
			code = SolidT::kVIB;
			break;
		case SolidT::kIsoVIBSimo:
			code = SolidT::kIsoVIBSimo;
			break;
		case SolidT::kIsoVIBOgden:
			code = SolidT::kIsoVIBOgden;
			break;
		case SolidT::kIsoVIBSimoJ2:
			code = SolidT::kIsoVIBSimoJ2;
			break;
		case SolidT::kABAQUS_BCJ:
			code = SolidT::kABAQUS_BCJ;
			break;
		default:
			cout << "\n operator>>SolidT::SolidT: unknown code: "
			<< i_code<< endl;
			throw eBadInputValue;	
	}
	return in;
}
