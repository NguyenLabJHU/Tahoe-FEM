/* $Id: SolidT.cpp,v 1.14 2002-06-20 01:18:58 thao Exp $ */
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
		case SolidT::kThermoViscoPlastic:
			code = SolidT::kThermoViscoPlastic;
			break;
		case SolidT::kPovirk2D:
			code = SolidT::kPovirk2D;
			break;
		case SolidT::kHyperEVP:
			code = SolidT::kHyperEVP;
			break;
		case SolidT::kBCJHypo:
			code = SolidT::kBCJHypo;
			break;
		case SolidT::kBCJHypoIsoDmgKE:
			code = SolidT::kBCJHypoIsoDmgKE;
			break;
		case SolidT::kBCJHypoIsoDmgYC:
			code = SolidT::kBCJHypoIsoDmgYC;
			break;
                case SolidT::kFDXtalElast:
                        code = SolidT::kFDXtalElast;
                        break;
		case SolidT::kLocXtalPlast:
			code = SolidT::kLocXtalPlast;
			break;
		case SolidT::kLocXtalPlast_C:
			code = SolidT::kLocXtalPlast_C;
			break;
		case SolidT::kGrdXtalPlast:
			code = SolidT::kGrdXtalPlast;
			break;
		case SolidT::kLocXtalPlastFp:
			code = SolidT::kLocXtalPlastFp;
			break;
		case SolidT::kLocXtalPlastFp_C:
			code = SolidT::kLocXtalPlastFp_C;
			break;
		case SolidT::kGrdXtalPlastFp:
			code = SolidT::kGrdXtalPlastFp;
			break;
		case SolidT::kOgdenViscVIB:
			code = SolidT::kOgdenViscVIB;
			break;
		case SolidT::kSSStandard:
			code = SolidT::kSSStandard;
			break;
		case SolidT::kFDStandard:
			code = SolidT::kFDStandard;
			break;
		case SolidT::kABAQUS_BCJ:
			code = SolidT::kABAQUS_BCJ;
			break;
		case SolidT::kABAQUS_VUMAT_BCJ:
			code = SolidT::kABAQUS_VUMAT_BCJ;
			break;
		case SolidT::kFCCEAM:
			code = SolidT::kFCCEAM;
			break;
		default:
			cout << "\n operator>>SolidT::SolidT: unknown code: "
			<< i_code<< endl;
			throw eBadInputValue;	
	}
	return in;
}
