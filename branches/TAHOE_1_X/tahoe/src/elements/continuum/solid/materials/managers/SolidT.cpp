/* $Id: SolidT.cpp,v 1.36.4.1 2004-08-26 04:23:01 thao Exp $ */
/* created: paklein (03/10/2001) */
#include "SolidT.h"

#include <iostream.h>
#include "ExceptionT.h"

namespace Tahoe {

/* stream extraction operator */ 
istream& operator>>(istream& in, SolidT::TypeT& code)
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
		case SolidT::kJ2SSKStV1D:
			code = SolidT::kJ2SSKStV1D;
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
		case SolidT::kMRSSKStV:
			code = SolidT::kMRSSKStV;
			break;
		case SolidT::kLJTr2D:
			code = SolidT::kLJTr2D;
			break;
		case SolidT::kHex2D:
			code = SolidT::kHex2D;
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
		case SolidT::kSSLinearVE:
			code = SolidT::kSSLinearVE;
			break;
		case SolidT::kSSLinearVEShear:
			code = SolidT::kSSLinearVEShear;
			break;
		case SolidT::kRGSplitVE:
			code = SolidT::kRGSplitVE;
			break;
		case SolidT::kChain1D:
			code = SolidT::kChain1D;
			break;
		case SolidT::kFossumSSIso:
			code = SolidT::kFossumSSIso;
			break;
		case SolidT::kFCC:
			code = SolidT::kFCC;
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
		case SolidT::kRGVIB:
			code = SolidT::kRGVIB;
			break;
		case SolidT::kRGSplit:
			code = SolidT::kRGSplit;
			break;
		case SolidT::kFDSVKStV:
			code = SolidT::kFDSVKStV;
			break;
		case SolidT::kSSSVKStV:
			code = SolidT::kSSSVKStV;
			break;
		case SolidT::kOgdenMat:
			code = SolidT::kOgdenMat;
			break;
		case SolidT::kSSJ2LinHard:
			code = SolidT::kSSJ2LinHard;
			break;
		case SolidT::kSSJ2LinHardplane:
			code = SolidT::kSSJ2LinHardplane;
			break;
		case SolidT::kLocJ2SSNlHard:
			code = SolidT::kLocJ2SSNlHard;
			break;
		case SolidT::kGrdJ2SSNlHard:
			code = SolidT::kGrdJ2SSNlHard;
			break;
		case SolidT::kGradJ2SS:
			code = SolidT::kGradJ2SS;
			break;
		case SolidT::kSIERRA_Hypoelastic:
			code = SolidT::kSIERRA_Hypoelastic;
			break;
		case SolidT::kSIERRA_Iso_Geomat:
			code = SolidT::kSIERRA_Iso_Geomat;
			break;
		case SolidT::kABAQUS_BCJ:
			code = SolidT::kABAQUS_BCJ;
			break;
		case SolidT::kABAQUS_BCJ_ISO:
			code = SolidT::kABAQUS_BCJ_ISO;
			break;
		case SolidT::kABAQUS_VUMAT_BCJ:
			code = SolidT::kABAQUS_VUMAT_BCJ;
			break;
		case SolidT::kFCCEAM:
			code = SolidT::kFCCEAM;
			break;
		default:
			ExceptionT::BadInputValue("operator>>SolidT::TypeT", "unknown code %d", i_code);
	}
	return in;
}

}
