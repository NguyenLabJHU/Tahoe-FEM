/* $Id: SolidT.cpp,v 1.37 2004-07-15 08:28:29 paklein Exp $ */
/* created: paklein (03/10/2001) */
#include "SolidT.h"
#include "ExceptionT.h"

using namespace Tahoe;

/* stream extraction operator */ 
SolidT::TypeT SolidT::int2TypeT(int i)
{
	switch (i)
	{
		case SolidT::kSSKStV:
			return SolidT::kSSKStV;

		case SolidT::kFDKStV:
			return SolidT::kFDKStV;

		case SolidT::kSSCubic:
			return SolidT::kSSCubic;

		case SolidT::kFDCubic:
			return SolidT::kFDCubic;

		case SolidT::kSimoIso:
			return SolidT::kSimoIso;

		case SolidT::kQuadLog:
			return SolidT::kQuadLog;

		case SolidT::kQuadLogOgden:
			return SolidT::kQuadLogOgden;

		case SolidT::kJ2SSKStV:
			return SolidT::kJ2SSKStV;

		case SolidT::kJ2SSKStV1D:
			return SolidT::kJ2SSKStV1D;

		case SolidT::kJ2Simo:
			return SolidT::kJ2Simo;

		case SolidT::kJ2QL:
			return SolidT::kJ2QL;

		case SolidT::kDPSSKStV:
			return SolidT::kDPSSKStV;

		case SolidT::kMRSSKStV:
			return SolidT::kMRSSKStV;

		case SolidT::kLJTr2D:
			return SolidT::kLJTr2D;

		case SolidT::kHex2D:
			return SolidT::kHex2D;

		case SolidT::kmodCauchyBornDC:
			return SolidT::kmodCauchyBornDC;

		case SolidT::kVIB:
			return SolidT::kVIB;

		case SolidT::kIsoVIBSimo:
			return SolidT::kIsoVIBSimo;

		case SolidT::kIsoVIBOgden:
			return SolidT::kIsoVIBOgden;

		case SolidT::kIsoVIBSimoJ2:
			return SolidT::kIsoVIBSimoJ2;

		case SolidT::kSSLinearVE:
			return SolidT::kSSLinearVE;

		case SolidT::kRGSplitVE:
			return SolidT::kRGSplitVE;

		case SolidT::kChain1D:
			return SolidT::kChain1D;

		case SolidT::kFossumSSIso:
			return SolidT::kFossumSSIso;

		case SolidT::kFCC:
			return SolidT::kFCC;

		case SolidT::kThermoViscoPlastic:
			return SolidT::kThermoViscoPlastic;

		case SolidT::kPovirk2D:
			return SolidT::kPovirk2D;

		case SolidT::kHyperEVP:
			return SolidT::kHyperEVP;

		case SolidT::kBCJHypo:
			return SolidT::kBCJHypo;

		case SolidT::kBCJHypoIsoDmgKE:
			return SolidT::kBCJHypoIsoDmgKE;

		case SolidT::kBCJHypoIsoDmgYC:
			return SolidT::kBCJHypoIsoDmgYC;

		case SolidT::kFDXtalElast:
			return SolidT::kFDXtalElast;

		case SolidT::kLocXtalPlast:
			return SolidT::kLocXtalPlast;

		case SolidT::kLocXtalPlast_C:
			return SolidT::kLocXtalPlast_C;

		case SolidT::kGrdXtalPlast:
			return SolidT::kGrdXtalPlast;

		case SolidT::kLocXtalPlastFp:
			return SolidT::kLocXtalPlastFp;

		case SolidT::kLocXtalPlastFp_C:
			return SolidT::kLocXtalPlastFp_C;

		case SolidT::kGrdXtalPlastFp:
			return SolidT::kGrdXtalPlastFp;

		case SolidT::kRGVIB:
			return SolidT::kRGVIB;

		case SolidT::kRGSplit:
			return SolidT::kRGSplit;

		case SolidT::kFDSVKStV:
			return SolidT::kFDSVKStV;

		case SolidT::kSSSVKStV:
			return SolidT::kSSSVKStV;

		case SolidT::kOgdenMat:
			return SolidT::kOgdenMat;

		case SolidT::kSSJ2LinHard:
			return SolidT::kSSJ2LinHard;

		case SolidT::kSSJ2LinHardplane:
			return SolidT::kSSJ2LinHardplane;

		case SolidT::kLocJ2SSNlHard:
			return SolidT::kLocJ2SSNlHard;

		case SolidT::kGrdJ2SSNlHard:
			return SolidT::kGrdJ2SSNlHard;

		case SolidT::kGradJ2SS:
			return SolidT::kGradJ2SS;

		case SolidT::kSIERRA_Hypoelastic:
			return SolidT::kSIERRA_Hypoelastic;

		case SolidT::kSIERRA_Iso_Geomat:
			return SolidT::kSIERRA_Iso_Geomat;

		case SolidT::kABAQUS_BCJ:
			return SolidT::kABAQUS_BCJ;

		case SolidT::kABAQUS_BCJ_ISO:
			return SolidT::kABAQUS_BCJ_ISO;

		case SolidT::kABAQUS_VUMAT_BCJ:
			return SolidT::kABAQUS_VUMAT_BCJ;

		case SolidT::kFCCEAM:
			return SolidT::kFCCEAM;

		default:
			ExceptionT::BadInputValue("SolidT::int2TypeT", "unknown code %d", i);
	}
	
	/* dummy */
	return kSSKStV;
}
