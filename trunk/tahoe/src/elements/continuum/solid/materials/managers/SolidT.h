/* $Id: SolidT.h,v 1.40 2004-07-20 23:21:30 rdorgan Exp $ */
/* created: paklein (03/10/2001) */
#ifndef _MATERIAL_T_H_
#define _MATERIAL_T_H_

namespace Tahoe {

class SolidT
{
public:

	/* solid material types */
	enum TypeT {
         kSSKStV = 1,			
         kFDKStV = 2,			
        kSSCubic = 3,			
        kFDCubic = 4,			
        kSimoIso = 5,
        kQuadLog = 6,
   kQuadLogOgden = 7,
       kJ2SSKStV = 8,
         kJ2Simo = 9,
           kJ2QL = 10,
       kDPSSKStV = 11,
         kLJTr2D = 12,
          kHex2D = 13,
         kFCCEAM = 14,
kmodCauchyBornDC = 15,
            kVIB = 16,
     kIsoVIBSimo = 17,
    kIsoVIBOgden = 18,
   kIsoVIBSimoJ2 = 19,
            kFCC = 20,
     kSSLinearVE = 21,
      kRGSplitVE = 22,
        kChain1D = 23,
    kFossumSSIso = 26,
       kMRSSKStV = 28,
kThermoViscoPlastic = 30,
       kPovirk2D = 31,
       kHyperEVP = 40,
        kBCJHypo = 45,
kBCJHypoIsoDmgKE = 46,
kBCJHypoIsoDmgYC = 47,
    kFDXtalElast = 49,
   kLocXtalPlast = 50,
 kLocXtalPlast_C = 51,
   kGrdXtalPlast = 52,
 kLocXtalPlastFp = 55,
kLocXtalPlastFp_C = 56,
 kGrdXtalPlastFp = 57,
          kRGVIB = 60,
	kRGSplit = 61,
       kFDSVKStV = 63,
       kSSSVKStV = 64,
       kOgdenMat = 65,
    kSSJ2LinHard = 67,
  kSSJ2LinHardplane = 66,
  kLocJ2SSNlHard = 70,
  kGrdJ2SSNlHard = 71,
   kGradJ2SSKStV = 72,
     kJ2SSKStV1D = 73,
     kABAQUS_BCJ = 80, /**< explicit integration */
 kABAQUS_BCJ_ISO = 81, /**< implicit integration */
kABAQUS_VUMAT_BCJ = 90,
kSIERRA_Hypoelastic = 110,
 kSIERRA_Iso_Geomat = 111
};

	/** convert integer to SolidT::TypeT */
	static TypeT int2TypeT(int i);
};

} /* namespace Tahoe */

#endif // _MATERIAL_T_H_
