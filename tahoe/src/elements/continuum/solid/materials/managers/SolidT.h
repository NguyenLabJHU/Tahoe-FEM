/* $Id: SolidT.h,v 1.38.4.2 2005-02-24 01:10:51 thao Exp $ */
/* created: paklein (03/10/2001) */
#ifndef _MATERIAL_T_H_
#define _MATERIAL_T_H_

/* forward declarations */
#include "ios_fwd_decl.h"

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
       kFDSVKStV = 63,
       kSSSVKStV = 64,
       kOgdenMat = 65,
    kSSJ2LinHard = 67,
  kSSJ2LinHardplane = 66,
kSSLinearVEShear = 68,
  kLocJ2SSNlHard = 70,
  kGrdJ2SSNlHard = 71,
       kGradJ2SS = 72,
     kJ2SSKStV1D = 73,
     kABAQUS_BCJ = 80, /**< explicit integration */
 kABAQUS_BCJ_ISO = 81, /**< implicit integration */
kABAQUS_VUMAT_BCJ = 90,
kSIERRA_Hypoelastic = 110,
 kSIERRA_Iso_Geomat = 111,
 kWLC = 120
};

	/* stream extraction operator */ 
	friend istream& operator>>(istream& in, SolidT::TypeT& code);

/* 2D types */
#if 0
const int kSWDC100         = 24;//improper CB material
const int kSWDC110         = 15;	//improper CB material
const int kD2VIB           = 23; // plane stress VIB + gradient terms
#endif

/* 3D types */
#if 0
const int kIsoVIB_X	    = 14; // remove
#endif
  
};

} // namespace Tahoe 
#endif // _MATERIAL_T_H_
