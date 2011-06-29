/* $Id: ThreeBodyT.cpp,v 1.2 2002-07-02 19:56:07 cjkimme Exp $ */
/* created: paklein (10/11/1997)                                          */
/* Base class for the 3 body contribution to the strain energy density    */

#include "ThreeBodyT.h"
#include "dMatrixT.h"
#include "iArray2DT.h"

/* parameters */

using namespace Tahoe;

const int kNumVars = 3; //number of arguments in Phi

/* constructor */
ThreeBodyT::ThreeBodyT(const dArrayT& lengths,
	const dArrayT& angles, const iArray2DT& bondpairs,
	const ThermalDilatationT* thermal):
	fLengths(lengths),
	fAngles(angles),
	fPairs(bondpairs),
	fPhi(fPairs.MajorDim()),
	fdPhi(fPairs.MajorDim(), kNumVars),
	fddPhi(fPairs.MajorDim(), kNumVars*kNumVars),
	fThermal(thermal)
{

}

/* Accessors */
const dArrayT& ThreeBodyT::Phi(void) const     { return(fPhi);   }
const dArray2DT& ThreeBodyT::dPhi(void) const  { return(fdPhi);  }
const dArray2DT& ThreeBodyT::ddPhi(void) const { return(fddPhi); }
