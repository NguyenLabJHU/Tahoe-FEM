/* $Id: TwoBodyT.cpp,v 1.1.1.1.10.1 2002-06-27 18:03:30 cjkimme Exp $ */
/* created: paklein (10/11/1997)                                          */

#include "TwoBodyT.h"

/* constructor */

using namespace Tahoe;

TwoBodyT::TwoBodyT(const dArrayT& lengths, const ThermalDilatationT* thermal):
	fLengths(lengths),
	fPhi(fLengths.Length()),
	fdPhi(fLengths.Length()),
	fddPhi(fLengths.Length()),
	fThermal(thermal)
{

}

/* Accessors */
const dArrayT& TwoBodyT::Phi(void) const  { return(fPhi);   }
const dArrayT& TwoBodyT::dPhi(void) const { return(fdPhi);  }
const dArrayT& TwoBodyT::ddPhi(void) const{ return(fddPhi); }
