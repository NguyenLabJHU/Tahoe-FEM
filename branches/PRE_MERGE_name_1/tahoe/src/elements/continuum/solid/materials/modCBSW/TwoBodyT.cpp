/* $Id: TwoBodyT.cpp,v 1.1.1.1 2001-01-29 08:20:26 paklein Exp $ */
/* created: paklein (10/11/1997)                                          */

#include "TwoBodyT.h"

/* constructor */
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
