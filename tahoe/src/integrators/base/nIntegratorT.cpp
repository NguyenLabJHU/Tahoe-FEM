/* $Id: nIntegratorT.cpp,v 1.6 2003-01-27 07:00:23 paklein Exp $ */
/* created: paklein (10/14/1996) */
#include "nIntegratorT.h"
#include "ExceptionT.h"

using namespace Tahoe;

/* constructor */
nIntegratorT::nIntegratorT(void) { }

/* destructor */
nIntegratorT::~nIntegratorT(void) { }

/* corrector. Maps ALL degrees of freedom forward. */
void nIntegratorT::Corrector(BasicFieldT& field, const dArray2DT& update)
{
#pragma unused(field)
#pragma unused(update)
#pragma message("nIntegratorT::Corrector: make me pure virtual")
ExceptionT::GeneralFail("nIntegratorT::Corrector", "not implemented");
}
