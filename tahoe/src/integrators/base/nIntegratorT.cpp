/* $Id: nIntegratorT.cpp,v 1.7 2003-05-20 10:27:06 paklein Exp $ */
/* created: paklein (10/14/1996) */
#include "nIntegratorT.h"
#include "ExceptionT.h"

using namespace Tahoe;

/* constructor */
nIntegratorT::nIntegratorT(void) { }

/* destructor */
nIntegratorT::~nIntegratorT(void) { }

/* register field with the integrator */
void nIntegratorT::Dimension(const BasicFieldT& field)
{
#pragma unused(field)
}

/* corrector. Maps ALL degrees of freedom forward. */
void nIntegratorT::Corrector(BasicFieldT& field, const dArray2DT& update)
{
#pragma unused(field)
#pragma unused(update)
#pragma message("nIntegratorT::Corrector: make me pure virtual")
ExceptionT::GeneralFail("nIntegratorT::Corrector", "not implemented");
}
