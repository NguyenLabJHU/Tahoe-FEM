/* $Id: nIntegratorT.cpp,v 1.8 2004-12-26 21:08:50 d-farrell2 Exp $ */
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
void nIntegratorT::Corrector(BasicFieldT& field, const dArray2DT& update, int fieldstart /*= 0*/, int fieldend /*= -1*/, int dummy /*= 0*/)
{
#pragma unused(field)
#pragma unused(update)
#pragma message("nIntegratorT::Corrector: make me pure virtual")
ExceptionT::GeneralFail("nIntegratorT::Corrector", "not implemented");
}
