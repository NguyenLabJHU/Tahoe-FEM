/* $Id: HookeanMatT.cpp,v 1.3 2002-07-02 19:55:38 cjkimme Exp $ */
/* created: paklein (06/09/1997)                                          */
/* Base class for all Hookean materials, defined as:                      */
/* 	stress_ij = moduli_ijkl strain_kl                                     */

#include "HookeanMatT.h"
#include "dSymMatrixT.h"

/* constructor */

using namespace Tahoe;

HookeanMatT::HookeanMatT(int nsd):
	fModulus(dSymMatrixT::NumValues(nsd))
{
	/* must set initialized later */
	fModulus =-1.0;
}

/* destructor */
HookeanMatT::~HookeanMatT(void)
{

}

/* initialization */
void HookeanMatT::Initialize(void)
{
	SetModulus(fModulus);
}

/***********************************************************************
* Protected
***********************************************************************/

/* symmetric stress */
void HookeanMatT::HookeanStress(const dSymMatrixT& strain, 
	dSymMatrixT& stress) const									
{
	/* symmetric rank-4 - rank-2 contraction */
	stress.A_ijkl_B_kl(fModulus, strain);
}								

/* strain energy density for the specified strain.
 * defined by:
 *
 *	w = 1/2 e_ij c_ijkl e_kl
 */
double HookeanMatT::HookeanEnergy(const dSymMatrixT& strain) const
{
	/* double contraction */
	return 0.5*strain.B_ij_A_ijkl_B_kl(fModulus);
}
