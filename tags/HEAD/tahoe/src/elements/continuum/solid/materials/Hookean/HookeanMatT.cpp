/* $Id: HookeanMatT.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/09/1997)                                          */
/* Base class for all Hookean materials, defined as:                      */
/* 	stress_ij = moduli_ijkl strain_kl                                     */

#include "HookeanMatT.h"
#include "dSymMatrixT.h"

/* constructor */
HookeanMatT::HookeanMatT(void) { }

/***********************************************************************
* Protected
***********************************************************************/

/* symmetric (2nd PK) stress */
void HookeanMatT::HookeanStress(const dMatrixT& moduli,
	const dSymMatrixT& strain, dSymMatrixT& stress)									
{
	/* symmetric rank-4 - rank-2 contraction */
	stress.A_ijkl_B_kl(moduli, strain);
}								

/* strain energy density for the specified strain.
* defined by:
*
*	w = 1/2 e_ij c_ijkl e_kl
*/
double HookeanMatT::HookeanEnergy(const dMatrixT& moduli,
	const dSymMatrixT& strain)
{
	/* double contraction */
	return 0.5*strain.B_ij_A_ijkl_B_kl(moduli);
}
