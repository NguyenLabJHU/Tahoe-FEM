/* $Id: HookeanMatT.cpp,v 1.4 2004-07-15 08:26:56 paklein Exp $ */
/* created: paklein (06/09/1997) */
#include "HookeanMatT.h"
#include "dSymMatrixT.h"

using namespace Tahoe;

/* constructor */
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

/* dimension */
void HookeanMatT::Dimension(int nsd)
{
	fModulus.Dimension(dSymMatrixT::NumValues(nsd));
	fModulus =-1.0;
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
