/* $Id: HookeanMatT.cpp,v 1.5 2004-07-22 21:09:32 paklein Exp $ */
/* created: paklein (06/09/1997) */
#include "HookeanMatT.h"
#include "dSymMatrixT.h"
#include "ParameterUtils.h"

using namespace Tahoe;

/* constructor */
HookeanMatT::HookeanMatT(int nsd):
	ParameterInterfaceT("Hookean"),
	fModulus(dSymMatrixT::NumValues(nsd))
{
	/* must set initialized later */
	fModulus =-1.0;
}

HookeanMatT::HookeanMatT(void):
	ParameterInterfaceT("Hookean")
{

}

/* destructor */
HookeanMatT::~HookeanMatT(void)
{

}

/* dimension */
void HookeanMatT::Dimension(int nsd)
{
	fModulus.Dimension(dSymMatrixT::NumValues(nsd));
	fModulus =-1.0;
}

/* information about subordinate parameter lists */
void HookeanMatT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* option to define all components of the modulus explicitly */
	sub_list.AddSub("full_modulus_matrix_choice", ParameterListT::ZeroOrOnce, true);
}

/* return the description of the given inline subordinate parameter list */
void HookeanMatT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "full_modulus_matrix_choice")
	{
		order = ParameterListT::Choice;

		/* 1D/2D/3D modulus matricies */	
		sub_lists.AddSub("modulus_1x1");
		sub_lists.AddSub("modulus_3x3");
		sub_lists.AddSub("modulus_6x6");
	}
	else /* inherited */
		ParameterInterfaceT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* HookeanMatT::NewSub(const StringT& name) const
{
	if (name == "modulus_1x1" || name == "modulus_3x3" || name == "modulus_6x6")
		return new MatrixParameterT(name, 'C');
	else /* inherited */
		return ParameterInterfaceT::NewSub(name);
}

/* accept parameter list */
void HookeanMatT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* look for modulus definition */
	const ParameterListT* modulus = list.ListChoice(*this, "full_modulus_matrix_choice");
	if (modulus) /* extract the matrix values */
	{
		int old_dim = fModulus.Rows();
		MatrixParameterT::Extract(*modulus, fModulus, 'C');
		if (old_dim > 0 && old_dim != fModulus.Rows())
			ExceptionT::GeneralFail("HookeanMatT::TakeParameterList", "\"%s\" should be dimension %dx%d",
				modulus->Name().Pointer(), old_dim, old_dim);
	}
	else /* derived class must handle initialization of the modulus */
		SetModulus(fModulus);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

void HookeanMatT::SetModulus(dMatrixT& modulus) {
#pragma unused(modulus)
	ExceptionT::GeneralFail("HookeanMatT::SetModulus", "must be overridden or input must include modulus");
}

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
