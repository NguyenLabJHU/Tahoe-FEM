/* $Id: VIB.cpp,v 1.14 2004-07-15 08:27:40 paklein Exp $ */
/* created: paklein (10/30/1997) */
#include "VIB.h"

#include "dSymMatrixT.h"
#include "C1FunctionT.h"

using namespace Tahoe;

/* constructors */
VIB::VIB(int nsd, int numstress, int nummoduli):
	ParameterInterfaceT("VIB_base"),
	fNumSD(nsd),
	fPotential(NULL),
	fNumStress(numstress),
	fNumModuli(nummoduli)
{

}

/* destructor */
VIB::~VIB(void) { delete fPotential; }

/* information about subordinate parameter lists */
void VIB::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* choice of parameters */
	sub_list.AddSub("VIB_potential_choice", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void VIB::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "VIB_potential_choice")
	{
		order = ParameterListT::Choice;
		sub_lists.AddSub("Lennard-Jones_6-12");
		sub_lists.AddSub("Smith-Ferrante");
		sub_lists.AddSub("Gao-Ji");
		sub_lists.AddSub("Gao-Ji_2");
		sub_lists.AddSub("Gao-Nguyen");
	}
	else /* inherited */
		ParameterInterfaceT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* VIB::NewSub(const StringT& name) const
{
	/* try C1FunctionT */
	C1FunctionT* C1 = C1FunctionT::New(name);
	if (C1)
		return C1;
	else /* inherited */
		return ParameterInterfaceT::NewSub(name);
}

/* accept parameter list */
void VIB::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* construct potential */
	const ParameterListT& potential = list.GetListChoice(*this, "VIB_potential_choice");
	fPotential = C1FunctionT::New(potential.Name());
	if (!fPotential) ExceptionT::GeneralFail("VIB::TakeParameterList", "could not construct potential");
	fPotential->TakeParameterList(potential);
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* allocate memory for all the tables */
void VIB::Dimension(int numbonds)
{
	/* length table */
	fLengths.Dimension(numbonds);

	/* potential tables */
	fU.Dimension(numbonds);
	fdU.Dimension(numbonds);
	fddU.Dimension(numbonds);

	/* jacobian table */
	fjacobian.Dimension(numbonds);

	/* STRESS angle tables - by associated stress component */
	fStressTable.Dimension(fNumStress, numbonds);
	  	
	/* MODULI angle tables - using Cauchy symmetry */ 	
	fModuliTable.Dimension(fNumModuli, numbonds);	
}
