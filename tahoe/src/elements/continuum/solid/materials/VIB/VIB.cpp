/* $Id: VIB.cpp,v 1.13.2.2 2004-07-12 16:06:14 paklein Exp $ */
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
#if 0
	/* set potential function */
	int potentialcode;
	in >> potentialcode;	
	switch(potentialcode)
	{
		case C1FunctionT::kLennardJones:
		{	
			double A, B;
			in >> A >> B;
			fPotential = new LennardJones612(A,B);
			break;
		}	
		case C1FunctionT::kSmithFerrante:
		{
			double A, B;
			in >> A >> B;		
			fPotential = new SmithFerrante(A,B);
			break;
		}
		case C1FunctionT::kGaoKlein:
		{
			double A, B, C;
			in >> A >> B >> C;		
			fPotential = new GaoKlein(A,B,C);
			break;
		}
		case C1FunctionT::kQuadraticPot:
		{
			double A, B;
			in >> A;		
	                in >> B;
			fPotential = new ParabolaT(A,B,1.0);
			break;
		}
		case C1FunctionT::kTriantafyllidis:
		{
			double A;
			in >> A;		
			fPotential = new Triantafyllidis(A);
			break;
		}
		case C1FunctionT::kGaoJi:
		{
			double A, B, C;
			in >> A >> B >> C;		
			fPotential = new GaoJi(A,B,C);
			break;
		}
		case C1FunctionT::kGaoJi2:
		{
			double A, B, C;
			in >> A >> B >> C;
			fPotential = new GaoJi2(A,B,C);
			break;
		}
		case C1FunctionT::kGaoVicky:
		{
			double A, B, C, D;
			in >> A >> B >> C >>D;
			fPotential = new GaoVicky(A,B,C,D);
			break;
		}
		case C1FunctionT::kSF2:
		{
			double A, B;
			in >> A >> B;		
			fPotential = new SF2(A,B);
			break;
		}
		default:		
			throw ExceptionT::kBadInputValue;	
	}
	if (!fPotential) throw ExceptionT::kOutOfMemory;
#endif
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
