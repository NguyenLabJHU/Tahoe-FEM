/* $Id: ModCB3DT.cpp,v 1.9 2004-07-15 08:28:36 paklein Exp $ */
/* created: paklein (10/14/1998) */
#include "ModCB3DT.h"

#include "ModCBSolverT.h"
#include "dMatrixT.h"

using namespace Tahoe;

/* material parameters */
const int kNSD  = 3;
const int kNDOF = 3;

const double sqrt2 = sqrt(2.0);
const double sqrt3 = sqrt(3.0);

/* plane codes - for crystal axes rotated wrt global axes*/
const int	kDCnatural 	= 0;
const int 	kDC110		= 1;
const int	kDC111		= 2;

/* constructor */
ModCB3DT::ModCB3DT(void):
	ParameterInterfaceT("Cauchy-Born_diamond"),
	fModCBSolver(NULL)
{

}

/* destructor */
ModCB3DT::~ModCB3DT(void) { delete fModCBSolver; }
 
/* information about subordinate parameter lists */
void ModCB3DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	NL_E_MatT::DefineSubs(sub_list);
	
	sub_list.AddSub("mod_Cauchy-Born_solver");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ModCB3DT::NewSub(const StringT& name) const
{
	if (name == "mod_Cauchy-Born_solver")
		return new ModCBSolverT(NULL);
	else /* inherited */
		return NL_E_MatT::NewSub(name);
}

/* accept parameter list */
void ModCB3DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	NL_E_MatT::TakeParameterList(list);

	/* dimension work space */
	fXsi.Dimension(kNDOF);
	fC.Dimension(kNSD);
	fPK2.Dimension(kNSD);
	
	/* construct Caucby-Born solver */
	fModCBSolver = new ModCBSolverT(fThermal);
	fModCBSolver->TakeParameterList(list.GetList("mod_Cauchy-Born_solver"));
}

/*************************************************************************
 * Protected
 *************************************************************************/

void ModCB3DT::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	/* compute moduli */
	fModCBSolver->SetModuli(fC, fXsi, moduli);
}

/* returns the stress corresponding to the given strain - the strain
* is the reduced index Green strain vector, while the stress is the
* reduced index 2nd Piola-Kirchhoff stress vector */
void ModCB3DT::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	/* compute stress */
	fModCBSolver->SetStress(fC, fXsi, fPK2);
	
	/* shape change */
	PK2.FromMatrix(fPK2);
}

/* returns the strain energy density for the specified strain */
double ModCB3DT::ComputeEnergyDensity(const dSymMatrixT& E)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	return fModCBSolver->StrainEnergyDensity(fC, fXsi);
}

/*************************************************************************
 * Private
 *************************************************************************/

/* compute the 3D stretch tensor from the 2D reduced index
* strain vector (assuming plane strain */
void ModCB3DT::StrainToStretch(const dSymMatrixT& E, dMatrixT& C)
{
	/* shape change */
	E.ToMatrix(C);

	/* convert */
	C *= 2.0;
	C.PlusIdentity();
}
