/* $Id: ModCB2DT.cpp,v 1.8.46.4 2004-06-25 01:30:28 paklein Exp $ */
/* created: paklein (05/31/1997) */
#include "ModCB2DT.h"

#include "ModCBSolverT.h"
#include "dMatrixT.h"

using namespace Tahoe;

/* material parameters */
const int knsd = 2;

/* constructor */
ModCB2DT::ModCB2DT(void):
	ParameterInterfaceT("Cauchy-Born_diamond_2D"),
	fModCBSolver(NULL)
{

}

/* destructor */
ModCB2DT::~ModCB2DT(void) { delete fModCBSolver; }

/* describe the parameters needed by the interface */
void ModCB2DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	NL_E_MatT::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStrain);
}

/* information about subordinate parameter lists */
void ModCB2DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	NL_E_MatT::DefineSubs(sub_list);
	
	sub_list.AddSub("mod_Cauchy-Born_solver");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ModCB2DT::NewSub(const StringT& list_name) const
{
	if (list_name == "mod_Cauchy-Born_solver")
		return new ModCBSolverT(NULL);
	else /* inherited */
		return NL_E_MatT::NewSub(list_name);
}

/* accept parameter list */
void ModCB2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	NL_E_MatT::TakeParameterList(list);

	/* dimension work space */
	fCij3D.Dimension(dSymMatrixT::NumValues(3));
	fXsi.Dimension(3); 
	fStretch3D.Dimension(3);
	fStretch2D.Dimension(2);
	fStress3D.Dimension(3);
	
	/* construct Caucby-Born solver */
	fModCBSolver = new ModCBSolverT(fThermal);
	fModCBSolver->TakeParameterList(list.GetList("mod_Cauchy-Born_solver"));
}

/*************************************************************************
 * Protected
 *************************************************************************/

void ModCB2DT::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* convert 2D Green strain to 3D stretch */
	StrainToStretch(E, fStretch3D);

	/* compute moduli */
	fModCBSolver->SetModuli(fStretch3D, fXsi, fCij3D);
	
	/* convert to 2D */
	moduli.Rank4ReduceFrom3D(fCij3D);
}

/*
* Returns the stress corresponding to the given strain - the strain
* is the reduced index Green strain vector, while the stress is the
* reduced index 2nd Piola-Kirchhoff stress vector.
*/
void ModCB2DT::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* convert 2D Green strain to 3D stretch */
	StrainToStretch(E, fStretch3D);

	/* compute stress */
	fModCBSolver->SetStress(fStretch3D, fXsi, fStress3D);
	
	/* 3D tensor2 -> 2D vector */
	PK2[0] = fStress3D(0,0);
	PK2[1] = fStress3D(1,1);
	PK2[2] = fStress3D(0,1);	
}

/* strain energy density for the specified strain */
double ModCB2DT::ComputeEnergyDensity(const dSymMatrixT& E)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* convert 2D Green strain to 3D stretch */
	StrainToStretch(E, fStretch3D);

	return( fModCBSolver->StrainEnergyDensity(fStretch3D,fXsi) );
}

/*************************************************************************
* Private
*************************************************************************/

/*
* Compute the 3D stretch tensor from the 2D reduced index
* strain vector (assuming plane strain)
*/
void ModCB2DT::StrainToStretch(const dSymMatrixT& strain2D,
	dMatrixT& stretch3D)
{
	stretch3D = 0.0;
	
	stretch3D(0,0) = 1.0 + 2.0*strain2D[0];
	stretch3D(1,1) = 1.0 + 2.0*strain2D[1];
	stretch3D(2,2) = 1.0;
	stretch3D(0,1) = stretch3D(1,0) = 2.0*strain2D[2];	
}