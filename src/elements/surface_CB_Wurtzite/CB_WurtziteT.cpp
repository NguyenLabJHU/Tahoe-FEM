/* $Id: CB_WurtziteT.cpp,v 1.1 2007-11-08 19:37:46 hspark Exp $ */
/* created: paklein (10/14/1998) */
#include "CB_WurtziteT.h"

#include "WurtziteSolverT.h"
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
CB_WurtziteT::CB_WurtziteT(void):
	ParameterInterfaceT("Wurtzite_CB"),
	fWurtziteSolver(NULL)
{

}

/* destructor */
CB_WurtziteT::~CB_WurtziteT(void) { delete fWurtziteSolver; }
 
/* information about subordinate parameter lists */
void CB_WurtziteT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	NL_E_MatT::DefineSubs(sub_list);
	
	sub_list.AddSub("Wurtzite_CB_solver");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* CB_WurtziteT::NewSub(const StringT& name) const
{
	if (name == "Wurtzite_CB_solver")
		return new WurtziteSolverT(NULL);
	else /* inherited */
		return NL_E_MatT::NewSub(name);
}

/* accept parameter list */
void CB_WurtziteT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	NL_E_MatT::TakeParameterList(list);

	/* dimension work space */
	fXsi.Dimension(kNDOF);
	fC.Dimension(kNSD);
	fPK2.Dimension(kNSD);
	
	/* construct Caucby-Born solver */
	fWurtziteSolver = new WurtziteSolverT(fThermal);
	fWurtziteSolver->TakeParameterList(list.GetList("Wurtzite_CB_solver"));
}

/*************************************************************************
 * Protected
 *************************************************************************/

void CB_WurtziteT::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{
	/* Temporary FD approximation to test new stress */
	//moduli = FSSolidMatT::c_ijkl();
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	/* compute moduli */
	fWurtziteSolver->SetModuli(fC, fXsi, moduli);
}

/* returns the stress corresponding to the given strain - the strain
* is the reduced index Green strain vector, while the stress is the
* reduced index 2nd Piola-Kirchhoff stress vector */
void CB_WurtziteT::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	/* compute stress */
	fWurtziteSolver->SetStress(fC, fXsi, fPK2);
	
	/* shape change */
	PK2.FromMatrix(fPK2);
}

/* returns the strain energy density for the specified strain */
double CB_WurtziteT::ComputeEnergyDensity(const dSymMatrixT& E)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	return fWurtziteSolver->StrainEnergyDensity(fC, fXsi);
}

/*************************************************************************
 * Private
 *************************************************************************/

/* compute the 3D stretch tensor from the 2D reduced index
* strain vector (assuming plane strain */
void CB_WurtziteT::StrainToStretch(const dSymMatrixT& E, dMatrixT& C)
{
	/* shape change */
	E.ToMatrix(C);

	/* convert */
	C *= 2.0;
	C.PlusIdentity();
}
