/* $Id: CB_TersoffDimerT.cpp,v 1.1 2007-11-09 16:53:55 hspark Exp $ */
/* created: paklein (10/14/1998) */
#include "CB_TersoffDimerT.h"

#include "TersoffDimerSolverT.h"
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
CB_TersoffDimerT::CB_TersoffDimerT(void):
	ParameterInterfaceT("TersoffDimer_CB"),
	fTersoffDimerSolver(NULL)
{

}

/* destructor */
CB_TersoffDimerT::~CB_TersoffDimerT(void) { delete fTersoffDimerSolver; }
 
/* information about subordinate parameter lists */
void CB_TersoffDimerT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	NL_E_MatT::DefineSubs(sub_list);
	
	sub_list.AddSub("TersoffDimer_CB_solver");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* CB_TersoffDimerT::NewSub(const StringT& name) const
{
	if (name == "TersoffDimer_CB_solver")
		return new TersoffDimerSolverT(NULL);
	else /* inherited */
		return NL_E_MatT::NewSub(name);
}

/* accept parameter list */
void CB_TersoffDimerT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	NL_E_MatT::TakeParameterList(list);

	/* dimension work space */
	fXsi.Dimension(kNDOF);
	fC.Dimension(kNSD);
	fPK2.Dimension(kNSD);
	
	/* construct Caucby-Born solver */
	fTersoffDimerSolver = new TersoffDimerSolverT(fThermal);
	fTersoffDimerSolver->TakeParameterList(list.GetList("TersoffDimer_CB_solver"));
}

/*************************************************************************
 * Protected
 *************************************************************************/

void CB_TersoffDimerT::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{
	/* Temporary FD approximation to test new stress */
	//moduli = FSSolidMatT::c_ijkl();
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	/* compute moduli */
	fTersoffDimerSolver->SetModuli(fC, fXsi, moduli);
}

/* returns the stress corresponding to the given strain - the strain
* is the reduced index Green strain vector, while the stress is the
* reduced index 2nd Piola-Kirchhoff stress vector */
void CB_TersoffDimerT::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	/* compute stress */
	fTersoffDimerSolver->SetStress(fC, fXsi, fPK2);
	
	/* shape change */
	PK2.FromMatrix(fPK2);
}

/* returns the strain energy density for the specified strain */
double CB_TersoffDimerT::ComputeEnergyDensity(const dSymMatrixT& E)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	return fTersoffDimerSolver->StrainEnergyDensity(fC, fXsi);
}

/*************************************************************************
 * Private
 *************************************************************************/

/* compute the 3D stretch tensor from the 2D reduced index
* strain vector (assuming plane strain */
void CB_TersoffDimerT::StrainToStretch(const dSymMatrixT& E, dMatrixT& C)
{
	/* shape change */
	E.ToMatrix(C);

	/* convert */
	C *= 2.0;
	C.PlusIdentity();
}
