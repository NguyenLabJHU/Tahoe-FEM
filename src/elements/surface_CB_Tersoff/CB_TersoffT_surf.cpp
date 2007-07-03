/* $Id: CB_TersoffT_surf.cpp,v 1.2 2007-07-03 17:47:23 hspark Exp $ */
/* created: paklein (10/14/1998) */
#include "CB_TersoffT_surf.h"

#include "TersoffSolverT_surf.h"
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
CB_TersoffT_surf::CB_TersoffT_surf(void):
	ParameterInterfaceT("Tersoff_CB_surf"),
	fSurfaceThickness(-1),
	fTersoffSolver_surf(NULL)
{

}

/* destructor */
CB_TersoffT_surf::~CB_TersoffT_surf(void) { delete fTersoffSolver_surf; }
 
/* information about subordinate parameter lists */
void CB_TersoffT_surf::DefineSubs(SubListT& sub_list) const
{
	cout << "CB_TersoffT_surf::DefineSubs" << endl;
	/* inherited */
	NL_E_MatT::DefineSubs(sub_list);
	
	sub_list.AddSub("Tersoff_CB_solver_surf");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* CB_TersoffT_surf::NewSub(const StringT& name) const
{
	cout << "CB_TersoffT_surf::NewSub" << endl;
	if (name == "Tersoff_CB_solver_surf")
		return new TersoffSolverT_surf(NULL,0);
	else /* inherited */
		return NL_E_MatT::NewSub(name);
}

/* accept parameter list */
void CB_TersoffT_surf::TakeParameterList(const ParameterListT& list)
{
	cout << "CB_TersoffT_surf::TakeParameterList" << endl;
	/* inherited */
	NL_E_MatT::TakeParameterList(list);

	/* dimension work space */
	fXsi.Dimension(kNDOF);
	fC.Dimension(kNSD);
	fPK2.Dimension(kNSD);
	
	cout << "Before construction" << endl;
	/* construct Surface Cauchy-Born solver */
	int normal_code = list.GetParameter("normal_code");
	fTersoffSolver_surf = new TersoffSolverT_surf(fThermal,normal_code);
	fTersoffSolver_surf->TakeParameterList(list.GetList("Tersoff_CB_solver_surf"));
	cout << "After construction" << endl;
	
	/* Get surface thickness */
	fSurfaceThickness = fTersoffSolver_surf->SurfaceThickness();
}

/*************************************************************************
 * Protected
 *************************************************************************/

void CB_TersoffT_surf::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{
	cout << "CB_TersoffT_surf::ComputeModuli" << endl;
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	/* compute moduli */
	fTersoffSolver_surf->SetModuli(fC, fXsi, moduli);
}

/* returns the stress corresponding to the given strain - the strain
* is the reduced index Green strain vector, while the stress is the
* reduced index 2nd Piola-Kirchhoff stress vector */
void CB_TersoffT_surf::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	/* compute stress */
	fTersoffSolver_surf->SetStress(fC, fXsi, fPK2);
	
	/* shape change */
	PK2.FromMatrix(fPK2);
}

/* returns the strain energy density for the specified strain */
double CB_TersoffT_surf::ComputeEnergyDensity(const dSymMatrixT& E)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	return fTersoffSolver_surf->StrainEnergyDensity(fC, fXsi);
}

/*************************************************************************
 * Private
 *************************************************************************/

/* compute the 3D stretch tensor from the 2D reduced index
* strain vector (assuming plane strain */
void CB_TersoffT_surf::StrainToStretch(const dSymMatrixT& E, dMatrixT& C)
{
	/* shape change */
	E.ToMatrix(C);

	/* convert */
	C *= 2.0;
	C.PlusIdentity();
}
