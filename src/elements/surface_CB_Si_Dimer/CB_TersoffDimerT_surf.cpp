/* $Id: CB_TersoffDimerT_surf.cpp,v 1.1 2007-11-09 16:53:55 hspark Exp $ */
/* created: paklein (10/14/1998) */
#include "CB_TersoffDimerT_surf.h"

#include "TersoffDimerSolverT_surf.h"
#include "dMatrixT.h"

using namespace Tahoe;

/* material parameters */
const int kNSD  = 3;
const int kNDOF = 4*kNSD; /* 4 lattice displacements */

const double sqrt2 = sqrt(2.0);
const double sqrt3 = sqrt(3.0);

/* plane codes - for crystal axes rotated wrt global axes*/
const int	kDCnatural 	= 0;
const int 	kDC110		= 1;
const int	kDC111		= 2;

/* constructor */
CB_TersoffDimerT_surf::CB_TersoffDimerT_surf(void):
	ParameterInterfaceT("TersoffDimer_CB_surf"),
	fSurfaceThickness(-1),
	fTersoffDimerSolver_surf(NULL)
{

}

/* destructor */
CB_TersoffDimerT_surf::~CB_TersoffDimerT_surf(void) { delete fTersoffDimerSolver_surf; }
 
/* information about subordinate parameter lists */
void CB_TersoffDimerT_surf::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	NL_E_MatT::DefineSubs(sub_list);
	
	sub_list.AddSub("TersoffDimer_CB_solver");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* CB_TersoffDimerT_surf::NewSub(const StringT& name) const
{
	if (name == "TersoffDimer_CB_solver")
		return new TersoffDimerSolverT_surf(NULL,0);
	else /* inherited */
		return NL_E_MatT::NewSub(name);
}

/* accept parameter list */
void CB_TersoffDimerT_surf::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	NL_E_MatT::TakeParameterList(list);

	/* dimension work space */
	fXsi.Dimension(kNDOF);
	fC.Dimension(kNSD);
	fPK2.Dimension(kNSD);
	
	/* construct Surface Cauchy-Born solver */
	int normal_code = list.GetParameter("normal_code");
	fTersoffDimerSolver_surf = new TersoffDimerSolverT_surf(fThermal,normal_code);
	fTersoffDimerSolver_surf->TakeParameterList(list.GetList("TersoffDimer_CB_solver"));
	
	/* Get surface thickness */
	fSurfaceThickness = fTersoffDimerSolver_surf->SurfaceThickness();
}

/* return the number of constitutive model output parameters */
int CB_TersoffDimerT_surf::NumOutputVariables(void) const { return kNDOF; }

/* return the labels for model output parameters */
void CB_TersoffDimerT_surf::OutputLabels(Tahoe::ArrayT<StringT>& labels) const
{
	int n_atom = kNDOF/kNSD;
	StringT xsi = "Xsi";
	labels.Dimension(kNDOF);
	int dex = 0;
	for (int i = 1; i <= n_atom; i++) {
		for (int j = 1; j <= kNSD; j++) {
			StringT label = xsi;
			label.Append(i);
			label.Append("_", j);
			labels[dex++] = label;
		}
	}
}

/* return material output variables */
void CB_TersoffDimerT_surf::ComputeOutput(Tahoe::dArrayT& output)
{
	/* set internal DOF's by calculating the stress */
	s_ij();

	/* copy */
	output = fXsi;
}

/*************************************************************************
 * Protected
 *************************************************************************/

void CB_TersoffDimerT_surf::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	/* compute moduli */
	fTersoffDimerSolver_surf->SetModuli(fC, fXsi, moduli);
}

/* returns the stress corresponding to the given strain - the strain
* is the reduced index Green strain vector, while the stress is the
* reduced index 2nd Piola-Kirchhoff stress vector */
void CB_TersoffDimerT_surf::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	/* compute stress */
	fTersoffDimerSolver_surf->SetStress(fC, fXsi, fPK2);
	
	/* shape change */
	PK2.FromMatrix(fPK2);
}

/* returns the strain energy density for the specified strain */
double CB_TersoffDimerT_surf::ComputeEnergyDensity(const dSymMatrixT& E)
{
	/* zero initial guess */
	fXsi = 0.0;

	/* E -> C */
	StrainToStretch(E, fC);

	return fTersoffDimerSolver_surf->StrainEnergyDensity(fC, fXsi);
}

/*************************************************************************
 * Private
 *************************************************************************/

/* compute the 3D stretch tensor from the 2D reduced index
* strain vector (assuming plane strain */
void CB_TersoffDimerT_surf::StrainToStretch(const dSymMatrixT& E, dMatrixT& C)
{
	/* shape change */
	E.ToMatrix(C);

	/* convert */
	C *= 2.0;
	C.PlusIdentity();
}
