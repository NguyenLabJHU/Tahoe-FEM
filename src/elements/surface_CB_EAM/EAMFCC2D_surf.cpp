/* $Id: EAMFCC2D_surf.cpp,v 1.1 2006-05-21 15:55:19 hspark Exp $ */
/* created: paklein (12/09/1996) */
#include "EAMFCC2D_surf.h"

#include "EAMFCC3DSym_surf.h"
#include "dMatrixT.h"

#include <math.h>

using namespace Tahoe;

/* material parameters */
const int knsd = 2;

/* constructor */
EAMFCC2D_surf::EAMFCC2D_surf(void):
	ParameterInterfaceT("FCC_EAM_2D"),
	fEAM(NULL)
{
	/* reset default */
	fConstraint = kPlaneStrain;
}

/* destructor */
EAMFCC2D_surf::~EAMFCC2D_surf(void) { delete fEAM; }

/* describe the parameters needed by the interface */
void EAMFCC2D_surf::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	NL_E_MatT::DefineSubs(sub_list);

	/* Cauchy-Born EAM parameters */
	sub_list.AddSub("FCC_EAM_Cauchy-Born");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* EAMFCC2D_surf::NewSub(const StringT& name) const
{
	if (name == "FCC_EAM_Cauchy-Born")
		return new EAMFCC3DSym_surf;
	else /* inherited */
		return NL_E_MatT::NewSub(name);
}

/* accept parameter list */
void EAMFCC2D_surf::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	NL_E_MatT::TakeParameterList(list);

	/* construct Cauchy-Born EAM solver */
	fEAM = new EAMFCC3DSym_surf;
	fEAM->TakeParameterList(list.GetList("FCC_EAM_Cauchy-Born"));
}

/*************************************************************************
 * Private
 *************************************************************************/

void EAMFCC2D_surf::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{
	/* EAM solver */
	fEAM->Moduli(moduli, E);
}

/* returns the stress corresponding to the given strain - the strain
* is the reduced index Green strain vector, while the stress is the
* reduced index 2nd Piola-Kirchhoff stress vector */
void EAMFCC2D_surf::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	/* EAM solver */
	fEAM->SetStress(E, PK2);
}

/* returns the strain energy density for the specified strain */
double EAMFCC2D_surf::ComputeEnergyDensity(const dSymMatrixT& E)
{
	/* EAM solver */
	return fEAM->EnergyDensity(E);
}
