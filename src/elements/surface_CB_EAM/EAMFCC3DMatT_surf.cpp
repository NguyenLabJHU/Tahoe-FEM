/* $Id: EAMFCC3DMatT_surf.cpp,v 1.2 2006-05-29 17:22:56 paklein Exp $ */
/* created: paklein (10/25/1998) */
#include "EAMFCC3DMatT_surf.h"

#include "EAMFCC3DSym_surf.h"
#include "dMatrixT.h"

#include <math.h>

using namespace Tahoe;

/* constructor */
EAMFCC3DMatT_surf::EAMFCC3DMatT_surf(void):
	ParameterInterfaceT("FCC_EAM"),
	fEAM(NULL)
{

}

/* destructor */
EAMFCC3DMatT_surf::~EAMFCC3DMatT_surf(void) { delete fEAM; }

/* describe the parameters needed by the interface */
void EAMFCC3DMatT_surf::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	NL_E_MatT::DefineParameters(list);

	/* number of neighbor shells */
	ParameterT n_shells(ParameterT::Integer, "shells");
	n_shells.AddLimit(1, LimitT::LowerInclusive);
	list.AddParameter(n_shells);
	
	/* surface normal */
	ParameterT normal(ParameterT::Integer, "normal_code");
	normal.AddLimit(0, LimitT::LowerInclusive);
	normal.AddLimit(5, LimitT::UpperInclusive);
	list.AddParameter(normal);
}

/* describe the parameters needed by the interface */
void EAMFCC3DMatT_surf::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	NL_E_MatT::DefineSubs(sub_list);

	/* Cauchy-Born EAM parameters */
	sub_list.AddSub("FCC_EAM_Cauchy-Born");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* EAMFCC3DMatT_surf::NewSub(const StringT& name) const
{
	if (name == "FCC_EAM_Cauchy-Born")
		return new EAMFCC3DSym_surf(0, 0);
	else /* inherited */
		return NL_E_MatT::NewSub(name);
}

/* accept parameter list */
void EAMFCC3DMatT_surf::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	NL_E_MatT::TakeParameterList(list);

	/* construct Cauchy-Born EAM solver */
	int shells = list.GetParameter("shells");
	int normal_code = list.GetParameter("normal_code");
	fEAM = new EAMFCC3DSym_surf(shells, normal_code);
	fEAM->TakeParameterList(list.GetList("FCC_EAM_Cauchy-Born"));
}

/*************************************************************************
 * Private
 *************************************************************************/

void EAMFCC3DMatT_surf::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{
	/* EAM solver */
	fEAM->Moduli(moduli, E);
}

/* returns the stress corresponding to the given strain - the strain
* is the reduced index Green strain vector, while the stress is the
* reduced index 2nd Piola-Kirchhoff stress vector */
void EAMFCC3DMatT_surf::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	/* EAM solver */
	fEAM->SetStress(E, PK2);
}

/* returns the strain energy density for the specified strain */
double EAMFCC3DMatT_surf::ComputeEnergyDensity(const dSymMatrixT& E)
{
	/* EAM solver */
	return fEAM->EnergyDensity(E);
}
