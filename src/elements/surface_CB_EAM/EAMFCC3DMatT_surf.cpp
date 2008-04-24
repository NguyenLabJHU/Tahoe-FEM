/* $Id: EAMFCC3DMatT_surf.cpp,v 1.7 2008-04-24 22:26:26 hspark Exp $ */
/* created: paklein (10/25/1998) */
#include "EAMFCC3DMatT_surf.h"

#include "EAMFCC3DSym_surf.h"
#include "dMatrixT.h"

#include <math.h>

using namespace Tahoe;

/* constructor */
EAMFCC3DMatT_surf::EAMFCC3DMatT_surf(void):
	ParameterInterfaceT("FCC_EAM_Surf"),
	fSurfaceThickness(-1),
	fAlpha(0.0),
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
	
	/* TEMP - GET SURFACE THICKNESS FROM EAMFCC3D_SURF */
	fSurfaceThickness = fEAM->SurfaceThickness();
	
	/* reset density from the atomistic parameters */
	fDensity = fEAM->Density();
	
	/* HSP ADDED 4/24/08 */
	fSS0 = FSSolidMatT::c_ijkl();
	/* Hopefully will return 0 strain stiffness since called initially */
	fAlpha = 0.5;
	fSS0*=fAlpha;
}

/*************************************************************************
 * Private
 *************************************************************************/

void EAMFCC3DMatT_surf::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli)
{
	/* Temporarily override for finite difference approximation */
	moduli = FSSolidMatT::c_ijkl();

	/* EAM solver */
	//fEAM->Moduli(moduli, E);
	
	/* Subtract off strain-dependent part */
	moduli-=fSS0;
}

/* returns the stress corresponding to the given strain - the strain
* is the reduced index Green strain vector, while the stress is the
* reduced index 2nd Piola-Kirchhoff stress vector */
void EAMFCC3DMatT_surf::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2)
{
	/* EAM solver */
	fEAM->SetStress(E, PK2);
	
	/* Subtract off strain-dependent part */
	dArrayT temp(6), temp2(6);
	dSymMatrixT product(3);
	/* Note that Miehe uses 11, 22, 33, 12, 23, 13 notation */
	temp[0] = E(0,0);	// epsilon 11
	temp[1] = E(1,1);	// epsilon 22
	temp[2] = E(2,2);	// epsilon 33
	temp[3] = E(0,1);	// epsilon 12
	temp[4] = E(1,2);	// epsilon 23
	temp[5] = E(0,2);	// epsilon 13
	
	/* Assuming that fSS0 is right (spatial?  material?) tangent modulus */
	fSS0.Multx(temp,temp2);
	product(0,0) = temp2[0];	// sigma 11
	product(1,1) = temp2[1];	// sigma 22
	product(2,2) = temp2[2];	// sigma 33
	product(1,2) = temp2[4];	// sigma 23
	product(0,2) = temp2[5];	// sigma 13
	product(0,1) = temp2[3];	// sigma 12
	product(1,0) = temp2[3];	// sigma 21
	product(2,0) = temp2[5];	// sigma 31
	product(2,1) = temp2[4];	// sigma 32	
	
	PK2-=product;
}

/* returns the strain energy density for the specified strain */
double EAMFCC3DMatT_surf::ComputeEnergyDensity(const dSymMatrixT& E)
{
	/* EAM solver */
	return fEAM->EnergyDensity(E);
}
