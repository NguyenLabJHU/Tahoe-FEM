/* $Id: ParticlePairT.cpp,v 1.4 2002-11-22 01:49:45 paklein Exp $ */
/* created: paklein (10/22/1996) */
#include "ParticlePairT.h"
#include "ifstreamT.h"
#include "eControllerT.h"

/* parameters */
const int kMemoryHeadRoom = 15; /* percent */

/* constructor */
ParticlePairT::ParticlePairT(const ElementSupportT& support, const FieldT& field):
	ParticleT(support, field),
	fReNeighborCounter(0),
	fNeighbors(kMemoryHeadRoom),
	fEqnos(kMemoryHeadRoom),
	fNeighborDistance(-1),
	fReNeighborIncr(-1)
{
	/* read parameters */
	ifstreamT& in = ElementSupport().Input();

	/* read parameters */
	in >> fNeighborDistance;
	in >> fReNeighborIncr;
	
	/* checks */
	if (fNeighborDistance < kSmall || fReNeighborIncr < 0)
		ExceptionT::BadInputValue("ParticlePairT::ParticlePairT");
}

/* initialization */
void ParticlePairT::Initialize(void)
{
	/* inherited */
	ParticleT::Initialize();
	
	/* set the neighborlists */
	SetConfiguration(true);
}

/* collecting element group equation numbers */
void ParticlePairT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_1)

	/* dimension equations array */
	fEqnos.Configure(fNeighbors, NumDOF());

	/* get local equations numbers */
	Field().SetLocalEqnos(fNeighbors, fEqnos);

	/* add to list of equation numbers */
	eq_2.Append(&fEqnos);
}

/* collecting element geometry connectivities */
void ParticlePairT::ConnectsX(AutoArrayT<const iArray2DT*>& connects) const
{
	/* NOTE: do not add anything to the geometry connectivity list */
#pragma unused(connects)
}

/* collecting element field connectivities */
void ParticlePairT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
#pragma unused(connects_1)
	connects_2.AppendUnique(&fNeighbors);
}

/* trigger reconfiguration */
GlobalT::RelaxCodeT ParticlePairT::RelaxSystem(void)
{
	/* reset neighbor lists */
	if (SetConfiguration())
		return GlobalT::kReEQ;
	else
		return GlobalT::kNoRelax;
}

/* close current time increment */
void ParticlePairT::CloseStep(void)
{
	/* inherited */
	ParticleT::CloseStep();
	
	/* increment counter */
	fReNeighborCounter++;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* form group contribution to the stiffness matrix */
void ParticlePairT::LHSDriver(void)
{
	/* time integration parameters */
	double constK = 0.0;
	double constM = 0.0;
	int formK = fController->FormK(constK);
	int formM = fController->FormM(constM);

//TEMP - no stiffness implemented
if (formK) ExceptionT::GeneralFail("ParticlePairT::LHSDriver", "stiffness not implemented");
	
}

/* form group contribution to the residual */
void ParticlePairT::RHSDriver(void)
{
	/* time integration parameters */
	double constMa = 0.0;
	double constKd = 0.0;
	int formMa = fController->FormMa(constMa);
	int formKd = fController->FormKd(constKd);


}

/* set neighborlists */
bool ParticlePairT::SetConfiguration(bool force)
{
	if (force || fReNeighborCounter == fReNeighborIncr)
	{
		/* reset neighbor lists */
		GenerateNeighborList(fGlobalTag, fNeighborDistance, false, fNeighbors);
	
		/* reset counter */
		fReNeighborCounter = 0;
		
		return true;
	}
	else return false;
}
