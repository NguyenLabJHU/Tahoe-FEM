/* $Id: ParticlePropertyT.cpp,v 1.11 2004-07-15 08:29:49 paklein Exp $ */
#include "ParticlePropertyT.h"
#include "ArrayT.h"

using namespace Tahoe;

namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<ParticlePropertyT*>::fByteCopy = true; 
DEFINE_TEMPLATE_STATIC const bool ArrayT<ParticlePropertyT>::fByteCopy = false; 
}

/* constructor */
ParticlePropertyT::ParticlePropertyT(void):
	ParameterInterfaceT("particle_property"),
	fMass(0.0),
	fRange(0.0),
	fNearestNeighbor(0.0)
{

}

/* describe the parameters needed by the interface */
void ParticlePropertyT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	ParameterT mass(fMass, "mass");
	mass.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(mass);
}

/* accept parameter list */
void ParticlePropertyT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	fMass = list.GetParameter("mass");
}
