/* $Id: ParticlePropertyT.cpp,v 1.10 2003-11-04 01:34:29 paklein Exp $ */
#include "ParticlePropertyT.h"
#include "ArrayT.h"
#include <iostream.h>

using namespace Tahoe;

namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<ParticlePropertyT*>::fByteCopy = true; 
DEFINE_TEMPLATE_STATIC const bool ArrayT<ParticlePropertyT>::fByteCopy = false; 
}

/* constructor */
ParticlePropertyT::ParticlePropertyT(void):
	ParameterInterfaceT("particle_property"),
	fMass(0),
	fRange(0)
{

}

/* write properties to output */
void ParticlePropertyT::Write(ostream& out) const
{
	out << " Mass. . . . . . . . . . . . . . . . . . . . . . = " << fMass << '\n';
	out << " Interaction range . . . . . . . . . . . . . . . = " << fRange << '\n';
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

namespace Tahoe {

/* stream extraction operator */
istream& operator>>(istream& in, ParticlePropertyT::TypeT& property)
{
	int i_property;
	in >> i_property;
	switch (i_property)
	{
		case ParticlePropertyT::kHarmonicPair:
			property = ParticlePropertyT::kHarmonicPair;
			break;
		case ParticlePropertyT::kLennardJonesPair:
			property = ParticlePropertyT::kLennardJonesPair;
			break;
		case ParticlePropertyT::kParadynPair:
			property = ParticlePropertyT::kParadynPair;
			break;
		case ParticlePropertyT::kParadynEAM:
			property = ParticlePropertyT::kParadynEAM;
			break;
		case ParticlePropertyT::kMatsuiPair:
			property = ParticlePropertyT::kMatsuiPair;
			break;	
		default:
			ExceptionT::BadInputValue("operator>>ParticlePropertyT::TypeT", 
				"unknown code: %d", i_property);
	}
	return in;
}

} /* namespace Tahoe */
