/* $Id: ParticlePropertyT.cpp,v 1.4 2003-04-09 20:22:28 cjkimme Exp $ */
#include "ParticlePropertyT.h"
#include "ArrayT.h"
#include <iostream.h>

using namespace Tahoe;

namespace Tahoe {
const bool ArrayT<ParticlePropertyT*>::fByteCopy = true; 
const bool ArrayT<ParticlePropertyT>::fByteCopy = false; 
}

/* constructor */
ParticlePropertyT::ParticlePropertyT(void):
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
		default:
			ExceptionT::BadInputValue("operator>>ParticlePropertyT::TypeT", 
				"unknown code: %d", i_property);
	}
	return in;
}

/* stream extraction operator */
istream& operator>>(istream& in, ParticlePropertyT::ThermostatT& property)
{
	int i_property;
	in >> i_property;
	switch (i_property)
	{
		case ParticlePropertyT::kFreeParticle:
			property = ParticlePropertyT::kFreeParticle;
			break;
		case ParticlePropertyT::kDampedParticle:
			property = ParticlePropertyT::kDampedParticle;
			break;
		case ParticlePropertyT::kLangevinParticle:
			property = ParticlePropertyT::kLangevinParticle;
			break;
		default:
			ExceptionT::BadInputValue("operator>>ParticlePropertyT::ThermostatT", 
				"unknown code: %d", i_property);
	}
	return in;
}

} /* namespace Tahoe */
