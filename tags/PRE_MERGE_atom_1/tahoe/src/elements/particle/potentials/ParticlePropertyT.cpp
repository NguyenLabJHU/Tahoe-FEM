/* $Id: ParticlePropertyT.cpp,v 1.1 2002-11-25 07:19:46 paklein Exp $ */
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
	out << " Interation range. . . . . . . . . . . . . . . . = " << fRange << '\n';
}
