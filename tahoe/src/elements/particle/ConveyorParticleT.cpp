/* $Id: ConveyorParticleT.cpp,v 1.1.2.1 2003-09-18 21:03:36 cjkimme Exp $ */
#include "ConveyorParticleT.h"

#include "ParticlePropertyT.h"

using namespace Tahoe;


/* constructors */
ConveyorParticleT::ConveyorParticleT(const ElementSupportT& support, const FieldT& field):
	ParticlePairT(support, field)
{
	// Do nothing 
}

/* destructor */
ConveyorParticleT::~ConveyorParticleT(void)
{
	// Do nothing
}

/* initialization */
void ConveyorParticleT::Initialize(void)
{
	const char caller[] = "ConveyorParticleT::Initialize";

	/* inherited */
	ParticlePairT::Initialize();
	
}

/* write restart data to the output stream */
void ConveyorParticleT::WriteRestart(ostream& out) const
{
	/* inherited */
	ParticlePairT::WriteRestart(out);
	
}

/* read restart data to the output stream */
void ConveyorParticleT::ReadRestart(istream& in)
{
	/* inherited */
	ParticlePairT::ReadRestart(in);
	
}

void ConveyorParticleT::CreateNoninteractingAtoms(iArrayT& topAtoms, iArrayT& bottomAtoms)
{
	fType_Actual.Dimension(fType.Length());
	fType.CopyInto(fType_Actual);
	
	fNumTypes_Actual = fNumTypes;
	fNumTypes += 2;

	int topType = fNumTypes_Actual;
	int bottomType = fNumTypes_Actual + 1;

	for (int i = 0; i < topAtoms.Length(); i++)
		fType[topAtoms[i]] = topType;
	for (int i = 0; i < bottomAtoms.Length(); i++)
		fType[bottomAtoms[i]] = bottomType;

	fPropertiesMap_Actual.Dimension(fNumTypes_Actual);
	fPropertiesMap_Actual = fPropertiesMap;
	
	fPropertiesMap.Dimension(fNumTypes);
	fPropertiesMap = fPropertiesMap_Actual(0,0); // WARNING :: ASSUMING ONLY 1 KIND OF INTERACTION
	fPropertiesMap(topType,bottomType) = fPropertiesMap(bottomType,topType) = ParticlePropertyT::kNull;

	return ;
}
