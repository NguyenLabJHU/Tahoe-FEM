/* $Id: ParticlePairT.cpp,v 1.5 2002-11-25 07:19:45 paklein Exp $ */
#include "ParticlePairT.h"
#include "PairPropertyT.h"
#include "fstreamT.h"
#include "eControllerT.h"

/* pair property types */
#include "LennardJonesPairT.h"
#include "HarmonicPairT.h"

/* parameters */
const int kMemoryHeadRoom = 15; /* percent */

/* constructor */
ParticlePairT::ParticlePairT(const ElementSupportT& support, const FieldT& field):
	ParticleT(support, field),
	fNeighbors(kMemoryHeadRoom),
	fEqnos(kMemoryHeadRoom),
	fNeighborDistance(-1)
{
	/* input stream */
	ifstreamT& in = ElementSupport().Input();

	/* read parameters */
	in >> fNeighborDistance;
	in >> fReNeighborIncr;

	/* checks */
	if (fNeighborDistance < kSmall || fReNeighborIncr < 0) ExceptionT::BadInputValue("ParticlePairT::ParticlePairT");
}

/* destructor */
ParticlePairT::~ParticlePairT(void)
{
	/* free properties list */
	for (int i = 0; i < fProperties.Length(); i++)
		delete fProperties[i];
}

/* initialization */
void ParticlePairT::Initialize(void)
{
	/* inherited */
	ParticleT::Initialize();
	
	/* verbose */
	if (ElementSupport().PrintInput())
	{
		ofstreamT& out = ElementSupport().Output();
		out << "\n Neighbor lists:\n\n";
		fNeighbors.WriteNumbered(out);
		out.flush();
	}
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

	/* assemble particle mass */
	if (formM) {

		/* collect mass per particle */
		dArrayT mass(fNumTypes);
		for (int i = 0; i < fNumTypes; i++)
			mass[i] = fProperties[fPropertiesMap(i,i)]->Mass();
		mass *= constM;
	
		AssembleParticleMass(mass);
	}
}

/* form group contribution to the residual */
void ParticlePairT::RHSDriver(void)
{
	/* time integration parameters */
	double constMa = 0.0;
	double constKd = 0.0;
	int formMa = fController->FormMa(constMa);
	int formKd = fController->FormKd(constKd);

//TEMP - interial force not implemented
if (formMa) ExceptionT::GeneralFail("ParticlePairT::RHSDriver", "inertial force not implemented");

	/* assembly information */
	const ElementSupportT& support = ElementSupport();
	int group = Group();
	int ndof = NumDOF();
	
	/* global coordinates */
	const dArray2DT& coords = support.CurrentCoordinates();

	/* pair properties function pointers */
	int current_property = -1;
	PairPropertyT::ForceFunction force_function = NULL;

	/* run through neighbor list */
	iArrayT eqs;
	iArrayT neighbors;
	dArrayT x_i, x_j, r_ij;
	for (int i = 0; i < fGlobalTag.Length(); i++)
	{
		/* global tag */
		int tag_i = fGlobalTag[i];
	
		/* initialize */
		fForce.Dimension(fEqnos.MinorDim(i));
		fEqnos.RowAlias(i, eqs);
		fNeighbors.RowAlias(i, neighbors);
		double* f_i = fForce.Pointer();
		double* f_j = fForce.Pointer(ndof);
		
		/* run though neighbors for one atom - first neighbor is self */
		coords.RowAlias(tag_i, x_i);
		for (int j = 1; j < neighbors.Length(); j++)
		{
			/* global tag */
			int tag_j = fGlobalTag[j];
			
			/* set pair property */
			int property = fPropertiesMap(tag_i, tag_j);
			if (property != current_property)
			{
				force_function = fProperties
				current_property = property;
			}
		
			/* global coordinates */
			coords.RowAlias(tag_j, x_j);
		
			/* connecting vector */
			r_ij.DiffOf(x_j, x_i);
			double r = r_ij.Magnitude();
			
		
		
		}

		/* assemble */
		support.AssembleRHS(group, fForce, eqs);
	}
}

/* set neighborlists */
void ParticlePairT::SetConfiguration(void)
{
	/* reset neighbor lists */
	GenerateNeighborList(fGlobalTag, fNeighborDistance, false, fNeighbors);
}

/* construct the list of properties from the given input stream */
void ParticlePairT::EchoProperties(ifstreamT& in, ofstreamT& out)
{
	/* read potentials */
	int num_potentials = -1;
	in >> num_potentials;
	fProperties.Dimension(num_potentials);
	prop_list = NULL;
	for (int i = 0; i < fProperties.Length(); i++)
	{
		ParticleT::PropertyT property;
		in >> property;
		switch (property)
		{
			case ParticleT::kHarmonicPair:
			{
				double R0, K;
				in >> R0 >> K;
				fProperties[i] = new HarmonicPairT(R0, K);
				break;
			}
			case ParticleT::kLennardJonesPair:
			{
				double eps, sigma, cut_off;
				in >> eps >> sigma >> cut_off;
				fProperties[i] = new LennardJonesPairT(eps, sigma, cut_off);
				break;
			}
			default:
				ExceptionT::BadInputValue("ParticlePairT::ReadProperties", 
					"unrecognized property type: %d", property);
		}
	}

	/* echo particle properties */
	out << "\n Particle properties:\n\n";
	out << " Number of properties. . . . . . . . . . . . . . = " << fProperties.Length() << '\n';
	for (int i = 0; i < fProperties.Length(); i++)
	{
		out << " Property: " << i+1 << '\n';
		fProperties[i]->Write(out);
	}
}
