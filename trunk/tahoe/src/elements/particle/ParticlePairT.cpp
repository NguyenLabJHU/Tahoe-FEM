/* $Id: ParticlePairT.cpp,v 1.7 2002-11-26 22:56:26 paklein Exp $ */
#include "ParticlePairT.h"
#include "PairPropertyT.h"
#include "fstreamT.h"
#include "eControllerT.h"
#include "InverseMapT.h"

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
	fNeighborDistance(-1),
	fForce_man(0, fForce)
{
	/* input stream */
	ifstreamT& in = ElementSupport().Input();

	/* read parameters */
	in >> fNeighborDistance;

	/* checks */
	if (fNeighborDistance < kSmall) 
		ExceptionT::BadInputValue("ParticlePairT::ParticlePairT");
}

/* destructor */
ParticlePairT::~ParticlePairT(void)
{
	/* free properties list */
	for (int i = 0; i < fProperties.Length(); i++)
		delete fProperties[i];
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

void ParticlePairT::WriteOutput(void)
{
	/* inherited */
	ParticleT::WriteOutput();
	
	int ndof = NumDOF();
	int num_output = ndof + 2; /* displacement + PE + KE */

	/* output arrays length fGlobalTag */
	int num_particles = fGlobalTag.Length();
	dArray2DT n_values(num_particles, num_output), e_values;
	n_values = 0.0;

	/* global coordinates */
	const dArray2DT& coords = ElementSupport().CurrentCoordinates();

	/* pair properties function pointers */
	int current_property = -1;
	PairPropertyT::EnergyFunction energy_function = NULL;
	
	/* the field */
	const FieldT& field = Field();
	const dArray2DT& displacement = field[0];
	const dArray2DT* velocities = NULL;
	if (field.Order() > 0) velocities = &(field[1]);

	/* collect mass per particle */
	dArrayT mass(fNumTypes);
	for (int i = 0; i < fNumTypes; i++)
		mass[i] = fProperties[fPropertiesMap(i,i)]->Mass();
	
	/* run through neighbor list */
	iArrayT neighbors;
	dArrayT x_i, x_j, r_ij(ndof), vec, values_i;
	for (int i = 0; i < fNeighbors.MajorDim(); i++)
	{
		/* row of neighbor list */
		fNeighbors.RowAlias(i, neighbors);

		/* tags */
		int   tag_i = neighbors[0]; /* self is 1st spot */
		int  type_i = fType[tag_i];		
		int local_i = fGlobalToLocal.Map(tag_i);
		
		/* values for particle i */
		n_values.RowAlias(local_i, values_i);
		
		/* displacements */
		vec.Set(ndof, values_i.Pointer());
		displacement.RowCopy(tag_i, vec);
		
		/* kinetic energy */
		if (velocities)
		{
			velocities->RowAlias(tag_i, vec);
			values_i[ndof+1] = 0.5*mass[type_i]*dArrayT::Dot(vec, vec);
		}
		
		/* run though neighbors for one atom - first neighbor is self
		 * to compute potential energy */
		coords.RowAlias(tag_i, x_i);
		for (int j = 1; j < neighbors.Length(); j++)
		{
			/* tags */
			int   tag_j = neighbors[j];
			int  type_j = fType[tag_j];		
			int local_j = fGlobalToLocal.Map(tag_j);
			
			/* set pair property (if not already set) */
			int property = fPropertiesMap(type_i, type_j);
			if (property != current_property)
			{
				energy_function = fProperties[property]->getEnergyFunction();
				current_property = property;
			}
		
			/* global coordinates */
			coords.RowAlias(tag_j, x_j);
		
			/* connecting vector */
			r_ij.DiffOf(x_j, x_i);
			double r = r_ij.Magnitude();
			
			/* split interaction energy */
			double uby2 = 0.5*energy_function(r, NULL, NULL);
			values_i[ndof] += uby2;
			n_values(local_j, ndof) += uby2;
		}
	}	

	/* send */
	ElementSupport().WriteOutput(fOutputID, n_values, e_values);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* generate labels for output data */
void ParticlePairT::GenerateOutputLabels(ArrayT<StringT>& labels) const
{
	if (NumDOF() > 3) ExceptionT::GeneralFail("ParticlePairT::GenerateOutputLabels");

	/* displacement labels */
	const char* disp[3] = {"D_X", "D_Y", "D_Z"};
	
	int num_labels =
		NumDOF() // displacements
		+ 2;     // PE and KE

	labels.Dimension(num_labels);
	int dex = 0;
	for (dex = 0; dex < NumDOF(); dex++)
		labels[dex] = disp[dex];
	labels[dex++] = "PE";
	labels[dex++] = "KE";
}

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
	dArrayT x_i, x_j, r_ij(ndof);
	for (int i = 0; i < fNeighbors.MajorDim(); i++)
	{
		/* row of neighbor list */
		fNeighbors.RowAlias(i, neighbors);

		/* type */
		int  tag_i = neighbors[0]; /* self is 1st spot */
		int type_i = fType[tag_i];
	
		/* initialize */
		fForce_man.SetLength(fEqnos.MinorDim(i), false);
		fForce = 0.0;
		double* f_i = fForce.Pointer();
		double* f_j = fForce.Pointer(ndof);
		
		/* run though neighbors for one atom - first neighbor is self */
		fEqnos.RowAlias(i, eqs);
		coords.RowAlias(tag_i, x_i);
		for (int j = 1; j < neighbors.Length(); j++)
		{
			/* global tag */
			int  tag_j = neighbors[j];
			int type_j = fType[tag_j];
			
			/* set pair property (if not already set) */
			int property = fPropertiesMap(type_i, type_j);
			if (property != current_property)
			{
				force_function = fProperties[property]->getForceFunction();
				current_property = property;
			}
		
			/* global coordinates */
			coords.RowAlias(tag_j, x_j);
		
			/* connecting vector */
			r_ij.DiffOf(x_j, x_i);
			double r = r_ij.Magnitude();
			
			/* interaction force */
			double f = force_function(r, NULL, NULL);
			double fbyr = formKd*f/r;
			for (int k = 0; k < ndof; k++)
			{
				double f_k = r_ij[k]*fbyr;
				f_i[k] += f_k;
				*f_j++ +=-f_k;
			}
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

	/* verbose */
	if (ElementSupport().PrintInput())
	{
		ofstreamT& out = ElementSupport().Output();
		out << "\n Neighbor lists (self as leading neighbor):\n";
		out << " Neighbor cut-off distance . . . . . . . . . . . = " << fNeighborDistance << '\n';
		out << setw(kIntWidth) << "row" << "  n..." << '\n';
		iArrayT tmp(fNeighbors.Length(), fNeighbors.Pointer());
		tmp++;
		fNeighbors.WriteNumbered(out);
		tmp--;
		out.flush();
	}
}

/* construct the list of properties from the given input stream */
void ParticlePairT::EchoProperties(ifstreamT& in, ofstreamT& out)
{
	/* read potentials */
	int num_potentials = -1;
	in >> num_potentials;
	fProperties.Dimension(num_potentials);
	fProperties = NULL;
	for (int i = 0; i < fProperties.Length(); i++)
	{
		ParticleT::PropertyT property;
		in >> property;
		switch (property)
		{
			case ParticleT::kHarmonicPair:
			{
				double mass, R0, K;
				in >> mass >> R0 >> K;
				fProperties[i] = new HarmonicPairT(mass, R0, K);
				break;
			}
			case ParticleT::kLennardJonesPair:
			{
				double mass, eps, sigma, cut_off;
				in >> mass >> eps >> sigma >> cut_off;
				fProperties[i] = new LennardJonesPairT(mass, eps, sigma, cut_off);
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
