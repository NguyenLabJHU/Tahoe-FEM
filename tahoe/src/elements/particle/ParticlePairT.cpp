/* $Id: ParticlePairT.cpp,v 1.13 2003-01-27 07:00:26 paklein Exp $ */
#include "ParticlePairT.h"
#include "PairPropertyT.h"
#include "fstreamT.h"
#include "eControllerT.h"
#include "InverseMapT.h"
#include "CommManagerT.h"

/* pair property types */
#include "LennardJonesPairT.h"
#include "HarmonicPairT.h"
#include "ParadynPairT.h"

using namespace Tahoe;

/* parameters */
const int kMemoryHeadRoom = 15; /* percent */

/* constructor */
ParticlePairT::ParticlePairT(const ElementSupportT& support, const FieldT& field):
	ParticleT(support, field),
	fNeighbors(kMemoryHeadRoom),
	fEqnos(kMemoryHeadRoom),
	fForce_list_man(0, fForce_list)
{

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

	/* muli-processor information */
	CommManagerT& comm_manager = ElementSupport().CommManager();
	const ArrayT<int>* proc_map = comm_manager.ProcessorMap();
	int rank = ElementSupport().Rank();

	/* dimensions */
	int ndof = NumDOF();
	int num_output = ndof + 2; /* displacement + PE + KE */

	/* output arrays length number of active nodes */
	int num_particles = fNeighbors.MajorDim();
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
		mass[i] = fPairProperties[fPropertiesMap(i,i)]->Mass();

	/* map from partition node index */
	const InverseMapT* inverse_map = fCommManager.PartitionNodes_inv();
	
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
		int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;
		
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
			
			/* set pair property (if not already set) */
			int property = fPropertiesMap(type_i, type_j);
			if (property != current_property)
			{
				energy_function = fPairProperties[property]->getEnergyFunction();
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
			
			/* second node may not be on processor */
			if (!proc_map || (*proc_map)[tag_j] == rank) {
				int local_j = (inverse_map) ? inverse_map->Map(tag_j) : tag_j;
				
				if (local_j < 0 || local_j >= n_values.MajorDim())
					cout << " out of range: " << local_j << '\n';
				else
					n_values(local_j, ndof) += uby2;

			}
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
void ParticlePairT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
	/* time integration parameters */
	double constK = 0.0;
	double constM = 0.0;
	int formK = fController->FormK(constK);
	int formM = fController->FormM(constM);

	/* assemble particle mass */
	if (formM) {

		/* collect mass per particle */
		dArrayT mass(fNumTypes);
		for (int i = 0; i < fNumTypes; i++)
			mass[i] = fPairProperties[fPropertiesMap(i,i)]->Mass();
		mass *= constM;
	
		AssembleParticleMass(mass);
	}
	
	/* assemble stiffness */
	if (formK)
	{
		if (sys_type != GlobalT::kDiagonal && ElementSupport().IterationNumber() == 1)
			cout << "\n ParticlePairT::LHSDriver: WARNING: LHS matrix is a diagonal approximation" << endl;
	
		/* assembly information */
		const ElementSupportT& support = ElementSupport();
		int group = Group();
		int ndof = NumDOF();
	
		/* global coordinates */
		const dArray2DT& coords = support.CurrentCoordinates();

		/* pair properties function pointers */
		int current_property = -1;
		PairPropertyT::ForceFunction force_function = NULL;
		PairPropertyT::StiffnessFunction stiffness_function = NULL;

		/* run through neighbor list */
		fForce = 0.0;
		iArrayT neighbors;
		dArrayT x_i, x_j, r_ij(ndof);
		for (int i = 0; i < fNeighbors.MajorDim(); i++)
		{
			/* row of neighbor list */
			fNeighbors.RowAlias(i, neighbors);

			/* type */
			int  tag_i = neighbors[0]; /* self is 1st spot */
			int type_i = fType[tag_i];
			double* k_i = fForce(tag_i);
		
			/* run though neighbors for one atom - first neighbor is self */
			coords.RowAlias(tag_i, x_i);
			for (int j = 1; j < neighbors.Length(); j++)
			{
				/* global tag */
				int  tag_j = neighbors[j];
				int type_j = fType[tag_j];
				double* k_j = fForce(tag_j);
			
				/* set pair property (if not already set) */
				int property = fPropertiesMap(type_i, type_j);
				if (property != current_property)
				{
					force_function = fPairProperties[property]->getForceFunction();
					stiffness_function = fPairProperties[property]->getStiffnessFunction();
					current_property = property;
				}
		
				/* global coordinates */
				coords.RowAlias(tag_j, x_j);
		
				/* connecting vector */
				r_ij.DiffOf(x_j, x_i);
				double r = r_ij.Magnitude();
			
				/* interaction functions */
				double F = force_function(r, NULL, NULL);
				double K = stiffness_function(r, NULL, NULL);
				K = (K < 0.0) ? 0.0 : K;

				double Fbyr = F/r;
				for (int k = 0; k < ndof; k++)
				{
					double r_k = r_ij[k]*r_ij[k]/r/r;
					double K_k = constK*(K*r_k + Fbyr*(1.0 - r_k));
					k_i[k] += K_k;
					k_j[k] += K_k;
				}
			}
		}

		/* assemble */
		support.AssembleLHS(group, fForce, Field().Equations());
	}
}

/* form group contribution to the residual */
void ParticlePairT::RHSDriver(void)
{
	int nsd = NumSD();
	if (nsd == 3)
		RHSDriver3D();
	else if (nsd == 2)
		RHSDriver2D();
	else
		ExceptionT::GeneralFail("ParticlePairT::RHSDriver");
}

void ParticlePairT::RHSDriver2D(void)
{
	/* function name */
	const char caller[] = "ParticlePairT::RHSDriver2D";

	/* check 2D */
	if (NumDOF() != 2) ExceptionT::GeneralFail(caller, "2D only: %d", NumDOF());

	/* time integration parameters */
	double constMa = 0.0;
	double constKd = 0.0;
	int formMa = fController->FormMa(constMa);
	int formKd = fController->FormKd(constKd);

	//TEMP - interial force not implemented
	if (formMa) ExceptionT::GeneralFail(caller, "inertial force not implemented");

	/* assembly information */
	const ElementSupportT& support = ElementSupport();
	int group = Group();
	int ndof = NumDOF();
	
	/* global coordinates */
	const dArray2DT& coords = support.CurrentCoordinates();

	/* pair properties function pointers */
	int current_property = -1;
	PairPropertyT::ForceFunction force_function = NULL;
	const double* Paradyn_table = NULL;
	double dr = 1.0;
	int row_size = 0, num_rows = 0;

	/* run through neighbor list */
	fForce = 0.0;
	iArrayT neighbors;
	for (int i = 0; i < fNeighbors.MajorDim(); i++)
	{
		/* row of neighbor list */
		fNeighbors.RowAlias(i, neighbors);

		/* type */
		int   tag_i = neighbors[0]; /* self is 1st spot */
		int  type_i = fType[tag_i];
		double* f_i = fForce(tag_i);
		double* x_i = coords(tag_i);
		
		/* run though neighbors for one atom - first neighbor is self */
		for (int j = 1; j < neighbors.Length(); j++)
		{
			/* global tag */
			int   tag_j = neighbors[j];
			int  type_j = fType[tag_j];
			double* f_j = fForce(tag_j);
			double* x_j = coords(tag_j);

			/* set pair property (if not already set) */
			int property = fPropertiesMap(type_i, type_j);
			if (property != current_property)
			{
				if (!fPairProperties[property]->getParadynTable(&Paradyn_table, dr, row_size, num_rows))
					force_function = fPairProperties[property]->getForceFunction();
				current_property = property;
			}
		
			/* connecting vector */
			double r_ij_0 = x_j[0] - x_i[0];
			double r_ij_1 = x_j[1] - x_i[1];
			double r = sqrt(r_ij_0*r_ij_0 + r_ij_1*r_ij_1);
			
			/* interaction force */
			double F;
			if (Paradyn_table)
			{
				double pp = r*dr;
				int kk = int(pp);
				int max_row = num_rows-2;
				kk = (kk < max_row) ? kk : max_row;
				pp -= kk;
				pp = (pp < 1.0) ? pp : 1.0;				
				const double* c = Paradyn_table + kk*row_size;
				F = c[4] + pp*(c[5] + pp*c[6]);
			}
			else
				F = force_function(r, NULL, NULL);
			double Fbyr = formKd*F/r;

			r_ij_0 *= Fbyr;
			f_i[0] += r_ij_0;
			f_j[0] +=-r_ij_0;

			r_ij_1 *= Fbyr;
			f_i[1] += r_ij_1;
			f_j[1] +=-r_ij_1;
		}
	}

	/* assemble */
	support.AssembleRHS(group, fForce, Field().Equations());
}

void ParticlePairT::RHSDriver3D(void)
{
	/* function name */
	const char caller[] = "ParticlePairT::RHSDriver3D";

	/* check 3D */
	if (NumDOF() != 3) ExceptionT::GeneralFail(caller, "3D only: %d", NumDOF());

	/* time integration parameters */
	double constMa = 0.0;
	double constKd = 0.0;
	int formMa = fController->FormMa(constMa);
	int formKd = fController->FormKd(constKd);

	//TEMP - interial force not implemented
	if (formMa) ExceptionT::GeneralFail(caller, "inertial force not implemented");

	/* assembly information */
	const ElementSupportT& support = ElementSupport();
	int group = Group();
	int ndof = NumDOF();
	
	/* global coordinates */
	const dArray2DT& coords = support.CurrentCoordinates();

	/* pair properties function pointers */
	int current_property = -1;
	PairPropertyT::ForceFunction force_function = NULL;
	const double* Paradyn_table = NULL;
	double dr = 1.0;
	int row_size = 0, num_rows = 0;

	/* run through neighbor list */
	fForce = 0.0;
	iArrayT neighbors;
	for (int i = 0; i < fNeighbors.MajorDim(); i++)
	{
		/* row of neighbor list */
		fNeighbors.RowAlias(i, neighbors);

		/* type */
		int   tag_i = neighbors[0]; /* self is 1st spot */
		int  type_i = fType[tag_i];
		double* f_i = fForce(tag_i);
		double* x_i = coords(tag_i);
		
		/* run though neighbors for one atom - first neighbor is self */
		for (int j = 1; j < neighbors.Length(); j++)
		{
			/* global tag */
			int   tag_j = neighbors[j];
			int  type_j = fType[tag_j];
			double* f_j = fForce(tag_j);
			double* x_j = coords(tag_j);

			/* set pair property (if not already set) */
			int property = fPropertiesMap(type_i, type_j);
			if (property != current_property)
			{
				if (!fPairProperties[property]->getParadynTable(&Paradyn_table, dr, row_size, num_rows))
					force_function = fPairProperties[property]->getForceFunction();
				current_property = property;
			}
		
			/* connecting vector */
			double r_ij_0 = x_j[0] - x_i[0];
			double r_ij_1 = x_j[1] - x_i[1];
			double r_ij_2 = x_j[2] - x_i[2];
			double r = sqrt(r_ij_0*r_ij_0 + r_ij_1*r_ij_1 + r_ij_2*r_ij_2);
			
			/* interaction force */
			double F;
			if (Paradyn_table)
			{
				double pp = r*dr;
				int kk = int(pp);
				int max_row = num_rows-2;
				kk = (kk < max_row) ? kk : max_row;
				pp -= kk;
				pp = (pp < 1.0) ? pp : 1.0;				
				const double* c = Paradyn_table + kk*row_size;
				F = c[4] + pp*(c[5] + pp*c[6]);
			}
			else
				F = force_function(r, NULL, NULL);
			double Fbyr = formKd*F/r;

			r_ij_0 *= Fbyr;
			f_i[0] += r_ij_0;
			f_j[0] +=-r_ij_0;

			r_ij_1 *= Fbyr;
			f_i[1] += r_ij_1;
			f_j[1] +=-r_ij_1;

			r_ij_2 *= Fbyr;
			f_i[2] += r_ij_2;
			f_j[2] +=-r_ij_2;
		}
	}

	/* assemble */
	support.AssembleRHS(group, fForce, Field().Equations());
}

/* set neighborlists */
void ParticlePairT::SetConfiguration(void)
{
	/* inherited */
	ParticleT::SetConfiguration();

	/* reset neighbor lists */
	const ArrayT<int>* part_nodes = fCommManager.PartitionNodes();
	GenerateNeighborList(part_nodes, fNeighborDistance, fNeighbors, false, true);
	
	ofstreamT& out = ElementSupport().Output();
	out << "\n Neighbor statistics:\n";
	out << " Total number of neighbors . . . . . . . . . . . = " << fNeighbors.Length() << '\n';
	out << " Minimum number of neighbors . . . . . . . . . . = " << fNeighbors.MinMinorDim(0) << '\n';
	out << " Maximum number of neighbors . . . . . . . . . . = " << fNeighbors.MaxMinorDim() << '\n';
	if (fNeighbors.MajorDim() > 0)
	out << " Average number of neighbors . . . . . . . . . . = " << double(fNeighbors.Length())/fNeighbors.MajorDim() << '\n';
	else
	out << " Average number of neighbors . . . . . . . . . . = " << 0 << '\n';

	/* verbose */
	if (ElementSupport().PrintInput())
	{
		out << " Neighbor lists (self as leading neighbor):\n";
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
	fPairProperties.Dimension(num_potentials);
	fPairProperties = NULL;
	for (int i = 0; i < fPairProperties.Length(); i++)
	{
		ParticleT::PropertyT property;
		in >> property;
		switch (property)
		{
			case ParticleT::kHarmonicPair:
			{
				double mass, R0, K;
				in >> mass >> R0 >> K;
				fPairProperties[i] = new HarmonicPairT(mass, R0, K);
				break;
			}
			case ParticleT::kLennardJonesPair:
			{
				double mass, eps, sigma, alpha;
				in >> mass >> eps >> sigma >> alpha;
				fPairProperties[i] = new LennardJonesPairT(mass, eps, sigma, alpha);
				break;
			}
			case ParticleT::kParadynPair:
			{
				StringT file;
				in >> file;
				file.ToNativePathName();

				StringT path;
				path.FilePath(in.filename());				
				file.Prepend(path);
			
				fPairProperties[i] = new ParadynPairT(file);
				break;
			}
			default:
				ExceptionT::BadInputValue("ParticlePairT::ReadProperties", 
					"unrecognized property type: %d", property);
		}
	}

	/* echo particle properties */
	out << "\n Particle properties:\n\n";
	out << " Number of properties. . . . . . . . . . . . . . = " << fPairProperties.Length() << '\n';
	for (int i = 0; i < fPairProperties.Length(); i++)
	{
		out << " Property: " << i+1 << '\n';
		fPairProperties[i]->Write(out);
	}
	
	/* copy into base class list */
	fParticleProperties.Dimension(fPairProperties.Length());
	for (int i = 0; i < fPairProperties.Length(); i++)
		fParticleProperties[i] = fPairProperties[i];
}
