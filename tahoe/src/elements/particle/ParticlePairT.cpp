/* $Id: ParticlePairT.cpp,v 1.44.6.1 2005-05-27 19:55:17 paklein Exp $ */
#include "ParticlePairT.h"

#include "PairPropertyT.h"
#include "eIntegratorT.h"
#include "InverseMapT.h"
#include "CommManagerT.h"
#include "dSPMatrixT.h"
#include "ModelManagerT.h"
#include "dSymMatrixT.h"
#include "dArray2DT.h"
#include "iGridManagerT.h"
#include "ParameterContainerT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"

#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>

/* pair property types */
#include "LennardJonesPairT.h"
#include "HarmonicPairT.h"
#include "ParadynPairT.h"
#include "MatsuiPairT.h"

using namespace Tahoe;

/* parameters */
const int kMemoryHeadRoom = 15; /* percent */
const int kNumOutput = 8;
static const char* OutputNames[kNumOutput] = {
	"displacement",
	"potential_energy",
	"kinetic_energy",
	"stress",
	"strain",
	"slip_vector",	
	"centrosymmetry",
	"coordination_number"
};

/* constructor */
ParticlePairT::ParticlePairT(const ElementSupportT& support):
	ParticleT(support),
	fNeighbors(kMemoryHeadRoom),
	fNearestNeighbors(kMemoryHeadRoom),
	fRefNearestNeighbors(kMemoryHeadRoom),
	fEqnos(kMemoryHeadRoom),
	fForce_list_man(0, fForce_list)
{
	SetName("particle_pair");
	fopen = false;
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
	const char caller[] = "ParticlePairT::WriteOutput";

	/* inherited */
	ParticleT::WriteOutput();

	/* muli-processor information */
	CommManagerT& comm_manager = ElementSupport().CommManager();
	const ArrayT<int>* proc_map = comm_manager.ProcessorMap();
	int rank = ElementSupport().Rank();

	/* dimensions */
	int ndof = NumDOF();
	int nstrs = dSymMatrixT::NumValues(ndof);

	/* number of output variables */
	iArrayT counts;
	SetOutputCount(fOutputFlags, counts);
	int num_output = counts.Sum();

	/* offsets to the different output values */
	iArrayT offsets(fOutputFlags.Length());
	offsets = 0;
	for (int i = 1; i < offsets.Length(); i++)
		offsets[i] = offsets[i-1] + counts[i-1];

	/* number of nodes */
	const ArrayT<int>* parition_nodes = comm_manager.PartitionNodes();
	int non = (parition_nodes) ? 
		parition_nodes->Length() : 
		ElementSupport().NumNodes();

	/* map from partition node index */
	const InverseMapT* inverse_map = comm_manager.PartitionNodes_inv();

	/* output arrays length number of active nodes */
	dArray2DT n_values(non, num_output), e_values;
	n_values = 0.0;

	/* global coordinates */
	const dArray2DT& coords = ElementSupport().CurrentCoordinates();

	/* pair properties function pointers */
	int current_property = -1;
	PairPropertyT::EnergyFunction energy_function = NULL;	
	PairPropertyT::ForceFunction force_function = NULL;

	/* the field */
	const FieldT& field = Field();
	const dArray2DT& displacement = field[0];
	const dArray2DT* velocities = NULL;
	if (field.Order() > 0) velocities = &(field[1]);

	/* atomic volume */
  	double V0 = 0.0;
  	if (NumSD() == 1)
    	V0 = fLatticeParameter;
	else if (NumSD() == 2)
		V0 = sqrt(3.0)*fLatticeParameter*fLatticeParameter/2.0; /* 2D hex */
	else /* 3D */
		V0 = fLatticeParameter*fLatticeParameter*fLatticeParameter/4.0; /* FCC */

	/* collect mass per particle */
	int num_types = fTypeNames.Length();
	dArrayT mass(num_types);
	for (int i = 0; i < num_types; i++)
		mass[i] = fPairProperties[fPropertiesMap(i,i)]->Mass();

	/* collect displacements and kinetic contribution to the virial */
	dArrayT vec, values_i;
	dSymMatrixT temp(ndof);
	for (int i = 0; i < non; i++)
	{
		/* particle ID */
		int   tag_i = (parition_nodes) ? (*parition_nodes)[i] : i;
		int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;
		int  type_i = fType[tag_i];

		/* values for particle i */
		n_values.RowAlias(local_i, values_i);

		/* displacements */
		if (fOutputFlags[kDisplacement]) {
			vec.Set(ndof, values_i.Pointer(offsets[kDisplacement]));
			displacement.RowCopy(tag_i, vec);
		}

		/* kinetic contribution to the virial */
		if (velocities && fOutputFlags[kStress]) {
			velocities->RowAlias(tag_i, vec);
			temp.Outer(vec);
			int index = offsets[kStress];
		 	for (int cc = 0; cc < nstrs; cc++)
		   		values_i[index++] = (fabs(V0) > kSmall) ? -mass[type_i]*temp[cc]/V0 : 0.0;
		}
	}
	
	/* run through neighbor list */
	iArrayT neighbors;
	dArrayT x_i, x_j, r_ij(ndof);
	dSymMatrixT vs_i(ndof);
	dArrayT SlipVector(ndof);
	dMatrixT Strain(ndof);
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

		/* initialize stress */
		if (fOutputFlags[kStress]) vs_i = 0.0;

		/* kinetic energy */
		if (velocities && fOutputFlags[kKE]) {
			velocities->RowAlias(tag_i, vec);
			values_i[offsets[kKE]] = 0.5*mass[type_i]*dArrayT::Dot(vec, vec);
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
			if (property != current_property) {
				energy_function = fPairProperties[property]->getEnergyFunction();
				force_function = fPairProperties[property]->getForceFunction();
				current_property = property;
			}
		
			/* global coordinates */
			coords.RowAlias(tag_j, x_j);

			/* connecting vector */
			r_ij.DiffOf(x_j, x_i);
			double r = r_ij.Magnitude();

			/* split interaction energy */
			double uby2 = 0.0;
			if (fOutputFlags[kPE]) {
				uby2 = 0.5*energy_function(r, NULL, NULL); 
				values_i[offsets[kPE]] += uby2;
			}

			/* stress contribution */			
			double Fbyr = 0.0;
			if (fOutputFlags[kStress])
			{
	      		/* interaction force */
				double F = force_function(r, NULL, NULL);
				Fbyr = F/r;

				temp.Outer(r_ij);
				vs_i.AddScaled(0.5*Fbyr, temp);
			}

			/* second node may not be on processor */
			if (!proc_map || (*proc_map)[tag_j] == rank) 
			{
				int local_j = (inverse_map) ? inverse_map->Map(tag_j) : tag_j;	
				if (local_j < 0 || local_j >= n_values.MajorDim())
					cout << caller << ": out of range: " << local_j << '\n';
				else {

					/* potential energy */
					if (fOutputFlags[kPE])
						n_values(local_j, offsets[kPE]) += uby2;
					
			 		/* accumulate into stress into array */
					if (fOutputFlags[kStress]) {
						int index = offsets[kStress];
		 				for (int cc = 0; cc < nstrs; cc++)
		   					n_values(local_j, index++) += (fabs(V0) > kSmall) ? 0.5*Fbyr*temp[cc]/V0 : 0.0;
		 			}
				}
			}
		}

		/* copy stress into array */
		if (fOutputFlags[kStress]) {
			int index = offsets[kStress];
			for (int cc = 0; cc < nstrs; cc++)
				values_i[index++] += (fabs(V0) > kSmall) ? vs_i[cc]/V0 : 0.0;
		}
	}

	/* calculate centrosymmetry parameter */
	dArrayT csp;
	if (fOutputFlags[kCS]) {
		csp.Dimension(non);
		Calc_CSP(fNearestNeighbors, csp);
	}
	
	// calculate coordination number
	iArrayT cnarray;
	if (fOutputFlags[kCN]) {
		cnarray.Dimension(non);
		Calc_CN(fNearestNeighbors, cnarray);
	}

	/* slip vector, stress, and strain */
	dArray2DT s_values;
	if (fOutputFlags[kSlipVector] || fOutputFlags[kStress] || fOutputFlags[kStrain]) 
	{
		/* work space */
		int num_s_vals = nstrs + 1 + ndof; /* {strain, J, slip vector} */
		s_values.Dimension(non, num_s_vals);
		s_values = 0.0;

	    /* flag for specifying Lagrangian (0) or Eulerian (1) strain */ 
	    const int kEulerLagr = 0;

		/* calculate slip vector and strain */
		Calc_Slip_and_Strain(s_values, fRefNearestNeighbors, kEulerLagr);
	}

	// combine strain, slip vector, centrosymmetry parameter and coordination number into n_values list
	for (int i = 0; i < fNeighbors.MajorDim(); i++)
	{
		/* row of neighbor list */
		fNeighbors.RowAlias(i, neighbors);

		/* tags */
		int   tag_i = neighbors[0]; /* self is 1st spot */
		int  type_i = fType[tag_i];
		int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;

		/* strain */
		if (fOutputFlags[kStrain])
			for (int is = 0; is < nstrs; is++)
				n_values(local_i, offsets[kStrain]+is) = s_values(local_i,is);

		/* stress */
		if (fOutputFlags[kStress])
		{
			/* recover J, the determinant of the deformation gradient, for atom i
			 * and divide stress values by it */
			double J = s_values(local_i,nstrs);
	    	for (int is = 0; is < nstrs; is++)
	    		if (fabs(J) > kSmall)
	    			n_values(local_i, offsets[kStress]+is) /= J;
	    		else
	    			n_values(local_i, offsets[kStress]+is) = 0.0;
	    }

		/* slip vector */
		if (fOutputFlags[kSlipVector])
			for (int n = 0; n < ndof; n++)
				n_values(local_i, offsets[kSlipVector]+n) = s_values(local_i, nstrs+1+n);

		/* centrosymmetry */
		if (fOutputFlags[kCS])
			n_values(local_i, offsets[kCS]) = csp[local_i];
		
		// coordination number
		if (fOutputFlags[kCN])
			n_values(local_i, offsets[kCN]) = cnarray[local_i];
	}

#if 0
	/* temporary to calculate crack propagation velocity */
	ifstreamT& in = ElementSupport().Input();
	ModelManagerT& model = ElementSupport().Model();
	const ArrayT<StringT> id_list = model.NodeSetIDs();
	iArrayT nodelist;
	dArray2DT partial;
	nodelist = model.NodeSet(id_list[id_list.Length()-1]); // want last nodeset
	const StringT& input_file = in.filename();
	fsummary_file.Root(input_file);
	fsummary_file.Append(".crack");
	double xcoord = coords(nodelist[0],0);
	double ydispcrit = .13;
	const double& time = ElementSupport().Time();
	for (int i = 0; i < nodelist.Length(); i++)
	{
	    int node = nodelist[i];
	    if (fabs(displacement(node,1)) >= ydispcrit)
	      xcoord = coords(node,0);
	}

	if (fopen)
	{
		fout.open_append(fsummary_file);
	  	fout.precision(13);
	  	fout << xcoord 
	  	     << setw(25) << time
	  	     << endl;
	}
	else
	{
	  	fout.open(fsummary_file);
		fopen = true;
	  	fout.precision(13);
	  	fout << "x-coordinate"
	  	     << setw(25) << "Time"
	  	     << endl;
	  	fout << xcoord 
	  	     << setw(25) << time
	  	     << endl;
	}
#endif

#if 0
	/* Temporary to calculate MD energy history and write to file */
	ifstreamT& in = ElementSupport().Input();
	ModelManagerT& model = ElementSupport().Model();
	const ArrayT<StringT> id_list = model.NodeSetIDs();
	iArrayT nodelist;
	dArray2DT partial;
	nodelist = model.NodeSet(id_list[id_list.Length()-1]);
	nodelist = model.NodeSet(id_list[3]);  // id_list[3]
	partial.Dimension(nodelist.Length(), n_values.MinorDim());
	partial.RowCollect(nodelist, n_values);
	const StringT& input_file = in.filename();
	fsummary_file.Root(input_file);
	fsummary_file2.Root(input_file);
	fsummary_file.Append(".sum");
	fsummary_file2.Append(".full");
	if (fopen)
	{
	        fout.open_append(fsummary_file);
			fout2.open_append(fsummary_file2);
			fout.precision(13);
			fout2.precision(13);
			fout << n_values.ColumnSum(3) 
				 << setw(25) << n_values.ColumnSum(2)
				 << setw(25) << n_values.ColumnSum(3) + n_values.ColumnSum(2)
				 << endl;
			fout2 << partial.ColumnSum(3) 
				 << setw(25) << partial.ColumnSum(2)
				 << setw(25) << partial.ColumnSum(3) + partial.ColumnSum(2)
				 << endl;
	}
	else
	{
			fout.open(fsummary_file);
			fout2.open(fsummary_file2);
			fopen = true;
			fout.precision(13);
			fout2.precision(13);
			fout << "Kinetic Energy"
				 << setw(25) << "Potential Energy"
				 << setw(25) << "Total Energy"
				 << endl;
			fout << n_values.ColumnSum(3) 
				 << setw(25) << n_values.ColumnSum(2)
				 << setw(25) << n_values.ColumnSum(3) + n_values.ColumnSum(2)
				 << endl;
			fout2 << "Kinetic Energy"
				 << setw(25) << "Potential Energy"
				 << setw(25) << "Total Energy"
				 << endl;
			fout2 << partial.ColumnSum(3) 
				 << setw(25) << partial.ColumnSum(2)
				 << setw(25) << partial.ColumnSum(3) + partial.ColumnSum(2)
				 << endl;
	}
#endif

	/* send */
	ElementSupport().WriteOutput(fOutputID, n_values, e_values);
}

/* write restart data to the output stream */
void ParticlePairT::WriteRestart(ofstreamT& out) const
{
	/* inherited */
	ParticleT::WriteRestart(out);

	/* open stream */	
	StringT file = out.filename();
	file.Append(".neighbor");
	ofstreamT my_out(file);
	if (!my_out.is_open())
		ExceptionT::GeneralFail("ParticlePairT::WriteRestart",
			"could not open file \"%s\"", file.Pointer());

	/* write neighbor lists */
	fNeighbors.Write(my_out); my_out << '\n';
	fNearestNeighbors.Write(my_out); my_out << '\n';
	fRefNearestNeighbors.Write(my_out); my_out << '\n';
}

/* read restart data to the output stream */
void ParticlePairT::ReadRestart(ifstreamT& in)
{
	/* inherited */
	ParticleT::ReadRestart(in);

	/* open stream */	
	StringT file = in.filename();
	file.Append(".neighbor");
	ifstreamT my_in(file);
	if (!my_in.is_open())
		ExceptionT::GeneralFail("ParticlePairT::ReadRestart",
			"could not open file \"%s\"", file.Pointer());

	/* read neighbor lists */
	fNeighbors.Read(my_in);
	fNearestNeighbors.Read(my_in);
	fRefNearestNeighbors.Read(my_in);
}

/* compute the part of the stiffness matrix */
void ParticlePairT::FormStiffness(const InverseMapT& col_to_col_eq_row_map,
	const iArray2DT& col_eq, dSPMatrixT& stiffness)
{
	const char caller[] = "ParticlePairT::FormStiffness";

	/* map should return -1 of out of range */
	if (col_to_col_eq_row_map.OutOfRange() != InverseMapT::MinusOne)
		ExceptionT::GeneralFail(caller, "inverse map out of range should return -1");

	/* assembly information */
	const ElementSupportT& support = ElementSupport();
	int group = Group();
	int ndof = NumDOF();
	fLHS.Dimension(2*ndof);
		
	/* global coordinates */
	const dArray2DT& coords = support.CurrentCoordinates();

	/* pair properties function pointers */
	int current_property = -1;
	PairPropertyT::ForceFunction force_function = NULL;
	PairPropertyT::StiffnessFunction stiffness_function = NULL;

	/* work space */
	dArrayT r_ij(NumDOF(), fRHS.Pointer());
	dArrayT r_ji(NumDOF(), fRHS.Pointer() + NumDOF());

	/* run through neighbor list */
	const iArray2DT& field_eqnos = Field().Equations();
	iArrayT row_eqnos, col_eqnos; 
	iArrayT neighbors;
	dArrayT x_i, x_j;
	for (int i = 0; i < fNeighbors.MajorDim(); i++)
	{
		/* row of neighbor list */
		fNeighbors.RowAlias(i, neighbors);

		/* type */
		int  tag_i = neighbors[0]; /* self is 1st spot */
		int type_i = fType[tag_i];
		
		/* particle equations */
		field_eqnos.RowAlias(tag_i, row_eqnos);

		/* run though neighbors for one atom - first neighbor is self */
		coords.RowAlias(tag_i, x_i);
		for (int j = 1; j < neighbors.Length(); j++)
		{
			/* global tag */
			int tag_j = neighbors[j];
			
			/* particle is a target column */
			int col_eq_index = col_to_col_eq_row_map.Map(tag_j);
			if (col_eq_index != -1)
			{
				/* more particle info */
				int type_j = fType[tag_j];

				/* particle equations */
				col_eq.RowAlias(col_eq_index, col_eqnos);

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
				r_ji.SetToScaled(-1.0, r_ij);

				/* interaction functions */
				double F = force_function(r, NULL, NULL);
				double K = stiffness_function(r, NULL, NULL);
				double Fbyr = F/r;

				/* 1st term */
				fLHS.Outer(fRHS, fRHS, (K - Fbyr)/r/r);

				/* 2nd term */				
				fLHS.AddScaled(Fbyr, fOneOne);

				/* assemble */
				for (int p = 0; p < row_eqnos.Length(); p++)
					for (int q = 0; q < col_eqnos.Length(); q++)
						stiffness.AddElement(row_eqnos[p]-1, col_eqnos[q]-1, fLHS(p,q));
			}
		}
	}
}

/* information about subordinate parameter lists */
void ParticlePairT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParticleT::DefineSubs(sub_list);

	/* interactions */
	sub_list.AddSub("pair_particle_interaction", ParameterListT::OnePlus);
	
	/* output */
	sub_list.AddSub("particle_pair_output", ParameterListT::ZeroOrOnce);
}

/* return the description of the given inline subordinate parameter list */
void ParticlePairT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "pair_property_choice")
	{
		order = ParameterListT::Choice;
		
		/* harmonic pair potential */
		sub_lists.AddSub("harmonic");

		/* Lennard-Jones 6/12 */
		sub_lists.AddSub("Lennard_Jones");

		/* Paradyn pair potential */
		sub_lists.AddSub("Paradyn_pair");

		/* Matsui pair potential */
		sub_lists.AddSub("Matsui");
	}
	else /* inherited */
		ParticleT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ParticlePairT::NewSub(const StringT& name) const
{
	/* try to construct potential */
	PairPropertyT* pair_property = PairPropertyT::New(name, &(ElementSupport()));
	if (pair_property)
		return pair_property;
	else if (name == "pair_particle_interaction")
	{
		ParameterContainerT* interactions = new ParameterContainerT(name);
		interactions->SetSubSource(this);

		/* particle type labels */
		interactions->AddParameter(ParameterT::Word, "label_1");
		interactions->AddParameter(ParameterT::Word, "label_2");
	
		/* properties choice list */
		interactions->AddSub("pair_property_choice", ParameterListT::Once, true);

		return interactions;
	}
	else if (name == "particle_pair_output")
	{
		ParameterContainerT* output = new ParameterContainerT(name);
		
		/* all true by default */
		for (int i = 0; i < kNumOutput; i++) {
			ParameterT var(ParameterT::Integer, OutputNames[i]);
			var.SetDefault(1);
			output->AddParameter(var, ParameterListT::ZeroOrOnce);
		}

		return output;
	}
	else /* inherited */
		return ParticleT::NewSub(name);
}

/* accept parameter list */
void ParticlePairT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParticleT::TakeParameterList(list);

	/* output variables */
	fOutputFlags.Dimension(kNumOutput);
	fOutputFlags = 0;
	const ParameterListT* output = list.List("particle_pair_output");
	if (output) 
	{
		/* set flags */
		for (int i = 0; i < kNumOutput; i++)
		{
			/* look for entry */
			const ParameterT* value = output->Parameter(OutputNames[i]);
			if (value) {
				int do_write = *value;
				if (do_write)
					fOutputFlags[i] = 1;
			}
		}
	}

	/* dimension */
	int ndof = NumDOF();
	fLHS.Dimension(2*ndof);
	fRHS.Dimension(2*ndof);

	/* constant matrix needed to calculate stiffness */
	fOneOne.Dimension(fLHS);
	dMatrixT one(ndof);
	one.Identity();
	fOneOne.SetBlock(0, 0, one);
	fOneOne.SetBlock(ndof, ndof, one);
	one *= -1;
	fOneOne.SetBlock(0, ndof, one);
	fOneOne.SetBlock(ndof, 0, one);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* return number of values for each output variable */
void ParticlePairT::SetOutputCount(const iArrayT& flags, iArrayT& counts) const
{
	/* dimension check */
	if (flags.Length() != kNumOutput)
		ExceptionT::SizeMismatch("ParticlePairT::SetOutputCount");
	
	/* initialize */
	counts.Dimension(flags.Length());
	counts = 0;

	/* set output flags */
	if (flags[kDisplacement]) counts[kDisplacement] = NumDOF();
	if (flags[kPE]) counts[kPE] = 1;
	if (flags[kKE]) counts[kKE] = 1;
	if (flags[kCS]) counts[kCS] = 1;
	if (flags[kCN]) counts[kCN] = 1;
	if (flags[kStress]) counts[kStress] = dSymMatrixT::NumValues(NumSD());
	if (flags[kStrain]) counts[kStrain] = dSymMatrixT::NumValues(NumSD());
	if (flags[kSlipVector]) counts[kSlipVector] = NumDOF();
}

/* generate labels for output data */
void ParticlePairT::GenerateOutputLabels(ArrayT<StringT>& labels) const
{
	const char caller[] = "ParticlePairT::GenerateOutputLabels";
	int ndof = NumDOF();
	if (ndof > 3) ExceptionT::GeneralFail(caller);

	/* number of output variables */
	iArrayT counts;
	SetOutputCount(fOutputFlags, counts);
	int num_output = counts.Sum();

	/* offsets to the different output values */
	iArrayT offsets(fOutputFlags.Length());
	offsets = 0;
	for (int i = 1; i < offsets.Length(); i++)
		offsets[i] = offsets[i-1] + counts[i-1];

	/* initialize */
	labels.Dimension(num_output);

	/* displacement labels */
	if (fOutputFlags[kDisplacement]) {
		const char* disp[3] = {"D_X", "D_Y", "D_Z"};
		int index = offsets[kDisplacement];
		for (int i = 0; i < ndof; i++)
			labels[index++] = disp[i];
	}

	/* potential energy */
	if (fOutputFlags[kPE])
		labels[offsets[kPE]] = "PE";

	/* kinetic energy */
	if (fOutputFlags[kKE])
		labels[offsets[kKE]] = "KE";
	
	/* centrosymmetry */
	if (fOutputFlags[kCS])
		labels[offsets[kCS]] = "CS";
	
	// coordination number
	if (fOutputFlags[kCN])
		labels[offsets[kCN]] = "CN";
	
	/* slip vector */
	if (fOutputFlags[kSlipVector]) {
		const char* SV[3] = {"SV_X", "SV_Y", "SV_Z"};
		int index = offsets[kSlipVector];
		for (int i = 0; i < ndof; i++)
			labels[index++] = SV[i];
	}
	
	/* stress */
	if (fOutputFlags[kStress]) {
		const char* s1D[1] = {"s11"};
		const char* s2D[3] = {"s11", "s22", "s12"};
		const char* s3D[6] = {"s11", "s22", "s33", "s23", "s13", "s12"};

		const char** slabels = NULL;
		if (ndof == 1) slabels = s1D;
		else if (ndof == 2) slabels = s2D;
		else if (ndof == 3) slabels = s3D;
		else ExceptionT::GeneralFail(caller);	
		
		int nstrs = dSymMatrixT::NumValues(ndof);
		int index = offsets[kStress];
		for (int i = 0; i < nstrs; i++)
			labels[index++] = slabels[i];
	}
	
	/* strain */
	if (fOutputFlags[kStrain]) {		
		const char* e1D[1] = {"e11"};
		const char* e2D[3] = {"e11", "e12", "e22"};
		const char* e3D[6] = {"e11", "e12", "e13", "e22", "e23", "e33"};

		const char** elabels = NULL;
		if (ndof == 1) elabels = e1D;
		else if (ndof == 2) elabels = e2D;
		else if (ndof == 3) elabels = e3D;
		else ExceptionT::GeneralFail(caller);	
		
		int nstrn = dSymMatrixT::NumValues(ndof);
		int index = offsets[kStrain];
		for (int i = 0; i < nstrn; i++)
			labels[index++] = elabels[i];
	}
}

/* form group contribution to the stiffness matrix */
void ParticlePairT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
	/* time integration parameters */
	double constK = 0.0;
	double constM = 0.0;
	int formK = fIntegrator->FormK(constK);
	int formM = fIntegrator->FormM(constM);

	/* assemble particle mass */
	if (formM) {

		/* collect mass per particle */
		int num_types = fTypeNames.Length();
		dArrayT mass(num_types);
		for (int i = 0; i < num_types; i++)
			mass[i] = fPairProperties[fPropertiesMap(i,i)]->Mass();
		mass *= constM;
	
		AssembleParticleMass(mass);
	}
	
	/* assemble diagonal stiffness */
	if (formK && sys_type == GlobalT::kDiagonal)
	{
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
	else if (formK)
	{
		/* assembly information */
		const ElementSupportT& support = ElementSupport();
		int group = Group();
		int ndof = NumDOF();
		fLHS.Dimension(2*ndof);
		
		/* global coordinates */
		const dArray2DT& coords = support.CurrentCoordinates();

		/* pair properties function pointers */
		int current_property = -1;
		PairPropertyT::ForceFunction force_function = NULL;
		PairPropertyT::StiffnessFunction stiffness_function = NULL;

		/* work space */
		dArrayT r_ij(NumDOF(), fRHS.Pointer());
		dArrayT r_ji(NumDOF(), fRHS.Pointer() + NumDOF());

		/* run through neighbor list */
		const iArray2DT& field_eqnos = Field().Equations();
		iArray2DT pair_eqnos(2, ndof); 
		iArrayT pair(2);
		iArrayT neighbors;
		dArrayT x_i, x_j;
		for (int i = 0; i < fNeighbors.MajorDim(); i++)
		{
			/* row of neighbor list */
			fNeighbors.RowAlias(i, neighbors);

			/* type */
			int  tag_i = neighbors[0]; /* self is 1st spot */
			int type_i = fType[tag_i];
			pair[0] = tag_i;
		
			/* run though neighbors for one atom - first neighbor is self */
			coords.RowAlias(tag_i, x_i);
			for (int j = 1; j < neighbors.Length(); j++)
			{
				/* global tag */
				int  tag_j = neighbors[j];
				int type_j = fType[tag_j];
				pair[1] = tag_j;
			
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
				r_ji.SetToScaled(-1.0, r_ij);
			
				/* interaction functions */
				double F = constK*force_function(r, NULL, NULL);
				double K = constK*stiffness_function(r, NULL, NULL);
				double Fbyr = F/r;

				/* 1st term */
				fLHS.Outer(fRHS, fRHS, (K - Fbyr)/r/r);
		
				/* 2nd term */
				fLHS.AddScaled(Fbyr, fOneOne);
				
				/* assemble */
				pair_eqnos.RowCollect(pair, field_eqnos);
				support.AssembleLHS(group, fLHS, pair_eqnos);
			}
		}
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
	else if (nsd == 1)
		RHSDriver1D();	
	else
		ExceptionT::GeneralFail("ParticlePairT::RHSDriver", "unsupported dimension %d", nsd);
		
	ApplyDamping(fNeighbors);
	
	/* assemble */
	ElementSupport().AssembleRHS(Group(), fForce, Field().Equations());
}

void ParticlePairT::RHSDriver1D(void)
{
	/* function name */
	const char caller[] = "ParticlePairT::RHSDriver1D";

	/* check 1D */
	if (NumDOF() != 1) ExceptionT::GeneralFail(caller, "1D only: %d", NumDOF());

	/* time integration parameters */
	double constMa = 0.0;
	double constKd = 0.0;
	int formMa = fIntegrator->FormMa(constMa);
	int formKd = fIntegrator->FormKd(constKd);

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
		const double* x_i = coords(tag_i);
		
		/* run though neighbors for one atom - first neighbor is self */
		for (int j = 1; j < neighbors.Length(); j++)
		{
			/* global tag */
			int   tag_j = neighbors[j];
			int  type_j = fType[tag_j];
			double* f_j = fForce(tag_j);
			const double* x_j = coords(tag_j);

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
			double r = sqrt(r_ij_0*r_ij_0);
			
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
		}
	}
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
	int formMa = fIntegrator->FormMa(constMa);
	int formKd = fIntegrator->FormKd(constKd);

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
		const double* x_i = coords(tag_i);
		
		/* run though neighbors for one atom - first neighbor is self */
		for (int j = 1; j < neighbors.Length(); j++)
		{
			/* global tag */
			int   tag_j = neighbors[j];
			int  type_j = fType[tag_j];
			double* f_j = fForce(tag_j);
			const double* x_j = coords(tag_j);

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

}

void ParticlePairT::RHSDriver3D(void)
{
	/* function name */
	const char caller[] = "ParticlePairT::RHSDriver3D";

// DEBUG
//CommunicatorT& fComm = ElementSupport().Communicator();
//fComm.Log(CommunicatorT::kUrgent, caller);


	/* check 3D */
	if (NumDOF() != 3) ExceptionT::GeneralFail(caller, "3D only: %d", NumDOF());

	/* time integration parameters */
	double constMa = 0.0;
	double constKd = 0.0;
	int formMa = fIntegrator->FormMa(constMa);
	int formKd = fIntegrator->FormKd(constKd);

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
		const double* x_i = coords(tag_i);
		
		/* run though neighbors for one atom - first neighbor is self */
		for (int j = 1; j < neighbors.Length(); j++)
		{
			/* global tag */
			int   tag_j = neighbors[j];
			int  type_j = fType[tag_j];
			double* f_j = fForce(tag_j);
			const double* x_j = coords(tag_j);

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
// DEBUG
//fComm.Log(CommunicatorT::kUrgent, caller);

}

/* set neighborlists */
void ParticlePairT::SetConfiguration(void)
{
	/* inherited */
	ParticleT::SetConfiguration();

	/* reset neighbor lists */
	CommManagerT& comm_manager = ElementSupport().CommManager();
	const ArrayT<int>* part_nodes = comm_manager.PartitionNodes();
	if (fActiveParticles) part_nodes = fActiveParticles;
	
	/* neighbor lists */	
	GenerateNeighborList(part_nodes, fNearestNeighborDistance, fNearestNeighbors, true, true);
	GenerateNeighborList(part_nodes, fNeighborDistance, fNeighbors, false, true);

	/* set the list of reference nearest neighbors */
	if (fRefNearestNeighbors.MajorDim() != fNearestNeighbors.MajorDim() &&
		(fOutputFlags[kSlipVector] || fOutputFlags[kStress] || fOutputFlags[kStrain]))
		SetRefNN(fNearestNeighbors, fRefNearestNeighbors);

	/* output stream */
	ofstreamT& out = ElementSupport().Output();

	/* write the search grid statistics */
	if (fGrid) fGrid->WriteStatistics(out);
	
	out << "\n Neighbor statistics:\n";
	out << " Total number of neighbors . . . . . . . . . . . = " << fNeighbors.Length() << '\n';
	out << " Minimum number of neighbors . . . . . . . . . . = " << fNeighbors.MinMinorDim(0) << '\n';
	out << " Maximum number of neighbors . . . . . . . . . . = " << fNeighbors.MaxMinorDim() << '\n';
	if (fNeighbors.MajorDim() > 0)
	out << " Average number of neighbors . . . . . . . . . . = " << double(fNeighbors.Length())/fNeighbors.MajorDim() << '\n';
	else
	out << " Average number of neighbors . . . . . . . . . . = " << 0 << '\n';

	/* verbose */
	if (ElementSupport().Logging() == GlobalT::kVerbose)
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

/* extract the properties information from the parameter list. See ParticleT::ExtractProperties */
void ParticlePairT::ExtractProperties(const ParameterListT& list, const ArrayT<StringT>& type_names,
	ArrayT<ParticlePropertyT*>& properties, nMatrixT<int>& properties_map)
{
	const char caller[] = "ParticlePairT::ExtractProperties";

	/* check number of interactions */
	int num_props = list.NumLists("pair_particle_interaction");
	int dim = 0;
	for (int i = 0; i < properties_map.Rows(); i++)
		dim += properties_map.Rows() - i;
	if (dim != num_props)
		ExceptionT::GeneralFail(caller, "%d types requires %d \"pair_particle_interaction\"",
			properties_map.Rows(), dim);

	/* read properties */
	fPairProperties.Dimension(num_props);
	for (int i = 0; i < num_props; i++) {

		const ParameterListT& interaction = list.GetList("pair_particle_interaction", i);
		
		/* type names */
		const StringT& label_1 = interaction.GetParameter("label_1");
		const StringT& label_2 = interaction.GetParameter("label_2");
		
		/* resolve index */
		int index_1 = -1, index_2 = -1;
		for (int j = 0; index_1 == -1 && j < type_names.Length(); j++)
			if (type_names[j] == label_1) index_1 = j;
		if (index_1 == -1) ExceptionT::GeneralFail(caller, "could not resolve index of \"%s\"", label_1.Pointer());

		for (int j = 0; index_2 == -1 && j < type_names.Length(); j++)
			if (type_names[j] == label_2) index_2 = j;
		if (index_2 == -1) ExceptionT::GeneralFail(caller, "could not resolve index of \"%s\"", label_2.Pointer());
		
		/* set properies map */
		if (properties_map(index_1, index_2) != -1)
			ExceptionT::GeneralFail(caller, "%s-%s interaction is already defined",
				label_1.Pointer(), label_2.Pointer());
		properties_map(index_1, index_2) = properties_map(index_2, index_1) = i; /* symmetric */
		
		/* read property */
		const ParameterListT& property = interaction.GetListChoice(*this, "pair_property_choice");
		PairPropertyT* pair_prop = PairPropertyT::New(property.Name(), &(ElementSupport()));
		if (!pair_prop) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", property.Name().Pointer());
		pair_prop->TakeParameterList(property);
		fPairProperties[i] = pair_prop;
	}
	
	/* copy */
	properties.Dimension(fPairProperties.Length());
	for (int i = 0; i < properties.Length(); i++)
		properties[i] = fPairProperties[i];
}
