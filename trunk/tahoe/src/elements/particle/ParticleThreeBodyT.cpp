/* $Id: ParticleThreeBodyT.cpp,v 1.2 2004-12-01 02:19:32 cjkimme Exp $ */
#include "ParticleThreeBodyT.h"

#include "ThreeBodyPropertyT.h"
#include "eIntegratorT.h"
#include "InverseMapT.h"
#include "CommManagerT.h"
#include "dSPMatrixT.h"
#include "ModelManagerT.h"
#include "dSymMatrixT.h"
#include "dArray2DT.h"
#include "iGridManagerT.h"
#include "ParameterContainerT.h"

#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>

/* pair property types */
#include "StillingerWeberT.h"

using namespace Tahoe;

/* parameters */
const int kMemoryHeadRoom = 15; /* percent */

/* constructor */
ParticleThreeBodyT::ParticleThreeBodyT(const ElementSupportT& support):
	ParticleT(support),
	fNeighbors(kMemoryHeadRoom),
	fNearestNeighbors(kMemoryHeadRoom),
	fRefNearestNeighbors(kMemoryHeadRoom),
	fEqnos(kMemoryHeadRoom),
	fForce_list_man(0, fForce_list)
{
	SetName("particle_three_body");
	fopen = false;
}

/* collecting element group equation numbers */
void ParticleThreeBodyT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
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
void ParticleThreeBodyT::ConnectsX(AutoArrayT<const iArray2DT*>& connects) const
{
	/* NOTE: do not add anything to the geometry connectivity list */
#pragma unused(connects)
}

/* collecting element field connectivities */
void ParticleThreeBodyT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
#pragma unused(connects_1)
	connects_2.AppendUnique(&fNeighbors);
}

void ParticleThreeBodyT::WriteOutput(void)
{
	const char caller[] = "ParticleThreeBodyT::WriteOutput";

	/* inherited */
	ParticleT::WriteOutput();

	/* muli-processor information */
	CommManagerT& comm_manager = ElementSupport().CommManager();
	const ArrayT<int>* proc_map = comm_manager.ProcessorMap();
	int rank = ElementSupport().Rank();

	/* dimensions */
	int ndof = NumDOF();
	int num_output = ndof + 2; /* displacement + PE + KE */

#ifndef NO_PARTICLE_STRESS_OUTPUT
	num_output++; /* includes centrosymmetry */
	num_output+=ndof; /* some more for slip vector */
#endif

	/* number of nodes */
	const ArrayT<int>* parition_nodes = comm_manager.PartitionNodes();
	int non = (parition_nodes) ? 
		parition_nodes->Length() : 
		ElementSupport().NumNodes();

	/* map from partition node index */
	const InverseMapT* inverse_map = comm_manager.PartitionNodes_inv();

#ifndef NO_PARTICLE_STRESS_OUTPUT
	dSymMatrixT vs_i(ndof), temp(ndof);
	int num_stresses = vs_i.NumValues(ndof);
	//dArray2DT vsvalues(non, num_stresses);
	num_output += num_stresses;
	num_output += num_stresses; //another for the strain
#endif

	/* output arrays length number of active nodes */
	dArray2DT n_values(non, num_output), e_values;
	n_values = 0.0;

	/* global coordinates */
	const dArray2DT& coords = ElementSupport().CurrentCoordinates();

	/* pair properties function pointers */
	int current_property = -1;
	PairPropertyT::EnergyFunction energy_function_2body = NULL;	
	PairPropertyT::ForceFunction force_function_2body = NULL;

	/* three body properties function pointers */
	ThreeBodyPropertyT::EnergyFunction energy_function_3body = NULL;
	ThreeBodyPropertyT::ForceFunction force_function_3body = NULL;

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
		mass[i] = fThreeBodyProperties[fPropertiesMap(i,i)]->Mass();

	/* collect displacements */
	dArrayT vec, values_i;
	for (int i = 0; i < non; i++) {
		int   tag_i = (parition_nodes) ? (*parition_nodes)[i] : i;
		int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;
		int  type_i = fType[tag_i];

		/* values for particle i */
		n_values.RowAlias(local_i, values_i);

		/* copy in */
		vec.Set(ndof, values_i.Pointer());
		displacement.RowCopy(tag_i, vec);

#ifndef NO_PARTICLE_STRESS_OUTPUT
		/* kinetic contribution to the virial */
		if (velocities) {
			velocities->RowAlias(tag_i, vec);
			temp.Outer(vec);
		 	for (int cc = 0; cc < num_stresses; cc++) {
				int ndex = ndof+2+cc;
		   		values_i[ndex] = (fabs(V0) > kSmall) ? -mass[type_i]*temp[cc]/V0 : 0.0;
		 	}
		}
#endif
	}
	
	/* run through neighbor list */
	iArrayT neighbors;
	dArrayT x_i, x_j, r_ij(ndof);

#ifndef NO_PARTICLE_STRESS_OUTPUT
	dArrayT SlipVector(ndof);
	dMatrixT Strain(ndof);
#endif
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

#ifndef NO_PARTICLE_STRESS_OUTPUT
		vs_i = 0.0;
#endif
		
		/* kinetic energy */
		if (velocities) {
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
			if (property != current_property) {
				energy_function_2body = fThreeBodyProperties[property]->getEnergyFunction();
				force_function_2body = fThreeBodyProperties[property]->getForceFunction();
				current_property = property;
			}
		
			/* global coordinates */
			coords.RowAlias(tag_j, x_j);

			/* connecting vector */
			r_ij.DiffOf(x_j, x_i);
			double r = r_ij.Magnitude();

			/* split interaction energy */
			double uby2 = 0.5*energy_function_2body(r, NULL, NULL); 
			values_i[ndof] += uby2;
			
	      	/* interaction force */
			double F = force_function_2body(r, NULL, NULL);
			double Fbyr = F/r;

#ifndef NO_PARTICLE_STRESS_OUTPUT
			temp.Outer(r_ij);
			vs_i.AddScaled(0.5*Fbyr, temp);
#endif

			/* second node may not be on processor */
			if (!proc_map || (*proc_map)[tag_j] == rank) 
			{
				int local_j = (inverse_map) ? inverse_map->Map(tag_j) : tag_j;	
				if (local_j < 0 || local_j >= n_values.MajorDim())
					cout << caller << ": out of range: " << local_j << '\n';
				else {

					/* potential energy */
					n_values(local_j, ndof) += uby2;

#ifndef NO_PARTICLE_STRESS_OUTPUT
			 		/* accumulate into stress into array */
		 			for (int cc = 0; cc < num_stresses; cc++) {
						int ndex = ndof+2+cc;
		   				n_values(local_j, ndex) += (fabs(V0) > kSmall) ? 0.5*Fbyr*temp[cc]/V0 : 0.0;
		 			}
#endif
				}
			}
		}


#ifndef NO_PARTICLE_STRESS_OUTPUT
		/* copy stress into array */
		for (int cc = 0; cc < num_stresses; cc++) {
		  int ndex = ndof+2+cc;
		  values_i[ndex] += (fabs(V0) > kSmall) ? vs_i[cc]/V0 : 0.0;
		}
#endif
	}
#ifndef NO_PARTICLE_STRESS_OUTPUT
    int num_s_vals = num_stresses+1+ndof+1;
    dArray2DT s_values(non,num_s_vals);
    s_values = 0.0;

    /* flag for specifying Lagrangian (0) or Eulerian (1) strain */ 
    const int kEulerLagr = 0;

	/* calculate slip vector and strain */
	Calc_Slip_and_Strain(s_values, fRefNearestNeighbors, kEulerLagr);

    /* calculate centrosymmetry parameter */
	Calc_CSP(s_values, fNearestNeighbors);

	/* combine strain, slip vector and centrosymmetry parameter into n_values list */
	for (int i = 0; i < fNeighbors.MajorDim(); i++)
	{
		/* row of neighbor list */
		fNeighbors.RowAlias(i, neighbors);

		/* tags */
		int   tag_i = neighbors[0]; /* self is 1st spot */
		int  type_i = fType[tag_i];
		int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;
									            
		int valuep = 0;
		for (int is = 0; is < num_stresses; is++)
			n_values(local_i, ndof+2+num_stresses+valuep++) = s_values(local_i,is);

		/* recover J, the determinant of the deformation gradient, for atom i
		 * and divide stress values by it */
		double J = s_values(local_i,num_stresses);
		for (int is = 0; is < num_stresses; is++) 
		  n_values(local_i,ndof+2+is) /= J;

		for (int n = 0; n < ndof; n++)
			n_values(local_i, ndof+2+num_stresses+num_stresses+n) = s_values(local_i, num_stresses+1+n);

		n_values(local_i, num_output-1) = s_values(local_i, num_s_vals-1);
	}
#endif

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

/* compute the part of the stiffness matrix */
void ParticleThreeBodyT::FormStiffness(const InverseMapT& col_to_col_eq_row_map,
	const iArray2DT& col_eq, dSPMatrixT& stiffness)
{
	const char caller[] = "ParticleThreeBodyT::FormStiffness";

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
	PairPropertyT::ForceFunction force_function_2body = NULL;
	PairPropertyT::StiffnessFunction stiffness_function_2body = NULL;

	/* three body properties function pointers */
	ThreeBodyPropertyT::ForceFunction force_function_3body = NULL;
	ThreeBodyPropertyT::StiffnessFunction stiffness_function_3body = NULL;

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
					force_function_2body = fThreeBodyProperties[property]->getForceFunction();
					stiffness_function_2body = fThreeBodyProperties[property]->getStiffnessFunction();
					current_property = property;
				}

				/* global coordinates */
				coords.RowAlias(tag_j, x_j);

				/* connecting vector */
				r_ij.DiffOf(x_j, x_i);
				double r = r_ij.Magnitude();
				r_ji.SetToScaled(-1.0, r_ij);

				/* interaction functions */
				double F = force_function_2body(r, NULL, NULL);
				double K = stiffness_function_2body(r, NULL, NULL);
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
void ParticleThreeBodyT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParticleT::DefineSubs(sub_list);

	/* interactions */
	sub_list.AddSub("three_body_particle_interaction", ParameterListT::OnePlus);
}

/* return the description of the given inline subordinate parameter list */
void ParticleThreeBodyT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "three_body_property_choice") {
		order = ParameterListT::Choice;
	 
		/* Stillinger-Weber */
		sub_lists.AddSub("Stillinger_Weber");
	} 
	else /* inherited */
		ParticleT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ParticleThreeBodyT::NewSub(const StringT& name) const
{
	/* try to construct potential */
	ThreeBodyPropertyT* three_body_property = ThreeBodyPropertyT::New(name, &(ElementSupport()));
	if (three_body_property)
	  	return three_body_property;
	else if (name == "three_body_particle_interaction") {
		ParameterContainerT* interactions = new ParameterContainerT(name);
		interactions->SetSubSource(this);
		
		/* particle type labels */
		interactions->AddParameter(ParameterT::Word, "label_1");
		interactions->AddParameter(ParameterT::Word, "label_2");
	    interactions->AddParameter(ParameterT::Word, "label_3");

		/* properties choice list */
		interactions->AddSub("three_body_property_choice", ParameterListT::Once, true);
		
		return interactions;
	}
	else /* inherited */
		return ParticleT::NewSub(name);
}

/* accept parameter list */
void ParticleThreeBodyT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParticleT::TakeParameterList(list);

	/* set the list of reference nearest neighbors */
	SetRefNN(fNearestNeighbors, fRefNearestNeighbors);

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

/* generate labels for output data */
void ParticleThreeBodyT::GenerateOutputLabels(ArrayT<StringT>& labels) const
{
  int ndof=NumDOF();
	if (ndof > 3) ExceptionT::GeneralFail("ParticleThreeBodyT::GenerateOutputLabels");

	/* displacement labels */
	const char* disp[3] = {"D_X", "D_Y", "D_Z"};
	const char* SV[3] = {"SV_X", "SV_Y", "SV_Z"};
	int num_labels =
		ndof // displacements
		+ 2;     // PE and KE

#ifndef NO_PARTICLE_STRESS_OUTPUT
	int num_stress=0;
	const char* stress[6];
	const char* strain[6];
	if (ndof==3){
	  num_stress=6;
	  stress[0]="s11";
	  stress[1]="s22";
	  stress[2]="s33";
	  stress[3]="s23";
	  stress[4]="s13";
	  stress[5]="s12";
	  }
	  else if (ndof==2) {
	   num_stress=3;
	  stress[0]="s11";
	  stress[1]="s22";
	  stress[2]="s12";
	  }
	  else if (ndof==1) {
	   num_stress=1;
	  stress[0] = "s11";
	  }
	if (ndof==3){
	  
	  strain[0]="e11";
	  strain[1]="e12";
	  strain[2]="e13";
	  strain[3]="e22";
	  strain[4]="e23";
	  strain[5]="e33";
	  }
	  else if (ndof==2) {
	   
	  strain[0]="e11";
	  strain[1]="e12";
	  strain[2]="e22";
	  }
	  else if (ndof==1) {

	  strain[0] = "e11";
	  }
	num_labels+=num_stress;
	num_labels++; //another label for the centrosymmetry
	num_labels+=num_stress; //another for the strain
	num_labels+=ndof; /*and another for the slip vector*/
#endif /* NO_PARTICLE_STRESS_OUTPUT */
	
	labels.Dimension(num_labels);
	int dex = 0;
	for (dex = 0; dex < NumDOF(); dex++)
		labels[dex] = disp[dex];
	labels[dex++] = "PE";
	labels[dex++] = "KE";

#ifndef NO_PARTICLE_STRESS_OUTPUT
	for (int ns =0 ; ns<num_stress; ns++)
	  labels[dex++]=stress[ns];
	for (int ns =0 ; ns<num_stress; ns++)
	  labels[dex++]=strain[ns];
	for (int i=0; i<ndof; i++)
	  labels[dex++]=SV[i];
	labels[dex++]= "CS";
#endif /* NO_PARTICLE_STRESS_OUTPUT */
}

/* form group contribution to the stiffness matrix */
void ParticleThreeBodyT::LHSDriver(GlobalT::SystemTypeT sys_type)
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
			mass[i] = fThreeBodyProperties[fPropertiesMap(i,i)]->Mass();
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
		ThreeBodyPropertyT::StiffnessFunction stiffness_function_3body = NULL;

		/* run through neighbor list */
		fForce = 0.0;
		iArrayT neighbors;
		dArrayT x_i, x_j, r_ij(ndof);
		const double *x_ip, *x_jp, *x_kp;
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
			x_ip = coords(tag_i);
			for (int j = 1; j < neighbors.Length(); j++)
			{
				/* global tag */
				int  tag_j = neighbors[j];
				if (tag_j < tag_i) {
					int type_j = fType[tag_j];
					double* k_j = fForce(tag_j);
				
					/* set pair property (if not already set) */
					int property = fPropertiesMap(type_i, type_j);
					if (property != current_property)
					{
						force_function = fThreeBodyProperties[property]->getForceFunction();
						stiffness_function = fThreeBodyProperties[property]->getStiffnessFunction();
						current_property = property;
					}
			
					/* global coordinates */
					coords.RowAlias(tag_j, x_j);
					x_jp = coords(tag_j);
			
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
					
					// additional loop over neighbors for 3-body terms
					for (int k = 1; k < neighbors.Length(); j++)
					{
						/* global tag */
						int   tag_k = neighbors[k];
						if (tag_k < tag_j) {
							int  type_k = fType[tag_k];
							double* f_k = fForce(tag_k);
							x_kp = coords(tag_k);

							/* set pair property (if not already set) */
							int property = fPropertiesMap(type_i, type_j);
							if (property != current_property)
							{
								stiffness_function_3body = fThreeBodyProperties[property]->getThreeBodyStiffnessFunction();
								current_property = property;
							}
						
							/* stiffness matrix */
							dMatrixT K_ijk(3*ndof);
							if (stiffness_function_3body(x_ip, x_jp, x_kp, K_ijk)) {
								// check signs and make sure I don't need factors of 1/3
								
							} 
						}
					}
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
		ThreeBodyPropertyT::StiffnessFunction stiffness_function_3body = NULL;

		/* work space */
		dArrayT r_ij(NumDOF(), fRHS.Pointer());
		dArrayT r_ji(NumDOF(), fRHS.Pointer() + NumDOF());

		/* run through neighbor list */
		const iArray2DT& field_eqnos = Field().Equations();
		iArray2DT pair_eqnos(2, ndof), triple_eqnos(3,ndof); 
		iArrayT pair(2), triple(3);
		iArrayT neighbors;
		dArrayT x_i, x_j;
		const double *x_ip, *x_jp, *x_kp;
		for (int i = 0; i < fNeighbors.MajorDim(); i++)
		{
			/* row of neighbor list */
			fNeighbors.RowAlias(i, neighbors);

			/* type */
			int  tag_i = neighbors[0]; /* self is 1st spot */
			int type_i = fType[tag_i];
			pair[0] = triple[0] = tag_i;
		
			/* run though neighbors for one atom - first neighbor is self */
			coords.RowAlias(tag_i, x_i);
			x_ip = coords(tag_i);
			for (int j = 1; j < neighbors.Length(); j++)
			{
				/* global tag */
				int  tag_j = neighbors[j];
				if (tag_j < tag_i) {
					int type_j = fType[tag_j];
					pair[1] = triple[1] = tag_j;
			
					/* set pair property (if not already set) */
					int property = fPropertiesMap(type_i, type_j);
					if (property != current_property)
					{
						force_function = fThreeBodyProperties[property]->getForceFunction();
						stiffness_function = fThreeBodyProperties[property]->getStiffnessFunction();
						current_property = property;
					}
			
					/* global coordinates */
					coords.RowAlias(tag_j, x_j);
					x_jp = coords(tag_j);
			
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
				
					// additional loop over neighbors for 3-body terms
					for (int k = 1; k < neighbors.Length(); j++)
					{
						/* global tag */
						int   tag_k = neighbors[k];
						if (tag_k < tag_j) {
							int  type_k = fType[tag_k];
							triple[2] = tag_k;
							double* f_k = fForce(tag_k);
							x_kp = coords(tag_k);

							/* set pair property (if not already set) */
							int property = fPropertiesMap(type_i, type_j);
							if (property != current_property)
							{
								stiffness_function_3body = fThreeBodyProperties[property]->getThreeBodyStiffnessFunction();
								current_property = property;
							}
						
							/* stiffness matrix */
							dMatrixT K_ijk(3*ndof);
							if (stiffness_function_3body(x_ip, x_jp, x_kp, K_ijk)) {
								/* assemble */
								triple_eqnos.RowCollect(triple, field_eqnos);
								support.AssembleLHS(group, K_ijk, pair_eqnos);
							} 
						}
					}
				}
			}
		}
	}
}

/* form group contribution to the residual */
void ParticleThreeBodyT::RHSDriver(void)
{
	int nsd = NumSD();
	if (nsd == 3)
		RHSDriver3D();
	else
		ExceptionT::GeneralFail("ParticleThreeBodyT::RHSDriver", "unsupported dimension %d", nsd);
		
	ApplyDamping(fNeighbors);
	
	/* assemble */
	ElementSupport().AssembleRHS(Group(), fForce, Field().Equations());
}

void ParticleThreeBodyT::RHSDriver3D(void)
{
	/* function name */
	const char caller[] = "ParticleThreeBodyT::RHSDriver3D";

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

	/* function pointers, etc. */
	int current_property = -1;
	PairPropertyT::ForceFunction force_function = NULL;
	ThreeBodyPropertyT::ForceFunction force_function_3body = NULL;
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
			if (tag_j < tag_i) {
				int  type_j = fType[tag_j];
				double* f_j = fForce(tag_j);
				const double* x_j = coords(tag_j);

				/* set pair property (if not already set) */
				int property = fPropertiesMap(type_i, type_j);
				if (property != current_property)
				{
					if (!fThreeBodyProperties[property]->getParadynTable(&Paradyn_table, dr, row_size, num_rows))
						force_function = fThreeBodyProperties[property]->getForceFunction();
					force_function_3body = fThreeBodyProperties[property]->getThreeBodyForceFunction();
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
				
				// additional loop over neighbors for 3-body terms
				for (int k = 1; k < neighbors.Length(); k++)
				{
					/* global tag */
					int   tag_k = neighbors[k];
					if (tag_k < tag_j) {
						int  type_k = fType[tag_k];
						double* f_k = fForce(tag_k);
						const double* x_k = coords(tag_k);
					
						/* connecting vector */
						double r_ik[3];
						r_ik[0] = x_k[0] - x_i[0];
						r_ik[1] = x_k[1] - x_i[1];
						r_ik[2] = x_k[2] - x_i[2];
						double rik = sqrt(r_ik[0]*r_ik[0] + r_ik[1]*r_ik[1] + r_ik[2]*r_ik[2]);
						
						/* interaction force */
						double f_ij[3], f_ik[3];
						if (force_function_3body(x_i, x_j, x_k, f_ij, f_ik)) {
							f_ij[0] *= formKd; 
							f_i[0] -= f_ij[0] + f_ik[0];
							f_j[0] += f_ij[0];
							f_k[0] += f_ik[0]; 

							f_ij[1] *= formKd; 
							f_i[1] -= f_ij[1] + f_ik[1];
							f_j[1] += f_ij[1];
							f_k[1] += f_ik[1]; 
						
							f_ij[0] *= formKd; 
							f_i[2] -= f_ij[2] + f_ik[2];
							f_j[2] += f_ij[2];
							f_k[2] += f_ik[2];
						} 
					}
				}
			}
		}
	}
}

/* set neighborlists */
void ParticleThreeBodyT::SetConfiguration(void)
{
	/* inherited */
	ParticleT::SetConfiguration();

	/* reset neighbor lists */
	CommManagerT& comm_manager = ElementSupport().CommManager();
	const ArrayT<int>* part_nodes = comm_manager.PartitionNodes();
	if (fActiveParticles) 
		part_nodes = fActiveParticles;
		
	GenerateNeighborList(part_nodes, fNearestNeighborDistance, fNearestNeighbors, true, true);
	// need full neighbor list for now
	GenerateNeighborList(part_nodes, fNeighborDistance, fNeighbors, true, true);

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

/* extract the properties information from the parameter list. See ParticleT::ExtractProperties */
void ParticleThreeBodyT::ExtractProperties(const ParameterListT& list, const ArrayT<StringT>& type_names,
	ArrayT<ParticlePropertyT*>& properties, nMatrixT<int>& properties_map)
{
	const char caller[] = "ParticleThreeBodyT::ExtractProperties";

	/* check number of interactions */
	int num_props = list.NumLists("three_body_particle_interaction");
	int dim = 0;
	for (int i = 0; i < properties_map.Rows(); i++)
		dim += properties_map.Rows() - i;
	if (dim != num_props)
		ExceptionT::GeneralFail(caller, "%d types requires %d \"three_body_particle_interaction\"",
			properties_map.Rows(), dim);

	/* read properties */
	fThreeBodyProperties.Dimension(num_props);
	for (int i = 0; i < num_props; i++) {

		const ParameterListT& interaction = list.GetList("three_body_particle_interaction", i);
		
		/* type names */
		const StringT& label_1 = interaction.GetParameter("label_1");
		const StringT& label_2 = interaction.GetParameter("label_2");
		const StringT& label_3 = interaction.GetParameter("label_3");
		
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
		const ParameterListT& property = interaction.GetListChoice(*this, "three_body_property_choice");
		ThreeBodyPropertyT* three_body_prop = ThreeBodyPropertyT::New(property.Name(), &(ElementSupport()));
		if (!three_body_prop) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", property.Name().Pointer());
		three_body_prop->TakeParameterList(property);
		fThreeBodyProperties[i] = three_body_prop;
	}
	
	/* copy */
	properties.Dimension(fThreeBodyProperties.Length());
	for (int i = 0; i < properties.Length(); i++)
		properties[i] = fThreeBodyProperties[i];
}
