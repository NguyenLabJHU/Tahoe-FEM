/* $Id: ParticleT.cpp,v 1.43 2004-07-22 08:23:04 paklein Exp $ */
#include "ParticleT.h"

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "eIntegratorT.h"
#include "OutputSetT.h"
#include "dArray2DT.h"
#include "ElementSupportT.h"
#include "ModelManagerT.h"
#include "iGridManagerT.h"
#include "iNodeT.h"
#include "ParticlePropertyT.h"
#include "RaggedArray2DT.h"
#include "CommManagerT.h"
#include "CommunicatorT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"

/* Thermostatting stuff */
#include "RandomNumberT.h"
#include "ScheduleT.h"
#include "ThermostatBaseT.h"
#include "GaussIsokineticT.h"
#include "LangevinT.h"
#include "NoseHooverT.h"
#include "RampedDampingT.h"

using namespace Tahoe;

/* class parameters */
/* parameters */
const int kAvgNodesPerCell = 20;
const int kMaxNumCells     =- 1; /* -1: no max */

/* constructors */
ParticleT::ParticleT(const ElementSupportT& support):
	ElementBaseT(support),
	fNeighborDistance(-1),
	fReNeighborDisp(-1),
	fReNeighborIncr(-1),
	fGrid(NULL),
	fReNeighborCounter(0),
	fDmax(0),
	fForce_man(fForce),
	fActiveParticles(NULL),
	fTypeMessageID(CommManagerT::kNULLMessageID),
	fLatticeParameter(-1.0)
{
	SetName("particle");

	/* set matrix format */
	fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
}

/* destructor */
ParticleT::~ParticleT(void)
{
	/* free search grid */
	delete fGrid;

	/* free properties list */
	for (int i = 0; i < fParticleProperties.Length(); i++)
		delete fParticleProperties[i];
		
	delete fActiveParticles;
	
	/* thermostats */
	for (int i = 0; i < fThermostats.Length(); i++)
		delete fThermostats[i];
}

/* form of tangent matrix */
GlobalT::SystemTypeT ParticleT::TangentType(void) const {
	return GlobalT::kSymmetric;
}

/* NOT implemented. Returns an zero force vector */
void ParticleT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
}

/* writing output */
void ParticleT::RegisterOutput(void)
{
	/* "point connectivities" needed for output */
	CommManagerT& comm_manager = ElementSupport().CommManager();
	const ArrayT<int>* parition_nodes = comm_manager.PartitionNodes();
	if (parition_nodes)
	{
		int num_nodes = parition_nodes->Length();
		fPointConnectivities.Alias(num_nodes, 1, parition_nodes->Pointer());
	}
	else /* ALL nodes */
	{
		fPointConnectivities.Dimension(ElementSupport().NumNodes(), 1);
		iArrayT tmp;
		tmp.Alias(fPointConnectivities);
		tmp.SetValueToPosition();				
	}

	/* block ID's */
	ArrayT<StringT> block_ID(fBlockData.Length());
	for (int i = 0; i < block_ID.Length(); i++)
		block_ID[i] = fBlockData[i].ID();

	/* get output labels (per node) */
	ArrayT<StringT> n_labels, e_labels;
	GenerateOutputLabels(n_labels);

	/* set output specifier */
	StringT set_ID;
	set_ID.Append(ElementSupport().ElementGroupNumber(this) + 1);
	OutputSetT output_set(GeometryT::kPoint, fPointConnectivities, n_labels, ChangingGeometry());
		
	/* register and get output ID */
	fOutputID = ElementSupport().RegisterOutput(output_set);
}

void ParticleT::WriteOutput(void)
{
	/* max distance traveled since last reneighboring */
	ofstreamT& out = ElementSupport().Output();
	out << "\n Maximum displacement since last re-neighboring. = " << fDmax << '\n';

	/* info about periodic boundaries */
	CommManagerT& comm_manager = ElementSupport().CommManager();
	const dArray2DT& periodic_bounds = comm_manager.PeriodicBoundaries();
	out << " Periodic bounds:\n";
	for (int i = 0; i < periodic_bounds.MajorDim(); i++)
		out << i+1 << ": {" << periodic_bounds(i,0) << ", " << periodic_bounds(i,1) << "}\n";
	
	/* reset connectivities */
	if (ChangingGeometry())
	{
		const ArrayT<int>* parition_nodes = comm_manager.PartitionNodes();
		if (parition_nodes)
		{
			int num_nodes = parition_nodes->Length();
			fPointConnectivities.Alias(num_nodes, 1, parition_nodes->Pointer());	
		}
		else
			ExceptionT::GeneralFail("ParticleT::WriteOutput", "expecting a partition nodes list");
	}
}

/* compute specified output parameter and send for smoothing */
void ParticleT::SendOutput(int kincode)
{
#pragma unused(kincode)
	//TEMP: for now, do nothing
}

/* trigger reconfiguration */
GlobalT::RelaxCodeT ParticleT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	/* multiprocessor support */
	CommManagerT& comm_manager = ElementSupport().CommManager();

	/* compute max distance traveled since last neighboring 
	 * (across all processes) */
	fDmax = comm_manager.Communicator().Max(MaxDisplacement());

	/* check damping regions */
	//fDampingCounters++;

	/* reset periodic bounds given stretching */
	bool has_moving = false;
	for (int i = 0; i < NumSD(); i++)
	{
		const ScheduleT* stretch = fStretchSchedule[i];
		if (stretch)
		{
			/* time during next solution step */
			double next_time = ElementSupport().Time() + ElementSupport().TimeStep();
		
			has_moving = true;
			double scale = stretch->Value(next_time);
			double x_min = scale*fPeriodicBounds(i,0);
			double x_max = scale*fPeriodicBounds(i,1);
	
			/* redefine bounds */
			comm_manager.SetPeriodicBoundaries(i, x_min, x_max);
		}
	}

	/* generate contact element data */
	fReNeighborCounter++;
	if (has_moving ||
	    (fReNeighborDisp > 0.0 && fDmax > fReNeighborDisp) || 
		(fReNeighborIncr > 0 && fReNeighborCounter >= fReNeighborIncr))
	{
		/* output stream */
		ofstreamT& out = ElementSupport().Output();
		if (fReNeighborDisp > 0.0 && fDmax > fReNeighborDisp)
			out << "\n ParticleT::RelaxSystem: max displacement since re-neighboring "
			    << fDmax << " > " << fReNeighborDisp << '\n';
		if (fReNeighborIncr > 0 && fReNeighborCounter >= fReNeighborIncr)
			out << "\n ParticleT::RelaxSystem: number of steps since re-neighboring "
			    << fReNeighborCounter << " >= " << fReNeighborIncr << '\n';
	
		/* (re-)set the neighborlists */
		SetConfiguration();

		/* reset counter */
		fReNeighborCounter = 0;
	
		return GlobalT::MaxPrecedence(relax, GlobalT::kReEQ);
	}
	else
		return relax;
}

/* write restart data to the output stream */
void ParticleT::WriteRestart(ostream& out) const
{
	/* write counter */
	out << fReNeighborCounter << '\n';
	
	for (int i = 0; i < fThermostats.Length(); i++)
		fThermostats[i]->WriteRestart(out);
}

/* read restart data to the output stream */
void ParticleT::ReadRestart(istream& in)
{
	/* read counter */
	in >> fReNeighborCounter;

	for (int i = 0; i < fThermostats.Length(); i++)
		fThermostats[i]->ReadRestart(in);
}

/* define the particles to skip */
void ParticleT::SetSkipParticles(const iArrayT& skip)
{
	if (skip.Length() == 0) {
		delete fActiveParticles;
		fActiveParticles = NULL;
	}
	else
	{
		int nnd = ElementSupport().NumNodes();
		iArrayT nodes_used(nnd);

		/* mark partition nodes as used */
		CommManagerT& comm_manager = ElementSupport().CommManager();
		const ArrayT<int>* part_nodes = comm_manager.PartitionNodes();
		if (part_nodes)
		{
			nodes_used = 0;
			int npn = part_nodes->Length();
			const int* p = part_nodes->Pointer();
			for (int i = 0; i < npn; i++)
				nodes_used[*p++] = 1;
		}
		else /* all are partition nodes */
			nodes_used = 1;

		/* mark nodes to skip */
		int nsn = skip.Length();
		const int* ps = skip.Pointer();
		for (int i = 0; i < nsn; i++)
			nodes_used[*ps++] = 0;
			
		
		/* collect active particles */	
		int nap = nodes_used.Count(1);
		if (!fActiveParticles)
			fActiveParticles = new AutoArrayT<int>;
		fActiveParticles->Dimension(nap);
		int dex = 0;
		for (int i = 0; i < nnd; i++)
			if (nodes_used[i] == 1)
				(*fActiveParticles)[dex++] = i;
	}
}

/* set neighborlists */
void ParticleT::SetConfiguration(void)
{
	/* set periodic boundary conditions */
	CommManagerT& comm_manager = ElementSupport().CommManager();
	comm_manager.EnforcePeriodicBoundaries(fPeriodicSkin);
	
	/* reset the types array */
	int nnd = ElementSupport().NumNodes();
	fType.Resize(nnd);
	
	/* exchange type information */
	if (fTypeMessageID == CommManagerT::kNULLMessageID)
		fTypeMessageID = ElementSupport().CommManager().Init_AllGather(MessageT::Integer, 1);
	iArray2DT type_wrapper(fType.Length(), 1, fType.Pointer());
	comm_manager.AllGather(fTypeMessageID, type_wrapper);
	
	/* resize working arrays */
	fForce_man.SetMajorDimension(nnd, false);

	/* collect current coordinates */
	const ArrayT<int>* part_nodes = comm_manager.PartitionNodes();
	const dArray2DT& curr_coords = ElementSupport().CurrentCoordinates();
	if (part_nodes)
	{
		/* collect */
		fReNeighborCoords.Dimension(part_nodes->Length(), curr_coords.MinorDim());
		fReNeighborCoords.RowCollect(*part_nodes, curr_coords);
	}
	else /* use ALL nodes */
		fReNeighborCoords = curr_coords;
}

/* contribution to the nodal residual forces */
const dArray2DT& ParticleT::InternalForce(int group)
{
	/* check */
	if (group != Group())
		ExceptionT::GeneralFail("ParticleT::InternalForce", 
			"expecting solver group %d not %d", Group(), group);
	return fForce;
}

/* add the element group's contribution to the lumped (scalar) mass of the given nodes */
void ParticleT::LumpedMass(const iArrayT& nodes, dArrayT& mass) const
{
	/* inherited */
	ElementBaseT::LumpedMass(nodes, mass);

	/* collect particle masses */
	for (int i = 0; i < nodes.Length(); i++)
	{
		int property_type = fType[nodes[i]];
		mass[i] += fParticleProperties[property_type]->Mass();
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

void ParticleT::ApplyDamping(const RaggedArray2DT<int>& fNeighbors)
{		
	if (fThermostats.Length() > 0)
	{
		const dArray2DT* velocities = NULL; 
     	if (Field().Order() > 0) // got velocities!
     	{
     		velocities = &(Field()[1]);
     		
     		for (int i = 0; i < fThermostats.Length(); i++)
				fThermostats[i]->ApplyDamping(fNeighbors,velocities,fForce,
										fType,fParticleProperties);
		}
	}
		
}

/* return true if connectivities are changing */
bool ParticleT::ChangingGeometry(void) const {
	return ElementSupport().CommManager().PartitionNodesChanging();
}

/* generate neighborlist */
void ParticleT::GenerateNeighborList(const ArrayT<int>* particle_tags, 
	double distance, RaggedArray2DT<int>& neighbors, 
	bool double_list, bool full_list)
{
	/* global coordinates */
	const dArray2DT& coords = ElementSupport().CurrentCoordinates();

	/* node to processor map */
	CommManagerT& comm_manager = ElementSupport().CommManager();
	const ArrayT<int>* n2p_map = comm_manager.ProcessorMap();

	/* construct grid (using all local nodes) */
	if (!fGrid) fGrid = new iGridManagerT(kAvgNodesPerCell, kMaxNumCells, coords, NULL);

	/* reset contents */
	fGrid->Reset();
	
	/* set up temp space */
	int init_num_neighbors = 6;
	int num_tags = (particle_tags) ? particle_tags->Length() : coords.MajorDim();
	int num_chunks = num_tags/250;
	num_chunks = (num_chunks < 1) ? 1 : num_chunks; 
	AutoFill2DT<int> auto_neighbors(num_tags, num_chunks, 20, init_num_neighbors);

	/* mark nodes owned by this processor, but skipped and still needed to
	 * ensure a full neighbor list of nodes in particle_tags */
	if (double_list) full_list = false;
	ArrayT<int> skipped;
	const ArrayT<int>* partition_nodes = comm_manager.PartitionNodes();
	int npn = (partition_nodes) ? partition_nodes->Length() : ElementSupport().NumNodes();
	if (full_list && particle_tags && npn > num_tags) {
		skipped.Dimension(npn);
		skipped = 1;
		for (int i = 0; i < particle_tags->Length(); i++)
			skipped[(*particle_tags)[i]] = 0; /* not skipped */
	}
	
	/* loop over tags */
	int nsd = coords.MinorDim();
	double distance2 = distance*distance;
	for (int i = 0; i < num_tags; i++)
	{
		/* this tag */
		int tag_i = (particle_tags) ? (*particle_tags)[i] : i;
		int  pr_i = (n2p_map) ? (*n2p_map)[tag_i] : 0;

		/* add self */
		auto_neighbors.Append(i, tag_i);
		
		/* gets points from grid */
		const double* coords_i = coords(tag_i);
		const AutoArrayT<iNodeT>& hits = fGrid->HitsInRegion(coords_i, distance);
		
		/* filter neighbors */
		for (int j = 0; j < hits.Length(); j++)
		{
			int tag_j = hits[j].Tag();
			int  pr_j = (n2p_map) ? (*n2p_map)[tag_j] : 0;
			
			if (double_list || /* double-linked neighbors */
				tag_j > tag_i || /* upper half */
				(full_list && tag_j < tag_i && skipped.Length() > 0 && skipped[tag_j]) || /* lower half + bonds to skipped nodes */
				(full_list && pr_i != pr_j && tag_j != tag_i)) /* lower half + off-processor bonds */
			{
				/* hit info */
				const double* coords_hit = hits[j].Coords();
			
				/* distance^2 */
				double d2 = 0.0;
				for (int k = 0; k < nsd; k++)
				{
					double dx = coords_i[k] - coords_hit[k];
					d2 += dx*dx;
				}
		
				/* it's a keeper */
				if (d2 <= distance2) 
					auto_neighbors.Append(i, tag_j);
			}
		}
	}
	
	/* copy/compress into return array */
	neighbors.Copy(auto_neighbors);
}

/* assemble particle mass matrix into LHS of global equation system */
void ParticleT::AssembleParticleMass(const dArrayT& mass)
{
	/* partition nodes */
	CommManagerT& comm_manager = ElementSupport().CommManager();
	const ArrayT<int>* part_nodes = comm_manager.PartitionNodes();
	
	fForce = 0.0;
	if (part_nodes) /* not all nodes */
	{
		for (int i = 0; i < part_nodes->Length(); i++)
		{
			int tag = (*part_nodes)[i];
		
			/* assemble into global array */
			fForce.SetRow(tag, mass[fType[tag]]);
		}
	}
	else /* all nodes */
	{
		int nnd = ElementSupport().NumNodes();
		for (int i = 0; i < nnd; i++)
			/* assemble into global array */
			fForce.SetRow(i, mass[fType[i]]);
	}
	
	/* assemble all */
	ElementSupport().AssembleLHS(Group(), fForce, Field().Equations());
}

/* return the maximum distance */
double ParticleT::MaxDisplacement(void) const
{
	const char caller[] = "ParticleT::MaxDisplacement";
	CommManagerT& comm_manager = ElementSupport().CommManager();
	const ArrayT<int>* part_nodes = comm_manager.PartitionNodes();
	const dArray2DT& curr_coords = ElementSupport().CurrentCoordinates();
	double dmax2 = 0.0;
	int nsd = curr_coords.MinorDim();
	const double *p_old = fReNeighborCoords.Pointer();
	if (part_nodes)
	{
		int nnd = part_nodes->Length();
		if (nnd != fReNeighborCoords.MajorDim()) ExceptionT::SizeMismatch(caller);
		if (nsd == 3)
		{
			for (int i = 0; i < nnd; i++)
			{
				const double* p_new = curr_coords((*part_nodes)[i]);
				double dx, d2 = 0.0;
				dx = (*p_old++) - (*p_new++);
				d2 += dx*dx;
				dx = (*p_old++) - (*p_new++);
				d2 += dx*dx;
				dx = (*p_old++) - (*p_new++);
				d2 += dx*dx;
				dmax2 = (d2 > dmax2) ? d2 : dmax2;
			}
		}
		else if (nsd == 2)
		{
			for (int i = 0; i < nnd; i++)
			{
				const double* p_new = curr_coords((*part_nodes)[i]);
				double dx, d2 = 0.0;
				dx = (*p_old++) - (*p_new++);
				d2 += dx*dx;
				dx = (*p_old++) - (*p_new++);
				d2 += dx*dx;
				dmax2 = (d2 > dmax2) ? d2 : dmax2;
			}		
		}
		else if (nsd == 1)
		{
			for (int i = 0; i < nnd; i++)
			{
				const double* p_new = curr_coords((*part_nodes)[i]);
				double dx, d2 = 0.0;
				dx = (*p_old++) - (*p_new++);
				d2 += dx*dx;
				dmax2 = (d2 > dmax2) ? d2 : dmax2;
			}		
		}
		else ExceptionT::GeneralFail(caller);
	}
	else /* use ALL nodes */
	{
		int nnd = curr_coords.MajorDim();
		const double *p_new = curr_coords.Pointer();
		if (nsd == 3)
		{
			for (int i = 0; i < nnd; i++)
			{
				double dx, d2 = 0.0;
				dx = (*p_old++) - (*p_new++);
				d2 += dx*dx;
				dx = (*p_old++) - (*p_new++);
				d2 += dx*dx;
				dx = (*p_old++) - (*p_new++);
				d2 += dx*dx;
				dmax2 = (d2 > dmax2) ? d2 : dmax2;
			}
		}
		else if (nsd == 2)
		{
			for (int i = 0; i < nnd; i++)
			{
				double dx, d2 = 0.0;
				dx = (*p_old++) - (*p_new++);
				d2 += dx*dx;
				dx = (*p_old++) - (*p_new++);
				d2 += dx*dx;
				dmax2 = (d2 > dmax2) ? d2 : dmax2;
			}		
		}
		else if (nsd == 1)
		{
			for (int i = 0; i < nnd; i++)
			{
				double dx, d2 = 0.0;
				dx = (*p_old++) - (*p_new++);
				d2 += dx*dx;
				dmax2 = (d2 > dmax2) ? d2 : dmax2;
			}		
		}
		else ExceptionT::GeneralFail(caller);
	}
	
	return sqrt(dmax2);
}

void ParticleT::SetDamping(const ParameterListT& list)
{
	const char caller[] = "ParticleT::SetDamping";

	AutoArrayT<ThermostatBaseT*> thermostats;
	int count = 0;
	const ParameterListT* thermostat_params = list.ListChoice(*this, "thermostats", count);
	while (thermostat_params) {
	
		/* construct new thermostat */
		ThermostatBaseT* thermostat = New_Thermostat(thermostat_params->Name(), true);

		/* initialize */
		thermostat->TakeParameterList(*thermostat_params);
		
		/* store */
		thermostats.Append(thermostat);
		
		/* look for another */	
		thermostat_params = list.ListChoice(*this, "thermostats", ++count);
	}

	/* free existing */
	if (fThermostats.Length() > 0) {
		for (int i = 0; i < fThermostats.Length(); i++)
			delete fThermostats[i];
		fThermostats.Dimension(0);
	}

	/* keep */
	fThermostats.Swap(thermostats);
}


/* describe the parameters needed by the interface */
void ParticleT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ElementBaseT::DefineParameters(list);

	ParameterT lattice_parameter(fLatticeParameter, "lattice_parameter");
	lattice_parameter.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(lattice_parameter);

	ParameterT max_neighbor_distance(fNeighborDistance, "max_neighbor_distance");
	max_neighbor_distance.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(max_neighbor_distance);

	ParameterT re_neighbor_disp(fReNeighborDisp, "re-neighbor_displacement");
	re_neighbor_disp.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(re_neighbor_disp);

	ParameterT re_neighbor_incr(fReNeighborIncr, "re-neighbor_increment");
	re_neighbor_incr.AddLimit(0, LimitT::LowerInclusive);
	re_neighbor_incr.SetDefault(0);
	list.AddParameter(re_neighbor_incr);
}

/* information about subordinate parameter lists */
void ParticleT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ElementBaseT::DefineSubs(sub_list);

	/* periodic boundary conditions */
	sub_list.AddSub("periodic_bc", ParameterListT::Any);

	/* thermostats - array of choices */
	sub_list.AddSub("thermostats", ParameterListT::Any, true);

	/* particle types */
	sub_list.AddSub("particle_type", ParameterListT::OnePlus);
}

/* return the description of the given inline subordinate parameter list */
void ParticleT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "thermostats")
	{
		order = ParameterListT::Choice;
		
		sub_lists.AddSub("velocity_damping");
		sub_lists.AddSub("ramped_damping");
		sub_lists.AddSub("Nose-Hoover");
		sub_lists.AddSub("Gauss_isokinetic");
		sub_lists.AddSub("Langevin");
	}
	else /* inherited */
		ElementBaseT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ParticleT::NewSub(const StringT& name) const
{
	/* try to construct thermostat */
	ThermostatBaseT* thermostat = New_Thermostat(name, false);
	if (thermostat)
		return thermostat;
	else if (name == "periodic_bc")
	{
		ParameterContainerT* pbc = new ParameterContainerT(name);
	
		pbc->AddParameter(ParameterT::Integer, "direction");
		pbc->AddParameter(ParameterT::Double, "x_min");
		pbc->AddParameter(ParameterT::Double, "x_max");
		pbc->AddParameter(ParameterT::Integer, "stretching_schedule", ParameterListT::ZeroOrOnce);

		return pbc;
	}
	else if (name == "particle_type")
	{
		ParameterContainerT* particle_type = new ParameterContainerT(name);
	
		particle_type->AddParameter(ParameterT::Word, "label");

		/* include all particles */
		ParameterT all_particles(ParameterT::Boolean, "all_particles");
		all_particles.SetDefault(true);
		particle_type->AddParameter(all_particles, ParameterListT::ZeroOrOnce);
		
		/* specify by node set IDs */
		particle_type->AddSub("node_ID_list", ParameterListT::ZeroOrOnce);

		/* specify by element block IDs */
		particle_type->AddSub("block_ID_list", ParameterListT::ZeroOrOnce);
		
		return particle_type;
	}
	else /* inherited */
		return ElementBaseT::NewSub(name);
}

/* accept parameter list */
void ParticleT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "ParticleT::TakeParameterList";

	/* inherited */
	ElementBaseT::TakeParameterList(list);

	/* extract parameters */
	fLatticeParameter = list.GetParameter("lattice_parameter");
	fNeighborDistance = list.GetParameter("max_neighbor_distance");
	fReNeighborDisp = list.GetParameter("re-neighbor_displacement");
	fReNeighborIncr = list.GetParameter("re-neighbor_increment");

	/* derived parameters */
	if (NumSD() == 1)
		fNearestNeighborDistance = fLatticeParameter*1.1;
    else if (NumSD() == 2)
		fNearestNeighborDistance = fLatticeParameter*1.1;
    else if (NumSD() == 3)
		fNearestNeighborDistance = fLatticeParameter*.78;
	else
		ExceptionT::GeneralFail(caller);
	fPeriodicSkin = fNeighborDistance;

	/* allocate work space */
	fForce_man.Dimension(ElementSupport().NumNodes(), NumDOF());

	/* periodic boundary conditions */
	fPeriodicBounds.Dimension(NumSD(), 2);
	fPeriodicBounds = 0.0;
	fPeriodicLengths.Dimension(NumSD());
	fPeriodicLengths= 0.0;
	fStretchSchedule.Dimension(NumSD());
	fStretchSchedule = NULL;
	int num_pbc = list.NumLists("periodic_bc");
	if (num_pbc > NumSD())
		ExceptionT::BadInputValue(caller, "expecting at most %d \"periodic_bc\" not %d",
			NumSD(), num_pbc);
	fhas_periodic = (num_pbc > 0) ? 1 : 0;
	for (int i = 0; i < num_pbc; i++) {
	
		/* periodic bc parameters */
		const ParameterListT& pbc = list.GetList("periodic_bc", i);
		
		/* extract parameters */
		int direction = pbc.GetParameter("direction"); direction--;
		double x_min = pbc.GetParameter("x_min");
		double x_max = pbc.GetParameter("x_max");
		const ScheduleT* schedule = NULL;
		const ParameterT* str_sched = pbc.Parameter("stretching_schedule");
		if (str_sched) {
			int sched_num = *str_sched;
			schedule = ElementSupport().Schedule(--sched_num);

			/* check - expecting f(0) = 1 */
			if (fabs(schedule->Value(0.0) - 1.0) > kSmall)
				ExceptionT::BadInputValue(caller, "schedule %d does not have value 1 at time 0", sched_num+1);
		}
		
		/* save */
		fPeriodicLengths[direction] = x_max-x_min;
		fPeriodicBounds(direction,0) = x_min;
		fPeriodicBounds(direction,1) = x_max;
		fStretchSchedule[direction] = schedule;

		/* send to CommManagerT */
		ElementSupport().CommManager().SetPeriodicBoundaries(direction, x_min, x_max);
	}

	/* resolve types */
	fType.Dimension(ElementSupport().NumNodes());
	fType = -1;
	int num_types = list.NumLists("particle_type");
	fTypeNames.Dimension(num_types);
	for (int i = 0; i < num_types; i++) {
		const ParameterListT& particle_type = list.GetList("particle_type", i);
		
		/* label */
		fTypeNames[i] = particle_type.GetParameter("label");
		for (int j = 0; j < i; j++)
			if (fTypeNames[j] == fTypeNames[i])
				ExceptionT::GeneralFail(caller, "types %d and %d have the same label \"%s\"",
					j+1, i+1, fTypeNames[i].Pointer());
					
		/* check "all" particles */
		const ParameterT* all_particles = particle_type.Parameter("all_particles");
		bool all = false;
		if (all_particles) {
			all = *all_particles;
			
			/* only one type can be "all" */
			if (all && num_types > 1)			
				ExceptionT::GeneralFail(caller, "\"all\" particles in \"%s\" conflicts with %d other types",
					fTypeNames[i].Pointer(), num_types-1);
					
			/* mark particles with type 0 */
			fType = 0;
		}
		
		/* look for node list */
		const ParameterListT* node_ID_list = particle_type.List("node_ID_list");
		if (node_ID_list) {
		
			/* collect id's */
			ArrayT<StringT> id_list;
			StringListT::Extract(*node_ID_list, id_list);
			if (all && id_list.Length() > 0)
				ExceptionT::GeneralFail(caller, "label \"%s\" is \"all\" and cannot include %d node ID's",
					fTypeNames[i].Pointer(), id_list.Length());

			/* access to the model database */
			ModelManagerT& model = ElementSupport().ModelManager();
			
			/* get id's */
			iArrayT tags;
			model.ManyNodeSets(id_list, tags);
			
			/* mark map */
			for (int j = 0; j < tags.Length(); j++)
				fType[tags[j]] = i;			
		}

		/* look for block ID list */
		const ParameterListT* block_ID_list = particle_type.List("block_ID_list");
		if (block_ID_list) {
		
			/* collect id's */
			ArrayT<StringT> id_list;
			StringListT::Extract(*block_ID_list, id_list);
			if (all && id_list.Length() > 0)
				ExceptionT::GeneralFail(caller, "label \"%s\" is \"all\" and cannot include %d block ID's",
					fTypeNames[i].Pointer(), id_list.Length());

			/* access to the model database */
			ModelManagerT& model = ElementSupport().ModelManager();
			
			/* read block information */
			for (int k = 0; k < id_list.Length(); k++) {
			
				const iArray2DT& tags = model.ElementGroup(id_list[k]);

				/* mark map - connectivities treated as 1D list */
				for (int j = 0; j < tags.Length(); j++)
					fType[tags[j]] = i;
			}
		}
		
		/* nothing declared */
		if (!all && !node_ID_list && !block_ID_list)
			ExceptionT::GeneralFail(caller, "label \"%s\" is empty", fTypeNames[i].Pointer());
	}

	/* check that all are typed */
	int not_marked_count = 0;
	for (int i = 0; i < fType.Length(); i++)
		if (fType[i] == -1) not_marked_count++;
	if (not_marked_count != 0)
		ExceptionT::BadInputValue(caller, "%d atoms not typed", not_marked_count);
	
	/* extract properties */
	fPropertiesMap.Dimension(num_types);
	fPropertiesMap = -1;
	ExtractProperties(list, fTypeNames, fParticleProperties, fPropertiesMap);

	/* set the neighborlists */
	SetConfiguration();
	
	/* construct thermostats */
	SetDamping(list);
}

/* return a new pair property or NULL if the name is invalid */
ThermostatBaseT* ParticleT::New_Thermostat(const StringT& name, bool throw_on_fail) const
{
	if (name == "velocity_damping")
		return new ThermostatBaseT(ElementSupport());
	else if (name == "ramped_damping")
		return new RampedDampingT(ElementSupport());
	else if (name == "Nose-Hoover")
		return new NoseHooverT(ElementSupport());
	else if (name == "Gauss_isokinetic")
		return new GaussIsokineticT(ElementSupport());
	else if (name == "Langevin")
		return new LangevinT(ElementSupport());
	else if (throw_on_fail)
		ExceptionT::GeneralFail("ParticleT::New_Thermostat",
			"unrecognized thermostat \"%s\"", name.Pointer());
	
	return NULL;
}

void ParticleT::Calc_Slip_and_Strain(dArray2DT &s_values, RaggedArray2DT<int> &RefNearestNeighbors, const int &kEulerLagr)
{
  int non = s_values.MajorDim();
  int num_s_vals = s_values.MinorDim();
  int ndof = NumDOF();
  int num_strains = num_s_vals - ndof - 1;
  iArrayT neighbors;
  dArrayT x_i(ndof), x_j(ndof), r_ij(ndof), R_ij(ndof), X_i(ndof), X_j(ndof);  
  dArrayT slipvector(ndof), svtemp(ndof);
  dMatrixT omega(ndof), eta(ndof), omegatemp(ndof), etatemp(ndof); 
  dMatrixT C_IJ(ndof), b_ij(ndof), F_iI(ndof), strain(ndof);
  dMatrixT etainverse(ndof), Id(ndof);
  int nslip;
  double J, svtol;

  /* set slip vector tolerance for either 1-D (chain), 
   * 2-D (triangular lattice), or 3-D (FCC lattice) systems */
	if (NumSD() == 1)
		svtol = 0.25*fNearestNeighborDistance;
	else if (NumSD() == 2)
		svtol = 0.5*fNearestNeighborDistance;
	else if (NumSD() == 3)
		svtol =  0.5*fNearestNeighborDistance/sqrt(3.0);
	else
		ExceptionT::BadInputValue("ParticleT::Calc_Slip_and_Strain");

  /* multi-processor information */
  CommManagerT& comm_manager = ElementSupport().CommManager();
  const ArrayT<int>* proc_map = comm_manager.ProcessorMap();
  int rank = ElementSupport().Rank();
  const InverseMapT* inverse_map = comm_manager.PartitionNodes_inv();
  const dArray2DT& refcoords = ElementSupport().InitialCoordinates();
  const dArray2DT& coords = ElementSupport().CurrentCoordinates();

  /* row of neighbor list */
  for (int i = 0; i < RefNearestNeighbors.MajorDim(); i++)
  {
   RefNearestNeighbors.RowAlias(i,neighbors);
   int tag_i = neighbors[0];
   int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;

   eta = 0.0; omega = 0.0; strain = 0.0; slipvector = 0.0; nslip = 0;
   J = 1.0;

   coords.RowAlias(tag_i,x_i);
   refcoords.RowAlias(tag_i,X_i);
   for (int j = 1; j<neighbors.Length(); j++)
   {
    int tag_j = neighbors[j];
    coords.RowAlias(tag_j,x_j);
    refcoords.RowAlias(tag_j,X_j);
    r_ij.DiffOf(x_j,x_i);
    R_ij.DiffOf(X_j,X_i);
    svtemp.DiffOf(r_ij,R_ij);
    if (fhas_periodic && (svtemp.Magnitude() > 0.5*fPeriodicLengths.Max()))
    {
     for (int m = 0; m < ndof;m++)
     {
      if (R_ij[m] > 0.5*fPeriodicLengths[m]) R_ij[m] -=fPeriodicLengths[m]; 
      if (R_ij[m] < -0.5*fPeriodicLengths[m]) R_ij[m] +=fPeriodicLengths[m]; 
      if (r_ij[m] > 0.5*fPeriodicLengths[m]) r_ij[m] -=fPeriodicLengths[m]; 
      if (r_ij[m] < -0.5*fPeriodicLengths[m]) r_ij[m] +=fPeriodicLengths[m]; 
     }
     svtemp.DiffOf(r_ij,R_ij);
    }
    slipvector -= svtemp; /* "-" sign is used so that slip is attributed in the correct direction */
    if (svtemp.Magnitude()>svtol) nslip += 1;

    omegatemp.Outer(r_ij,R_ij);
    etatemp.Outer(R_ij,R_ij);
    omega += omegatemp;
    eta += etatemp;
   } /* end of j loop */

   if (nslip>0) slipvector /= double(nslip);

   if (fabs(eta.Det())>kSmall)
   {
    etainverse = eta.Inverse();
    F_iI.MultAB(omega,etainverse);
    if (kEulerLagr)
    {
     b_ij.MultABT(F_iI,F_iI);
     double J2 = b_ij.Det();
     if (fabs(J2)>kSmall)
     {
      J = sqrt(J2);
      Id = 0.0;
      for (int m=0; m<ndof; m++) Id(m,m) = 1.0;
      strain.DiffOf(Id,b_ij.Inverse());
      strain *= 0.5;
     }
    }
    else
    {
     C_IJ.MultATB(F_iI,F_iI); 
     double J2 = C_IJ.Det();
     if (fabs(J2)>kSmall)
     {
      J = sqrt(J2);
      Id = 0.0;
      for (int m=0; m<ndof; m++) Id(m,m) = 1.0;
      strain.DiffOf(C_IJ,Id);
      strain *= 0.5;
     }
    }
   }

   /* put slip vector and strain info into global s_values array */
   int valuep = 0;
   for (int m = 0; m < ndof; m++)
   {
    for (int n = m; n < ndof; n++)
    {
     s_values(local_i,valuep++) = strain(m,n);
    }
   }

   s_values(local_i,valuep++) = J;

   for (int m = 0; m < ndof; m++)
   {
    s_values(local_i,valuep++) = slipvector[m];
   }
   //cout << i << "   " << strain(1,1) << endl;
  } /* end of i loop */
 
} 

int ParticleT::Combination(int n,int k)
  {
   int num_nk = 1; 
   int denom_nk = 1;
   int combo_nk;
   for (int i=(n-k+1); i<=n; i++) num_nk *= i;
   for (int i=1; i<=k; i++) denom_nk *=i;
   combo_nk = num_nk / denom_nk;
   return combo_nk;
  }

void ParticleT::Calc_CSP(dArray2DT &s_values, RaggedArray2DT<int> &NearestNeighbors)
{
	const char caller[] = "ParticleT::Calc_CSP";

  int num_s_vals = s_values.MinorDim();
  int ndof = NumDOF();
  iArrayT neighbors;
  dArrayT x_i(ndof), x_j(ndof), r_ij(ndof), rvec(ndof);  
  int ncspairs;

  /* set number of centrosymmetry pairs to be added up
   * according to the number of spatial dimensions */

  if (NumSD()==1)
   ncspairs = 2; 
  else if (NumSD()==2)
   ncspairs = 3; 
  else if (NumSD()==3)
   ncspairs = 6;
  else
   ExceptionT::BadInputValue(caller);

  /* multi-processor information */
  CommManagerT& comm_manager = ElementSupport().CommManager();
  const InverseMapT* inverse_map = comm_manager.PartitionNodes_inv();
  const dArray2DT& coords = ElementSupport().CurrentCoordinates();

  /* row of neighbor list */
  for (int i = 0; i < NearestNeighbors.MajorDim(); i++)
  {
   NearestNeighbors.RowAlias(i,neighbors);
   int tag_i = neighbors[0];
   int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;

   int nlen = neighbors.Length();
   dArray2DT neighdisp(nlen,ndof);
   int ncombos = Combination(nlen-1,2); 
   dArrayT ndsum(ncombos);

   coords.RowAlias(tag_i,x_i);
   for (int j = 1; j < nlen; j++)
   {
    int tag_j = neighbors[j];
    coords.RowAlias(tag_j,x_j);
    r_ij.DiffOf(x_j,x_i);
    if (fhas_periodic && (r_ij.Magnitude() > 0.5*fPeriodicLengths.Max()))
    {
     for (int m = 0; m < ndof; m++)
     {
      if (r_ij[m] > 0.5*fPeriodicLengths[m]) r_ij[m] -=fPeriodicLengths[m]; 
      if (r_ij[m] < -0.5*fPeriodicLengths[m]) r_ij[m] +=fPeriodicLengths[m]; 
     }
    }
    for (int m = 0; m < ndof; m++)
	{
	 neighdisp(j,m) = r_ij[m];
	}
   }  /* end of j loop */

   int icombos = 0;
   for (int j = 1; j < nlen-1; j++)
   {
	for (int k = j+1; k < nlen; k++)
	{
	  for (int m = 0; m < ndof; m++) 
	  {
	   rvec[m] = neighdisp(j,m) + neighdisp(k,m);
	  } /* end of m loop */
	  ndsum[icombos++] = pow(rvec.Magnitude(),2);
	} /* end of k loop */
   } /* end of j looop */
   if (icombos != ncombos) ExceptionT::SizeMismatch(caller);
   ndsum.SortAscending();
   double csp = 0.0;
   if (ncspairs >= ncombos) ncspairs = ncombos;
   for (int m = 0; m < ncspairs; m++) csp += ndsum[m];
   if (fabs(fLatticeParameter) > kSmall) csp /= pow(fLatticeParameter,2);  

   /* put centrosymmetry parameter into global s_values array */
   s_values(local_i, num_s_vals-1) = csp;
  //cout << i << "   " << csp << endl;
  } /* end of i loop */
 
}

void ParticleT::SetRefNN(RaggedArray2DT<int> &NearestNeighbors,RaggedArray2DT<int> &RefNearestNeighbors)
{
  int ndof = NumDOF();
  iArrayT neighbors;
  CommManagerT& comm_manager = ElementSupport().CommManager();

  /* copy NearestNeighbors list to a ReferenceNN list */
  RefNearestNeighbors = NearestNeighbors;

  const ArrayT<int>* GN = comm_manager.GhostNodes();
  const ArrayT<int>* NWG = comm_manager.NodesWithGhosts();

  if (NWG != NULL)
  {
   InverseMapT map;
   map.SetMap(*GN);
   map.SetOutOfRange(InverseMapT::MinusOne);

   /* looking at all atoms with nearest neighbors that are 
	* image/ghost atoms, replace the identity of those
	* ghost atoms with the ID of the actual atom it is a 
	* ghost of */
   for (int i = 0; i < RefNearestNeighbors.MajorDim(); i++)
   {
    RefNearestNeighbors.RowAlias(i,neighbors);
    for (int j = 1; j<neighbors.Length(); j++)
    {
     int tag_j = neighbors[j];
     int GN_num = map.Map(tag_j);
	 if (GN_num != -1)
	 {
	  int real_num = (*NWG)[GN_num];
	  //cout << i << "   " << tag_j << "   " << real_num << endl;
	  neighbors[j] = real_num;
	 }
    }
	//cout << i << "   " << neighbors.Length() << endl;
   }
  }
}

