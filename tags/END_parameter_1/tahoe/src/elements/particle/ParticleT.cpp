/* $Id: ParticleT.cpp,v 1.18 2003-04-22 01:23:15 cjkimme Exp $ */
#include "ParticleT.h"

#include "fstreamT.h"
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
ParticleT::ParticleT(const ElementSupportT& support, const FieldT& field):
	ElementBaseT(support, field),
	fNeighborDistance(-1),
	fReNeighborDisp(-1),
	fReNeighborIncr(-1),
	fNumTypes(-1),
	fGrid(NULL),
	fReNeighborCounter(0),
	fCommManager(support.CommManager()),
	fDmax(0),
	fForce_man(0, fForce, field.NumDOF()),
	fActiveParticles(NULL),
	fRandom(NULL)
{
	/* set matrix format */
	fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);

	/* read class parameters */
	ifstreamT& in = ElementSupport().Input();
	in >> fNeighborDistance
	   >> fReNeighborDisp
	   >> fReNeighborIncr;

	if (fNeighborDistance < 0) ExceptionT::GeneralFail();
	
	/* values < 0 mean ignore */
	fReNeighborDisp = (fReNeighborDisp < kSmall) ? -1 : fReNeighborDisp;
	fReNeighborIncr = (fReNeighborIncr <= 0) ? -1 : fReNeighborIncr;
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
		
	if (fRandom)
		delete fRandom;
	
}

/* initialization */
void ParticleT::Initialize(void)
{
	const char caller[] = "ParticleT::Initialize";

	/* inherited */
	ElementBaseT::Initialize();
	
	/* allocate work space */
	fForce_man.SetMajorDimension(ElementSupport().NumNodes(), false);

	/* write parameters */
	ofstreamT& out = ElementSupport().Output();
	int d_width = OutputWidth(out, &fNeighborDistance);
	out << " Neighbor cut-off distance . . . . . . . . . . . = " << fNeighborDistance << '\n';
	out << " Re-neighboring displacement trigger . . . . . . = " << fReNeighborDisp << '\n';
	out << " Re-neighboring interval . . . . . . . . . . . . = " << fReNeighborIncr << '\n';

	/* periodic boundary conditions */
	out << " Periodic boundary conditions:\n";
	out << setw(kIntWidth) << "dir"
	    << setw(d_width) << "min"
	    << setw(d_width) << "max" << '\n';
	ifstreamT& in = ElementSupport().Input();
	for (int i = 0; i < NumSD(); i++) {
	
		out << setw(kIntWidth) << i+1;
		int has_periodic = 0;
		in >> has_periodic;
		if (has_periodic) {
			double x_min = 0.0, x_max = 0.0;
			in >> x_min >> x_max;
			out << setw(d_width) << x_min << setw(d_width) << x_max << '\n';
			if (x_min > x_max)
				ExceptionT::BadInputValue(caller, "x_min > x_max: %g < %g", x_min, x_max);

			/* send to CommManagerT */
			ElementSupport().CommManager().SetPeriodicBoundaries(i, x_min, x_max);
		}
		else out << setw(d_width) << "-" << setw(d_width) << "-" << '\n';
	}
	
	/* read properties information */
	EchoProperties(in, out);

	/* set up communication of type information */
	fTypeMessageID = ElementSupport().CommManager().Init_AllGather(MessageT::Integer, 1);

	/* map of {type_a, type_b} -> potential number */
	fPropertiesMap.Dimension(fNumTypes);
	fPropertiesMap = -1;
	for (int i = 0; i < fPropertiesMap.Length(); i++)
	{
		int a, b, pot_num;
		in >> a >> b >> pot_num;
		a--; b--; /* offset */
		if (fPropertiesMap(a,b) != -1) 
			ExceptionT::BadInputValue(caller, "potential map(%d,%d) is already assigned to %d", a, b, fPropertiesMap(a,b));
		else if (pot_num < 1 && pot_num > fNumTypes)
			ExceptionT::BadInputValue(caller, "potential %d is out of range (%d,%d)", pot_num, 1, fNumTypes);
		else
			fPropertiesMap(a,b) = pot_num;
	}
	fPropertiesMap--; /* offset */

	/* set the neighborlists */
	SetConfiguration();
	
	EchoDamping(in, out);
	
}

/* form of tangent matrix */
GlobalT::SystemTypeT ParticleT::TangentType(void) const
{
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
	const ArrayT<int>* parition_nodes = fCommManager.PartitionNodes();
	if (parition_nodes)
	{
		int num_nodes = parition_nodes->Length();
		fPointConnectivities.Set(num_nodes, 1, parition_nodes->Pointer());
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
	
	/* reset connectivities */
	if (ChangingGeometry())
	{
		const ArrayT<int>* parition_nodes = fCommManager.PartitionNodes();
		if (parition_nodes)
		{
			int num_nodes = parition_nodes->Length();
			fPointConnectivities.Set(num_nodes, 1, parition_nodes->Pointer());	
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

	/* compute max distance traveled since last neighboring 
	 * (across all processes) */
	fDmax = fCommManager.Communicator().Max(MaxDisplacement());

	/* check damping regions */
	//fDampingCounters++;

	/* generate contact element data */
	fReNeighborCounter++;
	if ((fReNeighborDisp > 0.0 && fDmax > fReNeighborDisp) || 
		(fReNeighborIncr != -1 && fReNeighborCounter >= fReNeighborIncr))
	{
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
	
	if (fRandom != NULL)
		out << fRandom->RandSeed() << '\n';
	
	for (int i = 0; i < nThermostats; i++)
		fThermostats[i]->WriteRestart(out);
}

/* read restart data to the output stream */
void ParticleT::ReadRestart(istream& in)
{
	/* read counter */
	in >> fReNeighborCounter;
	
	if (fRandom != NULL)
	{
		long seed;
		in >> seed;
		fRandom->sRand(seed);
	}
	
	for (int i = 0; i < nThermostats; i++)
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
		const ArrayT<int>* part_nodes = fCommManager.PartitionNodes();
		if (part_nodes)
		{
			nodes_used = 0;
			int npn = part_nodes->Length();
			int*  p = part_nodes->Pointer();
			for (int i = 0; i < npn; i++)
				nodes_used[*p++] = 1;
		}
		else /* all are partition nodes */
			nodes_used = 1;

		/* mark nodes to skip */
		int nsn = skip.Length();
		int* ps = skip.Pointer();
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
	comm_manager.EnforcePeriodicBoundaries(fNeighborDistance);
	
	/* reset the types array */
	int nnd = ElementSupport().NumNodes();
	fType.Resize(nnd);
	
	/* exchange type information */
	iArray2DT type_wrapper(fType.Length(), 1, fType.Pointer());
	comm_manager.AllGather(fTypeMessageID, type_wrapper);
	
	/* resize working arrays */
	fForce_man.SetMajorDimension(nnd, false);

	/* collect current coordinates */
	const ArrayT<int>* part_nodes = fCommManager.PartitionNodes();
	const dArray2DT& curr_coords = ElementSupport().CurrentCoordinates();
	if (part_nodes)
	{
		/* collect */
		fReNeighborCoords.Dimension(part_nodes->Length(), curr_coords.MinorDim());
		fReNeighborCoords.RowCollect(*part_nodes, curr_coords);
	}
	else /* use ALL nodes */
		fReNeighborCoords = curr_coords;

	/* write the search grid statistics */
	if (fGrid) {
		ofstreamT& out = ElementSupport().Output();
		fGrid->WriteStatistics(out);
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

void ParticleT::ApplyDamping(const RaggedArray2DT<int>& fNeighbors)
{		
	if (QisDamped)
	{
		const dArray2DT* velocities = NULL; 
     	if (Field().Order() > 0) // got velocities!
     	{
     		velocities = &(Field()[1]);
     		
     		for (int i = 0; i < nThermostats; i++)
				fThermostats[i]->ApplyDamping(fNeighbors,velocities,fForce,
										fType,fParticleProperties);
		}
	}
		
}

/* return true if connectivities are changing */
bool ParticleT::ChangingGeometry(void) const
{
	return fCommManager.PartitionNodesChanging();
}

/* echo element connectivity data */
void ParticleT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
#pragma unused(out)
	const char caller[] = "ParticleT::EchoConnectivityData";
	
	/* particle types */
	fNumTypes = -1;
	in >> fNumTypes;
	if (fNumTypes < 1) ExceptionT::BadInputValue(caller, "must define at least one type");
	
	/* initialize type map */
	fType.Dimension(ElementSupport().NumNodes());
	fType = -1;

	int all_or_some = -99;
	in >> all_or_some; 
	if (all_or_some != 0 && all_or_some != 1) ExceptionT::BadInputValue(caller);
	if (fNumTypes > 1 && all_or_some == 0)
		ExceptionT::BadInputValue(caller, "atom types must be listed explicitly if there is more than 1 type");

	if (all_or_some == 0) /* ALL */
	{
		/* mark particle tags with type 0 */
		fType = 0;
	}
	else
	{
		/* access to the model database */
		ModelManagerT& model = ElementSupport().Model();

		/* read sets */
		for (int i = 0; i < fNumTypes; i++)
		{
			/* read node set ids */
			ArrayT<StringT> ids;
			model.NodeSetList(in, ids);
			iArrayT tags;
			model.ManyNodeSets(ids, tags);
	
			/* mark map */
			for (int j = 0; j < tags.Length(); j++)
				fType[tags[j]] = i;
		}
	}

	/* check that all are typed */
	int not_marked_count = 0;
	for (int i = 0; i < fType.Length(); i++)
		if (fType[i] == -1) not_marked_count++;
	if (not_marked_count != 0)
		ExceptionT::BadInputValue(caller, "%d atoms not typed", not_marked_count);
}

/* generate neighborlist */
void ParticleT::GenerateNeighborList(const ArrayT<int>* particle_tags, 
	double distance, RaggedArray2DT<int>& neighbors, 
	bool double_list, bool full_list)
{
	/* global coordinates */
	const dArray2DT& coords = ElementSupport().CurrentCoordinates();

	/* node to processor map */
	const ArrayT<int>* n2p_map = fCommManager.ProcessorMap();

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
	const ArrayT<int>* partition_nodes = fCommManager.PartitionNodes();
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
	const ArrayT<int>* part_nodes = fCommManager.PartitionNodes();
	
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
	const ArrayT<int>* part_nodes = fCommManager.PartitionNodes();
	const dArray2DT& curr_coords = ElementSupport().CurrentCoordinates();
	double dmax2 = 0.0;
	int nsd = curr_coords.MinorDim();
	double *p_old = fReNeighborCoords.Pointer();
	if (part_nodes)
	{
		int nnd = part_nodes->Length();
		if (nnd != fReNeighborCoords.MajorDim()) ExceptionT::SizeMismatch(caller);
		if (nsd == 3)
		{
			for (int i = 0; i < nnd; i++)
			{
				double* p_new = curr_coords((*part_nodes)[i]);
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
				double* p_new = curr_coords((*part_nodes)[i]);
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
				double* p_new = curr_coords((*part_nodes)[i]);
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
		double *p_new = curr_coords.Pointer();
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

void ParticleT::EchoDamping(ifstreamT& in, ofstreamT& out)
{
#pragma unused(out)
	
	const char caller[] = "ParticleT::EchoDamping";

	/* default flag assuming no damping/thermostatting present */
	QisDamped = false;
	
	/* flag for constructing a random number generator */
	bool QisRandom = false;
	
	in >> nThermostats;
	fThermostats.Dimension(nThermostats);
	
	for (int i = 0; i < nThermostats; i++)
	{
		bool QisLangevin = false;
	
		ThermostatBaseT::ThermostatT thermostat_i;
		in >> thermostat_i;
		switch (thermostat_i)
		{
			case ThermostatBaseT::kDamped:
			{
				QisDamped = true;
				
				fThermostats[i] = new ThermostatBaseT(in,
					ElementSupport().NumSD(),ElementSupport().TimeStep());
				
				break;
			}
			case ThermostatBaseT::kLangevin:
			{
				QisRandom = QisLangevin = true;
				QisDamped = true;
				
				fThermostats[i] = new LangevinT(in,
					ElementSupport().NumSD(),ElementSupport().TimeStep());
				
				break;
			}
			case ThermostatBaseT::kNoseHoover:
			{
				QisDamped = true;
				
				fThermostats[i] = new NoseHooverT(in,
					ElementSupport().NumSD(), ElementSupport().TimeStep());
				
				break;
			}
			case ThermostatBaseT::kGaussIsokinetic:
			{
				QisDamped = true;
				
				fThermostats[i] = new GaussIsokineticT(in,
					ElementSupport().NumSD(), ElementSupport().TimeStep());
				
				break;
			}
			case ThermostatBaseT::kRampedDamping:
			{
				QisDamped = true;
				
				fThermostats[i] = new RampedDampingT(in,
					ElementSupport().NumSD(), ElementSupport().TimeStep());
				
				break;
			}
			default:
			{
				ExceptionT::BadInputValue(caller,"Damping type does not exist or is not valid");
			}
		}
		int nodesOrRegion;
		in >> nodesOrRegion;
		switch (nodesOrRegion)
		{
			case ThermostatBaseT::kNodes:
			{
				if (thermostat_i == ThermostatBaseT::kRampedDamping)
					ExceptionT::BadInputValue(caller,"Ramped Damping requires spatial region");
				int all_or_some = -99;
				in >> all_or_some; 
				if (all_or_some != 0 && all_or_some != 1) ExceptionT::BadInputValue(caller);
		
				if (all_or_some == 0) /* ALL */
				{
					if (nThermostats != 1)
						ExceptionT::GeneralFail(caller,"For all particles damped, must have only 1 controller"); 
	 			
				}
				else
				{
					/* access to the model database */
					ModelManagerT& model = ElementSupport().Model();

					/* read node set ids */
					ArrayT<StringT> ids;
					model.NodeSetList(in, ids);
					iArrayT tags;
					model.ManyNodeSets(ids, fThermostats[i]->NodeList());	
				}
				
				break;
			}
			case ThermostatBaseT::kRegion:
			{
				ExceptionT::BadInputValue(caller,"Thermostatted region not implemented yet");
				break;
			}
			default:
			{
				ExceptionT::BadInputValue(caller,"Thermostat control type invalid");
			}
		}
		if (thermostat_i != ThermostatBaseT::kDamped)
		{
			int schedNum;
			double schedVal;
			in >> schedNum;
			schedNum--;
			in >> schedVal;
			const ScheduleT* sched = ElementSupport().Schedule(schedNum);
			if (!sched)
				ExceptionT::GeneralFail(caller,"Unable to get temperature schedule");
			fThermostats[i]->SetTemperatureSchedule(sched,schedVal);
		}
		if (QisLangevin)
		{
			if (fRandom == 0) // construct a singleton random number generator
			{	
				fRandom = new RandomNumberT(RandomNumberT::kParadynGaussian);
			}
			LangevinT* tmpPtr = dynamic_cast<LangevinT*>(fThermostats[i]); 
			
			if (!tmpPtr)
				ExceptionT::GeneralFail(caller,"Cannot send random number gen to thermostat");
			tmpPtr->SetRandNumGenerator(fRandom);
		}
	}
	
	if (QisRandom)
	{
		int randSeed;
		in >> randSeed;
		if (randSeed < -1) ExceptionT::BadInputValue(caller,"random seed must be >= -1");
		if (randSeed == -1) 
		{
			ExceptionT::GeneralFail(caller,"No machine-generated random seed available yet");
		}
		fRandom->sRand(randSeed);
	}
}


