/* $Id: ParticleT.cpp,v 1.41 2004-06-28 22:41:39 hspark Exp $ */

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
#include "ExceptionCodes.h"

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
	fDmax(0),
	fForce_man(0, fForce, field.NumDOF()),


	fActiveParticles(NULL),
	fRandom(NULL)
	
{
	SetName("particle");

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

	fPeriodicSkin = fNeighborDistance;
}

ParticleT::ParticleT(const ElementSupportT& support):
	ElementBaseT(support),
	fNeighborDistance(-1),
	fReNeighborDisp(-1),
	fReNeighborIncr(-1),
	fPeriodicSkin(-1),
	fNumTypes(-1),
	fGrid(NULL),
	fReNeighborCounter(0),
	fDmax(0),
	fActiveParticles(NULL),
	fRandom(NULL)
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

	fPeriodicBounds.Dimension(NumSD(), 2);
	fPeriodicBounds = 0.0;
	fPeriodicLengths.Dimension(NumSD());
	fPeriodicLengths=0.0;
	fStretchSchedule.Dimension(NumSD());
	fStretchSchedule = NULL;

	out << " Periodic boundary conditions:\n";
	out << setw(kIntWidth) << "dir"
	    << setw(d_width) << "min"
	    << setw(d_width) << "max" << '\n';
	ifstreamT& in = ElementSupport().Input();
	fhas_periodic=0;

	for (int i = 0; i < NumSD(); i++) {
	
		out << setw(kIntWidth) << i+1;

		 int has_periodic = 0;
		in >> has_periodic;
		
		if (has_periodic > 0) {
		  fhas_periodic= fhas_periodic|has_periodic;
			double x_min = 0.0, x_max = 0.0;
			in >> x_min >> x_max;
			out << setw(d_width) << x_min << setw(d_width) << x_max << '\n';
			if (x_min > x_max)
				ExceptionT::BadInputValue(caller, "x_min > x_max: %g < %g", x_min, x_max);

			fPeriodicLengths[i]=x_max-x_min;
			/* store */
			fPeriodicBounds(i,0) = x_min;
			fPeriodicBounds(i,1) = x_max;


			/* send to CommManagerT */
			ElementSupport().CommManager().SetPeriodicBoundaries(i, x_min, x_max);

			
			/* read stretch schedule */
			if (has_periodic > 1) {
				int schedule = -99;
				in >> schedule;
				schedule--;
				fStretchSchedule[i] = ElementSupport().Schedule(schedule);
				
				/* check - expecting f(0) = 1 */
				if (fabs(fStretchSchedule[i]->Value(0.0) - 1.0) > kSmall)
					ExceptionT::BadInputValue(caller, "schedule %d does not have value 1 at time 0", schedule+1);
			}

		}
		else out << setw(d_width) << "-" << setw(d_width) << "-" << '\n';
	}
	
	/* read properties information */
	EchoProperties(in, out);


	StringT key;
	in>>key;
	if (key =="lattice") in >> fLatticeParameter;
	else
	  {
	    in.rewind();
	    fLatticeParameter = 0.0; //defaults to nothing
	  }
	if (NumSD()==1)
	{
	 NearestNeighborDistance=fLatticeParameter*1.1;
	} 
    else if (NumSD()==2)
	{
	 NearestNeighborDistance=fLatticeParameter*1.1;
    }
    else if (NumSD()==3)
    { 
	 NearestNeighborDistance=fLatticeParameter*.78;
    }
	else
	{ 
	 throw eBadInputValue;
	}

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
		(fReNeighborIncr != -1 && fReNeighborCounter >= fReNeighborIncr))
	{
		/* output stream */
		ofstreamT& out = ElementSupport().Output();
		if (fReNeighborDisp > 0.0 && fDmax > fReNeighborDisp)
			out << "\n ParticleT::RelaxSystem: max displacement since re-neighboring "
			    << fDmax << " > " << fReNeighborDisp << '\n';
		if (fReNeighborIncr != -1 && fReNeighborCounter >= fReNeighborIncr)
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
	return ElementSupport().CommManager().PartitionNodesChanging();
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
	fThermostats = NULL;
	
	CommManagerT& comm_manager = ElementSupport().CommManager();
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
		char peekahead = in.next_char();
		int nodesOrRegion = atoi(&peekahead);
		switch (nodesOrRegion)
		{
			case ThermostatBaseT::kNodes:
			{
				in >> nodesOrRegion; // read in what I didn't above
				if (thermostat_i == ThermostatBaseT::kRampedDamping)
					ExceptionT::BadInputValue(caller,"Ramped Damping requires spatial region");
				
				fThermostats[i]->InitNodeSets(in, ElementSupport().Model());
				
				break;
			}
			case ThermostatBaseT::kRegion:
			{
				int nregions;
				in >> nregions;
				fThermostats[i]->InitRegion(in, ElementSupport().InitialCoordinates(),
											comm_manager.PartitionNodes());
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
			LangevinT* tmpPtr = TB_DYNAMIC_CAST(LangevinT*, fThermostats[i]); 
			
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


/* describe the parameters needed by the interface */
void ParticleT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ElementBaseT::DefineParameters(list);

	ParameterT neighbor_distance(fNeighborDistance, "neighbor_distance");
	neighbor_distance.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(neighbor_distance);

	ParameterT re_neighbor_disp(fReNeighborDisp, "re_neighbor_displacement");
	re_neighbor_disp.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(re_neighbor_disp, ParameterListT::ZeroOrOnce);



	ParameterT re_neighbor_incr(fReNeighborIncr, "re_neighbor_increment");
	re_neighbor_incr.AddLimit(0, LimitT::Lower);
	re_neighbor_incr.SetDefault(1);
	list.AddParameter(re_neighbor_incr, ParameterListT::ZeroOrOnce);
}

/* information about subordinate parameter lists */
void ParticleT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ElementBaseT::DefineSubs(sub_list);

	/* thermostats - array of choices */
	sub_list.AddSub("thermostats", ParameterListT::Any, true);
}

/* return the description of the given inline subordinate parameter list */
void ParticleT::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
	SubListT& sub_sub_list) const
{
	if (sub == "thermostats")
	{
		order = ParameterListT::Choice;
		
		sub_sub_list.AddSub("ramped_damping");
		sub_sub_list.AddSub("NoseHoover");
		sub_sub_list.AddSub("Gauss_isokinetic");
		sub_sub_list.AddSub("Langevin");
	}
	else /* inherited */
		ElementBaseT::DefineInlineSub(sub, order, sub_sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ParticleT::NewSub(const StringT& list_name) const
{
	/* try to construct thermostat */
	ThermostatBaseT* thermostat = New_Thermostat(list_name, false);
	if (thermostat)
		return thermostat;
	else /* inherited */
		return ElementBaseT::NewSub(list_name);
}

/* return a new pair property or NULL if the name is invalid */
ThermostatBaseT* ParticleT::New_Thermostat(const StringT& name, bool throw_on_fail) const
{
	if (name == "ramped_damping")
		return new RampedDampingT;
	else if (name == "NoseHoover")
		return new NoseHooverT;
	else if (name == "Gauss_isokinetic")
		return new GaussIsokineticT;
	else if (name == "Langevin")
		return new LangevinT;
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
  if (NumSD()==1)
  {
   svtol = 0.25*NearestNeighborDistance;
  }
  else if (NumSD()==2)
  {
   svtol = 0.5*NearestNeighborDistance;
  }
  else if (NumSD()==3)
  {
   svtol =  0.5*NearestNeighborDistance/sqrt(3.0);
  }
  else
  {
   throw eBadInputValue;
  } 
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
  int non = s_values.MajorDim();
  int num_s_vals = s_values.MinorDim();
  int ndof = NumDOF();
  iArrayT neighbors;
  dArrayT x_i(ndof), x_j(ndof), r_ij(ndof), rvec(ndof);  
  int ncspairs;

  /* set number of centrosymmetry pairs to be added up
   * according to the number of spatial dimensions */

  if (NumSD()==1)
  {
   ncspairs = 2; 
  }
  else if (NumSD()==2)
  {
   ncspairs = 3; 
  }
  else if (NumSD()==3)
  {
   ncspairs = 6;
  }
  else
  {
   throw eBadInputValue;
  }

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
   if (icombos != ncombos) throw eSizeMismatch;
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

