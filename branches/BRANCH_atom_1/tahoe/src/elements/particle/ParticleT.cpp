/* $Id: ParticleT.cpp,v 1.10.2.4 2003-01-11 01:15:25 paklein Exp $ */
#include "ParticleT.h"

#include "fstreamT.h"
#include "eControllerT.h"
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

using namespace Tahoe;

/* class parameters */
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
	fDmax(0)
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
}

/* initialization */
void ParticleT::Initialize(void)
{
	const char caller[] = "ParticleT::Initialize";

	/* inherited */
	ElementBaseT::Initialize();
	
	/* allocate work space */
	fForce.Dimension(ElementSupport().NumNodes(), NumDOF());

	/* write parameters */
	ofstreamT& out = ElementSupport().Output();
	out << " Neighbor cut-off distance . . . . . . . . . . . = " << fNeighborDistance << '\n';
	out << " Re-neighboring displacement trigger . . . . . . = " << fReNeighborDisp << '\n';
	out << " Re-neighboring interval . . . . . . . . . . . . = " << fReNeighborIncr << '\n';
	
	/* read properties information */
	ifstreamT& in = ElementSupport().Input();
	EchoProperties(in, out);

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

	/* write the search grid statistics */
	if (fGrid) {
		ofstreamT& out = ElementSupport().Output();
		fGrid->WriteStatistics(out);
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
}

/* read restart data to the output stream */
void ParticleT::ReadRestart(istream& in)
{
	/* read counter */
	in >> fReNeighborCounter;

	/* (re-)set configuration */
	SetConfiguration();	
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* return true if connectivities are changing */
bool ParticleT::ChangingGeometry(void) const
{
	return fCommManager.PartitionNodesChanging();
}

/* set neighborlists */
void ParticleT::SetConfiguration(void)
{
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
			
			if (double_list || tag_j > tag_i || (full_list && pr_i != pr_j && tag_j != tag_i))
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
				if (d2 <= distance2) auto_neighbors.Append(i, tag_j);
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

namespace Tahoe {

/** stream extraction operator */
istream& operator>>(istream& in, ParticleT::PropertyT& property)
{
	int i_property;
	in >> i_property;
	switch (i_property)
	{
		case ParticleT::kHarmonicPair:
			property = ParticleT::kHarmonicPair;
			break;
		case ParticleT::kLennardJonesPair:
			property = ParticleT::kLennardJonesPair;
			break;
		case ParticleT::kParadynPair:
			property = ParticleT::kParadynPair;
			break;
		default:
			ExceptionT::BadInputValue("operator>>ParticleT::PropertyT", 
				"unknown code: %d", i_property);
	}
	return in;
}

} /* namespace Tahoe */
