/* $Id: ParticleT.cpp,v 1.10 2002-12-04 06:34:18 paklein Exp $ */
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

using namespace Tahoe;

/* class parameters */
const int kAvgNodesPerCell = 20;
const int kMaxNumCells     =- 1; /* -1: no max */

/* constructors */
ParticleT::ParticleT(const ElementSupportT& support, const FieldT& field):
	ElementBaseT(support, field),
	fReNeighborIncr(-1),
	fNumTypes(-1),
	fGrid(NULL),
	fReNeighborCounter(0)
{
	/* set matrix format */
	fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);

	/* read class parameters */
	ifstreamT& in = ElementSupport().Input();
	in >> fReNeighborIncr;
	if (fReNeighborIncr < 0) ExceptionT::GeneralFail();
}

/* destructor */
ParticleT::~ParticleT(void)
{
	/* free search grid */
	delete fGrid;
}

/* initialization */
void ParticleT::Initialize(void)
{
	const char caller[] = "ParticleT::Initialize";

	/* inherited */
	ElementBaseT::Initialize();
	
	/* allocate work space */
	fForce.Dimension(ElementSupport().NumNodes(), NumDOF());
	
	/* read properties information */
	ifstreamT& in = ElementSupport().Input();
	ofstreamT& out = ElementSupport().Output();
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
	OutputSetT output_set(GeometryT::kPoint, fPointConnectivities, n_labels, true);
		
	/* register and get output ID */
	fOutputID = ElementSupport().RegisterOutput(output_set);
}

void ParticleT::WriteOutput(void)
{
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

	/* generate contact element data */
	fReNeighborCounter++;
	if (fReNeighborIncr != 0 && fReNeighborCounter >= fReNeighborIncr)
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
	/* assume that for a single processor calculation, the geometry
	 * is not changing */
	if (ElementSupport().Size() > 1)
		return true;
	else
		return false;
}

/* echo element connectivity data */
void ParticleT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	const char caller[] = "ParticleT::EchoConnectivityData";

	/* access to the model database */
	ModelManagerT& model = ElementSupport().Model();

	int all_or_some = -1;
	in >> all_or_some; 
	if (all_or_some != 0 && all_or_some != 1) ExceptionT::BadInputValue(caller);
	
	if (all_or_some == 0) /* ALL */
	{
		fID.Dimension(1);
		fID[0] = "ALL";
		
		/* make list of all nodes */
		fGlobalTag.Dimension(ElementSupport().NumNodes());
		fGlobalTag.SetValueToPosition();
	}
	else
	{
		/* read node set indexes */
		model.NodeSetList(in, fID);

		/* get tags */
		iArrayT tags;
		if (fID.Length() > 0) model.ManyNodeSets (fID, tags);
		
		/* remove duplicates */
		fGlobalTag.Union(tags);
	}

	/* "point connectivities" needed for output */
	fPointConnectivities.Set(fGlobalTag.Length(), 1, fGlobalTag.Pointer());
	
	/* set inverse map */
	fGlobalToLocal.SetMap(fGlobalTag);
	
	/* write connectivity info */
	out << " Number of particles . . . . . . . . . . . . . . = " << fGlobalTag.Length();
	if (all_or_some == 0) out << " (ALL)";
	out << endl;
	
	/* particle types */
	fNumTypes = -1;
	in >> fNumTypes;
	if (fNumTypes < 1) ExceptionT::BadInputValue(caller, "must define at least one type");
	
	/* initialize type map */
	fType.Dimension(ElementSupport().NumNodes());
	fType = -1;

	in >> all_or_some; 
	if (all_or_some != 0 && all_or_some != 1) ExceptionT::BadInputValue(caller);
	if (fNumTypes > 1 && all_or_some == 0)
		ExceptionT::BadInputValue(caller, "atom types must be listed explicitly if there is more than 1 type");

	if (all_or_some == 0) /* ALL */
	{
		/* mark particle tags with type 0 */
		for (int i = 0; i < fGlobalTag.Length(); i++)
			fType[fGlobalTag[i]] = 0;
	}
	else
	{
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
	for (int i = 0; i < fGlobalTag.Length(); i++)
		if (fType[fGlobalTag[i]] == -1) not_marked_count++;
	if (not_marked_count != 0)
		ExceptionT::BadInputValue(caller, "%d atoms not typed", not_marked_count);
}

/* generate neighborlist */
void ParticleT::GenerateNeighborList(const iArrayT& particle_tags, 
	double distance, bool double_list, RaggedArray2DT<int>& neighbors)
{
	/* global coordinates */
	const dArray2DT& coords = ElementSupport().CurrentCoordinates();

	/* construct grid */
	if (!fGrid) fGrid = new iGridManagerT(kAvgNodesPerCell, kMaxNumCells, coords, &particle_tags);

	/* reset contents */
	fGrid->Reset();
	
	/* set up temp space */
	int init_num_neighbors = 6;
	int num_tags = particle_tags.Length();
	int num_chunks = num_tags/250;
	num_chunks = (num_chunks < 1) ? 1 : num_chunks; 
	AutoFill2DT<int> auto_neighbors(num_tags, num_chunks, 20, init_num_neighbors);
	
	/* loop over tags */
	int nsd = coords.MinorDim();
	double distance2 = distance*distance;
	for (int i = 0; i < particle_tags.Length(); i++)
	{
		/* this tag */
		int tag_i = particle_tags[i];
	
		/* add self */
		auto_neighbors.Append(i, tag_i);
		
		/* gets points from grid */
		const double* coords_i = coords(tag_i);
		const AutoArrayT<iNodeT>& hits = fGrid->HitsInRegion(coords_i, distance);
		
		/* filter neighbors */
		for (int j = 0; j < hits.Length(); j++)
		{
			int tag_j = hits[j].Tag();
			
			if (tag_j > tag_i || double_list)
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
	/* loop over particles */
	fForce = 0.0;
	for (int i = 0; i < fGlobalTag.Length(); i++)
		/* assemble into global array */
		fForce.SetRow(fGlobalTag[i], mass[fType[i]]);
	
	/* assemble all */
	ElementSupport().AssembleLHS(Group(), fForce, Field().Equations());
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
