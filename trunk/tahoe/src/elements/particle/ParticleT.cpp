/* $Id: ParticleT.cpp,v 1.5 2002-11-22 01:49:45 paklein Exp $ */
#include "ParticleT.h"

#include "fstreamT.h"
#include "eControllerT.h"
#include "OutputSetT.h"
#include "dArray2DT.h"
#include "ElementSupportT.h"
#include "ModelManagerT.h"
#include "iGridManagerT.h"
#include "iNodeT.h"
#include "PotentialT.h"
#include "RaggedArray2DT.h"

using namespace Tahoe;

/* class parameters */
const int kAvgNodesPerCell = 20;
const int kMaxNumCells     =- 1; /* -1: no max */

/* constructors */
ParticleT::ParticleT(const ElementSupportT& support, const FieldT& field):
	ElementBaseT(support, field),
	fGrid(NULL)
{
	/* set matrix format */
	fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
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
	/* inherited */
	ElementBaseT::Initialize();
	
	/* dimension work space */
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
//TEMP
ExceptionT::Stop("ParticleT::WriteOutput", "not implemented");

	/* get list of nodes used by the group */
	iArrayT nodes_used;
	NodesUsed(nodes_used);

	/* temp space for group displacements */
	dArray2DT disp(nodes_used.Length(), NumDOF());
	
	/* collect group displacements */
	disp.RowCollect(nodes_used, Field()[0]);

	/* send */
	dArray2DT e_values;
	ElementSupport().WriteOutput(fOutputID, disp, e_values);
	
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
	int all_or_some = -1;
	in >> all_or_some; 
	if (all_or_some != 0 && all_or_some != 1) ExceptionT::BadInputValue("ParticleT::EchoConnectivityData");
	
	if (all_or_some == 0)
	{
		fID.Dimension(1);
		fID[0] = "ALL";
		
		/* make list of all nodes */
		fGlobalTag.Dimension(ElementSupport().NumNodes());
		fGlobalTag.SetValueToPosition();
	}
	else
	{
		/* access to the model database */
		ModelManagerT& model = ElementSupport().Model();

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
	
	/* write connectivity info */
	out << " Number of particles . . . . . . . . . . . . . . = " << fGlobalTag.Length();
	if (all_or_some == 0) out << " (ALL)";
	out << endl;
}

/* generate labels for output data */
void ParticleT::GenerateOutputLabels(ArrayT<StringT>& labels) const
{
	if (NumDOF() > 0) ExceptionT::GeneralFail("ParticleT::GenerateOutputLabels");

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

/* generate neighborlist */
void ParticleT::GenerateNeighborList(const iArrayT& particle_tags, 
	double distance, bool double_list, RaggedArray2DT<int>& neighbors)
{
	/* the global coordinate array */
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
