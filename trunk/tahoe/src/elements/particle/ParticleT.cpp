/* $Id: ParticleT.cpp,v 1.4 2002-11-21 01:11:14 paklein Exp $ */
#include "ParticleT.h"

#include "fstreamT.h"
#include "eControllerT.h"
#include "OutputSetT.h"
#include "dArray2DT.h"
#include "ElementSupportT.h"
#include "ModelManagerT.h"
#include "iGridManagerT.h"
#include "PotentialT.h"

using namespace Tahoe;

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
	const iArrayT& particle_type, const ArrayT<PotentialT*>& pots, 
	RaggedArray2DT<int>& neighbors, bool double_list)
{
	/* construct grid */
	if (!fGrid) {
	
	
	}

	/* reset contents */
	
	

}
