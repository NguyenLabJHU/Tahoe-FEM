/* $Id: SCNIMFT.cpp,v 1.17.2.2 2004-07-08 00:41:53 paklein Exp $ */
#include "SCNIMFT.h"

//#define VERIFY_B

#include "ArrayT.h"
#include "ofstreamT.h"
#include "ifstreamT.h"
#include "fstreamT.h"
#include "eIntegratorT.h"
#include "OutputSetT.h"
#include "ElementSupportT.h"
#include "ModelManagerT.h"
#include "iGridManagerT.h"
#include "iNodeT.h"
#include "CommManagerT.h"
#include "CommunicatorT.h"
#include "BasicFieldT.h"

#include "MeshFreeNodalShapeFunctionT.h"
#include "ContinuumMaterialT.h"
#include "SolidMaterialT.h"
#include "SolidMatSupportT.h"

/* materials lists */
#include "SSSolidMatList1DT.h"
#include "SSSolidMatList2DT.h"
#include "SSSolidMatList3DT.h"

#ifdef __QHULL__
#include "CompGeomT.h"
#endif

using namespace Tahoe;

/* constructors */
SCNIMFT::SCNIMFT(const ElementSupportT& support, const FieldT& field):
	ElementBaseT(support, field),
	fSD(ElementSupport().NumSD()),
	fMaterialList(NULL),
	fForce_man(0, fForce, field.NumDOF()),
	//fFakeGeometry(NULL),
	fVoronoi(NULL)
{
	SetName("mfparticle");

	/* set matrix format */
	fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
}

SCNIMFT::SCNIMFT(const ElementSupportT& support):
	ElementBaseT(support),
	//fFakeGeometry(NULL),
	fVoronoi(NULL)
{
	SetName("mfparticle");

	/* set matrix format */
	fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
}

/* destructor */
SCNIMFT::~SCNIMFT(void)
{
#ifdef __QHULL__		
	if (fVoronoi) 
	  delete fVoronoi;
#endif
	delete fMaterialList;
	delete fNodalShapes;
}

/* initialization */
void SCNIMFT::Initialize(void)
{
	const char caller[] = "SCNIMFT::Initialize";

	/* inherited */
	ElementBaseT::Initialize();

	/* re-dimension "element" force and stiffness contributions */
	fLHS.Dimension(fSD);
	
	/* allocate work space */
	fForce_man.SetMajorDimension(ElementSupport().NumNodes(), false);

	/* write parameters */
	ifstreamT& in = ElementSupport().Input();
	ostream& out = ElementSupport().Output();

	int qComputeVoronoiCell;
	in >> qComputeVoronoiCell;

	if (qComputeVoronoiCell)
	{
#ifndef __QHULL__
	        ExceptionT::GeneralFail(caller,"Requires the QHull library\n");
#else 

		// Do the heavy lifting for the Voronoi Diagram now
		fVoronoi = new CompGeomT(fDeloneVertices);
		fVoronoi->ComputeVoronoiDiagram();
		cout << " Computed\n";cout.flush();
		// Determine which cells are clipped by the boundary
		// Must be done before accessing data from the qhull library!!!
		fVoronoi->GenerateBoundaryCells(fBoundaryNodes, fBoundaryConnectivity,
			fBoundaryIsTriangulated);

		fVoronoiVertices.Alias(fVoronoi->VoronoiVertices());
		fVoronoiCells.Alias(fVoronoi->VoronoiCells()); 		
		fVoronoiFacetIndices.Alias(fVoronoi->VoronoiFacetIndices());

		// Data for integration over boundary of each Voronoi region
		fVoronoiFacetAreas.Alias(fVoronoi->VoronoiFacetAreas());
		fVoronoiFacetNormals.Alias(fVoronoi->VoronoiFacetNormals());
		fVoronoiCellVolumes.Alias(fVoronoi->VoronoiCellVolumes());

		fDeloneEdges.Alias(fVoronoi->DeloneEdges());
		fDualFacets.Alias(fVoronoi->DualFacets());
		fNumClippedFacets = fVoronoi->NumClippedFacets();
		fSelfDuals.Alias(fVoronoi->SelfDualFacets());
		fNumSelfDuals = fVoronoi->NumSelfDualFacets();

	
  		// Write output to file
		StringT vCellFile;
		in >> vCellFile;
		
		ofstreamT vout(vCellFile);

		if (vout.is_open())	
		{
			VoronoiDiagramToFile(vout);
			vout.close();
		}
		else 
  			cout  << " Unable to save data to file " << vCellFile << ". Ignoring error \n"; 
#endif
	} 
	else 
	{	// read in Voronoi information from a file
	  	StringT vCellFile;
		in >> vCellFile;
	    
		ifstreamT vin('#', vCellFile);

		if (!vin.is_open())
		  ExceptionT::GeneralFail(caller,"Unable to open file for reading");
	    
		VoronoiDiagramFromFile(vin);  
	    
		vin.close();
	}
	
	/* shape functions */
	/* only support single list of integration cells for now */
	if (fElementConnectivities.Length() > 1) {
	       ExceptionT::GeneralFail(caller,"Multiple ElementConnectivities not yet supported\n");
	}

	/* construct shape functions */
	fNodalShapes = new MeshFreeNodalShapeFunctionT(ElementSupport().NumSD(),
		ElementSupport().InitialCoordinates(), *fElementConnectivities[0], 
		fVoronoiVertices, ElementSupport().Input());
	if (!fNodalShapes) throw ExceptionT::kOutOfMemory;
	
	/* echo parameters */
	fNodalShapes->WriteParameters(ElementSupport().Output());
	
	/* initialize */
//	fNodalShapes->Initialize();

	/* MLS stuff */
	fNodalShapes->SetSupportSize();

	/* exchange nodal parameters (only Dmax for now) */
	const ArrayT<int>* p_nodes_in = ElementSupport().ExternalNodes();
	if (p_nodes_in)
	{
		/* skip MLS fit at external nodes */
		iArrayT nodes_in;
		nodes_in.Alias(*p_nodes_in);
		fNodalShapes->SetSkipNodes(nodes_in);
		
		/* exchange */
		CommManagerT& comm = ElementSupport().CommManager();

		/* send all */
		dArray2DT& nodal_params = fNodalShapes->NodalParameters();

		/* initialize the exchange */
		int id = comm.Init_AllGather(nodal_params);
		
		/* do the exchange */
		comm.AllGather(id, nodal_params);
		
		/* clear the communication */
		comm.Clear_AllGather(id);
	}
	
	/* set nodal neighborhoods */
	fNodalShapes->SetNeighborData();

	/* initialize support data */
/*	iArrayT surface_nodes;
	if (fAutoBorder) {
		ArrayT<StringT> IDs;
		ElementBlockIDs(IDs);
		ElementSupport().Model().SurfaceNodes(IDs, surface_nodes,
			fFakeGeometry);
	}
	MeshFreeFractureSupportT::InitSupport(
		ElementSupport().Input(), 
		ElementSupport().Output(),
		fElementCards, 
		surface_nodes, 
		NumDOF(), 
		ElementSupport().NumNodes(),
		&ElementSupport().Model()
	);*/
	
	//dArray2DT params(fNodes.Length(), 3);
	//params.SetColumn(0,0.);
	//params.SetColumn(1,0.);
	//params.SetColumn(2,0.);
	//fNodalShapes->SetNodalParameters(fNodes, params);
	
	/* final MLS initializations */
	fNodalShapes->WriteStatistics(ElementSupport().Output());
	
	/* initialize workspace for strain smoothing */
	ComputeBMatrices();	
	
	/** store shape functions at nodes */
	int nNodes = fNodes.Length();
	dArrayT nodalCoords;
	fNodalPhi.Dimension(nNodes);
	fNodalSupports.Dimension(nNodes);
	for (int i = 0; i < nNodes; i++)
	{
		nodalCoords.Set(fSD, fDeloneVertices(fNodes[i]));
	
		if (!fNodalShapes->SetFieldAt(nodalCoords, NULL)) // shift = 0 or not ?
			ExceptionT::GeneralFail("SCNIMFT::Initialize","Shape Function evaluation"
				"failed at Delone edge %d\n",fNodes[i]);		
		
		fNodalPhi[i].AppendArray(fNodalShapes->FieldAt().Length(),
									const_cast <double *> (fNodalShapes->FieldAt().Pointer()));
		fNodalSupports[i].AppendArray(fNodalShapes->Neighbors().Length(),
										const_cast <int *> (fNodalShapes->Neighbors().Pointer()));
	}
	
	/** Material Data */
	ReadMaterialData(in);
	
	//TEMP - only works for one material right now, else would have to check
	//       for the material active within the integration cell (element)
	/*if (HasActiveCracks() && fMaterialList->Length() != 1)
	{	
		cout << "\n MeshFreeSSSolidT::Initialize: can only have 1 material in the group\n";
		cout <<   "     with active cracks" << endl;
		throw ExceptionT::kBadInputValue;
	}*/

//TEMP - write nodal parameters
#if 0
	ostream& out = FEManager().Output();
	out << "\n MeshFreeFDSolidT::Initialize: d_max:\n";
	const dArray2DT& nodal_params = fMFShapes->NodalParameters();
	const iArrayT* node_map = FEManager().NodeMap();
	int d_width = out.precision() + kDoubleExtra;
	out << setw(kIntWidth) << "node"
	<< setw(nodal_params.MinorDim()*d_width) << "d_max" << '\n';
	if (node_map)
		for (int i = 0; i < nodal_params.MajorDim(); i++)
		{
			out << setw(kIntWidth) << (*node_map)[i] + 1;
			nodal_params.PrintRow(i, out);
			out << '\n';
		}
	else
		for (int j = 0; j < nodal_params.MajorDim(); j++)
		{
			out << setw(kIntWidth) << j + 1;
			nodal_params.PrintRow(j, out);
			out << '\n';
		}
#endif

//TEMP - U connectivities
#if 0
	out << "\n MeshFreeFDSolidT::Initialize: connects U: 0, 1,...\n";
	fElemNodesEX->WriteNumbered(out);
	out.flush();
#endif

//TEMP - trace node
#if 0
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		
	int hit_node = -1;
	if (size == 1)
		hit_node = 13607 - 1;
	else if (size == 2 && rank == 0)
		hit_node = 5420 - 1;
	if (hit_node > 0) TraceNode(ElementSupport().Output(), hit_node, *this);
#endif
}

/* form of tangent matrix */
GlobalT::SystemTypeT SCNIMFT::TangentType(void) const
{
	return GlobalT::kSymmetric;
}

/* NOT implemented. Returns an zero force vector */
void SCNIMFT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
}

/* writing output */
void SCNIMFT::RegisterOutput(void)
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

/* generate labels for output data */
void SCNIMFT::GenerateOutputLabels(ArrayT<StringT>& labels)
{
  	/* Reference Configuration */
	const char* ref[3] = {"X", "Y", "Z"};

	/* displacement labels */
	const char* disp[3] = {"D_X", "D_Y", "D_Z"};
	int num_labels = 2*fSD; // displacements
	int num_stress=0;

	const char* stress[6];
	const char* strain[6];
	
	stress[0] = "s11";
	stress[1] = "s22";
	strain[0] = "e11";
	strain[1] = "e22";
	num_stress = 3;
	
	if (fSD == 3)
	{
		num_stress=6;
		stress[2] = "s33";
	  	stress[3] = "s23";
	 	stress[4] = "s13";
	  	stress[5] = "s12";
	  	strain[2] = "e33";
	  	strain[3] = "e23";
	  	strain[4] = "e13";
	  	strain[5] = "e12";
	} 
	
	if (fSD == 2) 
	{
	  	stress[2] = "s12";
	  	strain[2] = "e12";
	}
		
	num_labels += 2 * num_stress + 1; 
	labels.Dimension(num_labels);
	int dex = 0;
	for (dex = 0; dex < fSD; dex++)
		labels[dex] = ref[dex];
	for (int ns = 0 ; ns < fSD; ns++)
	  	labels[dex++] = disp[ns];
	labels[dex++] = "mass";
	for (int ns = 0 ; ns < num_stress; ns++)
		labels[dex++] = strain[ns];
	for (int ns = 0 ; ns < num_stress; ns++)
		labels[dex++] = stress[ns];
}

void SCNIMFT::WriteOutput(void)
{
	//Nothin' to do -- moved to SS_SCNIMFT.cpp
}

/* compute specified output parameter(s) */
void SCNIMFT::SendOutput(int kincode)
{
#pragma unused(kincode)
	//TEMP: for now, do nothing
}

/* trigger reconfiguration */
GlobalT::RelaxCodeT SCNIMFT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	return relax;
}

/* write restart data to the output stream */
void SCNIMFT::WriteRestart(ostream& out) const
{
	ElementBaseT::WriteRestart(out);
	
}

/* read restart data to the output stream */
void SCNIMFT::ReadRestart(istream& in)
{
	ElementBaseT::ReadRestart(in);
}

/* contribution to the nodal residual forces */
//const dArray2DT& SCNIMFT::InternalForce(int group)
//{
	/* check */
//	if (group != Group())
//		ExceptionT::GeneralFail("SCNIMFT::InternalForce", 
//			"expecting solver group %d not %d", Group(), group);
//	return fForce;
//}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* return true if connectivities are changing */
bool SCNIMFT::ChangingGeometry(void) const
{
	return ElementSupport().CommManager().PartitionNodesChanging();
}

/* echo element connectivity data */
void SCNIMFT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
#pragma unused(out)
	const char caller[] = "SCNIMFT::EchoConnectivityData";
	
	/* access to the model database */
	ModelManagerT& model = ElementSupport().ModelManager();

	/* read node set ids */
	ArrayT<StringT> ids;
	model.NodeSetList(in, ids);
	model.ManyNodeSets(ids, fNodes);

	// temporarily, read in connectivies of the element group. I need this
	// to piggyback on the search that finds the support of each node
	iArrayT matnums;
	model.ElementBlockList(in, ids, matnums);
	if (ids.Length() != 1)
		ExceptionT::GeneralFail(caller,"Not helping my kludge by having more than one element block\n");
	
	fElementConnectivities.Dimension(1);
	model.ReadConnectivity(ids[0]);

	/* set pointer to connectivity list */
	fElementConnectivities[0] = model.ElementGroupPointer(ids[0]);
	
	/* get nodes and facets on the boundary */
	GeometryBaseT::CodeT facetType;
	iArrayT facet_numbers, element_numbers;
	model.SurfaceFacets(ids, facetType, fBoundaryConnectivity, fBoundaryNodes, 
		facet_numbers, element_numbers, NULL);
	fBoundaryIsTriangulated = (fSD == 2) ? (facetType == GeometryT::kLine) :
		(facetType == GeometryT::kTriangle);
	fBoundaryNodes.SortAscending();
	
	/* don't need this information */
	facet_numbers.Free();
	element_numbers.Free();
	    
	/* store block data  ASSUMING bth block is material number b*/
	//fBlockData[0].Set(new_id, 0, elem_count, 0); 
	  
	/* set up ElementCards and equation arrays */
	//fElementCards.Dimension(elem_count);
	//fEqnos.Dimension(1);
	//fEqnos[0].Dimension(elem_count,fSD);

	// Get nodal coordinates to use in Initialize
	fDeloneVertices.Dimension(fNodes.Length(), model.NumDimensions());
	fDeloneVertices.RowCollect(fNodes, model.Coordinates());
}

/* collecting element group equation numbers */
void SCNIMFT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_1)

	/* dimension equations array */
	fEqnos.Configure(nodeWorkSpace, NumDOF());

	/* get local equations numbers */
	Field().SetLocalEqnos(nodeWorkSpace, fEqnos);

	/* add to list of equation numbers */
	eq_2.Append(&fEqnos);
}

/* assemble particle mass matrix into LHS of global equation system */
void SCNIMFT::AssembleParticleMass(const double rho)
{
	
	fForce = 0.0;
	double* m = fForce.Pointer();
	int* nodes = fNodes.Pointer();
	double* volume = fVoronoiCellVolumes.Pointer();
	for (int i = 0; i < fNodes.Length(); i++)
	{
		for (int j = 0; j < fSD; j++)
			*m++ = *volume;
		volume++;
	}
	fForce *= rho;
	
	/* assemble all */
	ElementSupport().AssembleLHS(Group(), fForce, Field().Equations());
}

void SCNIMFT::ReadMaterialData(ifstreamT& in)
{
	const char caller[] = "SCNIMFT::ReadMaterialData";

	/* construct material list */
	int size;
	in >> size;
	fMaterialList = NewMaterialList(NumSD(), size);
	if (!fMaterialList) ExceptionT::OutOfMemory(caller);

	/* read */
	fMaterialList->ReadMaterialData(in);
	
	fMaterialNeeds.Dimension(fMaterialList->Length());
	for (int i = 0; i < fMaterialNeeds.Length(); i++)
	{
		/* allocate */
		ArrayT<bool>& needs = fMaterialNeeds[i];
		needs.Dimension(3);

		/* casts are safe since class contructs materials list */
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[i];
		SolidMaterialT* mat = (SolidMaterialT*) pcont_mat;

		/* collect needs */
		needs[0] = mat->NeedDisp();
		needs[1] = mat->NeedVel();
		needs[2] = mat->NeedLastDisp();
	}
	
	/* check range */
	/*for (int i = 0; i < fBlockData.Length(); i++)
		if (fBlockData[i].MaterialID() < 0 || fBlockData[i].MaterialID() >= size)
			ExceptionT::BadInputValue(caller, "material number %d for element block %d is out of range",
				fBlockData[i].MaterialID()+1, i+1);*/
}

int SCNIMFT::GlobalToLocalNumbering(iArrayT& nodes)
{
	if (!fNodes.Length())
		if (!nodes.Length())
			return 1;
		else 
			ExceptionT::GeneralFail("SCNIMFT::GlobalToLocalNumbering","No nodes exist\n");
	
	// Basic Idea: fNodes is sorted (it came from ModelManagerT::ManyNodeSets)
	// So, sort nodes with a key array and march down and compare. 
	iArrayT nodeMap(nodes.Length());
	nodeMap.SetValueToPosition();
	iArrayT nodeCopy(nodes.Length());
	nodeCopy = nodes;
	nodeMap.SortAscending(nodeCopy);

	// nodes[nodeMap[0]] is the smallest node in global numbering scheme
	// nodes[nodeMap[0]] should be that global node's position in fNodes
	int fNodesLen = fNodes.Length();
	int nodeMapLen = nodeMap.Length();
	int fNodesCtr = 0;
	int nodeMapCtr = 0;
	int *fNodesPtr = fNodes.Pointer();
	int *nodeMapPtr = nodeCopy.Pointer();
	
	if (*nodeMapPtr > *fNodesPtr)
		return 0;
	
	while (fNodesCtr < fNodesLen && nodeMapCtr < nodeMapLen) {
		while (*nodeMapPtr != *fNodesPtr && fNodesCtr < fNodesLen) {
			fNodesPtr++; 
			fNodesCtr++;
		}
		
		if (fNodesCtr != fNodesLen) {
			nodes[nodeMap[nodeMapCtr]] = fNodesCtr; // local numbering!
			nodeMapCtr++;
			nodeMapPtr++;
		}
	}
	
	if (nodeMapCtr != nodeMapLen)
		return 0;
	else
		return 1;
}

int SCNIMFT::GlobalToLocalNumbering(RaggedArray2DT<int>& nodes)
{
	iArrayT row_i;
	for (int i = 0; i < nodes.MajorDim(); i++) {
		row_i.Set(nodes.MinorDim(i), nodes(i)); 
		if (!GlobalToLocalNumbering(row_i))
			return 0;
	}
	
	return 1;
}

void SCNIMFT::InterpolatedFieldAtNodes(const iArrayT& nodes, dArray2DT& fieldAtNodes)
{
	/* displacements */
	const dArray2DT& u = Field()(0,0);
	dArrayT vec, values_i;
	for (int i = 0; i < nodes.Length(); i++) 
	{
		/* copy in */
		vec.Set(fSD, fieldAtNodes.Pointer() + i*fSD);
		vec = 0.;	
			
		LinkedListT<int>& supp_i = fNodalSupports[nodes[i]];
		LinkedListT<double>& phi_i = fNodalPhi[nodes[i]];
		supp_i.Top(); phi_i.Top();
		while (supp_i.Next() && phi_i.Next()) 
			vec.AddScaled(*(phi_i.CurrentValue()), u(*(supp_i.CurrentValue())));

	}

}

/** localNode is a local Number, so GlobalToLocalNumbering needs to have been called in whatever class 
  * calls this function. The node number returned in support are global. 
  */
void SCNIMFT::NodalSupportAndPhi(int localNode, LinkedListT<int>& support, LinkedListT<double>& phi)
{
#if __option(extended_errorcheck)
	if (localNode < 0 || localNode >= fNodes.Length()) throw ExceptionT::kSizeMismatch;
#endif

	support.Alias(fNodalSupports[localNode]);
	phi.Alias(fNodalPhi[localNode]);
}

void SCNIMFT::NodalSupportAndPhi(iArrayT& localNodes, RaggedArray2DT<int>& support, RaggedArray2DT<double>& phi)
{
	int nlnd = localNodes.Length();
	ArrayT<LinkedListT<int> > local_support(nlnd);
	ArrayT<LinkedListT<double> > local_phi(nlnd);
	
	for (int i = 0; i < nlnd; i++) {
		local_support[i].Alias(fNodalSupports[localNodes[i]]);
		local_phi[i].Alias(fNodalPhi[localNodes[i]]);
	}
	
	support.Configure(local_support);
	phi.Configure(local_phi);
	
	for (int i = 0; i < phi.MajorDim(); i++) {
		int* irow_i = support(i);
		double* drow_i = phi(i);
		LinkedListT<int>& ilist = local_support[i];
		LinkedListT<double>& dlist = local_phi[i];
		ilist.Top(); dlist.Top();
		while (ilist.Next() && dlist.Next()) {
			*irow_i++ = *(ilist.CurrentValue());
			*drow_i++ = *(dlist.CurrentValue());
		}
	}
}

int SCNIMFT::SupportSize(int localNode) 
{
	return fNodalPhi[localNode].Length();
}

void SCNIMFT::ComputeBMatrices(void)
{
	/* possible best implementation is to loop over all Delone edges
	 * and compute all the necessary values only once per Voronoi
	 * facet. This approach minimizes number of times that the support of
	 * an arbitrary point in space (the Voronoi facet centroid) has to be
	 * found.
	 */

	const RaggedArray2DT<int>& nodeSupport = fNodalShapes->NodeNeighbors();
	
	int nNodes = fNodes.Length();
	nodeWorkSpace.Dimension(nNodes);
	facetWorkSpace.Dimension(nNodes);
	
	dArrayT zeroFacet(3);
	zeroFacet = 0.0;
	for (int i = 0; i < nNodes; i++)
	{
		int l_supp_i = nodeSupport.MinorDim(i);
		iArrayT supp_i(l_supp_i);
		supp_i.Copy(nodeSupport(i));
		supp_i.SortAscending();
		nodeWorkSpace[i].AppendArray(l_supp_i, supp_i.Pointer());
		facetWorkSpace[i].AppendArray(l_supp_i, zeroFacet);
	}
	
	dArrayT facetCentroid(fSD), facetNormal(fSD), facetIntegral(fSD);
	double* currentB, *currentI;
	int n_0, n_1;
	bool traverseQ_0, traverseQ_1;
	int *next_0, *next_1;
	for (int i = 0; i < fDeloneEdges.MajorDim(); i++)
	{
		facetCentroid = 0.; 
		n_0 = fDeloneEdges(i,0);
		n_1 = fDeloneEdges(i,1);
		facetCentroid.Set(fSD, fDualFacetCentroids(i)); 
		facetNormal.DiffOf(fDeloneVertices(n_1), fDeloneVertices(n_0));
		facetNormal.UnitVector();
		
		if (!fNodalShapes->SetFieldAt(facetCentroid, NULL)) // shift = 0 or not ?
			ExceptionT::GeneralFail("SCNIMFT::ComputeBMatrices","Shape Function evaluation"
				"failed at Delone edge %d\n",i);
				
		const dArrayT& phiValues = fNodalShapes->FieldAt();			
		
		iArrayT centroid_cover(fNodalShapes->Neighbors());	
		int n_centroid_cover = centroid_cover.Length();	
		iArrayT centroid_cover_key(n_centroid_cover);
		centroid_cover_key.SetValueToPosition();
		centroid_cover_key.SortAscending(centroid_cover);

		LinkedListT<int>& supp_0 = nodeWorkSpace[n_0];
		LinkedListT<int>& supp_1 = nodeWorkSpace[n_1];
		LinkedListT< dArrayT >& bVectors_0 = facetWorkSpace[n_0];
		LinkedListT< dArrayT >& bVectors_1 = facetWorkSpace[n_1];
		int s_0 = -1;
		int s_1 = -1;
		/* Simultaneously loop over support of the two nodes that are endpoints of the
		 * current Delone edge and the nodes in the support of the midpoint of this
		 * edge. If a node covering the centroid is not in the support of n_0 or n_1,
		 * insert that covering node into the sorted list.
		 */
		 
		int* c = centroid_cover.Pointer();
		int* c_j = centroid_cover_key.Pointer();
		
		supp_0.Top(); bVectors_0.Top();
		supp_1.Top(); bVectors_1.Top();
		next_0 = supp_0.CurrentValue();
		next_1 = supp_1.CurrentValue();
		for (int j = 0; j < n_centroid_cover; j++, c++, c_j++)
		{
			facetIntegral = facetNormal;
			facetIntegral *= fDualAreas[i]*phiValues[*c_j];	
			
			if (next_0)
				traverseQ_0 = *next_0 <= *c;
			else
				traverseQ_0 = false;
					
			// advance supp_0 and supp_1 until they are greater than or equal to current node
			while (traverseQ_0 && supp_0.Next(s_0) && bVectors_0.Next())
			{
				next_0 = supp_0.PeekAhead(); 
				if (!next_0)
					traverseQ_0 = false;
				else
					if (*next_0 > *c)
						traverseQ_0 = false;
			}
				
			if (s_0 != *c) // means we're not at the end of the linked list
			{
				supp_0.InsertAtCurrent(*c);
				bVectors_0.InsertAtCurrent(zeroFacet);
				s_0 = *c;
				if (supp_0.AtTop()) // if we're inserting at the front, LinkedListT's behavior requires more work
				{
					supp_0.Next(); 
					bVectors_0.Next();
				}
			}
				
			currentI = facetIntegral.Pointer();
			currentB = bVectors_0.CurrentValue()->Pointer();
			for (int k = 0; k < fSD; k++)
				*currentB++ += *currentI++;
				
			if (next_1)
				traverseQ_1 = *next_1 <= *c;
			else
				traverseQ_1 = false;
				
			// advance supp_0 and supp_1 until they are greater than or equal to current node
			while (traverseQ_1 && supp_1.Next(s_1) && bVectors_1.Next())
			{
				next_1 = supp_1.PeekAhead(); 
				if (!next_1)
					traverseQ_1 = false;
				else
					if (*next_1 > *c)
						traverseQ_1 = false;
			}		
								
			if (s_1 != *c)
			{
				supp_1.InsertAtCurrent(*c);
				bVectors_1.InsertAtCurrent(zeroFacet);
				s_1 = *c;
				if (supp_1.AtTop()) // if we're inserting at the front, LinkedListT's behavior requires more work
				{
					supp_1.Next(); 
					bVectors_1.Next();
				}
			}
				 
			currentI = facetIntegral.Pointer();
			currentB =  bVectors_1.CurrentValue()->Pointer();
			for (int k = 0; k < fSD; k++)
				*currentB++ -= *currentI++; //NB change in sign; facet normal is inverted!
		}
	}
	
	/** Loop over remaining edges */
	for (int i = 0; i < fNonDeloneEdges.Length(); i++)
	{
		facetCentroid = 0.; 
		n_0 = fNonDeloneEdges[i];
		facetCentroid.Set(fSD, fNonDeloneCentroids(i)); 
		facetNormal.Set(fSD, fNonDeloneNormals(i));
		facetNormal.UnitVector();
		
		if (!fNodalShapes->SetFieldAt(facetCentroid, NULL)) // shift = 0 or not ?
			ExceptionT::GeneralFail("SCNIMFT::ComputeBMatrices","Shape Function evaluation"
				"failed at Delone edge %d\n",i);
				
		const dArrayT& phiValues = fNodalShapes->FieldAt();			
				
		iArrayT centroid_cover(fNodalShapes->Neighbors());	
		int n_centroid_cover = centroid_cover.Length();	
		iArrayT centroid_cover_key(n_centroid_cover);
		centroid_cover_key.SetValueToPosition();
		centroid_cover_key.SortAscending(centroid_cover);
		
		LinkedListT<int>& supp_0 = nodeWorkSpace[n_0];
		LinkedListT< dArrayT >& bVectors_0 = facetWorkSpace[n_0];
		int s_0;
		
		/* Merge support of the boundary node with covering of integration point
		 */
		int* c = centroid_cover.Pointer();
		int* c_j = centroid_cover_key.Pointer();
		
		supp_0.Top(); bVectors_0.Top();
		next_0 = supp_0.CurrentValue();
		for (int j = 0; j < n_centroid_cover; j++, c++, c_j++)
		{
			facetIntegral = facetNormal;
			facetIntegral *= fBoundaryIntegrationWeights[i]*phiValues[*c_j];		
		
			if (next_0)
				traverseQ_0 = *next_0 <= *c;
			else
				traverseQ_0 = false;
					
			// advance supp_0 and supp_1 until they are greater than or equal to current node
			while (traverseQ_0 && supp_0.Next(s_0) && bVectors_0.Next())
			{
				next_0 = supp_0.PeekAhead(); 
				if (!next_0)
					traverseQ_0 = false;
				else
					if (*next_0 > *c)
						traverseQ_0 = false;
			}
				
			if (s_0 != *c) // means we're not at the end of the linked list
			{
				supp_0.InsertAtCurrent(*c);
				bVectors_0.InsertAtCurrent(zeroFacet);
				s_0 = *c;
				if (supp_0.AtTop()) // if we're inserting at the front, LinkedListT's behavior requires more work
				{
					supp_0.Next(); 
					bVectors_0.Next();
				}
			}
				
			currentI = facetIntegral.Pointer();
			currentB =  bVectors_0.CurrentValue()->Pointer();
			for (int k = 0; k < fSD; k++)
				*currentB++ += *currentI++;
		}
	}
	
	// scale integrals by volumes of Voronoi cells
	dArrayT* currFacetIntegral;
	for (int i = 0; i < nNodes; i++)
	{
		LinkedListT<dArrayT>& bVectors_i = facetWorkSpace[i];
		LinkedListT<int>& nodes_i = nodeWorkSpace[i];
		bVectors_i.Top(); nodes_i.Top();
		while ((currFacetIntegral = bVectors_i.Next()))
			*currFacetIntegral *= 1./fVoronoiCellVolumes[i];
	}
	
	// move into more efficient storage for computation
	nodalCellSupports.Configure(nodeWorkSpace);
	bVectorArray.Configure(facetWorkSpace);
	
	for (int i = 0; i < nodalCellSupports.MajorDim(); i++) {
		int* irow_i = nodalCellSupports(i);
		dArrayT* drow_i = bVectorArray(i);
		LinkedListT<int>& ilist = nodeWorkSpace[i];
		LinkedListT<dArrayT>& dlist = facetWorkSpace[i];
		ilist.Top(); dlist.Top();
		while (ilist.Next() && dlist.Next()) {
			*irow_i++ = *(ilist.CurrentValue());
			*drow_i++ = *(dlist.CurrentValue());
		}
	}
}

void SCNIMFT::VoronoiDiagramToFile(ofstreamT& vout)
{
	
    int nVertices = fVoronoiVertices.MajorDim();
    int nSD = fVoronoiVertices.MinorDim();

    vout << fNodes.Length() << "\n";
    vout << nSD << "\n";
    vout << nVertices << "\n";

    // write out vertices of the VoronoiDiagram
    for (int i = 0; i < nVertices; i++)
    {
      	for (int j = 0; j < nSD; j++)
			vout << fVoronoiVertices(i,j) << " "; 
      	vout << "\n";
    }
    
    // write out Voronoi cells for each node
    for (int i = 0; i < fNodes.Length(); i++)
    {
		vout << i <<"\n";

		// number of vertices on its boundary
		vout << fVoronoiCells[i].Length() << " ";
		for (int j = 0; j < fVoronoiCells[i].Length(); j++)
		  vout << fVoronoiCells[i][j] << " ";
		vout << "\n";

		// number of facets
		vout << fVoronoiFacetIndices[i].Length() << "\n";
		
		// data for each facet
		for (int j = 0; j < fVoronoiFacetIndices[i].Length(); j++)
	  	{
	   		vout << j << " " << fVoronoiFacetIndices[i][j].Length() << " ";
	    	for (int k = 0; k < fVoronoiFacetIndices[i][j].Length(); k++)
	      		vout << fVoronoiFacetIndices[i][j][k] << " ";
	    
	    	// area of facet
	    	vout << fVoronoiFacetAreas[i][j] << " ";
	    
	    	// facet normal vector
	    	for (int k = 0; k < nSD; k++)
	      		vout << fVoronoiFacetNormals[i](j,k) << " ";
	    	vout << "\n";
	  	}
	
		// cell volume
		vout << fVoronoiCellVolumes[i] << "\n";
    }

    // write out Delone edge information
    int nDelone = fDeloneEdges.MajorDim();
    vout << nDelone << "\n"; // number of Delone edges
    for (int i = 0; i < nDelone; i++)
    {
    	vout << fDeloneEdges(i,0) << " " << fDeloneEdges(i,1) << "\n";
    }

    // write out areas and centroids of Voronoi facets dual to Delone edge
    // this assumes a 1-point integration scheme over the boundary of each node's cell
    if (fDualFacets.MajorDim() != nDelone)
    	ExceptionT::GeneralFail("SCNIMFT::VoronoiDiagramToFile","Dual edge/facet dimension mismatch\n");

    fDualAreas.Dimension(nDelone);
    fDualFacetCentroids.Dimension(nDelone, fSD);
    double *v1, *v2;
    for (int i = 0; i < fDualFacets.MajorDim(); i++) 
    {
    	v1 = fVoronoiVertices(fDualFacets(i,0));
    	v2 = fVoronoiVertices(fDualFacets(i,1));
    	fDualAreas[i] = sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1]));
    	vout << fDualAreas[i] << " ";
    	for (int j = 0; j < fSD; j++) {
			fDualFacetCentroids(i,j) = .5*(v1[j] + v2[j]);
			vout << fDualFacetCentroids(i,j) << " ";
    	}
    	vout << "\n";
    }
  
    // write out self-duals and allocate storage for them
    fNonDeloneEdges.Dimension(fNumSelfDuals);
    fNonDeloneNormals.Dimension(fNumSelfDuals, fSD);
    fNonDeloneCentroids.Dimension(fNumSelfDuals, fSD);
    fBoundaryIntegrationWeights.Dimension(fNumSelfDuals);
    vout << fNumSelfDuals << "\n";
    int ctr = 0;
    for (int i = 0; i < fSelfDuals.Length(); i++)
    	for (int j = 0; j < fSelfDuals[i].Length(); j++) 
    	{
			vout << fBoundaryNodes[i] << " ";
			fNonDeloneEdges[ctr++] = fBoundaryNodes[i];
    	}
    vout << "\n";

    // list of centroids of self-dual facets
    iArrayT* thisFacet;
    int thisFacetLength;
    dArrayT ptArray(fSD); // workspace for centroids
    double *pt = ptArray.Pointer(); 
    
    ctr = 0;
    for (int i = 0; i < fSelfDuals.Length(); i++)
    	for (int j = 0; j < fSelfDuals[i].Length(); j++) 
    	{
			thisFacet = &fVoronoiFacetIndices[fBoundaryNodes[i]][fSelfDuals[i][j]];
			thisFacetLength = thisFacet->Length();
			for (int l = 0; l < fSD; l++)
	 	 		pt[l] = 0.;
			for (int k = 0; k < thisFacetLength; k++) 
			{
	  			v1 = fVoronoiVertices(fVoronoiCells[fBoundaryNodes[i]][(*thisFacet)[k]]);
	  			for (int l = 0; l < fSD; l++)
	    			pt[l] += v1[l];
			}
			for (int l = 0; l < fSD; l++) 
			{
	  			pt[l] /= thisFacetLength;
	  			vout << pt[l] << " ";
	  			fNonDeloneCentroids(ctr,l) = pt[l];
			}
			ctr++;
			vout << "\n";
      	}

    // list of normals of self-dual facets
    ctr = 0;
    for (int i = 0; i < fSelfDuals.Length(); i++)
    	for (int j = 0; j < fSelfDuals[i].Length(); j++) 
    	{
			for (int k = 0; k < fSD; k++) 
			{
	  			fNonDeloneNormals(ctr,k) =  fVoronoiFacetNormals[fBoundaryNodes[i]](fSelfDuals[i][j],k);
	  			vout << fNonDeloneNormals(ctr,k) << " ";
			}
			vout << "\n";
			ctr++;
    	}

    // list of areas of self-dual facets
    ctr = 0;
    for (int i = 0; i < fSelfDuals.Length(); i++)
    	for (int j = 0; j < fSelfDuals[i].Length(); j++) 
    	{
			fBoundaryIntegrationWeights[ctr] =  fVoronoiFacetAreas[fBoundaryNodes[i]][fSelfDuals[i][j]];
			vout << fBoundaryIntegrationWeights[ctr++] << "\n";
      	}
}	
	
void SCNIMFT::VoronoiDiagramFromFile(ifstreamT& vin)
{
	const char caller[] = "SCNIMFT::VoronoiDiagramFromFile";	

    int nNodes, nSD, nVertices;
    vin >> nNodes;
    
    /* minor consistency checks */
    if (nNodes != fNodes.Length())
    {
		vin.close();
		ExceptionT::GeneralFail(caller,"Input Voronoi file does not match node number\n");
    }
    
    vin >> nSD;
    if (nSD != 2 )
    {
		vin.close();
		ExceptionT::GeneralFail(caller,"Input Voronoi file does not match SD\n");
    }
    
    vin >> nVertices;
    
    // allocate memory for Voronoi diagram data structures
    fVoronoiVertices.Dimension(nVertices, nSD);
    fVoronoiCells.Dimension(nNodes);
    fVoronoiFacetIndices.Dimension(nNodes);
    fVoronoiCellVolumes.Dimension(nNodes);
    fVoronoiFacetAreas.Dimension(nNodes);
    fVoronoiFacetNormals.Dimension(nNodes);

    for (int i = 0 ; i < nVertices; i++)
    	for (int j = 0; j < nSD; j++)
	  		vin >> fVoronoiVertices(i,j);

    for (int i = 0; i < nNodes; i++)
    {
		int itmp;
		vin >> itmp;
		if (itmp != i)
	  		ExceptionT::GeneralFail(caller,"Bad Input Voronoi file\n");

		// read in number of vertices on Voronoi cell i
		vin >> itmp;
		fVoronoiCells[i].Dimension(itmp);
		for (int j = 0; j < itmp; j++)
	  		vin >> fVoronoiCells[i][j];
	  		
		// number of facets
		int nFacets;
		vin >> nFacets;
		fVoronoiFacetIndices[i].Dimension(nFacets);
		fVoronoiFacetAreas[i].Dimension(nFacets);
		fVoronoiFacetNormals[i].Dimension(nFacets, nSD);
		for (int j = 0; j < nFacets; j++)
	  	{
	    	vin >> itmp;
	    	if (itmp != j)
	      		ExceptionT::GeneralFail(caller,"Bad Input Voronoi File\n");
	  
	    	// number of vertices in facet
	    	vin >> itmp;
	    	fVoronoiFacetIndices[i][j].Dimension(itmp);
	    
	    	for (int k = 0; k < itmp; k++)
				vin >> fVoronoiFacetIndices[i][j][k];
	    	
	    	vin >> fVoronoiFacetAreas[i][j];
	    	
	    	for (int k = 0; k < nSD; k++)
	      		vin >> fVoronoiFacetNormals[i](j,k);
        }
	
		vin >> fVoronoiCellVolumes[i];

	}
	
	/* Read in all DeloneEdges in or on the body */
	int nDelone;
	vin >> nDelone; 
	fDeloneEdges.Dimension(nDelone, 2);
	fDualAreas.Dimension(nDelone);
	fDualFacetCentroids.Dimension(nDelone, 2);
	
	for (int i = 0; i < nDelone; i++)
		vin >> fDeloneEdges(i,0) >> fDeloneEdges(i,1);
		
	for (int i = 0; i < nDelone; i++) {
		vin >> fDualAreas[i];
		for (int j = 0; j < fSD; j++)
			vin >> fDualFacetCentroids(i,j);
	}
		
	int nCentroids;
	vin >> nCentroids; // number of boundary facets (self-duals)
	fNonDeloneEdges.Dimension(nCentroids);
	fNonDeloneNormals.Dimension(nCentroids,fSD);
	fNonDeloneCentroids.Dimension(nCentroids,fSD);
	fBoundaryIntegrationWeights.Dimension(nCentroids);
	
	for (int i = 0; i < nCentroids; i++)
		vin >> fNonDeloneEdges[i];	

	for (int i = 0; i < nCentroids; i++)
		for (int j = 0; j < fSD; j++)
			vin >> fNonDeloneCentroids(i,j);
			
	for (int i = 0; i < nCentroids; i++)
		for (int j = 0; j < fSD; j++)
			vin >> fNonDeloneNormals(i,j);
			
	for (int i = 0; i < nCentroids; i++)
		vin >> fBoundaryIntegrationWeights[i];
}

// XML stuff below

/* describe the parameters needed by the interface */
void SCNIMFT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ElementBaseT::DefineParameters(list);

}

/* information about subordinate parameter lists */
void SCNIMFT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ElementBaseT::DefineSubs(sub_list);

}

/* return the description of the given inline subordinate parameter list */
void SCNIMFT::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
	SubListT& sub_sub_list) const
{
	/* inherited */
	ElementBaseT::DefineInlineSub(sub, order, sub_sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SCNIMFT::NewSub(const StringT& list_name) const
{
	/* inherited */
	return ElementBaseT::NewSub(list_name);
}
