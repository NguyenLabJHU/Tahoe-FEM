/* $Id: SCNIMFT.cpp,v 1.46 2005-01-19 17:46:18 cjkimme Exp $ */
#include "SCNIMFT.h"

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
#include "ParentDomainT.h"
#include "ParameterUtils.h"
#include "ParameterContainerT.h"
#include "InverseMapT.h"

#include "MeshFreeSupport2DT.h"
#include "MeshFreeSupport3DT.h"
#include "MeshFreeNodalShapeFunctionT.h"
#include "ContinuumMaterialT.h"
#include "SolidMaterialT.h"
#include "SolidMatSupportT.h"

#include "Traction_CardT.h"

/* materials lists */
#include "SSSolidMatList1DT.h"
#include "SSSolidMatList2DT.h"
#include "SSSolidMatList3DT.h"

#ifdef __QHULL__
#include "CompGeomT.h"
#endif

using namespace Tahoe;

/* uncomment this line for many, many boundary integration points;
 * --set the number of points in the routine compute B matrices
 */
//#define SEPARATE_BOUNDARY_INTEGRATION
/* uncomment this line for one point boundary integration rules over
 * voronoi facets. SEPARATE_BOUNDARY_INTEGRATION must be defined, too,
 * but its effects are overidden when EVALUATE_AT_NODES is defined
 */
//#define EVALUATE_AT_NODES
/* define method one for hard-coded boundary facet centroids for natural
 * BCs. 
 * define method two to use the same boundary integration as in 
 * ComputeBMatrices 
 */
//#define FIRST_WAY
#ifndef FIRST_WAY
#define SECOND_WAY
#endif
const int kNoTractionVector = -1;

/* constructors */
SCNIMFT::SCNIMFT(const ElementSupportT& support, const FieldT& field):
	ElementBaseT(support),
	fSD(ElementSupport().NumSD()),
	fMaterialList(NULL),
	fVoronoi(NULL),
	fNodalShapes(NULL),
	qComputeVoronoiCell(false),
	qJustVoronoiDiagram(false),
	fNumIP(1),
	vCellFile("voronoidiagram")
{
#pragma unused(field)

	SetName("mfparticle");

	/* set matrix format */
	fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
}

SCNIMFT::SCNIMFT(const ElementSupportT& support):
	ElementBaseT(support),
	fSD(0),
	fMaterialList(NULL),
	fVoronoi(NULL),
	fNodalShapes(NULL),
	qComputeVoronoiCell(false),
	qJustVoronoiDiagram(false),
	fNumIP(1),
	vCellFile("voronoidiagram")
{
	SetName("mfparticle");

	/* set matrix format */
	fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
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
void SCNIMFT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SCNIMFT::TakeParameterList";
	
	/* What dimension? */
	fSD = ElementSupport().NumSD();

	/* resolve path to Voronoi file */
	StringT path;
	path.FilePath(ElementSupport().InputFile());
	vCellFile = list.GetParameter("voronoi_file");
	vCellFile.ToNativePathName();
	vCellFile.Prepend(path);

	/* number of integration points used for surface integrals */
	fNumIP = list.GetParameter("num_ip");

	qComputeVoronoiCell = list.GetParameter("compute_voronoi");
	qJustVoronoiDiagram = list.GetParameter("just_voronoi_diagram");
	
	/* get parameters needed to construct shape functions */
	fMeshfreeParameters = list.ListChoice(*this, "meshfree_support_choice");
	
	/* access to the model database */
	ModelManagerT& model = ElementSupport().ModelManager();

	/* extract particle ID's */
	const ParameterListT& particle_ID_params = list.GetList("mf_particle_ID_list");
	ArrayT<StringT> particle_ID_list;
	StringListT::Extract(particle_ID_params, particle_ID_list);

	//get nodes from ModelManagerT
	model.ManyNodeSets(particle_ID_list, fNodes);

	/* inherited */
	ElementBaseT::TakeParameterList(list);

	/* re-dimension "element" force and stiffness contributions */
	fLHS.Dimension(fSD);
	
	/* allocate work space */
	fForce_man.SetWard(0, fForce, fSD);
	fForce_man.SetMajorDimension(ElementSupport().NumNodes(), false);

	/* write parameters */
	ostream& out = ElementSupport().Output();

	if (qJustVoronoiDiagram)
	  qComputeVoronoiCell = true;

	// Do the heavy lifting for the Voronoi Diagram now
	if (qComputeVoronoiCell) {
#ifndef __QHULL__
	        ExceptionT::GeneralFail(caller,"Requires the QHull library\n");
#else 

		fVoronoi = new CompGeomT(fDeloneVertices);
		fVoronoi->ComputeVoronoiDiagram(); 
		
		// Determine which cells are clipped by the boundary
		// Must be done before accessing data from the qhull library!!!
		fVoronoi->GenerateBoundaryCells(fBoundaryNodes, fBoundaryConnectivity,
			fBoundaryIsTriangulated);

		InitializeVoronoiData();
	
  		// Write output to file
		ofstreamT vout;
	        vout.open(vCellFile);

		if (vout.is_open())	 {
			VoronoiDiagramToFile(vout);
			vout.close();
			if (qJustVoronoiDiagram)
			  ExceptionT::GeneralFail(caller,"Thank you. Computation Successful.\n");
		} else {
  			cout  << " Unable to save data to file " << vCellFile << ". Ignoring error \n"; 
			if (qJustVoronoiDiagram)
			  ExceptionT::GeneralFail(caller,"Sorry. Unable to write to file.\n");
		}

	
#endif
	} 
	else  {	// read in Voronoi information from a file
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
	fNodalShapes = new MeshFreeNodalShapeFunctionT(fSD,
		ElementSupport().InitialCoordinates(), *fElementConnectivities[0], 
		fVoronoiVertices, *fMeshfreeParameters);
	if (!fNodalShapes) throw ExceptionT::kOutOfMemory;
	
	/* echo parameters */
	fNodalShapes->WriteParameters(ElementSupport().Output());

	/* MLS stuff */
	fNodalShapes->SetSupportSize();

	/* construct meshfree support before calling inherited method because
	 * support class needed to construct shape functions */
	//	fMFSupport = new MeshFreeSupport2DT;
	//      fMFSupport->TakeParameterList(list.GetList("meshfree_support_2D"));

	/* exchange nodal parameters (only Dmax for now) */
	const ArrayT<int>* p_nodes_in = ElementSupport().ExternalNodes();
	if (p_nodes_in) {

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
	
	/* final MLS initializations */
	fNodalShapes->WriteStatistics(ElementSupport().Output());
	
	/* initialize workspace for strain smoothing */
	ComputeBMatrices();	
	
	/** store shape functions at nodes */
	int nNodes = fNodes.Length();
	dArrayT nodalCoords;
	ArrayT< LinkedListT<double> > nodal_phi;
	ArrayT< LinkedListT<int> > nodal_supports;
	nodal_phi.Dimension(nNodes);
	nodal_supports.Dimension(nNodes);
	for (int i = 0; i < nNodes; i++) {
		nodalCoords.Set(fSD, fDeloneVertices(i));
	
		if (!fNodalShapes->SetFieldAt(nodalCoords, NULL)) // shift = 0 or not ?
			ExceptionT::GeneralFail("SCNIMFT::TakeParameterList","Shape Function evaluation"
				"failed at Delone edge %d\n",fNodes[i]);		
		
		nodal_phi[i].AppendArray(fNodalShapes->FieldAt().Length(),
									const_cast <double *> (fNodalShapes->FieldAt().Pointer()));
		nodal_supports[i].AppendArray(fNodalShapes->Neighbors().Length(),
										const_cast <int *> (fNodalShapes->Neighbors().Pointer()));
	}
	
	// move into RaggedArray2DT's
	// move into more efficient storage for computation
	fNodalPhi.Configure(nodal_phi);
	fNodalSupports.Configure(nodal_supports);
	
	if (nodal_supports.Length() != nodal_phi.Length())
		ExceptionT::GeneralFail(caller,"nodal support indices and shape function values do not match\n");
		
	for (int i = 0; i < nodal_supports.Length(); i++) {
		int* irow_i = fNodalSupports(i);
		double* drow_i = fNodalPhi(i);
		LinkedListT<int>& ilist = nodal_supports[i];
		LinkedListT<double>& dlist = nodal_phi[i];
		ilist.Top(); dlist.Top();
		while (ilist.Next() && dlist.Next()) {
			*irow_i++ = *(ilist.CurrentValue());
			*drow_i++ = *(dlist.CurrentValue());
		}
	}
#ifdef FIRST_WAY	
	/** do the same for the boundary facets */
	nNodes = fNonDeloneEdges.Length();
	dArrayT boundaryIPCoords(fSD);
	ArrayT< LinkedListT<double> > boundary_phi;
	ArrayT< LinkedListT<int> > boundary_supports;
	boundary_phi.Dimension(nNodes);
	boundary_supports.Dimension(nNodes);
	for (int i = 0; i < nNodes; i++) {
		double* v1 = fVoronoiVertices(fSelfDualFacets(i,0));
		double* v2 = fVoronoiVertices(fSelfDualFacets(i,1));
		for (int j = 0; j < fSD; j++)
			boundaryIPCoords[j] = v1[j] + v2[j];
		boundaryIPCoords /= double(fSD);
		
		if (!fNodalShapes->SetFieldAt(boundaryIPCoords, NULL)) // shift = 0 or not ?
			ExceptionT::GeneralFail("SCNIMFT::TakeParameterList","Shape Function evaluation"
				"failed at non-Delone edge %d\n",fNonDeloneEdges[i]);		
		
		boundary_phi[i].AppendArray(fNodalShapes->FieldAt().Length(),
									const_cast <double *> (fNodalShapes->FieldAt().Pointer()));
		boundary_supports[i].AppendArray(fNodalShapes->Neighbors().Length(),
										const_cast <int *> (fNodalShapes->Neighbors().Pointer()));
	}
	
	// move into RaggedArray2DT's
	// move into more efficient storage for computation
	fBoundaryPhi.Configure(boundary_phi);
	fBoundarySupports.Configure(boundary_supports);
	
	if (boundary_supports.Length() != boundary_phi.Length())
		ExceptionT::GeneralFail(caller,"nodal support indices and shape function values do not match\n");
		
	for (int i = 0; i < boundary_supports.Length(); i++) {
		double jw = fBoundaryIntegrationWeights[i];
		int* irow_i = fBoundarySupports(i);
		double* drow_i = fBoundaryPhi(i);
		LinkedListT<int>& ilist = boundary_supports[i];
		LinkedListT<double>& dlist = boundary_phi[i];
		ilist.Top(); dlist.Top();
		while (ilist.Next() && dlist.Next()) {
			*irow_i++ = *(ilist.CurrentValue());
			*drow_i++ = *(dlist.CurrentValue()) * jw;
		}
	}
#endif	
	
	/** Material Data */
	ParameterListT mat_params;
	CollectMaterialInfo(list, mat_params);
	fMaterialList = NewMaterialList(mat_params.Name(), mat_params.NumLists());
 
 	if (!fMaterialList)
	  ExceptionT::GeneralFail(caller,"could not construct material list \"%s\"", mat_params.Name().Pointer());
	fMaterialList->TakeParameterList(mat_params);
	
	/* body force */
	const ParameterListT* body_force = list.List("body_force");
	if (body_force) {
		int schedule = body_force->GetParameter("schedule");
		fBodySchedule = ElementSupport().Schedule(--schedule);

		/* body force vector */
		const ArrayT<ParameterListT>& body_force_vector = body_force->Lists();
		if (body_force_vector.Length() != NumDOF())
			ExceptionT::BadInputValue(caller, "body force is length %d not %d",
				body_force_vector.Length(), NumDOF());
		fBody.Dimension(NumDOF());
		for (int i = 0; i < fBody.Length(); i++)
			fBody[i] = body_force_vector[i].GetParameter("value");
	}
	
	/* extract natural boundary conditions */
	TakeNaturalBC(list);
}

/* extract natural boundary condition information */
void SCNIMFT::TakeNaturalBC(const ParameterListT& list)
{
	const char caller[] = "SCNIMFT::TakeNaturalBC";

	int num_natural_bc = list.NumLists("natural_bc");
	
	// allocate data structures
	fTractionVectors.Dimension(num_natural_bc, fSD == 2 ? 3 : 6);
	fTractionVectors = 0.;
	
	// 
	fTractionBoundaryCondition.Dimension(fBoundaryIntegrationWeights.Length());
	fTractionBoundaryCondition = kNoTractionVector;
	
	if (num_natural_bc > 0)
	{
		/* TEMP - turn on traction boundary condition for all boundary nodes */
		fTractionBoundaryCondition = 0;
		
		/* model manager */
		ModelManagerT& model = ElementSupport().ModelManager();
	
		/* temp space */
		ArrayT<StringT> block_ID(num_natural_bc);
	    ArrayT<iArray2DT> localsides(num_natural_bc);
	    iArrayT LTf(num_natural_bc);
	    ArrayT<Traction_CardT::CoordSystemT> coord_sys(num_natural_bc);
	    ArrayT<dArray2DT> values(num_natural_bc);

	    /* nodes on element facets */
	    iArrayT num_facet_nodes;
	    num_facet_nodes = fSD;
	    
	    /* loop over natural BC's */
	    int tot_num_sides = 0;
	    for (int i = 0; i < num_natural_bc; i++) 
	   	{
	    	const ParameterListT& natural_bc = list.GetList("natural_bc", i);
	    
	    	/* side set */
	    	const StringT& ss_ID = natural_bc.GetParameter("side_set_ID");
			localsides[i] = model.SideSet(ss_ID);
			int num_sides = localsides[i].MajorDim();
			tot_num_sides += num_sides;
			if (num_sides > 0)
			{
				block_ID[i] = model.SideSetGroupID(ss_ID);
				LTf[i] = natural_bc.GetParameter("schedule");
				coord_sys[i] = Traction_CardT::int2CoordSystemT(natural_bc.GetParameter("coordinate_system"));

				/* switch to elements numbering within the group */
				iArray2DT& side_set = localsides[i];
				iArrayT elems(num_sides);
				//side_set.ColumnCopy(0, elems);
				//BlockToGroupElementNumbers(elems, block_ID[i]);
				//side_set.SetColumn(0, elems);

				/* all facets in set must have the same number of nodes */
				//int num_nodes = num_facet_nodes[side_set(0,1)];
				//for (int f = 0; f < num_sides; f++)
				//	if (num_facet_nodes[side_set(f,1)] != num_nodes)
				//		ExceptionT::BadInputValue(caller, "faces side set \"%s\" have different numbers of nodes",
				//			ss_ID.Pointer());

				/* read traction nodal values */
				//dArray2DT& nodal_values = values[i];
				//nodal_values.Dimension(num_nodes, NumDOF());
				int num_traction_vectors = natural_bc.NumLists("DoubleList");
				//if (num_traction_vectors != 1 && num_traction_vectors != num_nodes)
				//	ExceptionT::GeneralFail(caller, "expecting 1 or %d vectors not %d",
				//		num_nodes, num_traction_vectors);
						
				/* constant over the face */
				if (num_traction_vectors == 1) {
					const ParameterListT& traction_vector = natural_bc.GetList("DoubleList");
					int dim = traction_vector.NumLists("Double");
						
					int minor_dim = fSD == 2 ? 3 : 6;
					if (dim != minor_dim)
						ExceptionT::GeneralFail(caller, "expecting traction vector length %d not %d",
							NumDOF(), dim);
					dArrayT t(minor_dim);
					for (int j = 0; j < minor_dim; j++)
						t[j] = traction_vector.GetList("Double", j).GetParameter("value");	
							
					fTractionVectors.SetRow(0, t);

					/* same for all face nodes */
					//for (int f = 0; f < NumDOF(); f++) {
					//	double t = traction_vector.GetList("Double", f).GetParameter("value");
					//	nodal_values.SetColumn(f, t);
					//}
				}
				else
				{
					ExceptionT::GeneralFail(caller,"Only constant traction over sideset implemented\n");
				
					/* read separate vector for each face node */
					//dArrayT t;
					//for (int f = 0; f < num_nodes; f++) {
						//const ParameterListT& traction_vector = natural_bc.GetList("DoubleList", f);
						//int dim = traction_vector.NumLists("Double");
						//if (dim != NumDOF())
							//ExceptionT::GeneralFail(caller, "expecting traction vector length %d not %d",
								//NumDOF(), dim);

						//nodal_values.RowAlias(f, t);
						//for (int j = 0; j < NumDOF(); j++)
							//t[j] = traction_vector.GetList("Double", j).GetParameter("value");
					}
				}
			}
	    }

		/* allocate all traction BC cards */
	    //fTractionList.Dimension(tot_num_sides);

	    /* correct numbering offset */
	    //LTf--;

		/* define traction cards */
		//if (tot_num_sides > 0)
		//{
			//iArrayT loc_node_nums;
			//int dex = 0;
			//for (int i = 0; i < num_natural_bc; i++)
			//{
				/* set traction BC cards */
				//iArray2DT& side_set = localsides[i];
				//int num_sides = side_set.MajorDim();
				//for (int j = 0; j < num_sides; j++)
				//{					
					/* get facet local node numbers */
					//We don't have a domain geometry. Can't call this.
					//fNodalShapes->NodesOnFacet(side_set(j, 1), loc_node_nums);
					//if (fSD != 2)
					//	ExceptionT::GeneralFail(caller,"Specialized to 2D so far\n");
					//loc_node_nums.Dimension(2);
					//loc_node_nums[0] = side_set(j,1);
					//if (side_set(j,1) != 3)
					//	loc_node_nums[1] = side_set(j,1) + 1;
					//else
					//	loc_node_nums[1] = 0;
					
					/* set and echo */
					//fTractionList[dex++].SetValues(ElementSupport(), side_set(j,0), side_set (j,1), LTf[i],
					//	 coord_sys[i], loc_node_nums, values[i]);
				//}
			//}
		//}

		/* check coordinate system specifications */
		//if (NumSD() != NumDOF())
		//	for (int i = 0; i < fTractionList.Length(); i++)
		//		if (fTractionList[i].CoordSystem() != Traction_CardT::kCartesian)
		//			ExceptionT::BadInputValue(caller, "coordinate system must be Cartesian if (nsd != ndof) for card %d", i+1);
	//}
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
	if (parition_nodes) {
		int num_nodes = parition_nodes->Length();
		fPointConnectivities.Alias(num_nodes, 1, parition_nodes->Pointer());
	}
	else { /* ALL nodes */
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
	
	if (fSD == 3) {
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
	
	if (fSD == 2) {
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

/* construct field */
void SCNIMFT::NodalDOFs(const iArrayT& nodes, dArray2DT& DOFs) const
{
	/* get local numbers */
	iArrayT nodes_local = nodes;
	if (!GlobalToLocalNumbering(nodes_local))
		ExceptionT::GeneralFail("SCNIMFT::NodalDOFs", "map to local numbering failed");

	/* compute field at nodes */
	InterpolatedFieldAtNodes(nodes_local, DOFs);
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
void SCNIMFT::DefineElements(const ArrayT<StringT>& block_ID, const ArrayT<int>& mat_index) 
{
#pragma unused(mat_index)

	const char caller[] = "SCNIMFT::DefineElements";
	
	/* access to the model database */
	ModelManagerT& model = ElementSupport().ModelManager();

	fElementConnectivities.Dimension(1);

	// NB THIS IS SPECIALIZED TO ONLY ONE ELEMENT BLOCK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	model.ReadConnectivity(block_ID[0]);

	/* set pointer to connectivity list */
	fElementConnectivities[0] = model.ElementGroupPointer(block_ID[0]);
	
	/* get nodes and facets on the boundary */
	GeometryBaseT::CodeT facetType;
	iArrayT facet_numbers, element_numbers;
	model.SurfaceFacets(block_ID, facetType, fBoundaryConnectivity, fBoundaryNodes, 
		facet_numbers, element_numbers, NULL);
	fBoundaryIsTriangulated = (fSD == 2) ? (facetType == GeometryT::kLine) :
		(facetType == GeometryT::kTriangle);
	fBoundaryNodes.SortAscending();
	
	/* write boundary nodes to output */
	if (ElementSupport().PrintInput()) {
		ofstreamT& out = ElementSupport().Output();
		fBoundaryNodes++;
		out << "\n " << caller << ": boundary nodes\n" << fBoundaryNodes.wrap(10) << endl;
		fBoundaryNodes--;
	}

	/* convert to local numbering for qhull */
	GlobalToLocalNumbering(fBoundaryNodes);
	iArrayT trickArray(fBoundaryConnectivity.Length(), fBoundaryConnectivity.Pointer());
	GlobalToLocalNumbering(trickArray);
	
	/* don't need this information */
	facet_numbers.Free();
	element_numbers.Free();

	// Get nodal coordinates to use in Initialize
	fDeloneVertices.Dimension(fNodes.Length(), fSD);
	fDeloneVertices.RowCollect(fNodes, model.Coordinates());

	/* set up element cards for state variable storage */
	fElementCards.Dimension(fNodes.Length()); /* one card per node */
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
	int* nodes = fNodes.Pointer();
	double* volume = fVoronoiCellVolumes.Pointer();
	for (int i = 0; i < fNodes.Length(); i++) {
		double* m = fForce(*nodes++);
		for (int j = 0; j < fSD; j++)
			*m++ = *volume;
		volume++;
	}
	fForce *= rho;
	
	/* assemble all */
	ElementSupport().AssembleLHS(Group(), fForce, Field().Equations());
}

/* contribution from natural BCs */
void SCNIMFT::RHSDriver(void) 
{
	// need a node to boundary facet map
	//fBoundaryIntegrationWeights and fNonDeloneEdges have weights and node numbers
	//I could also keep an fTractionBoundaryCondition which would have the number of
	//the traction boundary condition or -1 for none
	//that's the kludgiest way to get it going for now
	if (fTractionVectors.MajorDim()) {
	
		fForce = 0.;
	
		for (int i = 0; i < fNonDeloneEdges.Length(); i++) 
			if (fTractionBoundaryCondition[i] != kNoTractionVector) {
				dArrayT workspace(fSD == 2 ? 3 : 6); 
				dArrayT traction_vector(fSD);
				fTractionVectors.RowCopy(fTractionBoundaryCondition[i], workspace.Pointer());
				traction_vector[0] = workspace[0]*fNonDeloneNormals(i,0) + .5*workspace[2]*fNonDeloneNormals(i,1);
				traction_vector[1] = workspace[1]*fNonDeloneNormals(i,1) + .5*workspace[2]*fNonDeloneNormals(i,0);
 			
				int* supp_i = fBoundarySupports(i) ;
				double* phi_i = fBoundaryPhi(i);
				int n_supp = fBoundaryPhi.MinorDim(i);
				for (int j = 0; j < n_supp; j++) {
					double* fint = fForce(*supp_i++);
					for (int k = 0; k < fSD; k++) 
						*fint++ += traction_vector[k]*(*phi_i); 
					phi_i++;
				}
			}

		// fForce gets multiplied by constKd?

		ElementSupport().AssembleRHS(Group(),fForce,Field().Equations());	
	}
}

int SCNIMFT::GlobalToLocalNumbering(iArrayT& nodes) const
{
	InverseMapT inv_map;
	inv_map.SetMap(fNodes);
	for (int i = 0; i < nodes.Length(); i++)
	  nodes[i] = inv_map.Map(nodes[i]);
	
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

void SCNIMFT::InterpolatedFieldAtNodes(const iArrayT& nodes, dArray2DT& fieldAtNodes) const
{
	/* displacements */
	const dArray2DT& u = Field()(0,0);
	dArrayT vec, values_i;
	for (int i = 0; i < nodes.Length(); i++) {
		/* copy in */
		vec.Set(fSD, fieldAtNodes.Pointer() + i*fSD);
		vec = 0.;	
		
		int node_i = nodes[i];
		const int* nodal_supp = fNodalSupports(node_i);
		const double* phi_i = fNodalPhi(node_i);
		for (int j = 0; j < fNodalPhi.MinorDim(node_i); j++)
			vec.AddScaled(*phi_i++, u(*nodal_supp++));
	}
}

/** localNodes are local Numbers, so GlobalToLocalNumbering needs to have been called in whatever class 
  * calls this function. The node numbers returned in support are global. 
  */
void SCNIMFT::NodalSupportAndPhi(const iArrayT& localNodes, RaggedArray2DT<int>& support, 
	RaggedArray2DT<double>& phi) const
{
	int nlnd = localNodes.Length();
	iArrayT minorDims(localNodes.Length());
	for (int i = 0; i < nlnd; i++) 
		minorDims[i] = fNodalPhi.MinorDim(localNodes[i]);
		
	support.Configure(minorDims);
	phi.Configure(minorDims);
	
	const int *lndi = localNodes.Pointer();
	for (int i = 0; i < nlnd; i++) {
		support.SetRow(i, fNodalSupports(*lndi));
		phi.SetRow(i, fNodalPhi(*lndi++));
	}
}

int SCNIMFT::SupportSize(int localNode) const {
	return fNodalPhi.MinorDim(localNode);
}

void SCNIMFT::InitializeVoronoiData(void)
{
#ifdef __QHULL__
  // use qhull's data structures but make our own specialized versions
  CompGeomT::ConvexHullMap selfDuals;
  CompGeomT::VoronoiDiagramMap voronoiFacetIndices;
  ArrayT<dArrayT> voronoiFacetAreas;
  ArrayT<dArray2DT> voronoiFacetNormals;
  CompGeomT::ConvexHullMap voronoiCells;

  fVoronoiVertices.Alias(fVoronoi->VoronoiVertices());
  voronoiCells.Alias(fVoronoi->VoronoiCells()); 		
  voronoiFacetIndices.Alias(fVoronoi->VoronoiFacetIndices());
  
  // Data for integration over boundary of each Voronoi region
  voronoiFacetAreas.Alias(fVoronoi->VoronoiFacetAreas());
  voronoiFacetNormals.Alias(fVoronoi->VoronoiFacetNormals());
  fVoronoiCellVolumes.Alias(fVoronoi->VoronoiCellVolumes());
  
  fDeloneEdges.Alias(fVoronoi->DeloneEdges());
  fDualFacets.Alias(fVoronoi->DualFacets());
  selfDuals.Alias(fVoronoi->SelfDualFacets());
  int numSelfDuals = fVoronoi->NumSelfDualFacets();

  fNonDeloneEdges.Dimension(numSelfDuals);
  fNonDeloneNormals.Dimension(numSelfDuals, fSD);
  fBoundaryIntegrationWeights.Dimension(numSelfDuals);
  fSelfDualFacets.Dimension(numSelfDuals, 2);
  
  int ctr = 0;
  double *v1;
  
  // list of centroids of self-dual facets
  iArrayT* thisFacet;
  int thisFacetLength;
  dArrayT ptArray(fSD); // workspace for centroids
  double *pt = ptArray.Pointer(); 
  
  for (int i = 0; i < selfDuals.Length(); i++)
    	for (int j = 0; j < selfDuals[i].Length(); j++) {
      		fNonDeloneEdges[ctr] = fBoundaryNodes[i];
      		thisFacet = &voronoiFacetIndices[fBoundaryNodes[i]][selfDuals[i][j]];
      		thisFacetLength = thisFacet->Length();
      		ptArray = 0.;
      		for (int k = 0; k < thisFacetLength; k++) {
				v1 = fVoronoiVertices(voronoiCells[fBoundaryNodes[i]][(*thisFacet)[k]]);
				for (int l = 0; l < fSD; l++)
	  				pt[l] += v1[l];
				fSelfDualFacets(ctr,k) =  voronoiCells[fBoundaryNodes[i]][(*thisFacet)[k]]; 
      		}
      
      		for (int k = 0; k < fSD; k++) 
				fNonDeloneNormals(ctr,k) =  voronoiFacetNormals[fBoundaryNodes[i]](selfDuals[i][j],k);
      
      		fBoundaryIntegrationWeights[ctr] =  voronoiFacetAreas[fBoundaryNodes[i]][selfDuals[i][j]];
      		ctr++;
    	}
    	
  // Data for Axisymmetric mass matrix calculation
  int nVoronoiCells = voronoiCells.Length();
  fVoronoiCellCentroids.Dimension(nVoronoiCells, fSD);
  fVoronoiCellCentroids = 0.;
  for (int i = 0; i < nVoronoiCells; i++) {
  	iArrayT& cell_i = voronoiCells[i];
  	for (int j = 0; j < cell_i.Length(); j++) 
	  fVoronoiCellCentroids.AddToRowScaled(i,1.,fVoronoiVertices(cell_i[j]));
	fVoronoiCellCentroids.ScaleRow(i,1./double(cell_i.Length()));
  }
#endif
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
	for (int i = 0; i < nNodes; i++) {
		int l_supp_i = nodeSupport.MinorDim(i);
		iArrayT supp_i(l_supp_i);
		supp_i.Copy(nodeSupport(i));
		supp_i.SortAscending();
		nodeWorkSpace[i].AppendArray(l_supp_i, supp_i.Pointer());
		facetWorkSpace[i].AppendArray(l_supp_i, zeroFacet);
	}

	/* integration */
	int nfn = 2;
	int nsd = 2;
	ParentDomainT domain(GeometryT::kLine, fNumIP, nfn);
	domain.Initialize();
	LocalArrayT facet_coords(LocalArrayT::kInitCoords, nfn, nsd);
	facet_coords.SetGlobal(fVoronoiVertices);
	iArrayT keys;
	dArrayT ip_coords(nsd);
	dMatrixT jacobian(nsd, 1);
	const double* ip_weight = domain.Weight();

	dArrayT facetNormal(fSD), facetIntegral(fSD);
	double* currentB, *currentI;
	int n_0, n_1;
	bool traverseQ_0, traverseQ_1;
	int *next_0, *next_1;
	for (int i = 0; i < fDeloneEdges.MajorDim(); i++) {
		n_0 = fDeloneEdges(i,0);
		n_1 = fDeloneEdges(i,1); 
		
		facetNormal.DiffOf(fDeloneVertices(n_1), fDeloneVertices(n_0));
		facetNormal.UnitVector();

		/* copy face coordinates with local ordering */
		fDualFacets.RowAlias(i,keys);
		facet_coords.SetLocal(keys);
		for (int ii = 0; ii < fNumIP; ii++) {
		
			/* jacobian of the coordinate transformation */
			domain.DomainJacobian(facet_coords, ii, jacobian);
			double jw = ip_weight[ii]*domain.SurfaceJacobian(jacobian);

			/* integration point coordinates */
			domain.Interpolate(facet_coords, ip_coords, ii);			

			if (!fNodalShapes->SetFieldAt(ip_coords, NULL)) // shift = 0 or not ?
				ExceptionT::GeneralFail("SCNIMFT::ComputeBMatrices","Shape Function evaluation"
					"failed at Delone edge %d\n",i);
				
			const dArrayT& phiValues = fNodalShapes->FieldAt();			
		
			iArrayT ip_cover(fNodalShapes->Neighbors());	
			int n_ip_cover = ip_cover.Length();	
			iArrayT ip_cover_key(n_ip_cover);
			ip_cover_key.SetValueToPosition();
			ip_cover_key.SortAscending(ip_cover);

			LinkedListT<int>& supp_0 = nodeWorkSpace[n_0];
			LinkedListT<int>& supp_1 = nodeWorkSpace[n_1];
			LinkedListT< dArrayT >& bVectors_0 = facetWorkSpace[n_0];
			LinkedListT< dArrayT >& bVectors_1 = facetWorkSpace[n_1];
			int s_0 = -1;
			int s_1 = -1;
			/* Simultaneously loop over support of the two nodes that are endpoints of the
			 * current Delone edge and the nodes in the support of the midpoint of this
			 * edge. If a node covering the integration point is not in the support of n_0 
			 * or n_1, insert that covering node into the sorted list.
			 */
			 
			int* c = ip_cover.Pointer();
			int* c_j = ip_cover_key.Pointer();
			
			supp_0.Top(); bVectors_0.Top();
			supp_1.Top(); bVectors_1.Top();
			next_0 = supp_0.CurrentValue();
			next_1 = supp_1.CurrentValue();
			for (int j = 0; j < n_ip_cover; j++, c++, c_j++) {
		
				facetIntegral = facetNormal;
				facetIntegral *= phiValues[*c_j]*jw;	
				
				if (next_0)
					traverseQ_0 = *next_0 <= *c;
				else
					traverseQ_0 = false;
						
				// advance supp_0 and supp_1 until they are greater than or equal to current node
				while (traverseQ_0 && supp_0.Next(s_0) && bVectors_0.Next()) {
					next_0 = supp_0.PeekAhead(); 
					if (!next_0)
						traverseQ_0 = false;
					else
						if (*next_0 > *c)
							traverseQ_0 = false;
				}
					
				if (s_0 != *c) { // means we're not at the end of the linked list
					supp_0.InsertAtCurrent(*c);
					bVectors_0.InsertAtCurrent(zeroFacet);
					s_0 = *c;
					if (supp_0.AtTop()) { // if we're inserting at the front, LinkedListT's behavior requires more work
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
				while (traverseQ_1 && supp_1.Next(s_1) && bVectors_1.Next()) {
					next_1 = supp_1.PeekAhead(); 
					if (!next_1)
						traverseQ_1 = false;
					else
						if (*next_1 > *c)
							traverseQ_1 = false;
				}		
									
				if (s_1 != *c) {
					supp_1.InsertAtCurrent(*c);
					bVectors_1.InsertAtCurrent(zeroFacet);
					s_1 = *c;
					if (supp_1.AtTop()) { // if we're inserting at the front, LinkedListT's behavior requires more work
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
	}
#ifdef SECOND_WAY
	/** temporary storage for integration over the body boundary */
	ArrayT< LinkedListT<double> > boundary_phi(fNonDeloneEdges.Length());
	ArrayT< LinkedListT<int> > boundary_supports(fNonDeloneEdges.Length());
	
	double zero = 0.0;
	dArrayT boundaryIPCoord(fSD);
	for (int i = 0; i < fNonDeloneEdges.Length(); i++) {
		double* v1 = fVoronoiVertices(fSelfDualFacets(i,0));
		double* v2 = fVoronoiVertices(fSelfDualFacets(i,1));
		for (int j = 0; j < fSD; j++)
			boundaryIPCoord[j] = v1[j] + v2[j];
		boundaryIPCoord /= double(fSD);
		
		if (!fNodalShapes->SetFieldAt(boundaryIPCoord, NULL)) // shift = 0 or not ?
				ExceptionT::GeneralFail("SCNIMFT::ComputeBMatrices","Shape Function evaluation"
					"failed at Delone edge %d\n",i);
					
		const dArrayT& phiValues = fNodalShapes->FieldAt();			
					
		iArrayT supp_i(fNodalShapes->Neighbors());	
		int l_supp_i = supp_i.Length();
		supp_i.SortAscending();
		boundary_supports[i].AppendArray(l_supp_i, supp_i.Pointer());
		boundary_phi[i].AppendArray(l_supp_i, zero);
	}
#endif	
	/** Loop over remaining edges */
	for (int i = 0; i < fNonDeloneEdges.Length(); i++) {
		n_0 = fNonDeloneEdges[i];
		facetNormal.Set(fSD, fNonDeloneNormals(i));
		facetNormal.UnitVector();
		
#ifdef SEPARATE_BOUNDARY_INTEGRATION
		int num_bdry_pts = 1; // set this for trapezoidal rule
#ifdef EVALUATE_AT_NODES
		num_bdry_pts = 1;
#endif
		double jw = fBoundaryIntegrationWeights[i]/(num_bdry_pts);
		double dl = 1./double(num_bdry_pts);
		double* v1 = fVoronoiVertices(fSelfDualFacets(i,0));
		double* v2 = fVoronoiVertices(fSelfDualFacets(i,1));
		dArrayT ip_coord0(fSD,v1);
		dArrayT edgeVector(fSD);
		edgeVector.DiffOf(v2,v1);
#endif
		
		/* copy face coordinates with local ordering */
		fSelfDualFacets.RowAlias(i,keys);
		facet_coords.SetLocal(keys);
#ifndef SEPARATE_BOUNDARY_INTEGRATION
		for (int ii = 0; ii < fNumIP; ii++) {
			/* jacobian of the coordinate transformation */
			domain.DomainJacobian(facet_coords, ii, jacobian);
			double jw = ip_weight[ii]*domain.SurfaceJacobian(jacobian);

			/* integration point coordinates */
			domain.Interpolate(facet_coords, ip_coords, ii);	
#else
		for (int ii = 1; ii <= num_bdry_pts; ii++) {
#ifndef EVALUATE_AT_NODES
			ip_coords.SetToCombination(1.,ip_coord0,(ii-.5)*dl,edgeVector);
#else
			ip_coords.Set(fSD,fDeloneVertices(n_0));
#endif
#endif // SEPARATE_BOUNDARY_INTEGRATION				

			if (!fNodalShapes->SetFieldAt(ip_coords, NULL)) // shift = 0 or not ?
				ExceptionT::GeneralFail("SCNIMFT::ComputeBMatrices","Shape Function evaluation"
					"failed at Delone edge %d\n",i);
					
			const dArrayT& phiValues = fNodalShapes->FieldAt();			
					
			iArrayT ip_cover(fNodalShapes->Neighbors());	
			int n_ip_cover = ip_cover.Length();	
			iArrayT ip_cover_key(n_ip_cover);
			ip_cover_key.SetValueToPosition();
			ip_cover_key.SortAscending(ip_cover);
			
			LinkedListT<int>& supp_0 = nodeWorkSpace[n_0];
			LinkedListT< dArrayT >& bVectors_0 = facetWorkSpace[n_0];
			int s_0;
			
			/* Merge support of the boundary node with covering of integration point
			 */
			int* c = ip_cover.Pointer();
			int* c_j = ip_cover_key.Pointer();
			
			supp_0.Top(); bVectors_0.Top();
			next_0 = supp_0.CurrentValue();
			for (int j = 0; j < n_ip_cover; j++, c++, c_j++) {
				facetIntegral = facetNormal;
				facetIntegral *= phiValues[*c_j]*jw;		
			
				if (next_0)
					traverseQ_0 = *next_0 <= *c;
				else
					traverseQ_0 = false;
						
				// advance supp_0 and supp_1 until they are greater than or equal to current node
				while (traverseQ_0 && supp_0.Next(s_0) && bVectors_0.Next()) {
					next_0 = supp_0.PeekAhead(); 
					if (!next_0)
						traverseQ_0 = false;
					else
						if (*next_0 > *c)
							traverseQ_0 = false;
				}
					
				if (s_0 != *c) { // means we're not at the end of the linked list
					supp_0.InsertAtCurrent(*c);
					bVectors_0.InsertAtCurrent(zeroFacet);
					s_0 = *c;
					if (supp_0.AtTop()) { // if we're inserting at the front, LinkedListT's behavior requires more work
						supp_0.Next(); 
						bVectors_0.Next();
					}
				}
					
				currentI = facetIntegral.Pointer();
				currentB =  bVectors_0.CurrentValue()->Pointer();

				for (int k = 0; k < fSD; k++)
					*currentB++ += *currentI++;
			}
#ifdef SECOND_WAY			
			LinkedListT<int>& bsupp_0 = boundary_supports[i];
			LinkedListT< double >& phi_0 = boundary_phi[i];
			
			/* Merge support of the boundary facet with covering of integration point
			 */
			c = ip_cover.Pointer();
			c_j = ip_cover_key.Pointer();
			
			bsupp_0.Top(); phi_0.Top();
			next_0 = bsupp_0.CurrentValue();
			for (int j = 0; j < n_ip_cover; j++, c++, c_j++) {
				facetIntegral[0] = phiValues[*c_j]*jw;		
				if (next_0)
					traverseQ_0 = *next_0 <= *c;
				else
					traverseQ_0 = false;
						
				// advance supp_0 and supp_1 until they are greater than or equal to current node
				while (traverseQ_0 && bsupp_0.Next(s_0) && phi_0.Next()) {
					next_0 = supp_0.PeekAhead(); 
					if (!next_0)
						traverseQ_0 = false;
					else
						if (*next_0 > *c)
							traverseQ_0 = false;
				}
					
				if (s_0 != *c) { // means we're not at the end of the linked list
					bsupp_0.InsertAtCurrent(*c);
					phi_0.InsertAtCurrent(zero);
					s_0 = *c;
					if (supp_0.AtTop()) { // if we're inserting at the front, LinkedListT's behavior requires more work
						bsupp_0.Next(); 
						phi_0.Next();
					}
				}

				*(phi_0.CurrentValue()) += facetIntegral[0];
			}
#endif
		}	
	}
 
 	// scale integrals by volumes of Voronoi cells
	dArrayT* currFacetIntegral;
	for (int i = 0; i < nNodes; i++) {
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
#ifdef SECOND_WAY
	//do the same for the boundary integrals
	
	fBoundaryPhi.Configure(boundary_phi);
	fBoundarySupports.Configure(boundary_supports);
	
	for (int i = 0; i < boundary_supports.Length(); i++) {
		int* irow_i = fBoundarySupports(i);
		double* drow_i = fBoundaryPhi(i);
		LinkedListT<int>& ilist = boundary_supports[i];
		LinkedListT<double>& dlist = boundary_phi[i];
		ilist.Top(); dlist.Top();
		while (ilist.Next() && dlist.Next()) {
			*irow_i++ = *(ilist.CurrentValue());
			*drow_i++ = *(dlist.CurrentValue());
		}
	}

#endif
}

void SCNIMFT::VoronoiDiagramToFile(ofstreamT& vout)
{
	
    int nVertices = fVoronoiVertices.MajorDim();
    int nSD = fVoronoiVertices.MinorDim();

    vout << fNodes.Length() << "\n";
    vout << nSD << "\n";
    vout << nVertices << "\n";

    // write out vertices of the VoronoiDiagram
    for (int i = 0; i < nVertices; i++) {
      	for (int j = 0; j < nSD; j++)
			vout << fVoronoiVertices(i,j) << " "; 
      	vout << "\n";
    }
    
    // write out Voronoi cells for each node
    for (int i = 0; i < fNodes.Length(); i++) {
		vout << i <<"\n";

		for (int j = 0; j < fSD; j++)
			vout << fVoronoiCellCentroids(i,j) << " ";

		// cell volume
		vout << fVoronoiCellVolumes[i] << "\n";
    }

    // write out Delone edge information
    int nDelone = fDeloneEdges.MajorDim();
    vout << nDelone << "\n"; // number of Delone edges
    for (int i = 0; i < nDelone; i++) {
    	vout << fDeloneEdges(i,0) << " " << fDeloneEdges(i,1) << "\n";
    }

    // write out areas and centroids of Voronoi facets dual to Delone edge
    // this assumes a 1-point integration scheme over the boundary of each node's cell
    if (fDualFacets.MajorDim() != nDelone)
    	ExceptionT::GeneralFail("SCNIMFT::VoronoiDiagramToFile","Dual edge/facet dimension mismatch\n");

    for (int i = 0; i < fDualFacets.MajorDim(); i++) {
    	vout << fDualFacets(i,0) << " " << fDualFacets(i,1) << " ";
    	vout << "\n";
    }
  
    // write out self-duals and allocate storage for them
    // self-duals are facets on the body boundary dual to only 1 node in the body
    int numSelfDuals = fNonDeloneEdges.Length();
    vout << numSelfDuals << "\n";
    for (int i = 0; i < numSelfDuals; i++)
      vout << fNonDeloneEdges[i] << " ";
    vout << "\n";
    
    for (int i = 0; i < numSelfDuals; i++) { 
      for (int k = 0; k < fSD; k++)
	    vout << fSelfDualFacets(i,k) << " "; // SPECIALIZED TO 2D!!!
      vout << "\n";
    }

    // list of normals of self-dual facets
    for (int i = 0; i < numSelfDuals; i++) {
      for (int k = 0; k < fSD; k++) {
	    vout << fNonDeloneNormals(i,k) << " ";
      }
      vout << "\n";
    }

    // list of areas of self-dual facets
    for (int i = 0; i < numSelfDuals; i++)
      vout << fBoundaryIntegrationWeights[i] << "\n";
 
}	
	
void SCNIMFT::VoronoiDiagramFromFile(ifstreamT& vin)
{
	const char caller[] = "SCNIMFT::VoronoiDiagramFromFile";	

    int nNodes, nSD, nVertices;
    vin >> nNodes;
    
    /* minor consistency checks */
    if (nNodes != fNodes.Length()) {
		vin.close();
		ExceptionT::GeneralFail(caller,"Input Voronoi file does not match node number\n");
    }
    
    vin >> nSD;
    if (nSD != fSD) {
		vin.close();
		ExceptionT::GeneralFail(caller,"Input Voronoi file does not match SD\n");
    }
    
    vin >> nVertices;
    
    // allocate memory for Voronoi diagram data structures
    fVoronoiVertices.Dimension(nVertices, nSD);
    fVoronoiCellCentroids.Dimension(nNodes, nSD);
    fVoronoiCellVolumes.Dimension(nNodes);

    for (int i = 0 ; i < nVertices; i++)
    	for (int j = 0; j < nSD; j++)
	  		vin >> fVoronoiVertices(i,j);

    for (int i = 0; i < nNodes; i++) {
		int itmp;
		vin >> itmp;
		if (itmp != i)
	  		ExceptionT::GeneralFail(caller,"Bad Input Voronoi file\n");

	  	for (int j = 0; j < fSD; j++)
			vin >> fVoronoiCellCentroids(i,j);
	
		vin >> fVoronoiCellVolumes[i];

	}
	
	/* Read in all DeloneEdges in or on the body */
	int nDelone;
	vin >> nDelone; 
	fDeloneEdges.Dimension(nDelone, 2);
	fDualFacets.Dimension(nDelone, 2);
	
	for (int i = 0; i < nDelone; i++)
		vin >> fDeloneEdges(i,0) >> fDeloneEdges(i,1);
		
	for (int i = 0; i < nDelone; i++) {
	  vin >> fDualFacets(i,0) >> fDualFacets(i,1);
	}
		
	int nCentroids;
	vin >> nCentroids; // number of boundary facets (self-duals)
	fNonDeloneEdges.Dimension(nCentroids);
	fNonDeloneNormals.Dimension(nCentroids,fSD);
	fBoundaryIntegrationWeights.Dimension(nCentroids);

	fSelfDualFacets.Dimension(nCentroids, 2);
	
	for (int i = 0; i < nCentroids; i++)
		vin >> fNonDeloneEdges[i];	

	for (int i = 0; i < nCentroids; i++) {
	  vin >> fSelfDualFacets(i,0) >> fSelfDualFacets(i,1);
	}
			
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

	ParameterT compute_voronoi(qComputeVoronoiCell, "compute_voronoi");
	compute_voronoi.SetDefault(qComputeVoronoiCell);
	list.AddParameter(compute_voronoi);
	
	ParameterT voronoi_file(vCellFile, "voronoi_file");
	voronoi_file.SetDefault(vCellFile);
	list.AddParameter(voronoi_file);

	ParameterT num_ip(fNumIP, "num_ip");	
	num_ip.SetDefault(fNumIP);
	list.AddParameter(num_ip);

	ParameterT just_voronoi_diagram(qJustVoronoiDiagram,"just_voronoi_diagram");
	just_voronoi_diagram.SetDefault(qJustVoronoiDiagram);
	list.AddParameter(just_voronoi_diagram);
}

/* information about subordinate parameter lists */
void SCNIMFT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ElementBaseT::DefineSubs(sub_list);
	
	/* parameters for the meshfree support */
	sub_list.AddSub("meshfree_support_choice", ParameterListT::Once, true);

	/* list of node set ID's defining which nodes get integrated */
	sub_list.AddSub("mf_particle_ID_list", ParameterListT::OnePlus);
	
	/* optional body force */
	sub_list.AddSub("body_force", ParameterListT::ZeroOrOnce);
	
	/* tractions */
	sub_list.AddSub("natural_bc", ParameterListT::Any);
}

/* return the description of the given inline subordinate parameter list */
void SCNIMFT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	/* inherited */
	ElementBaseT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SCNIMFT::NewSub(const StringT& name) const
{
	SCNIMFT* non_const_this = const_cast<SCNIMFT*>(this);
	
	/* try material list */
	MaterialListT* material_list = non_const_this->NewMaterialList(name, 0);

	if (material_list)
		return material_list;
	else if (name == "meshfree_support_choice") {
	  ParameterContainerT* mf_choice = new ParameterContainerT(name);
	  mf_choice->SetSubSource(this);
	  mf_choice->SetListOrder(ParameterListT::Choice);
	  
	  mf_choice->AddSub("meshfree_support_2D");
	  mf_choice->AddSub("meshfree_support_3D");
	  
	  return mf_choice;
  	} else if (name == "meshfree_support_2D")
	  return new MeshFreeSupport2DT;	
	else if (name == "meshfree_support_3D")
	  return new MeshFreeSupport3DT;
  	else if (name == "body_force") { // body force
		ParameterContainerT* body_force = new ParameterContainerT(name);
	
		/* schedule number */
		body_force->AddParameter(ParameterT::Integer, "schedule");
	
		/* body force vector */
		body_force->AddSub("Double", ParameterListT::OnePlus); 		
		
		return body_force;
	} else if (name == "natural_bc") { /* traction bc */
		ParameterContainerT* natural_bc = new ParameterContainerT(name);

		natural_bc->AddParameter(ParameterT::Word, "side_set_ID");
		natural_bc->AddParameter(ParameterT::Integer, "schedule");

		ParameterT coord_sys(ParameterT::Enumeration, "coordinate_system");
		coord_sys.AddEnumeration("global", Traction_CardT::kCartesian);
		coord_sys.AddEnumeration( "local", Traction_CardT::kLocal);
		coord_sys.SetDefault(Traction_CardT::kCartesian);
		natural_bc->AddParameter(coord_sys);

		natural_bc->AddSub("DoubleList", ParameterListT::OnePlus); 		
		
		return natural_bc;
	} else /* inherited */
    	return ElementBaseT::NewSub(name);
}
