/* $Id: SCNIMFT.cpp,v 1.50 2005-01-26 20:21:00 cjkimme Exp $ */
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
#include "ParameterUtils.h"
#include "ParameterContainerT.h"
#include "InverseMapT.h"

#include "MeshFreeSupport2DT.h"
#include "MeshFreeSupport3DT.h"
#include "MeshFreeNodalShapeFunctionT.h"
#include "ContinuumMaterialT.h"
#include "SolidMaterialT.h"
#include "SolidMatSupportT.h"

/* cell geometries */
#include "CellGeometryT.h"
#include "VoronoiDiagramT.h"
#include "CellFromMeshT.h"

#include "Traction_CardT.h"

/* materials lists */
#include "SSSolidMatList1DT.h"
#include "SSSolidMatList2DT.h"
#include "SSSolidMatList3DT.h"

using namespace Tahoe;

const int kNoTractionVector = -1;

/* constructors */
SCNIMFT::SCNIMFT(const ElementSupportT& support, const FieldT& field):
	ElementBaseT(support),
	fSD(ElementSupport().NumSD()),
	fMaterialList(NULL),
	fNodalShapes(NULL),
	fCellGeometry(NULL),
	qIsAxisymmetric(false)
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
	fNodalShapes(NULL),
	fCellGeometry(NULL),
	qIsAxisymmetric(false)
{
	SetName("mfparticle");

	/* set matrix format */
	fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
}

/* destructor */
SCNIMFT::~SCNIMFT(void)
{
	
	delete fCellGeometry;
	delete fMaterialList;
	delete fNodalShapes;
}

/* initialization */
void SCNIMFT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SCNIMFT::TakeParameterList";
	
	/* What dimension? */
	fSD = ElementSupport().NumSD();
	
	/* get parameters needed to construct shape functions */
	fMeshfreeParameters = list.ListChoice(*this, "meshfree_support_choice");
	
	/* access to the model database */
	ModelManagerT& model = ElementSupport().ModelManager();

	if (fSD == 2) {
		/* extract particle ID's */
		const ParameterListT& particle_ID_params = list.GetList("mf_particle_ID_list");
		ArrayT<StringT> particle_ID_list;
		StringListT::Extract(particle_ID_params, particle_ID_list);

		//get nodes from ModelManagerT
		model.ManyNodeSets(particle_ID_list, fNodes);
	}
	
	// get the cell geometry parameters
	const ParameterListT* geometry_params = list.ListChoice(*this, "cell_geometry_choice");
	if (geometry_params->Name() == "voronoi_diagram")
		fCellGeometry = new VoronoiDiagramT(ElementSupport(), qIsAxisymmetric);
	else if (geometry_params->Name() == "cell_from_mesh")
		fCellGeometry = new CellFromMeshT(ElementSupport(), qIsAxisymmetric);
	else
		ExceptionT::GeneralFail(caller,"Cannot get valid cell geometry from input file\n");
	fCellGeometry->TakeParameterList(*geometry_params);
	fCellGeometry->SetNodalElements(this);
	
	/* inherited */
	ElementBaseT::TakeParameterList(list);

	/* re-dimension "element" force and stiffness contributions */
	fLHS.Dimension(fSD);
	
	/* allocate work space */
	fForce_man.SetWard(0, fForce, fSD);
	fForce_man.SetMajorDimension(ElementSupport().NumNodes(), false);

	/* write parameters */
	ostream& out = ElementSupport().Output();
	
	/* shape functions */
	/* only support single list of integration cells for now */
	if (fElementConnectivities.Length() > 1) {
	       ExceptionT::GeneralFail(caller,"Multiple ElementConnectivities not yet supported\n");
	}

	/* construct shape functions */
	// see what to pass in place of fNodalCoordinates
	fNodalShapes = new MeshFreeNodalShapeFunctionT(fSD,
		ElementSupport().InitialCoordinates(), *fElementConnectivities[0], 
		fNodalCoordinates, *fMeshfreeParameters);
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
	fCellGeometry->SetNodesAndShapes(fNodes, fNodalCoordinates, fNodalShapes);
	fCellGeometry->ComputeBMatrices(nodalCellSupports, bVectorArray, fCellVolumes, 
									fCellCentroids, circumferential_B);	
	
	/** store shape functions at nodes */
	int nNodes = fNodes.Length();
	dArrayT nodalCoords;
	ArrayT< LinkedListT<double> > nodal_phi;
	ArrayT< LinkedListT<int> > nodal_supports;
	nodal_phi.Dimension(nNodes);
	nodal_supports.Dimension(nNodes);
	for (int i = 0; i < nNodes; i++) {
		nodalCoords.Set(fSD, fNodalCoordinates(i));
	
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
	
	// store shape function information for boundary integration
	fCellGeometry->BoundaryShapeFunctions(fBoundaryPhi, fBoundarySupports, fBoundaryFacetNormals);
	
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
	fTractionBoundaryCondition.Dimension(fBoundaryFacetNormals.MajorDim());
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
	fPointConnectivities.Alias(fNodes.Length(), 1, fNodes.Pointer());

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

	const ArrayT<StringT>& disp_labels = Field().Labels();
	for (int ns = 0 ; ns < fSD; ns++)
	  	labels[dex++] = disp_labels[ns];

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

	// Get nodal coordinates 
	fNodalCoordinates.Dimension(fNodes.Length(), fSD);
	fNodalCoordinates.RowCollect(fNodes, model.Coordinates());

	/* set up element cards for state variable storage */
	fElementCards.Dimension(fNodes.Length()); /* one card per node */
	
	fCellGeometry->DefineElements(block_ID, mat_index);
}

/* collecting element group equation numbers */
void SCNIMFT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_1)

	/* dimension equations array */
	fEqnos.Configure(nodalCellSupports, NumDOF());

	/* get local equations numbers */
	Field().SetLocalEqnos(nodalCellSupports, fEqnos);

	/* add to list of equation numbers */
	eq_2.Append(&fEqnos);
}

/* assemble particle mass matrix into LHS of global equation system */
void SCNIMFT::AssembleParticleMass(const double rho)
{
	
	fForce = 0.0;
	int* nodes = fNodes.Pointer();
	double* volume = fCellVolumes.Pointer();
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
	//fBoundaryIntegrationWeights and fNonDeloneEdges have weights and node numbers
	//I could also keep an fTractionBoundaryCondition which would have the number of
	//the traction boundary condition or -1 for none
	//that's the kludgiest way to get it going for now
	if (fTractionVectors.MajorDim()) {
	
		fForce = 0.;
		int numBoundaryFacets = fBoundaryFacetNormals.MajorDim();
		for (int i = 0; i < numBoundaryFacets; i++) 
			if (fTractionBoundaryCondition[i] != kNoTractionVector) {
				dArrayT workspace(fSD == 2 ? 3 : 6); 
				dArrayT traction_vector(fSD);
				fTractionVectors.RowCopy(fTractionBoundaryCondition[i], workspace.Pointer());
				traction_vector[0] = workspace[0]*fBoundaryFacetNormals(i,0) + .5*workspace[2]*fBoundaryFacetNormals(i,1);
				traction_vector[1] = workspace[1]*fBoundaryFacetNormals(i,1) + .5*workspace[2]*fBoundaryFacetNormals(i,0);
 			
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
	
	/* parameters for the meshfree support */
	sub_list.AddSub("cell_geometry_choice", ParameterListT::Once, true);
	
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
	else if (name == "cell_geometry_choice") {
	  ParameterContainerT* cg_choice = new ParameterContainerT(name);
	  cg_choice->SetSubSource(this);
	  cg_choice->SetListOrder(ParameterListT::Choice);
	  
	  cg_choice->AddSub("voronoi_diagram");
	  cg_choice->AddSub("cell_from_mesh");
	  
	  return cg_choice;
  	} else if (name == "voronoi_diagram")
  		return new VoronoiDiagramT;
  	else if (name == "cell_from_mesh")
  		return new CellFromMeshT;
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
