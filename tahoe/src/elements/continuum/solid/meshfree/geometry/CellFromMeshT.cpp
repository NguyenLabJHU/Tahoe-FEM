/* $Id: CellFromMeshT.cpp,v 1.4 2005-01-27 17:52:45 paklein Exp $ */
#include "CellFromMeshT.h"

#include "ArrayT.h"
#include "dArrayT.h"
#include "LinkedListT.h"
#include "ModelManagerT.h"
#include "ShapeFunctionT.h"

using namespace Tahoe;

/* constructors */
CellFromMeshT::CellFromMeshT(const ElementSupportT& support, bool isAxisymmetric):
	CellGeometryT(support, isAxisymmetric)
{
	SetName("cell_from_mesh");
}

/* constructors */
CellFromMeshT::CellFromMeshT(void)
{
	SetName("cell_from_mesh");
}

void CellFromMeshT::ComputeBMatrices(RaggedArray2DT<int>& cellSupports, RaggedArray2DT<dArrayT>& bVectors,
	dArrayT& cellVolumes, dArray2DT& cellCentroids, RaggedArray2DT<double>& circumferential_B)
{
	const char caller[] = "CellFromMeshT::ComputeBMatrices";

#pragma unused(cellSupports)
#pragma unused(bVectors)
#pragma unused(cellVolumes)
#pragma unused(circumferential_B)

	/* For the Axisymmetric case, also calculates {Psi_I(X_L)/R)L,0.} */

	const RaggedArray2DT<int>& nodeSupport = fNodalShapes->NodeNeighbors();
	
	int nSD = fElementSupport->NumSD();
//	int nElements = fElementConnectivies.MajorDim();
	int nNodes = fNodalCoordinates->MajorDim();
//	int numFacets = numSD == 2 ? 4 : 6;
	
	ArrayT< LinkedListT<int> > nodeWorkSpace; 
	ArrayT< LinkedListT<dArrayT> > facetWorkSpace; 
	ArrayT< LinkedListT<double> > circumferentialWorkSpace;
	nodeWorkSpace.Dimension(nNodes);
	facetWorkSpace.Dimension(nNodes);
	if (qIsAxisymmetric) {
		circumferentialWorkSpace.Dimension(nNodes);
	}

	dArrayT zeroFacet(nSD);
	double zeroSingle = 0.;
	zeroFacet = 0.0;
	for (int i = 0; i < nNodes; i++) {
		int node_i = (*fNodes)[i];
		int l_supp_i = nodeSupport.MinorDim(node_i);
		iArrayT supp_i(l_supp_i);
		supp_i.Copy(nodeSupport(node_i));
		supp_i.SortAscending();
		nodeWorkSpace[i].AppendArray(l_supp_i, supp_i.Pointer());
		facetWorkSpace[i].AppendArray(l_supp_i, zeroFacet);
		if (qIsAxisymmetric)
			circumferentialWorkSpace[i].AppendArray(l_supp_i, zeroSingle);
	}

	/* initialize cell volumes */
	cellVolumes.Dimension(nNodes);
	cellVolumes = 0.0;
	cellCentroids.Dimension(nNodes, nSD);
	cellCentroids = 0.0;
	
	/* model information */
	ModelManagerT& model = ElementSupport().ModelManager();
	ArrayT<const iArray2DT*> connectivities;
	model.ElementGroupPointers(fBlockID, connectivities);
	int nElementNodes = connectivities[0]->MinorDim();

	/* determine cell geometry */
	GeometryT::CodeT cell_geometry = GeometryT::kNone;
	for (int i = 0; i < fBlockID.Length(); i++) {
		GeometryT::CodeT geometry = model.ElementGroupGeometry(fBlockID[i]);
		if (cell_geometry != GeometryT::kNone && cell_geometry != geometry)
			ExceptionT::GeneralFail(caller, "all cell geometries must be the same");
		cell_geometry = geometry;
	}

	/* shape functions over cell volume */
	LocalArrayT cell_coords(LocalArrayT::kInitCoords, nElementNodes, nSD);
	cell_coords.SetGlobal(ElementSupport().InitialCoordinates());
	ShapeFunctionT cell_shape(cell_geometry, 1, cell_coords);
	cell_shape.Initialize();
	const ParentDomainT& cell_parent_domain = cell_shape.ParentDomain();

	/* subdomain volume */
	GeometryT::CodeT sub_cell_geometry = cell_parent_domain.NodalSubDomainGeometry();
	int sub_cell_nodes = cell_parent_domain.NodalSubDomainNumPoints();
	LocalArrayT sub_cell_coords(LocalArrayT::kInitCoords, sub_cell_nodes, nSD);
	ShapeFunctionT sub_cell_shape(sub_cell_geometry, 1, sub_cell_coords);
	sub_cell_shape.Initialize();
	const ParentDomainT& sub_cell_parent_domain = sub_cell_shape.ParentDomain();

	/* subdomain boundaries */
	int n_faces = sub_cell_shape.NumFacets();
	ArrayT<GeometryT::CodeT> facet_geom(n_faces);
	iArrayT nfn(n_faces);
	sub_cell_shape.FacetGeometry(facet_geom, nfn);
	ParentDomainT sub_cell_face_domain(facet_geom[0], fNumIP, nfn[0]); /* assume all faces are the same */
	sub_cell_face_domain.Initialize();
	LocalArrayT facet_coords(LocalArrayT::kUnspecified, nfn[0], nSD);
	iArrayT facet_nodes(nfn[0]);
	const double* ip_weight = sub_cell_face_domain.Weight();

	//facet_coords.SetGlobal(fVoronoiVertices);
//	iArrayT keys;
	dArrayT ip_coords(nSD);
	dMatrixT jacobian(nSD, 1);
//	const double* ip_weight = domain.Weight();

	dMatrixT Q(nSD);
	iArrayT nodes_glb(nElementNodes), nodes_loc(nElementNodes);
	dArrayT facetNormal(nSD), facetIntegral(nSD);
	double* currentB, *currentI;
	int n_0, n_1;
	bool traverseQ_0, traverseQ_1;
	int *next_0, *next_1;
	for (int e = 0; e < connectivities.Length(); e++) /* loop over element blocks */
	{
		/* block connectivities */
		const iArray2DT& connects = *(connectivities[e]);
		int nElements = connects.MajorDim();
	
		for (int i = 0; i < nElements; i++) { /* loop over elements in the block */

			/* element nodes */
			connects.RowAlias(i, nodes_glb);
			nodes_loc = nodes_glb;
			fscnimft->GlobalToLocalNumbering(nodes_loc);

			/* collect coordinates over the current element */
			cell_coords.SetLocal(nodes_glb);

			//iArrayT& element_nodes = ; // global numbers of element nodes ? 
			//fscnimft->GlobalToLocalNumbering(element_nodes); // this function is available
		
			for (int k = 0; k < nElementNodes; k++) { /* loop over nodal subdomains */

				/* meshfree node number */
				n_0 = nodes_loc[k];

				/* subdomain coordinates */
				cell_parent_domain.NodalSubDomainCoordinates(cell_coords, k, sub_cell_coords);

				/* shape functions over the subdomain */
				sub_cell_shape.SetDerivatives();
				
				/* volumetric integration factors */
				const double* sub_cell_det = sub_cell_shape.IPDets();
				const double* sub_cell_wgt = sub_cell_shape.IPWeights();
				
				/* compute volume and centroid */
				sub_cell_shape.TopIP();
				while (sub_cell_shape.NextIP()) 
				{
					int ip = sub_cell_shape.CurrIP();
					double dv = sub_cell_det[ip]*sub_cell_wgt[ip];
				
					/* integrate nodal volume */
					cellVolumes[n_0] += dv;
					
					/* integrate centroid */
					sub_cell_shape.IPCoords(ip_coords);
					cellCentroids.AddToRowScaled(n_0, dv, ip_coords);
				}
				cellCentroids.ScaleRow(n_0, 1.0/cellVolumes[n_0]);

				/* loop over facets (internal and on the element boundary) for each node */
				for (int jj = 0; jj < n_faces; jj++) {
	
					/* collect coordinates over the facet */
					sub_cell_shape.NodesOnFacet(jj, facet_nodes);
					facet_coords.Collect(facet_nodes, sub_cell_coords);
	
					//facetNormal = ?????;
					//facetNormal.UnitVector();

					/* integrate over the face */
					int nIP = sub_cell_face_domain.NumIP();
					for (int ii = 0; ii < nIP; ii++) {
		
						/* jacobian of the coordinate transformation */
						sub_cell_face_domain.DomainJacobian(facet_coords, ii, jacobian);
						double jw = ip_weight[ii]*sub_cell_face_domain.SurfaceJacobian(jacobian, Q);
	
						/* surface normal - last column on transformation tensor */
						Q.CopyColumn(nSD-1, facetNormal);
	
						/* integration point coordinates */
						sub_cell_face_domain.Interpolate(facet_coords, ip_coords, ii);			

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
						LinkedListT< double > *circumf_0;
						if (qIsAxisymmetric) 
							circumf_0 = &circumferentialWorkSpace[n_0];
						int s_0 = -1;

						/* Loop over support of the current element node.*/
					 
						int* c = ip_cover.Pointer();
						int* c_j = ip_cover_key.Pointer();
			
						supp_0.Top(); bVectors_0.Top();
						next_0 = supp_0.CurrentValue();
						if (qIsAxisymmetric) 
							circumf_0->Top();
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
								if (qIsAxisymmetric) 
									circumf_0->InsertAtCurrent(0.);
								s_0 = *c;
								if (supp_0.AtTop()) { // if we're inserting at the front, LinkedListT's behavior requires more work
									supp_0.Next(); 
									bVectors_0.Next();
								}
							}
							
							currentI = facetIntegral.Pointer();
							currentB = bVectors_0.CurrentValue()->Pointer();
							for (int k = 0; k < nSD; k++)
								*currentB++ += *currentI++;
						}
					}
				}		
			}
		}
	}
 
 	// scale integrals by volumes of cells
	dArrayT* currFacetIntegral;
	for (int i = 0; i < nNodes; i++) {
		LinkedListT<dArrayT>& bVectors_i = facetWorkSpace[i];
		LinkedListT<int>& nodes_i = nodeWorkSpace[i];
		bVectors_i.Top(); nodes_i.Top();
		while ((currFacetIntegral = bVectors_i.Next()))
			*currFacetIntegral /= cellVolumes[i];
	}
	
	if (qIsAxisymmetric) {
		// calculate Psi/R terms. These are evaluated nodally, so this additional loop
		// is required
		dArrayT phis, nodal_init_coords;
		for (int i = 0; i < fNodalCoordinates->MajorDim(); i++) {
			nodal_init_coords.Alias(nSD, (*fNodalCoordinates)(i)); // This is the nodal coordinate.
			double R_i = nodal_init_coords[0];
			
			if (R_i >! kSmall)
			{
				if (!fNodalShapes->SetFieldAt(nodal_init_coords, NULL)) 
					ExceptionT::GeneralFail(caller,"Shape Function evaluation"
						"failed at node %d\n",i);
						
				const dArrayT& phiValues = fNodalShapes->FieldAt();	
				
				phis.Dimension(phiValues.Length());
				phis = phiValues;	
				phis /= R_i;
			}
			else
			{
				if (!fNodalShapes->SetDerivativesAt(nodal_init_coords))
					ExceptionT::GeneralFail(caller,"Shape Function derivate evaluation"
						"failed at node %d\n",i);
				
				const dArray2DT& DphiValues = fNodalShapes->DFieldAt();

				phis.Dimension(DphiValues.MajorDim());
				//Copy the first column of DphiValues, i.e. d phi / d R = lim_{R -> 0} phi/R
				phis.Copy(DphiValues.Pointer());
			} 		
			
			
			iArrayT node_cover(fNodalShapes->Neighbors());	
			int n_node_cover = node_cover.Length();	
			iArrayT node_cover_key(node_cover);
			node_cover_key.SetValueToPosition();
			node_cover_key.SortAscending(node_cover);

			LinkedListT<int>& supp_0 = nodeWorkSpace[i];
			LinkedListT< double >& circumf_0 = circumferentialWorkSpace[i];
			int s_0 = -1;
			/* Simultaneously loop over support of the two nodes that are endpoints of the
			 * current Delone edge and the nodes in the support of the midpoint of this
			 * edge. If a node covering the centroid is not in the support of n_0 or n_1,
			 * insert that covering node into the sorted list.
			 */
			 
			int* n = node_cover.Pointer();
			int* n_j = node_cover_key.Pointer();
			
			supp_0.Top(); circumf_0.Top();
			next_0 = supp_0.CurrentValue();
			for (int j = 0; j < n_node_cover; j++, n++, n_j++)
			{			
				if (next_0)
					traverseQ_0 = *next_0 <= *n;
				else
					ExceptionT::GeneralFail(caller,"Support list does not exist\n");
						
				// advance supp_0 until it is equal to current node
				// What we're really doing is skipping nodes in support of facet 
				// centroids (our smoothed strain integration points) that are 
				// not in support of the node
				while (traverseQ_0 && supp_0.Next(s_0) && circumf_0.Next())
				{
					next_0 = supp_0.PeekAhead(); 
					if (!next_0)
						traverseQ_0 = false;
					else
						if (*next_0 > *n)
							traverseQ_0 = false;
				}
					
				if (s_0 != *n) 
					ExceptionT::GeneralFail(caller,"Node %d in support of node %d but not in data\n",s_0,*n);
				if (supp_0.AtTop())
					circumf_0.Next();
				
				*(circumf_0.CurrentValue()) = phis[*n_j];
			}
		}
	}
	
	// move into more efficient storage for computation
	cellSupports.Configure(nodeWorkSpace);
	bVectors.Configure(facetWorkSpace);
	if (qIsAxisymmetric) 
		circumferential_B.Configure(circumferentialWorkSpace);
	
	for (int i = 0; i < cellSupports.MajorDim(); i++) {
		int* irow_i = cellSupports(i);
		dArrayT* drow_i = bVectors(i);
		LinkedListT<int>& ilist = nodeWorkSpace[i];
		LinkedListT<dArrayT>& dlist = facetWorkSpace[i];
		LinkedListT<double>* clist;
		double* crow_i;
		if (qIsAxisymmetric) {
			clist = &circumferentialWorkSpace[i];
			clist->Top();
			crow_i = circumferential_B(i);
		}
		ilist.Top(); dlist.Top();
		while (ilist.Next() && dlist.Next()) {
			*irow_i++ = *(ilist.CurrentValue());
			*drow_i++ = *(dlist.CurrentValue());
			if (qIsAxisymmetric) {
				clist->Next();
				*crow_i++ = *(clist->CurrentValue());
			}
		}
	}
}

void CellFromMeshT::BoundaryShapeFunctions(RaggedArray2DT<double>& phis, RaggedArray2DT<int>& supports, dArray2DT& normals)
{
#pragma unused(phis)
#pragma unused(supports)
#pragma unused(normals)
	// for traction BCs, these data structures are needed
	// phis are shape function values for nodes covering integration points on boundary facets
	// supports are the locally-numbered indices
	// normals are the facet normal vectors 
}

void CellFromMeshT::DefineElements(const ArrayT<StringT>& block_ID, const ArrayT<int>& mat_index) 
{
	const char caller[] = "CellFromMeshT::DefineElements";

	/* inherited */
	CellGeometryT::DefineElements(block_ID, mat_index);
	
	/* store block ID's */
	fBlockID = block_ID;
}
