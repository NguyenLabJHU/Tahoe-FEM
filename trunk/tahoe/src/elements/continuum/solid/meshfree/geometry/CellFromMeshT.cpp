/* $Id: CellFromMeshT.cpp,v 1.6 2005-01-28 02:41:41 paklein Exp $ */
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
	int nNodes = fNodalCoordinates->MajorDim();
	
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

	dArrayT ip_coords(nSD);
	dMatrixT jacobian(nSD, 1);

	dMatrixT Q(nSD);
	iArrayT nodes_glb(nElementNodes), nodes_loc(nElementNodes);
	dArrayT facetNormal(nSD), facetIntegral(nSD);
	int n_0;
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

				/* loop over facets (internal and on the element boundary) for each node */
				for (int jj = 0; jj < n_faces; jj++) {
	
					/* collect coordinates over the facet */
					sub_cell_shape.NodesOnFacet(jj, facet_nodes);
					facet_coords.Collect(facet_nodes, sub_cell_coords);

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
						const iArrayT& neighbors = fNodalShapes->Neighbors();	

						MergeFacetIntegral(n_0, jw, facetNormal, phiValues, neighbors);
					}
				}		
			}
		}
	}
	
	/* compute cell centroids */
	for (int e = 0; e < cellVolumes.Length(); e++)
		cellCentroids.ScaleRow(e, 1.0/cellVolumes[e]);
	
	ConfigureDataStructures(cellSupports, bVectors, circumferential_B, cellVolumes);
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
