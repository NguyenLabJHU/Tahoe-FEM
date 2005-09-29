/* $Id: CellFromMeshT.cpp,v 1.9 2005-09-29 19:11:39 jcmach Exp $ */
#include "CellFromMeshT.h"

#include "ArrayT.h"
#include "dArrayT.h"
#include "dSymMatrixT.h"
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
  dMatrixT jacobian(nSD, nSD-1); /* jacoban of surface mapping */

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
	if (!fscnimft->GlobalToLocalNumbering(nodes_loc))
	  ExceptionT::GeneralFail(caller, "list contains point that is not meshfree");

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

	      // n_0 (int): node number of node associated w/ subdomain k
	      // jw (double): jacobian-determinant * integration weight of point ii
	      // facetNormal (dArrayT(nsd)): normal vector at the facet
	      // phiValues (dArrayT(ip neighbors)): shape fcns at point ii
	      // neighbors (iArrayT(ip neighbors)): neighbors of point ii (order
	      //                                    corresponds to phiValues

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
//******************************************************************************************************************************************
void CellFromMeshT::ComputeBprimeMatricesSS(RaggedArray2DT<dMatrixT>& bprimeVectors, const RaggedArray2DT<int>& cellSupports,
					  const RaggedArray2DT<dArrayT>& bVectors, const dArrayT& cellVolumes,
					  const dArray2DT& cellCentroids, dArray2DT& Ymatrices)
{

  const char caller[] = "CellFromMeshT::ComputeBprimeMatricesSS";

  int nSD = fElementSupport->NumSD();
  int nNodes = fNodalCoordinates->MajorDim();
	
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
  ShapeFunctionT sub_cell_shape(sub_cell_geometry, 4, sub_cell_coords);
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

  dArrayT ip_coords(nSD), nd_coords(nSD);
  dMatrixT jacobian(nSD, nSD-1); /* jacoban of surface mapping */

  dMatrixT Q(nSD);
  iArrayT nodes_glb(nElementNodes), nodes_loc(nElementNodes);
  dArrayT facetNormal(nSD), facetIntegral(nSD);
  int n_0;

  /* Data structures used for calculation of bprime matrices */
  int size_sym = (nSD == 2 ? 3 : 6);
  Ymatrices.Dimension(nNodes,size_sym);
  Ymatrices = 0.;

  dArrayT Ymatrix(size_sym);
  Ymatrix = 0.;
  dSymMatrixT YmatrixInv(dSymMatrixT::int2DimensionT(nSD));
  YmatrixInv = 0.;

  dArrayT zvector(nSD);
  zvector = 0.;

  dArrayT yL(nSD), c_coords(nSD), DphiL(nSD);
  yL = 0.;
  c_coords = 0.;
  DphiL = 0.;

  dMatrixT DiagCmatrix(nSD), Dmatrix(nSD);
  DiagCmatrix = 0.;
  Dmatrix = 0.;

  RaggedArray2DT<dMatrixT> Cmatrices;
  iArrayT ncellSupports(nNodes), AllNodesGlb(nNodes), AllNodesLoc(nNodes);
  ncellSupports = 0;
  AllNodesGlb = 0;
  AllNodesLoc = 0;
  for (int i = 0; i < nNodes; i++)  
    AllNodesGlb[i] = (*fNodes)[i];
  AllNodesLoc = AllNodesGlb;
  if (!fscnimft->GlobalToLocalNumbering(AllNodesLoc))
    ExceptionT::GeneralFail(caller, "list contains point that is not meshfree");
  for (int i = 0; i < nNodes; i++) {
    int node_i = AllNodesLoc[i];
    ncellSupports[node_i] = cellSupports.MinorDim(node_i);
  }

  Cmatrices.Configure(ncellSupports);
  bprimeVectors.Configure(ncellSupports);

  dMatrixT zeroMatrix(nSD);
  zeroMatrix = 0.;
  Cmatrices = zeroMatrix;
  bprimeVectors = zeroMatrix;

  for (int e = 0; e < connectivities.Length(); e++) /* loop over element blocks */
    {
      /* block connectivities */
      const iArray2DT& connects = *(connectivities[e]);
      int nElements = connects.MajorDim();
	
      for (int i = 0; i < nElements; i++) { /* loop over elements in the block */

	/* element nodes */
	connects.RowAlias(i, nodes_glb);
	nodes_loc = nodes_glb;
	if (!fscnimft->GlobalToLocalNumbering(nodes_loc))
	  ExceptionT::GeneralFail(caller, "list contains point that is not meshfree");

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

	  /* compute Y matrices */
	  sub_cell_shape.TopIP();
	  while (sub_cell_shape.NextIP()) 
	    {
	      int ip = sub_cell_shape.CurrIP();

	      /* JacobianDet*IPweight */
	      double dv = sub_cell_det[ip]*sub_cell_wgt[ip];
					
	      /* IP coordinates */
	      sub_cell_shape.IPCoords(ip_coords);

	      /* Integrate Y matrices */
	      for (int kk = 0; kk < size_sym; kk++) {
		if (nSD == 2) {
		  switch (kk) {
		  case 0: 
		    Ymatrices(n_0,0) += (ip_coords[0] - cellCentroids(n_0,0))*(ip_coords[0] - cellCentroids(n_0,0))*dv;
		    break;
		  case 1:
		    Ymatrices(n_0,1) += (ip_coords[1] - cellCentroids(n_0,1))*(ip_coords[1] - cellCentroids(n_0,1))*dv;
		    break;
		  case 2:
		    Ymatrices(n_0,2) += (ip_coords[0] - cellCentroids(n_0,0))*(ip_coords[1] - cellCentroids(n_0,1))*dv;
		    break;
		  default:
		    ExceptionT::GeneralFail("SCNIMFT::ComputeBprimeMatricesSS","Failure integrating Ymatrices\n");
		    break;
		  }
		}
		else if (nSD == 3) {  
		  switch (kk) {
		  case 0: 
		    Ymatrices(n_0,0) += (ip_coords[0] - cellCentroids(n_0,0))*(ip_coords[0] - cellCentroids(n_0,0))*dv;
		    break;
		  case 1:
		    Ymatrices(n_0,1) += (ip_coords[1] - cellCentroids(n_0,1))*(ip_coords[1] - cellCentroids(n_0,1))*dv;
		    break;
		  case 2:
		    Ymatrices(n_0,2) += (ip_coords[2] - cellCentroids(n_0,2))*(ip_coords[2] - cellCentroids(n_0,2))*dv;
		    break;
		  case 3:
		    Ymatrices(n_0,3) += (ip_coords[1] - cellCentroids(n_0,1))*(ip_coords[2] - cellCentroids(n_0,2))*dv;
		    break;
		  case 4:
		    Ymatrices(n_0,4) += (ip_coords[0] - cellCentroids(n_0,0))*(ip_coords[2] - cellCentroids(n_0,2))*dv;
		    break;
		  case 5:
		    Ymatrices(n_0,5) += (ip_coords[0] - cellCentroids(n_0,0))*(ip_coords[1] - cellCentroids(n_0,1))*dv;
		    break;
		  default:
		    ExceptionT::GeneralFail("SCNIMFT::ComputeBprimeMatricesSS","Failure integrating Ymatrices\n");
		    break;
		  }
		}
		else
		  ExceptionT::GeneralFail("SCNIMFT::ComputeBprimeMatricesSS","Failure integrating Ymatrices\n)");
	      }	    
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

	      /* Map IP neighbors according to cell neighbors - figure out where IP neighbors
		 are in the current row of cellSupports */
	      int ncn = cellSupports.MinorDim(n_0);
	      iArrayT nSupport(ncn);
	      nSupport = 0;
	      cellSupports.RowAlias(n_0,nSupport);
	      int nIPn = neighbors.Length();
	      const int* IPn = neighbors.Pointer();
	      iArrayT mapn(nIPn);
	      mapn = 0;
	      int* pmapn = mapn.Pointer();
	      int testn = 0;
	      for (int nn = 0; nn < nIPn; nn++, IPn++, pmapn++) { 
		testn = nSupport.HasValue(*IPn,*pmapn);
		if (testn == 1) 
		  testn = 0;
		else {
		  ExceptionT::GeneralFail("SCNIMFT::ComputeBprimeMatricesSS",
					  "IP neighbor not found in cell neighbor list");
		  break;
		}
	      }

	      /* Loop over IP neighbors */
	      pmapn = mapn.Pointer();
	      for (int kk = 0; kk < nIPn; kk++, pmapn++) {
		dMatrixT& cI = Cmatrices(n_0,*pmapn);
		/* Compute boundary integral part of C matrix */
		for (int ii = 0; ii < nSD*nSD; ii++) {
		  if (nSD == 2) {
		    switch (ii) {
		    case 0: 
		      cI(0,0) += phiValues[kk]*facetNormal[0]*ip_coords[0]*jw;
		      break;
		    case 1:
		      cI(1,0) += phiValues[kk]*facetNormal[1]*ip_coords[0]*jw;
		      break;
		    case 2:
		      cI(0,1) += phiValues[kk]*facetNormal[0]*ip_coords[1]*jw;
		      break;
		    case 3:
		      cI(1,1) += phiValues[kk]*facetNormal[1]*ip_coords[1]*jw;
		      break;
		    default:
		      ExceptionT::GeneralFail("SCNIMFT::ComputeBprimeMatricesSS","Failure integrating Cmatrices\n");
		      break;
		    }
		  }
		  else if (nSD == 3) {  
		    switch (ii) {
		    case 0: 
		      cI(0,0) += phiValues[kk]*facetNormal[0]*ip_coords[0]*jw;
		      break;
		    case 1:
		      cI(1,0) += phiValues[kk]*facetNormal[1]*ip_coords[0]*jw;
		      break;
		    case 2:
		      cI(2,0) += phiValues[kk]*facetNormal[2]*ip_coords[0]*jw;
		      break;
		    case 3:
		      cI(0,1) += phiValues[kk]*facetNormal[0]*ip_coords[1]*jw;
		      break;
		    case 4:
		      cI(1,1) += phiValues[kk]*facetNormal[1]*ip_coords[1]*jw;
		      break;
		    case 5:
		      cI(2,1) += phiValues[kk]*facetNormal[2]*ip_coords[1]*jw;
		      break;
		    case 6:
		      cI(0,2) += phiValues[kk]*facetNormal[0]*ip_coords[2]*jw;
		      break;
		    case 7:
		      cI(1,2) += phiValues[kk]*facetNormal[1]*ip_coords[2]*jw;
		      break;
		    case 8:
		      cI(2,2) += phiValues[kk]*facetNormal[2]*ip_coords[2]*jw;
		      break; 
		    default:
		      ExceptionT::GeneralFail("SCNIMFT::ComputeBprimeMatricesSS","Failure integrating Cmatrices\n");
		      break;
		    }
		  }
		  else
		    ExceptionT::GeneralFail("SCNIMFT::ComputeBprimeMatricesSS","Failure integrating Cmatrices\n)");
		}
	      }
	    }
	  }		
	}
      }
    }
  

  /*  Subtract the other term from the Cmatrices, compute Dmatrices and zvectors for each cell (node) L */
  for (int ii = 0; ii < nNodes; ii++) {

    /* Local meshfree node number */
    n_0 = AllNodesLoc[ii];
    
    /* Compute YMatrixInv */
    Ymatrices.RowCopy(n_0,Ymatrix);
    double* pYm = Ymatrix.Pointer();
    YmatrixInv.Alias(dSymMatrixT::int2DimensionT(nSD),pYm);
    YmatrixInv = YmatrixInv.Inverse();

    /* Nodal coordinates */
    fNodalCoordinates->RowAlias(n_0,nd_coords);
    
    if (!fNodalShapes->SetFieldAt(nd_coords, NULL)) 
      ExceptionT::GeneralFail("SCNIMFT::ComputeBprimeMatricesSS","Shape Function evaluation failed at node %i\n",n_0);

    /* Get shape functions, shape function derivatives, and neighbor list at node L */
    const dArrayT& phiValuesL = fNodalShapes->FieldAt();
    const dArray2DT& DphiValuesL = fNodalShapes->DFieldAt();
    const iArrayT& neighborsL = fNodalShapes->Neighbors();

    /* Calculate yL = xL - x(bar)L */
    cellCentroids.RowAlias(n_0,c_coords);
    yL = nd_coords;
    yL -= c_coords;

    /* Compute zvector */
    YmatrixInv.Multx(c_coords,zvector,cellVolumes[n_0]);

    /* Loop over cell neighbors */
    /* If a cell neighbor is a node neighbor then compute diagonal term of Cmatrix */
    int ncn = cellSupports.MinorDim(n_0);
    const int* Cn = cellSupports(n_0);
    for (int kk = 0; kk < ncn; kk++, Cn++) {
      dMatrixT& cI = Cmatrices(n_0,kk);
      int knode;
      if (neighborsL.HasValue(*Cn,knode)) {
	DphiValuesL.ColumnCopy(knode,DphiL);
	DiagCmatrix = 0.;
	DiagCmatrix.Identity(cellVolumes[n_0]*(phiValuesL[knode] - DphiL.Dot(DphiL,yL)));
	/* Subtract Diagonal term from C matrix */
	cI -= DiagCmatrix;
      }
      /* Compute Dmatrix */
      cI.Transpose();
      Dmatrix.MultSymAB(YmatrixInv,cI);
      Dmatrix.Transpose();
      /* Compute bprimevector (really a matrix) */
      dMatrixT& bprimeI = bprimeVectors(n_0,kk);
      const dArrayT& bI = bVectors(n_0,kk);
      bprimeI.Outer(bI,zvector,-1.0);
      bprimeI += Dmatrix;
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

