/* $Id */
#include "VoronoiDiagramT.h"
 
#ifdef __QHULL__
#include "CompGeomT.h"
#endif

#include "ArrayT.h"
#include "ofstreamT.h"
#include "ifstreamT.h"
#include "fstreamT.h"
#include "dArrayT.h"
#include "LinkedListT.h"
#include "LocalArrayT.h"
#include "ParentDomainT.h"
#include "ModelManagerT.h"

#include "toolboxConstants.h"

/* uncomment this line for many, many boundary integration points;
 * --set the number of points in the routine computeBmatrices
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

using namespace Tahoe;

/* constructors */
VoronoiDiagramT::VoronoiDiagramT(const ElementSupportT& support, bool isAxisymmetric):
  CellGeometryT(support, isAxisymmetric),
  fVoronoi(NULL),
  qComputeVoronoiCell(false),
  qJustVoronoiDiagram(false),
  vCellFile("voronoidiagram")
{

  SetName("voronoi_diagram");

}

VoronoiDiagramT::VoronoiDiagramT():
  CellGeometryT(),
  fVoronoi(NULL),
  qComputeVoronoiCell(false),
  qJustVoronoiDiagram(false),
  vCellFile("voronoidiagram")
{

  SetName("voronoi_diagram");

}

/* destructor */
VoronoiDiagramT::~VoronoiDiagramT(void)
{
#ifdef __QHULL__		
  if (fVoronoi) 
    delete fVoronoi;
#endif
}

void VoronoiDiagramT::DefineElements(const ArrayT<StringT>& block_ID, const ArrayT<int>& mat_index) 
{

  const char caller[] = "VoronoiDiagramT::DefineElements";
	
  CellGeometryT::DefineElements(block_ID, mat_index);

  /* access to the model database */
  ModelManagerT& model = fElementSupport->ModelManager();
	
  /* get nodes and facets on the boundary */
  GeometryBaseT::CodeT facetType;
  iArrayT facet_numbers, element_numbers;
  model.SurfaceFacets(block_ID, facetType, fBoundaryConnectivity, fBoundaryNodes, 
		      facet_numbers, element_numbers, NULL);
  fBoundaryIsTriangulated = (facetType == GeometryT::kLine) ||
    (facetType == GeometryT::kTriangle);
  fBoundaryNodes.SortAscending();
	
  /* write boundary nodes to output */
  if (fElementSupport->PrintInput()) {
    ofstreamT& out = fElementSupport->Output();
    fBoundaryNodes++;
    out << "\n " << caller << ": boundary nodes\n" << fBoundaryNodes.wrap(10) << endl;
    fBoundaryNodes--;
  }

  /* convert to local numbering for qhull */
  fscnimft->GlobalToLocalNumbering(fBoundaryNodes);
  iArrayT trickArray(fBoundaryConnectivity.Length(), fBoundaryConnectivity.Pointer());
  fscnimft->GlobalToLocalNumbering(trickArray);
	
  /* don't need this information */
  facet_numbers.Free();
  element_numbers.Free();
}

void VoronoiDiagramT::ComputeBMatrices(RaggedArray2DT<int>& cellSupports, RaggedArray2DT<dArrayT>& bVectors,
				       dArrayT& cellVolumes, dArray2DT& cellCentroids, RaggedArray2DT<double>& circumferential_B)
{

  /* Here for the Axisymmetric case, also computes {Psi_I(X_L)/R)L,0.} 
   */

  const char caller[] = "VoronoiDiagramT::ComputeBMatrices";

  // Do the heavy lifting for the Voronoi Diagram now
  if (qComputeVoronoiCell) {
#ifndef __QHULL__
    ExceptionT::GeneralFail(caller,"Requires the QHull library\n");
#else 

    fVoronoi = new CompGeomT(*fNodalCoordinates);
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
	
  // centroid information is already here
  cellCentroids.Alias(fVoronoiCellCentroids);

  /* possible best implementation is to loop over all Delone edges
   * and compute all the necessary values only once per Voronoi
   * facet. This approach minimizes number of times that the support of
   * an arbitrary point in space (the Voronoi facet centroid) has to be
   * found.
   */
  const RaggedArray2DT<int>& nodeSupport = fNodalShapes->NodeNeighbors();
  cellVolumes.Alias(fVoronoiCellVolumes);
	
  int nNodes = fNodalCoordinates->MajorDim();
	
  nodeWorkSpace.Dimension(nNodes);
  facetWorkSpace.Dimension(nNodes);
  if (qIsAxisymmetric) {
    circumferentialWorkSpace.Dimension(nNodes);
  }
	
  int nsd = fNodalCoordinates->MinorDim();
  dArrayT zeroFacet(nsd);
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

  /* integration */
  int nfn = 2;
  ParentDomainT domain(GeometryT::kLine, fNumIP, nfn);
  domain.Initialize();
  LocalArrayT facet_coords(LocalArrayT::kInitCoords, nfn, nsd);
  facet_coords.SetGlobal(fVoronoiVertices);
  iArrayT keys;
  dArrayT ip_coords(nsd), phis;
  dMatrixT jacobian(nsd, 1);
  const double* ip_weight = domain.Weight();

  dArrayT facetNormal(nsd), facetIntegral(nsd);
  int n_0, n_1;
  for (int i = 0; i < fDeloneEdges.MajorDim(); i++) {
    n_0 = fDeloneEdges(i,0);
    n_1 = fDeloneEdges(i,1); 
		
    facetNormal.DiffOf((*fNodalCoordinates)(n_1), (*fNodalCoordinates)(n_0));
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
	ExceptionT::GeneralFail("VoronoiDiagramT::ComputeBMatrices","Shape Function evaluation"
				"failed at Delone edge %d\n",i);
				
      const dArrayT& phiValues = fNodalShapes->FieldAt();	
      const iArrayT& neighbors = fNodalShapes->Neighbors();		
			
      MergeFacetIntegral(n_0, jw, facetNormal, phiValues, neighbors);
			
      // this only needs to be done for loops over dual edges
      facetNormal *= -1.;
			
      MergeFacetIntegral(n_1, jw, facetNormal, phiValues, neighbors);
    }
  }

  /** temporary storage for integration over the body boundary */
  boundary_phi.Dimension(fNonDeloneEdges.Length());
  boundary_supports.Dimension(fNonDeloneEdges.Length());
	
  double zero = 0.0;
  dArrayT boundaryIPCoord(nsd);
  for (int i = 0; i < fNonDeloneEdges.Length(); i++) {
    double* v1 = fVoronoiVertices(fSelfDualFacets(i,0));
    double* v2 = fVoronoiVertices(fSelfDualFacets(i,1));
    for (int j = 0; j < nsd; j++)
      boundaryIPCoord[j] = v1[j] + v2[j];
    boundaryIPCoord /= double(nsd);
		
    if (!fNodalShapes->SetFieldAt(boundaryIPCoord, NULL)) // shift = 0 or not ?
      ExceptionT::GeneralFail("VoronoiDiagramT::ComputeBMatrices","Shape Function evaluation"
			      "failed at Delone edge %d\n",i);
					
    const dArrayT& phiValues = fNodalShapes->FieldAt();			
					
    iArrayT supp_i(fNodalShapes->Neighbors());	
    int l_supp_i = supp_i.Length();
    supp_i.SortAscending();
    boundary_supports[i].AppendArray(l_supp_i, supp_i.Pointer());
    boundary_phi[i].AppendArray(l_supp_i, zero);
  }
 
  /** Loop over remaining edges */
  for (int i = 0; i < fNonDeloneEdges.Length(); i++) {
    n_0 = fNonDeloneEdges[i];
    facetNormal.Set(nsd, fNonDeloneNormals(i));
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
    dArrayT ip_coord0(nsd,v1);
    dArrayT edgeVector(nsd);
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
	ip_coords.Set(nsd,fNodalCoordinates(n_0));
#endif
#endif // SEPARATE_BOUNDARY_INTEGRATION				

	if (!fNodalShapes->SetFieldAt(ip_coords, NULL)) // shift = 0 or not ?
	  ExceptionT::GeneralFail("VoronoiDiagramT::ComputeBMatrices","Shape Function evaluation"
				  "failed at Delone edge %d\n",i);
					
	const dArrayT& phiValues = fNodalShapes->FieldAt();
	phis.Dimension(phiValues.Length());
	phis = phiValues;
	phis *= jw;
					
	const iArrayT& neighbors = fNodalShapes->Neighbors();	
			
	MergeFacetIntegral(n_0, jw, facetNormal, phiValues, neighbors);
						
	/* Merge support of the boundary node with covering of integration point
	 */
	MergeNodalValues(i, phis, neighbors, boundary_supports, 
			 boundary_phi, true);
      }	
    }
	
    ConfigureDataStructures(cellSupports, bVectors, circumferential_B, cellVolumes);
	
  }


void VoronoiDiagramT::ComputeBprimeMatricesSS(RaggedArray2DT<dMatrixT>& bprimeVectors, const RaggedArray2DT<int>& cellSupports,
					  const RaggedArray2DT<dArrayT>& bVectors, const dArrayT& cellVolumes,
					  const dArray2DT& cellCentroids, dArray2DT& Ymatrices)
{
  ExceptionT::GeneralFail("VoronoiDiagramT::ComputeBprimeMatricesSS","Currently not available for Voronoi diagram - use cell from mesh");
}

void VoronoiDiagramT::VoronoiDiagramToFile(ofstreamT& vout)
{
	
  int nNodes = fNodalCoordinates->MajorDim();
  int nVertices = fVoronoiVertices.MajorDim();
  int nSD = fVoronoiVertices.MinorDim();

  vout << nNodes << "\n";
  vout << nSD << "\n";
  vout << nVertices << "\n";

  // write out vertices of the VoronoiDiagram
  for (int i = 0; i < nVertices; i++) {
    for (int j = 0; j < nSD; j++)
      vout << fVoronoiVertices(i,j) << " "; 
    vout << "\n";
  }
    
  // write out Voronoi cells for each node
  for (int i = 0; i < nNodes; i++) {
    vout << i <<"\n";

    for (int j = 0; j < nSD; j++)
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
    ExceptionT::GeneralFail("VornoiDiagramT::VoronoiDiagramToFile","Dual edge/facet dimension mismatch\n");

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
    for (int k = 0; k < nSD; k++)
      vout << fSelfDualFacets(i,k) << " "; // SPECIALIZED TO 2D!!!
    vout << "\n";
  }

  // list of normals of self-dual facets
  for (int i = 0; i < numSelfDuals; i++) {
    for (int k = 0; k < nSD; k++) {
      vout << fNonDeloneNormals(i,k) << " ";
    }
    vout << "\n";
  }

  // list of areas of self-dual facets
  for (int i = 0; i < numSelfDuals; i++)
    vout << fBoundaryIntegrationWeights[i] << "\n";
 
}	
	
void VoronoiDiagramT::VoronoiDiagramFromFile(ifstreamT& vin)
{
  const char caller[] = "VoronoiDiagramT::VoronoiDiagramFromFile";	

  int nNodes, nSD, nVertices;
  vin >> nNodes;
    
  /* minor consistency checks */
  if (nNodes != fNodalCoordinates->MajorDim()) {
    vin.close();
    ExceptionT::GeneralFail(caller,"Input Voronoi file does not match node number\n");
  }
    
  vin >> nSD;
  if (nSD != fNodalCoordinates->MinorDim()) {
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

    for (int j = 0; j < nSD; j++)
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
  fNonDeloneNormals.Dimension(nCentroids,nSD);
  fBoundaryIntegrationWeights.Dimension(nCentroids);

  fSelfDualFacets.Dimension(nCentroids, 2);
	
  for (int i = 0; i < nCentroids; i++)
    vin >> fNonDeloneEdges[i];	

  for (int i = 0; i < nCentroids; i++) {
    vin >> fSelfDualFacets(i,0) >> fSelfDualFacets(i,1);
  }
			
  for (int i = 0; i < nCentroids; i++)
    for (int j = 0; j < nSD; j++)
      vin >> fNonDeloneNormals(i,j);
			
  for (int i = 0; i < nCentroids; i++)
    vin >> fBoundaryIntegrationWeights[i];
}

void VoronoiDiagramT::InitializeVoronoiData(void)
{
#ifdef __QHULL__
  // use qhull's data structures but make our own specialized versions
  CompGeomT::ConvexHullMap selfDuals;
  CompGeomT::VoronoiDiagramMap voronoiFacetIndices;
  ArrayT<dArrayT> voronoiFacetAreas;
  ArrayT<dArray2DT> voronoiFacetNormals;
  CompGeomT::ConvexHullMap voronoiCells;

  fVoronoiVertices.Alias(fVoronoi->VoronoiVertices());
  int nsd = fVoronoiVertices.MinorDim();
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
  fNonDeloneNormals.Dimension(numSelfDuals, nsd);
  fBoundaryIntegrationWeights.Dimension(numSelfDuals);
  fSelfDualFacets.Dimension(numSelfDuals, 2);
  
  int ctr = 0;
  double *v1;
  
  // list of centroids of self-dual facets
  iArrayT* thisFacet;
  int thisFacetLength;
  dArrayT ptArray(nsd); // workspace for centroids
  double *pt = ptArray.Pointer(); 
  
  for (int i = 0; i < selfDuals.Length(); i++)
    for (int j = 0; j < selfDuals[i].Length(); j++) {
      fNonDeloneEdges[ctr] = fBoundaryNodes[i];
      thisFacet = &voronoiFacetIndices[fBoundaryNodes[i]][selfDuals[i][j]];
      thisFacetLength = thisFacet->Length();
      ptArray = 0.;
      for (int k = 0; k < thisFacetLength; k++) {
	v1 = fVoronoiVertices(voronoiCells[fBoundaryNodes[i]][(*thisFacet)[k]]);
	for (int l = 0; l < nsd; l++)
	  pt[l] += v1[l];
	fSelfDualFacets(ctr,k) =  voronoiCells[fBoundaryNodes[i]][(*thisFacet)[k]]; 
      }
      
      for (int k = 0; k < nsd; k++) 
	fNonDeloneNormals(ctr,k) =  voronoiFacetNormals[fBoundaryNodes[i]](selfDuals[i][j],k);
      
      fBoundaryIntegrationWeights[ctr] =  voronoiFacetAreas[fBoundaryNodes[i]][selfDuals[i][j]];
      ctr++;
    }
    	
  // Data for Axisymmetric mass matrix calculation
  int nVoronoiCells = voronoiCells.Length();
  fVoronoiCellCentroids.Dimension(nVoronoiCells, nsd);
  fVoronoiCellCentroids = 0.;
  for (int i = 0; i < nVoronoiCells; i++) {
    iArrayT& cell_i = voronoiCells[i];
    for (int j = 0; j < cell_i.Length(); j++) 
      fVoronoiCellCentroids.AddToRowScaled(i,1.,fVoronoiVertices(cell_i[j]));
    fVoronoiCellCentroids.ScaleRow(i,1./double(cell_i.Length()));
  }
#endif
}

void VoronoiDiagramT::BoundaryShapeFunctions(RaggedArray2DT<double>& phis, RaggedArray2DT<int>& supports, dArray2DT& normals)
{
  normals.Alias(fNonDeloneNormals);
  phis.Configure(boundary_phi);
  supports.Configure(boundary_supports);
	
  for (int i = 0; i < boundary_supports.Length(); i++) {
    int* irow_i = supports(i);
    double* drow_i = phis(i);
    LinkedListT<int>& ilist = boundary_supports[i];
    LinkedListT<double>& dlist = boundary_phi[i];
    ilist.Top(); dlist.Top();
    while (ilist.Next() && dlist.Next()) {
      *irow_i++ = *(ilist.CurrentValue());
      *drow_i++ = *(dlist.CurrentValue());
    }
  }
}

// xml stuff

/* initialization */
void VoronoiDiagramT::TakeParameterList(const ParameterListT& list)
{
  const char caller[] = "VoronoiDiagramT::TakeParameterList";
	
  CellGeometryT::TakeParameterList(list);

  /* resolve path to Voronoi file */
  StringT path;
  path.FilePath(fElementSupport->InputFile());
  vCellFile = list.GetParameter("voronoi_file");
  vCellFile.ToNativePathName();
  vCellFile.Prepend(path);

  qComputeVoronoiCell = list.GetParameter("compute_voronoi");
  qJustVoronoiDiagram = list.GetParameter("just_voronoi_diagram");

  if (qJustVoronoiDiagram) // override input error if we just nee the geometry
    qComputeVoronoiCell = true;

}

/* describe the parameters needed by the interface */
void VoronoiDiagramT::DefineParameters(ParameterListT& list) const
{
  /* inherited */
  CellGeometryT::DefineParameters(list);

  ParameterT compute_voronoi(qComputeVoronoiCell, "compute_voronoi");
  compute_voronoi.SetDefault(qComputeVoronoiCell);
  list.AddParameter(compute_voronoi);
	
  ParameterT voronoi_file(vCellFile, "voronoi_file");
  voronoi_file.SetDefault(vCellFile);
  list.AddParameter(voronoi_file);

  ParameterT just_voronoi_diagram(qJustVoronoiDiagram,"just_voronoi_diagram");
  just_voronoi_diagram.SetDefault(qJustVoronoiDiagram);
  list.AddParameter(just_voronoi_diagram);
}


/*void VoronoiDiagramT::FindCornersAndEdges(iArrayT& boundaryFacets, iArrayT& boundaryElements) {
  int offs = boundaryElements.Min();
  iArrayT elementCount(boundaryElements.Max()-offs+1);
  iArrayT keys(elementCount.Length());
  elementCount = 0;

  for (int i = 0; i < boundaryElements.Length(); i++)
  elementCount[offs + boundaryElements[i]]++;
 
  //keys.SetValueToPosition();
  //elementCount.SortAscending(keys);

  // BIG assumptions here
  bodyCorners.Dimension(elementCount.Count(3));
  bodyEdges.Dimension(elementCount.Count(2),3);

  // would like bodyCorners to have the node index
  // would like bodyEdges to have the edge index
  */
/*int i, i0;
  for (i = 0; i < elementCount.Length() && elementCount[keys[i]] < 2; i++)
  ; // skip interior elements and facets on the boundary
  i0 = i;
  for (; i < elementCount.Length() && elementCount[keys[i]] < 3; i++) 
  bodyEdges(i-i0,0) = keys[i] + offs; //element keys[i] + offs has a boundary edge
  i0 = i;
  for (; i < elementCount.Length(); i++) 
  bodyCorners[i - i0] = keys[i] + offs; // element keys[i] + offs has a corner node
  
  // would like to quickly know if an element has an edge or a corner (or both)
  InverseMapT edgeInverse, cornerInverse;
  edgeInverse.SetMap(bodyEdges);
  edgeInverse.SetOutOfRange(InverseMapT::MinusOne);
  cornerInverse.SetMap(bodyCorners);
  cornerInverse.SetOutOfRange(InverseMapT::MinusOne);

  for (int i = 0; i < boundaryElements.Length(); i++) {
  if (edgeInverse.Map(boundaryElements[i]) != -1) {
  ; // find the boundary edge
  }
  if (cornerInverse.Map(boundaryElements[i]) != -1) {
  ; // find the corner node
  }
  }      */

//}

