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

		fVoronoi = new CompGeomT(fNodalCoordinates);
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
	
	int nNodes = fNodalCoordinates.MajorDim();
	
	ArrayT< LinkedListT<int> > nodeWorkSpace; 
	ArrayT< LinkedListT<dArrayT> > facetWorkSpace; 
	ArrayT< LinkedListT<double> > circumferentialWorkSpace;
	nodeWorkSpace.Dimension(nNodes);
	facetWorkSpace.Dimension(nNodes);
	if (qIsAxisymmetric) {
		circumferentialWorkSpace.Dimension(nNodes);
	}
	
	dArrayT zeroFacet(3);
	double zeroSingle = 0.;
	zeroFacet = 0.0;
	for (int i = 0; i < nNodes; i++) {
		int l_supp_i = nodeSupport.MinorDim(i);
		iArrayT supp_i(l_supp_i);
		supp_i.Copy(nodeSupport(i));
		supp_i.SortAscending();
		nodeWorkSpace[i].AppendArray(l_supp_i, supp_i.Pointer());
		facetWorkSpace[i].AppendArray(l_supp_i, zeroFacet);
		if (qIsAxisymmetric)
			circumferentialWorkSpace[i].AppendArray(l_supp_i, zeroSingle);
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

	dArrayT facetNormal(nsd), facetIntegral(nsd);
	double* currentB, *currentI;
	int n_0, n_1;
	bool traverseQ_0, traverseQ_1;
	int *next_0, *next_1;
	for (int i = 0; i < fDeloneEdges.MajorDim(); i++) {
		n_0 = fDeloneEdges(i,0);
		n_1 = fDeloneEdges(i,1); 
		
		facetNormal.DiffOf(fNodalCoordinates(n_1), fNodalCoordinates(n_0));
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
		
			iArrayT ip_cover(fNodalShapes->Neighbors());	
			int n_ip_cover = ip_cover.Length();	
			iArrayT ip_cover_key(n_ip_cover);
			ip_cover_key.SetValueToPosition();
			ip_cover_key.SortAscending(ip_cover);

			LinkedListT<int>& supp_0 = nodeWorkSpace[n_0];
			LinkedListT<int>& supp_1 = nodeWorkSpace[n_1];
			LinkedListT< dArrayT >& bVectors_0 = facetWorkSpace[n_0];
			LinkedListT< dArrayT >& bVectors_1 = facetWorkSpace[n_1];
			LinkedListT< double > *circumf_0, *circumf_1;
			if (qIsAxisymmetric) {
				circumf_0 = &circumferentialWorkSpace[n_0];
				circumf_1 = &circumferentialWorkSpace[n_1];
			}	
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
			if (qIsAxisymmetric) {
				circumf_0->Top();
				circumf_1->Top();
			}
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
						if (qIsAxisymmetric)
							circumf_0->Next();
					}
				}
					
				currentI = facetIntegral.Pointer();
				currentB = bVectors_0.CurrentValue()->Pointer();
				for (int k = 0; k < nsd; k++)
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
					if (qIsAxisymmetric)
							circumf_1->InsertAtCurrent(0.);
					s_1 = *c;
					if (supp_1.AtTop()) { // if we're inserting at the front, LinkedListT's behavior requires more work
						supp_1.Next(); 
						bVectors_1.Next();
						if (qIsAxisymmetric)
							circumf_1->Next();
					}
				}
					 
				currentI = facetIntegral.Pointer();
				currentB =  bVectors_1.CurrentValue()->Pointer();
				for (int k = 0; k < nsd; k++)
					*currentB++ -= *currentI++; //NB change in sign; facet normal is inverted!
			}
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

				for (int k = 0; k < nsd; k++)
					*currentB++ += *currentI++;
			}
		
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

		}	
	}
 
 	// scale integrals by volumes of Voronoi cells
	dArrayT* currFacetIntegral;
	for (int i = 0; i < nNodes; i++) {
		LinkedListT<dArrayT>& bVectors_i = facetWorkSpace[i];
		LinkedListT<int>& nodes_i = nodeWorkSpace[i];
		bVectors_i.Top(); nodes_i.Top();
		while ((currFacetIntegral = bVectors_i.Next()))
			*currFacetIntegral *= 1./cellVolumes[i];
	}
	
	if (qIsAxisymmetric) {
		// calculate Psi/R terms. These are evaluated nodally, so this additional loop
		// is required
		dArrayT phis, nodal_init_coords;
		for (int i = 0; i < fNodalCoordinates.MajorDim(); i++) {
			nodal_init_coords.Set(nsd, fNodalCoordinates(i)); // This is the nodal coordinate.
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
				
				currentI = facetIntegral.Pointer();
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
			if (qIsAxisymmetric)
				*crow_i++ = *(clist->CurrentValue());
		}
	}
}

void VoronoiDiagramT::VoronoiDiagramToFile(ofstreamT& vout)
{
	
	int nNodes = fNodalCoordinates.MajorDim();
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
    if (nNodes != fNodalCoordinates.MajorDim()) {
		vin.close();
		ExceptionT::GeneralFail(caller,"Input Voronoi file does not match node number\n");
    }
    
    vin >> nSD;
    if (nSD != fNodalCoordinates.MinorDim()) {
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

