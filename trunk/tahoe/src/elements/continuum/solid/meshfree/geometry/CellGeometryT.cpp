/* $Id: CellGeometryT.cpp,v 1.7 2005-09-29 19:15:33 jcmach Exp $ */
#include "CellGeometryT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructors */
CellGeometryT::CellGeometryT(const ElementSupportT& support, bool isAxisymmetric):
	ParameterInterfaceT("cell_geometry"),
	fElementSupport(&support),
	fNodes(NULL),
	fNodalCoordinates(NULL),
	fNumIP(1),
	fscnimft(NULL),
	qIsAxisymmetric(isAxisymmetric)
{

}

/* constructors */
CellGeometryT::CellGeometryT(void):
	ParameterInterfaceT("cell_geometry"),
	fElementSupport(NULL),
	fNodes(NULL),
	fNodalCoordinates(NULL),
	fNumIP(1),
	fscnimft(NULL),
	qIsAxisymmetric(false)
{

}

void CellGeometryT::DefineElements(const ArrayT<StringT>& block_ID, const ArrayT<int>& mat_index) 
{
#pragma unused(block_ID)
#pragma unused(mat_index)
}

void CellGeometryT::SetNodesAndShapes(const iArrayT* nodes, const dArray2DT* nodal_coordinates, 
	MeshFreeNodalShapeFunctionT* nodalShapeFunctions)
{
	const char caller[] = "CellGeometryT::SetNodesAndShapes";

	fNodes = nodes;
	fNodalCoordinates = nodal_coordinates;
	if ((fNodes && !fNodalCoordinates) ||
	    (!fNodes && fNodalCoordinates))
	    ExceptionT::GeneralFail(caller, "nodes and coordinates must be set together");
	if (fNodes && fNodes->Length() != fNodalCoordinates->MajorDim())
		ExceptionT::SizeMismatch(caller);

	fNodalShapes = nodalShapeFunctions;
}

/* describe the parameters needed by the interface */
void CellGeometryT::DefineParameters(ParameterListT& list) const
{

	ParameterT num_ip(fNumIP, "num_ip");	
	num_ip.SetDefault(fNumIP);
	list.AddParameter(num_ip);

}

/* information about subordinate parameter lists */
void CellGeometryT::DefineSubs(SubListT& sub_list) const
{
#pragma unused(sub_list)
	
}

/* return the description of the given inline subordinate parameter list */
void CellGeometryT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
#pragma unused(name)
#pragma unused(order)
#pragma unused(sub_lists)
	
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* CellGeometryT::NewSub(const StringT& name) const
{
	return ParameterInterfaceT::NewSub(name);
}

/* initialization */
void CellGeometryT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "CellGeometryT::TakeParameterList";

	/* number of integration points used for surface integrals */
	fNumIP = list.GetParameter("num_ip");
}

void CellGeometryT::MergeFacetIntegral(int node_num, double weight, dArrayT& facetNormal, const dArrayT& phiValues,
						const iArrayT& neighbors) 
{
  /* If a node covering the integration point is not in the support of the node,
   * insert that covering node into the sorted list.
   */			 
  iArrayT cover(neighbors);	
  int n_cover = cover.Length();	
  iArrayT cover_key(n_cover);
  cover_key.SetValueToPosition();
  cover_key.SortAscending(cover);	

  LinkedListT<int>& supp = nodeWorkSpace[node_num];
  LinkedListT< dArrayT >& bVectors = facetWorkSpace[node_num];
  LinkedListT< double > *circumf;
  if (qIsAxisymmetric) 
    circumf = &circumferentialWorkSpace[node_num];
  int s = -1;
	
  supp.Top(); 
  bVectors.Top();
  int *next = supp.CurrentValue();
  if (qIsAxisymmetric) 
    circumf->Top();
		
  bool traverseQ;
  int* c = cover.Pointer();
  int* c_j = cover_key.Pointer();
  int nsd = facetNormal.Length();
  dArrayT facetIntegral(nsd), zeroFacet(nsd);
  zeroFacet = 0.;
  for (int j = 0; j < n_cover; j++, c++, c_j++) {
		
    facetIntegral = facetNormal;
    facetIntegral *= phiValues[*c_j]*weight;	
		
    if (next)
      traverseQ = *next <= *c;
    else
      traverseQ = false;
				
    // advance supp_0 and supp_1 until they are greater than or equal to current node
    while (traverseQ && supp.Next(s) && bVectors.Next()) {
      if (qIsAxisymmetric)
	circumf->Next();
      next = supp.PeekAhead(); 
      if (!next)
	traverseQ = false;
      else
	if (*next > *c)
	  traverseQ = false;
    }
			
    if (s != *c) { // means we're not at the end of the linked list
      supp.InsertAtCurrent(*c);
      bVectors.InsertAtCurrent(zeroFacet);
      if (qIsAxisymmetric) 
	circumf->InsertAtCurrent(0.);
      s = *c;
      if (supp.AtTop()) { // if we're inserting at the front, LinkedListT's behavior requires more work
	supp.Next(); 
	bVectors.Next();
	if (qIsAxisymmetric)
	  circumf->Next();
      }
    }
			
    double *currentI = facetIntegral.Pointer();
    double *currentB = bVectors.CurrentValue()->Pointer();
    for (int k = 0; k < nsd; k++)
      *currentB++ += *currentI++;
  }
}

void CellGeometryT::MergeNodalValues(int node_num, dArrayT& values,
						const iArrayT& neighbors, ArrayT< LinkedListT<int> >& suppWorkSpace, 
						ArrayT< LinkedListT<double> >& valWorkSpace, bool insertionQ) 
{
	/* If a node covering the integration point is not in the support of the node,
	 * insert that covering node into the sorted list.
	 */
	 
	iArrayT cover(neighbors);	
	int n_cover = cover.Length();	
	iArrayT cover_key(n_cover);
	cover_key.SetValueToPosition();
	cover_key.SortAscending(cover);			 

	LinkedListT<int>& supp = suppWorkSpace[node_num];
	LinkedListT<double>& vals = valWorkSpace[node_num]; 
	int s = -1;
	
	supp.Top(); 
	int *next = supp.CurrentValue();
	vals.Top();
		
	bool traverseQ;
	int* c = cover.Pointer();
	int* c_j = cover_key.Pointer();
	for (int j = 0; j < n_cover; j++, c++, c_j++) {
		
		if (next)
			traverseQ = *next <= *c;
		else
			traverseQ = false;
				
		// advance supp_0 and supp_1 until they are greater than or equal to current node
		while (traverseQ && supp.Next(s) && vals.Next()) {
			next = supp.PeekAhead(); 
			if (!next)
				traverseQ = false;
			else
				if (*next > *c)
					traverseQ = false;
		}
		
		if (insertionQ) {	
			if (s != *c) { // means we're not at the end of the linked list
				supp.InsertAtCurrent(*c);
				vals.InsertAtCurrent(0.0);
				s = *c;
				if (supp.AtTop()) { // if we're inserting at the front, LinkedListT's behavior requires more work
					supp.Next(); 
					vals.Next();
				}
			}
			(*vals.CurrentValue()) += values[*c_j];
		} else {
			if (s != *c)
				ExceptionT::GeneralFail("CellGeometryT::MergeNodalValues",
						"Node %d in support of node %d but not in data\n",s,*c);
			if (supp.AtTop())
				vals.Next();
			(*vals.CurrentValue()) = values[*c_j];
		}
	}
}

void CellGeometryT::ConfigureDataStructures(RaggedArray2DT<int>& cellSupports, RaggedArray2DT<dArrayT>& bVectors,
							RaggedArray2DT<double>& circumferential_B, dArrayT& cellVolumes)
{
	const char caller[] = "CellGeometry::FinishDataStuctures";

	int nNodes = fNodalCoordinates->MajorDim();
	
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
		int nsd = fNodalCoordinates->MinorDim();
		for (int i = 0; i < nNodes; i++) {
			nodal_init_coords.Alias(nsd, (*fNodalCoordinates)(i)); // This is the nodal coordinate.
			double R_i = nodal_init_coords[0];
			
			if (R_i > kSmall)
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
			
			const iArrayT& neighbors = fNodalShapes->Neighbors();	

			MergeNodalValues(i, phis, neighbors, nodeWorkSpace, 
				circumferentialWorkSpace, false);
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
	
	/* free linked list storage -- this might not need to happen if adapativity is an issue */
	for (int i = 0; i < nodeWorkSpace.Length(); i++) {
		nodeWorkSpace[i].Clear();
		facetWorkSpace[i].Clear();
		if (qIsAxisymmetric)
			circumferentialWorkSpace[i].Clear();
	}
}
