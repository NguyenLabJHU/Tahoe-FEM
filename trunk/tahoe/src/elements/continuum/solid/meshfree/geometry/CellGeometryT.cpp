/* $Id: CellGeometryT.cpp,v 1.3 2005-01-26 20:21:07 cjkimme Exp $ */
#include "CellGeometryT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructors */
CellGeometryT::CellGeometryT(const ElementSupportT& support, bool isAxisymmetric):
	ParameterInterfaceT("cell_geometry"),
	fElementSupport(&support),
	fNumIP(1),
	fscnimft(NULL),
	qIsAxisymmetric(isAxisymmetric)
{

}

/* constructors */
CellGeometryT::CellGeometryT(void):
	ParameterInterfaceT("cell_geometry"),
	fElementSupport(NULL),
	fNumIP(1),
	fscnimft(NULL),
	qIsAxisymmetric(false)
{

}

/* destructor */
CellGeometryT::~CellGeometryT(void)
{

}

void CellGeometryT::DefineElements(const ArrayT<StringT>& block_ID, const ArrayT<int>& mat_index) 
{
#pragma unused(block_ID)
#pragma unused(mat_index)
}

void CellGeometryT::SetNodesAndShapes(iArrayT& nodes, dArray2DT& nodal_coordinates, MeshFreeNodalShapeFunctionT* nodalShapeFunctions)
{
	fNodes.Alias(nodes);
	fNodalCoordinates.Alias(nodal_coordinates);
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
						iArrayT& ip_cover, iArrayT& ip_cover_key, ArrayT< LinkedListT<int> >& nodeWorkSpace, 
						ArrayT< LinkedListT<dArrayT> >& facetWorkSpace,
						ArrayT< LinkedListT<double> >& circumferentialWorkSpace) 
{
	/* If a node covering the integration point is not in the support of the node,
	 * insert that covering node into the sorted list.
	 */			 

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
	int n_cover = phiValues.Length();
	int* c = ip_cover.Pointer();
	int* c_j = ip_cover_key.Pointer();
	int nsd = facetNormal.Length();
	dArrayT facetIntegral(nsd), zeroFacet(3);
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
			//if (qIsAxisymmetric)
			//	circumf->Next();
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