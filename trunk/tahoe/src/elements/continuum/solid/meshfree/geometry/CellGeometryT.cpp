/* $Id: CellGeometryT.cpp,v 1.1 2005-01-25 02:24:04 cjkimme Exp $ */
#include "CellGeometryT.h"

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

void CellGeometryT::SetNodesAndShapes(dArray2DT& nodal_coordinates, MeshFreeNodalShapeFunctionT* nodalShapeFunctions)
{
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