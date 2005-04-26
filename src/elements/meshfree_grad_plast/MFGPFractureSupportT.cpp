/* $Id: MFGPFractureSupportT.cpp,v 1.1 2005-04-26 22:32:28 kyonten Exp $ */
#include "MFGPFractureSupportT.h"
//#include "MeshFreeFractureSupportT.h"

#include "StringT.h"
#include "ifstreamT.h"
#include "D3MeshFreeShapeFunctionT.h"
#include "ParameterContainerT.h"

using namespace Tahoe;

/* constructor */
MFGPFractureSupportT::MFGPFractureSupportT(void);


/* destructor */
//MFGPFractureSupportT::~MFGPFractureSupportT(void);

/* information about subordinate parameter lists */
//void MFGPFractureSupportT::DefineSubs(SubListT& sub_list) const
//{
	/* inherited */
	//MFGPElementSupportT::DefineSubs(sub_list);
	//MeshFreeFractureSupportT::DefineSubs(sub_list);
//}

/* a pointer to the ParameterInterfaceT of the given subordinate */
//ParameterInterfaceT* MFGPFractureSupportT::NewSub(const StringT& name) const
//{
	/* inherited */
	//return MFGPElementSupportT::NewSub(name);
	//return MeshFreeFractureSupportT::NewSub(name);
//}

/* accept parameter list */
//void MFGPFractureSupportT::TakeParameterList(const ParameterListT& list)
//{
	/* inherited */
	//MFGPElementSupportT::TakeParameterList(list);
	//MeshFreeFractureSupportT::TakeParameterList(list);

	/* collect pointers to parameter lists needed during InitSupport */
//}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* initialization */
void MFGPFractureSupportT::InitSupport(ostream& out, AutoArrayT<ElementCardT>& elem_cards_displ,
    AutoArrayT<ElementCardT>& elem_cards_plast, const iArrayT& surface_nodes, int numDOF_displ, 
    int numDOF_plast, int max_node_num, ModelManagerT* model)
{
	/* inherited */
	MFGPElementSupportT::InitSupport(out, elem_cards_displ, elem_cards_plast, surface_nodes, 
				numDOF_displ, numDOF_plast, max_node_num, model);
}

