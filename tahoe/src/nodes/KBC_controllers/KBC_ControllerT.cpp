/* $Id: KBC_ControllerT.cpp,v 1.12 2004-06-28 22:41:56 hspark Exp $ */
/* created: paklein (09/05/2000) */
#include "KBC_ControllerT.h"

#include "NodeManagerT.h"
#include "FEManagerT.h"
#include "ModelManagerT.h"
#include "ifstreamT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<KBC_ControllerT>::fByteCopy = false;
DEFINE_TEMPLATE_STATIC const bool ArrayT<KBC_ControllerT*>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
KBC_ControllerT::KBC_ControllerT(NodeManagerT& node_manager):
	ParameterInterfaceT("KBC_controller"),
	fNodeManager(node_manager)
{

}

/* destructor */
KBC_ControllerT::~KBC_ControllerT(void) { }

/* initialization */
void KBC_ControllerT::WriteParameters(ostream& out) const
{
#pragma unused(out)
}

void KBC_ControllerT::ReadRestart(istream& in)
{
#pragma unused(in)
}

void KBC_ControllerT::WriteRestart(ostream& out) const
{
#pragma unused(out)
}

/* returns true if the internal force has been changed since
* the last time step */
GlobalT::RelaxCodeT KBC_ControllerT::RelaxSystem(void)
{
	return GlobalT::kNoRelax;
}

/* output current configuration */
void KBC_ControllerT::WriteOutput(ostream& out) const
{
#pragma unused(out)
}

/**********************************************************************
* Protected
**********************************************************************/

/* read nodes from stream */
void KBC_ControllerT::ReadNodes(ifstreamT& in, ArrayT<StringT>& id_list,
	iArrayT& nodes) const
{
	/* top level */
	const FEManagerT& fe_man = fNodeManager.FEManager();
	ModelManagerT* model = fe_man.ModelManager();

	/* read node set indexes */
	model->NodeSetList (in, id_list);

	/* collect sets */
	model->ManyNodeSets(id_list, nodes);
}
