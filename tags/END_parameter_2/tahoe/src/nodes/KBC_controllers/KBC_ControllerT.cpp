/* $Id: KBC_ControllerT.cpp,v 1.11.6.3 2004-03-27 04:18:01 paklein Exp $ */
/* created: paklein (09/05/2000) */
#include "KBC_ControllerT.h"

#include "NodeManagerT.h"
#include "FEManagerT.h"
#include "ModelManagerT.h"
#include "fstreamT.h"

#include <string.h>

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<KBC_ControllerT>::fByteCopy = false;
DEFINE_TEMPLATE_STATIC const bool ArrayT<KBC_ControllerT*>::fByteCopy = true;
} /* namespace Tahoe */

/* converts strings to KBC_ControllerT::CodeT */
KBC_ControllerT::CodeT KBC_ControllerT::Code(const char* name)
{
	if (strcmp("K_field", name) == 0)
		return kK_Field;
	else if (strcmp("torsion", name) == 0)
		return kTorsion;
	else if (strcmp("mapped_nodes", name) == 0)
		return kMappedPeriodic;
	else
		return kNone;
}

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
#pragma message("delete me")

	/* top level */
	const FEManagerT& fe_man = fNodeManager.FEManager();
	ModelManagerT* model = fe_man.ModelManager();

	/* read node set indexes */
	model->NodeSetList (in, id_list);

	/* collect sets */
	model->ManyNodeSets(id_list, nodes);
}

/* read nodes from stream */
void KBC_ControllerT::GetNodes(const ArrayT<StringT>& id_list, iArrayT& nodes) const
{
	/* get the model */
	const FEManagerT& fe_man = fNodeManager.FEManager();
	ModelManagerT* model = fe_man.ModelManager();

	/* collect sets */
	model->ManyNodeSets(id_list, nodes);
}
