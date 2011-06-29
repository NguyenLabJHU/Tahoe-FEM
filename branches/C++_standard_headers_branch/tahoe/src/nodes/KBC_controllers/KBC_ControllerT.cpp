/* $Id: KBC_ControllerT.cpp,v 1.18 2009-05-26 12:29:41 tdnguye Exp $ */
/* created: paklein (09/05/2000) */
#include "KBC_ControllerT.h"
#include "BasicSupportT.h"
#include "ModelManagerT.h"
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
	if (strcmp("K-field", name) == 0)
		return kK_Field;
	else if (strcmp("bi-material_K-field", name) == 0)
		return kBimaterialK_Field;
//	else if (strcmp("K-field_3D", name) == 0)
//		return kK_Field_3D;
	else if (strcmp("torsion", name) == 0)
		return kTorsion;
	else if (strcmp("mapped_nodes", name) == 0)
		return kMappedPeriodic;
	else if (strcmp("scaled_velocity", name) == 0)
		return kScaledVelocityNodes;
	else if (strcmp("tied_nodes", name) == 0)
		return kTiedNodes;
	else if (strcmp("periodic_nodes", name) == 0)
		return kPeriodicNodes;
	else if (strcmp("conveyor", name) == 0)
		return kConveyor;
	else if (strcmp("symmetric_conveyor", name) == 0)
		return kConveyorSym;
//	else if (strcmp("angled_bc", name) == 0)
//		return kAngledBC;
	else
		return kNone;
}

/* constructor */
KBC_ControllerT::KBC_ControllerT(const BasicSupportT& support):
	ParameterInterfaceT("KBC_controller"),
	fSupport(support)
{

}

/* destructor */
KBC_ControllerT::~KBC_ControllerT(void) { }

void KBC_ControllerT::ReadRestart(ifstreamT& in)
{
#pragma unused(in)
}

void KBC_ControllerT::WriteRestart(ofstreamT& out) const
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

/* read nodes in node sets */
void KBC_ControllerT::GetNodes(const ArrayT<StringT>& id_list, iArrayT& nodes) const
{
	/* get the model */
	ModelManagerT& model = fSupport.ModelManager();

	/* collect sets */
	model.ManyNodeSets(id_list, nodes);
}
