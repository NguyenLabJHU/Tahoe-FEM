/* $Id: FieldSupportT.cpp,v 1.5 2004-07-15 08:31:09 paklein Exp $ */
#include "FieldSupportT.h"
#include "NodeManagerT.h"

using namespace Tahoe;

/* constructor */
FieldSupportT::FieldSupportT(void)
{

}

/* construct new KBC controller */
KBC_ControllerT* FieldSupportT::NewKBC_Controller(FieldT& field, int code) const {
	NodeManagerT& nodes = const_cast<NodeManagerT&>(NodeManager());
	return nodes.NewKBC_Controller(field, code);
}

FBC_ControllerT* FieldSupportT::NewFBC_Controller(int code) const {
	NodeManagerT& nodes = const_cast<NodeManagerT&>(NodeManager());
	return nodes.NewFBC_Controller(code);
}
