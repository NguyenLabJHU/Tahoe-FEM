/* $Id: FieldSupportT.cpp,v 1.4.16.3 2004-03-31 16:14:37 paklein Exp $ */
#include "FieldSupportT.h"
#include "NodeManagerT.h"

using namespace Tahoe;

/* constructor */
FieldSupportT::FieldSupportT(const FEManagerT& fe)
{
	SetFEManager(&fe);
}

/* construct new KBC controller */
KBC_ControllerT* FieldSupportT::NewKBC_Controller(FieldT& field, int code) const {
	NodeManagerT& nodes = const_cast<NodeManagerT&>(NodeManager());
	return nodes.NewKBC_Controller(field, code);
}

FBC_ControllerT* FieldSupportT::NewFBC_Controller(FieldT& field, int code) const {
	NodeManagerT& nodes = const_cast<NodeManagerT&>(NodeManager());
	return nodes.NewFBC_Controller(field, code);
}
