/* $Id: driveway.cpp,v 1.1.2.2 2003-05-03 09:08:27 paklein Exp $ */
#include "driveway.h"

driveway::driveway(void):
	ParameterInterfaceT("driveway"),
	length_(0.0),
	surface_(undefined)
{

}

void driveway::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	ParameterT the_length(ParameterT::Double, "length");
	the_length.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(the_length);

	ParameterT the_surface(ParameterT::Enumeration, "surface");
	the_surface.AddEnumeration("dirt", dirt);
	the_surface.AddEnumeration("gravel", gravel);
	the_surface.AddEnumeration("cobblestone", cobblestone);
	the_surface.AddEnumeration("asphalt", asphalt);
	the_surface.AddEnumeration("concrete", concrete);
	the_surface.SetDefault(asphalt);
	list.AddParameter(the_surface);
}

void driveway::SetParameters(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::SetParameters(list);

	list.GetParameter("length", length_);
	list.GetParameter("surface", enum2int<surface>(surface_));
}
