/* $Id: roof.cpp,v 1.1.2.2 2003-05-03 09:08:27 paklein Exp $ */
#include "roof.h"

roof::roof(void):
	ParameterInterfaceT("roof"),
	style_(undefined)
{

}

void roof::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	ParameterT the_style(ParameterT::Enumeration, "style");
	the_style.AddEnumeration("shingle", shingle);
	the_style.AddEnumeration("slate", slate);
	list.AddParameter(the_style);
}

void roof::SetParameters(const ParameterListT& list)
{
	list.GetParameter("style", enum2int<style>(style_));
}
