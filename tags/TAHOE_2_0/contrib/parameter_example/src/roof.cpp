/* $Id: roof.cpp,v 1.2 2003-05-04 22:49:50 paklein Exp $ */
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

void roof::TakeParameterList(const ParameterListT& list)
{
	list.GetParameter("style", enum2int<style>(style_));
}
