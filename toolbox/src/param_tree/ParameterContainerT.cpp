/* $Id: ParameterContainerT.cpp,v 1.5 2005-04-05 15:51:11 paklein Exp $ */
#include "ParameterContainerT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<ParameterContainerT>::fByteCopy = false;
DEFINE_TEMPLATE_STATIC const bool ArrayT<ParameterContainerT*>::fByteCopy = true;
}

/* constructor */
ParameterContainerT::ParameterContainerT(const StringT& name):
	ParameterInterfaceT(name),
	fSubSource(NULL)
{
	/* take default from ParameterListT */
	ParameterListT tmp;
	fListOrder = tmp.ListOrder();
	fInline = tmp.Inline();
}

ParameterContainerT::ParameterContainerT(void):
	ParameterInterfaceT("container"),
	fSubSource(NULL)
{
	/* take default from ParameterListT */
	ParameterListT tmp;
	fListOrder = tmp.ListOrder();
	fInline = tmp.Inline();
}

/* add parameter */
bool ParameterContainerT::AddParameter(const ParameterT& param, ParameterListT::OccurrenceT occur)
{
	/* add to the list */
	fParameters.Append(param);
	fParametersOccur.Append(occur);
	return true;
}

void ParameterContainerT::AddSub(const StringT& name, ParameterListT::OccurrenceT occur, 
	bool is_inline)
{
	/* order info */
	fs_or_c.Append('s');
	fs_or_c_index.Append(fSubs.Length());

	fSubs.AddSub(name, occur, is_inline);
}

void ParameterContainerT::AddSub(const SubListDescriptionT& sub)
{
	/* order info */
	fs_or_c.Append('s');
	fs_or_c_index.Append(fSubs.Length());

	fSubs.AddSub(sub);
}

void ParameterContainerT::AddSub(const ParameterContainerT& sub, ParameterListT::OccurrenceT occur, bool is_inline)
{
	/* order info */
	fs_or_c.Append('c');
	fs_or_c_index.Append(fContainers.Length());

	fContainers.Append(sub);
	fContainersOccur.Append(occur);
	fContainersInline.Append(is_inline);
}

/* set source for subs not defined by the container */
void ParameterContainerT::SetSubSource(const ParameterInterfaceT* sub_source) 
{
	if (sub_source == this)
		ExceptionT::GeneralFail("ParameterContainerT::SetSubSource", 
			"the sub source for %s cannot be self", Name().Pointer());
	fSubSource = sub_source; 
};

/* a pointer to the ParameterInterfaceT of the given subordinate*/
ParameterInterfaceT* ParameterContainerT::NewSub(const StringT& list_name) const
{
	/* inherited (get from self) */
	ParameterInterfaceT* sub = ParameterInterfaceT::NewSub(list_name);
	
	/* check list of containers */
	for (int i = 0; i < fContainers.Length(); i++)
		if (fContainers[i].Name() == list_name) {
			ParameterContainerT* container = new ParameterContainerT(fContainers[i]);
			return container;
		}

	/* get from sub source */
	if (!sub && fSubSource)
		sub = fSubSource->NewSub(list_name);
		
	return sub;
}

/* return the description of the given inline subordinate parameter list. */
void ParameterContainerT::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
	SubListT& sub_sub_list) const
{
	/* inherited (get from self) */
	ParameterInterfaceT::DefineInlineSub(sub, order, sub_sub_list);

	/* if list is empty, try sub source */
	if (fSubSource && sub_sub_list.Length() == 0)
		fSubSource->DefineInlineSub(sub, order, sub_sub_list);
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* describe the parameters needed by the interface. */
void ParameterContainerT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	/* set description */
	list.SetDescription(fDescription);

	/* register all parameters */
	for (int i = 0; i < fParameters.Length(); i++)
		list.AddParameter(fParameters[i], fParametersOccur[i]);
}

/* information about subordinate parameter lists */
void ParameterContainerT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* register all sublists */
	for (int i = 0; i < fs_or_c.Length(); i++) {
		int ii = fs_or_c_index[i];
		if (fs_or_c[i] == 's')
			sub_list.AddSub(fSubs[ii]);
		else
			sub_list.AddSub(fContainers[ii].Name(), fContainersOccur[ii], fContainersInline[ii]);
	}
}
