/* $Id: DotLine_FormatterT.cpp,v 1.3 2003-11-10 22:14:41 cjkimme Exp $ */
#include "DotLine_FormatterT.h"
#include "ParameterListT.h"
#include "ParameterT.h"

#include <iomanip.h>
#include <math.h>

using namespace Tahoe;
const int kLineLength = 50;

static int i_min(int a, int b) { return (a < b) ? a : b; };

DotLine_FormatterT::DotLine_FormatterT(void):
	fTabWidth(0),
	fTab(kLineLength-10),
	fDots(kLineLength+1)
{
	fDots[0] = fDots[kLineLength] = '\0';
	fTab.Fill(' ');
}

bool DotLine_FormatterT::InitParameterFile(ostream& out) const
{
	out << "# Generated by Tahoe::DotLine_FormatterT $Revision: 1.3 $\n";
	return true;
}

bool DotLine_FormatterT::WriteParameterList(ostream& out, const ParameterListT& list) const
{
	out << '\n';

	/* non-const this */
	DotLine_FormatterT* non_const_this = (DotLine_FormatterT*) this;

	/* set parameter path */
	int n_sep = 0;
	if (fPath.StringLength() > 0) {
		non_const_this->fPath.Append("::");
		n_sep = 2;
	}
	non_const_this->fPath.Append(list.Name());
	out << Tab() << "begin: " << fPath << '\n';

	/* name and description */
	if (list.Description().StringLength() > 0)
		out << Tab() << "description: " << list.Description() << '\n';
	
	/* parameters */
	const ArrayT<ParameterT>& params = list.Parameters();
	for (int i = 0; i < params.Length(); i++) {
		out << Tab() << params[i].Name() << Dots(params[i].Name()) << " = " << params[i] << '\n';
		if (params[i].Description().StringLength() > 0)
			out << Tab() << "description: " << params[i].Description() << '\n';
	}

	/* nested parameter lists */
	const ArrayT<ParameterListT>& nested_lists = list.Lists();
	for (int i = 0; i < nested_lists.Length(); i++)
	{
		TabOut();
		WriteParameterList(out, nested_lists[i]);
		TabIn();
	}

	/* restore the path */
	if (nested_lists.Length() > 0) out << '\n';
	out << Tab() << "end: " << fPath << '\n';
	non_const_this->fPath.Drop(-(list.Name().StringLength() + n_sep));
	return true;
}

bool DotLine_FormatterT::CloseParameterFile(ostream& out) const
{
#pragma unused(out)
	return true;
}

bool DotLine_FormatterT::InitDescriptionFile(ostream& out) const 
{
#pragma unused(out)
	return false;
}
bool DotLine_FormatterT::CloseDescriptionFile(ostream& out) const 
{
#pragma unused(out)

	return false;
}
bool DotLine_FormatterT::WriteDescription(ostream& out, const ParameterListT& list) const 
{
#pragma unused(out)
#pragma unused(list)

	return false;
}

/*************************************************************************
 * Private
 *************************************************************************/

/* return a row of dots which pads the given string */
const StringT& DotLine_FormatterT::Dots(const StringT& str) const
{
	DotLine_FormatterT* non_const_this = (DotLine_FormatterT*) this;

	int tab_length = Tab().StringLength();

	StringT& dots = non_const_this->fDots;
	char* p = dots.Pointer();
	int flip = int(fmod(str.StringLength(), 2.0));
	for (int i = str.StringLength()+tab_length; i < kLineLength; i++)
	{
		if (flip == 1) {
			*p++ = '.';
			flip = 0;
		} else {
			*p++ = ' ';
			flip = 1;
		}
	}
	
	/* terminate */
	*p = '\0';

	return fDots;
}

/* return a tabbing string */
const StringT& DotLine_FormatterT::Tab(void) const
{
	DotLine_FormatterT* non_const_this = (DotLine_FormatterT*) this;
	StringT& tab = non_const_this->fTab;

	/* set tab */
	tab.Fill(' ');
	int tab_space = i_min(Depth()*fTabWidth, tab.Length()-1);
	tab[tab_space] = '\0';
	
	return fTab;
}
