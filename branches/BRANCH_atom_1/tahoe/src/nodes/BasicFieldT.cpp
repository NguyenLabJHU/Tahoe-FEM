/* $Id: BasicFieldT.cpp,v 1.4 2002-09-20 16:14:06 paklein Exp $ */
#include "BasicFieldT.h"
#include "iArrayT.h"

/* constructor */

using namespace Tahoe;

BasicFieldT::BasicFieldT(const StringT& name, int ndof, int order):
	fName(name),
	fField(order+1),
	fEqnos(0, ndof)
{
	/* set default labels */
	fLabels.Dimension(ndof);
	for (int i = 0; i < fLabels.Length(); i++)
		fLabels[i].Append("D_", i+1);
}

/* set field labels */
void BasicFieldT::SetLabels(const ArrayT<StringT>& labels)
{
	fLabels = labels;
}

/* set number of nodes */
void BasicFieldT::Dimension(int nnd)
{
	/* number of degrees of freedom */
	int ndof = fEqnos.MinorDim();

	/* dimension field */
	for (int i = 0; i < fField.Length(); i++)
		fField[i].Dimension(nnd, ndof);
		
	/* allocate equations array */
	fEqnos.Dimension(nnd, ndof);
}

void BasicFieldT::WriteEquationNumbers(ostream& out, const iArrayT* node_map) const
{
	/* dimensions */
	int nnd  = NumNodes();
	int ndof = NumDOF();
	int columns = 1;
	
	/* print header */
	out << "\n Field: \"" << fName << "\"\n\n";
	for (int k = 0; k < columns; k++)
	{
		out << setw(kIntWidth) << "node";
		out << setw(kIntWidth) << "map";
		for (int j = 0; j < ndof; j++)
		{
			out << setw(kIntWidth - 2) << "d[";
			out << j+1 << "]";
		}
		out << "    ";
	}
	out << "\n\n";

	/* print data */
	int colcount = 0;
	for (int i = 0; i < nnd; i++)
	{
		out << setw(kIntWidth) << i+1;
		out << setw(kIntWidth) << ((node_map != NULL) ? (*node_map)[i]: i) + 1;
		for (int j = 0; j < ndof; j++)
			out << setw(kIntWidth) << fEqnos(i,j);
			
		/* wrap */
		if (++colcount == columns)
		{
			out << '\n';
			colcount = 0;
		}
		else
			out << "    ";	
	}
	
	if (colcount != 0) out << '\n'; //terminate mid-line
	out << '\n';
}
