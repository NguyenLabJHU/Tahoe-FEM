/* $Id: BasicFieldT.cpp,v 1.7 2003-08-14 05:34:21 paklein Exp $ */
#include "BasicFieldT.h"
#include "iArrayT.h"

using namespace Tahoe;

/* constructor */
BasicFieldT::BasicFieldT(void):
	fdArray2DGroup(0, false),
	fiArray2DGroup(0, false)
{

}

/* initialize field */
void BasicFieldT::Initialize(const StringT& name, int ndof, int order)
{
	fName = name;
	fField.Dimension(order+1);
	fEqnos.Dimension(0, ndof),
	fdArray2DGroup.SetMinorDimension(ndof);
	fiArray2DGroup.SetMinorDimension(ndof);

	/* set default labels */
	fLabels.Dimension(ndof);
	for (int i = 0; i < fLabels.Length(); i++)
		fLabels[i].Append("D_", i+1);

	/* register arrays */
	for (int i = 0; i < fField.Length(); i++)
		RegisterArray2D(fField[i]);
	RegisterArray2D(fEqnos);
}

/* set field labels */
void BasicFieldT::SetLabels(const ArrayT<StringT>& labels)
{
	fLabels = labels;
}

/* set number of nodes */
void BasicFieldT::Dimension(int nnd, bool copy_in)
{
	/* resize double arrays */
	fdArray2DGroup.SetMajorDimension(nnd, copy_in);
	
	/* resize integer arrays */
	fiArray2DGroup.SetMajorDimension(nnd, copy_in);
}

/* set all field values to 0.0 */
void BasicFieldT::Clear(void)
{
	for (int i = 0; i < fField.Length(); i++)
		fField[i] = 0.0;
}

void BasicFieldT::WriteEquationNumbers(ostream& out, const ArrayT<int>* node_map) const
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
