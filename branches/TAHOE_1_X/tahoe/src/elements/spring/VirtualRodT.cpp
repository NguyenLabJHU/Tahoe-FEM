/* $Id: VirtualRodT.cpp,v 1.7 2004-06-17 07:41:37 paklein Exp $ */
/* created: paklein (05/01/1997) */

#include "VirtualRodT.h"

#include <iomanip.h>

#include "ifstreamT.h"
#include "toolboxConstants.h"

/* decoding VElPair data */

using namespace Tahoe;

const int kBoundaryNode = 0;
const int kVirtualNode  = 1;
const int kActiveNode   = 2;

/* constructor */
VirtualRodT::VirtualRodT(const ElementSupportT& support, const FieldT& field): 
	UnConnectedRodT(support, field)
{

}

/* append element equations numbers to the list */
void VirtualRodT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_2)

	/* substitute in virtual node numbers */
	iArray2DT tempnodes; //temp copy to be modified
	tempnodes = *(fConnectivities[0]);
	
	/* swap nodes to get periodic local equation numbers */
	SwapVirtualNodes(tempnodes);
		//swapping here means the changes are not
		//accounted for in bandwidth reduction.

	/* set local equations numbers */
	Field().SetLocalEqnos(tempnodes, fEqnos[0]);

	/* add to list */
	eq_1.Append(&fEqnos[0]);
}

/***********************************************************************
* Protected
***********************************************************************/

/* element data */
void VirtualRodT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	/* inherited */
	UnConnectedRodT::EchoConnectivityData(in, out);

	/* virtual node triplets data */
	int numtriplets;
	in >> numtriplets;	
	out << " Number of virtual node triplets . . . . . . . . = " << numtriplets << "\n\n";
		
	if (numtriplets > 0)
	{
		/* memory */
		fVNodeTriplets.Dimension(numtriplets, 3);
	
		/* read data */
		fVNodeTriplets.ReadNumbered(in);
	
		/* echo to output */
		out << setw(kIntWidth) << "no."; 		
		out << setw(kIntWidth) << "n[1]"; 		
		out << setw(kIntWidth) << "n[2]"; 		
		out << setw(kIntWidth) << "n[2]" << "*" << '\n'; 		
	
		fVNodeTriplets.WriteNumbered(out);
	}
}

/***********************************************************************
* Private
***********************************************************************/

/* swaps the equation numbers for every occurence
* of the pair in this group */
void VirtualRodT::SwapVirtualNodes(iArray2DT& elnodelist) const
{
	/* shallow work space */
	iArrayT nodelist;
	
	int nen = NumElements();
	for (int i = 0; i < fVNodeTriplets.MajorDim(); i++)
	{
		/* warning flag */
		int found = 0;
		
		for (int el = 0; el < nen; el++)	
		{
			/* fetch local node list */
			elnodelist.RowAlias(el, nodelist);			

			if ( nodelist.HasValue( fVNodeTriplets(i,kBoundaryNode) ) &&
			     nodelist.HasValue( fVNodeTriplets(i,kVirtualNode ) ) )
			{
				found = 1;

				/* replace node numbers */
				int n_active  = fVNodeTriplets(i,kActiveNode);
				int n_virtual = fVNodeTriplets(i,kVirtualNode);
								
				nodelist.ChangeValue(n_virtual, n_active);			
			}
		}
			
		if (!found)
		{
			fVNodeTriplets.RowAlias(i,nodelist);
		
			cout << "\nVirtualRodT::SwapVirtualNodes: virtual node set:";
			cout << nodelist << "not found." << endl;
		}
	}
}

/* blind swap of all virtual node pairs */
void VirtualRodT::SwapVirtualNodes2(iArray2DT& elnodelist)
{
	iArrayT temp(elnodelist.Length(), elnodelist.Pointer());
	
	for (int i = 0; i < fVNodeTriplets.MajorDim(); i++)
	{
		/* replace node numbers */
		int n_active  = fVNodeTriplets(i,kActiveNode);
		int n_virtual = fVNodeTriplets(i,kVirtualNode);

		cout << "Count of: " << n_virtual << " is ";
		cout << temp.Count(n_virtual) << '\n';

		temp.ChangeValue(n_virtual, n_active);
	}
}
