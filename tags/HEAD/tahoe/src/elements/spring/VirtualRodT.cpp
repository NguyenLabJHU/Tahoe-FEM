/* $Id: VirtualRodT.cpp,v 1.1.1.1 2001-01-29 08:20:34 paklein Exp $ */
/* created: paklein (05/01/1997)                                          */
/* UnConnectedRodT plus virtual elements for periodic boundary            */
/* conditions.                                                            */

#include "VirtualRodT.h"

#include <iomanip.h>

#include "fstreamT.h"
#include "Constants.h"
#include "NodeManagerT.h"

/* decoding VElPair data */
const int kBoundaryNode = 0;
const int kVirtualNode  = 1;
const int kActiveNode   = 2;

/* constructor */
VirtualRodT::VirtualRodT(FEManagerT& fe_manager): UnConnectedRodT(fe_manager)
{

}

/* append element equations numbers to the list */
void VirtualRodT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_2)

	/* substitute in virtual node numbers */
	iArray2DT tempnodes; //temp copy to be modified
	tempnodes = fConnectivities;
	
	/* swap nodes to get periodic local equation numbers */
	SwapVirtualNodes(tempnodes);
		//swapping here means the changes are not
		//accounted for in bandwidth reduction.

	/* set local equations numbers */
	fNodes->SetLocalEqnos(tempnodes, fEqnos);

	/* add to list */
	eq_1.Append(&fEqnos);
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
		fVNodeTriplets.Allocate(numtriplets, 3);
	
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
	
	for (int i = 0; i < fVNodeTriplets.MajorDim(); i++)
	{
		/* warning flag */
		int found = 0;
		
		for (int el = 0; el < fNumElements; el++)	
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
