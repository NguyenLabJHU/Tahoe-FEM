/* $Id: VirtualSWDC.cpp,v 1.6 2004-06-17 07:41:37 paklein Exp $ */
/* created: paklein (05/05/1997) */
#include "VirtualSWDC.h"

#include <iomanip.h>

#include "ifstreamT.h"
#include "toolboxConstants.h"

/* fVNodePair decoding */

using namespace Tahoe;

const int kVirtualNode = 0;
const int kActiveNode  = 1;

/* constructor */
VirtualSWDC::VirtualSWDC(const ElementSupportT& support, const FieldT& field): 
	SWDiamondT(support, field)
{

}

/* append element equations numbers to the list */
void VirtualSWDC::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_2)

	/* set local 2-body equations numbers */
	Field().SetLocalEqnos(fNodes_2Body, fEqnos_2Body);

//	/* set periodic local equation numbers for 3-body */
//	iArray2DT temp_3Body; //temp copy to be modified
//	temp_3Body = fNodes_3Body;
//
//	/* swap nodes to get periodic local equation numbers */
//	SwapVirtualNodes(temp_3Body);
//		//swapping here means the changes are not
//		//accounted for in bandwidth reduction.
//
//	/* set local 2-body equations numbers */	
//	theNodes->SetLocalEqnos(temp_3Body, fEqnos_3Body);

	/* set local 3-body equations numbers */	
	Field().SetLocalEqnos(fPeriodicNodes_3Body, fEqnos_3Body);

	/* add to list */
	eq_1.Append(&fEqnos_3Body); //3 body is superset
		//still OK here since 2-body terms involving virtual nodes
		//don't contribute since the virtual nodes are inactive. (?)
}

/* appends group connectivities to the array */
void VirtualSWDC::ConnectsX(AutoArrayT<const iArray2DT*>&
	connects) const
{
	//need to return 2 connection sets, since periodic
	//conditions may leave some nodes unconnected
	connects.Append(&fPeriodicNodes_3Body);
	connects.Append(&fNodes_2Body);
}	
		
/***********************************************************************
* Protected
***********************************************************************/

void VirtualSWDC::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	/* inherited */
	SWDiamondT::EchoConnectivityData(in, out);

	/* virtual node pair data */
	int numvpairs;
	in >> numvpairs;
	out << " Number of virtual node pairs. . . . . . . . . . = " << numvpairs << "\n\n";
		
	if (numvpairs > 0)
	{
		/* memory */
		fVNodePairs.Dimension(numvpairs, 2);
	
		/* read data */
		fVNodePairs.ReadNumbered(in);
	
		/* echo to output */
		out << setw(kIntWidth) << "no."; 		
		out << setw(kIntWidth) << "node*"; 		
		out << setw(kIntWidth) << "node" << '\n'; 		
		fVNodePairs.WriteNumbered(out);
	}
	
	/* create periodic local node number lists */
	fPeriodicNodes_3Body = fNodes_3Body;
	//SwapVirtualNodes(fPeriodicNodes_3Body);
	SwapVirtualNodes2(fPeriodicNodes_3Body);
}

/***********************************************************************
* Private
***********************************************************************/

/* blind swap of all virtual node pairs */
void VirtualSWDC::SwapVirtualNodes(iArray2DT& elnodelist)
{
	iArrayT temp(elnodelist.Length(), elnodelist.Pointer());
	
	for (int i = 0; i < fVNodePairs.MajorDim(); i++)
	{
		/* replace node numbers */
		int n_active  = fVNodePairs(i,kActiveNode);
		int n_virtual = fVNodePairs(i,kVirtualNode);

		temp.ChangeValue(n_virtual, n_active);
	}
}

//can't do this since virtual nodes show up in
//interactions that are not periodic
void VirtualSWDC::SwapVirtualNodes2(iArray2DT& elnodelist)
{
	/* load local equation numbers - redundant */
	Field().SetLocalEqnos(fNodes_3Body, fEqnos_3Body);

	/* shallow work space */
	iArrayT localeqns, localnds;

	dArrayT X, Y, Z;

	/* loop over elements */
	for (int j = 0; j < fEqnos_3Body.MajorDim(); j++)
	{
		/* fetch element equations */
		fEqnos_3Body.RowAlias(j,localeqns);
			
		/* fetch node numbers */
		elnodelist.RowAlias(j,localnds);
	
		/* must have at least one active */
		//if (localeqns.Max() > 0)
		
		//middle node must be mapped and only
		//one slave node in the element
		//(or should be exactly one active?)
		//if (localeqns[4] == -99) &&
		//    (localeqns[1] != -98 || localeqns[7] != -98) )
		
		//get ref coords
		fLocX_3Body.SetLocal(localnds);
		
		X.Set(3, fLocX_3Body(0));
		Y.Set(3, fLocX_3Body(1));
		Z.Set(3, fLocX_3Body(2));
		
		if (X.Average() < 13.0 &&
		    Y.Average() < 13.0 &&
		    Z.Average() < 13.0 &&
		    localeqns.Min() < 0)
		{
			/* fetch node numbers */
			//elnodelist.RowAlias(j,localnds);
		
			/* blind swap all possible */
			for (int i = 0; i < fVNodePairs.MajorDim(); i++)
			{
				/* replace node numbers */
				int n_active  = fVNodePairs(i,kActiveNode);
				int n_virtual = fVNodePairs(i,kVirtualNode);
		
				/* safe swap */
				localnds.ChangeValue(n_virtual, n_active);
			}
		}
	}
}
