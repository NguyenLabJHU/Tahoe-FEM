/* $Id: UnConnectedRodT.cpp,v 1.2 2001-10-25 07:16:43 paklein Exp $ */
/* created: paklein (04/05/1997)                                          */

#include "UnConnectedRodT.h"

#include <iomanip.h>

#include "fstreamT.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "FindNeighborT.h"

/* constructor */
UnConnectedRodT::UnConnectedRodT(FEManagerT& fe_manager):
	RodT(fe_manager),
	fNumNodesUsed(0),
	fReconnectCount(0)
{
	/* read neighbor list parameters */
	fFEManager.Input() >> fReconnectInc >> fMaxNeighborCount >> fNeighborDist;

	/* checks */
	if (fMaxNeighborCount <  1  ) throw eBadInputValue;
	if (fNeighborDist     <= 0.0) throw eBadInputValue;
}

/* apply pre-conditions at the current time step.  Signal
* all listeners that the time has just been incremented */
void UnConnectedRodT::InitStep(void)
{
	/* inherited */
	RodT::InitStep();
	
	/* increment reconnection count */;
	fReconnectCount++;
}

/* resets to the last converged solution */
void UnConnectedRodT::ResetStep(void)
{
	/* pre-condition */
	fReconnectCount--;
	if (fReconnectCount < 0) throw eGeneralFail;
		// condition implies that equilibrium could not be
		// established with a reconnected system for which
		// there is no last converged solution to go back to
}

/* element level reconfiguration for the current solution */
GlobalT::RelaxCodeT UnConnectedRodT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	/* redetermine connectivities */
	if (++fReconnectCount == fReconnectInc)
	{
		if (fNumNodesUsed != -1) throw eGeneralFail;
			//not storing NodesUsed yet
			//so cannot reconnect.	 	
	
		/* re-connect - more neighbors and greater distance */
		FindNeighborT Connector(fNodes->CurrentCoordinates(), fMaxNeighborCount);
		Connector.GetNeighors(fConnectivities, fNeighborDist);

		/* reset local equation number lists */	
		ConfigureElementData();

		/* set block data */
		fBlockData(0, kBlockDim) = fConnectivities.MajorDim();
		
		/* reset count */
		fReconnectCount = 0;
		
		/* precedence */
		return GlobalT::MaxPrecedence(relax, GlobalT::kReEQRelax);
	}
	
	/* base class code falls through */
	return relax;
}

/***********************************************************************
* Protected
***********************************************************************/

/* print data */
void UnConnectedRodT::PrintControlData(ostream& out) const
{
	/* inherited */
	RodT::PrintControlData(out);
	
	out << " Reconnection increment. . . . . . . . . . . . . = " << fReconnectInc     << '\n';
	out << " Maximum number of neighbors . . . . . . . . . . = " << fMaxNeighborCount << '\n';
	out << " Neighbor cut-off distance . . . . . . . . . . . = " << fNeighborDist     << '\n';
}
	
/* element data */
void UnConnectedRodT::ReadMaterialData(ifstreamT& in)
{
	/* read element data */
	RodT::ReadMaterialData(in);
	
	/* should only be one material */
	if (fMaterialsList.Length() != 1) throw eGeneralFail;
	
	/* set current material (once) */
	fCurrMaterial = fMaterialsList[0];
}


void UnConnectedRodT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	in >> fNumNodesUsed;
	if (fNumNodesUsed != -1 && fNumNodesUsed < 1) throw eBadInputValue;

	/* read nodes used */
	if (fNumNodesUsed == -1) //use ALL nodes
	{
		/* connector */
		FindNeighborT Connector(fNodes->CurrentCoordinates(), fMaxNeighborCount);
	
		/* connect nodes - dimensions lists */
		Connector.GetNeighors(fConnectivities, fNeighborDist);
	}
	else                      //only use specified nodes
	{
		/* read specified nodes */
		iArrayT nodesused(fNumNodesUsed);
		in >> nodesused;
		
		/* echo data */
		out << "\n Number of search nodes. . . . . . . . . . . . . = ";
		out << nodesused.Length() << "\n\n";
		out << nodesused.wrap(5) << '\n';
		out << '\n';
		
		/* connector */
		nodesused--;
		FindNeighborT Connector(nodesused, fNodes->CurrentCoordinates(),
							    fMaxNeighborCount);
	
		/* connect nodes - dimensions lists */
		Connector.GetNeighors(fConnectivities, fNeighborDist);		
	}
			
	/* set element equation and node lists */
	ConfigureElementData();

	/* set block data */
	fBlockData.Allocate(1, kBlockDataSize);
	fBlockData(0, kStartNum) = 0;
	fBlockData(0, kBlockDim) = fConnectivities.MajorDim();
	fBlockData(0, kBlockMat) = 0;

	/* print connectivity data */
	WriteConnectivity(out);
}

/***********************************************************************
* Private
***********************************************************************/

/* call AFTER 2 and 3 body node lists are set */
void UnConnectedRodT::ConfigureElementData(void)
{
	fNumElements = fConnectivities.MajorDim();

	/* allocate memory */
	fElementCards.Allocate(fNumElements);
	fEqnos.Allocate(fNumElements, fNumElemEqnos);

	/* set 2 body element data */
	for (int i = 0; i < fNumElements; i++)	
	{
		/* node and equation numbers */			
		(fElementCards[i].NodesX()).Set(fNumElemNodes, fConnectivities(i) );		
		(fElementCards[i].Equations()).Set(fNumElemEqnos, fEqnos(i) );
		
		/* all have same material number */
		fElementCards[i].SetMaterialNumber(0);
	}
}
