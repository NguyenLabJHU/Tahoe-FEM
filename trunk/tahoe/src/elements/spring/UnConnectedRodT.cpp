/* $Id: UnConnectedRodT.cpp,v 1.15 2004-06-17 07:41:37 paklein Exp $ */
/* created: paklein (04/05/1997) */

#include "UnConnectedRodT.h"

#include <iomanip.h>

#include "ifstreamT.h"
#include "ModelManagerT.h"
#include "FindNeighborT.h"

/* constructor */

using namespace Tahoe;

UnConnectedRodT::UnConnectedRodT(const ElementSupportT& support, const FieldT& field):
	RodT(support, field),
	fNumNodesUsed(0),
	fReconnectCount(0)
{
	/* read neighbor list parameters */
	ElementSupport().Input() >> fReconnectInc >> fMaxNeighborCount >> fNeighborDist;

	/* checks */
	if (fMaxNeighborCount <  1  ) throw ExceptionT::kBadInputValue;
	if (fNeighborDist     <= 0.0) throw ExceptionT::kBadInputValue;
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
GlobalT::RelaxCodeT UnConnectedRodT::ResetStep(void)
{
	/* pre-condition */
	fReconnectCount--;
	if (fReconnectCount < 0) throw ExceptionT::kGeneralFail;
		// condition implies that equilibrium could not be
		// established with a reconnected system for which
		// there is no last converged solution to go back to

	/* inherited */
	return RodT::ResetStep();
}

/* element level reconfiguration for the current solution */
GlobalT::RelaxCodeT UnConnectedRodT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	/* redetermine connectivities */
	if (++fReconnectCount == fReconnectInc)
	{
		if (fNumNodesUsed != -1) throw ExceptionT::kGeneralFail;
			//not storing NodesUsed yet
			//so cannot reconnect.	 	
	
		/* re-connect - more neighbors and greater distance */
		iArray2DT rodconnects;
		FindNeighborT Connector(ElementSupport().CurrentCoordinates(), fMaxNeighborCount);
		Connector.GetNeighors(rodconnects, fNeighborDist);
		
		/* update model manager */
		ModelManagerT& model = ElementSupport().Model();
		model.UpdateElementGroup(fBlockData[0].ID(), rodconnects, true);
		fBlockData[0].SetDimension(rodconnects.MajorDim());

		/* reset local equation number lists */	
		ConfigureElementData();
		
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
	if (fMaterialsList.Length() != 1) throw ExceptionT::kGeneralFail;
	
	/* set current material (once) */
	fCurrMaterial = fMaterialsList[0];
}


void UnConnectedRodT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	in >> fNumNodesUsed;
	if (fNumNodesUsed != -1 && fNumNodesUsed < 1) throw ExceptionT::kBadInputValue;

	/* temp space */
	iArray2DT rodconnects;

	/* read nodes used */
	if (fNumNodesUsed == -1) //use ALL nodes
	{
		/* connector */
		FindNeighborT Connector(ElementSupport().CurrentCoordinates(), fMaxNeighborCount);
	
		/* connect nodes - dimensions lists */
		Connector.GetNeighors(rodconnects, fNeighborDist);
		
		/* list of all the nodes in the model */
		fGroupNodes.Dimension(ElementSupport().NumNodes());
		fGroupNodes.SetValueToPosition();
	}
	else //only use specified nodes
	{
		/* model information */
		ModelManagerT& model = ElementSupport().Model();
	
		/* model format specific */
		if (model.DatabaseFormat() == IOBaseT::kTahoe)
		{
			/* read specified nodes */
			fGroupNodes.Dimension(fNumNodesUsed);
			in >> fGroupNodes;
			fGroupNodes--;
		}
		else 
		{
			/* read node set ID's */
			ArrayT<StringT> node_ids(fNumNodesUsed);
			for (int i = 0; i < node_ids.Length(); i++)
				in >> node_ids[i];
		
			/* collect nodes */
			model.ManyNodeSets (node_ids, fGroupNodes);
		}	
	
		/* echo data */
		fGroupNodes++;
		out << "\n Number of search nodes. . . . . . . . . . . . . = ";
		out << fGroupNodes.Length() << "\n\n";
		out << fGroupNodes.wrap(5) << '\n';
		out << '\n';
		fGroupNodes--;
		
		/* connector */
		FindNeighborT Connector(fGroupNodes, 
			ElementSupport().CurrentCoordinates(), 
			fMaxNeighborCount);
	
		/* connect nodes - dimensions lists */
		Connector.GetNeighors(rodconnects, fNeighborDist);		
	}

	/* send connectivity data to ModelManagerT */
	ModelManagerT& model = ElementSupport().Model();
	StringT name("URod");
	name.Append(ElementSupport().ElementGroupNumber(this) + 1);
	GeometryT::CodeT code = GeometryT::kLine;
	model.RegisterElementGroup (name, rodconnects, code, true);

	/* set up fBlockData to store block ID */
	fBlockData.Dimension(1);
	fBlockData[0].Set(name, 0, rodconnects.MajorDim(), 0); // currently assume all interactions use potential 0
//	fNumElements = rodconnects.MajorDim();

	/* set up fConnectivities */
	fConnectivities.Allocate (1);
	fConnectivities[0] = model.ElementGroupPointer(name);
	
	/* set up base class equations array */
	fEqnos.Dimension(1);

	/* set element equation and node lists */
	ConfigureElementData();

	/* print connectivity data */
	WriteConnectivity(out);
}

/***********************************************************************
* Private
***********************************************************************/

/* call AFTER 2 and 3 body node lists are set */
void UnConnectedRodT::ConfigureElementData(void)
{
	/* base class equations arrays */
	const iArray2DT* connects = fConnectivities[0];
	iArray2DT& rod_eqnos = fEqnos[0];

	/* allocate memory */
	int nen = connects->MinorDim();
	int nel = connects->MajorDim();
	fElementCards.Dimension(nel);
	rod_eqnos.Dimension(nel, nen*NumDOF());

	/* set 2 body element data */
	int block_index = 0;
	for (int i = 0; i < nel; i++)	
	{
		/* element card */
		ElementCardT& card = fElementCards[i];
	
		/* node and equation numbers */			
		card.NodesX().Alias(NumElementNodes(), (*connects)(i));
		card.Equations().Alias(rod_eqnos.MinorDim(), rod_eqnos(i));
		
		/* material number */
		card.SetMaterialNumber(fBlockData[block_index].MaterialID());
		
		if (i == fBlockData[block_index].StartNumber() + 
		         fBlockData[block_index].Dimension() - 1)
			block_index++;
	}
}

/* print connectivity element data */
void UnConnectedRodT::PrintConnectivityData(ostream& out)
{
	out << " Number of 2 body interactions . . . . . . . . . = " 
	    << fConnectivities[0]->MajorDim() << '\n';

	/* 2-body connectivities */
	out << "\n Connectivities:\n\n";
	out << setw(kIntWidth) << "no."; 
	for (int i = 1; i <= 2; i++)
	{
		out << setw(kIntWidth - 2) << "n[";
		out << i << "]";
	}
	out << '\n';
	fConnectivities[0]->WriteNumbered(out);
}	
