/* $Id: MixedSWDiamondT.cpp,v 1.5 2004-06-17 07:41:37 paklein Exp $ */
/* created: paklein (03/22/1997) */

#include "MixedSWDiamondT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ifstreamT.h"
#include "FindNeighbor23T.h"
#include "ScheduleT.h"

/* parameters */

using namespace Tahoe;

const int kSWMaxNeighbors0 = 4;

/* constructor */
MixedSWDiamondT::MixedSWDiamondT(const ElementSupportT& support, const FieldT& field):
	SWDiamondT(support, field),
	fCurrMatType(-1)
{
	ElementSupport().Input() >> fLTfNum;
	fLTfPtr = ElementSupport().Schedule(fLTfNum);
}

/*
* Apply pre-conditions at the current time step.  Signal
* all listeners that the time has just been incremented.
*/
void MixedSWDiamondT::InitStep(void)
{
	/* inherited */
	SWDiamondT::InitStep();
	
	/* update variable material properties */
	double x   = fLTfPtr->Value();
	double xm1 = 1.0 - x; //assumes 0 ² x ² 1 during the run!!!!

	/* 3rd material is variable */
	fSWDataList[2].feps = xm1*fSWDataList[0].feps + x*fSWDataList[1].feps;
	fSWDataList[2].flambda = xm1*fSWDataList[0].flambda + x*fSWDataList[1].flambda;
	fSWDataList[2].fa = xm1*fSWDataList[0].fa + x*fSWDataList[1].fa;
}

/***********************************************************************
* Protected
***********************************************************************/

/*
* Print element group data.
*/
void MixedSWDiamondT::PrintControlData(ostream& out) const
{
	/* inherited */
	SWDiamondT::PrintControlData(out);

	out << " Material variation LTf number . . . . . . . . . = " << fLTfNum << '\n';
}

/* element data */
void MixedSWDiamondT::ReadMaterialData(ifstreamT& in)
{
	int numsets;
	
	in >> numsets;	if (numsets < 1) throw ExceptionT::kBadInputValue;

	fSWDataList.Dimension(numsets);

	for (int i = 0; i < numsets; i++)
		fSWDataList[i].Read(in);	
}

void MixedSWDiamondT::WriteMaterialData(ostream& out) const
{
	out << "\n Material Sets:\n";
	out << " Number of material sets . . . . . . . . . . . . = " << fSWDataList.Length() << '\n';

	for (int i = 0; i < fSWDataList.Length(); i++)
		fSWDataList[i].Write(out);
}

void MixedSWDiamondT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	int num_nodes_used;
	in >> num_nodes_used;
	if (num_nodes_used != -1 /* || num_nodes_used > 0 */) throw ExceptionT::kBadInputValue;
								//general case not implemented yet		

	/* neighbor distance assuming material 0*/
	double tolerance = 1.01*pow(2.0,1.0/6.0)*fSWDataList[0].fa;

	/* read nodes used */
	if (num_nodes_used == -1) //use ALL nodes
	{
		/* node type tags */
		EchoNodeTags(in, out);
	
		/* connector */
		FindNeighbor23T Connector(ElementSupport().CurrentCoordinates(), kSWMaxNeighbors0);
	
		/* connect nodes - dimensions lists */
		Connector.GetNeighors(fNodes_2Body, fNodes_3Body, tolerance);
	}
	else //only use specified nodes
	{
		throw ExceptionT::kGeneralFail;
	
		//see MixedSWDiamondT for sample code
	}
	
	/* set element equation and node lists */
	ConfigureElementData();

	/* print connectivity data */
	PrintConnectivityData(out);
}

/* element list increment */
bool MixedSWDiamondT::Next2Body(void)
{
	bool result = SWDiamondT::Next2Body();
	if (result)
	{
		iArrayT& nodes = CurrentElement().NodesX();

		if (fNodeTypes[nodes[0]] == fNodeTypes[nodes[1]] &&
		    fNodeTypes[nodes[1]] != fCurrMatType)
			CopyMaterialData( fNodeTypes[nodes[1]] );
		else if (fNodeTypes[nodes[0]] != fNodeTypes[nodes[1]]) //mix
			Mix2Body(fNodeTypes[nodes[0]], fNodeTypes[nodes[1]]);	
	}

	return result;
}

bool MixedSWDiamondT::Next3Body(void)
{
	bool result = SWDiamondT::Next3Body();
	if (result)
	{
		iArrayT& nodes = CurrentElement().NodesX();

		if (fNodeTypes[nodes[0]] == fNodeTypes[nodes[1]] &&
		    fNodeTypes[nodes[1]] == fNodeTypes[nodes[2]] &&
		    fNodeTypes[nodes[2]] != fCurrMatType)
			CopyMaterialData( fNodeTypes[nodes[2]] );
		else if (fNodeTypes[nodes[0]] != fNodeTypes[nodes[1]] ||
		         fNodeTypes[nodes[1]] != fNodeTypes[nodes[2]]) //mix
			Mix3Body(fNodeTypes[nodes[0]],
			         fNodeTypes[nodes[1]], //must be middle node
			         fNodeTypes[nodes[2]]);	
	}

	return result;
}

/***********************************************************************
* Private
***********************************************************************/

/* echo node type tags */
void MixedSWDiamondT::EchoNodeTags(istream& in, ostream& out)
{
	int numnodes = ElementSupport().NumNodes();
	int numtypes = fSWDataList.Length();

	/* allocate memory */
	fNodeTypes.Dimension(numnodes);

	/* header */
	out << "\n Node type tags : \n";

	/* echo data */
	int linecount = 0;
	for (int i = 0; i < numnodes; i++)
	{
		int nodenum;
		
		in >> nodenum;
		in >> fNodeTypes[nodenum];

		/* checks */
		if (fNodeTypes[nodenum] < 1 ||
		    fNodeTypes[nodenum] > numtypes) throw ExceptionT::kBadInputValue;
	
		out << setw(kIntWidth) << nodenum;
		out << setw(kIntWidth) << fNodeTypes[nodenum];
		
		if (++linecount == 4)
		{
			out << '\n';
			linecount = 0;
		}
		else
			out << "   ";
	}
	
	if (linecount != 0) out << '\n';
	
	/* correct offset */
	fNodeTypes += -1;
}

/* copy material properties from the specified set */
void MixedSWDiamondT::CopyMaterialData(int setnum)
{
	fCurrMatType = setnum;

	SWDataT& currmat = fSWDataList[fCurrMatType];

	/* unit scaling */
	feps   = currmat.feps;

	/* 2 body potential */
	fA     = currmat.fA;
	fB     = currmat.fB;
	fdelta = currmat.fdelta;
	 	
	/* 3 body potential */
	fgamma  = currmat.fgamma;
	flambda = currmat.flambda; 	
	frcut   = currmat.frcut;
	fa      = currmat.fa;	
}

void MixedSWDiamondT::Mix2Body(int m1, int m2)
{
	/* set constant values */
	CopyMaterialData(m1);

	/* mixture rules */
	feps    = sqrt(fSWDataList[m1].feps*fSWDataList[m2].feps);
	flambda = sqrt(fSWDataList[m1].flambda*fSWDataList[m2].flambda);
	fa      = 0.5*(fSWDataList[m1].fa + fSWDataList[m2].fa);

	fCurrMatType = -1;
}

void MixedSWDiamondT::Mix3Body(int m1, int m_mid, int m2)
{
	/* set constant values */
	CopyMaterialData(m1);

	/* mixture rules */
	feps    = sqrt(fSWDataList[m_mid].feps*
	               sqrt(fSWDataList[m1].feps*fSWDataList[m2].feps));
	flambda = sqrt(fSWDataList[m_mid].flambda*
	               sqrt(fSWDataList[m1].flambda*fSWDataList[m2].flambda));
	fa      = 0.5*(fSWDataList[m_mid].fa +
	               (0.5*(fSWDataList[m1].fa +
	                     fSWDataList[m2].fa)));

	fCurrMatType = -1;
}
