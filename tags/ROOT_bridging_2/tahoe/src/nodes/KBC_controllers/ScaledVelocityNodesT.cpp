/* $Id: ScaledVelocityNodesT.cpp,v 1.4 2003-05-06 19:59:44 cjkimme Exp $ */
#include "ScaledVelocityNodesT.h"
#include "NodeManagerT.h"
#include "FEManagerT.h"
#include "ModelManagerT.h"
#include "ifstreamT.h"
#include "RandomNumberT.h"
#include "iArrayT.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* constructor */
ScaledVelocityNodesT::ScaledVelocityNodesT(NodeManagerT& node_manager, BasicFieldT& field):
	KBC_ControllerT(node_manager),
	fField(field),
	fNodes(),
	fDummySchedule(1.0),
	fRandom(NULL),
	qFirstTime(false),
	qAllNodes(false),
	fIncs(0),
	fIncCt(0)
{
	// nein
}

ScaledVelocityNodesT::~ScaledVelocityNodesT(void)
{
   if (fRandom)
      delete fRandom;
}

void ScaledVelocityNodesT::WriteParameters(ostream& out) const
{
	/* inherited */
	KBC_ControllerT::WriteParameters(out);

}

void ScaledVelocityNodesT::Initialize(ifstreamT& in)
{

#pragma message("Add rescaled region")
	//int nodesOrRegion; 

	/* Get node sets */
	char nextun = in.next_char();
	int allOrSome;
	if (nextun == '-')
		in >> allOrSome;
	else
		allOrSome = atoi(&nextun);
	if (!allOrSome)
	{   // use all the nodes
		in >> allOrSome; // fast-forward the stream
		qAllNodes = true;
		fNodes.Dimension(fNodeManager.NumNodes());
		fNodes.SetValueToPosition();
	}
	else	
		if (allOrSome > 0)
		{   // Read in Node Sets and use them
			ReadNodes(in, fNodeIds, fNodes);
		}
		else
		{   // Read in Node Sets and use all nodes but theirs
			ArrayT<StringT> notTheseSets;
			int numNotSets = -allOrSome;
			numNotSets = abs(numNotSets);
			notTheseSets.Dimension(numNotSets);
			ModelManagerT* model = fNodeManager.FEManager().ModelManager();
			for (int i=0; i < numNotSets; i++)
			{
	  			StringT& name = notTheseSets[i];
	  			in >> name;
	  			int index = model->NodeSetIndex(name);
	  			if (index < 0) {
	  				cout << "\n ScaledVelocityT::Initialize: error retrieving node set " << name << endl;
	  				throw ExceptionT::kDatabaseFail;
	  			}
			}
			// get all the nodes we don't want
			iArrayT fNotNodes;
			model->ManyNodeSets(notTheseSets, fNotNodes);
			// get all the nodes
			fNodes.Dimension(model->NumNodes());
			fNodes.SetValueToPosition();
			// Take the complement
			for (int i = 0; i < fNotNodes.Length(); i++)
				fNodes[fNotNodes[i]] = -fNodes[fNotNodes[i]] - 1; // if the node is to be deleted, make it < 0
			fNodes.SortDescending();
			fNodes.Resize(fNodes.Length() - fNotNodes.Length());
		}
	
	in >> fMass; 
	if (fMass < kSmall) ExceptionT::BadInputValue("ScaledVelocityNodesT::Initialize","mass must be positive");

	int BCnotIC;
	in >> BCnotIC;
	
	if (BCnotIC)
	{
		in >> fnumTempSchedule >> fTempScale; fnumTempSchedule--;
		fTempSchedule = fNodeManager.Schedule(fnumTempSchedule);	
		if (!fTempSchedule) throw ExceptionT::kBadInputValue;
	
		in >> fIncs;
		qIConly = false;
	}
	else
	{
		qIConly = true;
		in >> fT_0;
	}

	// parameters for initial velocity distribution random num gen
	int rseed;
	in >> rseed;
	
	fRandom = new RandomNumberT(RandomNumberT::kParadynGaussian);
	if (!fRandom)
		ExceptionT::GeneralFail("ScaledVelocityNodesT::Initialize","Cannot create random number gen");
	fRandom->sRand(rseed);

}

void ScaledVelocityNodesT::InitStep(void)
{
	/* really bad, this */
	if (qIConly)
	{
		if (!qFirstTime)
			qFirstTime = true;
		else
			qFirstTime = false;
	}
	else
	{
		fIncCt++;
		if (fIncCt == fIncs) // time to rescale velocities
		{
			fIncCt = 0;
			SetBCCards();
		}
	}
					
	/* inherited */
	KBC_ControllerT::InitStep();
	
	if (qIConly && !qFirstTime)
		fKBC_Cards.Dimension(0);

}


void ScaledVelocityNodesT::InitialCondition(void)
{
	if (!qIConly)
	{
		fT_0 = fTempScale*fTempSchedule->Value();
	}
	
	/* get processor number */
	int np = fNodeManager.Rank();
	const ArrayT<int>* pMap = fNodeManager.ProcessorMap();
	
	/* number of scaled nodes */
	int n_scaled = 0;
	int ndof = fField.NumDOF();
	
	/* workspace to generate velocities */
	dArray2DT velocities(fNodeManager.NumNodes(), ndof); 
	dArrayT vCOM(ndof);
	
	/* generate gaussian dist of random vels */
	fRandom->RandomArray(velocities);
	
	double tKE = 0.;
	vCOM = 0.;
	//double totalMass = fMass*n_scaled;
	for (int i = 0; i < fNodes.Length(); i++)
	{
		int node_i = fNodes[i];
		if (!pMap || (*pMap)[node_i] == np)
		{
			double* v_i = velocities(node_i);
		
	    	for (int j = 0; j < ndof; j++)
			{	
				tKE += (*v_i)*(*v_i);
				vCOM[j] += *v_i++;
			} 
			
			n_scaled++;
		}
	}
	vCOM /= n_scaled;
	
	for (int i = 0; i < fNodes.Length(); i++)
		for (int j = 0; j < ndof; j++)
			if (!pMap || (*pMap)[fNodes[i]] == np)
				velocities(fNodes[i],j) -= vCOM[j];
	
	/* adjust KE to COM frame  and convert it to a temperature */
	for (int j = 0; j < ndof; j++)
		tKE -= n_scaled*vCOM[j]*vCOM[j];
	tKE *= fMass/n_scaled/ndof/fkB;
	
	/* set the temperature */
	velocities *= sqrt(fT_0/tKE);

	/* generate BC cards */
	fKBC_Cards.Dimension(n_scaled*ndof);
	KBC_CardT* pcard = fKBC_Cards.Pointer();
	for (int i = 0; i < fNodes.Length(); i++)
	{
		int node_i = fNodes[i];
		
		if (!pMap || (*pMap)[node_i] == np)
		{
			double* v_i = velocities(node_i);	
			
	    	for (int j = 0; j < ndof; j++)
			{	
				/* set values */
				pcard->SetValues(node_i, j, KBC_CardT::kVel, 0, *v_i++);

				/* dummy schedule */
				pcard->SetSchedule(&fDummySchedule);
				pcard++;
			}
		} 
	}

}

/**********************************************************************
 * Protected
 **********************************************************************/

/* initialize the current step */
void ScaledVelocityNodesT::SetBCCards(void)
{
	/* get processor number */
	int np = fNodeManager.Rank();
	const ArrayT<int>* pMap = fNodeManager.ProcessorMap();

	/* number of scaled nodes */
	int n_scaled = 0; 
	int ndof = fField.NumDOF();

	/* grab the velocities */
	const dArray2DT* velocities = NULL;
	if (fField.Order() > 0)
		velocities = &fField[1];		
	if (!velocities)
		ExceptionT::GeneralFail("ScaledVelocityNodesT::SetBCCards","Cannot get velocity field ");
		
	/* 	assume uniform mass for now */

	/* calculate CM velocity and temperature */
	dArrayT vCOM(ndof);
	double tKE = 0.; // total kinetic energy
	double vscale; // scale velocity by this after subtracting off vCOM
	if (fNodes.Length() > 0)
	{
		vCOM = 0.;
		//double totalMass = fMass*n_scaled;
		for (int i = 0; i < n_scaled; i++)
		{
			int node_i = fNodes[i];
			
			if (!pMap || (*pMap)[node_i] == np)
			{
				double* v_i = (*velocities)(node_i);
		
	    		for (int j = 0; j < ndof; j++)
				{	
					tKE += (*v_i)*(*v_i);
					vCOM[j] += *v_i++;
				} 
				
				n_scaled++;
			}
		}
		vCOM /= n_scaled;
		
		/* adjust KE to COM frame */
		for (int j = 0; j < ndof; j++)
			tKE -= n_scaled*vCOM[j]*vCOM[j];
		tKE *= fMass;
		
		/* want new KE to be ndof/2*n_scaled * kT */
		vscale = sqrt(ndof*n_scaled*fkB*fTempScale*(fTempSchedule->Value())/tKE);
		
		/* generate BC cards */
		fKBC_Cards.Dimension(n_scaled*ndof);
		KBC_CardT* pcard = fKBC_Cards.Pointer();
		for (int i = 0; i < n_scaled; i++)
		{
			int node_i = fNodes[i];
			if (!pMap || (*pMap)[node_i] == np)
			{
				double* v_i = (*velocities)(node_i);	
				
		    	for (int j = 0; j < ndof; j++)
				{	
					/* set values */
					pcard->SetValues(node_i, j, KBC_CardT::kVel, 0, (*v_i++-vCOM[j])*vscale);

					/* dummy schedule */
					pcard->SetSchedule(&fDummySchedule);
					pcard++;
				}
			} 
		}
	}
}

