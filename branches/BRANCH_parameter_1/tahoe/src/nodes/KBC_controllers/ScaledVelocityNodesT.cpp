/* $Id: ScaledVelocityNodesT.cpp,v 1.1 2003-04-24 20:40:24 cjkimme Exp $ */
#include "ScaledVelocityNodesT.h"
#include "NodeManagerT.h"
#include "ifstreamT.h"
#include "RandomNumberT.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* constructor */
ScaledVelocityNodesT::ScaledVelocityNodesT(NodeManagerT& node_manager, BasicFieldT& field):
	KBC_ControllerT(node_manager),
	fField(field),
	fSteps(0),
	fNodes(),
	fDummySchedule(1.0),
	fRandom(NULL),
	qFirstTime(false)
{
	// nein
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
	ReadNodes(in, fNodeIds, fNodes);
	
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
		if (!qFirstTime)
			qFirstTime = true;
		else
			qFirstTime = false;
			
	/* inherited */
	KBC_ControllerT::InitStep();
	
	if (!qFirstTime)
		fKBC_Cards.Dimension(0);

}


void ScaledVelocityNodesT::InitialCondition(void)
{
	if (!qIConly)
	{
		fT_0 = fTempSchedule->Value();
	}
	
	/* number of scaled nodes */
	int n_scaled = fNodes.Length();
	int ndof = fField.NumDOF();
	fKBC_Cards.Dimension(n_scaled*ndof);
	
	/* workspace to generate velocities */
	dArray2DT velocities(n_scaled, ndof);
	dArrayT vCOM(ndof);
	
	/* generate gaussian dist of random vels */
	fRandom->RandomArray(velocities);
	velocities *= sqrt(fT_0*fkB/fMass);
	
	/* get centre of mass v */
	for (int j = 0; j < ndof; j++)
		vCOM[j] = velocities.ColumnSum(j)/n_scaled;
	/* subtract it off and calculate resultant total kinetic energy */
	double tKE = 0.;
	for (int i = 0; i < n_scaled; i++)
		for (int j = 0; j < ndof; j++)
		{	
			velocities(i,j) -= vCOM[j];
			tKE += velocities(i,j)*velocities(i,j);
		}
	tKE *= fMass/n_scaled/ndof/fkB;
	
	velocities *= sqrt(fT_0/tKE);
	
	/* grab the velocities */
//	const dArray2DT* velocities = NULL;
//	if (fField.Order() > 0)
//		velocities = &fField[1];		
//	if (!velocities)
//		ExceptionT::GeneralFail("ScaledVelocityNodesT::InitialCondition","Cannot get velocity field ");
	
	/* modify them */
//	dArray2DT* nonconvs = const_cast<dArray2DT*>(velocities);
	
	/* generate BC cards */
	KBC_CardT* pcard = fKBC_Cards.Pointer();
	for (int i = 0; i < n_scaled; i++)
	{
		double* v_j = velocities.Pointer();	
			
	    for (int j = 0; j < ndof; j++)
		{	
			/* set values */
			pcard->SetValues(fNodes[i], j, KBC_CardT::kVel, 0, *v_j++);

			/* dummy schedule */
			pcard->SetSchedule(&fDummySchedule);
			pcard++;
		} 
	}

}

/**********************************************************************
 * Protected
 **********************************************************************/

/* initialize the current step */
void ScaledVelocityNodesT::SetBCCards(void)
{
	/* number of scaled nodes */
	int n_scaled = fNodes.Length();
	int ndof = fField.NumDOF();
	fKBC_Cards.Dimension(n_scaled*ndof);

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
	if (n_scaled > 0)
	{
		vCOM = 0.;
		//double totalMass = fMass*n_scaled;
		for (int i = 0; i < n_scaled; i++)
		{
			double* v_j = (*velocities)(fNodes[i]);
		
	    	for (int j = 0; j < ndof; j++)
			{	
				tKE += (*v_j)*(*v_j);
				vCOM[j] += *v_j++;
			} 
		}
		
		/* adjust KE to COM frame */
		for (int j = 0; j < n_scaled; j++)
			tKE -= n_scaled*vCOM[j]*vCOM[j];
		tKE *= fMass;
		
		/* want new KE to be ndof/2*n_scaled * kT */
		vscale = sqrt(ndof*n_scaled*fkB*fTempScale*(fTempSchedule->Value())/tKE);
		
		/* generate BC cards */
		KBC_CardT* pcard = fKBC_Cards.Pointer();
		for (int i = 0; i < n_scaled; i++)
		{
			double* v_j = (*velocities)(fNodes[i]);	
			
	    	for (int j = 0; j < ndof; j++)
			{	
				/* set values */
				pcard->SetValues(fNodes[i], j, KBC_CardT::kVel, 0, (*v_j++-vCOM[j])*vscale);

				/* dummy schedule */
				pcard->SetSchedule(&fDummySchedule);
				pcard++;
			} 
		}
	}
}

