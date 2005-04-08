/* $Id: ConveyorAdhesiveT.cpp,v 1.1.2.1 2005-04-08 00:50:24 thao Exp $ */
#include "NodeManagerT.h"
#include "FEManagerT.h"
#include "ModelManagerT.h"
#include "LocalArrayT.h"
#include "ElementBaseT.h"
#include "KBC_PrescribedT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"
#include "ConveyorAdhesiveT.h"
#include "ContinuumElementT.h"
#include "CSEAnisoT.h"

using namespace Tahoe;

/* parameters */
const double TipNoise = 1.0e-06;

/* constructor */
ConveyorAdhesiveT::ConveyorAdhesiveT(NodeManagerT& node_manager, FieldT& field):
	ConveyorT(node_manager,field) {}

/* set to initial conditions */
void ConveyorAdhesiveT::InitialCondition(void)
{
	/* inherited */
	ConveyorT::InitialCondition();
}

/* initialize data - called immediately after construction */
void ConveyorAdhesiveT::Initialize(ifstreamT& in)
{
        const char caller[] = "ConveyorAdhesiveT::Initialize";

	/* only 2D for now */
	int nsd = fNodeManager.NumSD();
	if (nsd != 2) ExceptionT::GeneralFail(caller, "only tested for 2D: %d", nsd);
	
	/* prescribed dimensions */
	in >> fMeshRepeatLength;
	in >> fWindowShiftDistance;
	in >> fRightMinSpacing;
	in >> fTrackingInterval;
	if (fMeshRepeatLength < kSmall) ExceptionT::BadInputValue(caller, "%g < 0", fMeshRepeatLength);
	if (fTrackingInterval < 0) ExceptionT::BadInputValue(caller);
	
	/* boundary stretching */
	in >> fULBC_Code;
	in >> fULBC_Value;
	in >> fULBC_ScheduleNumber; fULBC_ScheduleNumber--;
	fULBC_Schedule = fNodeManager.Schedule(fULBC_ScheduleNumber);
	if (!fULBC_Schedule) ExceptionT::BadInputValue(caller, "could not resolve schedule %d", fULBC_ScheduleNumber+1);
	
	/* initial crack tip position */
	in >> fTipElementGroup
	   >> fTipX_0
	   >> fTipY_0
	   >> fTipOutputCode
	   >> fTipColumnNum
	   >> fTipThreshold;
	fTipElementGroup--;
	fTipColumnNum--;
	
	/* damping */
	in >> fDampingWidth
	   >> fDampingCoefficient;
	if (fDampingWidth < 0.0 || fDampingCoefficient < 0.0) ExceptionT::BadInputValue(caller, "improper damping");
	
	/* read lower/upper boundary nodes */
	ArrayT<StringT> id_list;
	ReadNodes(in, id_list, fBottomNodes);
	ReadNodes(in, id_list, fTopNodes);
	
	fKBC_Cards.Dimension(nsd*(fBottomNodes.Length() + fTopNodes.Length()));
//	fKBC_Cards.Dimension(fBottomNodes.Length() + fTopNodes.Length());
	int node = 0;

	for (int i = 0; i < fBottomNodes.Length(); i++) {
		KBC_CardT& card = fKBC_Cards[node++];
		card.SetValues(fBottomNodes[i], 1,KBC_CardT::kFix,0,0);
		card.SetSchedule(fULBC_Schedule);
	}
	for (int i = 0; i < fTopNodes.Length(); i++) {
		KBC_CardT& card = fKBC_Cards[node++];
		card.SetValues(fTopNodes[i], 1, fULBC_Code, fULBC_ScheduleNumber, fULBC_Value);
		card.SetSchedule(fULBC_Schedule);
	}
	
	/* set stretching tangent cards */
	for (int i = 0; i < fBottomNodes.Length(); i++) {
		KBC_CardT& card = fKBC_Cards[node++];
		card.SetValues(fBottomNodes[i], 0, KBC_CardT::kFix, 0, 0);
	}
	for (int i = 0; i < fTopNodes.Length(); i++) {
		KBC_CardT& card = fKBC_Cards[node++];
		card.SetValues(fTopNodes[i], 0, KBC_CardT::kFix, 0, 0);
	} 

     
	/* find boundaries */
	const dArray2DT& init_coords = fNodeManager.InitialCoordinates();
	dArrayT X2(init_coords.MajorDim());
	init_coords.ColumnCopy(0, X2);
	X2.MinMax(fX_Left, fX_Right);
	X2.Free();
	
        /* set the periodic distance */
        fX_PeriodicLength = fX_Right - fX_Left + fMeshRepeatLength;
        fWidthDeadZone = fMeshRepeatLength*1.5;
        if (fWidthDeadZone > fRightMinSpacing) ExceptionT::GeneralFail(caller);

        /* open file for tracking information */
        StringT file;
        file.Root(in.filename());
        file.Append(".tracking");
        fTrackingOutput.open(file);


        /* create controller for the right and left edge of the domain */
        fRightEdge = new KBC_PrescribedT(fNodeManager);
        fField.AddKBCController(fRightEdge);

        /*find the right edge*/
        int nnd = init_coords.MajorDim();
        iAutoArrayT rightnodes(0);
        const double* px = init_coords.Pointer();
        for (int i = 0; i < nnd; i++)
        {
                if (fabs(*px - fX_Right) < kSmall)
                        rightnodes.Append(i);
                px += nsd;
        }

        /*fix the right edge*/
        ArrayT<KBC_CardT>& cards = fRightEdge->KBC_Cards();
        cards.Dimension(rightnodes.Length());
        for (int i=0; i< cards.Length(); i++) {
                KBC_CardT& card = cards[i];
                card.SetValues(rightnodes[i], 0, KBC_CardT::kFix, 0, 0.0);
        }

        fShiftedNodes.Dimension(rightnodes);
        rightnodes.CopyInto(fShiftedNodes);
}

void ConveyorAdhesiveT::WriteParameters(ostream& out) const
{
	/* inherited */
	KBC_ControllerT::WriteParameters(out);
}

void ConveyorAdhesiveT::Reset(void)
{
	/* inherited */
	KBC_ControllerT::Reset();

	/* reset system */
	fTrackingCount--;
	fTrackingPoint = fTrackingPoint_last;
	fX_Left = fX_Left_last;
	fX_Right = fX_Right_last;
}

/* open time interva; */
void ConveyorAdhesiveT::InitStep(void)
{
	/* inherited */
	ConveyorT::InitStep();
}

/* computing residual force */
void ConveyorAdhesiveT::FormRHS(void)
{
	/* inherited */
	ConveyorT::FormRHS();
}

/* apply the update to the solution. Does nothing by default. */
void ConveyorAdhesiveT::Update(const dArrayT& update)
{
	/* inherited */
	ConveyorT::Update(update);
}

/* signal that the solution has been found */
void ConveyorAdhesiveT::CloseStep(void)
{
	/* inherited */
	KBC_ControllerT::CloseStep();

	/* report tracking point */
	if (fTrackingCount == fTrackingInterval) {
	
		/* boundary displacement */
		dArray2DT& u_field = fField[0];

		/* nodes to get reference stretch */
		double uY_top = u_field(fTopNodes[0], 1);
		
		const dArray2DT& initial_coords = fNodeManager.InitialCoordinates();
		for (int i = 0; i < fShiftedNodes.Length(); i++) {
		  if ((fabs(initial_coords(fShiftedNodes[i],1) - fTipY_0) < kSmall) && (u_field(fShiftedNodes[i],1) > kSmall))
		    fuY_interface = u_field(fShiftedNodes[i],1);
		}

		fTrackingOutput << fNodeManager.FEManager().Time() << ' ' 
		                << fTrackingPoint << ' '
						<< fuY_interface << ' '
		                << uY_top <<'\n';
		fTrackingCount = 0;
	}

	/* update history */
	fTrackingPoint_last = fTrackingPoint;
        fX_Left_last = fX_Left;
        fX_Right_last = fX_Right;
}

/* returns true if the internal force has been changed since
 * the last time step */
GlobalT::RelaxCodeT ConveyorAdhesiveT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = KBC_ControllerT::RelaxSystem();

	/* don't track this time */
	if (fTrackingCount != fTrackingInterval) return relax;
	
	/* check to see if focus needs to be reset */
	fTrackingPoint = TrackPoint(kRightMost, fTipThreshold);
	bool initiation = fabs(fTrackingPoint_last - fX_Left) < kSmall;

	if ( (fTrackingPoint - fTrackingPoint_last < fWindowShiftDistance*1.1) || initiation)
	{
		if (SetSystemFocus(fTrackingPoint)) 
		{
			fDampingReset = true;
			relax = GlobalT::MaxPrecedence(relax, GlobalT::kReEQ);
		
			/* message */
			cout << "\n ConveyorT::RelaxSystem: setting system focus = " << fTrackingPoint << endl;
		}
		return relax;
	}
	else 
	{
		cout << "\n Crack growth "<< fTrackingPoint-fTrackingPoint_last << " exceeds window-shift-distance.";
		cout << "\nfTrackingPoint: "<<fTrackingPoint<<"\nfTrackingPoint_last: "<<fTrackingPoint_last;
		relax = GlobalT::MaxPrecedence(relax, GlobalT::kFailReset);
		return relax;
	}
}

void ConveyorAdhesiveT::ReadRestart(ifstreamT& in)
{
	/*inheritde*/
	ConveyorT::ReadRestart(in);	
}

void ConveyorAdhesiveT::WriteRestart(ofstreamT& out) const
{
        /* inherited */
        ConveyorT::WriteRestart(out);
}

/**********************************************************************
 * Protected
 **********************************************************************/

/* reset system to new center */
bool ConveyorAdhesiveT::SetSystemFocus(double focus)
{
	/* no need to shift the window */
	bool no_shift = fX_Right - focus > fRightMinSpacing;
	
	if (no_shift) return false;

	/* shift window */
	fX_Left  += fWindowShiftDistance;
	fX_Right += fWindowShiftDistance;

	/* model information */
	const FEManagerT& fe = fNodeManager.FEManager();
	ModelManagerT* model = fe.ModelManager();

	/* reference coordinates */
	const dArray2DT& initial_coords = fNodeManager.InitialCoordinates();

	/* fields */
	dArray2DT& u_field = fField[0];
	dArray2DT* Du_field = NULL;
	if (fField.Order() > 0) Du_field = &(fField[1]);
	dArray2DT* DDu_field = NULL;
	if (fField.Order() > 1) DDu_field = &(fField[2]);

	/* nodes to get reference stretch */
	double  Y_top = initial_coords(fTopNodes[0], 1);
	double uY_top = u_field(fTopNodes[0], 1);
	double duY_dY = (uY_top-fuY_interface)/(Y_top-fTipY_0);

	/* has damping */
	bool has_damping = (fabs(fDampingWidth) > kSmall && fabs(fDampingCoefficient) > kSmall) ? true : false;

	/* shift reference coordinates and correct fields */
	int nnd = initial_coords.MajorDim();
	int nsd = initial_coords.MinorDim();
	const double* px = initial_coords.Pointer();
	dArrayT new_coords(nsd);
	fShiftedNodes.Dimension(0);
	fDampingNodes.Dimension(0);
	for (int i = 0; i < nnd; i++)
	{
		/* node outside the window */
		if (*px < fX_Left-fMeshRepeatLength/10.0)
		{
			/* store */
			fShiftedNodes.Append(i);
		
			/* shift reference coordinates */
			new_coords[0] = initial_coords(i,0) + fX_PeriodicLength;
			new_coords[1] = initial_coords(i,1);
			model->UpdateNode(new_coords, i);

			/* correct displacements */
			u_field(i,0) = 0.0;
			u_field(i,1) = fuY_interface + duY_dY*(initial_coords(i,1) - fTipY_0); /* interpolate between fixed boundary and interface */
			
			/* zero higher order components */
			if (Du_field) {
				(*Du_field)(i,0) = 0.0;
				(*Du_field)(i,1) = 0.0;
			}
			if (DDu_field) {
				(*DDu_field)(i,0) = 0.0;
				(*DDu_field)(i,1) = 0.0;
			}
		}
		
		/* check for damping */
		if (has_damping)
			if (*px - fX_Left < fDampingWidth || fX_Right - *px < fDampingWidth)
				fDampingNodes.Append(i);
		
		px += nsd;
	}
	fNodeManager.UpdateCurrentCoordinates();
	
	/* reset cards for right edge */
	ArrayT<KBC_CardT>& cards = fRightEdge->KBC_Cards();
	cards.Dimension(fShiftedNodes.Length());
	for (int i = 0; i < cards.Length(); i++) {
		KBC_CardT& card = cards[i];
		card.SetValues(fShiftedNodes[i], 0, KBC_CardT::kFix, 0, 0.0);
	}

	/* mark elements linking left to right edge as inactive */
	MarkElements();

	/* signal change */
	return true;
}

