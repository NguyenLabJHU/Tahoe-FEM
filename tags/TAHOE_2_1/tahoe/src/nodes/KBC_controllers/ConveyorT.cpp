/* $Id: ConveyorT.cpp,v 1.5 2004-07-22 08:42:46 paklein Exp $ */
#include "ConveyorT.h"
#include "NodeManagerT.h"
#include "FEManagerT.h"
#include "ModelManagerT.h"
#include "LocalArrayT.h"
#include "ElementBaseT.h"
#include "KBC_ControllerT.h"

using namespace Tahoe;

/* parameters */
const double TipNoise = 1.0e-06;

/* constructor */
ConveyorT::ConveyorT(const BasicSupportT& support, FieldT& field):
	KBC_ControllerT(support),
	fField(field),
	fULBC_Value(0.0),
	fULBC_Code(KBC_CardT::kFix),
	fULBC_ScheduleNumber(-1),
	fULBC_Schedule(NULL),
	fTrackingInterval(-1),
	fRightEdge(NULL),
	fDampingWidth(0.0),
	fDampingCoefficient(0.0),
	fDampingReset(true)
{
	
}

/* set to initial conditions */
void ConveyorT::InitialCondition(void)
{
	/* reset */
	fTrackingCount = 0;
	fTrackingPoint = fX_Left;

	/* mark elements linking left to right edge as inactive */
	MarkElements();	
}

#if 0
/* initialize data - called immediately after construction */
void ConveyorT::Initialize(ifstreamT& in)
{
	const char caller[] = "ConveyorT::Initialize";

	/* only 2D for now */
	int nsd = fSupport.NumSD();
	if (nsd != 2) ExceptionT::GeneralFail(caller, "only tested for 2D: %d", nsd);

	/* prescribed dimensions */
	in >> fMeshRepeatLength;
	in >> fWindowShiftDistance;
	in >> fRightMinSpacing;
	in >> fTrackingInterval;
	if (fMeshRepeatLength < kSmall) ExceptionT::BadInputValue(caller, "%g < 0", fMeshRepeatLength);
	if (fTrackingInterval < 0) ExceptionT::BadInputValue(caller);

	/* boundary stretching */
	int schedule_number = -1;
	in >> fULBC_Code;
	in >> fULBC_Value;
	in >> fULBC_ScheduleNumber; fULBC_ScheduleNumber--;
	fULBC_Schedule = fSupport.Schedule(fULBC_ScheduleNumber);	
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

	/* set stretching BC cards */
	fKBC_Cards.Dimension(nsd*(fBottomNodes.Length() + fTopNodes.Length()));
	int node = 0;
	double valueby2 = fULBC_Value/2.0;
	for (int i = 0; i < fBottomNodes.Length(); i++) {
		KBC_CardT& card = fKBC_Cards[node++];
		card.SetValues(fBottomNodes[i], 1, fULBC_Code, fULBC_Schedule, -valueby2);
	}
	for (int i = 0; i < fTopNodes.Length(); i++) {
		KBC_CardT& card = fKBC_Cards[node++];
		card.SetValues(fTopNodes[i], 1, fULBC_Code, fULBC_Schedule, valueby2);
		}

	/* set stretching tangent cards */
	for (int i = 0; i < fBottomNodes.Length(); i++) {
		KBC_CardT& card = fKBC_Cards[node++];
		card.SetValues(fBottomNodes[i], 0, KBC_CardT::kFix, NULL, 0);
	}
	for (int i = 0; i < fTopNodes.Length(); i++) {
		KBC_CardT& card = fKBC_Cards[node++];
		card.SetValues(fTopNodes[i], 0, KBC_CardT::kFix, NULL, 0);
	}

	/* find boundaries */
	const dArray2DT& init_coords = fSupport.InitialCoordinates();
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
	
	/* create controller for the right edge of the domain */
	fRightEdge = new KBC_ControllerT(fSupport);
	fField.AddKBCController(fRightEdge);
}
#endif

void ConveyorT::Reset(void)
{
	/* inherited */
	KBC_ControllerT::Reset();

	/* reset system */
	fTrackingCount--;
	fTrackingPoint = fTrackingPoint_last;
}

/* open time interva; */
void ConveyorT::InitStep(void)
{
	/* inherited */
	KBC_ControllerT::InitStep();
	fTrackingCount++;
	
	/* create pre-crack */
	if (fSupport.FEManager().Time() < kSmall)
		CreatePrecrack();
}

/* computing residual force */
void ConveyorT::FormRHS(void)
{
	/* inherited */
	KBC_ControllerT::FormRHS();

	/* apply damping */
	if (!fDampingReset && fDampingNodes.Length()) {
	
		/* collect nodal velocities */
		fDampingForce.RowCollect(fDampingNodes, fField[1]);
	
		/* convert velocity to force */
		for (int i = 0; i < fDampingForce.Length(); i++)
			fDampingForce[i] = -fDampingForce[i]*fDampingCoeff[i];

		/* assemble */
		const FEManagerT& fe_man = fSupport.FEManager();
		fe_man.AssembleRHS(fField.Group(), fDampingForce, fDampingEqnos);
	}
}

/* apply the update to the solution. Does nothing by default. */
void ConveyorT::Update(const dArrayT& update)
{
	/* inherited */
	KBC_ControllerT::Update(update);

	/* reset damping force coefficients */
	if (fDampingReset) {
	
		/* field dimension */
		int ndof = fField.NumDOF();
		int ndn = fDampingNodes.Length();
		
		/* dimension workspace */
		fDampingForce.Dimension(ndn, ndof);
		fDampingCoeff.Dimension(ndn, ndof);
		fDampingEqnos.Dimension(ndn, ndof);
	
		/* collect LHS components for damped nodes */
		if (ndn > 0) {
			/* collect equation numbers */
			iArrayT tags;
			tags.Alias(fDampingNodes);
			fField.SetLocalEqnos(tags, fDampingEqnos);

			/* collect nodal masses */
			const FEManagerT& fe = fSupport.FEManager();
			dArrayT diagonals;
			diagonals.Alias(fDampingCoeff);
			fe.DisassembleLHSDiagonal(fField.Group(), diagonals, fDampingEqnos);

			/* convert to mass proportionate damping coefficient */
			for (int i = 0; i < fDampingEqnos.Length(); i++)
				if (fDampingEqnos[i] > 0 && fabs(fDampingCoeff[i]) > kSmall)
					fDampingCoeff[i] = fDampingCoefficient/fDampingCoeff[i];
				else
					fDampingCoeff[i] = 0.0;
		}
		
		/* reset flag */
		fDampingReset = false;
	}
}

/* signal that the solution has been found */
void ConveyorT::CloseStep(void)
{
	/* inherited */
	KBC_ControllerT::CloseStep();

	/* report tracking point */
	if (fTrackingCount == fTrackingInterval) {
	
		/* boundary displacement */
		dArray2DT& u_field = fField[0];

		/* nodes to get reference stretch */
		double uY_bottom = u_field(fBottomNodes[0], 1);
		double uY_top = u_field(fTopNodes[0], 1);
		
		fTrackingOutput << fSupport.Time() << ' ' 
		                << fTrackingPoint << ' '
		                << uY_top - uY_bottom <<'\n';
		fTrackingCount = 0;
	}

	/* update history */
	fTrackingPoint_last = fTrackingPoint;
}

/* returns true if the internal force has been changed since
 * the last time step */
GlobalT::RelaxCodeT ConveyorT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = KBC_ControllerT::RelaxSystem();

	/* don't track this time */
	if (fTrackingCount != fTrackingInterval) return relax;
	
	/* check to see if focus needs to be reset */
	fTrackingPoint = TrackPoint(kRightMost, fTipThreshold);
	if (SetSystemFocus(fTrackingPoint)) 
	{
		fDampingReset = true;
		relax = GlobalT::MaxPrecedence(relax, GlobalT::kReEQ);
		
		/* message */
		cout << "\n ConveyorT::RelaxSystem: setting system focus = " << fTrackingPoint << endl;
	}
	return relax;
}

/**********************************************************************
 * Protected
 **********************************************************************/

/* locate new tracking point */
double ConveyorT::TrackPoint(TrackingTypeT tracking_type, double threshold)
{
	const char caller[] = "ConveyorT::TrackPoint";

	/* near tip element group */
	const FEManagerT& fe_man = fSupport.FEManager();
	ElementBaseT* neartip_group = fe_man.ElementGroup(fTipElementGroup);
	if (!neartip_group)
		ExceptionT::GeneralFail(caller, "could not resolve near tip element group number %d", fTipElementGroup+1);

	/* signal to accumulate nodal values */
	neartip_group->SendOutput(fTipOutputCode);

	/* nodes */
	NodeManagerT& node_manager = fSupport.NodeManager();

	/* tracking type */
	if (tracking_type == kMax)
	{
		/* find the node with max opening stress */
		int maxrow;
		double maxval;
		node_manager.MaxInColumn(fTipColumnNum, maxrow, maxval);
		if (maxrow == -1) ExceptionT::GeneralFail(caller);

		/* tracking point x-coordinate */
		if (maxval > threshold)
			return (fSupport.InitialCoordinates())(maxrow,0);
		else /* default to left edge of domain */
			return fX_Left;
	}
	else if (tracking_type == kRightMost)
	{
		/* reference coordinates */
		const dArray2DT& ref_coords = node_manager.InitialCoordinates();
	
		/* find right most value greater than threshold */
		double right_most = fX_Left;
		dArrayT values;
		node_manager.TopAverage();
		int node = node_manager.NextAverageRow(values);
		while (node != -1)
		{
			if (values[fTipColumnNum] > threshold && ref_coords(node,0) > right_most)
				right_most = ref_coords(node,0);
		
			node = node_manager.NextAverageRow(values);
		}

		return right_most;
	}
	else 
		ExceptionT::GeneralFail(caller, "unknown tracking type %d", tracking_type);
	
	return fX_Left;
}

/* reset system to new center */
bool ConveyorT::SetSystemFocus(double focus)
{
	/* no need to shift the window */
	if (fX_Right - focus > fRightMinSpacing) return false;

	/* shift window */
	fX_Left  += fWindowShiftDistance;
	fX_Right += fWindowShiftDistance;

	/* model information */
	ModelManagerT& model = fSupport.ModelManager();

	/* reference coordinates */
	const dArray2DT& initial_coords = fSupport.InitialCoordinates();

	/* fields */
	dArray2DT& u_field = fField[0];
	dArray2DT* Du_field = NULL;
	if (fField.Order() > 0) Du_field = &(fField[1]);
	dArray2DT* DDu_field = NULL;
	if (fField.Order() > 1) DDu_field = &(fField[2]);

	/* nodes to get reference stretch */
	double  Y_bottom = initial_coords(fBottomNodes[0], 1);
	double uY_bottom = u_field(fBottomNodes[0], 1);
	double  Y_top = initial_coords(fTopNodes[0], 1);
	double uY_top = u_field(fTopNodes[0], 1);
	double duY_dY = (uY_top - uY_bottom)/(Y_top - Y_bottom);

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
		if (*px < fX_Left)
		{
			/* store */
			fShiftedNodes.Append(i);
		
			/* shift reference coordinates */
			new_coords[0] = initial_coords(i,0) + fX_PeriodicLength;
			new_coords[1] = initial_coords(i,1);
			model.UpdateNode(new_coords, i);

			/* correct displacements */
			u_field(i,0) = 0.0;
			u_field(i,1) = uY_bottom + duY_dY*(initial_coords(i,1) - Y_bottom); /* interpolate between lower and upper boundary */
			
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
	fSupport.NodeManager().UpdateCurrentCoordinates();
	
	/* reset cards for right edge */
	ArrayT<KBC_CardT>& cards = fRightEdge->KBC_Cards();
	cards.Dimension(fShiftedNodes.Length());
	for (int i = 0; i < cards.Length(); i++) {
		KBC_CardT& card = cards[i];
		card.SetValues(fShiftedNodes[i], 0, KBC_CardT::kFix, NULL, 0.0);
	}

	/* mark elements linking left to right edge as inactive */
	MarkElements();

	/* signal change */
	return true;
}

/* mark elements linking left to right edge as inactive */
void ConveyorT::MarkElements(void)
{
	/* system information */
	const FEManagerT& fe = fSupport.FEManager();
	const dArray2DT& current_coords = fSupport.CurrentCoordinates();

	/* zone to activate/deactivate elements */
	double X_R_on  = fX_Right - fRightMinSpacing;

	/* mark dead elements (this group only?) */
	int num_element_groups = fe.NumElementGroups();
	for (int i = 0; i < num_element_groups; i++)
	{
		/* element group */
		ElementBaseT* element_group = fe.ElementGroup(i);
		int nel = element_group->NumElements();
		int nen = element_group->NumElementNodes();

		/* local coordinate array */
		LocalArrayT curr_coords(LocalArrayT::kCurrCoords, nen, 2);
		curr_coords.SetGlobal(current_coords);
	
		//TEMP
		if (!element_group->InGroup(fField.Group())) ExceptionT::GeneralFail();
	
		ArrayT<ElementBaseT::StatusT> status(nel);
		for (int j = 0; j < nel; j++)
		{
			/* element info */
			ElementCardT& card = element_group->ElementCard(j);
		
			/* copy status */
			status[j] = (ElementBaseT::StatusT) card.Flag();
		
			/* collect local coordinates */
			curr_coords.SetLocal(card.NodesX());
			
			/* set flag based on x-location and size */
			const double* px = curr_coords(0);
			double x_min = *px;
			double x_max = *px;
			px++;
			for (int k = 1; k < nen; k++) {
				/* bounds */
				x_min = (*px < x_min) ? *px : x_min;
				x_max = (*px > x_max) ? *px : x_max;
				px++;
			}
			
			/* goes end to end */
			if (x_max - x_min > fX_PeriodicLength/2.0)
				status[j] = ElementBaseT::kOFF;
			else if (x_min > X_R_on) /* in reactivation zone */
				status[j] = ElementBaseT::kON;
		}
		
		/* reset element status */
		element_group->SetStatus(status);
	}
}

/* deactivate elements to create a pre-crack */
void ConveyorT::CreatePrecrack(void)
{
	/* element group */
	const FEManagerT& fe = fSupport.FEManager();
	const dArray2DT& current_coords = fSupport.CurrentCoordinates();
	ElementBaseT* element_group = fe.ElementGroup(fTipElementGroup);
	if (!element_group) ExceptionT::GeneralFail();
	int nel = element_group->NumElements();
	int nen = element_group->NumElementNodes();

	/* local coordinate array */
	LocalArrayT curr_coords(LocalArrayT::kCurrCoords, nen, 2);
	curr_coords.SetGlobal(current_coords);
	
	//TEMP
	if (!element_group->InGroup(fField.Group())) ExceptionT::GeneralFail();
	
	ArrayT<ElementBaseT::StatusT> status(nel);
	for (int j = 0; j < nel; j++)
	{
		/* element info */
		ElementCardT& card = element_group->ElementCard(j);
		
		/* copy status */
		status[j] = (ElementBaseT::StatusT) card.Flag();
		
		/* collect local coordinates */
		curr_coords.SetLocal(card.NodesX());
			
		/* element lies "behind" the initial tip position */
		bool X_check_OK = true;
		bool Y_check_OK = true;
		const double* px = curr_coords(0);
		const double* py = curr_coords(1);
		for (int k = 0; X_check_OK && Y_check_OK && k < nen; k++) {

			/* "behind" tip */
			X_check_OK = *px++ < fTipX_0;
			Y_check_OK = fabs(*py++ - fTipY_0) < kSmall;
		}
			
		/* goes end to end */
		if (X_check_OK && Y_check_OK) status[j] = ElementBaseT::kOFF;
	}
		
	/* reset element status */
	element_group->SetStatus(status);

}
