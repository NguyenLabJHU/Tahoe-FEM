/* $Id: ConveyorSymT.cpp,v 1.1 2004-12-29 16:00:58 thao Exp $ */
#include "ConveyorSymT.h"
#include "NodeManagerT.h"
#include "FEManagerT.h"
#include "ModelManagerT.h"
#include "LocalArrayT.h"
#include "ElementBaseT.h"
#include "KBC_ControllerT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "ifstreamT.h"
#include "ofstreamT.h"

using namespace Tahoe;

/* constructor */
ConveyorSymT::ConveyorSymT(const BasicSupportT& support, FieldT& field):
  ConveyorT(support, field)
{
  SetName("symmetric_conveyor");
}

/* set to initial conditions */
void ConveyorSymT::InitialCondition(void)
{
  /* inherited */
  ConveyorT::InitialCondition();
}

void ConveyorSymT::Reset(void)
{
  /* inherited */
  ConveyorT::Reset();
}

/* open time interva; */
void ConveyorSymT::InitStep(void)
{
  /* inherited */
  ConveyorT::InitStep();
}

/* computing residual force */
void ConveyorSymT::FormRHS(void)
{
  /* inherited */
  ConveyorT::FormRHS();
}

/* apply the update to the solution. Does nothing by default. */
void ConveyorSymT::Update(const dArrayT& update)
{
  /* inherited */
  ConveyorT::Update(update);
}

/* signal that the solution has been found */
void ConveyorSymT::CloseStep(void)
{
  /* inherited */
  KBC_ControllerT::CloseStep();

  /* report tracking point */
  if (fTrackingCount == fTrackingInterval) {

    /* boundary displacement */
    dArray2DT& u_field = fField[0];

    /* nodes to get reference stretch */
    double uY_top = u_field(fTopNodes[0], 1);

    const dArray2DT& initial_coords = fSupport.InitialCoordinates();
    for (int i = 0; i < fShiftedNodes.Length(); i++) {
      if (fabs(initial_coords(fShiftedNodes[i],1) - fTipY_0) < kSmall)
	fuY_interface = u_field(fShiftedNodes[i],1);
    }

    fTrackingOutput << fSupport.FEManager().Time() << ' '
		    << fTrackingPoint << ' '
		    << fuY_interface << ' '
		    << 2.0*uY_top <<'\n';
    fTrackingCount = 0;
  
  }

  /* update history */
  fX_Left_last = fX_Left;
  fX_Right_last = fX_Right;
  fTrackingPoint_last = fTrackingPoint;
}

/* returns true if the internal force has been changed since
 * the last time step */
GlobalT::RelaxCodeT ConveyorSymT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = KBC_ControllerT::RelaxSystem();

	/* don't track this time */
	if (fTrackingCount != fTrackingInterval) return relax;
	
	/* check to see if focus needs to be reset */
	fTrackingPoint = TrackPoint(fTrackingType, fTipThreshold);
	if (SetSystemFocus(fTrackingPoint)) 
	{
		fDampingReset = true;
		relax = GlobalT::MaxPrecedence(relax, GlobalT::kReEQ);
		
		/* message */
		cout << "\n ConveyorT::RelaxSystem: {time, focus} = {" 
		     << fSupport.Time() << ", " << fTrackingPoint << "}" << endl;
	}
	return relax;
}

void ConveyorSymT::ReadRestart(ifstreamT& in)
{
	/* inherited */
	ConveyorT::ReadRestart(in);
}

void ConveyorSymT::WriteRestart(ofstreamT& out) const
{
  /* inherited */
  ConveyorT::WriteRestart(out);
}

/* describe the parameters needed by the interface */
void ConveyorSymT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ConveyorT::DefineParameters(list); 
}

/* information about subordinate parameter lists */
void ConveyorSymT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	KBC_ControllerT::DefineSubs(sub_list);

	/* tip trackng information */
	sub_list.AddSub("focus_tracking_method", ParameterListT::Once, true);

	/* KBC for upper surfaces */
	sub_list.AddSub("upper_kinematic_BC");

	/* initial tip coordinates */
	sub_list.AddSub("initial_focus_coordinates");}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ConveyorSymT::NewSub(const StringT& name) const
{
	const char caller[] = "ConveyorSymT::NewSub";

	if (name == "focus_tracking_method")
	{
		ParameterContainerT* method = new ParameterContainerT(name);
		method->SetListOrder(ParameterListT::Choice);
	
		/* right most position above threshold value */
		ParameterContainerT right_most("focus_right_most");
		right_most.AddParameter(ParameterT::Double, "threshold");
		method->AddSub(right_most);

		/* left most position above threshold value */
		ParameterContainerT left_most("focus_left_most");
		left_most.AddParameter(ParameterT::Double, "threshold");
		method->AddSub(left_most);
		
		/* at maximum (above threshold) */
		ParameterContainerT max("focus_at_maximum");
		max.AddParameter(ParameterT::Double, "threshold");
		method->AddSub(max);	
	
		return method;
	}
	else if (name == "initial_focus_coordinates")
		return new VectorParameterT(name, 2, 'x');
	else if (name == "upper_kinematic_BC")
	{
		ParameterContainerT* kbc = new ParameterContainerT(name);
		
		ParameterT BC_type(ParameterT::Enumeration, "type");
		BC_type.AddEnumeration("u", 0);
		BC_type.AddEnumeration("D_u", 1);
		BC_type.AddEnumeration("DD_u", 2);
		BC_type.AddEnumeration("D3_u", 3);
		BC_type.AddEnumeration("D4_u", 4);
		BC_type.SetDefault(0);
		kbc->AddParameter(BC_type);
		ParameterT schedule(ParameterT::Integer, "schedule");
		schedule.SetDefault(0);
		kbc->AddParameter(schedule);
		ParameterT value(ParameterT::Double, "value");
		value.SetDefault(0.0);
		kbc->AddParameter(value);

		/* the nodes */
		kbc->AddSub("upper_ID_list");

		return kbc;
	}
	else /* inherited */
		return KBC_ControllerT::NewSub(name);}

/* accept parameter list */
void ConveyorSymT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "ConveyorSymT::TakeParameterList";

	/* only 2D for now */
	int nsd = fSupport.NumSD();
	if (nsd != 2) ExceptionT::GeneralFail(caller, "only tested for 2D: %d", nsd);

	/* inherited */
	KBC_ControllerT::TakeParameterList(list);

	/* dimensions */
	fMeshRepeatLength = list.GetParameter("mesh_repeat_length");
	fWindowShiftDistance = list.GetParameter("window_shift");
	fRightMinSpacing = list.GetParameter("min_right_space");

	/* boundary stretching */
	const ParameterListT& ul_kbc = list.GetList("upper_kinematic_BC");
	int i_code = ul_kbc.GetParameter("type");
	fULBC_Code = KBC_CardT::int2CodeT(i_code + 1);
	fULBC_Value = ul_kbc.GetParameter("value");
	fULBC_ScheduleNumber = ul_kbc.GetParameter("schedule"); fULBC_ScheduleNumber--;
	fULBC_Schedule = fSupport.Schedule(fULBC_ScheduleNumber);	
	if (!fULBC_Schedule) ExceptionT::BadInputValue(caller, "could not resolve schedule %d", fULBC_ScheduleNumber+1);
	
	ArrayT<StringT> id_list;
	const ParameterListT& upper_nodes = ul_kbc.GetList("upper_ID_list");
	StringListT::Extract(upper_nodes,  id_list);	
	GetNodes(id_list, fTopNodes);

	/* tracking */
	fTrackingInterval = list.GetParameter("focus_tracking_increment");
	const ParameterListT& init_focus = list.GetList("initial_focus_coordinates");
	fTipX_0 = init_focus.GetParameter("x_1");
	fTipY_0 = init_focus.GetParameter("x_2");
	fTipElementGroup = list.GetParameter("focus_element_group"); fTipElementGroup--;
	fTipOutputVariable = list.GetParameter("focus_output_variable");

	/* resolve tracking method */
	const ParameterListT& tracking_method = list.GetListChoice(*this, "focus_tracking_method");
	if (tracking_method.Name() == "focus_at_maximum") {
		fTrackingType = kMax;
		fTipThreshold = tracking_method.GetParameter("threshold");
	}
	else if (tracking_method.Name() == "focus_right_most") {
		fTrackingType = kRightMost;
		fTipThreshold = tracking_method.GetParameter("threshold");	
	}
	else if (tracking_method.Name() == "focus_left_most") {
		fTrackingType = kLeftMost;
		fTipThreshold = tracking_method.GetParameter("threshold");	
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized tracking method \"%s\"", tracking_method.Name().Pointer());

	/* damping */
	fDampingWidth = list.GetParameter("damping_width");
	fDampingCoefficient = list.GetParameter("damping_coefficient");

	/* set stretching BC cards */
	fKBC_Cards.Dimension(nsd*fTopNodes.Length());
	int node = 0;
	double valueby2 = fULBC_Value/2.0;
	for (int i = 0; i < fTopNodes.Length(); i++) {
		KBC_CardT& card = fKBC_Cards[node++];
		card.SetValues(fTopNodes[i], 1, fULBC_Code, fULBC_Schedule, valueby2);
	}
	for (int i = 0; i < fTopNodes.Length(); i++) {
		KBC_CardT& card = fKBC_Cards[node++];
		card.SetValues(fTopNodes[i], 0, KBC_CardT::kFix, 0, 0);
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
	file.Root(fSupport.InputFile());
	file.Append(".tracking");
	fTrackingOutput.open(file);

	/* create controller for the right edge of the domain */
	fRightEdge = new KBC_ControllerT(fSupport);
	fField.AddKBCController(fRightEdge);

	/*find the right edge*/
	int nnd = init_coords.MajorDim();
	iAutoArrayT rightnodes(0);
	const double* px = init_coords.Pointer();
	//double rightmost = TrackPoint(kRightMost,kSmall);
	/* find and store right edge */
	for (int i = 0; i < nnd; i++)
	  {
	    //if (fabs(*px - rightmost) < kSmall) rightnodes.Append(i);
	    if (fabs(*px - fX_Right) < kSmall) rightnodes.Append(i);
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

/**********************************************************************
 * Protected
 **********************************************************************/

/* reset system to new center */
bool ConveyorSymT::SetSystemFocus(double focus)
{
	/* no need to shift the window */
	if (fX_Right - focus > fRightMinSpacing) return false;

        /* reference coordinates */
        const dArray2DT& initial_coords = fSupport.InitialCoordinates();

        /* fields */
        dArray2DT& u_field = fField[0];
        dArray2DT* Du_field = NULL;
        if (fField.Order() > 0) Du_field = &(fField[1]);
        dArray2DT* DDu_field = NULL;
        if (fField.Order() > 1) DDu_field = &(fField[2]);


	/* shift window */
	fX_Left  += fWindowShiftDistance;
	fX_Right += fWindowShiftDistance;

	/* model information */
	ModelManagerT& model = fSupport.ModelManager();

        /* nodes to get reference stretch */
        double  Y_top = initial_coords(fTopNodes[0], 1);
        double uY_top = u_field(fTopNodes[0], 1);
        double duY_dY = (uY_top-fuY_interface)/(Y_top-fTipY_0);

	/* has damping */
	bool has_damping = (fabs(fDampingWidth) > kSmall && fabs(fDampingCoefficient) > kSmall) ? true : false;
	cout << "\nDisplacement of interface: " << fuY_interface;

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
		if (*px < fX_Left - fMeshRepeatLength/10.0)
		{ 
			/* store */
			fShiftedNodes.Append(i);
		
			/* shift reference coordinates */
			new_coords[0] = initial_coords(i,0) + fX_PeriodicLength;
			new_coords[1] = initial_coords(i,1);
			model.UpdateNode(new_coords, i);

			/* correct displacements */
			u_field(i,0) = 0.0;
                        u_field(i,1) = fuY_interface + duY_dY*(initial_coords(i,1) - fTipY_0);										  
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

