/* $Id: ConveyorT.cpp,v 1.13 2005-01-19 08:56:44 paklein Exp $ */
#include "ConveyorT.h"
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
#include "ContinuumElementT.h"
#include "CSEAnisoT.h"
#include "CommunicatorT.h"

using namespace Tahoe;

#define _FIX_RIGHT_EDGE_ 1

/* constructor */
ConveyorT::ConveyorT(const BasicSupportT& support, FieldT& field):
	KBC_ControllerT(support),
	fField(field),
	fMeshRepeatLength(0.0),
	fWindowShiftDistance(0.0),
	fRightMinSpacing(0.0),
	fTrackingInterval(1),
	fULBC_Value(0.0),
	fULBC_Code(KBC_CardT::kFix),
	fULBC_ScheduleNumber(-1),
	fULBC_Schedule(NULL),
	fRightEdge(NULL),
	fDampingWidth(0.0),
	fDampingCoefficient(0.0),
	fDampingReset(true),
	fTipOutputCode(-1),
	fTipColumnNum(-1),
	fUy_node_upper(-1),
	fUy_node_lower(-1),
	fNumSamples(20)
{
	SetName("conveyor");
}

/* set to initial conditions */
void ConveyorT::InitialCondition(void)
{
	/* reset */
	fTrackingCount = 0;
	fTrackingPoint = fX_Left;
	fTrackingPoint_last = fTrackingPoint;
	fX_Left_last = fX_Left;
	fX_Right_last = fX_Right;

	/* sort nodes into upper and lower groups - could not be done during 
	 * ConveyorT::TakeParameterList because it requires element group information
	 * that is not available then */
	const dArray2DT& init_coords = fSupport.InitialCoordinates();	 
	iAutoArrayT rightnodes_U(25), rightnodes_L(25);
	AutoArrayT<double> rightnodes_Y_U(25), rightnodes_Y_L(25);
	for (int i = 0; i < fShiftedNodesU.Length(); i++)
	{
		int nd = fShiftedNodesU[i];
		if (UpperLower(nd) == 1) {
	    	rightnodes_U.Append(nd);
	    	rightnodes_Y_U.Append(init_coords(nd,1));
	    }
	    else {
	    	rightnodes_L.Append(nd);
	    	rightnodes_Y_L.Append(init_coords(nd,1));
	    }
	}	
	fShiftedNodesU.Dimension(rightnodes_U);
	rightnodes_U.CopyInto(fShiftedNodesU);
	fShiftedNodesL.Dimension(rightnodes_L);
	rightnodes_L.CopyInto(fShiftedNodesL);

	/* sort in ascending y-coordinates */
	iArrayT shifted_nodes_tmp;
	shifted_nodes_tmp.Alias(fShiftedNodesU);
	shifted_nodes_tmp.SortAscending(rightnodes_Y_U);
	shifted_nodes_tmp.Alias(fShiftedNodesL);
	shifted_nodes_tmp.SortAscending(rightnodes_Y_L);

	/* mark elements linking left to right edge as inactive */
	MarkElements();	
}

void ConveyorT::Reset(void)
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
void ConveyorT::InitStep(void)
{
	/* inherited */
	KBC_ControllerT::InitStep();
	fTrackingCount++;

	/* update running average of interface compliance */
	dArray2DT& u_field = fField[0];
	if (fUy_node_upper != -1) {
		fUy_samples_upper.Push(u_field(fUy_node_upper,1));
		fUy_samples_upper.Resize(fNumSamples);
	}
	if (fUy_node_lower != -1) {
		fUy_samples_lower.Push(u_field(fUy_node_lower,1));
		fUy_samples_lower.Resize(fNumSamples);
	}

#ifdef _FIX_RIGHT_EDGE_
	/* coordinates */
	const dArray2DT& initial_coords = fSupport.InitialCoordinates();
	const dArray2DT& current_coords = fSupport.CurrentCoordinates();

	/* redistribute y-displacements on the right edge */	
	double Y_top = (   fTopNodes.Length() > 0) ? initial_coords(   fTopNodes[0], 1) : fTipY_0;
	double Y_bot = (fBottomNodes.Length() > 0) ? initial_coords(fBottomNodes[0], 1) : fTipY_0;
	double Y_U_LB = (fUy_node_upper != -1) ? initial_coords(fUy_node_upper, 1) : fTipY_0;
	double Y_L_UB = (fUy_node_lower != -1) ? initial_coords(fUy_node_lower, 1) : fTipY_0;

	double uY_top = (   fTopNodes.Length() > 0) ? u_field(   fTopNodes[0], 1) : 0.0;
	double uY_bot = (fBottomNodes.Length() > 0) ? u_field(fBottomNodes[0], 1) : 0.0;
	double uY_U_LB = 0.0;
	double uY_L_UB = 0.0;
	if (fUy_node_upper != -1) /* use running average */{
		dArrayT tmp;
		tmp.Alias(fUy_samples_upper);
		uY_U_LB = tmp.Sum()/fNumSamples;
	}
	if (fUy_node_lower != -1) /* use running average */{
		dArrayT tmp;
		tmp.Alias(fUy_samples_lower);
		uY_L_UB = tmp.Sum()/fNumSamples;
	}

	/* displacement gradient */
	double dY;
	dY = Y_top - Y_U_LB;
	double duY_dY_U = (fabs(dY) > kSmall) ? (uY_top - uY_U_LB)/dY : 0.0;
	dY = Y_L_UB - Y_bot;
	double duY_dY_L = (fabs(dY) > kSmall) ? (uY_L_UB - uY_bot)/dY : 0.0;

	/* reset cards for right edge */
	ArrayT<KBC_CardT>& cards = fRightEdge->KBC_Cards();
	KBC_CardT* card = fRightEdge->KBC_Cards().Pointer();
	for (int i = 0; i < fShiftedNodesU.Length(); i++) 
	{
		int nd = fShiftedNodesU[i];
		
		/* fix x-displacement */
		card->SetValues(nd, 0, KBC_CardT::kFix, NULL, 0.0);
		card++;

		/* interpolate y-displacement between interface and upper boundary */
		double uy = uY_U_LB + duY_dY_U*(initial_coords(nd,1) - Y_U_LB); 
		card->SetValues(nd, 1, KBC_CardT::kDsp, NULL, uy);
		card++;
	}
	for (int i = 0; i < fShiftedNodesL.Length(); i++) 
	{
		int nd = fShiftedNodesL[i];
		
		/* fix x-displacement */
		card->SetValues(nd, 0, KBC_CardT::kFix, NULL, 0.0);
		card++;

		/* interpolate y-displacement between lower boundary and interface */
		double uy = uY_bot + duY_dY_L*(initial_coords(nd,1) - Y_bot); 
		card->SetValues(nd, 1, KBC_CardT::kDsp, NULL, uy);
		card++;
	}
#endif

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

//TEMP
		bool print_force = false;
		if (print_force) {
			ofstreamT& out = fSupport.Output();
			int prec = out.precision();
			out.precision(12);
			out << "\n ConveyorT::FormRHS: time = " << fSupport.Time() << "\ndamping force = \n" << fDampingForce << '\n';
			out.precision(prec);
		}
	}
}

/* apply the update to the solution. Does nothing by default. */
void ConveyorT::Update(const dArrayT& update)
{
	/* inherited */
	KBC_ControllerT::Update(update);

	/* reset damping force coefficients - reset must be done here to ensure LHS is up to date */
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
	if (fSupport.Rank() == 0 && fTrackingCount == fTrackingInterval) {
	
		/* boundary displacement */
		dArray2DT& u_field = fField[0];

		/* nodes to get reference stretch */
		double uY_bot = (fBottomNodes.Length() > 0) ? u_field(fBottomNodes[0], 1) : 0.0;
		double uY_top = (   fTopNodes.Length() > 0) ? u_field(   fTopNodes[0], 1) : 0.0;
		
		fTrackingOutput << fSupport.Time() << ' ' 
		                << fTrackingPoint << ' '
		                << uY_top - uY_bot <<'\n';
		fTrackingCount = 0;
	}

	/* update history */
	fX_Left_last = fX_Left;
	fX_Right_last = fX_Right;
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

void ConveyorT::ReadRestart(ifstreamT& in)
{
	/* inherited */
	KBC_ControllerT::ReadRestart(in);

	/* external file */
	StringT file = in.filename();
	file.Append(".", fField.FieldName());
	file.Append(".", Name());
	ifstreamT my_in(file);
	if (!my_in.is_open()) 
		ExceptionT::GeneralFail("ConveyorT::ReadRestart", "could not open file \"%s\"",
			file.Pointer());

	/* read dimensions */
	my_in >> fTrackingCount >> fTrackingPoint >> fX_Left >> fX_Right;

	/* read node lists */
	iArrayT tmp;
	int num_shifted = -1;
	my_in >> num_shifted;
	if (num_shifted > 0) {
		tmp.Dimension(num_shifted);
		my_in >> tmp;
		fShiftedNodesU = tmp;
	}
	my_in >> num_shifted;
	if (num_shifted > 0) {
		tmp.Dimension(num_shifted);
		my_in >> tmp;
		fShiftedNodesL = tmp;
	}

#if _FIX_RIGHT_EDGE_
	num_shifted = fShiftedNodesU.Length() + fShiftedNodesL.Length();
	/* reset cards for right edge*/
	const dArray2DT& u_field = fField[0];
	fRightEdge->KBC_Cards().Dimension(2*num_shifted);
	KBC_CardT* card = fRightEdge->KBC_Cards().Pointer();
	for (int i = 0; i < fShiftedNodesU.Length(); i++) {
		int nd = fShiftedNodesU[i];
		card->SetValues(nd, 0, KBC_CardT::kFix, NULL, 0.0);
		card++;
		card->SetValues(nd, 1, KBC_CardT::kDsp, NULL, u_field(nd,1));
		card++;
	}
	for (int i = 0; i < fShiftedNodesL.Length(); i++) {
		int nd = fShiftedNodesL[i];
		card->SetValues(nd, 0, KBC_CardT::kFix, NULL, 0.0);
		card++;
		card->SetValues(nd, 1, KBC_CardT::kDsp, NULL, u_field(nd,1));
		card++;
	}
#endif

	int reset, num_damped = -1;
	my_in >> reset;
	fDampingReset = (reset) ? true : false;
	my_in >> num_damped;
	tmp.Dimension(num_damped);
	my_in >> tmp;
	fDampingNodes = tmp;

	int ndof = fField.NumDOF();
	int ndn = fDampingNodes.Length();
	fDampingForce.Dimension(ndn, ndof);
	fDampingCoeff.Dimension(ndn, ndof);
	fDampingEqnos.Dimension(ndn, ndof);
	my_in >> fDampingCoeff >> fDampingEqnos;

	/* interface displacement history */
	my_in >> fUy_node_upper >> fUy_node_lower;
	dArrayT samples_tmp;
	samples_tmp.Alias(fUy_samples_upper);
	my_in >> samples_tmp;
	samples_tmp.Alias(fUy_samples_lower);
	my_in >> samples_tmp;

	/* shifted reference coordinates */
	ModelManagerT& model = fSupport.ModelManager();
	dArray2DT shifted_init_coords(model.NumNodes(), model.NumDimensions());
	my_in >> shifted_init_coords;
	model.UpdateNodes(shifted_init_coords, false);
	fSupport.NodeManager().UpdateCurrentCoordinates();

	/* initialize history */
	fTrackingPoint_last = fTrackingPoint;
	fX_Left_last = fX_Left;
	fX_Right_last = fX_Right;

//TEMP - need to reset the equation system because it was set before
//       reading the restart files and does not reflect the equations
//       the system had when the restart was written because the kbc's
//       on right edge nodes where not in place at the time.
FEManagerT& fe_man = const_cast<FEManagerT&>(fSupport.FEManager());
fe_man.SetEquationSystem(fField.Group(), 0); // what about the equation start shift?
}

void ConveyorT::WriteRestart(ofstreamT& out) const
{
	/* inherited */
	KBC_ControllerT::WriteRestart(out);

	/* external file */
	StringT file = out.filename();
	file.Append(".", fField.FieldName());
	file.Append(".", Name());
	ofstreamT my_out(file);
	my_out.precision(out.precision());

	/* write dimensions */
	my_out << fTrackingCount << '\n'
	    << fTrackingPoint << '\n'
	    << fX_Left << '\n'
	    << fX_Right << '\n';
	
	/* nodes on right edge */
	iArrayT tmp;
	tmp.Alias(fShiftedNodesU);
	my_out << tmp.Length() << '\n' << tmp.wrap(10) << '\n';
	tmp.Alias(fShiftedNodesL);
	my_out << tmp.Length() << '\n' << tmp.wrap(10) << '\n';
	
	/* damping */
	my_out << ((fDampingReset) ? 1 : 0) << '\n';
	tmp.Alias(fDampingNodes);
	my_out << tmp.Length() << '\n' << tmp.wrap(10) << '\n';
	my_out << fDampingCoeff << '\n';
	my_out << fDampingEqnos << '\n';

	/* interface displacement history */
	my_out << fUy_node_upper << ' ' << fUy_node_lower << '\n';
	dArrayT samples_tmp;
	samples_tmp.Alias(fUy_samples_upper);
	my_out << samples_tmp << '\n';
	samples_tmp.Alias(fUy_samples_lower);
	my_out << samples_tmp << '\n';

	/* write the modified reference coordinates */
	my_out << fSupport.InitialCoordinates() << '\n';
}

/* describe the parameters needed by the interface */
void ConveyorT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	KBC_ControllerT::DefineParameters(list);

	/* bound */
	LimitT zero_bound(0.0, LimitT::Lower);
	LimitT zero_bound_inc(0.0, LimitT::LowerInclusive);

	/* dimensions */
	ParameterT repeat_length(fMeshRepeatLength, "mesh_repeat_length");
	repeat_length.AddLimit(zero_bound);
	list.AddParameter(repeat_length);

	ParameterT window_shift(fWindowShiftDistance, "window_shift");
	window_shift.AddLimit(zero_bound);
	list.AddParameter(window_shift);

	ParameterT min_right_space(fRightMinSpacing, "min_right_space");
	min_right_space.AddLimit(zero_bound);
	list.AddParameter(min_right_space);

	ParameterT num_samples(fNumSamples, "interface_displacement_samples");
	num_samples.AddLimit(1, LimitT::LowerInclusive);
	num_samples.SetDefault(fNumSamples);
	list.AddParameter(num_samples, ParameterListT::ZeroOrOnce);

	/* tip tracking */
	ParameterT tracking_increment(fTrackingInterval, "focus_tracking_increment");
	tracking_increment.SetDefault(fTrackingInterval);
	list.AddParameter(tracking_increment);
	list.AddParameter(ParameterT::Integer, "focus_element_group");
	list.AddParameter(ParameterT::Word, "focus_output_variable");

	/* boundary damping */
	ParameterT damping_width(fDampingWidth, "damping_width");
	damping_width.AddLimit(zero_bound_inc);
	list.AddParameter(damping_width);

	ParameterT damping_coeff(fDampingCoefficient, "damping_coefficient");
	damping_coeff.AddLimit(zero_bound_inc);
	list.AddParameter(damping_coeff);
}

/* information about subordinate parameter lists */
void ConveyorT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	KBC_ControllerT::DefineSubs(sub_list);

	/* tip trackng information */
	sub_list.AddSub("focus_tracking_method", ParameterListT::Once, true);

	/* KBC for lower and upper surfaces */
	sub_list.AddSub("lower_upper_kinematic_BC");

	/* initial tip coordinates */
	sub_list.AddSub("initial_focus_coordinates");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ConveyorT::NewSub(const StringT& name) const
{
	const char caller[] = "ConveyorT::NewSub";

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
	else if (name == "lower_upper_kinematic_BC")
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
		kbc->AddSub("lower_ID_list", ParameterListT::ZeroOrOnce);
		kbc->AddSub("upper_ID_list", ParameterListT::ZeroOrOnce);

		return kbc;
	}
	else /* inherited */
		return KBC_ControllerT::NewSub(name);
}

/* accept parameter list */
void ConveyorT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "ConveyorT::TakeParameterList";

	/* only 2D for now */
	int nsd = fSupport.NumSD();
	if (nsd != 2) ExceptionT::GeneralFail(caller, "only tested for 2D: %d", nsd);

	/* inherited */
	KBC_ControllerT::TakeParameterList(list);

	/* dimensions */
	fMeshRepeatLength = list.GetParameter("mesh_repeat_length");
	fWindowShiftDistance = list.GetParameter("window_shift");
	fRightMinSpacing = list.GetParameter("min_right_space");

	const ParameterT* num_samples = list.Parameter("interface_displacement_samples");
	if (num_samples)
		fNumSamples = *num_samples;

	/* boundary stretching */
	const ParameterListT& ul_kbc = list.GetList("lower_upper_kinematic_BC");
	int i_code = ul_kbc.GetParameter("type");
	fULBC_Code = KBC_CardT::int2CodeT(i_code + 1);
	fULBC_Value = ul_kbc.GetParameter("value");
	fULBC_ScheduleNumber = ul_kbc.GetParameter("schedule"); fULBC_ScheduleNumber--;
	fULBC_Schedule = fSupport.Schedule(fULBC_ScheduleNumber);	
	if (!fULBC_Schedule) ExceptionT::BadInputValue(caller, "could not resolve schedule %d", fULBC_ScheduleNumber+1);
	
	ArrayT<StringT> id_list;
	const ParameterListT* lower_nodes = ul_kbc.List("lower_ID_list");
	if (lower_nodes) {
		StringListT::Extract(*lower_nodes, id_list);	
		GetNodes(id_list, fBottomNodes);
	}
	const ParameterListT* upper_nodes = ul_kbc.List("upper_ID_list");
	if (upper_nodes) {
		StringListT::Extract(*upper_nodes,  id_list);	
		GetNodes(id_list, fTopNodes);
	}
	if (!upper_nodes && !lower_nodes)
		ExceptionT::GeneralFail(caller, "expecting at least \"lower_ID_list\" or \"lower_ID_list\"");

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
		card.SetValues(fBottomNodes[i], 0, KBC_CardT::kFix, NULL, 0.0);
	}
	for (int i = 0; i < fTopNodes.Length(); i++) {
		KBC_CardT& card = fKBC_Cards[node++];
		card.SetValues(fTopNodes[i], 0, KBC_CardT::kFix, NULL, 0.0);
	}

	/* find boundaries */
	const dArray2DT& init_coords = fSupport.InitialCoordinates();
	dArrayT X2(init_coords.MajorDim());
	init_coords.ColumnCopy(0, X2);
	X2.MinMax(fX_Left, fX_Right);
	X2.Free();
	
	/* exchange bounds */
	fX_Left  = fSupport.Communicator().Min(fX_Left );
	fX_Right = fSupport.Communicator().Max(fX_Right);

	/* check */
	if (fX_Right - fTipX_0 < fRightMinSpacing)
		ExceptionT::GeneralFail(caller, "initial tip position %g violates \"min_right_space\" %g",
			fTipX_0, fRightMinSpacing);
	
	/* set the periodic distance */
	fX_PeriodicLength = fX_Right - fX_Left + fMeshRepeatLength;
	fWidthDeadZone = fMeshRepeatLength*1.5;
	if (fWidthDeadZone > fRightMinSpacing) ExceptionT::GeneralFail(caller);

	/* open file for tracking information */
	if (fSupport.Rank() == 0) {
		StringT file;
		file.Root(fSupport.InputFile());
		file.Append(".tracking");
		fTrackingOutput.open(file);
	}

	/* create controller for the right edge of the domain */
	fRightEdge = new KBC_ControllerT(fSupport);
	fField.AddKBCController(fRightEdge);

	/* find and store right edge - in ConveyorT::InitialCondition, these will be sorted
	 * into fShiftedNodesU and fShiftedNodesU, but that can't be done now because it requires
	 * element group information that's not available now */
	int nnd = init_coords.MajorDim();
	iAutoArrayT rightnodes(25);
	const double* px = init_coords.Pointer();
	for (int i = 0; i < nnd; i++) {
	    if (fabs(*px - fX_Right) < kSmall) 
	    	if (!fTopNodes.HasValue(i) && !fBottomNodes.HasValue(i)) /* do not include nodes on the lower/upper edges */
				rightnodes.Append(i);
	    px += nsd;
	}
	fShiftedNodesU.Dimension(rightnodes);
	rightnodes.CopyInto(fShiftedNodesU);

#if _FIX_RIGHT_EDGE_
	/* fix the right edge */
	int num_shifted = fShiftedNodesU.Length();
	ArrayT<KBC_CardT>& cards = fRightEdge->KBC_Cards();
	cards.Dimension(2*num_shifted);
	KBC_CardT* card = cards.Pointer();
	for (int i = 0; i < fShiftedNodesU.Length(); i++) {
		int nd = fShiftedNodesU[i];
		card->SetValues(nd, 0, KBC_CardT::kFix, NULL, 0.0);
		card++;
		card->SetValues(nd, 1, KBC_CardT::kDsp, NULL, 0.0);
		card++;
	}
#endif

	/* initialize running average */
	fUy_samples_upper.Dimension(fNumSamples);
	fUy_samples_upper = 0.0;
	fUy_samples_lower.Dimension(fNumSamples);
	fUy_samples_lower = 0.0;
}

/**********************************************************************
 * Protected
 **********************************************************************/

/* locate new tracking point */
double ConveyorT::TrackPoint(TrackingTypeT tracking_type, double threshold)
{
	const char caller[] = "ConveyorT::TrackPoint";

	/* near tip element group */
	ElementBaseT& neartip_group = fSupport.ElementGroup(fTipElementGroup);	
	
	/* resolve tracking variable */
	if (fTipOutputCode == -1 || fTipColumnNum == -1) {

		/* try to resolve output variable */
		neartip_group.ResolveOutputVariable(fTipOutputVariable, fTipOutputCode, fTipColumnNum);
		
		/* check */
		if (fTipOutputCode == -1 || fTipColumnNum == -1)
			ExceptionT::GeneralFail(caller, 
				"could not resolve output variable \"%s\" in element group %d",
					fTipOutputVariable.Pointer(), fTipElementGroup+1);
	}

	/* signal to accumulate nodal values */
	neartip_group.SendOutput(fTipOutputCode);

	/* nodes */
	NodeManagerT& node_manager = fSupport.NodeManager();

	/* tracking type */
	if (tracking_type == kMax)
	{
//TEMP
if (fSupport.Size() != 1)
	ExceptionT::GeneralFail(caller, "tracking method %d not parallelized", tracking_type);

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

		/* exchange result */
		right_most = fSupport.Communicator().Max(right_most);

		return right_most;
	}
	else if (tracking_type == kLeftMost)
	{
		/* reference coordinates */
		const dArray2DT& ref_coords = node_manager.InitialCoordinates();
	
		/* find left most value greater than threshold */
		double left_most = fX_Right;
		dArrayT values;
		node_manager.TopAverage();
		int node = node_manager.NextAverageRow(values);
		while (node != -1)
		{
			if (values[fTipColumnNum] > threshold && ref_coords(node,0) < left_most)
				left_most = ref_coords(node,0);
		
			node = node_manager.NextAverageRow(values);
		}

		/* exchange result */
		left_most = fSupport.Communicator().Min(left_most);

		return left_most;
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

//TEMP - not parallelized
if (fSupport.Size() > 1)
	ExceptionT::GeneralFail("ConveyorT::SetSystemFocus", "not parallelized");

	/* shift window */
	fX_Left  += fWindowShiftDistance;
	fX_Right += fWindowShiftDistance;

	/* model information */
	ModelManagerT& model = fSupport.ModelManager();

	/* coordinates */
	const dArray2DT& initial_coords = fSupport.InitialCoordinates();
	const dArray2DT& current_coords = fSupport.CurrentCoordinates();

	/* fields */
	dArray2DT& u_field = fField[0];
	dArray2DT* Du_field = NULL;
	if (fField.Order() > 0) Du_field = &(fField[1]);
	dArray2DT* DDu_field = NULL;
	if (fField.Order() > 1) DDu_field = &(fField[2]);
	
	/* has boundary conditions */
	bool has_top = fTopNodes.Length() > 0;
	bool has_bot = fBottomNodes.Length() > 0;

#if 0
	/* reset x-displacement functions */
	dArray2DT points;
	points.Dimension(upper_nodes.Length(), 2);
	for (int i = 0; i < upper_nodes.Length(); i++) {
		int nd = upper_nodes[i];
		points(i, 0) = initial_coords(nd, 1);
		points(i, 1) = u_field(nd, 0);
	}
	fUx_upper.SetPoints(points);	

	points.Dimension(lower_nodes.Length(), 2);
	for (int i = 0; i < lower_nodes.Length(); i++) {
		int nd = lower_nodes[i];
		points(i, 0) = initial_coords(nd, 1);
		points(i, 1) = u_field(nd, 0);
	}
	fUx_lower.SetPoints(points);
#endif

	/* nodes to get reference stretch */
	double Y_top = (has_top) ? initial_coords(   fTopNodes[0], 1) : fTipY_0;
	double Y_bot = (has_bot) ? initial_coords(fBottomNodes[0], 1) : fTipY_0;
	double Y_U_LB = (has_top) ? initial_coords(fShiftedNodesU[0], 1) : fTipY_0;
	double Y_L_UB = (has_bot) ? initial_coords(fShiftedNodesL.Last(), 1) : fTipY_0;

	double uY_top = (has_top) ? u_field(   fTopNodes[0], 1) : 0.0;
	double uY_bot = (has_bot) ? u_field(fBottomNodes[0], 1) : 0.0;
	double uY_U_LB = (has_top) ? u_field(fShiftedNodesU[0], 1) : 0.0;
	double uY_L_UB = (has_bot) ? u_field(fShiftedNodesL.Last(), 1) : 0.0;
	if (fUy_node_upper != -1) /* use running average */{
		dArrayT tmp;
		tmp.Alias(fUy_samples_upper);
		uY_U_LB = tmp.Sum()/fNumSamples;
	}
	if (fUy_node_lower != -1) /* use running average */{
		dArrayT tmp;
		tmp.Alias(fUy_samples_lower);
		uY_L_UB = tmp.Sum()/fNumSamples;
	}

	/* reset tracking node */
	fUy_node_upper = (has_top) ? fShiftedNodesU[0] : -1;
	fUy_node_lower = (has_bot) ? fShiftedNodesL.Last() : -1;

	/* displacement gradient */
	double dY;
	dY = Y_top - Y_U_LB;
	double duY_dY_U = (fabs(dY) > kSmall) ? (uY_top - uY_U_LB)/dY : 0.0;
	dY = Y_L_UB - Y_bot;
	double duY_dY_L = (fabs(dY) > kSmall) ? (uY_L_UB - uY_bot)/dY : 0.0;

	/* has damping */
	bool has_damping = (fabs(fDampingWidth) > kSmall && fabs(fDampingCoefficient) > kSmall) ? true : false;

	/* shift reference coordinates and correct fields */
	int nnd = initial_coords.MajorDim();
	int nsd = initial_coords.MinorDim();
	const double* px = initial_coords.Pointer();
	dArrayT new_coords(nsd);
	fDampingNodes.Dimension(0);
	fShiftedNodesU.Dimension(0);
	fShiftedNodesL.Dimension(0);
	AutoArrayT<double> shifter_nodes_Y_U(25), shifter_nodes_Y_L(25);
	for (int i = 0; i < nnd; i++)
	{
		/* node outside the window */
		if (*px < fX_Left - fMeshRepeatLength/10.0)
		{
			/* above or below the cleavage plane */
			int upper_lower = UpperLower(i);
		
			/* store - do not include nodes on the upper and lower surface */
			if (!fTopNodes.HasValue(i) && !fBottomNodes.HasValue(i)) {
				if (upper_lower == 1) {
					fShiftedNodesU.Append(i);
					shifter_nodes_Y_U.Append(px[1]);
				} else {
					fShiftedNodesL.Append(i);
					shifter_nodes_Y_L.Append(px[1]);
				}
			}
		
			/* shift reference coordinates */
			new_coords[0] = initial_coords(i,0) + fX_PeriodicLength;
			new_coords[1] = initial_coords(i,1);
			model.UpdateNode(new_coords, i);

			/* correct displacements */
			u_field(i,0) = 0.0;
			if (upper_lower == 1)
				u_field(i,1) = uY_U_LB + duY_dY_U*(initial_coords(i,1) - Y_U_LB); /* interpolate between interface and upper boundary */			
			else
				u_field(i,1) = uY_bot + duY_dY_L*(initial_coords(i,1) - Y_bot); /* interpolate between lower boundary and interface */

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

	/* sort shifted nodes in ascending y-coordinate */
	iArrayT right_nodes_tmp;
	right_nodes_tmp.Alias(fShiftedNodesU);
	right_nodes_tmp.SortAscending(shifter_nodes_Y_U);
	right_nodes_tmp.Alias(fShiftedNodesL);
	right_nodes_tmp.SortAscending(shifter_nodes_Y_L);

#if _FIX_RIGHT_EDGE_
	/* reset cards for right edge */
	int num_shifted = fShiftedNodesU.Length() + fShiftedNodesL.Length();
	ArrayT<KBC_CardT>& cards = fRightEdge->KBC_Cards();
	cards.Dimension(2*num_shifted);
	KBC_CardT* card = cards.Pointer();
	for (int i = 0; i < fShiftedNodesU.Length(); i++) {
		int nd = fShiftedNodesU[i];
		card->SetValues(nd, 0, KBC_CardT::kFix, NULL, 0.0);
		card++;
		card->SetValues(nd, 1, KBC_CardT::kDsp, NULL, u_field(nd,1));
		card++;
	}
	for (int i = 0; i < fShiftedNodesL.Length(); i++) {
		int nd = fShiftedNodesL[i];
		card->SetValues(nd, 0, KBC_CardT::kFix, NULL, 0.0);
		card++;
		card->SetValues(nd, 1, KBC_CardT::kDsp, NULL, u_field(nd,1));
		card++;
	}
#endif

	/* mark elements linking left to right edge as inactive */
	MarkElements();

	/* signal change */
	return true;
}

/* mark elements linking left to right edge as inactive */
void ConveyorT::MarkElements(void)
{
	const char caller[] = "ConveyorT::MarkElements";

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
		if (!element_group->InGroup(fField.Group())) ExceptionT::GeneralFail(caller);
	
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
		
		ContinuumElementT* cont_elem = TB_DYNAMIC_CAST(ContinuumElementT*, element_group);
		if (!cont_elem) {
			CSEAnisoT* cse_aniso_elem = TB_DYNAMIC_CAST(CSEAnisoT*, element_group);
			if (!cse_aniso_elem)
				ExceptionT::GeneralFail(caller,  "could not cast element group %d to ContinuumElementT", 
					element_group+1);
			else cse_aniso_elem->SetStatus(status);
		}
		else cont_elem->SetStatus(status);
 	}
}

/* deactivate elements to create a pre-crack */
void ConveyorT::CreatePrecrack(void)
{
	const char caller[] = "ConveyorT::CreatePrecrack";

	/* element group */
	const FEManagerT& fe = fSupport.FEManager();
	const dArray2DT& current_coords = fSupport.CurrentCoordinates();
	ElementBaseT* element_group = fe.ElementGroup(fTipElementGroup);
	if (!element_group) ExceptionT::GeneralFail(caller);
	int nel = element_group->NumElements();
	int nen = element_group->NumElementNodes();

	/* local coordinate array */
	LocalArrayT curr_coords(LocalArrayT::kCurrCoords, nen, 2);
	curr_coords.SetGlobal(current_coords);
	
	//TEMP
	if (!element_group->InGroup(fField.Group())) ExceptionT::GeneralFail(caller);
	
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
		bool X_check_OK = false;
		bool Y_check_OK = true;
		const double* px = curr_coords(0);
		const double* py = curr_coords(1);
		for (int k = 0; Y_check_OK && k < nen; k++) {
		
			/* "behind" tip */
			X_check_OK = (*px++ < fTipX_0) ? true : X_check_OK; /* just one node needed */
			Y_check_OK = fabs(*py++ - fTipY_0) < kSmall; /* all nodes needed */
		}
			
		/* goes end to end */
		if (X_check_OK && Y_check_OK) status[j] = ElementBaseT::kOFF;
	}
		
	/* reset element status */
	element_group->SetStatus(status);
}

/* resolve the given node number into the area above or below the crack plane */
int ConveyorT::UpperLower(int node) const
{
	const dArray2DT& init_coords = fSupport.InitialCoordinates();
	double y_node = init_coords(node, 1);
	if (fabs(y_node - fTipY_0) < kSmall) /* on the crack plane */
	{
		/* run through element groups */
		const FEManagerT& fe = fSupport.FEManager();
		int num_element_groups = fe.NumElementGroups();
		for (int i = 0; i < num_element_groups; i++)
		{
			/* element group */
			ElementBaseT* element_group = fe.ElementGroup(i);
			int nel = element_group->NumElements();
			int nen = element_group->NumElementNodes();

			/* look for connected nodes in upper/lower region */
			for (int j = 0; j < nel; j++) {
				const iArrayT& nodes = element_group->ElementCard(j).NodesX();
				if (nodes.HasValue(node)) {
					for (int k = 0; k < nodes.Length(); k++)
					{
						double y_node_k = init_coords(nodes[k], 1);
						if (fabs(y_node_k - fTipY_0) > kSmall) /* not on the crack plane */ {
							if (y_node_k > fTipY_0)
								return 1;
							else
								return -1;
						}
					}
				}
			}
		}
		
		/* should not fall through */
		ExceptionT::GeneralFail("ConveyorT::UpperLower", "could not place node %d", node+1);
		
		return 0;
	}
	else if (y_node > fTipY_0)
		return 1;
	else
		return -1;
}
