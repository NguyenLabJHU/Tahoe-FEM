/* $Id: PenaltyRegionT.cpp,v 1.15.4.1 2004-03-31 16:20:14 paklein Exp $ */
/* created: paklein (04/30/1998) */
#include "PenaltyRegionT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include <ctype.h>

#include "toolboxConstants.h"
#include "GlobalT.h"
#include "FEManagerT.h"
#include "ModelManagerT.h"
#include "CommManagerT.h"
#include "fstreamT.h"
#include "ScheduleT.h"
#include "eIntegratorT.h"
#include "IOBaseT.h"
#include "OutputSetT.h"
#include "CommunicatorT.h"

#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "FieldSupportT.h"

using namespace Tahoe;

const double Pi = acos(-1.0);

#if 0
/* constructor */
PenaltyRegionT::PenaltyRegionT(FEManagerT& fe_manager,
	int group,
	const iArray2DT& eqnos,
	const dArray2DT& coords,
	const dArray2DT& disp,
	const dArray2DT* vels):
	FBC_ControllerT(fe_manager, group),
	
	/* references to NodeManagerT data */
	rEqnos(eqnos),
	rCoords(coords),
	rDisp(disp),
	pVels(vels),
	
	/* wall parameters */
	fx0(rCoords.MinorDim()),
	fv0(rCoords.MinorDim()),
	fMass(0.0),
	fLTf(NULL),

	/* state variables */
	fx(rCoords.MinorDim()),
	fv(rCoords.MinorDim()),
	fxlast(rCoords.MinorDim()),
	fvlast(rCoords.MinorDim()),
	fOutputID(-1)
{

}

PenaltyRegionT::PenaltyRegionT(void):
	fMass(0.0),
	fLTf(NULL),
	fOutputID(-1)
{

}
#endif

#if 0
/* input processing */
void PenaltyRegionT::EchoData(ifstreamT& in, ostream &out)
{
	/* echo parameters */
	in >> fx0;
	in >> fv0;
	in >> fk;
	in >> fSlow;

	/* motion control */
	int numLTf;
	if (fSlow == kSchedule)
	{
		in >> numLTf;
		numLTf--;
		fLTf = fFEManager.Schedule(numLTf);
	}
	else if (fSlow == kImpulse)
		in >> fMass;

	out << "\n P e n a l t y   R e g i o n   P a r a m e t e r s :\n\n";
	out << " Initial position. . . . . . . . . . . . . . . . =\n" << fx0 << '\n';
	out << " Initial velocity. . . . . . . . . . . . . . . . =\n" << fv0 << '\n';
	out << " Penalty stiffness . . . . . . . . . . . . . . . = " << fk << '\n';
	out << " Momentum option . . . . . . . . . . . . . . . . = " << fSlow << '\n';
	out << "    eq. " << kConstantVelocity << ", constant velocity\n";
	out << "    eq. " << kImpulse          << ", slow with contact impulse\n";
	out << "    eq. " << kSchedule         << ", velocity load time function\n";
	if (fSlow == kSchedule)
		out << " Velocity load time function . . . . . . . . . . = " << numLTf << endl;
	else if (fSlow == kImpulse)
		out << " Mass. . . . . . . . . . . . . . . . . . . . . . = " << fMass << '\n';

	/* checks */
	if (fk <= 0.0) throw ExceptionT::kBadInputValue;
	if (fSlow != kImpulse && fSlow != kConstantVelocity && fSlow != kSchedule) throw ExceptionT::kBadInputValue;
	if (fSlow == kImpulse && fMass <= 0.0) throw ExceptionT::kBadInputValue;
		
	/* read contact nodes */
	ModelManagerT* model = fFEManager.ModelManager ();

	/* read node set indexes */
	ArrayT<StringT> ns_ID;
	model->NodeSetList (in, ns_ID);

	if (ns_ID.Length() > 0)
	  {
	    model->ManyNodeSets (ns_ID, fContactNodes);
	    fNumContactNodes = fContactNodes.Length();
	  }
	else
	  fNumContactNodes = 0;
	
	/* remove "external" nodes */
	CommManagerT* comm = fFEManager.CommManager();
	const ArrayT<int>* p_map = comm->ProcessorMap();
	if (fNumContactNodes > 0 && p_map)
	{
		/* wrap it */
		iArrayT processor;
		processor.Alias(*p_map);
	
		/* count processor nodes */
		int rank = comm->Rank();
		int num_local_nodes = 0;
		for (int i = 0; i < fNumContactNodes; i++)
			if (processor[fContactNodes[i]] == rank)
				num_local_nodes++;
		
		/* remove off-processor nodes */
		if (num_local_nodes != fNumContactNodes)
		{
			/* report */
			out << " Number of external contact nodes (removed). . . = "
			    << fNumContactNodes - num_local_nodes << '\n';

			/* collect processor nodes */
			iArrayT nodes_temp(num_local_nodes);
			int index = 0;
			for (int i = 0; i < fNumContactNodes; i++)
				if (processor[fContactNodes[i]] == rank)
					nodes_temp[index++] = fContactNodes[i];
		
			/* reset values */
			fContactNodes.Swap(nodes_temp);
			fNumContactNodes = fContactNodes.Length();
		}
	}
	
	/* write contact nodes */
	out << " Number of contact nodes . . . . . . . . . . . . = "
	    << fNumContactNodes << endl;	
	if (fFEManager.PrintInput()) {
		fContactNodes++;
		out << fContactNodes.wrap(6) << '\n';
		fContactNodes--;				
	}	
}
#endif

/* initialize data */
void PenaltyRegionT::Initialize(void)
{
#if 0
	/* allocate memory for equation numbers */
	int numDOF = rEqnos.MinorDim();
	fContactEqnos.Dimension(fNumContactNodes*numDOF);
	
	/* allocate memory for force vector */
	fContactForce2D.Dimension(fNumContactNodes,numDOF);
	fContactForce.Set(fNumContactNodes*numDOF, fContactForce2D.Pointer());
	fContactForce2D = 0.0; // will be generate impulse at ApplyPreSolve

	/* allocate space to store gaps (for output) */
	fGap.Dimension(fNumContactNodes);
	fGap = 0.0;
#endif
}

/* form of tangent matrix */
GlobalT::SystemTypeT PenaltyRegionT::TangentType(void) const
{
	/* no tangent */
	return GlobalT::kDiagonal;
}

void PenaltyRegionT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* inherited */
	FBC_ControllerT::Equations(eq_1, eq_2);

	/* NOTE - just collect equation numbers, don't append
	 * since contact nodes only act on self */

	/* collect contact equation numbers */
	const iArray2DT& eqnos = Field().Equations();
	iArray2DT eqtemp(fContactNodes.Length(), eqnos.MinorDim(), fContactEqnos.Pointer());
	
	/* set current equation number */
	eqtemp.RowCollect(fContactNodes, eqnos);
}

void PenaltyRegionT::InitialCondition(void)
{
	/* initialize position and velocity */
	fx = fx0;
	fv = fv0;
	
	/* initialize history */
	fxlast = fx;
	fvlast = fv;
}

void PenaltyRegionT::ReadRestart(istream& in)
{
	/* inherited */
	FBC_ControllerT::ReadRestart(in);

	in >> fx;     // position
	in >> fv;     // velocity
	in >> fxlast; // last converged position
	in >> fvlast; // last converged velocity
}

void PenaltyRegionT::WriteRestart(ostream& out) const
{
	/* inherited */
	FBC_ControllerT::WriteRestart(out);

	out << fx << '\n';     // position
	out << fv << '\n';     // velocity
	out << fxlast << '\n'; // last converged position
	out << fvlast << '\n'; // last converged velocity
}

/* compute the nodal contribution to the residual force vector */
void PenaltyRegionT::ApplyRHS(void)
{
	double constKd = 0.0;
	int formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* recompute contact forces */
	ComputeContactForce(constKd);

	/* assemble */
	FieldSupport().AssembleRHS(fGroup, fContactForce, fContactEqnos);
}

/* apply kinematic boundary conditions */
void PenaltyRegionT::InitStep(void)
{
	/* the time step */
	double time_step = FieldSupport().TimeStep();

	/* compute impulse acting on region */
	if (fSlow == kImpulse)
		for (int i = 0; i < fContactForce2D.MinorDim(); i++)
			fv[i] -= time_step*fContactForce2D.ColumnSum(i)/fMass;
	else if (fSlow == kSchedule)
		fv.SetToScaled(fLTf->Value(), fv0);
	
	/* compute new position */
	fx.AddScaled(time_step, fv);
}

/* finalize step */
void PenaltyRegionT::CloseStep(void)
{
	/* store position and velocity */
	fxlast = fx;
	fvlast = fv;
}

/* reset to the last known solution */
void PenaltyRegionT::Reset(void)
{
	/* restore position and velocity */
	fx = fxlast;
	fv = fvlast;
}

/* register data for output */
void PenaltyRegionT::RegisterOutput(void)
{
	/* initialize connectivities */
	fContactNodes2D.Alias(fContactNodes.Length(), 1, fContactNodes.Pointer());
	
	/* output labels */
	const iArray2DT& eqnos = Field().Equations();
	int ndof = eqnos.MinorDim();
	int num_output = ndof + /* displacements */
	                 ndof + /* force */
	                 1;     /* penetration depth */
	ArrayT<StringT> n_labels(num_output);
	const char* d_labels[] = {"D_X", "D_Y", "D_Z"};
	const char* f_labels[] = {"F_X", "F_Y", "F_Z"};
	int index = 0;
	for (int i = 0; i < ndof; i++)
		n_labels[index++] = d_labels[i];
	for (int i = 0; i < ndof; i++)
		n_labels[index++] = f_labels[i];
	n_labels[index] = "h";
	
	/* register output */
	OutputSetT output_set(GeometryT::kPoint, fContactNodes2D, n_labels);
	fOutputID = FieldSupport().RegisterOutput(output_set);
}

/* writing results */
void PenaltyRegionT::WriteOutput(ostream& out) const
{
	int d_width = out.precision() + kDoubleExtra;

	/* mp support */
	CommunicatorT& comm = FieldSupport().Communicator();

	/* maximum penetration */
	double h_max = 0.0;
	if (fContactNodes.Length() > 0)
		h_max = fGap.Min();

	out << "\n P e n a l t y   R e g i o n   D a t a :\n\n";
	out << " Local maximum penetration. . . . . . . . . =\n"
	    << setw(d_width) << h_max << '\n';
	out << " Global maximum penetration . . . . . . . . =\n"
	    << setw(d_width) << comm.Min(h_max) << '\n';
	out << " Position . . . . . . . . . . . . . . . . . =\n" << fx << '\n';
	out << " Velocity . . . . . . . . . . . . . . . . . =\n" << fv << '\n';
	
	/* compute contact force */
	dArrayT loc_sum(fContactForce2D.MinorDim());
	for (int i = 0; i < fContactForce2D.MinorDim(); i++)
		loc_sum[i] = -fContactForce2D.ColumnSum(i);
	dArrayT global_sum(loc_sum.Length());
	comm.Sum(loc_sum, global_sum);

	/* write output */
	out << " Local contact force. . . . . . . . . . . . =\n";
	for (int i = 0; i < loc_sum.Length(); i++)
		out << setw(kDoubleWidth) << loc_sum[i] << '\n';

	out << " Glocal contact force . . . . . . . . . . . =\n";
	for (int i = 0; i < global_sum.Length(); i++)
		out << setw(kDoubleWidth) << global_sum[i] << '\n';

	/* collect output data */
	int ndof = fContactForce2D.MinorDim();
	int num_output = ndof + /* displacements */
	                 ndof + /* force */
	                 1;     /* penetration depth */
	dArray2DT n_values(fContactNodes.Length(), num_output);
	dArray2DT n_disp(fContactNodes.Length(), ndof);

	/* collect displacements */
	int index = 0;
	const dArray2DT& disp = Field()[0];
	n_disp.RowCollect(fContactNodes, disp);
	for (int i = 0; i < ndof; i++)
		n_values.ColumnCopy(index++, n_disp, i);

	/* collect the forces */
	for (int i = 0; i < ndof; i++)
		n_values.ColumnCopy(index++, fContactForce2D, i);	

	/* collect gaps */
	n_values.SetColumn(index, fGap);

	/* send output */
	dArray2DT e_values;
	FieldSupport().WriteOutput(fOutputID, n_values, e_values);
}

/* describe the parameters needed by the interface */
void PenaltyRegionT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FBC_ControllerT::DefineParameters(list);

	ParameterT k(fk, "stiffness");
	k.AddLimit("0.0", LimitT::LowerInclusive);
	list.AddParameter(k);
}

/* information about subordinate parameter lists */
void PenaltyRegionT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FBC_ControllerT::DefineSubs(sub_list);

	/* motion control */
	sub_list.AddSub("motion_control_choice", ParameterListT::Once, true);	

	/* nodes */
	sub_list.AddSub("node_ID_list");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* PenaltyRegionT::NewSub(const StringT& list_name) const
{
	if (list_name == "motion_control_choice") {

		ParameterContainerT* motion = new ParameterContainerT(list_name);
		motion->SetListOrder(ParameterListT::Choice);

		/* options */	
		motion->AddSub(ParameterContainerT("velocity_constant"));

		ParameterContainerT v_sched("velocity_schedule");
		v_sched.AddParameter(ParameterT::Integer, "schedule");
		motion->AddSub(v_sched);

		ParameterContainerT v_slow("integrate_impulse");
		v_slow.AddParameter(ParameterT::Double, "mass");
		motion->AddSub(v_slow);
		
		return motion;
	}
	if (list_name == "node_ID_list")
		return new StringListT(list_name);
	else /* inherited */
		return FBC_ControllerT::NewSub(list_name);
}

/* accept parameter list */
void PenaltyRegionT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "PenaltyRegionT::TakeParameterList";

	/* inherited */
	FBC_ControllerT::TakeParameterList(list);

	fk = list.GetParameter("stiffness");
	
	/* motion control */
	fSlow = kConstantVelocity;
	const ParameterListT* motion = list.ResolveListChoice(*this, "motion_control_choice");
	if (!motion) ExceptionT::GeneralFail(caller, "expecting \"motion_control_choice\"");
	if (motion->Name() == "velocity_schedule") {
		fSlow = kSchedule;
		int schedule = motion->GetParameter("schedule");
		schedule--;
		fLTf = FieldSupport().Schedule(schedule);
		if (!fLTf) ExceptionT::BadInputValue(caller, "could not resolve schedule %d", schedule+1);
	}
	else if (motion->Name() == "integrate_impulse") {
		fSlow = kImpulse;
		fMass = motion->GetParameter("mass");
	}

	/* affected nodes */
	ArrayT<StringT> ns_ID;
	StringListT::Extract(list.GetList("node_ID_list"), ns_ID);
	if (ns_ID.Length() > 0) {
		ModelManagerT& model = FieldSupport().ModelManager();
		model.ManyNodeSets(ns_ID, fContactNodes);
		fNumContactNodes = fContactNodes.Length();
	}
	else
	  fNumContactNodes = 0;
	
	/* remove "external" nodes */
	ofstreamT& out = FieldSupport().Output();
	CommManagerT& comm = FieldSupport().CommManager();
	const ArrayT<int>* p_map = comm.	ProcessorMap();
	if (fNumContactNodes > 0 && p_map)
	{
		/* wrap it */
		iArrayT processor;
		processor.Alias(*p_map);
	
		/* count processor nodes */
		int rank = comm.Rank();
		int num_local_nodes = 0;
		for (int i = 0; i < fNumContactNodes; i++)
			if (processor[fContactNodes[i]] == rank)
				num_local_nodes++;
		
		/* remove off-processor nodes */
		if (num_local_nodes != fNumContactNodes)
		{
			/* report */
			out << " Number of external contact nodes (removed). . . = "
			    << fNumContactNodes - num_local_nodes << '\n';

			/* collect processor nodes */
			iArrayT nodes_temp(num_local_nodes);
			int index = 0;
			for (int i = 0; i < fNumContactNodes; i++)
				if (processor[fContactNodes[i]] == rank)
					nodes_temp[index++] = fContactNodes[i];
		
			/* reset values */
			fContactNodes.Swap(nodes_temp);
			fNumContactNodes = fContactNodes.Length();
		}
	}
	
	/* write contact nodes */
	out << " Number of contact nodes . . . . . . . . . . . . = "
	    << fNumContactNodes << endl;	
	if (FieldSupport().PrintInput()) {
		fContactNodes++;
		out << fContactNodes.wrap(6) << '\n';
		fContactNodes--;				
	}	

	/* allocate memory for equation numbers */
	const iArray2DT& eqnos = Field().Equations();
	int numDOF = eqnos.MinorDim();
	fContactEqnos.Dimension(fNumContactNodes*numDOF);
	
	/* allocate memory for force vector */
	fContactForce2D.Dimension(fNumContactNodes,numDOF);
	fContactForce.Set(fNumContactNodes*numDOF, fContactForce2D.Pointer());
	fContactForce2D = 0.0; // will be generate impulse at ApplyPreSolve

	/* allocate space to store gaps (for output) */
	fGap.Dimension(fNumContactNodes);
	fGap = 0.0;
}
