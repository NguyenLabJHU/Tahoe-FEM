/* $Id: PenaltyRegionT.cpp,v 1.12 2003-08-23 20:08:57 paklein Exp $ */
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

using namespace Tahoe;

const double Pi = acos(-1.0);

/* constructor */
PenaltyRegionT::PenaltyRegionT(FEManagerT& fe_manager,
	int group,
	const iArray2DT& eqnos,
	const dArray2DT& coords,
	const dArray2DT* vels):
	FBC_ControllerT(fe_manager, group),
	
	/* references to NodeManagerT data */
	rEqnos(eqnos),
	rCoords(coords),
	pVels(vels),
	
	/* wall parameters */
	fx0(rCoords.MinorDim()),
	fv0(rCoords.MinorDim()),
	fMass(0.0),
	fLTf(NULL),

	/* state variables */
	fh_max(0.0),
	fx(rCoords.MinorDim()),
	fv(rCoords.MinorDim()),
	fxlast(rCoords.MinorDim()),
	fvlast(rCoords.MinorDim())
{

}

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
	if (fFEManager.PrintInput())
		switch (fFEManager.OutputFormat())
		{
			case IOBaseT::kTahoe:
			case IOBaseT::kTecPlot:
			case IOBaseT::kEnSight:
			case IOBaseT::kEnSightBinary:
			case IOBaseT::kExodusII:
			
				
				fContactNodes++;
				out << fContactNodes.wrap(6) << '\n';
				fContactNodes--;				
				break;
	
			default:

				cout << "\n PenaltyRegionT::EchoData: unsupported output format: ";
				cout << fFEManager.OutputFormat() << endl;
				throw ExceptionT::kGeneralFail;
		}
}

/* initialize data */
void PenaltyRegionT::Initialize(void)
{
	/* allocate memory for equation numbers */
	int numDOF = rEqnos.MinorDim();
	fContactEqnos.Dimension(fNumContactNodes*numDOF);
	
	/* allocate memory for force vector */
	fContactForce2D.Dimension(fNumContactNodes,numDOF);
	fContactForce.Set(fNumContactNodes*numDOF, fContactForce2D.Pointer());
	fContactForce2D = 0.0; // will be generate impulse at ApplyPreSolve
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
	iArray2DT eqtemp(fContactNodes.Length(),
	                 rEqnos.MinorDim(),
	                 fContactEqnos.Pointer()); //alias to fContactEqnos
	
	/* set current equation number */
	eqtemp.RowCollect(fContactNodes, rEqnos);
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
	fFEManager.AssembleRHS(fGroup, fContactForce, fContactEqnos);
}

/* apply kinematic boundary conditions */
void PenaltyRegionT::InitStep(void)
{
	/* compute impulse acting on region */
	if (fSlow == 1)
		for (int i = 0; i < rCoords.MinorDim(); i++)
			fv[i] -= fFEManager.TimeStep()*fContactForce2D.ColumnSum(i)/fMass;
	else if (fSlow == 2)
		fv.SetToScaled(fLTf->Value(), fv0);
	
	/* compute new position */
	fx.AddScaled(fFEManager.TimeStep(), fv);
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

/* writing results */
void PenaltyRegionT::WriteOutput(ostream& out) const
{
	int d_width = out.precision() + kDoubleExtra;

	out << "\n P e n a l t y   R e g i o n   D a t a :\n\n";
	out << " Maximum penetration. . . . . . . . . . . . =\n"
	    << setw(d_width) << fh_max << '\n';
	out << " Position . . . . . . . . . . . . . . . . . =\n" << fx << '\n';
	out << " Velocity . . . . . . . . . . . . . . . . . =\n" << fv << '\n';
	
	/* contact force */
	out << " Contact force:\n";
	for (int i = 0; i < rCoords.MinorDim(); i++)
		out << setw(kDoubleWidth) << -fContactForce2D.ColumnSum(i) << '\n';
}
