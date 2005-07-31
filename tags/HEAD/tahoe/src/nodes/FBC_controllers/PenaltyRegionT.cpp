/* $Id: PenaltyRegionT.cpp,v 1.1.1.1 2001-01-29 08:20:40 paklein Exp $ */
/* created: paklein (04/30/1998)                                          */
/* base class for moving rigid, penalty regions. contact nodes            */
/* that enter the region are expelled by a quadratic penetration          */
/* potential. derived classes are responsilble for computing              */
/* the penetration depth and reaction for based on the geometry           */
/* of the region.                                                         */

#include "PenaltyRegionT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include <ctype.h>

#include "Constants.h"
#include "GlobalT.h"
#include "FEManagerT.h"
#include "ExodusT.h"
#include "ModelFileT.h"
#include "fstreamT.h"
#include "LoadTime.h"
#include "eControllerT.h"
#include "IOBaseT.h"

const double Pi = acos(-1.0);

/* constructor */
PenaltyRegionT::PenaltyRegionT(FEManagerT& fe_manager,
	const iArray2DT& eqnos,
	const dArray2DT& coords,
	const dArray2DT* vels):
	FBC_ControllerT(fe_manager),
	
	/* references to NodeManagerT data */
	rEqnos(eqnos),
	rCoords(coords),
	pVels(vels),
	
	/* wall parameters */
	fx0(rCoords.MinorDim()),
	fv0(rCoords.MinorDim()),
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
	in >> fMass;

	/* read load time information */
	int numLTf;
	if (fSlow == 2)
	{
		in >> numLTf;
		numLTf--;
		fLTf = fFEManager.GetLTfPtr(numLTf);
	}

	out << "\n P e n a l t y   R e g i o n   P a r a m e t e r s :\n\n";
	out << " Initial position. . . . . . . . . . . . . . . . =\n" << fx0 << '\n';
	out << " Initial velocity. . . . . . . . . . . . . . . . =\n" << fv0 << '\n';
	out << " Penalty stiffness . . . . . . . . . . . . . . . = " << fk << '\n';
	out << " Momentum option . . . . . . . . . . . . . . . . = " << fSlow << '\n';
	out << "    eq. 0, constant velocity\n";
	out << "    eq. 1, calculate contact impulse\n";
	out << "    eq. 2, velocity load time function\n";
	out << " Mass. . . . . . . . . . . . . . . . . . . . . . = " << fMass << '\n';
	if (fSlow == 2)
		out << " Velocity load time function . . . . . . . . . . = " << numLTf << endl;

	/* checks */
	if (fk <= 0.0) throw eBadInputValue;
	if (fSlow != 1 && fSlow != 0 && fSlow != 2) throw eBadInputValue;
	if (fSlow == 1 && fMass <= 0.0) throw eBadInputValue;
		
	/* read contact nodes */
	switch (fFEManager.InputFormat())
	{
		case IOBaseT::kTahoe:
		{
			ifstreamT tmp;
			ifstreamT& in2 = fFEManager.OpenExternal(in, tmp, out, true,
				"PenaltyRegionT::EchoData: could not open file");

			in2 >> fNumContactNodes;
			fContactNodes.Allocate(fNumContactNodes);
			in2 >> fContactNodes;
			break;
		}
		case IOBaseT::kTahoeII:
		{
			/* number of node sets */
			int num_sets;
			in >> num_sets;
			out << " Number of node set ID's: " << num_sets << endl;
			if (num_sets > 0)
			{
				/* open database */
				ModelFileT model_file;
				model_file.OpenRead(fFEManager.ModelFile());

				/* echo set ID's */
				iArrayT ID_list(num_sets);
				in >> ID_list;
				out << ID_list.wrap(10) << '\n';
				
				/* collect */
				if (model_file.GetNodeSets(ID_list, fContactNodes) !=
				    ModelFileT::kOK) throw eBadInputValue;
				
				/* dimension */
				fNumContactNodes = fContactNodes.Length();
				
			}
			else
				fNumContactNodes = 0;
				
			break;
		}
		case IOBaseT::kExodusII:
		{
			/* number of node sets */
			int num_sets;
			in >> num_sets;
			out << " Number of node set ID's: " << num_sets << endl;
			if (num_sets > 0)
			{
				/* echo set ID's */
				iArrayT ID_list(num_sets);
				in >> ID_list;
				out << ID_list.wrap(10) << '\n';

				/* open database */
				ExodusT database(out);
				database.OpenRead(fFEManager.ModelFile());
				
				/* read collect all nodes in sets */
				database.ReadNodeSets(ID_list, fContactNodes);
				
				/* dimension */
				fNumContactNodes = fContactNodes.Length();				
			}
			else
				fNumContactNodes = 0;
				
			break;
		}
		default:

			cout << "\n PenaltyRegionT::EchoData: unsupported input format: ";
			cout << fFEManager.InputFormat() << endl;
			throw eGeneralFail;
	}

	/* internal numbering */
	fContactNodes--;
	
	/* remove "external" nodes */
	iArrayT nodes_in;
	fFEManager.IncomingNodes(nodes_in);
	if (fNumContactNodes > 0 && nodes_in.Length() > 0)
	{
		/* brute force search */
		iArrayT is_external(fNumContactNodes);
		for (int i = 0; i < fNumContactNodes; i++)
			is_external[i] = nodes_in.HasValue(fContactNodes[i]);
	
		/* found external nodes */
		if (is_external.HasValue(1))
		{
			int num_external = is_external.Count(1);
			out << " Number of external contact nodes (removed). . . = "
			    << num_external << '\n';
		
			int num_contact_nodes = fNumContactNodes - num_external;
			iArrayT nodes_temp(num_contact_nodes);
			int count = 0;
			for (int i = 0; i < fNumContactNodes; i++)
				if (!is_external[i])
					nodes_temp[count++] = fContactNodes[i];
		
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
				throw eGeneralFail;
		}
}

/* initialize data */
void PenaltyRegionT::Initialize(void)
{
	/* allocate memory for equation numbers */
	int numDOF = rEqnos.MinorDim();
	fContactEqnos.Allocate(fNumContactNodes*numDOF);
	
	/* allocate memory for force vector */
	fContactForce2D.Allocate(fNumContactNodes,numDOF);
	fContactForce.Set(fNumContactNodes*numDOF, fContactForce2D.Pointer());
	fContactForce2D = 0.0; // will be generate impulse at ApplyPreSolve
}

void PenaltyRegionT::Reinitialize(void)
{
	// do nothing
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
	int formKd = fController->FormKd(constKd);
	if (!formKd) return;

	/* recompute contact forces */
	ComputeContactForce(constKd);

	/* assemble */
	fFEManager.AssembleRHS(fContactForce, fContactEqnos);
}

/* apply kinematic boundary conditions */
void PenaltyRegionT::InitStep(void)
{
	/* compute impulse acting on region */
	if (fSlow == 1)
		for (int i = 0; i < rCoords.MinorDim(); i++)
			fv[i] -= fFEManager.TimeStep()*fContactForce2D.ColumnSum(i)/fMass;
	else if (fSlow == 2)
		fv.SetToScaled(fLTf->LoadFactor(), fv0);
	
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
