/* $Id: AugLagSphereT.cpp,v 1.1.1.1 2001-01-29 08:20:40 paklein Exp $ */
/* created: paklein (03/24/1999)                                          */

#include "AugLagSphereT.h"

#include <iostream.h>
#include <iomanip.h>

#include "Constants.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "XDOF_ManagerT.h"
#include "eControllerT.h"

/* parameters */
const int kNumAugLagDOF = 1;

/* constructor */
AugLagSphereT::AugLagSphereT(FEManagerT& fe_manager, XDOF_ManagerT* XDOF_nodes,
	const iArray2DT& eqnos,
	const dArray2DT& coords,
	const dArray2DT* vels):
	PenaltySphereT(fe_manager, eqnos, coords, vels),
	fXDOF_Nodes(XDOF_nodes)
{
	/* (re-)dimension the tangent matrix */
	fLHS.Allocate(rEqnos.MinorDim() + 1); // additional DOF
}

/* initialize data */
void AugLagSphereT::Initialize(void)
{
	/* inherited */
	PenaltySphereT::Initialize();
	
	/* set dimensions */
	int numDOF = rEqnos.MinorDim() + 1; // additional DOF
	fContactEqnos.Allocate(fNumContactNodes*numDOF);
	fContactEqnos2D.Set(fNumContactNodes, numDOF, fContactEqnos.Pointer());
	
	/* allocate memory for force vector */
	fContactForce2D.Allocate(fNumContactNodes, numDOF);
	fContactForce.Set(fNumContactNodes*numDOF, fContactForce2D.Pointer());
	fContactForce2D = 0.0;

	/* register with node manager - sets initial fContactDOFtags */
	fXDOF_Nodes->Register(this, kNumAugLagDOF);	
}

void AugLagSphereT::SetEquationNumbers(void)
{
// don't need to set FBC destinations as with the base class because
// the class collects the global equation numbers during SendEqnsToSolver()
}

/* append element equations numbers to the list */
void AugLagSphereT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_2)

	/* dimensions */
	int ndof_u = rCoords.MinorDim();

	/* collect displacement DOF's */
	iArray2DT disp_eq(fContactNodes.Length(), ndof_u);
	NodeManagerT* nodemanager = fFEManager.NodeManager();
	nodemanager->SetLocalEqnos(fContactNodes, disp_eq);

	int eq_col = 0;
	iArrayT eq_temp(fContactNodes.Length());

	/* displacement equations */
	for (int i = 0; i < ndof_u; i++)
	{
		disp_eq.ColumnCopy(i, eq_temp);
		fContactEqnos2D.SetColumn(eq_col++, eq_temp);
	}

	/* constraint equations */
	const iArray2DT& auglageqs = fXDOF_Nodes->XDOF_Eqnos(this);
	for (int j = 0; j < auglageqs.MinorDim(); j++)
	{
		auglageqs.ColumnCopy(j, eq_temp);
		fContactEqnos2D.SetColumn(eq_col++, eq_temp);
	}

	/* send to solver */
	eq_1.Append(&fContactEqnos2D);
}

void AugLagSphereT::Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
#pragma unused(connects_2)
	connects_1.Append(&fContactTags);
}

void AugLagSphereT::ReadRestart(istream& in)
{
	/* inherited */
	PenaltySphereT::ReadRestart(in);

	in >> fLastDOF; // previous solution
}

void AugLagSphereT::WriteRestart(ostream& out) const
{
	/* inherited */
	PenaltySphereT::WriteRestart(out);

	out << fLastDOF; // previous solution
}

void AugLagSphereT::CloseStep(void)
{
	/* inherited */
	PenaltySphereT::CloseStep();

	/* store last converged DOF array */
	dArrayT constraints;
	constraints.Alias(fXDOF_Nodes->XDOF(this));
	fLastDOF = constraints;
}

/* restore the DOF values to the last converged solution */
void AugLagSphereT::ResetDOF(dArray2DT& DOF) const
{
	dArrayT constraints;
	constraints.Alias(DOF);
	constraints = fLastDOF;
}

/* tangent term */
void AugLagSphereT::ApplyLHS(void)
{
	/* time integration */
	double constK = 0.0;
	int formK = fController->FormK(constK);
	if (!formK) return;

	/* get current values of constraints */
	const dArray2DT& constr = fXDOF_Nodes->XDOF(this);
	const dArrayT force(constr.MajorDim(), constr.Pointer());

	/* workspace */
	int nsd = rCoords.MinorDim();
	dArrayT norm(nsd);
	dArrayT vec;
	dMatrixT ULblock(nsd);
	dMatrixT mat;

	/* node by node */
	for (int i = 0; i < fNumContactNodes; i++)
	{
		/* initialize */
		fLHS = 0.0;
	
		/* gap and augmented Lagrangian */
		double v = fDistances[i];
		double h = v - fRadius;
		double g = force[i] + fk*h;

		/* contact */
		if (g <= 0.0)
		{
			/* unit gap vector (from the force) */
			vec.Set(nsd, fContactForce2D(i));
			norm.SetToScaled(-1.0/g, vec);

			/* the long way */
			ULblock.Outer(norm, norm);
			ULblock *= (fk - g/v);
			ULblock.PlusIdentity(g/v);

			/* assemble element matrix */
			fLHS.SetBlock(0, 0, ULblock);
			
			mat.Set(1, nsd, norm.Pointer());
			fLHS.SetBlock(nsd, 0, mat);

			mat.Set(nsd, 1, norm.Pointer());
			fLHS.SetBlock(0, nsd, mat);
		}
		/* gap */
		else
		{
			/* augmented Lagrangian DOF */
			int dex = fLHS.Rows() - 1;
			fLHS(dex,dex) = -1.0/fk;							
		}

		/* time integration coefficient */
		fLHS *= constK;
		
		/* send to global equations */
		fContactEqnos2D.RowAlias(i,fi_sh);
		fFEManager.AssembleLHS(fLHS,fi_sh);
	}	
}

/* returns the array for the DOF tags needed for the current config */
iArrayT& AugLagSphereT::SetDOFTags(void)
{
// NOTE: this would be the place to determine the contact configuration
//       and collect the list of active nodes

	/* ALL constraints ALWAYS active */
	fContactDOFtags.Allocate(fContactNodes.Length());

	return fContactDOFtags;
}

const iArrayT& AugLagSphereT::DOFTags(void) const
{
	return fContactDOFtags;
}

/* generate nodal connectivities - does nothing here */
void AugLagSphereT::GenerateElementData(void)
{
	/* allocate space */
	fContactTags.Allocate(fContactNodes.Length(), 2);
	
	/* collect tags - {contact node, DOF tag} */
	fContactTags.SetColumn(0, fContactNodes);
	fContactTags.SetColumn(1, fContactDOFtags);
}

/* return the contact elements */
const iArray2DT& AugLagSphereT::DOFConnects(void) const
{
	return fContactTags;
}

/* returns 1 if group needs to reconfigure DOF's, else 0 */
int AugLagSphereT::Reconfigure(void) { return 0; }

/**********************************************************************
* Private
**********************************************************************/

/* accumulate the contact force vector fContactForce */
void AugLagSphereT::ComputeContactForce(double kforce)
{
	/* dimensions */
	int ndof_u = rCoords.MinorDim();
	int ndof   = fContactForce2D.MinorDim();

	/* initialize */
	fContactForce2D = 0.0;	

	/* get current values of constraints */
	const dArray2DT& constr = fXDOF_Nodes->XDOF(this);
	const dArrayT force(constr.MajorDim(), constr.Pointer());

	dArrayT f_u;
	for (int i = 0; i < fNumContactNodes; i++)
	{
		/* displacement DOF's */
		f_u.Set(ndof_u, fContactForce2D(i));
	
		/* center to striker */
		rCoords.RowCopy(fContactNodes[i], fv_OP);
		fv_OP -= fx;
		
		/* augmented Lagrangian multiplier */
		double v = fv_OP.Magnitude();
		double h = v - fRadius;
		double g = force[i] + fk*h;
	
		/* contact */
		if (g <= 0.0)
		{
			/* displace DOF's */
			f_u.SetToScaled(-g*kforce/v, fv_OP);

			/* augmented Lagrangian DOF */
			fContactForce2D(i, ndof - 1) = -h*kforce;
		}
		/* no contact */
		else
		{
			/* grad_disp contribution */
			f_u = 0.0;

			/* augmented Lagrangian DOF */
			fContactForce2D(i, ndof - 1) = force[i]*kforce/fk;			
		}

		//NOTE: This contact force is the negative of the element
		//      force in Heegaard.

		/* store */
		fDistances[i] = v;
	}
}
