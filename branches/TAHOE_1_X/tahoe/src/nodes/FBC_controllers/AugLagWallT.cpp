/* $Id: AugLagWallT.cpp,v 1.11 2003-10-04 19:14:05 paklein Exp $ */
#include "AugLagWallT.h"

#include <iostream.h>
#include <iomanip.h>

#include "toolboxConstants.h"
#include "FEManagerT.h"
#include "XDOF_ManagerT.h"
#include "eIntegratorT.h"
#include "FieldT.h"

using namespace Tahoe;

/* parameters */
const int kNumAugLagDOF = 1;

/* constructor */
AugLagWallT::AugLagWallT(FEManagerT& fe_manager, XDOF_ManagerT* XDOF_nodes,
	const FieldT& field, const dArray2DT& coords, const dArray2DT& disp):
	PenaltyWallT(fe_manager, 
		field.Group(),
		field.Equations(), 
		coords,
		disp,
		(field.Order() > 0) ? &(field[1]): NULL),
	fXDOF_Nodes(XDOF_nodes),
	fField(field)
{
	/* (re-)dimension the tangent matrix */
	fLHS.Dimension(rEqnos.MinorDim() + 1); // additional DOF
}

/* initialize data */
void AugLagWallT::Initialize(void)
{
	/* inherited */
	PenaltyWallT::Initialize();
	
	/* set dimensions */
	int numDOF = rEqnos.MinorDim() + 1; // additional DOF
	fContactEqnos.Dimension(fNumContactNodes*numDOF);
	fContactEqnos2D.Set(fNumContactNodes, numDOF, fContactEqnos.Pointer());
	fFloatingDOF.Dimension(fNumContactNodes);
	fFloatingDOF = 0;
	
	/* allocate memory for force vector */
	fContactForce2D.Dimension(fNumContactNodes, numDOF);
	fContactForce.Set(fNumContactNodes*numDOF, fContactForce2D.Pointer());
	fContactForce2D = 0.0;

	/* register with node manager - sets initial fContactDOFtags */
	iArrayT set_dims(1);
	set_dims = kNumAugLagDOF;
	fXDOF_Nodes->XDOF_Register(this, set_dims);	
}

void AugLagWallT::SetEquationNumbers(void)
{
// don't need to set FBC destinations as with the base class because
// the class collects the global equation numbers during SendEqnsToSolver()
}

/* append element equations numbers to the list */
void AugLagWallT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_2)

	/* dimensions */
	int ndof_u = rCoords.MinorDim();

	/* collect displacement DOF's */
	iArray2DT disp_eq(fContactNodes.Length(), ndof_u);
	fField.SetLocalEqnos(fContactNodes, disp_eq);

	int eq_col = 0;
	iArrayT eq_temp(fContactNodes.Length());

	/* displacement equations */
	fFloatingDOF = 0;
	for (int i = 0; i < ndof_u; i++)
	{
		disp_eq.ColumnCopy(i, eq_temp);
		fContactEqnos2D.SetColumn(eq_col++, eq_temp);

		/* check for floating DOF's */
		for (int j = 0; j < eq_temp.Length(); j++)
			if (eq_temp[j] < 1)
				fFloatingDOF[j] = 1; /* mark */
	}

	/* warning */
	if (fFloatingDOF.HasValue(1))
		cout << "\n AugLagWallT::Equations: node with constraint has prescribed DOF\n" 
		     <<   "     Stiffness may be approximate." << endl;	

	/* constraint equations */
	const iArray2DT& auglageqs = fXDOF_Nodes->XDOF_Eqnos(this, 0);
	for (int j = 0; j < auglageqs.MinorDim(); j++)
	{
		auglageqs.ColumnCopy(j, eq_temp);
		fContactEqnos2D.SetColumn(eq_col++, eq_temp);
	}

	/* send to solver */
	eq_1.Append(&fContactEqnos2D);
}

void AugLagWallT::Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
	AutoArrayT<const iArray2DT*>& equivalent_nodes) const
{
#pragma unused(connects_2)
#pragma unused(equivalent_nodes)
	connects_1.Append(&fContactTags);
}

void AugLagWallT::ReadRestart(istream& in)
{
	/* inherited */
	PenaltyWallT::ReadRestart(in);

	in >> fLastDOF; // previous solution
}

void AugLagWallT::WriteRestart(ostream& out) const
{
	/* inherited */
	PenaltyWallT::WriteRestart(out);

	out << fLastDOF; // previous solution
}

void AugLagWallT::CloseStep(void)
{
	/* inherited */
	PenaltyWallT::CloseStep();

	/* store last converged DOF array */
	dArrayT constraints;
	constraints.Alias(fXDOF_Nodes->XDOF(this, 0));
	fLastDOF = constraints;
}

/* restore the DOF values to the last converged solution */
void AugLagWallT::ResetDOF(dArray2DT& DOF, int tag_set) const
{
#pragma unused (tag_set)

	dArrayT constraints;
	constraints.Alias(DOF);
	constraints = fLastDOF;
}

/* tangent term */
void AugLagWallT::ApplyLHS(GlobalT::SystemTypeT sys_type)
{
#pragma unused (sys_type)

	/* time integration */
	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	/* get current values of constraints */
	const dArray2DT& constr = fXDOF_Nodes->XDOF(this, 0);
	const dArrayT force(constr.MajorDim(), constr.Pointer());

	/* workspace */
	int nsd = rCoords.MinorDim();
	dArrayT vec;
	dMatrixT ULblock(nsd);
	dMatrixT mat;

	/* node by node */
	for (int i = 0; i < fNumContactNodes; i++)
	{
		/* initialize */
		fLHS = 0.0;
	
		/* gap and augmented Lagrangian */
		double h = fp_i.DotRow(i, fnormal); /* fp_i set during ComputeContactForce */
		double g = force[i] + fk*h;

		/* contact */
		if (g <= 0.0)
		{
			/* the long way */
			ULblock.Outer(fnormal, fnormal);
			ULblock *= fk;

			/* assemble element matrix */
			fLHS.SetBlock(0, 0, ULblock);
			
			mat.Set(1, nsd, fnormal.Pointer());
			fLHS.SetBlock(nsd, 0, mat);

			mat.Set(nsd, 1, fnormal.Pointer());
			fLHS.SetBlock(0, nsd, mat);

			/* augmented Lagrangian DOF */
			if (fFloatingDOF[i] && fabs(g) < kSmall) {
				int dex = fLHS.Rows() - 1;
				fLHS(dex,dex) = -1.0/fk;							
			}
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
		fFEManager.AssembleLHS(fGroup, fLHS, fi_sh);
	}	
}

/* returns the array for the DOF tags needed for the current config */
void AugLagWallT::SetDOFTags(void)
{
// NOTE: this would be the place to determine the contact configuration
//       and collect the list of active nodes

	/* ALL constraints ALWAYS active */
	fContactDOFtags.Dimension(fContactNodes.Length());
}

iArrayT& AugLagWallT::DOFTags(int tag_set)
{
#pragma unused (tag_set)
	return fContactDOFtags;
}

/* generate nodal connectivities - does nothing here */
void AugLagWallT::GenerateElementData(void)
{
	/* allocate space */
	fContactTags.Dimension(fContactNodes.Length(), 2);
	
	/* collect tags - {contact node, DOF tag} */
	fContactTags.SetColumn(0, fContactNodes);
	fContactTags.SetColumn(1, fContactDOFtags);
}

/* return the contact elements */
const iArray2DT& AugLagWallT::DOFConnects(int tag_set) const
{
#pragma unused (tag_set)
	return fContactTags;
}

/* returns 1 if group needs to reconfigure DOF's, else 0 */
int AugLagWallT::Reconfigure(void) { return 0; }

/* return the equation group */
int AugLagWallT::Group(void) const { return fField.Group(); };

/**********************************************************************
* Private
**********************************************************************/

/* accumulate the contact force vector fContactForce */
void AugLagWallT::ComputeContactForce(double kforce)
{
	/* dimensions */
	int ndof_u = rCoords.MinorDim();
	int ndof   = fContactForce2D.MinorDim();

	/* initialize */
	fContactForce2D = 0.0;	

	/* get current values of constraints */
	const dArray2DT& constr = fXDOF_Nodes->XDOF(this, 0);
	const dArrayT force(constr.MajorDim(), constr.Pointer());

	/* compute relative positions */
	fp_i.RowCollect(fContactNodes, rCoords);
	for (int j = 0; j < fNumContactNodes; j++)
		fp_i.AddToRowScaled(j, -1.0, fx);

	dArrayT f_u;
	for (int i = 0; i < fNumContactNodes; i++)
	{
		/* displacement DOF's */
		f_u.Set(ndof_u, fContactForce2D(i));
	
		/* augmented Lagrangian multiplier */
		double h = fp_i.DotRow(i, fnormal);
		double g = force[i] + fk*h;
	
		/* contact */
		if (g <= 0.0)
		{
			/* displace DOF's */
			f_u.SetToScaled(-g*kforce, fnormal);

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
		
		/* store gap */
		fGap[i] = h;
	}
}
