/* $Id: AugLagWallT.cpp,v 1.16 2004-10-14 20:25:02 paklein Exp $ */
#include "AugLagWallT.h"
#include "FieldT.h"
#include "eIntegratorT.h"
#include "FieldSupportT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "XDOF_ManagerT.h"

using namespace Tahoe;

/* parameters */
const int kNumAugLagDOF = 1;

/* constructor */
AugLagWallT::AugLagWallT(void):
	fUzawa(false),
	fPrimalIterations(-1),
	fPenetrationTolerance(-1.0),
	fRecomputeForce(false)
{
	SetName("wall_augmented_Lagrangian");
}

void AugLagWallT::SetEquationNumbers(void)
{
// don't need to set FBC destinations when not solving by Uzawa because
// the class collects the global equation numbers during SendEqnsToSolver()
}

/* append element equations numbers to the list */
void AugLagWallT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_2)

	/* Uzawa algorithm has local update */
	if (fUzawa)
		/* inherited */
		PenaltyWallT::Equations(eq_1, eq_2);
	else
	{
		/* dimensions */
		int ndof_u = Field().NumDOF();

		/* collect displacement DOF's */
		iArray2DT disp_eq(fContactNodes.Length(), ndof_u);
		Field().SetLocalEqnos(fContactNodes, disp_eq);

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
		const iArray2DT& auglageqs = FieldSupport().XDOF_Manager().XDOF_Eqnos(this, 0);
		for (int j = 0; j < auglageqs.MinorDim(); j++)
		{
			auglageqs.ColumnCopy(j, eq_temp);
			fContactEqnos2D.SetColumn(eq_col++, eq_temp);
		}

		/* send to solver */
		eq_1.Append(&fContactEqnos2D);
	}
}

void AugLagWallT::Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
	AutoArrayT<const iArray2DT*>& equivalent_nodes) const
{
#pragma unused(connects_2)
#pragma unused(equivalent_nodes)
	if (!fUzawa) connects_1.Append(&fContactTags);
}

void AugLagWallT::ReadRestart(istream& in)
{
	/* inherited */
	PenaltyWallT::ReadRestart(in);

	/* previous solution */
	if (fUzawa)
		in >> fDOF;
	else
		in >> fLastDOF;
}

void AugLagWallT::WriteRestart(ostream& out) const
{
	/* inherited */
	PenaltyWallT::WriteRestart(out);

	/* previous solution */
	if (fUzawa)
		out << fDOF;
	else
		out << fLastDOF;
}

void AugLagWallT::InitStep(void)
{
	/* inherited */
	PenaltyWallT::InitStep();

	/* store solution */
	if (fUzawa) fLastDOF = fDOF;
}

void AugLagWallT::CloseStep(void)
{
	/* inherited */
	PenaltyWallT::CloseStep();

	/* store last converged DOF array */
	if (!fUzawa) {
		dArrayT constraints;
		constraints.Alias(FieldSupport().XDOF_Manager().XDOF(this, 0));
		fLastDOF = constraints;
	}
}

/* update constrain forces */
GlobalT::RelaxCodeT AugLagWallT::RelaxSystem(void)
{
	GlobalT::RelaxCodeT relax = PenaltyWallT::RelaxSystem();
	
	/* check penetration tolerance */
	if (fUzawa)
	{
		/* compute relative positions */
		const dArray2DT& coords = FieldSupport().CurrentCoordinates();
		fp_i.RowCollect(fContactNodes, coords);
		for (int j = 0; j < fNumContactNodes; j++)
			fp_i.AddToRowScaled(j, -1.0, fx);

		/* evaluate constraints */
		double penetration_norm = 0.0;
		for (int i = 0; i < fNumContactNodes; i++)
		{
			/* augmented Lagrangian multiplier */
			double h = fp_i.DotRow(i, fnormal);
			double g = fDOF[i] + fk*h;
	
			/* contact */
			if (g <= 0.0) penetration_norm += h*h;

			/* update stored gap */
			fGap[i] = h;
		}
		
		/* check convergence */
		if (sqrt(penetration_norm) > fPenetrationTolerance) {
			relax = GlobalT::MaxPrecedence(relax, GlobalT::kRelax);
			fRecomputeForce = true;
		}
		else
			relax = GlobalT::MaxPrecedence(relax, GlobalT::kNoRelax);
	}

	/* return */
	return relax;
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

	/* no stiffness contribution with Uzawa */
	if (fUzawa) return;

	/* time integration */
	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	/* get current values of constraints */
	const dArray2DT& constr = FieldSupport().XDOF_Manager().XDOF(this, 0);
	const dArrayT force(constr.MajorDim(), constr.Pointer());

	/* workspace */
	int nsd = FieldSupport().NumSD();
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
			/* upper-left block */
			ULblock.Outer(fnormal, fnormal, fk);

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
		FieldSupport().AssembleLHS(fGroup, fLHS, fi_sh);
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
int AugLagWallT::Group(void) const { return Field().Group(); };

/* information about subordinate parameter lists */
void AugLagWallT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	PenaltyWallT::DefineSubs(sub_list);

	/* parameters for Uzawa's method */
	sub_list.AddSub("Uzawa_method", ParameterListT::ZeroOrOnce);
}
	
/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* AugLagWallT::NewSub(const StringT& name) const
{
	if (name == "Uzawa_method")
	{
		ParameterContainerT* Uzawa = new ParameterContainerT(name);
	
		ParameterT primal_its(fPrimalIterations, "primal_iterations");
		primal_its.AddLimit(1, LimitT::LowerInclusive);
		Uzawa->AddParameter(primal_its);

		ParameterT pen_tolerance(fPenetrationTolerance, "penetration_tolerance");
		pen_tolerance.AddLimit(0, LimitT::Lower);
		Uzawa->AddParameter(pen_tolerance);		
	
		return Uzawa;
	}
	else /* inherited */
		return PenaltyWallT::NewSub(name);
}

/* accept parameter list */
void AugLagWallT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	PenaltyWallT::TakeParameterList(list);

	/* look for Uzawa parameters */
	const ParameterListT* Uzawa = list.List("Uzawa_method");
	if (Uzawa)
	{
		fUzawa = true;
		fPrimalIterations = Uzawa->GetParameter("primal_iterations");
		fPenetrationTolerance = Uzawa->GetParameter("penetration_tolerance");
	}
	else
		fUzawa = false;

	/* do Uzawa iterations or solve concurrently */
	if (!fUzawa)
	{
		/* (re-)dimension the tangent matrix */
		int ndof = Field().NumDOF() + 1; // additional DOF
		fLHS.Dimension(ndof); 
	
		/* set dimensions */
		fContactEqnos.Dimension(fNumContactNodes*ndof);
		fContactEqnos2D.Set(fNumContactNodes, ndof, fContactEqnos.Pointer());
		fFloatingDOF.Dimension(fNumContactNodes);
		fFloatingDOF = 0;
	
		/* allocate memory for force vector */
		fContactForce2D.Dimension(fNumContactNodes, ndof);
		fContactForce.Set(fNumContactNodes*ndof, fContactForce2D.Pointer());
		fContactForce2D = 0.0;

		/* register with node manager - sets initial fContactDOFtags */
		iArrayT set_dims(1);
		set_dims = kNumAugLagDOF;
		FieldSupport().XDOF_Manager().XDOF_Register(this, set_dims);	
	}
	else /* allocated space for multipliers */ {
		fDOF.Dimension(fNumContactNodes);
		fDOF = 0.0;
	}
}

/**********************************************************************
* Private
**********************************************************************/

/* accumulate the contact force vector fContactForce */
void AugLagWallT::ComputeContactForce(double kforce)
{
	/* Uzawa */
	if (fUzawa)
	{
		/* check for update */
		int iter = FieldSupport().IterationNumber();
		if (!fRecomputeForce && (iter != -1 && (iter+1)%fPrimalIterations != 0)) return;
		fRecomputeForce = false;
	
		/* dimensions */
		int ndof = Field().NumDOF();

		/* initialize */
		fContactForce2D = 0.0;	

		/* compute relative positions */
		const dArray2DT& coords = FieldSupport().CurrentCoordinates();
		fp_i.RowCollect(fContactNodes, coords);
		for (int j = 0; j < fNumContactNodes; j++)
			fp_i.AddToRowScaled(j, -1.0, fx);

		dArrayT f_u;
		for (int i = 0; i < fNumContactNodes; i++)
		{
			/* displacement DOF's */
			f_u.Set(ndof, fContactForce2D(i));
	
			/* augmented Lagrangian multiplier */
			double h = fp_i.DotRow(i, fnormal);
			double g = fDOF[i] + fk*h;
	
			/* contact */
			if (g <= 0.0) {
				f_u.SetToScaled(-g*kforce, fnormal);
				fDOF[i] = g;
			}
			else
				fDOF[i] = 0.0;
		
			/* store gap */
			fGap[i] = h;
		}
	}
	else /* solve primal and dual simultaneously */
	{
		/* dimensions */
		int ndof_u = Field().NumDOF();
		int ndof   = fContactForce2D.MinorDim();

		/* initialize */
		fContactForce2D = 0.0;	

		/* get current values of constraints */
		const dArray2DT& constr = FieldSupport().XDOF_Manager().XDOF(this, 0);
		const dArrayT force(constr.MajorDim(), constr.Pointer());
		const dArray2DT& coords = FieldSupport().CurrentCoordinates();
	
		/* compute relative positions */
		fp_i.RowCollect(fContactNodes, coords);
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
}
