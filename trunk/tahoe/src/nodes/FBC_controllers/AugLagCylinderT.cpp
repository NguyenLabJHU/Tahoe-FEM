/* $Id: AugLagCylinderT.cpp,v 1.1 2004-09-15 15:38:41 paklein Exp $ */
#include "AugLagCylinderT.h"
#include "FieldT.h"
#include "eIntegratorT.h"
#include "FieldSupportT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "XDOF_ManagerT.h"

/* for Uzawa parameters */
#include "AugLagWallT.h"

using namespace Tahoe;

/* parameters */
const int kNumAugLagDOF = 1;

/* constructor */
AugLagCylinderT::AugLagCylinderT(void):
	fUzawa(false),
	fPrimalIterations(-1),
	fPenetrationTolerance(-1.0),
	fRecomputeForce(false)
{
	SetName("cylinder_augmented_Lagrangian");
}

void AugLagCylinderT::SetEquationNumbers(void)
{
// don't need to set FBC destinations when not solving by Uzawa because
// the class collects the global equation numbers during SendEqnsToSolver()
}

/* append element equations numbers to the list */
void AugLagCylinderT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_2)

	/* Uzawa algorithm has local update */
	if (fUzawa)
		/* inherited */
		PenaltyCylinderT::Equations(eq_1, eq_2);
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
			cout << "\n AugLagCylinderT::Equations: node with constraint has prescribed DOF\n" 
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

void AugLagCylinderT::Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
	AutoArrayT<const iArray2DT*>& equivalent_nodes) const
{
#pragma unused(connects_2)
#pragma unused(equivalent_nodes)
	if (!fUzawa) connects_1.Append(&fContactTags);
}

void AugLagCylinderT::ReadRestart(istream& in)
{
	/* inherited */
	PenaltyCylinderT::ReadRestart(in);

	/* previous solution */
	if (fUzawa)
		in >> fDOF;
	else
		in >> fLastDOF;
}

void AugLagCylinderT::WriteRestart(ostream& out) const
{
	/* inherited */
	PenaltyCylinderT::WriteRestart(out);

	/* previous solution */
	if (fUzawa)
		out << fDOF;
	else
		out << fLastDOF;
}

void AugLagCylinderT::InitStep(void)
{
	/* inherited */
	PenaltyCylinderT::InitStep();

	/* store solution */
	if (fUzawa) fLastDOF = fDOF;
}

void AugLagCylinderT::CloseStep(void)
{
	/* inherited */
	PenaltyCylinderT::CloseStep();

	/* store last converged DOF array */
	if (!fUzawa) {
		dArrayT constraints;
		constraints.Alias(FieldSupport().XDOF_Manager().XDOF(this, 0));
		fLastDOF = constraints;
	}
}

/* update constrain forces */
GlobalT::RelaxCodeT AugLagCylinderT::RelaxSystem(void)
{
	GlobalT::RelaxCodeT relax = PenaltyCylinderT::RelaxSystem();
	
	/* check penetration tolerance */
	if (fUzawa)
	{
		/* evaluate constraints */
		const dArray2DT& coords = FieldSupport().CurrentCoordinates();
		double penetration_norm = 0.0;
		for (int i = 0; i < fNumContactNodes; i++)
		{
			/* center to striker */
			coords.RowCopy(fContactNodes[i], fv_OP);
			fv_OP -= fx;

			/* vector in radial direction */
			fR.SetToCombination(1.0, fv_OP, -dArrayT::Dot(fv_OP, fDirection), fDirection);

			/* augmented Lagrangian multiplier */
			double dist = fR.Magnitude();
			double h = dist - fRadius;
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
void AugLagCylinderT::ResetDOF(dArray2DT& DOF, int tag_set) const
{
#pragma unused (tag_set)

	dArrayT constraints;
	constraints.Alias(DOF);
	constraints = fLastDOF;
}

/* returns the array for the DOF tags needed for the current config */
void AugLagCylinderT::SetDOFTags(void)
{
// NOTE: this would be the place to determine the contact configuration
//       and collect the list of active nodes

	/* ALL constraints ALWAYS active */
	fContactDOFtags.Dimension(fContactNodes.Length());
}

iArrayT& AugLagCylinderT::DOFTags(int tag_set)
{
#pragma unused (tag_set)
	return fContactDOFtags;
}

/* generate nodal connectivities - does nothing here */
void AugLagCylinderT::GenerateElementData(void)
{
	/* allocate space */
	fContactTags.Dimension(fContactNodes.Length(), 2);
	
	/* collect tags - {contact node, DOF tag} */
	fContactTags.SetColumn(0, fContactNodes);
	fContactTags.SetColumn(1, fContactDOFtags);
}

/* return the contact elements */
const iArray2DT& AugLagCylinderT::DOFConnects(int tag_set) const
{
#pragma unused (tag_set)
	return fContactTags;
}

/* returns 1 if group needs to reconfigure DOF's, else 0 */
int AugLagCylinderT::Reconfigure(void) { return 0; }

/* return the equation group */
int AugLagCylinderT::Group(void) const { return Field().Group(); };

/* tangent term */
void AugLagCylinderT::ApplyLHS(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)

	/* no stiffness contribution with Uzawa */
	if (fUzawa) return;

	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	/* equations */
	const iArray2DT& eqnos = Field().Equations();

	/* support class */
	const FieldSupportT& support = FieldSupport();

	/* get current values of constraints */
	const dArray2DT& constr = support.XDOF_Manager().XDOF(this, 0);
	const dArrayT force(constr.MajorDim(), constr.Pointer());

	/* workspace */
	int nsd = support.NumSD();
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
		double gap  = fGap[i];
		double dist = gap + fRadius;
		double g = force[i] + fk*gap;		

		/* active */
		if (g <= 0.0)
		{
			/* get force vector (has normal direction) */
			fContactForce2D.RowAlias(i, fd_sh);
			
			/* unit normal */
			norm.SetToScaled(-1.0/g, fd_sh);

			/* upper-left block */			
			double dPhi = gap*fk;
			ULblock.Outer(fd_sh, fd_sh, constK*((fk/dPhi) - (1.0/dist))/dPhi);
			ULblock.PlusIdentity(constK*dPhi/dist);
			ULblock.Outer(fDirection, fDirection, -constK*dPhi/dist, dMatrixT::kAccumulate);		
			fLHS.SetBlock(0, 0, ULblock);

			mat.Alias(1, nsd, norm.Pointer());
			fLHS.SetBlock(nsd, 0, mat);

			mat.Alias(nsd, 1, norm.Pointer());
			fLHS.SetBlock(0, nsd, mat);

			/* augmented Lagrangian DOF */
			if (fFloatingDOF[i] && fabs(g) < kSmall) {
				int dex = fLHS.Rows() - 1;
				fLHS(dex,dex) = -1.0/fk;							
			}
		}
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

/* information about subordinate parameter lists */
void AugLagCylinderT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	PenaltyCylinderT::DefineSubs(sub_list);

	/* direction */
	sub_list.AddSub("Uzawa_method");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* AugLagCylinderT::NewSub(const StringT& name) const
{
	if (name == "Uzawa_method")
	{
		/* parameters defined by wall */
		AugLagWallT wall;
		return wall.NewSub(name);
	}
	else /* inherited */
		return PenaltyCylinderT::NewSub(name);
}

/* accept parameter list */
void AugLagCylinderT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "AugLagCylinderT::TakeParameterList";

	/* inherited */
	PenaltyCylinderT::TakeParameterList(list);

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
 * Protected
 **********************************************************************/

/* compute the nodal contribution to the residual force vector */
void AugLagCylinderT::ComputeContactForce(double kforce)
{
	/* Uzawa */
	if (fUzawa)
	{

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

		/* loop over strikers */
		const dArray2DT& coords = FieldSupport().CurrentCoordinates();
		dArrayT f_u;
		for (int i = 0; i < fNumContactNodes; i++)
		{
			/* displacement DOF's */
			f_u.Set(ndof_u, fContactForce2D(i));

			/* center to striker */
			coords.RowCopy(fContactNodes[i], fv_OP);
			fv_OP -= fx;

			/* vector in radial direction */
			fR.SetToCombination(1.0, fv_OP, -dArrayT::Dot(fv_OP, fDirection), fDirection);
		
			/* penetration */
			double dist = fR.Magnitude();
			double pen  = dist - fRadius;
			double g = force[i] + fk*pen;			
			if (g <= 0.0)
			{
				/* displace DOF's */
				f_u.SetToScaled(-g*kforce/dist, fR);

				/* augmented Lagrangian DOF */
				fContactForce2D(i, ndof - 1) = -pen*kforce;
			}
			else
			{
				/* grad_disp contribution */
				f_u = 0.0;

				/* augmented Lagrangian DOF */
				fContactForce2D(i, ndof - 1) = force[i]*kforce/fk;			
			}

			/* store */
			fGap[i] = pen;
		}
	}
}
