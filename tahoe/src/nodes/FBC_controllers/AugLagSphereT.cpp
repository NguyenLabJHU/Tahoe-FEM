/* $Id: AugLagSphereT.cpp,v 1.16.14.1 2005-07-02 22:50:38 paklein Exp $ */
/* created: paklein (03/24/1999) */
#include "AugLagSphereT.h"
#include "FieldT.h"
#include "eIntegratorT.h"
#include "FieldSupportT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "XDOF_ManagerT.h"
#include "AugLagWallT.h"

using namespace Tahoe;

/* parameters */
const int kNumAugLagDOF = 1;

/* constructor */
AugLagSphereT::AugLagSphereT(void):
	fUzawa(false),
	fPrimalIterations(-1),
	fPenetrationTolerance(-1.0),
	fRecomputeForce(false),
	fIterationi(-2)
{
	SetName("sphere_augmented_Lagrangian");
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
	/* Uzawa algorithm has local update */
	if (fUzawa)
		/* inherited */
		PenaltySphereT::Equations(eq_1, eq_2);
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
		for (int i = 0; i < ndof_u; i++)
		{
			disp_eq.ColumnCopy(i, eq_temp);
			fContactEqnos2D.SetColumn(eq_col++, eq_temp);
		}

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

void AugLagSphereT::Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
	AutoArrayT<const iArray2DT*>& equivalent_nodes) const
{
#pragma unused(connects_2)
#pragma unused(equivalent_nodes)
	connects_1.Append(&fContactTags);
}

void AugLagSphereT::ReadRestart(istream& in)
{
	/* inherited */
	PenaltySphereT::ReadRestart(in);

	/* previous solution */
	if (fUzawa)
		in >> fDOF;
	else
		in >> fLastDOF;
}

void AugLagSphereT::WriteRestart(ostream& out) const
{
	/* inherited */
	PenaltySphereT::WriteRestart(out);

	/* previous solution */
	if (fUzawa)
		out << fDOF;
	else
		out << fLastDOF;
}

void AugLagSphereT::InitStep(void)
{
	/* inherited */
	PenaltySphereT::InitStep();

	/* store solution */
	if (fUzawa) {
		fLastDOF = fDOF;
		fIterationi = -2;
	}
}

void AugLagSphereT::CloseStep(void)
{
	/* inherited */
	PenaltySphereT::CloseStep();

	/* store last converged DOF array */
	if (!fUzawa) {
		dArrayT constraints;
		constraints.Alias(FieldSupport().XDOF_Manager().XDOF(this, 0));
		fLastDOF = constraints;
	}
}

/* update constrain forces */
GlobalT::RelaxCodeT AugLagSphereT::RelaxSystem(void)
{
	GlobalT::RelaxCodeT relax = PenaltySphereT::RelaxSystem();
	
	/* check penetration tolerance */
	if (fUzawa)
	{
		/* evaluate constraints */
		const dArray2DT& coords = FieldSupport().CurrentCoordinates();
		double penetration_norm = 0.0;
		for (int i = 0; i < fContactNodes.Length(); i++)
		{
			/* center to striker */
			coords.RowCopy(fContactNodes[i], fv_OP);
			fv_OP -= fx;

			/* center to striker */
			coords.RowCopy(fContactNodes[i], fv_OP);
			fv_OP -= fx;
		
			/* augmented Lagrangian multiplier */
			double v = fv_OP.Magnitude();
			double h = v - fRadius;
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
void AugLagSphereT::ResetDOF(dArray2DT& DOF, int tag_set) const
{
#pragma unused (tag_set)

	dArrayT constraints;
	constraints.Alias(DOF);
	constraints = fLastDOF;
}

/* tangent term */
void AugLagSphereT::ApplyLHS(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)

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
	dArrayT norm(nsd);
	dArrayT vec;
	dMatrixT ULblock(nsd);
	dMatrixT mat;

	/* node by node */
	for (int i = 0; i < fContactNodes.Length(); i++)
	{
		/* initialize */
		fLHS = 0.0;
	
		/* gap and augmented Lagrangian */
		double h = fGap[i];
		double v = h + fRadius;
		double g = force[i] + fk*h;

		/* contact */
		if (g <= 0.0)
		{
			/* unit gap vector (from the force) */
			vec.Alias(nsd, fContactForce2D(i));			
			norm.SetToScaled(-1.0/g, vec);

			/* the long way */
			ULblock.Outer(norm, norm, fk - g/v);
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
		FieldSupport().AssembleLHS(fGroup, fLHS,fi_sh);
	}	
}

/* returns the array for the DOF tags needed for the current config */
void AugLagSphereT::SetDOFTags(void)
{
// NOTE: this would be the place to determine the contact configuration
//       and collect the list of active nodes

	/* ALL constraints ALWAYS active */
	fContactDOFtags.Dimension(fContactNodes.Length());
}

iArrayT& AugLagSphereT::DOFTags(int tag_set)
{
#pragma unused (tag_set)
	return fContactDOFtags;
}

/* generate nodal connectivities - does nothing here */
void AugLagSphereT::GenerateElementData(void)
{
	/* allocate space */
	fContactTags.Dimension(fContactNodes.Length(), 2);
	
	/* collect tags - {contact node, DOF tag} */
	fContactTags.SetColumn(0, fContactNodes);
	fContactTags.SetColumn(1, fContactDOFtags);
}

/* return the contact elements */
const iArray2DT& AugLagSphereT::DOFConnects(int tag_set) const
{
#pragma unused (tag_set)
	return fContactTags;
}

/* returns 1 if group needs to reconfigure DOF's, else 0 */
int AugLagSphereT::Reconfigure(void) { return 0; }

/* return the equation group */
int AugLagSphereT::Group(void) const { return Field().Group(); };

/* information about subordinate parameter lists */
void AugLagSphereT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	PenaltySphereT::DefineSubs(sub_list);

	/* direction */
	sub_list.AddSub("Uzawa_method", ParameterListT::ZeroOrOnce);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* AugLagSphereT::NewSub(const StringT& name) const
{
	if (name == "Uzawa_method")
	{
		/* parameters defined by wall */
		AugLagWallT wall;
		return wall.NewSub(name);
	}
	else /* inherited */
		return PenaltySphereT::NewSub(name);
}

/* accept parameter list */
void AugLagSphereT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	PenaltySphereT::TakeParameterList(list);

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

ExceptionT::GeneralFail("AugLagSphereT::TakeParameterList", "not up to date");
#if 0
	/* do Uzawa iterations or solve concurrently */
	if (!fUzawa)
	{
		/* (re-)dimension the tangent matrix */
		int ndof = Field().NumDOF() + 1; // additional DOF
		fLHS.Dimension(ndof); 

		/* set dimensions */
		fContactEqnos.Dimension(fNumContactNodes*ndof);
		fContactEqnos2D.Set(fNumContactNodes, ndof, fContactEqnos.Pointer());
	
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
#endif
}

/**********************************************************************
* Private
**********************************************************************/

/* accumulate the contact force vector fContactForce */
void AugLagSphereT::ComputeContactForce(double kforce)
{
	/* Uzawa */
	if (fUzawa)
	{
		/* check for update */
		int iter = FieldSupport().IterationNumber();
		if (!fRecomputeForce && (iter != -1 && (iter+1)%fPrimalIterations != 0)) return;
		fRecomputeForce = false;
	
		/* store last iterate */
		if (fIterationi != iter) {
			fIterationi = iter;
			fDOFi = fDOF;
		}
	
		/* dimensions */
		int ndof = Field().NumDOF();

		/* initialize */
		fContactForce2D = 0.0;	

		/* loop over strikers */
		const dArray2DT& coords = FieldSupport().CurrentCoordinates();
		dArrayT f_u;
		for (int i = 0; i < fContactNodes.Length(); i++)
		{
			/* displacement DOF's */
			f_u.Set(ndof, fContactForce2D(i));
	
			/* center to striker */
			coords.RowCopy(fContactNodes[i], fv_OP);
			fv_OP -= fx;
		
			/* augmented Lagrangian multiplier */
			double v = fv_OP.Magnitude();
			double h = v - fRadius;
			double g = fDOFi[i] + fk*h;
	
			/* contact */
			if (g <= 0.0) {
				f_u.SetToScaled(-g*kforce/v, fv_OP);
				fDOF[i] = g;
			}
			/* no contact */
			else
				fDOF[i] = 0.0;

			/* store */
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
		dArrayT f_u;
		for (int i = 0; i < fContactNodes.Length(); i++)
		{
			/* displacement DOF's */
			f_u.Set(ndof_u, fContactForce2D(i));
	
			/* center to striker */
			coords.RowCopy(fContactNodes[i], fv_OP);
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
			fGap[i] = h;
		}
	}
}
