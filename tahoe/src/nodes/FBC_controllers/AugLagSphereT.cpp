/* $Id: AugLagSphereT.cpp,v 1.11.28.1 2004-04-08 07:33:51 paklein Exp $ */
/* created: paklein (03/24/1999) */
#include "AugLagSphereT.h"
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
AugLagSphereT::AugLagSphereT(void)
{
	SetName("sphere_augmented_Lagrangian");
}

#if 0
/* initialize data */
void AugLagSphereT::Initialize(void)
{
	/* inherited */
	PenaltySphereT::Initialize();
	
	/* set dimensions */
	int numDOF = rEqnos.MinorDim() + 1; // additional DOF
	fContactEqnos.Dimension(fNumContactNodes*numDOF);
	fContactEqnos2D.Set(fNumContactNodes, numDOF, fContactEqnos.Pointer());
	
	/* allocate memory for force vector */
	fContactForce2D.Dimension(fNumContactNodes, numDOF);
	fContactForce.Set(fNumContactNodes*numDOF, fContactForce2D.Pointer());
	fContactForce2D = 0.0;

	/* register with node manager - sets initial fContactDOFtags */
	iArrayT set_dims(1);
	set_dims = kNumAugLagDOF;
	fXDOF_Nodes->XDOF_Register(this, set_dims);	
}
#endif

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
	constraints.Alias(FieldSupport().XDOF_Manager().XDOF(this, 0));
	fLastDOF = constraints;
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
	for (int i = 0; i < fNumContactNodes; i++)
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

/* accept parameter list */
void AugLagSphereT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	PenaltySphereT::TakeParameterList(list);

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

/**********************************************************************
* Private
**********************************************************************/

/* accumulate the contact force vector fContactForce */
void AugLagSphereT::ComputeContactForce(double kforce)
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
	for (int i = 0; i < fNumContactNodes; i++)
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
