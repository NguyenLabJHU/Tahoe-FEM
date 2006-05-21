/* $Id: RigidCSEAnisoT.cpp,v 1.5 2006-05-21 17:47:59 paklein Exp $ */
#include "RigidCSEAnisoT.h"

#include "XDOF_ManagerT.h"
#include "ifstreamT.h"
#include "eIntegratorT.h"
#include "SurfaceShapeT.h"

using namespace Tahoe;

/* constructor */
RigidCSEAnisoT::RigidCSEAnisoT(const ElementSupportT& support):
	CSEAnisoT(support),
	fr(0.0),
	fDisp(LocalArrayT::kDisp),
	fLastDisp(LocalArrayT::kLastDisp),
	fConstraintShapes(NULL)
{
	SetName("rigid_anisotropic_CSE");
}

/* destructor */
RigidCSEAnisoT::~RigidCSEAnisoT(void)
{
	if (fConstraintShapes != fShapes)
		delete fConstraintShapes;
}

/* append element equations numbers to the list */
void RigidCSEAnisoT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	const char caller[] = "RigidCSEAnisoT::Equations";

	/* inherited */
	CSEAnisoT::Equations(eq_1, eq_2);

	/* source for XDOF information */
	XDOF_ManagerT& xdof_manager = ElementSupport().XDOF_Manager();

	/* resize the work space */
	int num_constraints = fConstraintXDOFTags.Length();
	fXDOFEqnos_man.SetMajorDimension(num_constraints, false);

	/* get the equations associate with the XDOF */
	const iArray2DT& xdof_eqnos = xdof_manager.XDOF_Eqnos(this, 0);
	if (xdof_eqnos.Length() != fConstraintXDOFTags.Length())
		ExceptionT::GeneralFail(caller, "expecting %d xdof equations not %d",
			fConstraintXDOFTags.Length(), xdof_eqnos.Length());
	
	/* collect equations for active constraints */
	int constraint_dex = 0;
	int nel = fConstraintStatus.MajorDim();
	int nex = fConstraintStatus.MinorDim(); /* number of constraints per element */
	int ndof = NumDOF();
	ArrayT<char> status;
	iArrayT constraint_eqnos;
	for (int i = 0; i < nel; i++)
	{
		fConstraintStatus.RowAlias(i, status);
		const iArrayT& disp_eq = ElementCard(i).Equations();
		int dof = 0;
		for (int j = 0; j < nex; j++)
		{
			/* collect equations */
			if (status[j] == kActive)
			{
				fXDOFEqnos.RowAlias(constraint_dex, constraint_eqnos);

				/* displacement equations */
				int n_disp_eq = constraint_eqnos.Length() - 1;
				int disp_eq_dex = dof; 
				for (int k = 0; k < n_disp_eq; k++) {
					constraint_eqnos[k] = disp_eq[disp_eq_dex];
					disp_eq_dex += ndof;
				}

				/* last equation is the constraint */
				constraint_eqnos.Last() = xdof_eqnos[constraint_dex];
			
				/* next */
				constraint_dex++;
			}
			
			/* cycle through nodal dof */
			if (++dof == ndof)
				dof = 0;
		}
	}

	/* add to list */
	eq_1.Append(&fXDOFEqnos);	
}

/* close current time increment */
void RigidCSEAnisoT::CloseStep(void)
{
	/* inherited */
	CSEAnisoT::CloseStep();

	/* update status history */
	fConstraintStatus_last = fConstraintStatus;

	/* update constraint history */
	XDOF_ManagerT& xdof_manager = ElementSupport().XDOF_Manager();
	const dArray2DT& constraints = xdof_manager.XDOF(this, 0);
	fConstraints_last.Dimension(constraints.Length());
	fConstraints_last = constraints;
}

/* restore the state */
GlobalT::RelaxCodeT RigidCSEAnisoT::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = CSEAnisoT::ResetStep();

	/* reset the status */
	fConstraintStatus = fConstraintStatus_last;

	/* Note: the XDOF_ManagerT handles calls to reset the values of
	 *       the constraints. */

	return relax;
}

/* determine the number of constraints needed */
void RigidCSEAnisoT::SetDOFTags(void)
{
	/* count number of active constraints */
	int count = 0;
	char* pa = fConstraintStatus.Pointer();
	int  len = fConstraintStatus.Length();
	for (int i = 0; i < len; i++)
		if (*pa++ == kActive)
			count++;

	/* resize tags array */
	fConstraintXDOFTags_man.SetLength(count, false);
}

/* returns 1 if group needs to reconfigure DOF's, else 0 */
int RigidCSEAnisoT::Reconfigure(void)
{
#if 0
	//TEMP - free one of the constraints
	double time = ElementSupport().Time();
	if (time > 0.29 && time < 0.31) {
		fConstraintStatus(0,1) = kFree;
		return 1;
	}
	else if (time > 0.69 && time < 0.71) {
		fConstraintStatus(0,3) = kFree;
		return 1;
	}
	else
#endif
		return 0;
}

/* collecting element connectivities for the field */
void RigidCSEAnisoT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1, 
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
	/* inherited */
	CSEAnisoT::ConnectsU(connects_1, connects_2);

	/* add connectivities with constraints */
	connects_1.AppendUnique(&fXDOFConnectivities);
}

iArrayT& RigidCSEAnisoT::DOFTags(int tag_set)
{
#if __option(extended_errorcheck)
	/* check */
	if (tag_set != 0)
		ExceptionT::OutOfRange("AugLagContact2DT::DOFTags", "expecting tag set 0: %d", tag_set);
#else
#pragma unused(tag_set)
#endif
	return fConstraintXDOFTags;
}

/* return the contact elements */
const iArray2DT& RigidCSEAnisoT::DOFConnects(int tag_set) const
{
#if __option(extended_errorcheck)
	/* check */
	if (tag_set != 0)
		ExceptionT::OutOfRange("AugLagContact2DT::DOFConnects", "expecting tag set 0: %d", tag_set);
#else
#pragma unused(tag_set)
#endif
	return fXDOFConnectivities;
}

/* restore the DOF values to the last converged solution */
void RigidCSEAnisoT::ResetDOF(dArray2DT& DOF, int tag_set) const
{
#if __option(extended_errorcheck)
	/* check */
	if (tag_set != 0)
		ExceptionT::OutOfRange("AugLagContact2DT::ResetDOF", "expecting tag set 0: %d", tag_set);
#else
#pragma unused(tag_set)
#endif

	int curr_dof_count = 0;
	int last_dof_count = 0;
	int num_status = fConstraintStatus.Length();
	for (int i = 0; i < num_status; i++)
	{
		char curr_status = fConstraintStatus[i];
		char last_status = fConstraintStatus_last[i];
	
		if (curr_status == kActive) {
			if (last_status == kActive)
				DOF[curr_dof_count] = fConstraints_last[last_dof_count];
			else
				DOF[curr_dof_count] = 0.0;
			curr_dof_count++;
		}
		
		if (last_status == kActive)
			last_dof_count++;
	}
}

/* generate nodal connectivities */
void RigidCSEAnisoT::GenerateElementData(void)
{
	/* resize the work space */
	int num_constraints = fConstraintXDOFTags.Length();
	fXDOFConnectivities_man.SetMajorDimension(num_constraints, false);

	/* collect nodes and DOF tags */
	int num_element_constraints = fConstraintStatus.MinorDim();
	int tag = 0;
	int nel = NumElements();
	iArrayT row;
	for (int i = 0; i < nel; i++)
	{
		ElementCardT& element_card = ElementCard(i);
		const iArrayT& nodes_U = element_card.NodesU();

		char* constraint_status = fConstraintStatus(i);
		for (int j = 0; j < num_element_constraints; j++)
			if (*constraint_status++ == kActive)
			{
				/* collect connectivities */
				fXDOFConnectivities.RowAlias(tag, row);
				row.CopyIn(0, nodes_U);
				row.Last() = fConstraintXDOFTags[tag];
		
				/* next */
				tag++;
			}
	}
}

/* describe the parameters needed by the interface */
void RigidCSEAnisoT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	CSEAnisoT::DefineParameters(list);

	/* regularization */
	ParameterT regularization(fr, "regularization");
	regularization.SetDefault(fr);
	regularization.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(regularization);
}

/* accept parameter list */
void RigidCSEAnisoT::TakeParameterList(const ParameterListT& list)
{
	/* regularization */
	fr = list.GetParameter("regularization");

	/* inherited */
	CSEAnisoT::TakeParameterList(list);

	/* construct surface shape functions (1 integration point) */
//	fConstraintShapes = new SurfaceShapeT(fGeometryCode, 1, NumElementNodes(), NumDOF(), fLocInitCoords1);
//	if (!fConstraintShapes) throw ExceptionT::kOutOfMemory;
//	fConstraintShapes->Initialize();
	fConstraintShapes = fShapes;

	/* set local arrays */
	int nen = NumElementNodes();
	int ndof = NumDOF();
	fDisp.Dimension(nen, ndof);
	fLastDisp.Dimension(nen, ndof);
	Field().RegisterLocal(fDisp);
	Field().RegisterLocal(fLastDisp);

	/* flags array */
	fConstraintStatus.Dimension(NumElements(), NumDOF()*fConstraintShapes->NumIP());
	fConstraintStatus = kActive;

	/* reset base class parameters */
	int neq = NumElementNodes() + 1; // 1 additional dof

	/* dynamic work space managers */
	fConstraintXDOFTags_man.SetWard(0, fConstraintXDOFTags);
	fXDOFConnectivities_man.SetWard(0, fXDOFConnectivities, NumElementNodes() + 1);		
	fXDOFEqnos_man.SetWard(0, fXDOFEqnos, NumElementNodes() + 1);

	/* register with node manager - sets initial fContactDOFtags */
	iArrayT xdof_tags(1);
	xdof_tags = 1;
	ElementSupport().XDOF_Manager().XDOF_Register(this, xdof_tags);
}

/* force vector */
void RigidCSEAnisoT::RHSDriver(void)
{
	const char caller[] = "RigidCSEAnisoT::RHSDriver";

	/* time-integration parameters */
	double constKd = 1.0;
	int formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* source for XDOF information */
	XDOF_ManagerT& xdof_manager = ElementSupport().XDOF_Manager();

	/* current values of the constraints */
	const dArray2DT& constraints = xdof_manager.XDOF(this, 0);
	if (constraints.Length() != fConstraintXDOFTags.Length())
		ExceptionT::GeneralFail(caller, "expecting %d constraints not %d",
			fConstraintXDOFTags.Length(), constraints.Length());

	/* jumps */
	int ndof = NumDOF();
	dArrayT jump(ndof), last_jump(ndof), jump_Na;

	dArrayT rhs(fXDOFEqnos.MinorDim());
	dArrayT rhs_u(fXDOFEqnos.MinorDim() - 1, rhs.Pointer());
	int constraint_dex = 0;
	ArrayT<char> status;
	iArrayT constraint_eqnos;
	Top();
	while (NextElement())
	{
		/* collect displacements */
		SetLocalU(fDisp);
		SetLocalU(fLastDisp);
		
		/* constraint status over the element */
		fConstraintStatus.RowAlias(CurrElementNumber(), status);

		/* loop over integration points */
		int element_constraint = 0;
		fConstraintShapes->TopIP();
		while (fConstraintShapes->NextIP())
		{
			/* compute displacement jumps */
			fConstraintShapes->InterpolateJump(fDisp, jump);
			fConstraintShapes->InterpolateJump(fLastDisp, last_jump);

			/* loop over dof */
			for (int i = 0; i < ndof; i++)
				if (status[element_constraint++] == kActive)
				{
					/* constraint and multiplier */
					double h = jump[i] - last_jump[i];
					double f = constraints[constraint_dex];

					/* dh/du */
					fConstraintShapes->JumpShapes(jump_Na);
				
					/* residual of the displacement equation */
					rhs_u.SetToScaled(-(f + fr*h), jump_Na);
					
					/* residual of the constraint equation */
					rhs.Last() = -h;

					/* get equation numbers */
					fXDOFEqnos.RowAlias(constraint_dex, constraint_eqnos);
	
					/* assemble */
					ElementSupport().AssembleRHS(Group(), rhs, constraint_eqnos);
			
					/* next constraint */
					constraint_dex++;
				}
		}
	}
}

/* tangent matrix */
void RigidCSEAnisoT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)

	/* time-integration parameters */
	double constK = 1.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	/* jumps */
	int ndof = NumDOF();

	/* workspace */
	dArrayT vec_nee(fXDOFEqnos.MinorDim()), jump_Na;
	vec_nee = 0.0;
	ElementMatrixT lhs(fXDOFEqnos.MinorDim(), ElementMatrixT::kSymmetric);

	int constraint_dex = 0;
	ArrayT<char> status;
	iArrayT constraint_eqnos;
	Top();
	while (NextElement())
	{
		/* constraint status over the element */
		fConstraintStatus.RowAlias(CurrElementNumber(), status);

		/* loop over integration points */
		int element_constraint = 0;
		fConstraintShapes->TopIP();
		while (fConstraintShapes->NextIP())
		{
			/* loop over dof */
			for (int i = 0; i < ndof; i++)
				if (status[element_constraint++] == kActive)
				{
					/* dh/du -> into vec_nee */
					fConstraintShapes->JumpShapes(jump_Na);
					vec_nee.CopyIn(0, jump_Na);
				
					/* displacement part */
					lhs.Outer(vec_nee, vec_nee, fr);

					/* cross derivative */
					int dex = fXDOFEqnos.MinorDim() - 1;
					lhs.SetRow(dex, vec_nee);
					lhs.SetCol(dex, vec_nee);

					/* get equation numbers */
					fXDOFEqnos.RowAlias(constraint_dex, constraint_eqnos);
	
					/* assemble */
					ElementSupport().AssembleLHS(Group(), lhs, constraint_eqnos);
			
					/* next constraint */
					constraint_dex++;
				}
		}
	}
}
