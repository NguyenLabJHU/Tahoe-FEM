/* $Id: LinearSolver.cpp,v 1.1.1.1 2001-01-29 08:20:34 paklein Exp $ */
/* created: paklein (05/30/1996)                                          */

#include "LinearSolver.h"
#include "FEManagerT.h"

/* constructor */
LinearSolver::LinearSolver(FEManagerT& fe_manager):
	SolverT(fe_manager),
	fFormLHS(1)
{

}

/* signal new reconfigured system */
void LinearSolver::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* inherited */
	SolverT::Initialize(tot_num_eq, loc_num_eq, start_eq);
	
	/* flag to reform LHS */
	fFormLHS = 1;

#if 0
	int analysis_code = fFEManager.Analysis();
	if (analysis_code != GlobalT::kLinExpDynamic &&
		analysis_code != GlobalT::kNLExpDynamic)
		fFormLHS = 1;
#endif
//NOTE: these checks were added because explicit dynamics with
//      contact was reforming the mass matrix whenever the contact
//      configuration was changed. until this state is more clearly
//      defined {ReEQ, Relax, ReEQRelax, ????}
//
//  Need a flag to say that the force configured:
//    ReEQ_LHS, ReEQ_RHS
}

/* solve for the current time sequence */
void LinearSolver::Run(void)
{
	/* single-step update loop */
	while (Step())
	{
		try
		{
			/* initialize */
			fRHS = 0.0;
			
			/* apply kinematic BC's */
			fFEManager.InitStep();

			/* form the residual force vector */
			fFEManager.FormRHS();
					
			/* solve equation system */
			if (fFormLHS)
			{
				/* initialize */
				fLHS->Clear();
	
				/* form the stiffness matrix */
				fFEManager.FormLHS();
				
				/* flag not to reform */
				fFormLHS = 0;
			}

			/* determine update vector */
			fLHS->Solve(fRHS);

			/* update displacements */
			fFEManager.Update(fRHS);		
			
			/* relaxation */
			GlobalT::RelaxCodeT relaxcode = fFEManager.RelaxSystem();
				
			/* relax for configuration change */
			if (relaxcode == GlobalT::kRelax) fFormLHS = 1;
				//NOTE: NLSolver calls "fFEManager.Reinitialize()". Should this happen
				//      here, too? For statics, should also reset the structure of
				//      global stiffness matrix, but since EFG only breaks connections
				//      and doesn't make new ones, this should be OK for now. PAK (03/04/99)
			
			/* trigger set of new equations */
			if (relaxcode == GlobalT::kReEQ ||
			    relaxcode == GlobalT::kReEQRelax)
				fFEManager.Reinitialize();
				
			/* finalize */
			fFEManager.CloseStep();
		}

		catch (int code) { fFEManager.HandleException(code); }
	}
}	
