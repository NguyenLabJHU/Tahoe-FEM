/* $Id: ExpCD_DRSolver.cpp,v 1.12 2004-10-14 20:25:07 paklein Exp $ */
/* created: paklein (08/19/1998) */
#include "ExpCD_DRSolver.h"

#include <iostream.h>
#include <math.h>

#include "ifstreamT.h"
#include "FEManagerT.h"
#include "DiagonalMatrixT.h"

using namespace Tahoe;

/* constructor */
ExpCD_DRSolver::ExpCD_DRSolver(FEManagerT& fe_manager, int group):
	SolverT(fe_manager, group),
	fMaxIterations(-1),
	fTolerance(0.0),
	fdt(1.0)
{
ExceptionT::GeneralFail("ExpCD_DRSolver::ExpCD_DRSolver", "out of date");
#if 0
	/* check matrix type */
	if (fMatrixType != kDiagonalMatrix)
	{
		cout << "\n ExpCD_DRSolver::ExpCD_DRSolver: expecting matrix type: "
		     << kDiagonalMatrix << endl;
		throw ExceptionT::kBadInputValue;
	}
	else
	{
		DiagonalMatrixT* diagonal_matrix = TB_DYNAMIC_CAST(DiagonalMatrixT*, fLHS);
		if (!fLHS)
		{
			cout << "\n ExpCD_DRSolver::ExpCD_DRSolver: cast of global matrix to\n"
			     <<   "    diagonal matrix failed" << endl;
			throw ExceptionT::kGeneralFail;
		}

		/* reset assembly mode */
		diagonal_matrix->SetAssemblyMode(DiagonalMatrixT::kAbsRowSum);
	}
	
	ifstreamT& in = fFEManager.Input();
	
	/* read parameters */
	in >> fMaxIterations;
	in >> fTolerance;
	in >> fMass_scaling;
	in >> fDamp_scaling;

	/* convergence history output stream */
	int nodenum, dofnum;
	in >> nodenum >> dofnum;
	
	/* correct offset */
	nodenum--; dofnum--;
	
	if (nodenum > -1)
	{
//NOTE - not updated for generalized multi-field. Would need to know
//       which field to query
//		fOutputDOF = fFEManager.GlobalEquationNumber(nodenum, dofnum);
cout << "\n ExpCD_DRSolver::ExpCD_DRSolver: not update for multi-field" << endl;
throw ExceptionT::kGeneralFail;

		StringT outname("dof");
		outname.Append(".", fOutputDOF);
		outname.Append(".out");
		
		/* open output stream */
		fhist_out.open(outname);
	}
	else
		fOutputDOF = 0;

	/* print parameters */
	ostream& out = fFEManager.Output();	
	out << "\n O p t i m i z a t i o n   P a r a m e t e r s :\n\n";
	out << " Maximum number of iterations. . . . . . . . . . = " << fMaxIterations << '\n';
	out << " Convergence tolerance . . . . . . . . . . . . . = " << fTolerance     << '\n';
	out << " Mass scaling. . . . . . . . . . . . . . . . . . = " << fMass_scaling  << '\n';
	out << " Damping scaling . . . . . . . . . . . . . . . . = " << fDamp_scaling  << '\n';
	out << " Output DOF. . . . . . . . . . . . . . . . . . . = ";
	if (fOutputDOF > 0)
		out << fOutputDOF;
	else
		out << "<none>";
	out << '\n';

	/* checks */
	if (fMaxIterations < 0) throw ExceptionT::kBadInputValue;
	if (fTolerance < 0.0 || fTolerance > 1.0) throw ExceptionT::kBadInputValue;
	if (fMass_scaling < kSmall) throw ExceptionT::kBadInputValue;
	if (fDamp_scaling < kSmall) throw ExceptionT::kBadInputValue;
#endif
}

/* (re-)configure the global equation system */
void ExpCD_DRSolver::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* inherited */
	SolverT::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* dimension */
	fDis.Dimension(loc_num_eq);
	fVel.Dimension(loc_num_eq);
	fAcc.Dimension(loc_num_eq);
}

/* solve the system over the current time increment */
SolverT::SolutionStatusT ExpCD_DRSolver::Solve(int num_iterations)
{
	try {

	/* reset iteration count */
	fNumIteration = -1;

	/* "initial conditions" for DR */
	fDis = 0.0;
	fVel = 0.0;
	fAcc = 0.0;
	
	/* generate pseudo-mass matrix */
	if (fFEManager.StepNumber() == 0) SetMass();
	
	/* form the residual force vector */
	fRHS = 0.0;
	fFEManager.FormRHS(Group());	
	double error = fRHS.Magnitude();
			
	/* loop on error */
	SolutionStatusT solutionflag = ExitIteration(error);
	while (solutionflag == kContinue && 
		(num_iterations == -1 || fNumIteration < num_iterations))
	{
		/* explicit central difference time stepping */
		error = SolveAndForm();
		solutionflag = ExitIteration(error);

		/* convergence history */
		if (fOutputDOF > 0)
		{
			fhist_out << fVel[fOutputDOF - 1] << '\t';
			fhist_out << error/fError0      << '\n';
		}
	}

	/* found solution */
	if (solutionflag == kConverged)
	{
		/* relaxation */
		GlobalT::RelaxCodeT relaxcode = fFEManager.RelaxSystem(Group());
				
		/* reset global equations */
		if (relaxcode == GlobalT::kReEQ ||
			relaxcode == GlobalT::kReEQRelax)
			fFEManager.SetEquationSystem(Group());
					
		/* new equilibrium */					
		if (relaxcode == GlobalT::kRelax ||
			relaxcode == GlobalT::kReEQRelax)
			Relax();
	} 
			
	return solutionflag;
	} /* end try */

	/* abnormal */
	catch (ExceptionT::CodeT code) { 
		return kFailed;
	}
}

/*************************************************************************
* Protected
*************************************************************************/

/* returns 1 if the iteration loop should be left, otherwise
* returns 0.  The iteration loop can be exited for the
* following reasons:
*
*	(1) error meets the convergence criteria
*	(2) the iteration limit has been exceeded
*	(3) the error is diverging
*
* For (2) and (3), the load increment will be cut and the
* iteration re-entered with the next Step() call */
SolverT::SolutionStatusT ExpCD_DRSolver::ExitIteration(double error)
{
	/* CONVERGENCE MOVIES */
//TEMP
//	if (fKinePrintInc > 0 && (fNumIteration+1)%fKinePrintInc == 0)
//		fFEManager.PrintKinematic(fOutput,fTime);

	++fNumIteration;
	
	/* first pass */
	if (fNumIteration == 0)
	{
		cout <<   " Absolute error = " << error << '\n';
		cout <<   " Relative error :\n\n";
		cout << setw(kIntWidth)    << fNumIteration;

		fError0 = error;
		
		/* hit on first try */
		if (fError0 < kSmall)
		{
			cout << setw(kDoubleWidth) << error << endl;
			return kConverged;
		}
		else
		{
			cout << setw(kDoubleWidth) << 1.0 << endl;
			return kContinue;
		}
	}
	/* iteration limit hit */
	else if (fNumIteration >= fMaxIterations)
	{
		cout << "\n ExpCD_DRSolver::ExitIteration: max iterations hit\n" << endl;
		return kFailed;
	}
	/* interpret error */
	else
	{
		double relerror = error/fError0;
		cout << setw(kIntWidth)    << fNumIteration;
		cout << setw(kDoubleWidth) << relerror << endl;

		/* converged */
		if (relerror < fTolerance || error < 10.0*kSmall)
		{
			cout << '\n';
			
			fFEManager.Output() << " Converged at time = " << fFEManager.Time() << '\n';
					
			return kConverged;
		}
		/* continue iterations */
		else
			return kContinue;
	}
}

SolverT::SolutionStatusT ExpCD_DRSolver::ExitRelaxation(double error)
{
	++fNumIteration;
	
	/* first pass */
	if (fNumIteration == 0)
	{
		cout <<   " Absolute error = " << error << '\n';
		cout <<   " Relative error :\n\n";

		fError0 = error;

		/* hit on first try */
		if (fError0 < kSmall)
		{
			cout << '\t' << error << endl;
			return kConverged;
		}
		else
		{
			cout << '\t' << 1.0 << endl;
			return kContinue;
		}
	}
	/* iteration limit hit */
	else if (fNumIteration >= fMaxIterations)
	{
		cout << "\n ExpCD_DRSolver::ExitRelaxation: max iterations hit\n" << endl;
		return kFailed;
	}
	/* interpret error */
	else
	{
		double relerror = error/fError0;
		cout << '\t' << relerror << endl;

		/* converged */
		if (relerror < fTolerance || error < 10.0*kSmall)
		{
			cout << '\n';
			fFEManager.Output() << " Converged at time = " << fFEManager.Time() << '\n';
					
			return kConverged;
		}
		/* continue iterations */
		else
			return kContinue;
	}
}

/* form and solve the equation system */
double ExpCD_DRSolver::SolveAndForm(void)
{
	/* set damping parameter */
	double c = SetDamping();

	/* build dynamic residual */
	double* pres = fRHS.Pointer();
	double* pvel = fVel.Pointer();
	double* pmas = pMass->Pointer();
	double* pacc = fAcc.Pointer();	
	int   numeqs = fRHS.Length();
	for (int i = 0; i < numeqs; i++)
	{
		*pres += (-(*pacc)/(*pmas)) - c*(*pvel);	
		pres++; pvel++; pmas++; pacc++;
	}
		
	/* get a_n+1 */
	if (!fLHS->Solve(fRHS)) throw ExceptionT::kBadJacobianDet;

	/* updates */
	fDis.SetToCombination(fdt, fVel, 0.5*fdt*fdt, fAcc);
	fVel.AddCombination(0.5*fdt, fAcc, 0.5*fdt, fRHS);
	fAcc = fRHS;

	/* incremental update */
	fFEManager.Update(Group(), fDis);

	/* compute new residual */
	fNumIteration++;
	fRHS = 0.0;
	fFEManager.FormRHS(Group());

	return fRHS.Magnitude();
}

/* relax system */
void ExpCD_DRSolver::Relax(int newtancount)
{	
	cout << " Relaxation:" << '\n';

	/* reset iteration count */
	fNumIteration = -1;
		
	/* residual loop */
	try {
	
		int count = newtancount - 1;
	
		/* form the residual force vector */
		fRHS = 0.0;
		fFEManager.FormRHS(Group());
		double error = fRHS.Magnitude();	
		
		while ( !ExitRelaxation(error) )
		{
			if (++count == newtancount) count = 0;
			
			error = SolveAndForm();
		}	
	}
			
	catch (ExceptionT::CodeT ErrorCode)
	{
		cout << "\nExpCD_DRSolver::Run: encountered exception on relaxation:\n";
		throw ExceptionT::kGeneralFail;
	}
}

/*************************************************************************
* Private
*************************************************************************/

/* set pseudo-mass */
void ExpCD_DRSolver::SetMass(void)
{
	/* get diagonal stiffness */
	fLHS->Clear();
	fFEManager.FormLHS(Group(), GlobalT::kDiagonal);
	
	/* get the matrix - should be safe */
	DiagonalMatrixT* lhs = (DiagonalMatrixT*) fLHS;
	dArrayT&  massmatrix = lhs->TheMatrix();

	/* should all be positive */
	if (massmatrix.Min() < kSmall) throw ExceptionT::kGeneralFail;
	
	/* store the diagonal stiffness matrix (for optimal damping calc) */
	fK_0 = massmatrix;

	/* stable mass matrix - P.Underwood (2.40) */
	massmatrix *= 0.25*fMass_scaling*1.1*1.1*fdt*fdt;

	/* keep a pointer to the mass array */
	pMass = &massmatrix;
}

double ExpCD_DRSolver::SetDamping(void)
{
	/* dominant mode estimation by Papadrakakis */
	if (fNumIteration < 2)
		return 0.0; // no damping on first iteration
	else
	{
		/* Underwood optimal damping (3.2) */
		double* pdisp = fDis.Pointer();
		double* pmass = pMass->Pointer();
		double* pK_0  = fK_0.Pointer();
		
		int length = fDis.Length();
		double num = 0.0;
		double den = 0.0;
		for (int i = 0; i < length; i++)
		{
			double d_sqr = (*pdisp)*(*pdisp);
		
			num += d_sqr*(*pK_0);
			den += d_sqr*(*pmass);
			
			pdisp++; pmass++; pK_0++;
		}
		
		return fDamp_scaling*2.0*sqrt(num/den);
	}
}
