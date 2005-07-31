/* $Id: ExpCD_DRSolver.cpp,v 1.1.1.1 2001-01-29 08:20:34 paklein Exp $ */
/* created: paklein (08/19/1998)                                          */

#include "ExpCD_DRSolver.h"

#include <iostream.h>
#include <math.h>

#include "fstreamT.h"
#include "Constants.h"
#include "ExceptionCodes.h"
#include "FEManagerT.h"
#include "DiagonalMatrixT.h"

/* iteration status flags */
const int kContinue  = 0;
const int kConverged = 1;
const int kFailed    = 2;

/* constructor */
ExpCD_DRSolver::ExpCD_DRSolver(FEManagerT& fe_manager):
	SolverT(fe_manager),
	fMaxIterations(-1),
	fTolerance(0.0),
	fdt(1.0)
{
	/* check matrix type */
	if (fMatrixType != kDiagonalMatrix)
	{
		cout << "\n ExpCD_DRSolver::ExpCD_DRSolver: expecting matrix type: "
		     << kDiagonalMatrix << endl;
		throw eBadInputValue;
	}
	else
	{
#ifdef __NO_RTTI__
		DiagonalMatrixT* diagonal_matrix = (DiagonalMatrixT*) fLHS;
#else	
		DiagonalMatrixT* diagonal_matrix = dynamic_cast<DiagonalMatrixT*>(fLHS);
		if (!fLHS)
		{
			cout << "\n ExpCD_DRSolver::ExpCD_DRSolver: cast of global matrix to\n"
			     <<   "    diagonal matrix failed" << endl;
			throw eGeneralFail;
		}
#endif
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
		fOutputDOF = fFEManager.GlobalEquationNumber(nodenum, dofnum);

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
	if (fMaxIterations < 0) throw eBadInputValue;
	if (fTolerance < 0.0 || fTolerance > 1.0) throw eBadInputValue;
	if (fMass_scaling < kSmall) throw eBadInputValue;
	if (fDamp_scaling < kSmall) throw eBadInputValue;
}

/* (re-)configure the global equation system */
void ExpCD_DRSolver::Initialize(int tot_num_eq, int loc_num_eq, int start_eq)
{
	/* inherited */
	SolverT::Initialize(tot_num_eq, loc_num_eq, start_eq);

	/* dimension */
	fDis.Allocate(loc_num_eq);
	fVel.Allocate(loc_num_eq);
	fAcc.Allocate(loc_num_eq);
}

/* generate the solution for the current time sequence */
void ExpCD_DRSolver::Run(void)
{
	/* generate pseudo-mass matrix */
	SetMass();

	/* solve displacements for quasi-static load sequence */
	while (Step())
	{			
		/* residual loop */
		try {
	
			/* apply kinematic BC's */
			fFEManager.InitStep();
		
			/* form the residual force vector */
			fRHS = 0.0;
			fFEManager.FormRHS();	
			double error = fRHS.Magnitude();
			
			/* loop on error */
			int solutionflag = ExitIteration(error);
			while(solutionflag == kContinue)
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
				GlobalT::RelaxCodeT relaxcode = fFEManager.RelaxSystem();
				
				/* reset global equations */
				if (relaxcode == GlobalT::kReEQ ||
				    relaxcode == GlobalT::kReEQRelax)
					fFEManager.Reinitialize();
					
				/* new equilibrium */					
				if (relaxcode == GlobalT::kRelax ||
				    relaxcode == GlobalT::kReEQRelax)
					Relax();
			}
			
			/* finalize */
			fFEManager.CloseStep();
		}

		catch (int code) { fFEManager.HandleException(code); }		
	}
}

/*************************************************************************
* Protected
*************************************************************************/

/* advance to next load step. Returns 0 if there are no more
* steps. Overload to add class dependent initializations */
int ExpCD_DRSolver::Step(void)
{
	/* reset iteration count */
	fNumIteration = -1;

	/* "initial conditions" for DR */
	fDis = 0.0;
	fVel = 0.0;
	fAcc = 0.0;
		
	/* inherited */
	return SolverT::Step();
}

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
int ExpCD_DRSolver::ExitIteration(double error)
{
	/* CONVERGENCE MOVIES */
//TEMP
//	if (fKinePrintInc > 0 && fmod(fNumIteration+1,fKinePrintInc) == 0)
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

		fFEManager.ResetStep();
		fFEManager.DecreaseLoadStep();
			
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

int ExpCD_DRSolver::ExitRelaxation(double error)
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

		fFEManager.ResetStep();
		fFEManager.DecreaseLoadStep();

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
	fLHS->Solve(fRHS);

	/* updates */
	fDis.SetToCombination(fdt, fVel, 0.5*fdt*fdt, fAcc);
	fVel.AddCombination(0.5*fdt, fAcc, 0.5*fdt, fRHS);
	fAcc = fRHS;

	/* incremental update */
	fFEManager.Update(fDis);

	/* compute new residual */
	fRHS = 0.0;
	fFEManager.FormRHS();

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
		fFEManager.FormRHS();
		double error = fRHS.Magnitude();	
		
		while ( !ExitRelaxation(error) )
		{
			if (++count == newtancount) count = 0;
			
			error = SolveAndForm();
		}	
	}
			
	catch (int ErrorCode)
	{
		cout << "\nExpCD_DRSolver::Run: encountered exception on relaxation:\n";
		throw eGeneralFail;
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
	fFEManager.FormLHS();
	
	/* get the matrix - should be safe */
	DiagonalMatrixT* lhs = (DiagonalMatrixT*) fLHS;
	dArrayT&  massmatrix = lhs->TheMatrix();

	/* should all be positive */
	if (massmatrix.Min() < kSmall) throw eGeneralFail;
	
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
