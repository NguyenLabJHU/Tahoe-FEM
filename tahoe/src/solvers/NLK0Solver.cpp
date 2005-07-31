/* $Id: NLK0Solver.cpp,v 1.1.1.1 2001-01-29 08:20:34 paklein Exp $ */
/* created: paklein (10/01/1996)                                          */
/* Solver that forms and factorizes the tangent matrix only once          */
/* (for now)                                                              */

#include "NLK0Solver.h"
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include "Constants.h"
#include "ExceptionCodes.h"
#include "FEManagerT.h"
#include "fstreamT.h"

/* line search parameters */
const double ks_max_factor = 5.0; //normally 1.0

/* constructor */
NLK0Solver::NLK0Solver(FEManagerT& fe_manager):
	NLSolver(fe_manager),
	fFormTangent(1),
	fLastTangent(fe_manager.Output(), fLHS->CheckCode())
{
#ifdef __NO_RTTI__
	pCCSLHS = NULL; //not allowed for now
#else
	pCCSLHS = dynamic_cast<CCSMatrixT*>(fLHS);
#endif
	if (!pCCSLHS) throw eGeneralFail;
}

/* generate the solution for the current time sequence */
void NLK0Solver::RunK0(void) //not used
{
	/* form tangent for inital state */
	fLHS->Clear();
	fFEManager.FormLHS();
	
	/* solve displacements for quasi-static load sequence */
	while ( Step() )
	{					
		/* residual loop */
		try {

			/* apply kinematic BC's - set displacements and update geometry */
			fFEManager.InitStep();
	
			/* form the residual force vector */
			fRHS = 0.0;
			fFEManager.FormRHS();
			double error = fRHS.Magnitude();	
			
			/* loop on error */
			while( !ExitIteration(error) )
			{
				/* solve equation system */
				fLHS->Solve(fRHS);
			
				/* update system */
				fFEManager.Update(fRHS);
				//LineSearchUpdate(1.0);
								
				/* compute new residual */
				fRHS = 0.0;
				fFEManager.FormRHS();
				
				/* compute new error */
				error = fRHS.Magnitude();
			}

			/* relaxation if elements deleted */
			if ( 0 && fFEManager.RelaxSystem() )
				Relax(); //disabled since this reforms the tangent
				
			/* finalize */
			fFEManager.CloseStep();
		}

		catch (int code) { fFEManager.HandleException(code); }
	}
}

/***********************************************************************
* Protected
***********************************************************************/

/* form and solve - returns the magnitude of the residual */
double NLK0Solver::SolveAndForm(bool junk)
{		
#pragma unused(junk)

	/* form the stiffness matrix */
	if (fFormTangent)
	{
		pCCSLHS->Clear();
		fFEManager.FormLHS();

		/* solve equation system */
		fUpdate = fRHS;
		pCCSLHS->Solve(fRHS);
	
		/* check for positive definiteness */
		if ( pCCSLHS->HasNegativePivot() )
		{
			fFormTangent = 0; //no switching back, yet
		
			/* restore RHS */
			fRHS = fUpdate;
		}
		else /* store tangent */
			fLastTangent = (*pCCSLHS);
	}
	
	/* use last positive definite tangent */
	if (!fFormTangent)
		fLastTangent.Solve(fRHS);
			 		
	/* update system */
	fFEManager.Update(fRHS);
								
	/* compute new residual */
	fRHS = 0.0;
	fFEManager.FormRHS();

	/* combine residual magnitude with update magnitude */
	/* e = a1 |R| + a2 |delta_d|                        */
	//not implemented!
			
	return( fRHS.Magnitude() );
}

/***********************************************************************
* Private
***********************************************************************/

/* perform linesearching to find the best step length */
void NLK0Solver::InitLineSearch(void)
{
	fUpdate = fRHS;
	s_current = 0.0;
	fRecursionDepth = 0;
}

void NLK0Solver::TermLineSearch(void)
{
	if (s_current != 0.0)
	{
		/* set fRHS = delta_d */
		fRHS = fUpdate;
		
		/* restore to ith state */
		fUpdate *= (-s_current);
		fFEManager.Update(fUpdate);
		s_current = 0.0; //insurance	
	}
}
void NLK0Solver::LineSearchUpdate(double s_upper)
{
	int d_width = fFEManager.Output().precision() + kDoubleExtra;

	/* initialize */
	InitLineSearch();

	//TEMP - plot G(s) ( 0.0 < s < s_upper )
	int PrintStuff = 0;
	if (PrintStuff)
	{
		fFEManager.Output() << "\n************\n";
		
		/* print 50 points within the range */
		for (int i = 0; i < 51; i++)
		{
			double value = s_upper*i/50.0;
			fFEManager.Output() << setw(d_width) << value;
			fFEManager.Output() << setw(d_width) << GValue(value);
			fFEManager.Output() << setw(d_width) << fRHS.Magnitude() << '\n';
		}
		
		fFEManager.Output() << "\n************\n";
		
		int QuitNow = 0;
		if (QuitNow) throw eGeneralFail;
	}
		
	/* initial bounds */
	double s_a = 0.0;
	double s_b = s_upper;

	double G_b = GValue(s_b);
	
	/* quick exit - full step */
	if (fabs(G_b)*100 < kSmall)
	{
		/* print step size */
		cout << setw(d_width) << s_current;
		
		return;
	}
	
	double G_a = GValue(s_a);	
	double G_0 = ( fabs(G_a) > fabs(G_b) ) ? G_b : G_a;

	/* minimize G with secant method */		
	int count = 0;
	double G_new;
	int GiveUp = 0;
	do {
		double m = (G_a - G_b)/(s_a - s_b);
		double b = G_b - m*s_b;
		double s_new = -b/m;
		
		/* limits Newton step - jump out */
		if (s_new > s_upper*ks_max_factor)
		{
			double z = s_upper*ks_max_factor;
			GValue(z);
			GiveUp = 1;
			break;
		}
		else if (s_new < -s_upper*ks_max_factor)		
		{
			double z =-s_upper*ks_max_factor;
			GValue(z);
			GiveUp = 1;
			break;
		}
				
		G_new = GValue(s_new);
		if ( fabs(G_a) > fabs(G_new) && fabs(G_a) > fabs(G_b) )
		{
			G_a = G_new;
			s_a = s_new;
			GiveUp = 0;
		}
		else if ( fabs(G_b) > fabs(G_new) && fabs(G_b) > fabs(G_a) )
		{
			G_b = G_new;
			s_b = s_new;
			GiveUp = 0;
		}
		else
		{
			/* G_a and G_b don't bracket zero */
			if (G_b*G_a > 0)
			{
				if (G_a*G_new < 0)
				{
					G_a = G_new;
					s_a = s_new;
				}
				else if (G_b*G_new < 0)
				{
					G_b = G_new;
					s_b = s_new;
				}
				else /* exit */
					GiveUp = 1;
			}
			else /* exit */
				GiveUp = 1;
		
		}
		
	} while (fabs(G_new)*100 > kSmall &&
	         //fabs(G_0/G_new) < 10000 &&
	         fabs(G_0/G_new) < 1.0e10 && // "exact" line search
	         ++count < 15 &&
	         !GiveUp);
	
	/* recursive reset for jump out */
	if (GiveUp == 2)
		LineSearchUpdate(s_upper/2.0);		
	
	/* print step size */
	cout << setw(d_width) << s_current;
	
	/* assume line search failed */
//	if (count == 10)
//		GValue(1.0);
	
	/* try parabolic fit if line search fails */
//	if (fRHS.Magnitude() > upper)
//	{
//		TermLineSearch();
//	
//		ParabolicUpdate(upper);
//	}	
	
//	cout << "      NLK0Solver::LineSearchUpdate: step: " << s_current << ": repetitions: " << count;	
} 	

/* return the line search weight function for the given step size.
* The degrees of freedom:
*
*	G = R(d_i+1).delta_d */
double NLK0Solver::GValue(double& step)
{
	/* scale update vector */
	fRHS = fUpdate;
	
	/* scale accounting for current step update */
	fRHS *= (step - s_current);
	s_current = step;
	
	/* compute residual */
	fFEManager.Update(fRHS);
	fRHS = 0.0;
	try { fFEManager.FormRHS(); }

	/* Jacobian error during form RHS */
/*
	catch (int ErrorCode)
	{
		if (ErrorCode == eBadJacobianDet)
			cout << "\nNLK0Solver::GValue: Jacobian exception\n" << endl;
		else
			cout << "\nNLK0Solver::GValue: unknown exception\n" << endl;
		throw eGeneralFail;
	}
*/

	catch (int ErrorCode)
	{
		if (ErrorCode == eBadJacobianDet)
		{
			cout << "\nNLK0Solver::GValue: Jacobian exception\n" << endl;
			throw eGeneralFail;
		
			if (++fRecursionDepth > 8)
			{
				cout << "\nNLK0Solver::GValue: max recursion depth reached\n" << endl;
				throw ;
			}
			else
			{
				step *= -0.5;
				return(GValue(step));
			}
		}
		else
		{
			cout << "\nNLK0Solver::GValue: unknown exception\n" << endl;
			throw eGeneralFail;
		}
	}

	return InnerProduct(fUpdate,fRHS);
}
