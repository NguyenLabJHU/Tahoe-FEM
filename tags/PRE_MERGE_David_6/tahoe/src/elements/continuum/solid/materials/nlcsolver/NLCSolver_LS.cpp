/* $Id: NLCSolver_LS.cpp,v 1.8 2003-01-27 07:00:28 paklein Exp $ */
#include "NLCSolver_LS.h"
#include "ExceptionT.h"

// some macros
#ifndef max
static int max(int i1, int i2) {return i1 >= i2 ? i1 : i2;};
static double max(double d1, double d2) {return d1 >= d2 ? d1 : d2;};
#endif

using namespace Tahoe;

NLCSolver_LS::NLCSolver_LS(const int dim, const NewtonMethodBase& method)
  : NLCSolver(dim, method), 
    fAlpha (1.e-4)
{ }

NLCSolver_LS::~NLCSolver_LS()
{ }

void NLCSolver_LS::ResetMemberData()
{
  // inherited
  NLCSolver::ResetMemberData();

  // reset member data
  fPrevF = 0.;
}

void NLCSolver_LS::SetAlpha(double alpha)
{
	if (alpha >= 0.5) ExceptionT::GeneralFail("NLCSolver_LS::SetAlpha", "alpha >= 0.5: %g", alpha);
	fAlpha = alpha;
}

void NLCSolver_LS::ComputeTrialPoint(dArrayT& X)
{
  double tempLambda = 0.e0;

  // get newton step, if needed
  NLCSolver::ComputeTrialPoint(X);

  // try full newton step
  if (fRejectionCount == 0)
    {
      fSlope = dArrayT::Dot(fGradF, fNewtonStep);
      double relLength = 0.;
      for (int i = 0; i < fDim; i++)
	relLength = max(relLength, fabs(fNewtonStep[i]/max(fabs(X[i]), fTypX[i])));
      fMinLambda = fStepTol / relLength;
      fLambda = 1.;
      if (NLCS_MESSAGES) {
       	cout << "NLCSolver_LS::ComputeTrialPoint: trying full Newton step" << endl;
	cout << "  fMinLambda = " << fMinLambda << endl;
      }
    }
  // first backtrack, use quadratic fit
  else if (fRejectionCount == 1)
    {
      tempLambda = - fSlope / (2.* (fF - fLastAcceptedF - fSlope));
      if (NLCS_MESSAGES) {
       	cout << "NLCSolver_LS::ComputeTrialPoint: doing quad backtrack" << endl;
	cout << "  lambda = " << tempLambda << endl;
      }
    }
  // second backtrack and on, use cubic fit
  else  
    {
      double l2  = 1. / fLambda / fLambda;
      double lp2 = 1. / fPrevLambda / fPrevLambda;
      double f1   = fF - fLastAcceptedF - fLambda * fSlope;
      double f2  = fPrevF - fLastAcceptedF - fPrevLambda * fSlope; 
      double d   = 1. / (fLambda - fPrevLambda);
      double c1  = (l2 * f1 - lp2 * f2) * d;
      double c2  = (-fPrevLambda * l2 * f1 + fLambda * lp2 * f2) * d;
      if (c1 == 0)
	tempLambda = - fSlope / (2. * c2);
      else
	tempLambda = (-c2 + sqrt(c2 * c2 - 3. * c1 * fSlope)) / (3. * c1);
      if (tempLambda > 0.5 * fLambda) tempLambda = 0.5 * fLambda;
      if (NLCS_MESSAGES) {
       	cout << "NLCSolver_LS::ComputeTrialPoint: doing cubic backtrack " << endl;
	cout << "  lambda = " << tempLambda << endl;
      }
    }

  // check if value of tempLambda is NaN or <= 1.e-150 
//  if ( is_NaN(tempLambda) || (fabs(tempLambda) <= 1.e-150 && fRejectionCount > 0) )
  //if ( is_NaN(tempLambda) )
	if (fabs(tempLambda) <= 1.e-150 && fRejectionCount > 0)
    {
      //writeMessage("NLCSolver_LS::ComputeTrialPoint: lambda = NaN");
      if (NLCS_MESSAGES) {
         cout << "* tempLambda = * = " << tempLambda << endl
              << "  fLambda        = " << fLambda   << endl
              << "  fPrevLambda    = " << fPrevLambda << endl
              << "  fF             = " << fF        << endl
              << "  fPrevF         = " << fPrevF        << endl
              << "  fLastAcceptedF = " << fLastAcceptedF << endl
              << "  fSlope         = " << fSlope    << endl
              << "  fTypF          = " << fTypF     << endl
              << "  fNewtonStep    = " << endl << fNewtonStep << endl 
              << "  fGradF         = " << endl << fGradF      << endl
              << "  X              = " << endl << X           << endl
              << "  fLastAcceptedX = " << endl << fLastAcceptedX << endl
              << "  fTypX          = " << endl << fTypX       << endl
              << "  fRHS           = " << endl << fRHS        << endl
              << "  fLHS           = " << endl << fLHS        << endl;
      }
		ExceptionT::GeneralFail("NLCSolver_LS::ComputeTrialPoint", "lambda is NaN or < 1.e-150");
    }

  // after first backtrack, be prepared for cubic backtrack
  if (fRejectionCount > 0)
    {
      fPrevLambda = fLambda;
      fPrevF = fF;
      fLambda = max(tempLambda, 0.1*fLambda);
    }
  
  // compute trial solution
  X.SetToCombination(1., fLastAcceptedX, fLambda, fNewtonStep);
}

void NLCSolver_LS::TestTrialPoint(dArrayT& trialX)
{
	const char caller[] = "NLCSolver_LS::TestTrialPoint";
  // check value of Lambda
  if (fLambda < fMinLambda) 
    {
      if (NLCS_MESSAGES)
        ExceptionT::GeneralFail(caller, "lambda < minlambda: %g < %g", fLambda, fMinLambda);
    }
  
  // check sufficient function decrease
  if (fF <= fLastAcceptedF + fAlpha * fLambda * fSlope)
    {
      AcceptTrialPoint(trialX);
      if (fLambda == 1 && fNewtonStepLength > 0.99 * fMaxStep)
	fMaxStepTaken = true;
      else
	fMaxStepTaken = false;
      if (NLCS_MESSAGES)
	cout << caller << ": trial point accepted " << endl;
    }
  else
    {
      RejectTrialPoint(trialX);
      if (NLCS_MESSAGES)
	cout << caller << ": trial point rejected " << endl;
    }
}











