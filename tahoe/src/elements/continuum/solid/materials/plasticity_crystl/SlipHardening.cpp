/*
  File: SlipHardening.cpp
*/

#include "PolyCrystalMatT.h"
#include "SlipHardening.h"
#include "NLCSolver.h"
#include "NLCSolver_LS.h"
#include "Utils.h"
#include "ifstreamT.h" 


using namespace Tahoe;

SlipHardening::SlipHardening(PolyCrystalMatT& poly):
  fPolyXtal  (poly),
  fdt        (poly.TimeStep()),
  fTau       (poly.GetResolvedShearStress()),
  fDGamma    (poly.GetIncrSlipShearStrain()),
  fSolver    (NULL), 
  fSolverPtr ( new SolverWrapperHard(*this) )
{ 
  // if needed, these allocations are overridden in derived classes
  fTauIso.Dimension(1);
  fTauKin.Dimension(1);
  fTauInc.Dimension(1);

  // other allocations are carried out in derived classes
}

SlipHardening::~SlipHardening() { delete fSolver; }

void SlipHardening::SetHardeningSolver(ifstreamT& in, int numvar)
{
  int ksolve;
  in >> ksolve;

  switch(ksolve)
    {
      // Newton's method + line search
    case NLCSolver::kNLCSolver_LS: 
      fSolver = new NLCSolver_LS(numvar);
      break;

      /* // Newton's method + trust region
    case NLCSolver::kNLCSolver_TR:
      fSolver = new NLCSolver_TR(numvar);
      break; */

    default:
      throwRunTimeError("SlipHardening::SetHardeningSolver: Bad ksolve");
    }
  if (!fSolver) throwMemoryError("SlipHardening::SetHardeningSolver");

  // modify some default values of solver
  int maxiter;
  in >> maxiter;
  fSolver->SetMaxIterations(maxiter);

  double functol;
  in >> functol;
  fSolver->SetFuncTol(functol);

  double gradtol;
  in >> gradtol;
  fSolver->SetGradTol(gradtol);
}

// default implementation returns nondirectional hardening stress
const dArrayT& SlipHardening::ComputeHardQnts() { return fTauIso; }
const double SlipHardening::ComputeHardQnts(int is) { return fTauIso[is]; }

const double SlipHardening::IsoHardeningStress(int is) const { return fTauIso[is]; }
const double SlipHardening::KinHardeningStress(int is) const { return fTauKin[is]; }
