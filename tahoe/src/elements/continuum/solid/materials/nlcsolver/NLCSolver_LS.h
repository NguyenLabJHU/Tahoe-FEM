/*
  File: NLCSolver_LS.h
*/

#ifndef _NLC_SOLVER_LS_H_
#define _NLC_SOLVER_LS_H_

#include "NLCSolver.h"
#include "NewtonMethod.h"


namespace Tahoe {

class dArrayT;

class NLCSolver_LS: public NLCSolver
{
 public:

  // constructor
  NLCSolver_LS(const int dim, const NewtonMethodBase& method = NewtonMethod());

  // destuctor
  ~NLCSolver_LS();
  
  // reset bool variables and counters
  virtual void ResetMemberData();

  // trial solution computations
  virtual void ComputeTrialPoint(dArrayT& X);
  virtual void TestTrialPoint(dArrayT& X);

  // variable in function decrease testing
  void SetAlpha(double alpha);

 private:

  double fLambda;
  double fMinLambda;
  double fAlpha;
  double fSlope;
  double fPrevF;
  double fPrevLambda;
};

} // namespace Tahoe 
#endif  /* _NLC_SOLVER_LS_H_ */
