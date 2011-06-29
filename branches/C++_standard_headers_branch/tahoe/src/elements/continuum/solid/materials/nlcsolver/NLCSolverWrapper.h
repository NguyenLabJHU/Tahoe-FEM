/*
  File: NLCSolverWrapper.h
*/

#ifndef _NLC_SOLVER_WRAPPER_H_
#define _NLC_SOLVER_WRAPPER_H_


namespace Tahoe {

class dArrayT;
class dMatrixT;
class PolyCrystalMatT;
class SlipHardening;
class EVPFDBaseT;

/**
   abstract base class to wrap non-related 
   objects for using the NLCSolver class 
*/

class NLCSolverWrapper
{
 public:
  virtual ~NLCSolverWrapper() = 0;

  virtual void FormRHS(dArrayT& x, dArrayT& rhs) = 0;

  virtual void FormLHS(dArrayT& x, dMatrixT& lhs) = 0;
};

/*
  NLCSolver Wrapper for Polycrystal objects 
*/

class SolverWrapperPoly: public NLCSolverWrapper
{
 public:
  SolverWrapperPoly(PolyCrystalMatT& poly);

  ~SolverWrapperPoly();

  virtual void FormRHS(dArrayT& x, dArrayT& rhs);

  virtual void FormLHS(dArrayT& x, dMatrixT& lhs);

 private:
  PolyCrystalMatT& fpoly;
};

/*
  NLCSolver Wrapper for SlipHardening objects 
*/

class SolverWrapperHard: public NLCSolverWrapper
{
 public:
  SolverWrapperHard(SlipHardening& hard);
  
  ~SolverWrapperHard();

  virtual void FormRHS(dArrayT& x, dArrayT& rhs);

  virtual void FormLHS(dArrayT& x, dMatrixT& lhs);

 private:
  SlipHardening& fhard;
};

/*
  NLCSolver Wrapper for EVPFDBaseT objects 
*/

class SolverWrapperEVPBase: public NLCSolverWrapper
{
 public:
  SolverWrapperEVPBase(EVPFDBaseT& evp);
  
  ~SolverWrapperEVPBase();

  virtual void FormRHS(dArrayT& x, dArrayT& rhs);

  virtual void FormLHS(dArrayT& x, dMatrixT& lhs);

 private:
  EVPFDBaseT& fevp;
};

} // namespace Tahoe 
#endif /* _NLC_SOLVER_WRAPPER_H_ */
