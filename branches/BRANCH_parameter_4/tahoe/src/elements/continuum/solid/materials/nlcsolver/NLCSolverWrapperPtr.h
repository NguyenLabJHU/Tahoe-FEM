/*
  File: NLCSolverWrapperPtr.h
*/

#ifndef _NLC_SOLVER_WRAPPER_PTR_H_
#define _NLC_SOLVER_WRAPPER_PTR_H_

#include "NLCSolverWrapper.h"


namespace Tahoe {

class dArrayT;
class dMatrixT;

/** handle class for NLCSolverWrapper objects */

class NLCSolverWrapperPtr
{
 public:
  NLCSolverWrapperPtr(NLCSolverWrapper* ptr);

  NLCSolverWrapperPtr(const NLCSolverWrapperPtr& other);

  ~NLCSolverWrapperPtr();

  NLCSolverWrapperPtr& operator=(const NLCSolverWrapperPtr& other);

  void FormRHS(dArrayT& x, dArrayT& rhs);

  void FormLHS(dArrayT& x, dMatrixT& lhs);

 private:
  int* refCount_;
  NLCSolverWrapper* ptr_;
};

inline void NLCSolverWrapperPtr::FormRHS(dArrayT& x, dArrayT& rhs)
{ ptr_->FormRHS(x, rhs); }

inline void NLCSolverWrapperPtr::FormLHS(dArrayT& x, dMatrixT& lhs)
{ ptr_->FormLHS(x, lhs); }
 
} // namespace Tahoe 
#endif /* _NLC_SOLVER_WRAPPER_PTR_H_ */
