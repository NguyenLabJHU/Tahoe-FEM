/*
  File NLCSolverWrapperPtr.cpp
*/

#include "NLCSolverWrapperPtr.h"
#include "ExceptionT.h"


using namespace Tahoe;

NLCSolverWrapperPtr::NLCSolverWrapperPtr(NLCSolverWrapper* ptr)
  : refCount_(0), ptr_(ptr)
{
  refCount_ = new int;
  if (refCount_ == 0) throw ExceptionT::kOutOfMemory;
  *refCount_ = 1;
}

NLCSolverWrapperPtr::NLCSolverWrapperPtr(const NLCSolverWrapperPtr& other)
  : refCount_(other.refCount_),
    ptr_(other.ptr_)
{
  if (refCount_ != 0) (*refCount_)++;
}

NLCSolverWrapperPtr::~NLCSolverWrapperPtr()
{
  if (refCount_ != 0 && --(*refCount_) == 0)
    {
      delete ptr_;
      ptr_ = 0;
      delete refCount_;
      refCount_ = 0;
    }
}

NLCSolverWrapperPtr& NLCSolverWrapperPtr::operator=(const NLCSolverWrapperPtr& other)
{
  if (ptr_ != other.ptr_)
    {
      if (refCount_ != 0 && --(*refCount_) == 0)
	{
	  delete ptr_;
	  delete refCount_;
	}
      ptr_ = other.ptr_;
      refCount_ = other.refCount_;
      if (refCount_ != 0) ++(*refCount_);
    }
  return *this;
}
