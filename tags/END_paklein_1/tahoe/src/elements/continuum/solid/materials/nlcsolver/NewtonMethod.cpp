/* 
   File: NewtonMethod.cpp
*/

#include "NewtonMethod.h"
#include "ExceptionT.h"


using namespace Tahoe;

NewtonMethod::NewtonMethod() { }

NewtonMethodBase* NewtonMethod::clone() const
{
  NewtonMethodBase* rtn = new NewtonMethod(*this);
  if (!rtn) throw ExceptionT::kOutOfMemory;
  return rtn;
}
