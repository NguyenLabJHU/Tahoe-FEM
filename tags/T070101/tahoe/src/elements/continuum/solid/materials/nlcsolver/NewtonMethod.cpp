/* 
   File: NewtonMethod.cpp
*/

#include "NewtonMethod.h"
#include "ExceptionCodes.h"

NewtonMethod::NewtonMethod() { }

NewtonMethodBase* NewtonMethod::clone() const
{
  NewtonMethodBase* rtn = new NewtonMethod(*this);
  if (!rtn) throw eOutOfMemory;
  return rtn;
}