/*
  File: NewtonMethod.h
*/

#ifndef _NEWTON_METHOD_H_
#define _NEWTON_METHOD_H_

#include "NewtonMethodBase.h"


namespace Tahoe {

class NewtonMethod : public NewtonMethodBase
{
 public:

  // constructor
  NewtonMethod();

  // virtual constructor
  virtual NewtonMethodBase* clone() const;
};

} // namespace Tahoe 
#endif /* _NEWTON_METHOD_H_ */
