/*
  File: NewtonMethodBase.h
*/

#ifndef _NEWTON_METHOD_BASE_H_
#define _NEWTON_METHOD_BASE_H_

class NLCSolver;
class dArrayT;

class NewtonMethodBase
{
 public:

  // copy constructor
  NewtonMethodBase();

  // solve Ax = b
  void GetNewtonStep(NLCSolver& nlcsolve, dArrayT& step) const;

  // make a copy of *this and returns a pointer to it  
  virtual NewtonMethodBase* clone() const = 0;
};

#endif /* _NEWTON_METHOD_BASE_H_ */
