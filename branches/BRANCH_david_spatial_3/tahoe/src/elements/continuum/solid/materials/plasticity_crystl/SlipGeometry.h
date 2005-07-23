/*
  File: SlipGeometry.h
*/

#ifndef _SLIP_GEOMETRY_H_
#define _SLIP_GEOMETRY_H_

#include <iostream.h>

#include "ArrayT.h"
#include "dArrayT.h"
#include "dMatrixT.h"

namespace Tahoe {

/* Base Class */
class SlipGeometry
{
 public:
  // enum variable to select crystal structure
  enum CrystalType { kFCC=1, kBCC=2, kHCP=3, kFCC24=4 };

  // constructor
  SlipGeometry(int numslip);

  // destructor
  virtual ~SlipGeometry();

  // accesors
  int NumSlip() const;
  const ArrayT<dArrayT>& GetSlipVectorS() const;
  const ArrayT<dArrayT>& GetSlipVectorM() const;
  const ArrayT<dMatrixT>& GetSchmidtTensor() const;

  // print crystal structure type
   virtual void PrintName(ostream& out) const = 0;

 protected:
  // compute slip quantities
  void InitializeSlipQnts();

 private:
  // normalize slip system vectors
  void NormalizeSlipVectors();    

  // compute Schmidt tensor s(x)m
  void SetSchmidtTensor();

 protected:
  // number of slip systems
  int fNumSlip;

  // slip system vectors
  ArrayT<dArrayT> fVecS;
  ArrayT<dArrayT> fVecM;

  // Schmidt tensor: fVecS(x)fVecM
  ArrayT<dMatrixT> fZ;
};

inline int SlipGeometry::NumSlip() const { return fNumSlip; } 
inline const ArrayT<dArrayT>& SlipGeometry::GetSlipVectorS() const
  { return fVecS; }
inline const ArrayT<dArrayT>& SlipGeometry::GetSlipVectorM() const
  { return fVecM; }
inline const ArrayT<dMatrixT>& SlipGeometry::GetSchmidtTensor() const
  { return fZ; }

/* Derived Classes */

class FCCGeometry: public SlipGeometry
{
 public:
 // constructor
  FCCGeometry();

  // destructor
  ~FCCGeometry(); 
  
  // print crystal structure name
  virtual void PrintName(ostream& out) const;

 private:
  // build FCC slip system vectors from Miller indices
  void SetSlipVectors();
};

class BCCGeometry: public SlipGeometry
{
 public:
  // constructor
  BCCGeometry();

  // destructor
  ~BCCGeometry();

  // print crystal structure name
  virtual void PrintName(ostream& out) const;

 private:
  // build BCC slip system vectors from Miller indices
  void SetSlipVectors();
};

class HCPGeometry: public SlipGeometry
{
 public:
  // constructor
  HCPGeometry();

  // destructor
  ~HCPGeometry();

  // print crystal structure name
  virtual void PrintName(ostream& out) const;

 private:
  // build HCP slip system vectors from Miller indices
  void SetSlipVectors();
};

class FCCGeometry24: public SlipGeometry
{
 public:
 // constructor
  FCCGeometry24();

  // destructor
  ~FCCGeometry24();

  // print crystal structure name
  virtual void PrintName(ostream& out) const;

 private:
  // build FCC24 slip system vectors from Miller indices
  void SetSlipVectors();
};

} // namespace Tahoe 
#endif /* _SLIP_GEOMETRY_H_ */

