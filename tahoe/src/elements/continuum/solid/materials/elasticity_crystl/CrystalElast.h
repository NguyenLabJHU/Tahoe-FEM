/* $Id: CrystalElast.h,v 1.2 2001-08-20 22:15:40 rdorgan Exp $ */
/*
  File: CrystalElast.h
*/
#ifndef _CRYSTAL_ELAST_H_
#define _CRYSTAL_ELAST_H_

#include "FDHookeanMatT.h"

#include "GlobalT.h"
#include "ArrayT.h"
#include "dArrayT.h"
#include "Array2DT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "ifstreamT.h" 

#include "LocalArrayT.h"

class CrystalElastLat;
class CrystalElastMat;

class ElasticT;

class CrystalElast : public FDHookeanMatT
{
 public:
  // constructor
  CrystalElast(ifstreamT& in, const FiniteStrainT& element);

  // destructor
  virtual ~CrystalElast();

  // allocate space/initialize crystal arrays (all)
  virtual bool NeedsInitialization() const;
  virtual void Initialize();

  // some methods to set/initialize member data
  virtual void InitializeCrystalVariables() = 0;
  virtual int  NumVariablesPerElement() = 0;

  // general accesors
  ifstreamT& Input_x();
  const int NumGrain() const;
  const dArrayT& MaterialProperties() const;

  // print data read
  virtual void Print(ostream& out) const;

 protected:

     /* set (material) tangent modulus */
     virtual void SetModulus(dMatrixT& modulus);

  // print name
  virtual void PrintName(ostream& out) const;

  // allocate all elements at once
  void AllocateElements();

 private:
  // read lattice orientation data, construct array fEuler
  void SetLatticeOrientation() ;

  // crystal elasticity parameters
  void SetCrystalElasticity();

  // solver for nonlinear constitutive equations
  void SetConstitutiveSolver();

 protected:
  // initial and final temperature
  double fInit_Temp_DegC;
  double fTemp_DegC;

  // status flag
  const GlobalT::StateT& fStatus;

  // stream for crystal input data
  ifstreamT fInput;

  // number crystals at each IP
  int fNumGrain;

  // control code for supporting classes
  int fODFCode;

  // pointers to supporting classes
  CrystalElastLat* fCrystalElastLat;
  CrystalElastMat* fCrystalElastMat;

  // material properties
  dArrayT fMatProp;

  // array for Euler angles at integration point
  ArrayT<dArrayT> fangles;

  // huge temporary array to hold all euler angles
  ArrayT<Array2DT<dArrayT> > fEuler; 
};

  // general needed accesors
  inline ifstreamT& CrystalElast::Input_x() { return fInput; }
  inline const int CrystalElast::NumGrain() const { return fNumGrain; }
  inline const dArrayT& CrystalElast::MaterialProperties() const { return fMatProp; }

#endif /* _CRYSTAL_ELAST_H_ */
