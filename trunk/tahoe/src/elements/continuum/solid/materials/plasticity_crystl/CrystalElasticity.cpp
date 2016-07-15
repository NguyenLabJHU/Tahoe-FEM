/*
  File: CrystalElasticity.cpp
*/

#include "CrystalElasticity.h"
#include "PolyCrystalMatT.h"

#include <iostream>
#include "dMatrixT.h"
#include "dArrayT.h"

#include "Utils.h"

// base class

using namespace Tahoe;

CrystalElasticity::CrystalElasticity() 
{
  // initialize general stiffness coefficients
    fC11 = fC12 = fC44 = 0.;
}

CrystalElasticity::~CrystalElasticity() { }

void CrystalElasticity::ElasticityProps(dArrayT& matprop) const
{
  // material properties used in computations
  matprop[0] = fC44;                      // mu 
  matprop[1] = fC12;                      // lambda
  matprop[2] = 0.5*(2.0*fC44+fC12-fC11);  // beta (=0 for isotropic)
}

void CrystalElasticity::ComputeModuli(dMatrixT& moduli) const
{
  // initialize moduli
  moduli = 0.0;

  // 3D symmetric Cij reduced index matrix
  moduli(2,2) = moduli(1,1) = moduli(0,0) = fC11;
  moduli(1,2) = moduli(0,1) = moduli(0,2) = fC12;
  moduli(5,5) = moduli(4,4) = moduli(3,3) = fC44;

  // symmetric
  moduli.CopySymmetric();
}

bool CrystalElasticity::IsIsotropic() const { return false; }

// derived class: IsotropicCrystalElast
IsotropicCrystalElast::IsotropicCrystalElast(PolyCrystalMatT& poly)
{
  // input file
  ifstreamT& in = poly.Input_x();

  // read material constants
  in >> fYoung;	
  in >> fPoisson;

  // check input values
  if (fYoung < 0.0) 
    throwRunTimeError("IsotropicCrystalElast: fYoung < 0");
  if (fPoisson > 0.5 || fPoisson < -1.0) 
    throwRunTimeError("IsotropicCrystalElast: 0.5 < fPoisson < -1");

  // general stiffness coefficients
  fC44 = 0.5*fYoung/(1.0 + fPoisson);              // mu 
  fC12 = 2.0*fC44*fPoisson/(1.0 - 2.0*fPoisson);   // lambda
  fC11 = 2.0*fC44 + fC12;                          // alpha
}

IsotropicCrystalElast::~IsotropicCrystalElast() { }

bool IsotropicCrystalElast::IsIsotropic() const { return true; }

//derived class: CubicCrystalElast
CubicCrystalElast::CubicCrystalElast(PolyCrystalMatT& poly)
{
  // input file
  ifstreamT& in = poly.Input_x();

  // read material constants
  in >> fC11 >> fC12 >> fC44;	

  // check input values
  bool test = ( fC44 > 0.0 && fC11 > fabs(fC12) && (fC11+2.*fC12) > 0.0 );
  if (!test) 
    throwRunTimeError("CubicCrystalElast: Bad values of fC11, fC12, fC44");
}  

CubicCrystalElast::~CubicCrystalElast() { }

//derived class: PECrystalElast
PECrystalElast::PECrystalElast(PolyCrystalMatT& poly)
{
    // initialize general stiffness coefficients
    fC11 = fC22 = fC33 = 0.;  //E11, E22, E33
    fC12  = fC13 = fC23 = 0.; // nu12, nu13, nu23
    fC44=fC55= fC66 = 0.; //G23, G13, G12
    

    // input file
    ifstreamT& in = poly.Input_x();
    
    // read material constants
    in >> fC11 >> fC22 >> fC33 >> fC12 >> fC13 >> fC23 >> fC44 >> fC55 >> fC66;
    /*convention: 12 -> 1122, 44-> 2323, 55-> 1313, 66-> 1212*/
    
}

PECrystalElast::~PECrystalElast() { }

void PECrystalElast::ElasticityProps(dArrayT& matprop) const
{
    // material properties used in computations
    matprop[0] = fC44;                      // mu
    matprop[1] = fC12;                      // lambda
    matprop[2] = 0.5*(2.0*fC44+fC12-fC11);  // beta (=0 for isotropic)
}

void PECrystalElast::ComputeModuli(dMatrixT& moduli) const
{
    // initialize moduli
    moduli = 0.0;
    
    // 3D symmetric Sij compliance index matrix
    moduli(0,0) = 1.0/fC11;
    moduli(0,1) = -fC12/fC22; // First lines are Aurelie's original.
//    moduli(0,1) = -fC12/fC11;
    
    moduli(1,1) = 1.0/fC22;
    moduli(1,2) = -fC23/fC33;
//    moduli(1,2) = -fC23/fC22;
    
    moduli(2,2) = 1.0/fC33;
    moduli(0,2) = -fC13/fC11;
//    moduli(0,2) = -fC13/fC33;
    
    moduli(3,3) = 1.0/fC44;
    moduli(4,4) = 1.0/fC55;
    moduli(5,5) = 1.0/fC66;
    // symmetric
    moduli.CopySymmetric();
    //3D stiffness matrix.
    moduli.Inverse();
//    cout << "\nmoduli: "<< moduli;
}

bool PECrystalElast::IsIsotropic() const { return false; }

