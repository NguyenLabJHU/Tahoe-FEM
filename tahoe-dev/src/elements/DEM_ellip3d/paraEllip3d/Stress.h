///////////////////////////////////////////////////////////////////////////////////////////////////////
//                                   Code: ParaEllip3d-CFD                                           //
//                                 Author: Dr. Beichuan Yan                                          //
//                                  Email: beichuan.yan@colorado.edu                                 //
//                              Institute: University of Colorado Boulder                            //
///////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef STRESS_STRAIN

#ifndef STRESS_H
#define STRESS_H

#include "realtypes.h"
#include <iostream>
#include <iomanip>
#include <boost/mpi.hpp>

namespace dem {
  
  // This class is only for parallel gathering and printing, because Eigen cannot be used by Boost C++.
  // Matrix operations should use Eigen in Assembly.cpp
  class Stress {
    
  public:
    Stress();
    void setZero();
    void print(std::ostream &ofs) const;
    
  public:
    REAL coord[3];
    REAL density;
    REAL voidRatio;

    REAL stress[6];           // symmetric
    REAL stressRate[6];       // symmetric

    REAL OldroStressRate[9];  // unsymmetric
    REAL TruesStressRate[9];  // unsymmetric

    REAL deformGradient[9];   // unsymmetric
    REAL rotation[9];         // unsymmetric
    REAL stretch[6];          // symmetric
    REAL greenStrain[6];      // symmetric
    REAL eulerStrain[6];      // symmetric
    REAL Jacobian;

    REAL velocityGradient[9]; // unsymmetric
    REAL rateOfDeform[6];     // symmetric
    REAL spin[3];             // skew-symmetric

    REAL stressEigenValue[3];
    REAL stressEigenVector[9];
    REAL stressRateEigenValue[3];
    REAL stressRateEigenVector[9];
    REAL rateOfDeformEigenValue[3];
    REAL rateOfDeformEigenVector[9];

    REAL unitVec[3];
    REAL angle;

    REAL EDot[6];             // symmetric
    REAL eDot[6];             // symmetric
    REAL JDot;
    
    REAL OldroEigenValue[3];
    REAL OldroEigenVector[9];
    REAL TruesEigenValue[3];
    REAL TruesEigenVector[9];
    REAL EDotEigenValue[3];
    REAL EDotEigenVector[9];
    REAL eDotEigenValue[3];
    REAL eDotEigenVector[9];

    REAL fabric[6];           // symmetric

    REAL norm[12];            // norm of 12 tensors

    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & coord;
      ar & density;
      ar & voidRatio;
      ar & stress;
      ar & stressRate;
      ar & OldroStressRate;
      ar & TruesStressRate;
      ar & deformGradient;
      ar & rotation;
      ar & stretch;
      ar & greenStrain;
      ar & eulerStrain;
      ar & Jacobian;
      ar & velocityGradient;
      ar & rateOfDeform;
      ar & spin;
      ar & stressEigenValue;
      ar & stressEigenVector;
      ar & stressRateEigenValue;
      ar & stressRateEigenVector;
      ar & rateOfDeformEigenValue;
      ar & rateOfDeformEigenVector;
      ar & unitVec;
      ar & angle;
      ar & EDot;
      ar & eDot;
      ar & JDot;
      ar & OldroEigenValue;
      ar & OldroEigenVector;
      ar & TruesEigenValue;
      ar & TruesEigenVector;
      ar & EDotEigenValue;
      ar & EDotEigenVector;
      ar & eDotEigenValue;
      ar & eDotEigenVector;
      ar & fabric;
      ar & norm;
    }
    
  };
  
} // namespace dem

#endif

#endif
