///////////////////////////////////////////////////////////////////////////////////////////////////////
//                                   Code: ParaEllip3d-CFD                                           //
//                                 Author: Dr. Beichuan Yan                                          //
//                                  Email: beichuan.yan@colorado.edu                                 //
//                              Institute: University of Colorado Boulder                            //
///////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef STRESS_STRAIN

#ifndef TETRA_H
#define TETRA_H

#include "Vec.h"
#include "Particle.h"
#include "realtypes.h"
#include <map>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

typedef Eigen::Matrix<double, 4, 3> Matrix43d;
typedef Eigen::Matrix<double, 3, 4> Matrix34d;
typedef Eigen::Matrix<double, 3, 1> Matrix31d;

namespace dem {

  class Tetra {

  private:
    int m, n, i, j;    // node number of a tetrahedron
    std::map<int, Particle*> ptclMap; // map between node number and particle
    Matrix43d matrixB; // "strain" matrix, based on shape function N(x), x is current coordinates
    Eigen::Matrix3d  matrix_l;   // velocity gradient, dvdx, using shape function
    Eigen::Matrix3d  matrixF;    // deformation gradient, directly resolved
    Eigen::Matrix3d  matrixFDot; // rate of deformation gradient, directly resolved
    //Eigen::Matrix3d  matrix_l2;  // velocity gradient, dvdx, directly resolved, exactly the same

  public:
    Tetra()
    : m(0), n(0), i(0), j(0) {
      matrixB.setZero();
      matrix_l.setZero();
      matrixF.setZero();
      matrixFDot.setZero();
      //matrix_l2.setZero();
    }
    
    Tetra(int d1, int d2, int d3, int d4, Particle *p1, Particle *p2, Particle *p3, Particle *p4)
    : m(d1), n(d2), i(d3), j(d4) {
      ptclMap[m] = p1;
      ptclMap[n] = p2;
      ptclMap[i] = p3;
      ptclMap[j] = p4;
      // should not compute the tensor here because we need to set node order first
      matrixB.setZero();
      matrix_l.setZero();
      matrixF.setZero();
      matrixFDot.setZero();
      //matrix_l2.setZero();
    }
    
    int getM() const {return m;}
    int getN() const {return n;}
    int getI() const {return i;}
    int getJ() const {return j;}

    // need a const version?
    Particle* getParticleM() {return ptclMap[m];}
    Particle* getParticleN() {return ptclMap[n];}
    Particle* getParticleI() {return ptclMap[i];}
    Particle* getParticleJ() {return ptclMap[j];}

    Matrix43d getMatrixB() const {return matrixB;}
    Eigen::Matrix3d  getMatrix_l() const {return matrix_l;}
    Eigen::Matrix3d  getMatrixF() const {return matrixF;}
    Eigen::Matrix3d  getMatrixFDot() const {return matrixFDot;}

    Vec  getCentroid();
    REAL getVolume();
    REAL getInitVolume();
    void setNodeOrderCalcMatrix();
    void calcMatrixB();
    void calcMatrix_l();
    //void calcMatrix_l2();
    void calcMatrixF();
  };

} // namespace dem

#endif

#endif
