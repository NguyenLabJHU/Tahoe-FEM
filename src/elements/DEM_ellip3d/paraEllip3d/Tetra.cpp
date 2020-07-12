///////////////////////////////////////////////////////////////////////////////////////////////////////
//                                   Code: ParaEllip3d-CFD                                           //
//                                 Author: Dr. Beichuan Yan                                          //
//                                  Email: beichuan.yan@colorado.edu                                 //
//                              Institute: University of Colorado Boulder                            //
///////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef STRESS_STRAIN

#include "Tetra.h"
#include "const.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>

namespace dem {

  REAL Tetra::getInitVolume() {
    REAL x1, y1, z1;
    REAL x2, y2, z2;
    REAL x3, y3, z3;
    REAL x4, y4, z4;
    x1 = ptclMap[m]->getInitPos().getX();
    y1 = ptclMap[m]->getInitPos().getY();
    z1 = ptclMap[m]->getInitPos().getZ();
    x2 = ptclMap[n]->getInitPos().getX();
    y2 = ptclMap[n]->getInitPos().getY();
    z2 = ptclMap[n]->getInitPos().getZ();
    x3 = ptclMap[i]->getInitPos().getX();
    y3 = ptclMap[i]->getInitPos().getY();
    z3 = ptclMap[i]->getInitPos().getZ();
    x4 = ptclMap[j]->getInitPos().getX();
    y4 = ptclMap[j]->getInitPos().getY();
    z4 = ptclMap[j]->getInitPos().getZ();

    return (x2*y3*z4-x2*y4*z3-x3*y2*z4+x3*y4*z2+x4*y2*z3-x4*y3*z2
	    -x1*y3*z4+x1*y4*z3+x3*y1*z4-x3*y4*z1-x4*y1*z3+x4*y3*z1
	    +x1*y2*z4-x1*y4*z2-x2*y1*z4+x2*y4*z1+x4*y1*z2-x4*y2*z1
	    -x1*y2*z3+x1*y3*z2+x2*y1*z3-x2*y3*z1-x3*y1*z2+x3*y2*z1)/6.0;
  }


  REAL Tetra::getVolume() {
    REAL x1, y1, z1;
    REAL x2, y2, z2;
    REAL x3, y3, z3;
    REAL x4, y4, z4;
    x1 = ptclMap[m]->getCurrPos().getX();
    y1 = ptclMap[m]->getCurrPos().getY();
    z1 = ptclMap[m]->getCurrPos().getZ();
    x2 = ptclMap[n]->getCurrPos().getX();
    y2 = ptclMap[n]->getCurrPos().getY();
    z2 = ptclMap[n]->getCurrPos().getZ();
    x3 = ptclMap[i]->getCurrPos().getX();
    y3 = ptclMap[i]->getCurrPos().getY();
    z3 = ptclMap[i]->getCurrPos().getZ();
    x4 = ptclMap[j]->getCurrPos().getX();
    y4 = ptclMap[j]->getCurrPos().getY();
    z4 = ptclMap[j]->getCurrPos().getZ();

    return (x2*y3*z4-x2*y4*z3-x3*y2*z4+x3*y4*z2+x4*y2*z3-x4*y3*z2
	    -x1*y3*z4+x1*y4*z3+x3*y1*z4-x3*y4*z1-x4*y1*z3+x4*y3*z1
	    +x1*y2*z4-x1*y4*z2-x2*y1*z4+x2*y4*z1+x4*y1*z2-x4*y2*z1
	    -x1*y2*z3+x1*y3*z2+x2*y1*z3-x2*y3*z1-x3*y1*z2+x3*y2*z1)/6.0;
  }


  Vec Tetra::getCentroid() {
    REAL x1, y1, z1;
    REAL x2, y2, z2;
    REAL x3, y3, z3;
    REAL x4, y4, z4;
    x1 = ptclMap[m]->getCurrPos().getX();
    y1 = ptclMap[m]->getCurrPos().getY();
    z1 = ptclMap[m]->getCurrPos().getZ();
    x2 = ptclMap[n]->getCurrPos().getX();
    y2 = ptclMap[n]->getCurrPos().getY();
    z2 = ptclMap[n]->getCurrPos().getZ();
    x3 = ptclMap[i]->getCurrPos().getX();
    y3 = ptclMap[i]->getCurrPos().getY();
    z3 = ptclMap[i]->getCurrPos().getZ();
    x4 = ptclMap[j]->getCurrPos().getX();
    y4 = ptclMap[j]->getCurrPos().getY();
    z4 = ptclMap[j]->getCurrPos().getZ();
    
    return(Vec((x1+x2+x3+x4)/4, (y1+y2+y3+y4)/4, (z1+z2+z3+z4)/4));
  }


  void Tetra::setNodeOrderCalcMatrix() {
    if(getVolume() < 0) { // swap node 2 and node 3
      /*
      std::cout << "tetra nodes are not numbered counter-clockwise, volume=" << getVolume() << " iteration=" << iteration << std::endl
		<< std::setw(OWID) << m << std::setw(OWID) << n << std::setw(OWID) << i << std::setw(OWID) << j << std::endl
		<< std::setw(OWID) << ptclMap[m]->getCurrPos().getX() << std::setw(OWID) << ptclMap[m]->getCurrPos().getY() << std::setw(OWID) << ptclMap[m]->getCurrPos().getZ() << std::endl
		<< std::setw(OWID) << ptclMap[n]->getCurrPos().getX() << std::setw(OWID) << ptclMap[n]->getCurrPos().getY() << std::setw(OWID) << ptclMap[n]->getCurrPos().getZ() << std::endl
		<< std::setw(OWID) << ptclMap[i]->getCurrPos().getX() << std::setw(OWID) << ptclMap[i]->getCurrPos().getY() << std::setw(OWID) << ptclMap[i]->getCurrPos().getZ() << std::endl
		<< std::setw(OWID) << ptclMap[j]->getCurrPos().getX() << std::setw(OWID) << ptclMap[j]->getCurrPos().getY() << std::setw(OWID) << ptclMap[j]->getCurrPos().getZ() << std::endl;
      */
      int iTmp = n;
      n = i;
      i = iTmp;
    }

    calcMatrixB();
    calcMatrix_l();
    calcMatrixF();
    //calcMatrix_l2();

    /*
    std::cout << "iteration=" << iteration 
	      << " /////////////////////////////////////////" << std::endl
	      << "matrixF=" << std::endl << matrixF << std::endl << std::endl
	      << "matrixFDot=" << std::endl << matrixFDot << std::endl << std::endl
	      << "matrix_l=" << std::endl << matrix_l << std::endl << std::endl
	      << "matrix_l2=" << std::endl << matrix_l2 << std::endl << std::endl;
    */
  }


  void Tetra::calcMatrixB() { // based on current configuration x, not initial configuration X, B=B(x)
    REAL x1, y1, z1;
    REAL x2, y2, z2;
    REAL x3, y3, z3;
    REAL x4, y4, z4;
    x1 = ptclMap[m]->getCurrPos().getX();
    y1 = ptclMap[m]->getCurrPos().getY();
    z1 = ptclMap[m]->getCurrPos().getZ();
    x2 = ptclMap[n]->getCurrPos().getX();
    y2 = ptclMap[n]->getCurrPos().getY();
    z2 = ptclMap[n]->getCurrPos().getZ();
    x3 = ptclMap[i]->getCurrPos().getX();
    y3 = ptclMap[i]->getCurrPos().getY();
    z3 = ptclMap[i]->getCurrPos().getZ();
    x4 = ptclMap[j]->getCurrPos().getX();
    y4 = ptclMap[j]->getCurrPos().getY();
    z4 = ptclMap[j]->getCurrPos().getZ();

    // intermediate variables
    REAL a1, a2, a3, a4;
    REAL b1, b2, b3, b4;
    REAL c1, c2, c3, c4;
    a1 = y2*(z4-z3)-y3*(z4-z2)+y4*(z3-z2);
    a2 = -y1*(z4-z3)+y3*(z4-z1)-y4*(z3-z1);
    a3 = y1*(z4-z2)-y2*(z4-z1)+y4*(z2-z1);
    a4 = -y1*(z3-z2)+y2*(z3-z1)-y3*(z2-z1);

    b1 = -x2*(z4-z3)+x3*(z4-z2)-x4*(z3-z2);
    b2 = x1*(z4-z3)-x3*(z4-z1)+x4*(z3-z1);
    b3 = -x1*(z4-z2)+x2*(z4-z1)-x4*(z2-z1);
    b4 = x1*(z3-z2)-x2*(z3-z1)+x3*(z2-z1);

    c1 = x2*(y4-y3)-x3*(y4-y2)+x4*(y3-y2);
    c2 = -x1*(y4-y3)+x3*(y4-y1)-x4*(y3-y1);
    c3 = x1*(y4-y2)-x2*(y4-y1)+x4*(y2-y1);
    c4 = -x1*(y3-y2)+x2*(y3-y1)-x3*(y2-y1);

    matrixB(0,0) = a1;
    matrixB(1,0) = a2;
    matrixB(2,0) = a3;
    matrixB(3,0) = a4;

    matrixB(0,1) = b1;
    matrixB(1,1) = b2;
    matrixB(2,1) = b3;
    matrixB(3,1) = b4;
	
    matrixB(0,2) = c1;
    matrixB(1,2) = c2;
    matrixB(2,2) = c3;
    matrixB(3,2) = c4;

    // current tetrahedron could have coplanar/collinear vertices
    if (getVolume() == 0)
      matrixB.setZero();
    else
      matrixB /= (6.0 * getVolume());
  }


  void Tetra::calcMatrix_l() {
    Matrix31d v1(3,1), v2(3,1), v3(3,1), v4(3,1);

    v1 << ptclMap[m]->getCurrVeloc().getX(), ptclMap[m]->getCurrVeloc().getY(), ptclMap[m]->getCurrVeloc().getZ();
    v2 << ptclMap[n]->getCurrVeloc().getX(), ptclMap[n]->getCurrVeloc().getY(), ptclMap[n]->getCurrVeloc().getZ();
    v3 << ptclMap[i]->getCurrVeloc().getX(), ptclMap[i]->getCurrVeloc().getY(), ptclMap[i]->getCurrVeloc().getZ();
    v4 << ptclMap[j]->getCurrVeloc().getX(), ptclMap[j]->getCurrVeloc().getY(), ptclMap[j]->getCurrVeloc().getZ();

    Matrix34d velocity;
    velocity(0,0) = v1(0);
    velocity(1,0) = v1(1);
    velocity(2,0) = v1(2);

    velocity(0,1) = v2(0);
    velocity(1,1) = v2(1);
    velocity(2,1) = v2(2);

    velocity(0,2) = v3(0);
    velocity(1,2) = v3(1);
    velocity(2,2) = v3(2);

    velocity(0,3) = v4(0);
    velocity(1,3) = v4(1);
    velocity(2,3) = v4(2);

    matrix_l = velocity * matrixB; // both velocity and matrixB are current configurations
  }


  void Tetra::calcMatrixF() { // dx = F dX, a linear solver for F based on any three independent edges of a tetrahedron
    // three edge vectors: jm, jn, ji
    Eigen::Matrix3d dXdYdZ; // in Ax=b, this is A
    Eigen::Vector3d dxdxdx, dydydy, dzdzdz; // in Ax=b, this is b
    Eigen::Vector3d pdx, pdy, pdz; // in Ax=b, this is x; partial differential
    dXdYdZ.setZero();
    dxdxdx.setZero();
    dydydy.setZero();
    dzdzdz.setZero();
    pdx.setZero();
    pdy.setZero();
    pdz.setZero();

    // for dX
    Vec JM = ptclMap[m]->getInitPos() - ptclMap[j]->getInitPos();
    Vec JN = ptclMap[n]->getInitPos() - ptclMap[j]->getInitPos();
    Vec JI = ptclMap[i]->getInitPos() - ptclMap[j]->getInitPos();

    dXdYdZ << JM.getX(), JM.getY(), JM.getZ(),
              JN.getX(), JN.getY(), JN.getZ(),
              JI.getX(), JI.getY(), JI.getZ();

    /*
    if (fabs(dXdYdZ.determinant()) < EPS ) {
      std::cout << "iteration=" << iteration << " /////////////////////////////////////////" << std::endl;
      std::cout<<  "m="<<m<<"\n init=";ptclMap[m]->getInitPos().print(std::cout);std::cout<<"\n curr=";ptclMap[m]->getCurrPos().print(std::cout);
      std::cout<<"\nn="<<n<<"\n init=";ptclMap[n]->getInitPos().print(std::cout);std::cout<<"\n curr=";ptclMap[n]->getCurrPos().print(std::cout);
      std::cout<<"\ni="<<i<<"\n init=";ptclMap[i]->getInitPos().print(std::cout);std::cout<<"\n curr=";ptclMap[i]->getCurrPos().print(std::cout);
      std::cout<<"\nj="<<j<<"\n init=";ptclMap[j]->getInitPos().print(std::cout);std::cout<<"\n curr=";ptclMap[j]->getCurrPos().print(std::cout);
      std::cout<<" volume="<<getVolume()<<" initVolume="<<getInitVolume()<<" det=" << dXdYdZ.determinant();
      std::cout<< std::endl;
    }
    */

    // if initial tetrahedron have coplanar/collinear vertices, the determinant is close to zero, and
    // this condition does not work well (V=1/6|det|) and may generate "nan".
    if (fabs(dXdYdZ.determinant()) > EPS) { 
      Eigen::Matrix3d inv =  dXdYdZ.inverse(); // more efficient than LU decomposition for very small matrix

      // for dx
      Vec jm = ptclMap[m]->getCurrPos() - ptclMap[j]->getCurrPos();
      Vec jn = ptclMap[n]->getCurrPos() - ptclMap[j]->getCurrPos();
      Vec ji = ptclMap[i]->getCurrPos() - ptclMap[j]->getCurrPos();

      dxdxdx << jm.getX(), jn.getX(), ji.getX();
      pdx = inv * dxdxdx;

      dydydy << jm.getY(), jn.getY(), ji.getY();
      pdy = inv * dydydy;

      dzdzdz << jm.getZ(), jn.getZ(), ji.getZ();
      pdz = inv * dzdzdz;

      matrixF << pdx.transpose(), pdy.transpose(), pdz.transpose();

      // for dv
      jm = ptclMap[m]->getCurrVeloc() - ptclMap[j]->getCurrVeloc();
      jn = ptclMap[n]->getCurrVeloc() - ptclMap[j]->getCurrVeloc();
      ji = ptclMap[i]->getCurrVeloc() - ptclMap[j]->getCurrVeloc();

      dxdxdx << jm.getX(), jn.getX(), ji.getX();
      pdx = inv * dxdxdx;

      dydydy << jm.getY(), jn.getY(), ji.getY();
      pdy = inv * dydydy;

      dzdzdz << jm.getZ(), jn.getZ(), ji.getZ();
      pdz = inv * dzdzdz;

      matrixFDot << pdx.transpose(), pdy.transpose(), pdz.transpose();
    }

  }

  /*
  void Tetra::calcMatrix_l2() {
    // three edge vectors: jm, jn, ji
    Eigen::Matrix3d dxdydz; // in Ax=b, this is A
    Eigen::Vector3d dvxdvxdvx, dvydvydvy, dvzdvzdvz; // in Ax=b, this is b
    Eigen::Vector3d pdx, pdy, pdz; // in Ax=b, this is x; partial differential
    dxdydz.setZero();
    dvxdvxdvx.setZero();
    dvydvydvy.setZero();
    dvzdvzdvz.setZero();
    pdx.setZero();
    pdy.setZero();
    pdz.setZero();

    // for dx
    Vec jm = ptclMap[m]->getCurrPos() - ptclMap[j]->getCurrPos();
    Vec jn = ptclMap[n]->getCurrPos() - ptclMap[j]->getCurrPos();
    Vec ji = ptclMap[i]->getCurrPos() - ptclMap[j]->getCurrPos();

    dxdydz << jm.getX(), jm.getY(), jm.getZ(),
              jn.getX(), jn.getY(), jn.getZ(),
              ji.getX(), ji.getY(), ji.getZ();

    Eigen::Matrix3d inv =  dxdydz.inverse(); // more efficient than LU decomposition for very small matrix

    // for dvx
    Vec dvjm = ptclMap[m]->getCurrVeloc() - ptclMap[j]->getCurrVeloc();
    Vec dvjn = ptclMap[n]->getCurrVeloc() - ptclMap[j]->getCurrVeloc();
    Vec dvji = ptclMap[i]->getCurrVeloc() - ptclMap[j]->getCurrVeloc();
      
    dvxdvxdvx << dvjm.getX(), dvjn.getX(), dvji.getX();
    pdx = inv * dvxdvxdvx;

    dvydvydvy << dvjm.getY(), dvjn.getY(), dvji.getY();
    pdy = inv * dvydvydvy;

    dvzdvzdvz << dvjm.getZ(), dvjn.getZ(), dvji.getZ();
    pdz = inv * dvzdvzdvz;

    matrix_l2 << pdx.transpose(), pdy.transpose(), pdz.transpose();
  }
  */
}

#endif
