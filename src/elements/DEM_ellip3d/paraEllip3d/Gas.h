#ifndef GAS_H
#define GAS_H

#include "MPIFrame.h"
#include "Parameter.h"
#include "realtypes.h"
#include "Vec.h"
#include "Gradation.h"
#include "Rectangle.h"
#include "Boundary.h"
#include "Particle.h"
#include <cstddef>
#include <vector>
#include <utility>
#include <map>
#include <boost/mpi.hpp>


namespace dem {
  
  class IJK {
  public:
    std::size_t i;
    std::size_t j;
    std::size_t k;

    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & i;
      ar & j;
      ar & k;
    }
    
  public:
    IJK() : i(0), j(0), k(0) {}
    IJK(int _i, int _j, int _k) : i(_i), j(_j), k(_k) {}
  };

  struct cmpIJK {
    bool operator()(const IJK &v1, const IJK &v2) const {

      // compare in i, j, k order, important!
      if (v1.i < v2.i) {
	return true;
      } else if (v1.i == v2.i) {

	if (v1.j < v2.j) {
	  return true;
	} else if (v1.j == v2.j) {

	  if (v1.k < v2.k) {
	    return true;
	  } else {
	    return false;
	  }

	} else {
	  return false;
	}

      } else {
	return false;
      }

    }
  }; 

  // class for MPI transmission
  // must transmit conserved variables across discontinuity!
  // with no discontiuity, it is OK transmit primitive variables.
  class GasVar { 
  public:
    REAL density;
    Vec  momentum;
    REAL energy;
    Vec  velocity;
    REAL pressure;
    Vec  coords;

    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & density;
      ar & momentum;
      ar & energy;
      ar & velocity;
      ar & pressure;
      ar & coords;
    }

  public:
  GasVar() : density(0), momentum(0), energy(0), velocity(0), pressure(0), coords(0) {}
  GasVar(REAL d, Vec mo, REAL en, Vec ve, REAL pr, Vec co) 
      : density(d), momentum(mo), energy(en), velocity(ve), pressure(pr), coords(co) {}
  };

  class ValCoord { // class to store a Vec variable and its coordinates
  public:
    Vec  val;
    Vec  coords;

    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & val;
      ar & coords;
    }

  public:
    ValCoord() : val(0), coords(0) {}
    ValCoord(Vec v1, Vec v2) 
      : val(v1), coords(v2) {}
  };

  typedef std::map<IJK, GasVar, cmpIJK> IJK_GasVar_Map;
  typedef std::map<IJK, ValCoord, cmpIJK> IJK_arrayPenalForce_Map;

  class Bound {
  public:
    std::size_t lowX, lowY, lowZ;
    std::size_t uppX, uppY, uppZ;
    std::size_t sizeX, sizeY, sizeZ;

  public:
    Bound() : lowX(0), lowY(0), lowZ(0), uppX(0), uppY(0), uppZ(0), sizeX(0), sizeY(0), sizeZ(0) {}
    Bound(std::size_t lowx, std::size_t lowy, std::size_t lowz, std::size_t uppx, std::size_t uppy, std::size_t uppz) 
      : lowX(lowx), lowY(lowy), lowZ(lowz), 
        uppX(uppx), uppY(uppy), uppZ(uppz), 
        sizeX(uppx - lowx + 1), sizeY(uppy - lowy + 1), sizeZ(uppz - lowz + 1) {}
  };

  class Gas {

    typedef std::vector< std::vector< std::vector <REAL> > > Array3D;
    typedef std::vector< std::vector< std::vector <std::vector<REAL> > > > Array4D;
    typedef std::vector< std::vector< std::vector <std::vector< std::vector<REAL>  > > > > Array5D;

  private:

    MPIFrame mpi;
    
    static const REAL Rs; // specific gas constant

    // only for all gas domain
    std::size_t allGridNx;// allGridNx = total cell centers = parts + two boundary points in x direction
    std::size_t allGridNy;// allGridNy = total cell centers = parts + two boundary points in y direction
    std::size_t allGridNz;// allGridNz = total cell centers = parts + two boundary points in z direction

    // for each gas partition
    std::size_t gridNx;   // gridNx = cell centers = parts + ghost boundary points
    std::size_t gridNy;
    std::size_t gridNz;

    Bound boundPrn; // CFD printing bound
    Bound boundCup; // DEM-CFD coupling bound
    Bound boundCfd; // CFD resolving bound
    Bound boundGod; // God resolving bound

    std::size_t ptclGrid; // approximate grids accross particle in each dimension
    std::size_t haloGrid;// ghost layer thickness for DEM particles

    std::size_t haloGridX;
    std::size_t haloGridY;
    std::size_t haloGridZ;

    REAL gridDx;          // grid size in x direction
    REAL gridDy;          // grid size in y direction
    REAL gridDz;          // grid size in z direction
    REAL x1F, x2F, y1F, y2F, z1F, z2F; // fluid domain

    REAL Cd;           // drag coefficient
    REAL porosity;     // particle porosity as porous media
    REAL Cdi;          // fictitious drag coefficient inside porous media
    REAL RK;           // Runge-Kutta scheme
    REAL CFL;          // Courant-Friedrichs-Lewy condition
    REAL gama;         // ratio of specific heat capacity of air
    REAL arrayBC[6];   // boundary condition

    int  leftType;     // type of left part
    REAL x1L, x2L, y1L, y2L, z1L, z2L; // left part of fluid domain
    REAL x0L;          // center of left part of sphere, x-coordinate
    REAL y0L;          // center of left part of sphere, y-coordinate
    REAL z0L;          // center of left part of sphere, z-coordinate
    REAL r0L;          // radius of left part

    REAL rhoR, uR, pR; // known for Rankine-Hugoniot conditions (RHC)
    REAL MachShock;    // shock Mach number, known for RHC
    REAL MachL;        // Mach number for left part
    REAL rhoL, uL, pL; // unknown for RHC
    REAL shockSpeed;   // shock/discontinuity speed
    REAL rhoBL, uBL, pBL; // below left part

    std::size_t nDim, nVar, nInteg, varDen, varEng, varPrs, varMsk;
    std::size_t varMom[3], varVel[3];

    bool negPrsDen;    // handle cases of negative pressure or density

    Array4D arrayU;
    Array4D arrayURota;  // copy of arrayU in rotation.
    Array4D arrayUN;     // for Runge-Kutta 
    Array4D arrayUStar;  // for Runge-Kutta 
    Array4D arrayUStar2; // for Runge-Kutta 
    // 4-dimensional, defined at cell centers
    // gridNx, gridNy, gridNz, nVar
    // (a) fixed:
    // arrayU[i][j][k][0]: varDen
    // arrayU[i][j][k][1]: varMom[0], in x direction
    // arrayU[i][j][k][2]: varMom[1], in y direction
    // arrayU[i][j][k][3]: varMom[2], in z direction
    // arrayU[i][j][k][4]: varEng    // total energy per unit volume, E = rho * (1/2*V^2 + e), NOT total specific energy
    // arrayU[i][j][k][5]: varVel[0], in x direction, corresponding to u in 3D Euler equations.
    // arrayU[i][j][k][6]: varVel[1], in y direction, corresponding to v in 3D Euler equations.
    // arrayU[i][j][k][7]: varVel[2], in z direction, corresponding to w in 3D Euler equations.
    // arrayU[i][j][k][8]: varPrs
    // (b) extended:
    // arrayU[i][j][k][9]: varMsk

    Array4D arrayGridCoord; 
    // fluid grid coordinates, 4-dimensional
    // gridNx, gridNy, gridNz, nDim
    // arrayGridCoord[i][j][k][0]: coordX
    // arrayGridCoord[i][j][k][1]: coordY
    // arrayGridCoord[i][j][k][2]: coordZ

    Array4D arrayPenalForce;
    // fluid grid forces, 4-dimensional
    // gridNx, gridNy, gridNz, nDim
    // arrayPenalForce[i][j][k][0]: forceX
    // arrayPenalForce[i][j][k][1]: forceY
    // arrayPenalForce[i][j][k][2]: forceZ

    Array4D arrayPressureForce;
    // fluid grid forces, 4-dimensional
    // gridNx, gridNy, gridNz, nDim
    // arrayPressureForce[i][j][k][0]: forceX
    // arrayPressureForce[i][j][k][1]: forceY
    // arrayPressureForce[i][j][k][2]: forceZ

    Array4D arrayFlux;
    // 4-dimensional, defined at cell centers
    // gridNx, gridNy, gridNz, nInteg
    // arrayFlux[i][j][k][0]: varDen
    // arrayFlux[i][j][k][1]: varMom[0]
    // arrayFlux[i][j][k][2]: varMom[1]
    // arrayFlux[i][j][k][3]: varMom[2]
    // arrayFlux[i][j][k][4]: varEng

    Array5D arrayGodFlux; 
    // 5-dimensional, defined at cell faces
    // gridNx-1, gridNy-1, gridNz-1, nInteg, nDim
    // arrayGodFlux[i][j][k][0]: varDen
    // arrayGodFlux[i][j][k][1]: varMom[0]
    // arrayGodFlux[i][j][k][2]: varMom[1]
    // arrayGodFlux[i][j][k][3]: varMom[2]
    // arrayGodFlux[i][j][k][4]: varEng

    Array4D arrayGodFluxTmp;
    // 4-dimensional
    // gridNx-1, gridNy-1, gridNz-1, nInteg

    Array3D arrayH; 
    // Enthalpy, 3-dimensional, total specific enthalpy, NOT static specific enthalpy
    // gridNx, gridNy, gridNz

    Array3D arraySoundSpeed; 
    // speed of sound, 3-dimensional
    // gridNx, gridNy, gridNz

    std::vector<std::size_t> printPtcls;

  public:  
    void setMPI(MPIFrame &m) {mpi = m;}
    void commu6();
    void commu20();
    void backCommu20();
    void findGasInRectangle(IJK &start, IJK &end, std::vector<GasVar> &foundGas);
    void findPenalInRectangle(IJK &start, IJK &end, std::vector<ValCoord> &foundPenal);

  public:
    Gas() {}
    void initParameter(Gradation &gradation);
    void initPureGasParameter();
    void initSharedParameter();
    void printSharedParameter();

    void allocArray();
    void initializePureGas();
    void initialize();
    void initialCondition();
    void calcTimeStep();
    void RankineHugoniot();
    void initGhostPoints();
    void soundSpeed();
    void enthalpy();
    void flux(std::size_t, std::vector<Particle *> &ptcls);
    void LaxFrieScheme(REAL UL[], REAL UR[], REAL FL[], REAL FR[], std::size_t iDim, std::size_t it, std::size_t jt, std::size_t kt);
    void LaxWendScheme(REAL UL[], REAL UR[], REAL FL[], REAL FR[], std::size_t iDim, std::size_t it, std::size_t jt, std::size_t kt);
    void exactSolver(REAL uL[], REAL uR[], REAL relaCoord, std::size_t idim, std::size_t i, std::size_t j, std::size_t k);
    void RoeSolver(REAL uL[], REAL uR[], REAL FL[], REAL FR[], REAL HL, REAL HR, std::size_t idim, std::size_t i, std::size_t j, std::size_t k);
    void HllcSolver(REAL uL[], REAL uR[], REAL FL[], REAL FR[], REAL HL, REAL HR, std::size_t idim, std::size_t i, std::size_t j, std::size_t k);
    void HlleSolver(REAL uL[], REAL uR[], REAL FL[], REAL FR[], REAL HL, REAL HR, std::size_t idim, std::size_t i, std::size_t j, std::size_t k);
    void UtoW(); // U - conserved; W - primitive
    void WtoU();
    void rotateIJK(std::vector<Particle *> &ptcls);
    void inteStep1(std::vector<Particle *> &ptcls);
    void RK2InteStep2(std::vector<Particle *> &ptcls);
    void RK3InteStep2(std::vector<Particle *> &ptcls);
    void RK3InteStep3(std::vector<Particle *> &ptcls);

    void coordToGlobalIndex(REAL x, REAL y, REAL z, std::size_t &i, std::size_t &j, std::size_t &k);
    void coordToGlobalIndex(Vec v, std::size_t &i, std::size_t &j, std::size_t &k);
    void coordToGlobalIndex(Vec v, IJK &t);
    void localIndexToGlobal(IJK &local, IJK &global);

    void getPtclInfo(std::vector<Particle *> &ptcls, Gradation &gradation);
    void runOneStep(std::vector<Particle *> &ptcls);
    void calcPtclForce(std::vector<Particle *> &ptcls);
    void penalize(std::vector<Particle *> &ptcls);
    void plot(const char *, std::size_t);
    void checkMomentum(std::vector<Particle *> &ptcl);  

  private:
    void guessPressure(REAL uL[], REAL uR[], REAL &pInit); 
    void evalF(REAL &f, REAL &fd, REAL &p, REAL &dk, REAL &pk, REAL &ck);
    void findPrsVel(REAL uL[], REAL uR[], REAL &p, REAL &u);
    void sampling(REAL uL[], REAL uR[], const REAL pStar, const REAL uStar, const REAL s, REAL &d, REAL &u, REAL &p);
  };
  
} // name space dem
#endif
