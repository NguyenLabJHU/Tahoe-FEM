#ifndef FLUID_H
#define FLUID_H

#include "Parameter.h"
#include "realtypes.h"
#include "Vec.h"
#include "Gradation.h"
#include "Rectangle.h"
#include "Boundary.h"
#include "Particle.h"
#include <cstddef>
#include <valarray>

namespace dem {
  
  class Fluid {
    typedef std::valarray< std::valarray< std::valarray <REAL> > > Array3D;
    typedef std::valarray< std::valarray< std::valarray <std::valarray<REAL> > > > Array4D;
    typedef std::valarray< std::valarray< std::valarray <std::valarray< std::valarray<REAL>  > > > > Array5D;

  private:
    static const REAL Rs  = 287.06; // specific gas constant

    std::size_t gridNx;   // gridNx = total cell centers = parts + two boundary points in x direction
    std::size_t gridNy;   // gridNy = total cell centers = parts + two boundary points in y direction
    std::size_t gridNz;   // gridNz = total cell centers = parts + two boundary points in z direction
    std::size_t ptclGrid; // approximate grids accross particle in each dimension
    REAL gridDx;          // grid size in x direction
    REAL gridDy;          // grid size in y direction
    REAL gridDz;          // grid size in z direction
    REAL x1F, x2F, y1F, y2F, z1F, z2F; // fluid domain

    REAL Cd;           // drag coefficient
    REAL porosity;     // particle porosity as porous media
    REAL Cdi;          // fictitious drag coefficient inside porous media
    REAL velMod;       // velocity correction coefficient in total enthalpy as porous media
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
    Array4D arrayURota; // copy of arrayU in rotation.
    Array4D arrayUPrev; // copy of arrayU for previous step used in Runge-Kutta
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
    // arrayGridCoord[i][j][k][0]: coorgridDx
    // arrayGridCoord[i][j][k][1]: coorgridDy
    // arrayGridCoord[i][j][k][2]: coorgridDz

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
    Array5D arrayGodFluxStep2; 
    Array5D arrayGodFluxStep3; 
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
    Fluid() {}
    
    void initParameter(Rectangle &container, Gradation &gradation);
    void initialize();
    void initialCondition();
    void calcTimeStep();
    void RankineHugoniot();
    void initGhostPoints();
    void soundSpeed();
    void enthalpy();
    void flux(std::size_t, std::vector<Particle *> &ptcls);
    void exactSolver(REAL uL[], REAL uR[], REAL relaCoord, std::size_t idim, std::size_t i, std::size_t j, std::size_t k);
    void RoeSolver(REAL uL[], REAL uR[], REAL FL[], REAL FR[], REAL HL, REAL HR, std::size_t idim, std::size_t i, std::size_t j, std::size_t k);
    void HllcSolver(REAL uL[], REAL uR[], REAL FL[], REAL FR[], REAL HL, REAL HR, std::size_t idim, std::size_t i, std::size_t j, std::size_t k);
    void HlleSolver(REAL uL[], REAL uR[], REAL FL[], REAL FR[], REAL HL, REAL HR, std::size_t idim, std::size_t i, std::size_t j, std::size_t k);
    void UtoW(); // U - conserved; W - primitive
    void WtoU();
    void rotateIJK(std::vector<Particle *> &ptcls);
    void inteStep1(std::vector<Particle *> &ptcls);
    void inteStep2(std::vector<Particle *> &ptcls);
    void inteStep3(std::vector<Particle *> &ptcls);

    void getPtclInfo(std::vector<Particle *> &ptcls);
    void runOneStep(std::vector<Particle *> &ptcls);
    void calcPtclForce(std::vector<Particle *> &ptcls);
    void penalize(std::vector<Particle *> &ptcls);
    void plot(const char *) const;
    void checkMomentum(std::vector<Particle *> &ptcl);   

    void exactGuessPressure(REAL uL[], REAL uR[], REAL &pInit); 
    void exactEvalF(REAL &f, REAL &fd, REAL &p, REAL &dk, REAL &pk, REAL &ck);
    void exactFindPrsVel(REAL uL[], REAL uR[], REAL &p, REAL &u);
    void exactSampling(REAL uL[], REAL uR[], const REAL pStar, const REAL uStar, const REAL s, REAL &d, REAL &u, REAL &p);
  };
  
} // name space dem
#endif
