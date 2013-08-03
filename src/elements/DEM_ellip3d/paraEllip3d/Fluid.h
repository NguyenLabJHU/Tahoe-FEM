#ifndef FLUID_H
#define FLUID_H

#include "Parameter.h"
#include "realtypes.h"
#include "Vec.h"
#include "Gradation.h"
#include "Rectangle.h"
#include "Boundary.h"
#include "Particle.h"
#include <valarray>

namespace dem {
  
  class Fluid {
    typedef std::valarray< std::valarray< std::valarray <REAL> > > Array3D;
    typedef std::valarray< std::valarray< std::valarray <std::valarray<REAL> > > > Array4D;
    typedef std::valarray< std::valarray< std::valarray <std::valarray< std::valarray<REAL>  > > > > Array5D;

  private:
    REAL minX;
    REAL minY;
    REAL minZ;
    REAL maxX;
    REAL maxY;
    REAL maxZ;

    std::size_t nx; // nx = total cellcenters = parts + two boundary points in x direction
    std::size_t ny; // ny = total cellcenters = parts + two boundary points in y direction
    std::size_t nz; // nz = total cellcenters = parts + two boundary points in z direction
    REAL dx;
    REAL dy;
    REAL dz;
    REAL dt;

    REAL CFL;          // Courant Fredics Levy
    REAL gamma;
    bool reflecting;   // 0 - non-reflecting; 1 - reflecting
    REAL rhoL, uL, pL; // known
    REAL rhoR, uR, pR; // only uR is known
    REAL shockSpeed;   // unknown
    REAL z0;           // initial discontinuity plane in Z direction

    std::size_t n_dim, n_var, n_integ, var_den, var_eng, var_prs, var_msk;
    std::size_t var_mom[3], var_vel[3];

    Array4D arrayGridCoord; 
    // fluid compute grid coordinates, 4-dimensional
    // nx, ny, nz, n_dim
    // arrayMesh[i][j][k][0]: coord_x
    // arrayMesh[i][j][k][1]: coord_y
    // arrayMesh[i][j][k][2]: coord_z

    Array4D arrayU;
    Array4D arrayUtmp;
    // 4-dimensional
    // nx, ny, nz, n_var
    // (a) fixed:
    // arrayU[i][j][k][0]: var_den
    // arrayU[i][j][k][1]: var_mom[1]
    // arrayU[i][j][k][2]: var_mom[2]
    // arrayU[i][j][k][3]: var_mom[3]
    // arrayU[i][j][k][4]: var_eng
    // arrayU[i][j][k][5]: var_vel[1]
    // arrayU[i][j][k][6]: var_vel[2]
    // arrayU[i][j][k][7]: var_vel[3]
    // arrayU[i][j][k][8]: var_prs
    // (b) extended:
    // arrayU[i][j][k][9]: var_msk

    Array4D arrayFlux;
    // 4-dimensional
    // nx, ny, nz, n_integ
    // arrayFlux[i][j][k][0]: var_den
    // arrayFlux[i][j][k][1]: var_mom[1]
    // arrayFlux[i][j][k][2]: var_mom[2]
    // arrayFlux[i][j][k][3]: var_mom[3]
    // arrayFlux[i][j][k][4]: var_eng

    Array5D arrayRoeFlux; 
    // 5-dimensional
    // nx-1, ny-1, nz-1, n_integ, n_dim
    // arrayRoeFlux[i][j][k][0]: var_den
    // arrayRoeFlux[i][j][k][1]: var_mom[1]
    // arrayRoeFlux[i][j][k][2]: var_mom[2]
    // arrayRoeFlux[i][j][k][3]: var_mom[3]
    // arrayRoeFlux[i][j][k][4]: var_eng

    Array4D arrayRoeFluxTmp;
    // 4-dimensional
    // nx-1, ny-1, nz-1, n_integ

    Array3D arrayH; 
    // Enthalpy, 3-dimensional
    // nx, ny, nz

    Array3D arraySoundSpeed; 
    // speed of sound, 3-dimensional
    // nx, ny, nz

  public:
    Fluid() {}
    
    void initParameter(Rectangle &container, Gradation &gradation);
    void initialize();
    void initialCondition();
    REAL timeStep();
    void RankineHugoniot();
    void addGhostPoints();
    void soundSpeed();
    void enthalpy();
    void flux();
    void RoeFlux(REAL uL[], REAL uR[], REAL FL[], REAL FR[], REAL HL, REAL HR, std::size_t idim,  std::size_t i, std::size_t j, std::size_t k);
    void UtoW(); // U - integrated; W - primitive
    void WtoU();

    void getParticleInfo(std::vector<Particle *> &ptcls);
    void runOneStep();
    void calcParticleForce(std::vector<Particle *> &ptcls);
    void plot(const char *) const;
    
  };
  
} // name space dem
#endif
