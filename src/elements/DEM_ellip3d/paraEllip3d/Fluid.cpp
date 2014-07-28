#include "Fluid.h"
#include "const.h"
#include <cmath>
#include <algorithm>

namespace dem {

  const REAL Fluid::Rs;

  void Fluid::initParameter(Rectangle &container, Gradation &gradation) {

    std::size_t ptclGrid = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["ptclGrid"]);
    Cd   = dem::Parameter::getSingleton().parameter["Cd"];
    porosity = dem::Parameter::getSingleton().parameter["porosity"];
    RK = dem::Parameter::getSingleton().parameter["RK"];
    CFL = dem::Parameter::getSingleton().parameter["CFL"];
    gamma = dem::Parameter::getSingleton().parameter["airGamma"];
    rhoR = dem::Parameter::getSingleton().parameter["rightDensity"];
    pR   = dem::Parameter::getSingleton().parameter["rightPressure"];
    uR   = dem::Parameter::getSingleton().parameter["rightVelocity"];
    
    leftType = static_cast<int> (dem::Parameter::getSingleton().parameter["leftType"]);    
    REAL x1F = dem::Parameter::getSingleton().parameter["x1F"];
    REAL y1F = dem::Parameter::getSingleton().parameter["y1F"];
    REAL z1F = dem::Parameter::getSingleton().parameter["z1F"];
    REAL x2F = dem::Parameter::getSingleton().parameter["x2F"];
    REAL y2F = dem::Parameter::getSingleton().parameter["y2F"];
    REAL z2F = dem::Parameter::getSingleton().parameter["z2F"];
    arrayBC[0] = dem::Parameter::getSingleton().parameter["x1Reflecting"];
    arrayBC[1] = dem::Parameter::getSingleton().parameter["x2Reflecting"];
    arrayBC[2] = dem::Parameter::getSingleton().parameter["y1Reflecting"];
    arrayBC[3] = dem::Parameter::getSingleton().parameter["y2Reflecting"];
    arrayBC[4] = dem::Parameter::getSingleton().parameter["z1Reflecting"];
    arrayBC[5] = dem::Parameter::getSingleton().parameter["z2Reflecting"];

    if (leftType == 1) {
      z2L = dem::Parameter::getSingleton().parameter["z2L"];
      Mach= dem::Parameter::getSingleton().parameter["shockMach"];
    }      
    else if (leftType == 2) {
      z2L = dem::Parameter::getSingleton().parameter["z2L"];
      rhoL= dem::Parameter::getSingleton().parameter["leftDensity"];
      pL  = dem::Parameter::getSingleton().parameter["leftPressure"];
      uL  = dem::Parameter::getSingleton().parameter["leftVelocity"];
    }
    else if (leftType == 3) {
      x1L = dem::Parameter::getSingleton().parameter["x1L"];
      x2L = dem::Parameter::getSingleton().parameter["x2L"];
      y1L = dem::Parameter::getSingleton().parameter["y1L"];
      y2L = dem::Parameter::getSingleton().parameter["y2L"];
      z1L = dem::Parameter::getSingleton().parameter["z1L"];
      z2L = dem::Parameter::getSingleton().parameter["z2L"];
      rhoL= dem::Parameter::getSingleton().parameter["leftDensity"];
      pL  = dem::Parameter::getSingleton().parameter["leftPressure"];
      uL  = dem::Parameter::getSingleton().parameter["leftVelocity"];
    } 
    else if (leftType == 4) {
      x0L = dem::Parameter::getSingleton().parameter["x0L"];
      y0L = dem::Parameter::getSingleton().parameter["y0L"];
      z0L = dem::Parameter::getSingleton().parameter["z0L"];
      r0L = dem::Parameter::getSingleton().parameter["r0L"];
      rhoL= dem::Parameter::getSingleton().parameter["leftDensity"];
      pL  = dem::Parameter::getSingleton().parameter["leftPressure"];
      uL  = dem::Parameter::getSingleton().parameter["leftVelocity"];
    }
    else if (leftType == 5) {
      x1L = dem::Parameter::getSingleton().parameter["x1L"];
      x2L = dem::Parameter::getSingleton().parameter["x2L"];
      y1L = dem::Parameter::getSingleton().parameter["y1L"];
      y2L = dem::Parameter::getSingleton().parameter["y2L"];
      z1L = dem::Parameter::getSingleton().parameter["z1L"];
      z2L = dem::Parameter::getSingleton().parameter["z2L"];
      rhoL= dem::Parameter::getSingleton().parameter["leftDensity"];
      pL  = dem::Parameter::getSingleton().parameter["leftPressure"];
      uL  = dem::Parameter::getSingleton().parameter["leftVelocity"];
      rhoBL= dem::Parameter::getSingleton().parameter["belowLeftDensity"];
      pBL  = dem::Parameter::getSingleton().parameter["belowLeftPressure"];
      uBL  = dem::Parameter::getSingleton().parameter["belowLeftVelocity"];
    }    

    REAL minR = gradation.getPtclMinRadius();
    dx = (minR * 2) / ptclGrid;
    dy = dx;
    dz = dx;
    nx = static_cast<std::size_t> (ceil((x2F - x1F) / dx));
    ny = static_cast<std::size_t> (ceil((y2F - y1F) / dy));
    nz = static_cast<std::size_t> (ceil((z2F - z1F) / dz));

    dx = (x2F - x1F) / nx;
    dy = (y2F - y1F) / ny;
    dz = (z2F - z1F) / nz;

    nx += 2;
    ny += 2;
    nz += 2;

    // fixed
    nDim = 3;
    nVar = 0; 
    nInteg = 0;

    varDen = nVar++; nInteg++;
    varMom[0] = nVar++; nInteg++;
    varMom[1] = nVar++; nInteg++;
    varMom[2] = nVar++; nInteg++;
    varEng = nVar++; nInteg++;

    varVel[0] = nVar++;
    varVel[1] = nVar++;
    varVel[2] = nVar++;
    varPrs = nVar++; 

    // extended
    varMsk = nVar++;

    // print
    debugInf << std::setw(OWID) << "ptclGrid" << std::setw(OWID) << ptclGrid << std::endl;
    debugInf << std::setw(OWID) << "Cd" << std::setw(OWID) << Cd << std::endl;
    debugInf << std::setw(OWID) << "porosity" << std::setw(OWID) << porosity << std::endl;
    debugInf << std::setw(OWID) << "Runge-Kutta" << std::setw(OWID) << (int) RK << std::endl;
    debugInf << std::setw(OWID) << "CFL" << std::setw(OWID) << CFL << std::endl;
    debugInf << std::setw(OWID) << "gamma" << std::setw(OWID) << gamma << std::endl;
    debugInf << std::setw(OWID) << "rhoR" << std::setw(OWID) << rhoR << std::endl;
    debugInf << std::setw(OWID) << "pR" << std::setw(OWID) << pR << std::endl;
    debugInf << std::setw(OWID) << "uR" << std::setw(OWID) << uR << std::endl;

    debugInf << std::setw(OWID) << "gridSize" << std::setw(OWID) << dx << std::endl;
    debugInf << std::setw(OWID) << "gridX" << std::setw(OWID) << nx << std::endl;
    debugInf << std::setw(OWID) << "gridY" << std::setw(OWID) << ny << std::endl;
    debugInf << std::setw(OWID) << "gridZ" << std::setw(OWID) << nz << std::endl;

    debugInf << std::setw(OWID) << "leftType" << std::setw(OWID) << leftType << std::endl;
    debugInf << std::setw(OWID) << "x1F" << std::setw(OWID) << x1F << std::endl;
    debugInf << std::setw(OWID) << "y1F" << std::setw(OWID) << y1F << std::endl;
    debugInf << std::setw(OWID) << "z1F" << std::setw(OWID) << z1F << std::endl;
    debugInf << std::setw(OWID) << "x2F" << std::setw(OWID) << x2F << std::endl;
    debugInf << std::setw(OWID) << "y2F" << std::setw(OWID) << y2F << std::endl;
    debugInf << std::setw(OWID) << "z2F" << std::setw(OWID) << z2F << std::endl;
    debugInf << std::setw(OWID) << "x1Rflecting" << std::setw(OWID) << (int) arrayBC[0] << std::endl;
    debugInf << std::setw(OWID) << "x2Rflecting" << std::setw(OWID) << (int) arrayBC[1] << std::endl;
    debugInf << std::setw(OWID) << "y1Rflecting" << std::setw(OWID) << (int) arrayBC[2] << std::endl;
    debugInf << std::setw(OWID) << "y2Rflecting" << std::setw(OWID) << (int) arrayBC[3] << std::endl;
    debugInf << std::setw(OWID) << "z1Rflecting" << std::setw(OWID) << (int) arrayBC[4] << std::endl;
    debugInf << std::setw(OWID) << "z2Rflecting" << std::setw(OWID) << (int) arrayBC[5] << std::endl;

    if (leftType == 1) {
      debugInf << std::setw(OWID) << "z2L" << std::setw(OWID) << z2L << std::endl;
      debugInf << std::setw(OWID) << "shockMach" << std::setw(OWID) << Mach << std::endl;
    }
    else if (leftType == 2) {
      debugInf << std::setw(OWID) << "z2L" << std::setw(OWID) << z2L << std::endl;
      debugInf << std::setw(OWID) << "rhoL" << std::setw(OWID) << rhoL << std::endl;
      debugInf << std::setw(OWID) << "pL" << std::setw(OWID) << pL << std::endl; 
      debugInf << std::setw(OWID) << "uL" << std::setw(OWID) << uL << std::endl;      
    } else if (leftType == 3) {
      debugInf << std::setw(OWID) << "x1L" << std::setw(OWID) << x1L << std::endl;
      debugInf << std::setw(OWID) << "x2L" << std::setw(OWID) << x2L << std::endl;
      debugInf << std::setw(OWID) << "y1L" << std::setw(OWID) << y1L << std::endl;
      debugInf << std::setw(OWID) << "y2L" << std::setw(OWID) << y2L << std::endl;
      debugInf << std::setw(OWID) << "z1L" << std::setw(OWID) << z1L << std::endl;
      debugInf << std::setw(OWID) << "z2L" << std::setw(OWID) << z2L << std::endl;
      debugInf << std::setw(OWID) << "rhoL" << std::setw(OWID) << rhoL << std::endl;
      debugInf << std::setw(OWID) << "pL" << std::setw(OWID) << pL << std::endl;  
      debugInf << std::setw(OWID) << "uL" << std::setw(OWID) << uL << std::endl;     
    } else if(leftType == 4) {
      debugInf << std::setw(OWID) << "x0L" << std::setw(OWID) << x0L << std::endl;
      debugInf << std::setw(OWID) << "y0L" << std::setw(OWID) << y0L << std::endl;
      debugInf << std::setw(OWID) << "z0L" << std::setw(OWID) << z0L << std::endl;
      debugInf << std::setw(OWID) << "r0L" << std::setw(OWID) << r0L << std::endl;
      debugInf << std::setw(OWID) << "rhoL" << std::setw(OWID) << rhoL << std::endl;
      debugInf << std::setw(OWID) << "pL" << std::setw(OWID) << pL << std::endl;    
      debugInf << std::setw(OWID) << "uL" << std::setw(OWID) << uL << std::endl;
    } else if (leftType == 5) {
      debugInf << std::setw(OWID) << "x1L" << std::setw(OWID) << x1L << std::endl;
      debugInf << std::setw(OWID) << "x2L" << std::setw(OWID) << x2L << std::endl;
      debugInf << std::setw(OWID) << "y1L" << std::setw(OWID) << y1L << std::endl;
      debugInf << std::setw(OWID) << "y2L" << std::setw(OWID) << y2L << std::endl;
      debugInf << std::setw(OWID) << "z1L" << std::setw(OWID) << z1L << std::endl;
      debugInf << std::setw(OWID) << "z2L" << std::setw(OWID) << z2L << std::endl;
      debugInf << std::setw(OWID) << "rhoL" << std::setw(OWID) << rhoL << std::endl;
      debugInf << std::setw(OWID) << "pL" << std::setw(OWID) << pL << std::endl;  
      debugInf << std::setw(OWID) << "uL" << std::setw(OWID) << uL << std::endl; 
      debugInf << std::setw(OWID) << "rhoBL" << std::setw(OWID) << rhoBL << std::endl;
      debugInf << std::setw(OWID) << "pBL" << std::setw(OWID) << pBL << std::endl;  
      debugInf << std::setw(OWID) << "uBL" << std::setw(OWID) << uBL << std::endl; 
    }

    /*
    debugInf << "nVar " << nVar << std::endl;
    debugInf << "nInteg " << nInteg << std::endl;
    debugInf << "varDen " << varDen  << std::endl;    
    debugInf << "varMom[0] " << varMom[0] << std::endl;    
    debugInf << "varMom[1] " << varMom[1] << std::endl;
    debugInf << "varMom[2] " << varMom[2] << std::endl;    
    debugInf << "varEng " << varEng  << std::endl;    
    debugInf << "varVel[0] " << varVel[0] << std::endl;    
    debugInf << "varVel[1] " << varVel[1] << std::endl;    
    debugInf << "varVel[2] " << varVel[2] << std::endl;    
    debugInf << "varPrs " << varPrs  << std::endl;    
    debugInf << "varMsk " << varMsk  << std::endl;   
    */ 

    // nx, ny, nz, nDim
    arrayGridCoord.resize(nx);
    for (std::size_t i = 0; i < arrayGridCoord.size(); ++i) {
      arrayGridCoord[i].resize(ny);
      for (std::size_t j = 0; j < arrayGridCoord[i].size(); ++j) {
	arrayGridCoord[i][j].resize(nz);
	for (std::size_t k = 0; k < arrayGridCoord[i][j].size(); ++k) 
	  arrayGridCoord[i][j][k].resize(nDim);
      }
    }

    // coordinates
    for (std::size_t i = 0; i < arrayGridCoord.size(); ++i)
      for (std::size_t j = 0; j < arrayGridCoord[i].size(); ++j)
	for (std::size_t k = 0; k < arrayGridCoord[i][j].size(); ++k) {
	  arrayGridCoord[i][j][k][0] = (x1F - dx/2) + i * dx;
	  arrayGridCoord[i][j][k][1] = (y1F - dy/2) + j * dy;
	  arrayGridCoord[i][j][k][2] = (z1F - dz/2) + k * dz;
	}

    // nx, ny, nz, nDim
    arrayPenalForce.resize(nx);
    for (std::size_t i = 0; i < arrayPenalForce.size(); ++i) {
      arrayPenalForce[i].resize(ny);
      for (std::size_t j = 0; j < arrayPenalForce[i].size(); ++j) {
	arrayPenalForce[i][j].resize(nz);
	for (std::size_t k = 0; k < arrayPenalForce[i][j].size(); ++k) 
	  arrayPenalForce[i][j][k].resize(nDim);
      }
    }

    // nx, ny, nz, nDim
    arrayPressureForce.resize(nx);
    for (std::size_t i = 0; i < arrayPressureForce.size(); ++i) {
      arrayPressureForce[i].resize(ny);
      for (std::size_t j = 0; j < arrayPressureForce[i].size(); ++j) {
	arrayPressureForce[i][j].resize(nz);
	for (std::size_t k = 0; k < arrayPressureForce[i][j].size(); ++k) 
	  arrayPressureForce[i][j][k].resize(nDim);
      }
    }

    // nx, ny, nz, nVar
    arrayU.resize(nx);
    for (std::size_t i = 0; i < arrayU.size(); ++i) {
      arrayU[i].resize(ny);
      for (std::size_t j = 0; j < arrayU[i].size(); ++j) {
	arrayU[i][j].resize(nz);
	for (std::size_t k = 0; k < arrayU[i][j].size(); ++k) 
	  arrayU[i][j][k].resize(nVar);
      }
    }

    // nx, ny, nz, nVar
    arrayUtmp.resize(nx);
    for (std::size_t i = 0; i < arrayUtmp.size(); ++i) {
      arrayUtmp[i].resize(ny);
      for (std::size_t j = 0; j < arrayUtmp[i].size(); ++j) {
	arrayUtmp[i][j].resize(nz);
	for (std::size_t k = 0; k < arrayUtmp[i][j].size(); ++k) 
	  arrayUtmp[i][j][k].resize(nVar);
      }
    }

    // nx, ny, nz, nInteg
    arrayFlux.resize(nx);
    for (std::size_t i = 0; i < arrayFlux.size(); ++i) {
      arrayFlux[i].resize(ny);
      for (std::size_t j = 0; j < arrayFlux[i].size(); ++j) {
	arrayFlux[i][j].resize(nz);
	for (std::size_t k = 0; k < arrayFlux[i][j].size(); ++k) 
	  arrayFlux[i][j][k].resize(nInteg);
      }
    }

    // nx-1, ny-1, nz-1, nInteg, nDim
    arrayRoeFlux.resize(nx-1);
    for (std::size_t i = 0; i < arrayRoeFlux.size(); ++i) {
      arrayRoeFlux[i].resize(ny-1);
      for (std::size_t j = 0; j < arrayRoeFlux[i].size(); ++j) {
	arrayRoeFlux[i][j].resize(nz-1);
	for (std::size_t k = 0; k < arrayRoeFlux[i][j].size(); ++k) {
	  arrayRoeFlux[i][j][k].resize(nInteg);
	  for (std::size_t m = 0; m < arrayRoeFlux[i][j][k].size(); ++m)
	    arrayRoeFlux[i][j][k][m].resize(nDim);
	}
      }
    }

    if (RK >= 1) {
      // nx-1, ny-1, nz-1, nInteg, nDim
      arrayRoeFluxStep2.resize(nx-1);
      for (std::size_t i = 0; i < arrayRoeFluxStep2.size(); ++i) {
	arrayRoeFluxStep2[i].resize(ny-1);
	for (std::size_t j = 0; j < arrayRoeFluxStep2[i].size(); ++j) {
	  arrayRoeFluxStep2[i][j].resize(nz-1);
	  for (std::size_t k = 0; k < arrayRoeFluxStep2[i][j].size(); ++k) {
	    arrayRoeFluxStep2[i][j][k].resize(nInteg);
	    for (std::size_t m = 0; m < arrayRoeFluxStep2[i][j][k].size(); ++m)
	      arrayRoeFluxStep2[i][j][k][m].resize(nDim);
	  }
	}
      }
    }

    if (RK == 2) {
      // nx-1, ny-1, nz-1, nInteg, nDim
      arrayRoeFluxStep3.resize(nx-1);
      for (std::size_t i = 0; i < arrayRoeFluxStep3.size(); ++i) {
	arrayRoeFluxStep3[i].resize(ny-1);
	for (std::size_t j = 0; j < arrayRoeFluxStep3[i].size(); ++j) {
	  arrayRoeFluxStep3[i][j].resize(nz-1);
	  for (std::size_t k = 0; k < arrayRoeFluxStep3[i][j].size(); ++k) {
	    arrayRoeFluxStep3[i][j][k].resize(nInteg);
	    for (std::size_t m = 0; m < arrayRoeFluxStep3[i][j][k].size(); ++m)
	      arrayRoeFluxStep3[i][j][k][m].resize(nDim);
	  }
	}
      }
    }

    // nx-1, ny-1, nz-1, nInteg
    arrayRoeFluxTmp.resize(nx-1);
    for (std::size_t i = 0; i < arrayRoeFluxTmp.size(); ++i) {
      arrayRoeFluxTmp[i].resize(ny-1);
      for (std::size_t j = 0; j < arrayRoeFluxTmp[i].size(); ++j) {
	arrayRoeFluxTmp[i][j].resize(nz-1);
	for (std::size_t k = 0; k < arrayRoeFluxTmp[i][j].size(); ++k) {
	  arrayRoeFluxTmp[i][j][k].resize(nInteg);
	}
      }
    }
    
    // nx, ny, nz
    arrayH.resize(nx);
    for (std::size_t i = 0; i < arrayH.size(); ++i) {
      arrayH[i].resize(ny);
      for (std::size_t j = 0; j < arrayH[i].size(); ++j)
	arrayH[i][j].resize(nz);
    }

    // nx, ny, nz
    arraySoundSpeed.resize(nx);
    for (std::size_t i = 0; i < arraySoundSpeed.size(); ++i) {
      arraySoundSpeed[i].resize(ny);
      for (std::size_t j = 0; j < arraySoundSpeed[i].size(); ++j)
	arraySoundSpeed[i][j].resize(nz);
    }

  }

  void Fluid::initialize() {
    RankineHugoniot();
    initialCondition(); 
    soundSpeed(); // for printing Mach number
  }

  void Fluid::runOneStep(std::vector<Particle *> &ptcls) {
    inteStep1(ptcls);

    if (RK >= 1) {
      arrayRoeFluxStep2 = arrayRoeFlux;
      inteStep2(ptcls);
      if (RK == 2) {
	arrayRoeFluxStep3 = arrayRoeFlux;    
	inteStep3(ptcls);
      }
    }
  }

  void Fluid::inteStep1(std::vector<Particle *> &ptcls) { 
    addGhostPoints();
    soundSpeed();
    timeStep = std::min(timeStep, calcTimeStep());
    debugInf << std::setw(OWID) << iteration 
	     << std::setw(OWID) << timeStep 
	     << std::setw(OWID) << timeAccrued << std::endl;
    enthalpy();
    rotateIJK(ptcls);

    // update conserved variables at the next time step
    for (std::size_t i = 1; i < nx - 1 ; ++i)
      for (std::size_t j = 1; j < ny - 1; ++j)
	for (std::size_t k = 1; k < nz - 1; ++k) {
	  for (std::size_t m = 0; m < nInteg; ++m) {
	    arrayU[i][j][k][m] -= (   timeStep / dx * (arrayRoeFlux[i][j][k][m][0] - arrayRoeFlux[i-1][j][k][m][0])
				    + timeStep / dy * (arrayRoeFlux[i][j][k][m][1] - arrayRoeFlux[i][j-1][k][m][1])
				    + timeStep / dz * (arrayRoeFlux[i][j][k][m][2] - arrayRoeFlux[i][j][k-1][m][2]) );
	    //if (arrayU[i][j][k][varMsk] == 1) arrayU[i][j][k][m] /= porosity;
	  }
	}

    // calculate primitive after finding conserved variables
    UtoW(); 
  }
  
  void Fluid::inteStep2(std::vector<Particle *> &ptcls) { 
    addGhostPoints();
    soundSpeed();
    enthalpy();
    rotateIJK(ptcls);

    // update conserved variables at the next time step
    for (std::size_t i = 1; i < nx - 1 ; ++i)
      for (std::size_t j = 1; j < ny - 1; ++j)
	for (std::size_t k = 1; k < nz - 1; ++k) {
	  for (std::size_t m = 0; m < nInteg; ++m) {
	    arrayU[i][j][k][m] -= (   timeStep / (2*RK*dx) * (arrayRoeFlux[i][j][k][m][0] - arrayRoeFlux[i-1][j][k][m][0] 
						        + (arrayRoeFluxStep2[i][j][k][m][0] - arrayRoeFluxStep2[i-1][j][k][m][0]) )
				    + timeStep / (2*RK*dy) * (arrayRoeFlux[i][j][k][m][1] - arrayRoeFlux[i][j-1][k][m][1] 
						        + (arrayRoeFluxStep2[i][j][k][m][1] - arrayRoeFluxStep2[i][j-1][k][m][1]) )
				    + timeStep / (2*RK*dz) * (arrayRoeFlux[i][j][k][m][2] - arrayRoeFlux[i][j][k-1][m][2] 
						        + (arrayRoeFluxStep2[i][j][k][m][2] - arrayRoeFluxStep2[i][j][k-1][m][2])) );
	    //if (arrayU[i][j][k][varMsk] == 1) arrayU[i][j][k][m] /= porosity;
	  }
	}

    // calculate primitive after finding conserved variables
    UtoW(); 
  }

  void Fluid::inteStep3(std::vector<Particle *> &ptcls) { 
    addGhostPoints();
    soundSpeed();
    enthalpy();
    rotateIJK(ptcls);

    // update conserved variables at the next time step
    for (std::size_t i = 1; i < nx - 1 ; ++i)
      for (std::size_t j = 1; j < ny - 1; ++j)
	for (std::size_t k = 1; k < nz - 1; ++k) {
	  for (std::size_t m = 0; m < nInteg; ++m) {
	    arrayU[i][j][k][m] -= (   timeStep / (6*dx) * (arrayRoeFlux[i][j][k][m][0] - arrayRoeFlux[i-1][j][k][m][0]
						         +(arrayRoeFluxStep2[i][j][k][m][0] - arrayRoeFluxStep2[i-1][j][k][m][0])
						      + 4*(arrayRoeFluxStep3[i][j][k][m][0] - arrayRoeFluxStep3[i-1][j][k][m][0]))
				      
				    + timeStep / (6*dy) * (arrayRoeFlux[i][j][k][m][1] - arrayRoeFlux[i][j-1][k][m][1]
						        + (arrayRoeFluxStep2[i][j][k][m][1] - arrayRoeFluxStep2[i][j-1][k][m][1])
						      + 4*(arrayRoeFluxStep3[i][j][k][m][1] - arrayRoeFluxStep3[i][j-1][k][m][1]))

				    + timeStep / (6*dz) * (arrayRoeFlux[i][j][k][m][2] - arrayRoeFlux[i][j][k-1][m][2]
						        + (arrayRoeFluxStep2[i][j][k][m][2] - arrayRoeFluxStep2[i][j][k-1][m][2])
						     + 4*( arrayRoeFluxStep3[i][j][k][m][2] - arrayRoeFluxStep3[i][j][k-1][m][2])) );
	    //if (arrayU[i][j][k][varMsk] == 1) arrayU[i][j][k][m] /= porosity;
	  }
	}

    // calculate primitive after finding conserved variables
    UtoW(); 
  }

  void Fluid::rotateIJK(std::vector<Particle *> &ptcls) {
    std::size_t id[3][3] = {{0,1,2},{1,0,2},{2,1,0}};

    // for x, y, z directions
    for (std::size_t idim = 0; idim < nDim; ++idim) {
      arrayUtmp = arrayU; // must have same rank and extent

      // switch components
      for (std::size_t i = 0; i < nx; ++i)
	for (std::size_t j = 0; j < ny; ++j)
	  for (std::size_t k = 0; k < nz; ++k) 
	    for (std::size_t jdim = 0; jdim < nDim; ++jdim) {
	      arrayUtmp[i][j][k][  varMom[jdim]  ] = arrayU[i][j][k][  varMom[id[idim][jdim]]  ];
	      arrayUtmp[i][j][k][  varVel[jdim]  ] = arrayU[i][j][k][  varVel[id[idim][jdim]]  ];
	    }

      flux(idim, ptcls); // variables defined at cell centers

      // for local Riemann problem
      for (std::size_t i = 0; i < nx - 1; ++i) { // variables defined at cell faces
	for (std::size_t j = 0; j < ny - 1; ++j) {
	  for (std::size_t k = 0; k < nz -1; ++k) {
	    std::size_t IL[3] = {i, j, k};
	    std::size_t IR[3] = {i, j, k};
	    IR[idim] += 1;
	    REAL uL[9], uR[9], FL[5], FR[5], HL, HR; // local variable only
	    HL = arrayH[IL[0]] [IL[1]] [IL[2]];
	    HR = arrayH[IR[0]] [IR[1]] [IR[2]];
	    for (std::size_t m = 0; m < nVar; ++m) {
	      uL[m] = arrayUtmp[IL[0]] [IL[1]] [IL[2]] [m];
	      uR[m] = arrayUtmp[IR[0]] [IR[1]] [IR[2]] [m];
	    }	
	    for (std::size_t m = 0; m < nInteg; ++m) {
	      FL[m] = arrayFlux[IL[0]] [IL[1]] [IL[2]] [m];
	      FR[m] = arrayFlux[IR[0]] [IR[1]] [IR[2]] [m];
	    }    
	    RoeFlux(uL, uR, FL, FR, HL, HR, idim, i, j, k);
	  }
	}
      }

      for (std::size_t i = 0; i < nx -1; ++i)
	for (std::size_t j = 0; j < ny -1; ++j)
	  for (std::size_t k = 0; k < nz -1; ++k)
	    for (std::size_t m = 0; m < nInteg; ++m)
	      arrayRoeFluxTmp[i][j][k][m] = arrayRoeFlux[i][j][k][m][idim];

      // switch components back for consistency with u
      for (std::size_t i = 0; i < nx - 1; ++i)
	for (std::size_t j = 0; j < ny - 1; ++j)
	  for (std::size_t k = 0; k < nz -1; ++k)
	    for (std::size_t m = 0; m < nDim; ++m)
	      arrayRoeFlux[i][j][k][varMom[m]][idim] = arrayRoeFluxTmp[i][j][k][ varMom[id[idim][m]] ];

    } // end of for x, y, z directions

  }

  void Fluid::penalize(std::vector<Particle *> &ptcls) {
    // 2nd implementation: higher efficiency
    // for cells that are enclosed by particle volumes
    for (std::vector<Particle *>::const_iterator it = ptcls.begin(); it != ptcls.end(); ++it) {
      std::vector< std::vector<REAL> > fluidGrid = (*it)->getFluidGrid();
      for (std::size_t iter = 0; iter < fluidGrid.size(); ++iter) {

	std::size_t i = static_cast<std::size_t> (fluidGrid[iter][0]);
	std::size_t j = static_cast<std::size_t> (fluidGrid[iter][1]);
	std::size_t k = static_cast<std::size_t> (fluidGrid[iter][2]);

	// calculate momentum and density quotient before modification
	bool inGrid = (i > 0 && i < nx-1 && j > 0 && j < ny-1 && k > 0 && k < nz-1 );
	REAL momQuot[3], denQuot[3], u0[3];
	if (inGrid) {
	  momQuot[0] = (arrayU[i+1][j][k][varMom[0]] - arrayU[i-1][j][k][varMom[0]]) / (2*dx);
	  momQuot[1] = (arrayU[i][j+1][k][varMom[1]] - arrayU[i][j-1][k][varMom[1]]) / (2*dy);
	  momQuot[2] = (arrayU[i][j][k+1][varMom[2]] - arrayU[i][j][k-1][varMom[2]]) / (2*dz);

	  denQuot[0] = (arrayU[i+1][j][k][varDen] - arrayU[i-1][j][k][varDen]) / (2*dx);
	  denQuot[1] = (arrayU[i][j+1][k][varDen] - arrayU[i][j-1][k][varDen]) / (2*dy);
	  denQuot[2] = (arrayU[i][j][k+1][varDen] - arrayU[i][j][k-1][varDen]) / (2*dz);

	  REAL coordX = arrayGridCoord[i][j][k][0];
	  REAL coordY = arrayGridCoord[i][j][k][1];
	  REAL coordZ = arrayGridCoord[i][j][k][2];
	  Vec dist = Vec(coordX, coordY, coordZ) - (*it)->getCurrPos();
	  Vec omgar = (*it)->getCurrOmga() % dist; // w X r = omga % dist, where % is overloaded as cross product
	  u0[0] = (*it)->getCurrVeloc().getX() + omgar.getX(); 
	  u0[1] = (*it)->getCurrVeloc().getY() + omgar.getY(); 
	  u0[2] = (*it)->getCurrVeloc().getZ() + omgar.getZ();
	}

	// 1. momentum penalization
	for (std::size_t m = 0; m < nDim; ++m) {
	  // momentum penalization
	  arrayU[i][j][k][varMom[m]] -= arrayU[i][j][k][varMsk] * arrayPenalForce[i][j][k][m] * timeStep;
	  // influence of momentum penalization on energy
	  arrayU[i][j][k][varEng]    -= arrayU[i][j][k][varMsk] * arrayPenalForce[i][j][k][m] * arrayU[i][j][k][varVel[m]] * timeStep;
	}

	// 2. mass penalization
	// mass penalization and particle velocity term are incorporated by changes in flux()
	// influence of mass penalization (1-1.0/porosity)*momQuot[m] and particle velocity term u0[m]/porosity*denQuot[m] on energy
	if (inGrid) {
	  for (std::size_t m = 0; m < nDim; ++m)
	    arrayU[i][j][k][varEng]  += arrayU[i][j][k][varMsk] * (-0.5*pow(arrayU[i][j][k][varVel[m]],2))
                                         * ( (1-1.0/porosity)*momQuot[m] + u0[m]/porosity*denQuot[m] ) * timeStep;
	}
      }
    }

    /*
    // 1st implementation: lower efficiency
    for (std::size_t i = 0; i < nx; ++i)
      for (std::size_t j = 0; j < ny; ++j)
	for (std::size_t k = 0; k < nz; ++k) {

	  // calculate momentum quotient before modification
	  bool inGrid = (i > 0 && i < nx-1 && j > 0 && j < ny-1 && k > 0 && k < nz-1 );
	  REAL momQuot[3];
	  if (inGrid) {
	    momQuot[0] = (arrayU[i+1][j][k][varMom[0]] - arrayU[i-1][j][k][varMom[0]]) / (2*dx);
	    momQuot[1] = (arrayU[i][j+1][k][varMom[1]] - arrayU[i][j-1][k][varMom[1]]) / (2*dy);
	    momQuot[2] = (arrayU[i][j][k+1][varMom[2]] - arrayU[i][j][k-1][varMom[2]]) / (2*dz);
	  }

	  // penalization
	  for (std::size_t m = 0; m < nDim; ++m) {
	    // momentum penalization
	    arrayU[i][j][k][varMom[m]] -= arrayU[i][j][k][varMsk] * arrayPenalForce[i][j][k][m] * timeStep;
	    // energy penalization
	    arrayU[i][j][k][varEng]    -= arrayU[i][j][k][varMsk] * arrayPenalForce[i][j][k][m] * arrayU[i][j][k][varVel[m]] * timeStep;
	  }

	  // influence of mass penalization on energy
	  if (inGrid) {
	    for (std::size_t m = 0; m < nDim; ++m)
	      arrayU[i][j][k][varEng]  += arrayU[i][j][k][varMsk] * (0.5*pow(arrayU[i][j][k][varVel[m]],2)*(1.0/porosity-1)) * momQuot[m] * timeStep;
	  }
	}
    */
  }
  
  void Fluid::addGhostPoints() {
    // non-reflecting BCs
    for (std::size_t j = 1; j < ny - 1; ++j)
      for (std::size_t k = 1; k < nz - 1; ++k)
	for (std::size_t m = 0; m < nVar; ++m) {
	  arrayU[0][j][k][m]    = arrayU[1][j][k][m]; 
	  arrayU[nx-1][j][k][m] = arrayU[nx-2][j][k][m]; 
	}

    for (std::size_t i = 1; i < nx - 1; ++i)
      for (std::size_t k = 1; k < nz -1; ++k)
	for (std::size_t m = 0; m < nVar; ++m) {
	  arrayU[i][0][k][m]    = arrayU[i][1][k][m]; 
	  arrayU[i][ny-1][k][m] = arrayU[i][ny-2][k][m]; 
	}

    for (std::size_t i = 1; i < nx - 1; ++i)
      for (std::size_t j = 1; j < ny - 1; ++j)
	for (std::size_t m = 0; m < nVar; ++m) {
	  arrayU[i][j][0][m]    = arrayU[i][j][1][m]; 
	  arrayU[i][j][nz-1][m] = arrayU[i][j][nz-2][m]; 
	}

    // reflecting BCs
    bool reflecting = false;
    for (std::size_t it = 0; it < 6; ++it) {
      if (arrayBC[it] > 0) {
	reflecting = true;
	break;
      }
    }

    if (reflecting) {
      for (std::size_t j = 1; j < ny - 1; ++j)
	for (std::size_t k = 1; k < nz - 1; ++k)
	  for (std::size_t m = 0; m < 1; ++m) { // x-direction
	    arrayU[0][j][k][varMom[m]]    *= (1-2*arrayBC[0]); 
	    arrayU[nx-1][j][k][varMom[m]] *= (1-2*arrayBC[1]); 
	    arrayU[0][j][k][varVel[m]]    *= (1-2*arrayBC[0]); 
	    arrayU[nx-1][j][k][varVel[m]] *= (1-2*arrayBC[1]); 
	  }

      for (std::size_t i = 1; i < nx - 1; ++i)
	for (std::size_t k = 1; k < nz - 1; ++k)
	  for (std::size_t m = 1; m < 2; ++m) { // y-direction
	    arrayU[i][0][k][varMom[m]]    *= (1-2*arrayBC[2]); 
	    arrayU[i][ny-1][k][varMom[m]] *= (1-2*arrayBC[3]);
	    arrayU[i][0][k][varVel[m]]    *= (1-2*arrayBC[2]); 
	    arrayU[i][ny-1][k][varVel[m]] *= (1-2*arrayBC[3]);  
	  }

      for (std::size_t i = 1; i < nx - 1; ++i)
	for (std::size_t j = 1; j < ny - 1; ++j)
	  for (std::size_t m = 2; m < 3; ++m) { // z-direction
	    arrayU[i][j][0][varMom[m]]    *= (1-2*arrayBC[4]); 
	    arrayU[i][j][nz-1][varMom[m]] *= (1-2*arrayBC[5]); 
	    arrayU[i][j][0][varVel[m]]    *= (1-2*arrayBC[4]); 
	    arrayU[i][j][nz-1][varVel[m]] *= (1-2*arrayBC[5]); 
	  }
    }
    
  }
  
  REAL Fluid::calcTimeStep() {
    std::valarray<REAL> allGrid(nx * ny * nz);
    for (std::size_t i = 0; i < nx ; ++i)
      for (std::size_t j = 0; j < ny; ++j)
	for (std::size_t k = 0; k < nz; ++k)
	  allGrid[i + j * nx + k * nx * ny] = fabs(arrayU[i][j][k][varVel[0]]) + arraySoundSpeed[i][j][k];
    REAL sx = allGrid.max();

    for (std::size_t i = 0; i < nx ; ++i)
      for (std::size_t j = 0; j < ny; ++j)
	for (std::size_t k = 0; k < nz; ++k)
	  allGrid[i + j * nx + k * nx * ny] = fabs(arrayU[i][j][k][varVel[1]]) + arraySoundSpeed[i][j][k];
    REAL sy = allGrid.max();

    for (std::size_t i = 0; i < nx ; ++i)
      for (std::size_t j = 0; j < ny; ++j)
	for (std::size_t k = 0; k < nz; ++k)
	  allGrid[i + j * nx + k * nx * ny] = fabs(arrayU[i][j][k][varVel[2]]) + arraySoundSpeed[i][j][k];
    REAL sz = allGrid.max();

    std::valarray<REAL> dtMin(3);
    dtMin[0] = dx/sx;
    dtMin[1] = dy/sy;
    dtMin[2] = dz/sz;
    
    return CFL * dtMin.min();
  }

  void Fluid::soundSpeed() {
    for (std::size_t i = 0; i < nx ; ++i)
      for (std::size_t j = 0; j < ny; ++j)
	for (std::size_t k = 0; k < nz; ++k)
	  arraySoundSpeed[i][j][k] = sqrt(gamma * arrayU[i][j][k][varPrs] / arrayU[i][j][k][varDen]);
  }

  // total enthalphy, not static enthalpy
  void Fluid::enthalpy() {
    for (std::size_t i = 0; i < nx ; ++i)
      for (std::size_t j = 0; j < ny; ++j)
	for (std::size_t k = 0; k < nz; ++k)
	  arrayH[i][j][k] = (arrayU[i][j][k][varEng] + arrayU[i][j][k][varPrs]) / arrayU[i][j][k][varDen];
  }

  void Fluid::initialCondition() {
    if (leftType == 1 || leftType == 2) { // normal shock w/ and w/o Rankine-Hugoniot conditions
      for (std::size_t i = 0; i < nx; ++i)
	for (std::size_t j = 0; j < ny; ++j)
	  for (std::size_t k = 0; k < nz; ++k) {
	    if (arrayGridCoord[i][j][k][2] <= z2L) {
	      arrayU[i][j][k][varDen] = rhoL;
	      arrayU[i][j][k][varPrs] = pL;
	      arrayU[i][j][k][varVel[0]] = 0;
	      arrayU[i][j][k][varVel[1]] = 0;
	      arrayU[i][j][k][varVel[2]] = uL;
	    } else {
	      arrayU[i][j][k][varDen] = rhoR;
	      arrayU[i][j][k][varPrs] = pR;
	      arrayU[i][j][k][varVel[0]] = 0;
	      arrayU[i][j][k][varVel[1]] = 0;
	      arrayU[i][j][k][varVel[2]] = uR;
	    }
	  }
    }
    else if (leftType == 3) { // normal shock in x, y, z directions
      for (std::size_t i = 0; i < nx; ++i)
	for (std::size_t j = 0; j < ny; ++j)
	  for (std::size_t k = 0; k < nz; ++k) {
	    if ( arrayGridCoord[i][j][k][2] >= z1L && arrayGridCoord[i][j][k][2] <= z2L &&
		 arrayGridCoord[i][j][k][0] >= x1L && arrayGridCoord[i][j][k][0] <= x2L &&
		 arrayGridCoord[i][j][k][1] >= y1L && arrayGridCoord[i][j][k][1] <= y2L) {
	      arrayU[i][j][k][varDen] = rhoL;
	      arrayU[i][j][k][varPrs] = pL;
	      arrayU[i][j][k][varVel[0]] = 0;
	      arrayU[i][j][k][varVel[1]] = 0;
	      arrayU[i][j][k][varVel[2]] = 0;
	    } else {
	      arrayU[i][j][k][varDen] = rhoR;
	      arrayU[i][j][k][varPrs] = pR;
	      arrayU[i][j][k][varVel[0]] = 0;
	      arrayU[i][j][k][varVel[1]] = 0;
	      arrayU[i][j][k][varVel[2]] = 0;
	    }
	  }
    } else if (leftType == 4) { // spherical shock
      for (std::size_t i = 0; i < nx; ++i)
	for (std::size_t j = 0; j < ny; ++j)
	  for (std::size_t k = 0; k < nz; ++k) {
	    REAL radius = sqrt(pow(arrayGridCoord[i][j][k][0]-x0L,2) + pow(arrayGridCoord[i][j][k][1]-y0L,2) + pow(arrayGridCoord[i][j][k][2]-z0L,2));
	    if ( radius <= r0L) {
	      arrayU[i][j][k][varDen] = rhoL;
	      arrayU[i][j][k][varPrs] = pL;
	      arrayU[i][j][k][varVel[0]] = uL*(arrayGridCoord[i][j][k][0]-x0L)/radius;
	      arrayU[i][j][k][varVel[1]] = uL*(arrayGridCoord[i][j][k][1]-y0L)/radius;
	      arrayU[i][j][k][varVel[2]] = uL*(arrayGridCoord[i][j][k][2]-z0L)/radius;
	    } else {
	      arrayU[i][j][k][varDen] = rhoR;
	      arrayU[i][j][k][varPrs] = pR;
	      arrayU[i][j][k][varVel[0]] = 0;
	      arrayU[i][j][k][varVel[1]] = 0;
	      arrayU[i][j][k][varVel[2]] = 0;
	    }
	  }
    } else if (leftType == 5) { // normal shock z directions, three initial zones
      for (std::size_t i = 0; i < nx; ++i)
	for (std::size_t j = 0; j < ny; ++j)
	  for (std::size_t k = 0; k < nz; ++k) {
	    if ( arrayGridCoord[i][j][k][2] >= z1L && arrayGridCoord[i][j][k][2] <= z2L ) {
	      arrayU[i][j][k][varDen] = rhoL;
	      arrayU[i][j][k][varPrs] = pL;
	      arrayU[i][j][k][varVel[0]] = 0;
	      arrayU[i][j][k][varVel[1]] = 0;
	      arrayU[i][j][k][varVel[2]] = 0;
	    } else if ( arrayGridCoord[i][j][k][2] > z2L) {
	      arrayU[i][j][k][varDen] = rhoR;
	      arrayU[i][j][k][varPrs] = pR;
	      arrayU[i][j][k][varVel[0]] = 0;
	      arrayU[i][j][k][varVel[1]] = 0;
	      arrayU[i][j][k][varVel[2]] = 0;
	    } else if ( arrayGridCoord[i][j][k][2] < z1L) {
	      arrayU[i][j][k][varDen] = rhoBL;
	      arrayU[i][j][k][varPrs] = pBL;
	      arrayU[i][j][k][varVel[0]] = 0;
	      arrayU[i][j][k][varVel[1]] = 0;
	      arrayU[i][j][k][varVel[2]] = 0;
	    }
	  }
    }

    WtoU();
  }

  void Fluid::RankineHugoniot() { // Rankine-Hugoniot conditions
    if (leftType == 1) {
      shockSpeed = Mach*sqrt(gamma*pR/rhoR);
      rhoL = ( pow(rhoR*(shockSpeed-uR),2)*(1+gamma) ) / ( rhoR*pow(shockSpeed-uR,2)*(gamma-1) + 2*pR*gamma);
      pL   = (pR*(1-gamma)+2*rhoR*pow(shockSpeed-uR,2)) / (1+gamma);
      uL   = ( rhoR*(shockSpeed-uR)*(2*shockSpeed + uR*(gamma-1)) - 2*pR*gamma ) / (rhoR * (shockSpeed-uR) * (1+gamma));
      debugInf << std::setw(OWID) << "shockSpeed" << std::setw(OWID) << shockSpeed << std::endl;
      debugInf << std::setw(OWID) << "rhoL" << std::setw(OWID) << rhoL << std::endl;
      debugInf << std::setw(OWID) << "pL" << std::setw(OWID) << pL << std::endl;
      debugInf << std::setw(OWID) << "uL" << std::setw(OWID) << uL << std::endl << std::endl;
    } else {
      Mach = uL / sqrt(gamma*pL/rhoL);
      debugInf << std::setw(OWID) << "MachL" << std::setw(OWID) << Mach << std::endl << std::endl;
    }
  }

  void Fluid::flux(std::size_t idim, std::vector<Particle *> &ptcls) {
    for (std::size_t i = 0; i < nx; ++i)
      for (std::size_t j = 0; j < ny; ++j)
	for (std::size_t k = 0; k < nz; ++k) {
	  arrayFlux[i][j][k][varDen]    = arrayUtmp[i][j][k][varDen] * arrayUtmp[i][j][k][varVel[0]]; // rho*u
	  arrayFlux[i][j][k][varMom[0]] = arrayUtmp[i][j][k][varDen] * pow(arrayUtmp[i][j][k][varVel[0]],2) + arrayUtmp[i][j][k][varPrs]; // rho*u^2 + p
	  arrayFlux[i][j][k][varMom[1]] = arrayUtmp[i][j][k][varDen] * arrayUtmp[i][j][k][varVel[0]] * arrayUtmp[i][j][k][varVel[1]]; // rho*u*v
	  arrayFlux[i][j][k][varMom[2]] = arrayUtmp[i][j][k][varDen] * arrayUtmp[i][j][k][varVel[0]] * arrayUtmp[i][j][k][varVel[2]]; // rho*u*w
	  arrayFlux[i][j][k][varEng]    = arrayUtmp[i][j][k][varVel[0]] * (arrayUtmp[i][j][k][varEng] + arrayUtmp[i][j][k][varPrs]);  // u*(E + p)
	}  

    // for cells that are enclosed by particle volumes
    for (std::vector<Particle *>::const_iterator it = ptcls.begin(); it != ptcls.end(); ++it) {
      std::vector< std::vector<REAL> > fluidGrid = (*it)->getFluidGrid();
      for (std::size_t iter = 0; iter < fluidGrid.size(); ++iter) {
	std::size_t i = static_cast<std::size_t> (fluidGrid[iter][0]);
	std::size_t j = static_cast<std::size_t> (fluidGrid[iter][1]);
	std::size_t k = static_cast<std::size_t> (fluidGrid[iter][2]);
	REAL u0[3];
	REAL coordX = arrayGridCoord[i][j][k][0];
	REAL coordY = arrayGridCoord[i][j][k][1];
	REAL coordZ = arrayGridCoord[i][j][k][2];
	Vec dist = Vec(coordX, coordY, coordZ) - (*it)->getCurrPos();
	Vec omgar = (*it)->getCurrOmga() % dist; // w X r = omga % dist, where % is overloaded as cross product
	u0[0] = (*it)->getCurrVeloc().getX() + omgar.getX(); 
	u0[1] = (*it)->getCurrVeloc().getY() + omgar.getY(); 
	u0[2] = (*it)->getCurrVeloc().getZ() + omgar.getZ();

	// continuity equation modified on porosity
	arrayFlux[i][j][k][varDen]    = 1.0/porosity * arrayUtmp[i][j][k][varDen] * (arrayUtmp[i][j][k][varVel[0]] - u0[idim]) ; // rho*(u-u0) / porosity

	// momentum equations modified on porosity
	arrayFlux[i][j][k][varMom[0]] = 1.0/porosity * arrayUtmp[i][j][k][varDen] * pow(arrayUtmp[i][j][k][varVel[0]],2) + arrayUtmp[i][j][k][varPrs]; // rho*u^2 / porosity + p
	arrayFlux[i][j][k][varMom[1]] = 1.0/porosity * arrayUtmp[i][j][k][varDen] * arrayUtmp[i][j][k][varVel[0]] * arrayUtmp[i][j][k][varVel[1]]; // rho*u*v / porosity
	arrayFlux[i][j][k][varMom[2]] = 1.0/porosity * arrayUtmp[i][j][k][varDen] * arrayUtmp[i][j][k][varVel[0]] * arrayUtmp[i][j][k][varVel[2]]; // rho*u*w / porosity

	// momentum equations modified on porosity
	//arrayFlux[i][j][k][varEng]    = 1.0/porosity * arrayUtmp[i][j][k][varVel[0]] * (arrayUtmp[i][j][k][varEng] + arrayUtmp[i][j][k][varPrs]);  // u*(E + p) / porosity
      }
    }
  }

  void Fluid::RoeFlux(REAL uL[], REAL uR[], REAL FL[], REAL FR[], REAL HL, REAL HR, std::size_t idim, std::size_t it, std::size_t jt, std::size_t kt) {

    // it, jt, kt defined at cell faces
    if (uL[varDen] < 0 || uR[varDen] < 0)
      debugInf << std::setw(OWID) << "RoeFlux:den<0"
	       << std::setw(5) << it
	       << std::setw(5) << jt
	       << std::setw(5) << kt
	       << std::setw(5) << idim
	       << std::setw(5) << (int) uL[varMsk]
	       << std::setw(5) << (int) uR[varMsk]
	       << std::setw(OWID) << uL[varDen]
	       << std::setw(OWID) << uR[varDen]
	       << std::setw(OWID) << uL[varEng]
	       << std::setw(OWID) << uR[varEng]
	       << std::endl;

    REAL avgRho =  sqrt(uL[varDen]*uR[varDen]);
    REAL avgH   = (sqrt(uL[varDen])*HL + sqrt(uR[varDen])*HR)/(sqrt(uL[varDen]) + sqrt(uR[varDen]));
    REAL avgU   = (sqrt(uL[varDen])*uL[varVel[0]] + sqrt(uR[varDen])*uR[varVel[0]])/(sqrt(uL[varDen]) + sqrt(uR[varDen]));
    REAL avgV   = (sqrt(uL[varDen])*uL[varVel[1]] + sqrt(uR[varDen])*uR[varVel[1]])/(sqrt(uL[varDen]) + sqrt(uR[varDen]));
    REAL avgW   = (sqrt(uL[varDen])*uL[varVel[2]] + sqrt(uR[varDen])*uR[varVel[2]])/(sqrt(uL[varDen]) + sqrt(uR[varDen]));
    REAL avgh   = avgH - 0.5*(avgU*avgU + avgV*avgV + avgW*avgW); // static specific enthalpy
    REAL avgSoundSpeed = sqrt((gamma-1)*avgh);

    if (avgh <= 0) {
      debugInf << std::setw(OWID) << "RoeFlux:avgh<=0"
	       << std::setw(5) << it
	       << std::setw(5) << jt
	       << std::setw(5) << kt
	       << std::setw(5) << idim
	       << std::setw(5) << (int) uL[varMsk]
	       << std::setw(5) << (int) uR[varMsk]
	       << std::setw(OWID) << HL
	       << std::setw(OWID) << HR
	       << std::setw(OWID) << uL[varDen]
	       << std::setw(OWID) << uR[varDen]
	       << std::setw(OWID) << uL[varVel[0]]
	       << std::setw(OWID) << uR[varVel[0]]
	       << std::setw(OWID) << uL[varVel[1]]
	       << std::setw(OWID) << uR[varVel[1]]
	       << std::setw(OWID) << uL[varVel[2]]
	       << std::setw(OWID) << uR[varVel[2]]
	       << std::setw(OWID) << uL[varEng]
	       << std::setw(OWID) << uR[varEng]
	       << std::setw(OWID) << avgH
	       << std::setw(OWID) << avgU
	       << std::setw(OWID) << avgV
	       << std::setw(OWID) << avgW
	       << std::setw(OWID) << avgSoundSpeed
	       << std::setw(OWID) << uL[varPrs]
	       << std::setw(OWID) << uR[varPrs]
	       << std::endl;
    }

    REAL eigen[5];
    eigen[varDen]    = avgU - avgSoundSpeed;
    eigen[varMom[0]] = eigen[varMom[1]] = eigen[varMom[2]] = avgU;
    eigen[varEng]    = avgU + avgSoundSpeed;

    REAL avgWaveStr[5], du[9];
    for (std::size_t i = 0; i < nVar; ++i)
      du[i] = uR[i] - uL[i];

    avgWaveStr[varDen]    = (du[varPrs] - avgRho*avgSoundSpeed*du[varVel[0]]) / (2*avgSoundSpeed*avgSoundSpeed);
    avgWaveStr[varMom[0]] = du[varDen] - du[varPrs]/(avgSoundSpeed*avgSoundSpeed);
    avgWaveStr[varMom[1]] = avgRho * du[varVel[1]];
    avgWaveStr[varMom[2]] = avgRho * du[varVel[2]];
    avgWaveStr[varEng]    = (du[varPrs] + avgRho*avgSoundSpeed*du[varVel[0]]) / (2*avgSoundSpeed*avgSoundSpeed);
    
    REAL avgK[5][5]; // right eigenvectors
    avgK[varDen][varDen]    = 1;
    avgK[varMom[0]][varDen] = avgU - avgSoundSpeed;
    avgK[varMom[1]][varDen] = avgV;
    avgK[varMom[2]][varDen] = avgW;
    avgK[varEng][varDen]    = avgH - avgU * avgSoundSpeed;

    avgK[varDen][varMom[0]]    = 1;
    avgK[varMom[0]][varMom[0]] = avgU;
    avgK[varMom[1]][varMom[0]] = avgV;
    avgK[varMom[2]][varMom[0]] = avgW;
    avgK[varEng][varMom[0]]    = 0.5 * (avgU*avgU + avgV*avgV + avgW*avgW);

    avgK[varDen][varMom[1]]    = 0;
    avgK[varMom[0]][varMom[1]] = 0;
    avgK[varMom[1]][varMom[1]] = 1;
    avgK[varMom[2]][varMom[1]] = 0;
    avgK[varEng][varMom[1]]    = avgV;

    avgK[varDen][varMom[2]]    = 0;
    avgK[varMom[0]][varMom[2]] = 0;
    avgK[varMom[1]][varMom[2]] = 0;
    avgK[varMom[2]][varMom[2]] = 1;
    avgK[varEng][varMom[2]]    = avgW;

    avgK[varDen][varEng]    = 1;
    avgK[varMom[0]][varEng] = avgU + avgSoundSpeed;
    avgK[varMom[1]][varEng] = avgV;
    avgK[varMom[2]][varEng] = avgW;
    avgK[varEng][varEng]    = avgH + avgU * avgSoundSpeed;

    REAL RF[5];
    for (std::size_t i = 0; i < nInteg; ++i)
      RF[i] = 0.5*(FL[i] + FR[i]);

    for (std::size_t ie = 0; ie < nInteg; ++ie) {
      for (std::size_t je = 0; je < nInteg; ++je)
	RF[ie] -= 0.5 * avgWaveStr[je] * fabs(eigen[je]) * avgK[ie][je];
      arrayRoeFlux[it][jt][kt][ie][idim] = RF[ie];
    }
  }

  void Fluid::UtoW() {// converting conserved variables into primitive
    for (std::size_t i = 0; i < nx; ++i)
      for (std::size_t j = 0; j < ny; ++j)
	for (std::size_t k = 0; k < nz; ++k) {

	  for (std::size_t m = 0; m < nDim; ++m)
	    arrayU[i][j][k][varVel[m]] = arrayU[i][j][k][varMom[m]] / arrayU[i][j][k][varDen];

	  arrayU[i][j][k][varPrs] = 0;
	  for (std::size_t m = 0; m < nDim; ++m)
	    arrayU[i][j][k][varPrs] += arrayU[i][j][k][varDen] * pow(arrayU[i][j][k][varVel[m]],2)/2 ;
	  arrayU[i][j][k][varPrs] = (arrayU[i][j][k][varEng] - arrayU[i][j][k][varPrs]) * (gamma-1);
	}
  }

  void Fluid::WtoU() { // converting primitive variables into conserved
    for (std::size_t i = 0; i < nx; ++i)
      for (std::size_t j = 0; j < ny; ++j)
	for (std::size_t k = 0; k < nz; ++k) {
	  for (std::size_t m = 0; m < nDim; ++m)
	    arrayU[i][j][k][varMom[m]] = arrayU[i][j][k][varDen] * arrayU[i][j][k][varVel[m]];

	  arrayU[i][j][k][varEng] = 0;
	  for (std::size_t m = 0; m < nDim; ++m)
	    arrayU[i][j][k][varEng] += arrayU[i][j][k][varDen] * pow(arrayU[i][j][k][varVel[m]],2)/2 ;
	  arrayU[i][j][k][varEng] += arrayU[i][j][k][varPrs] / (gamma-1);
	}
  }

  void Fluid::getParticleInfo(std::vector<Particle *> &ptcls) {
    for (std::vector<Particle*>::const_iterator it = ptcls.begin(); it != ptcls.end(); ++it)
      (*it)->clearFluidGrid();

    // 0 ~ (n-1), including boundaries
    for (std::size_t i = 0; i < arrayGridCoord.size() ; ++i)
      for (std::size_t j = 0; j <  arrayGridCoord[i].size(); ++j)
	for (std::size_t k = 0; k <  arrayGridCoord[i][j].size(); ++k) {

	  arrayU[i][j][k][varMsk] = 0;
	  REAL coordX = arrayGridCoord[i][j][k][0];
	  REAL coordY = arrayGridCoord[i][j][k][1];
	  REAL coordZ = arrayGridCoord[i][j][k][2];

	  for (std::vector<Particle*>::iterator it = ptcls.begin(); it != ptcls.end(); ++it)
	    if ( (*it)->surfaceError(Vec(coordX, coordY, coordZ)) <= 0 ) { // inside particle surface
	      arrayU[i][j][k][varMsk] = 1; 
	      (*it)->recordFluidGrid(i, j, k);
	    }
	}
  }

  void Fluid::calcParticleForce(std::vector<Particle *> &ptcls, std::ofstream &ofs) {
    // must clear forces each loop, otherwise Fluid::plot prints wrong values;
    // but Fluid::penalize works OK since it uses masks.
    for (std::size_t i = 0; i < nx ; ++i)
      for (std::size_t j = 0; j < ny; ++j)
	for (std::size_t k = 0; k < nz; ++k)
	  for (std::size_t m = 0; m < nDim; ++m) {
	    arrayPenalForce[i][j][k][m] = 0;
	    arrayPressureForce[i][j][k][m] = 0;
	  }

    for (std::vector<Particle *>::const_iterator it = ptcls.begin(); it != ptcls.end(); ++it) {
      REAL etaBx = 8.0/3.0 * (*it)->getA() / Cd; // local direction x (i.e. a)
      REAL etaBy = 8.0/3.0 * (*it)->getB() / Cd; // local direction y (i.e. b)
      REAL etaBz = 8.0/3.0 * (*it)->getC() / Cd; // local direction z (i.e. c)

      Vec penalForce  = 0, presForce  = 0;
      Vec penalMoment = 0, presMoment = 0;
      REAL avgDen = 0, avgVel = 0, avgPrs = 0;
      std::vector< std::vector<REAL> > fluidGrid = (*it)->getFluidGrid();
      for (std::size_t iter = 0; iter < fluidGrid.size(); ++iter) {
	std::size_t i = static_cast<std::size_t> (fluidGrid[iter][0]);
	std::size_t j = static_cast<std::size_t> (fluidGrid[iter][1]);
	std::size_t k = static_cast<std::size_t> (fluidGrid[iter][2]);
	avgDen += arrayU[i][j][k][varDen];
	avgVel += sqrt( pow(arrayU[i][j][k][varVel[0]],2) + pow(arrayU[i][j][k][varVel[1]],2) + pow(arrayU[i][j][k][varVel[2]],2) );
	avgPrs += arrayU[i][j][k][varPrs];

	REAL coordX = arrayGridCoord[i][j][k][0];
	REAL coordY = arrayGridCoord[i][j][k][1];
	REAL coordZ = arrayGridCoord[i][j][k][2];

	REAL uxFluid = arrayU[i][j][k][varVel[0]];
	REAL uyFluid = arrayU[i][j][k][varVel[1]];
	REAL uzFluid = arrayU[i][j][k][varVel[2]];

	Vec dist = Vec(coordX, coordY, coordZ) - (*it)->getCurrPos();
	Vec omgar = (*it)->getCurrOmga() % dist; // w X r = omga % dist, where % is overloaded as cross product

	REAL ux = (*it)->getCurrVeloc().getX() + omgar.getX(); 
	REAL uy = (*it)->getCurrVeloc().getY() + omgar.getY(); 
	REAL uz = (*it)->getCurrVeloc().getZ() + omgar.getZ();

	// principal axis decomposition
	Vec globalDelta = Vec(fabs(uxFluid - ux)*(uxFluid - ux), fabs(uyFluid - uy)*(uyFluid - uy), fabs(uzFluid - uz)*(uzFluid - uz));
	Vec localDelta = (*it)->globalToLocal(globalDelta);
	Vec localPenal, globalPenal;
	// localDelta needs to project in local frame in order to calculate local penalization forces
	localPenal.setX(arrayU[i][j][k][varDen] * localDelta.getX() / etaBx);
	localPenal.setY(arrayU[i][j][k][varDen] * localDelta.getY() / etaBy);
	localPenal.setZ(arrayU[i][j][k][varDen] * localDelta.getZ() / etaBz);
	globalPenal = (*it)->localToGlobal(localPenal);
	// one grid could have multiple particles intruded, +=, not =
	arrayPenalForce[i][j][k][0] += globalPenal.getX(); 
	arrayPenalForce[i][j][k][1] += globalPenal.getY();
	arrayPenalForce[i][j][k][2] += globalPenal.getZ();

	// restrict pressure gradient grids
 	if (i > 0 && i < nx-1 && j > 0 && j < ny-1 && k > 0 && k < nz-1 ) { // do not use (i-1) for std::size_t because (i-1) is postive when i=0
	  arrayPressureForce[i][j][k][0] = -(arrayU[i+1][j][k][varPrs] - arrayU[i-1][j][k][varPrs])/(2*dx);
	  arrayPressureForce[i][j][k][1] = -(arrayU[i][j+1][k][varPrs] - arrayU[i][j-1][k][varPrs])/(2*dy);
	  arrayPressureForce[i][j][k][2] = -(arrayU[i][j][k+1][varPrs] - arrayU[i][j][k-1][varPrs])/(2*dz);
	}

	penalForce += Vec(arrayPenalForce[i][j][k][0], arrayPenalForce[i][j][k][1], arrayPenalForce[i][j][k][2]);
	presForce  += Vec(arrayPressureForce[i][j][k][0], arrayPressureForce[i][j][k][1], arrayPressureForce[i][j][k][2]);

	// r X F,  % is overloaded as cross product
	penalMoment += dist % Vec(arrayPenalForce[i][j][k][0], arrayPenalForce[i][j][k][1], arrayPenalForce[i][j][k][2]);
	presMoment  += dist % Vec(arrayPressureForce[i][j][k][0], arrayPressureForce[i][j][k][1], arrayPressureForce[i][j][k][2]);
      } // end of fluidGrid loop

      avgDen /= fluidGrid.size();
      avgVel /= fluidGrid.size();
      avgPrs /= fluidGrid.size();

      penalForce *= dx*dy*dz;
      presForce  *= dx*dy*dz;
      (*it)->addForce(penalForce);
      (*it)->addForce(presForce);

      penalMoment *= dx*dy*dz;
      presMoment  *= dx*dy*dz;
      (*it)->addMoment(penalMoment);
      (*it)->addMoment(presMoment);

      if ((*it)->getId() == 1) {
	ofs << std::setw(OWID) << iteration
	    << std::setw(OWID) << timeAccrued

	    << std::setw(OWID) << penalForce.getX()
	    << std::setw(OWID) << penalForce.getY()
	    << std::setw(OWID) << penalForce.getZ()
	    << std::setw(OWID) << presForce.getX()
	    << std::setw(OWID) << presForce.getY()
	    << std::setw(OWID) << presForce.getZ()

	    << std::setw(OWID) << penalMoment.getX()
	    << std::setw(OWID) << penalMoment.getY()
	    << std::setw(OWID) << penalMoment.getZ()
	    << std::setw(OWID) << presMoment.getX()
	    << std::setw(OWID) << presMoment.getY()
	    << std::setw(OWID) << presMoment.getZ()

	    << std::setw(OWID) << (*it)->getAccel().getX()
	    << std::setw(OWID) << (*it)->getAccel().getY()
	    << std::setw(OWID) << (*it)->getAccel().getZ()

	    << std::setw(OWID) << (*it)->getCurrVeloc().getX()
	    << std::setw(OWID) << (*it)->getCurrVeloc().getY()
	    << std::setw(OWID) << (*it)->getCurrVeloc().getZ()

	    << std::setw(OWID) << avgDen
	    << std::setw(OWID) << avgVel
	    << std::setw(OWID) << avgPrs

	    << std::endl;    
      } 
    } // end of particle loop  
  }

  void Fluid::plot(const char *str) const {
    std::ofstream ofs(str);
    if(!ofs) { std::cout << "stream error: Fluid::plot" << std::endl; exit(-1); }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    
    ofs	<< std::setw(OWID) << "VARIABLES = \"x\""
	<< std::setw(OWID) << "\"y\""
	<< std::setw(OWID) << "\"z\""
	<< std::setw(OWID) << "\"Mach\""
	<< std::setw(OWID) << "\"density\""
	<< std::setw(OWID) << "\"momentumX\""
	<< std::setw(OWID) << "\"momentumY\""
	<< std::setw(OWID) << "\"momentumZ\""
	<< std::setw(OWID) << "\"energy\""
	<< std::setw(OWID) << "\"velocityX\""
	<< std::setw(OWID) << "\"velocityY\""
	<< std::setw(OWID) << "\"velocityZ\""
	<< std::setw(OWID) << "\"pressure\""
	<< std::setw(OWID) << "\"temperature\""
	<< std::setw(OWID) << "\"mask\""
	<< std::setw(OWID) << "\"penalFx\""
	<< std::setw(OWID) << "\"penalFy\""
	<< std::setw(OWID) << "\"penalFz\""
	<< std::setw(OWID) << "\"pressureFx\""
	<< std::setw(OWID) << "\"pressureFy\""
	<< std::setw(OWID) << "\"pressureFz\""
	<< std::endl;

    ofs << "ZONE I=" << nx -2
	<< ", J=" << ny -2
	<< ", K=" << nz -2
	<< ", DATAPACKING=POINT"
	<< std::endl;

    for (std::size_t k = 1; k < nz - 1; ++k)
      for (std::size_t j = 1; j < ny - 1; ++j)
	for (std::size_t i = 1; i < nx -1; ++i) {
	  ofs << std::setw(OWID) << arrayGridCoord[i][j][k][0]
	      << std::setw(OWID) << arrayGridCoord[i][j][k][1]
	      << std::setw(OWID) << arrayGridCoord[i][j][k][2]
	      << std::setw(OWID) << vfabs( Vec(arrayU[i][j][k][varVel[0]], arrayU[i][j][k][varVel[1]], arrayU[i][j][k][varVel[2]]) ) / arraySoundSpeed[i][j][k]
	      << std::setw(OWID) << arrayU[i][j][k][varDen]
	      << std::setw(OWID) << arrayU[i][j][k][varMom[0]]
	      << std::setw(OWID) << arrayU[i][j][k][varMom[1]]
	      << std::setw(OWID) << arrayU[i][j][k][varMom[2]]
	      << std::setw(OWID) << arrayU[i][j][k][varEng]
	      << std::setw(OWID) << arrayU[i][j][k][varVel[0]]
	      << std::setw(OWID) << arrayU[i][j][k][varVel[1]]
	      << std::setw(OWID) << arrayU[i][j][k][varVel[2]]
	      << std::setw(OWID) << arrayU[i][j][k][varPrs]
	      << std::setw(OWID) << arrayU[i][j][k][varPrs]/(Rs*arrayU[i][j][k][varDen])
	      << std::setw(OWID) << arrayU[i][j][k][varMsk]  
	      << std::setw(OWID) << arrayPenalForce[i][j][k][0] 
	      << std::setw(OWID) << arrayPenalForce[i][j][k][1] 
	      << std::setw(OWID) << arrayPenalForce[i][j][k][2]
	      << std::setw(OWID) << arrayPressureForce[i][j][k][0] 
	      << std::setw(OWID) << arrayPressureForce[i][j][k][1] 
	      << std::setw(OWID) << arrayPressureForce[i][j][k][2] 
	      << std::endl;
	}

    ofs.close();
  }
  
} // name space dem
