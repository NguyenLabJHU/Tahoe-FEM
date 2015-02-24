#include "Fluid.h"
#include "const.h"
#include <cmath>
#include <algorithm>

namespace dem {

  char *combineStr(char *cstr, const char *str, std::size_t num, std::size_t width) {
    std::string obj(str);
    std::stringstream ss;
    ss << std::setw(width) << std::setfill('0') << std::right << num;
    obj += ss.str();
    return strcpy( cstr, obj.c_str() );
  }

  const REAL Fluid::Rs;

  void Fluid::initParameter(Rectangle &container, Gradation &gradation) {

    ptclGrid = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["ptclGrid"]);
    Cd   = dem::Parameter::getSingleton().parameter["Cd"];
    porosity = dem::Parameter::getSingleton().parameter["porosity"];
    Cdi = dem::Parameter::getSingleton().parameter["Cdi"];
    RK = dem::Parameter::getSingleton().parameter["RK"];
    CFL = dem::Parameter::getSingleton().parameter["CFL"];
    gama = dem::Parameter::getSingleton().parameter["airGamma"];
    rhoR = dem::Parameter::getSingleton().parameter["rightDensity"];
    pR   = dem::Parameter::getSingleton().parameter["rightPressure"];
    uR   = dem::Parameter::getSingleton().parameter["rightVelocity"];
    
    leftType = static_cast<int> (dem::Parameter::getSingleton().parameter["leftType"]);    
    x1F = dem::Parameter::getSingleton().parameter["x1F"];
    y1F = dem::Parameter::getSingleton().parameter["y1F"];
    z1F = dem::Parameter::getSingleton().parameter["z1F"];
    x2F = dem::Parameter::getSingleton().parameter["x2F"];
    y2F = dem::Parameter::getSingleton().parameter["y2F"];
    z2F = dem::Parameter::getSingleton().parameter["z2F"];
    arrayBC[0] = dem::Parameter::getSingleton().parameter["x1Reflecting"];
    arrayBC[1] = dem::Parameter::getSingleton().parameter["x2Reflecting"];
    arrayBC[2] = dem::Parameter::getSingleton().parameter["y1Reflecting"];
    arrayBC[3] = dem::Parameter::getSingleton().parameter["y2Reflecting"];
    arrayBC[4] = dem::Parameter::getSingleton().parameter["z1Reflecting"];
    arrayBC[5] = dem::Parameter::getSingleton().parameter["z2Reflecting"];

    if (leftType == 1) {
      z2L = dem::Parameter::getSingleton().parameter["z2L"];
      MachShock= dem::Parameter::getSingleton().parameter["shockMach"];
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

    printPtcls = dem::Parameter::getSingleton().cfdPrintPtcls;

    REAL minR = gradation.getPtclMinRadius();
    gridDz = (minR * 2) / ptclGrid;
    gridNz = static_cast<std::size_t> (ceil((z2F - z1F) / gridDz));
    gridDz = (z2F - z1F) / gridNz;
    gridDx = gridDz;
    gridDy = gridDz;
    gridNx = static_cast<std::size_t> (ceil((x2F - x1F) / gridDx));
    gridNy = static_cast<std::size_t> (ceil((y2F - y1F) / gridDy));

    gridNx += 2; // add two boundary cells
    gridNy += 2;
    gridNz += 2;

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
    debugInf << std::setw(OWID) << "Cdi" << std::setw(OWID) << Cdi << std::endl;
    debugInf << std::setw(OWID) << "Runge-Kutta" << std::setw(OWID) << (int) RK << std::endl;
    debugInf << std::setw(OWID) << "CFL" << std::setw(OWID) << CFL << std::endl;
    debugInf << std::setw(OWID) << "gama" << std::setw(OWID) << gama << std::endl;
    debugInf << std::setw(OWID) << "rhoR" << std::setw(OWID) << rhoR << std::endl;
    debugInf << std::setw(OWID) << "pR" << std::setw(OWID) << pR << std::endl;
    debugInf << std::setw(OWID) << "uR" << std::setw(OWID) << uR << std::endl;

    debugInf << std::setw(OWID) << "gridDx" << std::setw(OWID) << gridDx << std::endl;
    debugInf << std::setw(OWID) << "gridDy" << std::setw(OWID) << gridDy << std::endl;
    debugInf << std::setw(OWID) << "gridDz" << std::setw(OWID) << gridDz << std::endl;
    debugInf << std::setw(OWID) << "gridNx" << std::setw(OWID) << gridNx << std::endl;
    debugInf << std::setw(OWID) << "gridNy" << std::setw(OWID) << gridNy << std::endl;
    debugInf << std::setw(OWID) << "gridNz" << std::setw(OWID) << gridNz << std::endl;

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
    debugInf << std::setw(OWID) << "printPtclNum" << std::setw(OWID) << printPtcls.size() << std::endl;

    if (leftType == 1) {
      debugInf << std::setw(OWID) << "z2L" << std::setw(OWID) << z2L << std::endl;
      debugInf << std::setw(OWID) << "shockMach" << std::setw(OWID) << MachShock << std::endl;
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

    REAL aNum = (z2F - z1F) / gridDz;
    if (aNum != ceil(aNum))
      debugInf << "fluid domain in z direction is not precisely grided!" << std::endl;

    aNum = (x2F - x1F) / gridDx;
    if (aNum != ceil(aNum))
      debugInf << "fluid domain in x direction is not precisely grided!" << std::endl;
    aNum = (y2F - y1F) / gridDy;
    if (aNum != ceil(aNum))
      debugInf << "fluid domain in y direction is not precisely grided!" << std::endl;

    /*
    debugInf << "nVar " << nVar << std::endl;
    debugInf << "nInteg " << nInteg << std::endl;
    debugInf << "varDen " << varDen << std::endl;    
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

    // gridNx, gridNy, gridNz, nDim
    arrayGridCoord.resize(gridNx);
    for (std::size_t i = 0; i < arrayGridCoord.size(); ++i) {
      arrayGridCoord[i].resize(gridNy);
      for (std::size_t j = 0; j < arrayGridCoord[i].size(); ++j) {
	arrayGridCoord[i][j].resize(gridNz);
	for (std::size_t k = 0; k < arrayGridCoord[i][j].size(); ++k) 
	  arrayGridCoord[i][j][k].resize(nDim);
      }
    }
    // coordinates
    for (std::size_t i = 0; i < arrayGridCoord.size(); ++i)
      for (std::size_t j = 0; j < arrayGridCoord[i].size(); ++j)
	for (std::size_t k = 0; k < arrayGridCoord[i][j].size(); ++k) {
	  arrayGridCoord[i][j][k][0] = (x1F - gridDx/2) + i * gridDx;
	  arrayGridCoord[i][j][k][1] = (y1F - gridDy/2) + j * gridDy;
	  arrayGridCoord[i][j][k][2] = (z1F - gridDz/2) + k * gridDz;
	}

    // gridNx, gridNy, gridNz, nDim
    arrayPenalForce.resize(gridNx);
    for (std::size_t i = 0; i < arrayPenalForce.size(); ++i) {
      arrayPenalForce[i].resize(gridNy);
      for (std::size_t j = 0; j < arrayPenalForce[i].size(); ++j) {
	arrayPenalForce[i][j].resize(gridNz);
	for (std::size_t k = 0; k < arrayPenalForce[i][j].size(); ++k) 
	  arrayPenalForce[i][j][k].resize(nDim);
      }
    }

    // gridNx, gridNy, gridNz, nDim
    arrayPressureForce.resize(gridNx);
    for (std::size_t i = 0; i < arrayPressureForce.size(); ++i) {
      arrayPressureForce[i].resize(gridNy);
      for (std::size_t j = 0; j < arrayPressureForce[i].size(); ++j) {
	arrayPressureForce[i][j].resize(gridNz);
	for (std::size_t k = 0; k < arrayPressureForce[i][j].size(); ++k) 
	  arrayPressureForce[i][j][k].resize(nDim);
      }
    }

    // gridNx, gridNy, gridNz, nVar
    arrayU.resize(gridNx);
    for (std::size_t i = 0; i < arrayU.size(); ++i) {
      arrayU[i].resize(gridNy);
      for (std::size_t j = 0; j < arrayU[i].size(); ++j) {
	arrayU[i][j].resize(gridNz);
	for (std::size_t k = 0; k < arrayU[i][j].size(); ++k) 
	  arrayU[i][j][k].resize(nVar);
      }
    }

    // gridNx, gridNy, gridNz, nVar
    arrayURota.resize(gridNx);
    for (std::size_t i = 0; i < arrayURota.size(); ++i) {
      arrayURota[i].resize(gridNy);
      for (std::size_t j = 0; j < arrayURota[i].size(); ++j) {
	arrayURota[i][j].resize(gridNz);
	for (std::size_t k = 0; k < arrayURota[i][j].size(); ++k) 
	  arrayURota[i][j][k].resize(nVar);
      }
    }

    // gridNx, gridNy, gridNz, nVar
    arrayUPrev.resize(gridNx);
    for (std::size_t i = 0; i < arrayUPrev.size(); ++i) {
      arrayUPrev[i].resize(gridNy);
      for (std::size_t j = 0; j < arrayUPrev[i].size(); ++j) {
	arrayUPrev[i][j].resize(gridNz);
	for (std::size_t k = 0; k < arrayUPrev[i][j].size(); ++k) 
	  arrayUPrev[i][j][k].resize(nVar);
      }
    }

    // gridNx, gridNy, gridNz, nInteg
    arrayFlux.resize(gridNx);
    for (std::size_t i = 0; i < arrayFlux.size(); ++i) {
      arrayFlux[i].resize(gridNy);
      for (std::size_t j = 0; j < arrayFlux[i].size(); ++j) {
	arrayFlux[i][j].resize(gridNz);
	for (std::size_t k = 0; k < arrayFlux[i][j].size(); ++k) 
	  arrayFlux[i][j][k].resize(nInteg);
      }
    }

    // gridNx-1, gridNy-1, gridNz-1, nInteg, nDim
    arrayGodFlux.resize(gridNx-1);
    for (std::size_t i = 0; i < arrayGodFlux.size(); ++i) {
      arrayGodFlux[i].resize(gridNy-1);
      for (std::size_t j = 0; j < arrayGodFlux[i].size(); ++j) {
	arrayGodFlux[i][j].resize(gridNz-1);
	for (std::size_t k = 0; k < arrayGodFlux[i][j].size(); ++k) {
	  arrayGodFlux[i][j][k].resize(nInteg);
	  for (std::size_t m = 0; m < arrayGodFlux[i][j][k].size(); ++m)
	    arrayGodFlux[i][j][k][m].resize(nDim);
	}
      }
    }

    if (RK >= 1) {
      // gridNx-1, gridNy-1, gridNz-1, nInteg, nDim
      arrayGodFluxStep2.resize(gridNx-1);
      for (std::size_t i = 0; i < arrayGodFluxStep2.size(); ++i) {
	arrayGodFluxStep2[i].resize(gridNy-1);
	for (std::size_t j = 0; j < arrayGodFluxStep2[i].size(); ++j) {
	  arrayGodFluxStep2[i][j].resize(gridNz-1);
	  for (std::size_t k = 0; k < arrayGodFluxStep2[i][j].size(); ++k) {
	    arrayGodFluxStep2[i][j][k].resize(nInteg);
	    for (std::size_t m = 0; m < arrayGodFluxStep2[i][j][k].size(); ++m)
	      arrayGodFluxStep2[i][j][k][m].resize(nDim);
	  }
	}
      }
    }

    if (RK == 2) {
      // gridNx-1, gridNy-1, gridNz-1, nInteg, nDim
      arrayGodFluxStep3.resize(gridNx-1);
      for (std::size_t i = 0; i < arrayGodFluxStep3.size(); ++i) {
	arrayGodFluxStep3[i].resize(gridNy-1);
	for (std::size_t j = 0; j < arrayGodFluxStep3[i].size(); ++j) {
	  arrayGodFluxStep3[i][j].resize(gridNz-1);
	  for (std::size_t k = 0; k < arrayGodFluxStep3[i][j].size(); ++k) {
	    arrayGodFluxStep3[i][j][k].resize(nInteg);
	    for (std::size_t m = 0; m < arrayGodFluxStep3[i][j][k].size(); ++m)
	      arrayGodFluxStep3[i][j][k][m].resize(nDim);
	  }
	}
      }
    }

    // gridNx-1, gridNy-1, gridNz-1, nInteg
    arrayGodFluxTmp.resize(gridNx-1);
    for (std::size_t i = 0; i < arrayGodFluxTmp.size(); ++i) {
      arrayGodFluxTmp[i].resize(gridNy-1);
      for (std::size_t j = 0; j < arrayGodFluxTmp[i].size(); ++j) {
	arrayGodFluxTmp[i][j].resize(gridNz-1);
	for (std::size_t k = 0; k < arrayGodFluxTmp[i][j].size(); ++k) {
	  arrayGodFluxTmp[i][j][k].resize(nInteg);
	}
      }
    }
    
    // gridNx, gridNy, gridNz
    arrayH.resize(gridNx);
    for (std::size_t i = 0; i < arrayH.size(); ++i) {
      arrayH[i].resize(gridNy);
      for (std::size_t j = 0; j < arrayH[i].size(); ++j)
	arrayH[i][j].resize(gridNz);
    }

    // gridNx, gridNy, gridNz
    arraySoundSpeed.resize(gridNx);
    for (std::size_t i = 0; i < arraySoundSpeed.size(); ++i) {
      arraySoundSpeed[i].resize(gridNy);
      for (std::size_t j = 0; j < arraySoundSpeed[i].size(); ++j)
	arraySoundSpeed[i][j].resize(gridNz);
    }
  }

  void Fluid::initialize() {
    negPrsDen = false;
    RankineHugoniot();
    initialCondition(); 
    soundSpeed(); // for printing Mach number
    debugInf << std::setw(OWID) << "iteration" 
	     << std::setw(OWID) << "timeStep" 
	     << std::setw(OWID) << "timeAccrued" /*
	     << std::setw(OWID) << "uZMax"
	     << std::setw(OWID) << "uZMin"
	     << std::setw(OWID) << "soundSpeedMax"
	     << std::setw(OWID) << "(|uZ|+a)Max" */;
  }

  void Fluid::runOneStep(std::vector<Particle *> &ptcls) {
    if (!negPrsDen) arrayUPrev = arrayU; 
    else arrayU = arrayUPrev;
  label1: inteStep1(ptcls); // arrayU updated
    debugInf << std::setw(OWID) << " inteStep1";

    if (RK >= 1) {
      if (!negPrsDen) arrayUPrev = arrayU; 
      else { arrayU = arrayUPrev; goto label1;}
      arrayGodFluxStep2 = arrayGodFlux;
      inteStep2(ptcls); // arrayU updated
      debugInf << std::setw(OWID) << " inteStep2";      

      if (RK == 2) {
	if (!negPrsDen) arrayUPrev = arrayU; 
	else { arrayU = arrayUPrev; goto label1;}
	arrayGodFluxStep3 = arrayGodFlux;    
	inteStep3(ptcls); // arrayU updated
	debugInf << std::setw(OWID) << " inteStep3";
      }
    }
  }

  void Fluid::inteStep1(std::vector<Particle *> &ptcls) { 
    initGhostPoints();
    soundSpeed();
    calcTimeStep();
    enthalpy();
    rotateIJK(ptcls);

    // update conserved variables at the next time step
    for (std::size_t i = 1; i < gridNx - 1 ; ++i)
      for (std::size_t j = 1; j < gridNy - 1; ++j)
	for (std::size_t k = 1; k < gridNz - 1; ++k)
	  for (std::size_t m = 0; m < nInteg; ++m)
	    arrayU[i][j][k][m] -= (   timeStep / gridDx * (arrayGodFlux[i][j][k][m][0] - arrayGodFlux[i-1][j][k][m][0])
				    + timeStep / gridDy * (arrayGodFlux[i][j][k][m][1] - arrayGodFlux[i][j-1][k][m][1])
				    + timeStep / gridDz * (arrayGodFlux[i][j][k][m][2] - arrayGodFlux[i][j][k-1][m][2]) );
    UtoW();
  }
  
  void Fluid::inteStep2(std::vector<Particle *> &ptcls) { 
    initGhostPoints();
    soundSpeed();
    enthalpy();
    rotateIJK(ptcls);

    // update conserved variables at the next time step
    for (std::size_t i = 1; i < gridNx - 1 ; ++i)
      for (std::size_t j = 1; j < gridNy - 1; ++j)
	for (std::size_t k = 1; k < gridNz - 1; ++k)
	  for (std::size_t m = 0; m < nInteg; ++m)
	    arrayU[i][j][k][m] -= (   timeStep / (2*RK*gridDx) * (arrayGodFlux[i][j][k][m][0] - arrayGodFlux[i-1][j][k][m][0] 
						        + (arrayGodFluxStep2[i][j][k][m][0] - arrayGodFluxStep2[i-1][j][k][m][0]) )
				    + timeStep / (2*RK*gridDy) * (arrayGodFlux[i][j][k][m][1] - arrayGodFlux[i][j-1][k][m][1] 
						        + (arrayGodFluxStep2[i][j][k][m][1] - arrayGodFluxStep2[i][j-1][k][m][1]) )
				    + timeStep / (2*RK*gridDz) * (arrayGodFlux[i][j][k][m][2] - arrayGodFlux[i][j][k-1][m][2] 
						        + (arrayGodFluxStep2[i][j][k][m][2] - arrayGodFluxStep2[i][j][k-1][m][2])) );
    UtoW();
  }

  void Fluid::inteStep3(std::vector<Particle *> &ptcls) { 
    initGhostPoints();
    soundSpeed();
    enthalpy();
    rotateIJK(ptcls);

    // update conserved variables at the next time step
    for (std::size_t i = 1; i < gridNx - 1 ; ++i)
      for (std::size_t j = 1; j < gridNy - 1; ++j)
	for (std::size_t k = 1; k < gridNz - 1; ++k)
	  for (std::size_t m = 0; m < nInteg; ++m)
	    arrayU[i][j][k][m] -= (   timeStep / (6*gridDx) * (arrayGodFlux[i][j][k][m][0] - arrayGodFlux[i-1][j][k][m][0]
						         +(arrayGodFluxStep2[i][j][k][m][0] - arrayGodFluxStep2[i-1][j][k][m][0])
						      + 4*(arrayGodFluxStep3[i][j][k][m][0] - arrayGodFluxStep3[i-1][j][k][m][0]))
				      
				    + timeStep / (6*gridDy) * (arrayGodFlux[i][j][k][m][1] - arrayGodFlux[i][j-1][k][m][1]
						        + (arrayGodFluxStep2[i][j][k][m][1] - arrayGodFluxStep2[i][j-1][k][m][1])
						      + 4*(arrayGodFluxStep3[i][j][k][m][1] - arrayGodFluxStep3[i][j-1][k][m][1]))

				    + timeStep / (6*gridDz) * (arrayGodFlux[i][j][k][m][2] - arrayGodFlux[i][j][k-1][m][2]
						        + (arrayGodFluxStep2[i][j][k][m][2] - arrayGodFluxStep2[i][j][k-1][m][2])
						     + 4*( arrayGodFluxStep3[i][j][k][m][2] - arrayGodFluxStep3[i][j][k-1][m][2])) );
    UtoW();
  }

  void Fluid::rotateIJK(std::vector<Particle *> &ptcls) {
    std::size_t id[3][3] = {{0,1,2},{1,0,2},{2,1,0}};

    // for x, y, z sweeps: 
    // iDim=0 for x sweep
    // iDim=1 for y sweep
    // iDim=2 for z sweep
    for (std::size_t iDim = 0; iDim < nDim; ++iDim) {
      arrayURota = arrayU; // must have same rank and extent
      
      // switch components
      for (std::size_t i = 0; i < gridNx; ++i)
	for (std::size_t j = 0; j < gridNy; ++j)
	  for (std::size_t k = 0; k < gridNz; ++k) 
	    for (std::size_t jdim = 0; jdim < nDim; ++jdim) {
	      arrayURota[i][j][k][  varMom[jdim]  ] = arrayU[i][j][k][  varMom[id[iDim][jdim]]  ];
	      arrayURota[i][j][k][  varVel[jdim]  ] = arrayU[i][j][k][  varVel[id[iDim][jdim]]  ];
	    }

      flux(iDim, ptcls); // variables defined at cell centers

      // for local Riemann problem
      for (std::size_t i = 0; i < gridNx - 1; ++i) { // variables defined at cell faces
	for (std::size_t j = 0; j < gridNy - 1; ++j) {
	  for (std::size_t k = 0; k < gridNz -1; ++k) {
	    std::size_t IL[3] = {i, j, k};
	    std::size_t IR[3] = {i, j, k};
	    IR[iDim] += 1;
	    REAL UL[9], UR[9], FL[5], FR[5], HL, HR; // local variable only
	    HL = arrayH[IL[0]] [IL[1]] [IL[2]];
	    HR = arrayH[IR[0]] [IR[1]] [IR[2]];
	    for (std::size_t m = 0; m < nVar; ++m) {
	      UL[m] = arrayURota[IL[0]] [IL[1]] [IL[2]] [m];
	      UR[m] = arrayURota[IR[0]] [IR[1]] [IR[2]] [m];
	    }	
	    for (std::size_t m = 0; m < nInteg; ++m) {
	      FL[m] = arrayFlux[IL[0]] [IL[1]] [IL[2]] [m];
	      FR[m] = arrayFlux[IR[0]] [IR[1]] [IR[2]] [m];
	    }    

	    REAL avgH   = (sqrt(UL[varDen])*HL + sqrt(UR[varDen])*HR)/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
	    REAL avgU   = (sqrt(UL[varDen])*UL[varVel[0]] + sqrt(UR[varDen])*UR[varVel[0]])/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
	    REAL avgV   = (sqrt(UL[varDen])*UL[varVel[1]] + sqrt(UR[varDen])*UR[varVel[1]])/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
	    REAL avgW   = (sqrt(UL[varDen])*UL[varVel[2]] + sqrt(UR[varDen])*UR[varVel[2]])/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
	    REAL avgh   = avgH - 0.5*(avgU*avgU + avgV*avgV + avgW*avgW); // static specific enthalpy
	    REAL pStar;
	    guessPressure(UL, UR, pStar);

	    int solver = static_cast<int> (dem::Parameter::getSingleton().parameter["solver"]);
	    switch (solver) {
	    case -1:
	      // exactSolver itself can provide exact solution to the full domain, if run separately.
	      // herein exactSolver is used at each discritized face to give exact solution to local Riemann problem, for the purpose 
	      // of comparing with other solvers.
	      exactSolver(UL, UR, 0, iDim, i, j, k); // 0 = x/t = s
	      break;
	    case 0: 
	      RoeSolver(UL, UR, FL, FR, HL, HR, iDim, i, j, k);
	      break;
	    case 1:
	      HllcSolver(UL, UR, FL, FR, HL, HR, iDim, i, j, k);
	      break;
	    case 2:
	      HlleSolver(UL, UR, FL, FR, HL, HR, iDim, i, j, k);
	      break;
	    case 3:
	      if (pStar <= UL[varPrs] || pStar <= UR[varPrs] || avgh <= 0 || negPrsDen) // first 2 expressions cover rarefaction and contact waves
		HllcSolver(UL, UR, FL, FR, HL, HR, iDim, i, j, k);
	      else // use Roe solver for shock waves
		RoeSolver(UL, UR, FL, FR, HL, HR, iDim, i, j, k); 
	      break;
	    case 4:
	      if (pStar <= UL[varPrs] || pStar <= UR[varPrs] || avgh <= 0 || negPrsDen) // first 2 expressions cover rarefaction and contact waves
		HlleSolver(UL, UR, FL, FR, HL, HR, iDim, i, j, k);
	      else // use Roe solver for shock waves
		RoeSolver(UL, UR, FL, FR, HL, HR, iDim, i, j, k);
	      break;
	    case 5:
	      if (pStar <= UL[varPrs] || pStar <= UR[varPrs] || avgh <= 0 || negPrsDen) // first 2 expressions cover rarefaction and contact waves
		exactSolver(UL, UR, 0, iDim, i, j, k);
	      else // use Roe solver for shock waves
		RoeSolver(UL, UR, FL, FR, HL, HR, iDim, i, j, k);
	      break;
	    }

	  }
	}
      }

      for (std::size_t i = 0; i < gridNx -1; ++i)
	for (std::size_t j = 0; j < gridNy -1; ++j)
	  for (std::size_t k = 0; k < gridNz -1; ++k)
	    for (std::size_t m = 0; m < nInteg; ++m)
	      arrayGodFluxTmp[i][j][k][m] = arrayGodFlux[i][j][k][m][iDim];

      // switch components back for consistency with u
      for (std::size_t i = 0; i < gridNx - 1; ++i)
	for (std::size_t j = 0; j < gridNy - 1; ++j)
	  for (std::size_t k = 0; k < gridNz -1; ++k)
	    for (std::size_t m = 0; m < nDim; ++m)
	      arrayGodFlux[i][j][k][varMom[m]][iDim] = arrayGodFluxTmp[i][j][k][ varMom[id[iDim][m]] ];

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
	bool inGrid = (i > 0 && i < gridNx-1 && j > 0 && j < gridNy-1 && k > 0 && k < gridNz-1 );
	REAL momQuot[3], denQuot[3], u0[3];
	if (inGrid) {
	  momQuot[0] = (arrayU[i+1][j][k][varMom[0]] - arrayU[i-1][j][k][varMom[0]]) / (2*gridDx);
	  momQuot[1] = (arrayU[i][j+1][k][varMom[1]] - arrayU[i][j-1][k][varMom[1]]) / (2*gridDy);
	  momQuot[2] = (arrayU[i][j][k+1][varMom[2]] - arrayU[i][j][k-1][varMom[2]]) / (2*gridDz);

	  denQuot[0] = (arrayU[i+1][j][k][varDen] - arrayU[i-1][j][k][varDen]) / (2*gridDx);
	  denQuot[1] = (arrayU[i][j+1][k][varDen] - arrayU[i][j-1][k][varDen]) / (2*gridDy);
	  denQuot[2] = (arrayU[i][j][k+1][varDen] - arrayU[i][j][k-1][varDen]) / (2*gridDz);

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
	  // a. momentum penalization
	  arrayU[i][j][k][varMom[m]] -= arrayU[i][j][k][varMsk] * arrayPenalForce[i][j][k][m] * timeStep * (Cdi/Cd+1) ;
	  // b. influence of momentum penalization on energy
	  arrayU[i][j][k][varEng]    -= arrayU[i][j][k][varMsk] * arrayPenalForce[i][j][k][m] * arrayU[i][j][k][varVel[m]] * timeStep * (Cdi/Cd+1);
	}

	// 2. porosity correction
	//   a. porosity corrections are incorporated by changes in flux()

	//   b. influence of porosity correction (1-1.0/porosity)*momQuot[m] and particle velocity term u0[m]/porosity*denQuot[m] on energy
	if (inGrid) {
	  for (std::size_t m = 0; m < nDim; ++m)
	    arrayU[i][j][k][varEng]  += arrayU[i][j][k][varMsk] * (-0.5*pow(arrayU[i][j][k][varVel[m]],2))
                                         * ( (1-1.0/porosity)*momQuot[m] + u0[m]/porosity*denQuot[m] ) * timeStep;
	}

      }
    }
  }

  /*
  // 1st implementation: lower efficiency
  for (std::size_t i = 0; i < gridNx; ++i)
    for (std::size_t j = 0; j < gridNy; ++j)
      for (std::size_t k = 0; k < gridNz; ++k) {

	// calculate momentum quotient before modification
	bool inGrid = (i > 0 && i < gridNx-1 && j > 0 && j < gridNy-1 && k > 0 && k < gridNz-1 );
	REAL momQuot[3];
	if (inGrid) {
	  momQuot[0] = (arrayU[i+1][j][k][varMom[0]] - arrayU[i-1][j][k][varMom[0]]) / (2*gridDx);
	  momQuot[1] = (arrayU[i][j+1][k][varMom[1]] - arrayU[i][j-1][k][varMom[1]]) / (2*gridDy);
	  momQuot[2] = (arrayU[i][j][k+1][varMom[2]] - arrayU[i][j][k-1][varMom[2]]) / (2*gridDz);
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
  
  void Fluid::initGhostPoints() {
    // non-reflecting BCs
    for (std::size_t j = 1; j < gridNy - 1; ++j)
      for (std::size_t k = 1; k < gridNz - 1; ++k)
	for (std::size_t m = 0; m < nVar; ++m) {
	  arrayU[0][j][k][m]    = arrayU[1][j][k][m]; 
	  arrayU[gridNx-1][j][k][m] = arrayU[gridNx-2][j][k][m]; 
	}

    for (std::size_t i = 1; i < gridNx - 1; ++i)
      for (std::size_t k = 1; k < gridNz -1; ++k)
	for (std::size_t m = 0; m < nVar; ++m) {
	  arrayU[i][0][k][m]    = arrayU[i][1][k][m]; 
	  arrayU[i][gridNy-1][k][m] = arrayU[i][gridNy-2][k][m]; 
	}

    for (std::size_t i = 1; i < gridNx - 1; ++i)
      for (std::size_t j = 1; j < gridNy - 1; ++j)
	for (std::size_t m = 0; m < nVar; ++m) {
	  arrayU[i][j][0][m]    = arrayU[i][j][1][m]; 
	  arrayU[i][j][gridNz-1][m] = arrayU[i][j][gridNz-2][m]; 
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
      for (std::size_t j = 1; j < gridNy - 1; ++j)
	for (std::size_t k = 1; k < gridNz - 1; ++k)
	  for (std::size_t m = 0; m < 1; ++m) { // x-direction
	    arrayU[0][j][k][varMom[m]]    *= (1-2*arrayBC[0]); 
	    arrayU[gridNx-1][j][k][varMom[m]] *= (1-2*arrayBC[1]); 
	    arrayU[0][j][k][varVel[m]]    *= (1-2*arrayBC[0]); 
	    arrayU[gridNx-1][j][k][varVel[m]] *= (1-2*arrayBC[1]); 
	  }

      for (std::size_t i = 1; i < gridNx - 1; ++i)
	for (std::size_t k = 1; k < gridNz - 1; ++k)
	  for (std::size_t m = 1; m < 2; ++m) { // y-direction
	    arrayU[i][0][k][varMom[m]]    *= (1-2*arrayBC[2]); 
	    arrayU[i][gridNy-1][k][varMom[m]] *= (1-2*arrayBC[3]);
	    arrayU[i][0][k][varVel[m]]    *= (1-2*arrayBC[2]); 
	    arrayU[i][gridNy-1][k][varVel[m]] *= (1-2*arrayBC[3]);  
	  }

      for (std::size_t i = 1; i < gridNx - 1; ++i)
	for (std::size_t j = 1; j < gridNy - 1; ++j)
	  for (std::size_t m = 2; m < 3; ++m) { // z-direction
	    arrayU[i][j][0][varMom[m]]    *= (1-2*arrayBC[4]); 
	    arrayU[i][j][gridNz-1][varMom[m]] *= (1-2*arrayBC[5]); 
	    arrayU[i][j][0][varVel[m]]    *= (1-2*arrayBC[4]); 
	    arrayU[i][j][gridNz-1][varVel[m]] *= (1-2*arrayBC[5]); 
	  }
    }  
  }
  
  void Fluid::calcTimeStep() {
    std::valarray<REAL> gridX(gridNx * gridNy * gridNz);
    std::valarray<REAL> gridY(gridNx * gridNy * gridNz);
    std::valarray<REAL> gridZ(gridNx * gridNy * gridNz);

    std::valarray<REAL> velZ(gridNx * gridNy * gridNz);
    std::valarray<REAL> sound(gridNx * gridNy * gridNz);

    for (std::size_t i = 0; i < gridNx ; ++i)
      for (std::size_t j = 0; j < gridNy; ++j)
	for (std::size_t k = 0; k < gridNz; ++k) {
	  gridX[i + j * gridNx + k * gridNx * gridNy] = fabs(arrayU[i][j][k][varVel[0]]) + arraySoundSpeed[i][j][k];
	  gridY[i + j * gridNx + k * gridNx * gridNy] = fabs(arrayU[i][j][k][varVel[1]]) + arraySoundSpeed[i][j][k];
	  gridZ[i + j * gridNx + k * gridNx * gridNy] = fabs(arrayU[i][j][k][varVel[2]]) + arraySoundSpeed[i][j][k];

	  velZ[i + j * gridNx + k * gridNx * gridNy]  = arrayU[i][j][k][varVel[2]];
	  sound[i + j * gridNx + k * gridNx * gridNy] = arraySoundSpeed[i][j][k];
	}

    std::valarray<REAL> dtMin(3);
    dtMin[0] = gridDx / gridX.max();
    dtMin[1] = gridDy / gridY.max();
    dtMin[2] = gridDz / gridZ.max();
    
    timeStep = std::min(timeStep, CFL * dtMin.min());
    timeAccrued += timeStep;
    debugInf << std::endl
	     << std::setw(OWID) << iteration 
	     << std::setw(OWID) << timeStep 
	     << std::setw(OWID) << timeAccrued /*
	     << std::setw(OWID) << velZ.max()
	     << std::setw(OWID) << velZ.min()
	     << std::setw(OWID) << sound.max()
	     << std::setw(OWID) << gridZ.max() */;	    
  }

  void Fluid::soundSpeed() {
    for (std::size_t i = 0; i < gridNx ; ++i)
      for (std::size_t j = 0; j < gridNy; ++j)
	for (std::size_t k = 0; k < gridNz; ++k)
	  arraySoundSpeed[i][j][k] = sqrt(gama * arrayU[i][j][k][varPrs] / arrayU[i][j][k][varDen]);
  }

  // total specific enthalphy, NOT static specific enthalpy
  void Fluid::enthalpy() {
    for (std::size_t i = 0; i < gridNx ; ++i)
      for (std::size_t j = 0; j < gridNy; ++j)
	for (std::size_t k = 0; k < gridNz; ++k)
	  arrayH[i][j][k] = (arrayU[i][j][k][varEng] + arrayU[i][j][k][varPrs]) / arrayU[i][j][k][varDen];
  }

  void Fluid::initialCondition() {
    if (leftType == 1 || leftType == 2) { // normal shock w/ and w/o Rankine-Hugoniot conditions
      for (std::size_t i = 0; i < gridNx; ++i)
	for (std::size_t j = 0; j < gridNy; ++j)
	  for (std::size_t k = 0; k < gridNz; ++k) {
	    if (arrayGridCoord[i][j][k][2] <= z2L) {
	      arrayU[i][j][k][varDen] = rhoL;
	      arrayU[i][j][k][varPrs] = pL;
	      arrayU[i][j][k][varVel[0]] = 0;
	      arrayU[i][j][k][varVel[1]] = 0;
	      arrayU[i][j][k][varVel[2]] = uL;
	      if (arrayGridCoord[i][j][k][2] == z2L) arrayU[i][j][k][varVel[2]] = 0; // for 123 problem
	    } else if (arrayGridCoord[i][j][k][2] > z2L) {
	      arrayU[i][j][k][varDen] = rhoR;
	      arrayU[i][j][k][varPrs] = pR;
	      arrayU[i][j][k][varVel[0]] = 0;
	      arrayU[i][j][k][varVel[1]] = 0;
	      arrayU[i][j][k][varVel[2]] = uR;
	    } 
	  }
    }
    else if (leftType == 3) { // normal shock in x, y, z directions
      for (std::size_t i = 0; i < gridNx; ++i)
	for (std::size_t j = 0; j < gridNy; ++j)
	  for (std::size_t k = 0; k < gridNz; ++k) {
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
      for (std::size_t i = 0; i < gridNx; ++i)
	for (std::size_t j = 0; j < gridNy; ++j)
	  for (std::size_t k = 0; k < gridNz; ++k) {
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
      for (std::size_t i = 0; i < gridNx; ++i)
	for (std::size_t j = 0; j < gridNy; ++j)
	  for (std::size_t k = 0; k < gridNz; ++k) {
	    if ( arrayGridCoord[i][j][k][2] >= z1L && arrayGridCoord[i][j][k][2] <= z2L ) {
	      arrayU[i][j][k][varDen] = rhoL;
	      arrayU[i][j][k][varPrs] = pL;
	      arrayU[i][j][k][varVel[0]] = 0;
	      arrayU[i][j][k][varVel[1]] = 0;
	      arrayU[i][j][k][varVel[2]] = uL;
	    } else if ( arrayGridCoord[i][j][k][2] > z2L) {
	      arrayU[i][j][k][varDen] = rhoR;
	      arrayU[i][j][k][varPrs] = pR;
	      arrayU[i][j][k][varVel[0]] = 0;
	      arrayU[i][j][k][varVel[1]] = 0;
	      arrayU[i][j][k][varVel[2]] = uR;
	    } else if ( arrayGridCoord[i][j][k][2] < z1L) {
	      arrayU[i][j][k][varDen] = rhoBL;
	      arrayU[i][j][k][varPrs] = pBL;
	      arrayU[i][j][k][varVel[0]] = 0;
	      arrayU[i][j][k][varVel[1]] = 0;
	      arrayU[i][j][k][varVel[2]] = uBL;
	    }
	  }
    }

    WtoU();
  }

  void Fluid::RankineHugoniot() { // Rankine-Hugoniot conditions
    if (leftType == 1) {
      shockSpeed = MachShock*sqrt(gama*pR/rhoR);
      rhoL = ( pow(rhoR*(shockSpeed-uR),2)*(1+gama) ) / ( rhoR*pow(shockSpeed-uR,2)*(gama-1) + 2*pR*gama);
      pL   = (pR*(1-gama)+2*rhoR*pow(shockSpeed-uR,2)) / (1+gama);
      uL   = ( rhoR*(shockSpeed-uR)*(2*shockSpeed + uR*(gama-1)) - 2*pR*gama ) / (rhoR * (shockSpeed-uR) * (1+gama));
      debugInf << std::setw(OWID) << "shockSpeed" << std::setw(OWID) << shockSpeed << std::endl;
      debugInf << std::setw(OWID) << "rhoL" << std::setw(OWID) << rhoL << std::endl;
      debugInf << std::setw(OWID) << "pL" << std::setw(OWID) << pL << std::endl;
      debugInf << std::setw(OWID) << "uL" << std::setw(OWID) << uL << std::endl;
    } 

    MachL = uL / sqrt(gama*pL/rhoL);
    debugInf << std::setw(OWID) << "MachL" << std::setw(OWID) << MachL << std::endl << std::endl;
  }

  void Fluid::flux(std::size_t iDim, std::vector<Particle *> &ptcls) {
    for (std::size_t i = 0; i < gridNx; ++i)
      for (std::size_t j = 0; j < gridNy; ++j)
	for (std::size_t k = 0; k < gridNz; ++k) {
	  arrayFlux[i][j][k][varDen]    = arrayURota[i][j][k][varDen] * arrayURota[i][j][k][varVel[0]]; // rho*u
	  arrayFlux[i][j][k][varMom[0]] = arrayURota[i][j][k][varDen] * pow(arrayURota[i][j][k][varVel[0]],2) + arrayURota[i][j][k][varPrs]; // rho*u^2 + p
	  arrayFlux[i][j][k][varMom[1]] = arrayURota[i][j][k][varDen] * arrayURota[i][j][k][varVel[0]] * arrayURota[i][j][k][varVel[1]]; // rho*u*v
	  arrayFlux[i][j][k][varMom[2]] = arrayURota[i][j][k][varDen] * arrayURota[i][j][k][varVel[0]] * arrayURota[i][j][k][varVel[2]]; // rho*u*w
	  arrayFlux[i][j][k][varEng]    = arrayURota[i][j][k][varVel[0]] * (arrayURota[i][j][k][varEng] + arrayURota[i][j][k][varPrs]);  // u*(E + p)
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

	// all 5 equations are modified in terms of porosity and Darcy's velocity for high-porosity material
	// continuity equation modified on porosity
	arrayFlux[i][j][k][varDen]    = 1.0/porosity * arrayURota[i][j][k][varDen] * (arrayURota[i][j][k][varVel[0]] - u0[iDim]) ; // rho*(u-u0) / porosity

	// momentum equations modified on porosity
	arrayFlux[i][j][k][varMom[0]] = 1.0/porosity * arrayURota[i][j][k][varDen] * pow(arrayURota[i][j][k][varVel[0]],2) + arrayURota[i][j][k][varPrs]; // rho*u^2 / porosity + p*porosity
	arrayFlux[i][j][k][varMom[1]] = 1.0/porosity * arrayURota[i][j][k][varDen] * arrayURota[i][j][k][varVel[0]] * arrayURota[i][j][k][varVel[1]]; // (rho*u)*v / porosity
	arrayFlux[i][j][k][varMom[2]] = 1.0/porosity * arrayURota[i][j][k][varDen] * arrayURota[i][j][k][varVel[0]] * arrayURota[i][j][k][varVel[2]]; // (rho*u)*w / porosity

	// energy equation modified on porosity
	//arrayFlux[i][j][k][varEng]    = 1.0/porosity * arrayURota[i][j][k][varVel[0]] * (arrayURota[i][j][k][varEng] + arrayURota[i][j][k][varPrs]);  // u*(E + p) / porosity
      }
    }
  }

  void Fluid::exactSolver(REAL UL[], REAL UR[], REAL relaCoord, std::size_t iDim, std::size_t it, std::size_t jt, std::size_t kt) {
    REAL pStar, uStar;
    findPrsVel(UL, UR, pStar, uStar);

    REAL s = relaCoord/timeAccrued;
    REAL d, u, p;
    sampling(UL, UR, pStar, uStar, s, d, u, p);
    //if (iDim==2 && it==3 && jt==3 && kt==250) debugInf << " relaCoord=" << relaCoord << " pstar=" << pStar << " uStar=" << uStar << " den=" << d << " u=" << u << " p=" << p; 
    arrayGodFlux[it][jt][kt][varDen][iDim] = d*u;
    arrayGodFlux[it][jt][kt][varMom[0]][iDim] = d*u*u+p;
    REAL v = uStar >= s ? UL[varVel[1]] : UR[varVel[1]]; // Equ. (4.115) of Toro
    REAL w = uStar >= s ? UL[varVel[2]] : UR[varVel[2]];
    arrayGodFlux[it][jt][kt][varMom[1]][iDim] = d*u*v;
    arrayGodFlux[it][jt][kt][varMom[2]][iDim] = d*u*w;
    arrayGodFlux[it][jt][kt][varEng][iDim] = u*(p*gama/(gama-1)+0.5*d*(u*u+v*v+w*w)); // u*(E+p)
  }

  void Fluid::RoeSolver(REAL UL[], REAL UR[], REAL FL[], REAL FR[], REAL HL, REAL HR, std::size_t iDim, std::size_t it, std::size_t jt, std::size_t kt) {
    // it, jt, kt defined at cell faces
    REAL avgRho =  sqrt(UL[varDen]*UR[varDen]);
    REAL avgH   = (sqrt(UL[varDen])*HL + sqrt(UR[varDen])*HR)/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
    REAL avgU   = (sqrt(UL[varDen])*UL[varVel[0]] + sqrt(UR[varDen])*UR[varVel[0]])/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
    REAL avgV   = (sqrt(UL[varDen])*UL[varVel[1]] + sqrt(UR[varDen])*UR[varVel[1]])/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
    REAL avgW   = (sqrt(UL[varDen])*UL[varVel[2]] + sqrt(UR[varDen])*UR[varVel[2]])/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
    REAL avgh   = avgH - 0.5*(avgU*avgU + avgV*avgV + avgW*avgW); // static specific enthalpy 
    REAL avgSoundSpeed = sqrt((gama-1)*avgh);
    
    REAL eigen[5];
    eigen[varDen]    = avgU - avgSoundSpeed;
    eigen[varMom[0]] = eigen[varMom[1]] = eigen[varMom[2]] = avgU;
    eigen[varEng]    = avgU + avgSoundSpeed;

    REAL avgWaveStr[5], du[9];
    for (std::size_t i = 0; i < nVar; ++i)
      du[i] = UR[i] - UL[i];

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

    // entropy fix: use HartenHyman Entropy Fix and Roe-Averaged States
    // for left transonic rarefaction
    REAL denStar = UL[varDen] + avgWaveStr[varDen];
    REAL uStar = ( UL[varDen]*UL[varVel[0]] + avgWaveStr[varDen]*(avgU-avgSoundSpeed) ) / denStar;
    REAL pStar = (gama-1)* (UL[varEng] + avgWaveStr[varDen]*(avgH-avgU*avgSoundSpeed) - 0.5*denStar*uStar*uStar );
    REAL aStarL= sqrt(gama*pStar/denStar);
    REAL aL = sqrt(gama*UL[varPrs]/UL[varDen]);
    REAL eigenL = UL[varVel[0]] - aL;
    REAL eigenR = uStar - aStarL;
    if (eigenL < 0 && eigenR > 0)
      eigen[varDen] = eigenL * (eigenR-eigen[varDen])/(eigenR-eigenL);
    // for right transonic rarefaction
    denStar = UR[varDen] - avgWaveStr[varEng];
    uStar = ( UR[varDen]*UR[varVel[0]] - avgWaveStr[varEng]*(avgU+avgSoundSpeed) ) / denStar;
    pStar = (gama-1)* (UR[varEng] - avgWaveStr[varEng]*(avgH+avgU*avgSoundSpeed) - 0.5*denStar*uStar*uStar );
    REAL aStarR = sqrt(gama*pStar/denStar);
    REAL aR = sqrt(gama*UR[varPrs]/UR[varDen]);
    eigenL = uStar + aStarR;
    eigenR = UR[varVel[0]] + aR;
    if (eigenL < 0 && eigenR > 0)
      eigen[varEng] = eigenR * (eigen[varEng]-eigenL)/(eigenR-eigenL);
    // end of entropy fix

    /*
    // entropy fix: use HartenHyman Entropy Fix and PrimitiveVariable Riemann Solver (PVRS)
    REAL aBar = 0.5*(sqrt(gama*UL[varPrs]/UL[varDen])+sqrt(gama*UR[varPrs]/UR[varDen]));
    REAL denBar = 0.5*(UL[varDen]+UR[varDen]);
    REAL pStar = 0.5*(UL[varPrs]+UR[varPrs]) + 0.5*(UL[varVel[0]]-UR[varVel[0]])/(denBar*aBar);
    REAL uStar = 0.5*(UL[varVel[0]]+UR[varVel[0]]) + 0.5*(UL[varPrs]-UR[varPrs])/(denBar*aBar);
    REAL denStarL = UL[varDen] + (UL[varVel[0]]-uStar)*denBar/aBar;
    REAL denStarR = UR[varDen] + (uStar-UR[varVel[0]])*denBar/aBar;

    // for left transonic rarefaction
    REAL eigenL = UL[varVel[0]] - sqrt(gama*UL[varPrs]/UL[varDen]);
    REAL eigenR = uStar - sqrt(gama*pStar/denStarL);
    if (eigenL < 0 && eigenR > 0)
      eigen[varDen] = eigenL * (eigenR-eigen[varDen])/(eigenR-eigenL);
    // for right transonic rarefaction
    eigenL = uStar + sqrt(gama*pStar/denStarR);
    eigenR = UR[varVel[0]] + sqrt(gama*UR[varPrs]/UR[varDen]);
    if (eigenL < 0 && eigenR > 0)
      eigen[varEng] = eigenR * (eigen[varEng]-eigenL)/(eigenR-eigenL);
    // end of entropy fix
    */

    /*
    // entropy fix: use HartenHyman Entropy Fix and TwoRarefaction Riemann Solver (TRRS)
    // it is supposed to work for 123 problem, but it actually does not.
    REAL z = 0.5*(gama-1)/gama;
    REAL aL = sqrt(gama*UL[varPrs]/UL[varDen]);
    REAL aR = sqrt(gama*UR[varPrs]/UR[varDen]);
    REAL pStar = pow((aL+aR-0.5*(gama-1)*(UR[varVel[0]]-UL[varVel[0]]))/(aL/pow(UL[varPrs],z)+aR/pow(UR[varPrs],z)), 1/z);
    // for left transonic rarefaction
    REAL aStarL = aL*pow(pStar/UL[varPrs],z);
    REAL uStar = UL[varVel[0]] + 2/(gama-1)*(aL-aStarL);
    REAL eigenL = UL[varVel[0]] - aL;
    REAL eigenR = uStar - aStarL;
    if (eigenL < 0 && eigenR > 0)
      eigen[varDen] = eigenL * (eigenR-eigen[varDen])/(eigenR-eigenL);
    // for right transonic rarefaction
    REAL aStarR = aR*pow(pStar/UR[varPrs],z);
    uStar = UR[varVel[0]] + 2/(gama-1)*(aStarR-aR);
    eigenL = uStar + aStarR;
    eigenR = UR[varVel[0]] + aR;
    if (eigenL < 0 && eigenR > 0)
      eigen[varEng] = eigenR * (eigen[varEng]-eigenL)/(eigenR-eigenL);
    // end of entropy fix
    */

    REAL RF[5];
    for (std::size_t i = 0; i < nInteg; ++i)
      RF[i] = 0.5*(FL[i] + FR[i]);

    for (std::size_t ie = 0; ie < nInteg; ++ie) {
      for (std::size_t je = 0; je < nInteg; ++je) {
	RF[ie] -= 0.5 * avgWaveStr[je] * fabs(eigen[je]) * avgK[ie][je];
      }
      arrayGodFlux[it][jt][kt][ie][iDim] = RF[ie];
    }
  }

  void Fluid::HlleSolver(REAL UL[], REAL UR[], REAL FL[], REAL FR[], REAL HL, REAL HR, std::size_t iDim, std::size_t it, std::size_t jt, std::size_t kt) {
    REAL avgH   = (sqrt(UL[varDen])*HL + sqrt(UR[varDen])*HR)/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
    REAL avgU   = (sqrt(UL[varDen])*UL[varVel[0]] + sqrt(UR[varDen])*UR[varVel[0]])/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
    REAL avgV   = (sqrt(UL[varDen])*UL[varVel[1]] + sqrt(UR[varDen])*UR[varVel[1]])/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
    REAL avgW   = (sqrt(UL[varDen])*UL[varVel[2]] + sqrt(UR[varDen])*UR[varVel[2]])/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
    REAL avgh   = avgH - 0.5*(avgU*avgU + avgV*avgV + avgW*avgW); // static specific enthalpy 
    REAL avgSoundSpeed = sqrt((gama-1)*avgh);

    REAL aL = sqrt(gama*UL[varPrs]/UL[varDen]);
    REAL aR = sqrt(gama*UR[varPrs]/UR[varDen]);
    REAL SL = std::min(avgU-avgSoundSpeed, UL[varVel[0]]-aL); // this is the core of HLLE, positively conservative
    REAL SR = std::max(avgU+avgSoundSpeed, UR[varVel[0]]+aR);
    REAL uBar = (SR+SL)/2;
    REAL cBar = (SR-SL)/2;

    /*
    REAL dBar = (sqrt(UL[varDen])*aL*aL+sqrt(UR[varDen])*aR*aR)/(sqrt(UL[varDen])+sqrt(UR[varDen])) 
                + (0.5*sqrt(UL[varDen]*UR[varDen])/pow(sqrt(UL[varDen])+sqrt(UR[varDen]),2))*pow(UR[varVel[0]]-UL[varVel[0]],2);
    REAL SL = avgU - dBar; // the fastest signal velocities, avgU is uBar in Toro's book
    REAL SR = avgU + dBar;

    REAL avgR[5]; // eigenvector for contact wave, not obtained from Roe's Jacobian
    avgR[varDen]    = 1;
    avgR[varMom[0]] = uBar;
    avgR[varMom[1]] = 0;
    avgR[varMom[2]] = 0;
    avgR[varEng]    = 0.5 * (uBar*uBar); // what about 3D?

    /*
    REAL l[5] = {1-(gama-1)/2*uBar*uBar/(cBar*cBar), (gama-1)*uBar/(cBar*cBar), 0, 0, (1-gama)/(cBar*cBar)}; // what about 3D?
    REAL eta = 0;
    for (std::size_t i = 0; i < nInteg; ++i)
      eta += l[i] * (UR[i] - UL[i]);
    REAL bNeg = std::min(SL, .0);
    REAL bPos = std::max(SR, .0);
    REAL delta = 1/(timeStep*(cBar+fabs(uBar)));
    */
    if (SL >= 0) {
      for (std::size_t ie = 0; ie < nInteg; ++ie)
	arrayGodFlux[it][jt][kt][ie][iDim] = FL[ie];
    } else if (SL < 0 && SR > 0) {
      for (std::size_t ie = 0; ie < nInteg; ++ie)
	arrayGodFlux[it][jt][kt][ie][iDim] = (SR*FL[ie]-SL*FR[ie]+SL*SR*(UR[ie]-UL[ie]))/(SR-SL);      
    }else if (SR <= 0) {
      for (std::size_t ie = 0; ie < nInteg; ++ie)
	arrayGodFlux[it][jt][kt][ie][iDim] = FR[ie];
    }
    /*
    // modification to HLL
    for (std::size_t ie = 0; ie < nInteg; ++ie)
      arrayGodFlux[it][jt][kt][ie][iDim] -= bPos*bNeg*timeStep/2*delta*eta*avgR[ie];
    */
  }

  void Fluid::HllcSolver(REAL UL[], REAL UR[], REAL FL[], REAL FR[], REAL HL, REAL HR, std::size_t iDim, std::size_t it, std::size_t jt, std::size_t kt) {
    REAL avgU = (sqrt(UL[varDen])*UL[varVel[0]] + sqrt(UR[varDen])*UR[varVel[0]])/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
    REAL aL = sqrt(gama*UL[varPrs]/UL[varDen]);
    REAL aR = sqrt(gama*UR[varPrs]/UR[varDen]);

    ///*
    // HLLC part
    // Pressure-based estimate by Toro
    REAL pStar=std::max(0.0, 0.5*(UL[varPrs]+UR[varPrs])-0.5*(UR[varVel[0]]-UL[varVel[0]])*0.25*(UL[varDen]+UR[varDen])*(aL+aR));

    // Pressure-based estimate by TRRS
    //REAL z = 0.5*(gama-1)/gama;
    //REAL pStar = pow((aL+aR-0.5*(gama-1)*(UR[varVel[0]]-UL[varVel[0]]))/(aL/pow(UL[varPrs],z)+aR/pow(UR[varPrs],z)), 1/z);

    REAL qL, qR;
    if (pStar <= UL[varPrs])
      qL = 1;
    else 
      qL = sqrt(1+(gama+1)/(2*gama)*(pStar/UL[varPrs]-1));
    if (pStar <= UR[varPrs])
      qR = 1;
    else 
      qR = sqrt(1+(gama+1)/(2*gama)*(pStar/UR[varPrs]-1));
    REAL SL = UL[varVel[0]] - aL*qL; // the fastest signal velocities
    REAL SR = UR[varVel[0]] + aR*qR;
    // end of HLLC part
    //*/

    /*
    // Hlle part, this part does not work for 123 problem.
    REAL dBar = (sqrt(UL[varDen])*aL*aL+sqrt(UR[varDen])*aR*aR)/(sqrt(UL[varDen])+sqrt(UR[varDen])) 
                + (0.5*sqrt(UL[varDen]*UR[varDen])/pow(sqrt(UL[varDen])+sqrt(UR[varDen]),2))*pow(UR[varVel[0]]-UL[varVel[0]],2);
    REAL SL = avgU - dBar; // the fastest signal velocities, avgU is uBar in Toro's book
    REAL SR = avgU + dBar;
    // end of Hlle part
    */

    REAL SStar = (UR[varPrs]-UL[varPrs]+UL[varDen]*UL[varVel[0]]*(SL-UL[varVel[0]])-UR[varDen]*UR[varVel[0]]*(SR-UR[varVel[0]])) / 
                 (UL[varDen]*(SL-UL[varVel[0]])-UR[varDen]*(SR-UR[varVel[0]]));
    REAL DStar[5] = {0, 1, 0, 0, SStar};

    if (SL >= 0) {
      for (std::size_t ie = 0; ie < nInteg; ++ie)
	arrayGodFlux[it][jt][kt][ie][iDim] = FL[ie];
    } else if (SL < 0 && SStar >= 0) {
      REAL uStarL[5];
      uStarL[0] = 1;
      uStarL[1] = SStar;
      uStarL[2] = UL[varVel[1]];
      uStarL[3] = UL[varVel[2]];
      uStarL[4] = UL[varEng]/UL[varDen] + (SStar-UL[varVel[0]])*(SStar+UL[varPrs]/UL[varDen]/(SL-UL[varVel[0]]));
      for (std::size_t i = 0; i < nInteg; ++i)
	uStarL[i] *= UL[varDen]*(SL-UL[varVel[0]])/(SL-SStar);

      //REAL pLR = 0.5*(UL[varPrs]+UR[varPrs] + UL[varDen]*(SL-UL[varVel[0]])*(SStar-UL[varVel[0]]) + UR[varDen]*(SR-UR[varVel[0]])*(SStar-UR[varVel[0]]));
      for (std::size_t ie = 0; ie < nInteg; ++ie) {
	arrayGodFlux[it][jt][kt][ie][iDim] = FL[ie] + SL*(uStarL[ie]-UL[ie]);
	// variant 1:
	//arrayGodFlux[it][jt][kt][ie][iDim] = (SStar*(SL*UL[ie]-FL[ie]) + SL*(UL[varPrs]+UL[varDen]*(SL-UL[varVel[0]])*(SStar-UL[varVel[0]]))*DStar[ie]) / (SL-SStar);
	// variant 2:
	//arrayGodFlux[it][jt][kt][ie][iDim] = (SStar*(SL*UL[ie]-FL[ie]) + SL*pLR*DStar[ie]) / (SL-SStar);	
      }
    } else if (SStar < 0 && SR > 0) {
      REAL uStarR[5];
      uStarR[0] = 1;
      uStarR[1] = SStar;
      uStarR[2] = UR[varVel[1]];
      uStarR[3] = UR[varVel[2]];
      uStarR[4] = UR[varEng]/UR[varDen] + (SStar-UR[varVel[0]])*(SStar+UR[varPrs]/UR[varDen]/(SR-UR[varVel[0]]));
      for (std::size_t i = 0; i < nInteg; ++i)
	uStarR[i] *= UR[varDen]*(SR-UR[varVel[0]])/(SR-SStar);

      //REAL pLR = 0.5*(UL[varPrs]+UR[varPrs] + UL[varDen]*(SL-UL[varVel[0]])*(SStar-UL[varVel[0]]) + UR[varDen]*(SR-UR[varVel[0]])*(SStar-UR[varVel[0]]));
      for (std::size_t ie = 0; ie < nInteg; ++ie) {
	arrayGodFlux[it][jt][kt][ie][iDim] = FR[ie] + SR*(uStarR[ie]-UR[ie]);
	// variant 1:
	//arrayGodFlux[it][jt][kt][ie][iDim] = (SStar*(SR*UR[ie]-FR[ie]) + SR*(UR[varPrs]+UR[varDen]*(SR-UR[varVel[0]])*(SStar-UR[varVel[0]]))*DStar[ie]) / (SR-SStar);
	// variant 2:
	//arrayGodFlux[it][jt][kt][ie][iDim] = (SStar*(SR*UR[ie]-FR[ie]) + SR*pLR*DStar[ie]) / (SR-SStar);	
      }
    } else if (SR <= 0) {
      for (std::size_t ie = 0; ie < nInteg; ++ie)
	arrayGodFlux[it][jt][kt][ie][iDim] = FR[ie];
    }
  }

  void Fluid::UtoW() { // converting conserved variables into primitive
    negPrsDen = false;
    
    for (std::size_t i = 0; i < gridNx && !negPrsDen; ++i)  // stop if negPrsDen
      for (std::size_t j = 0; j < gridNy && !negPrsDen; ++j)  // stop if negPrsDen
	for (std::size_t k = 0; k < gridNz && !negPrsDen; ++k) {  // stop if negPrsDen

	  for (std::size_t m = 0; m < nDim; ++m)
	    arrayU[i][j][k][varVel[m]] = arrayU[i][j][k][varMom[m]] / arrayU[i][j][k][varDen];

	  arrayU[i][j][k][varPrs] = 0;
	  for (std::size_t m = 0; m < nDim; ++m) {
	    arrayU[i][j][k][varPrs] += pow(arrayU[i][j][k][varVel[m]],2)/2 ;
	  }
	  arrayU[i][j][k][varPrs] = (arrayU[i][j][k][varEng] - arrayU[i][j][k][varDen]*arrayU[i][j][k][varPrs]) * (gama-1);

	  if (arrayU[i][j][k][varPrs] <= 0) {
	    debugInf << std::setw(OWID) << " UtoW:prs<0";
	    negPrsDen = true;
	  }
	  if (arrayU[i][j][k][varDen] <= 0) {
	    debugInf << std::setw(OWID) << " UtoW:den<0";
	    negPrsDen = true;
	  }
	}
  }

  void Fluid::WtoU() { // converting primitive variables into conserved
    for (std::size_t i = 0; i < gridNx; ++i)
      for (std::size_t j = 0; j < gridNy; ++j)
	for (std::size_t k = 0; k < gridNz; ++k) {
	  for (std::size_t m = 0; m < nDim; ++m)
	    arrayU[i][j][k][varMom[m]] = arrayU[i][j][k][varDen] * arrayU[i][j][k][varVel[m]];

	  arrayU[i][j][k][varEng] = 0;
	  for (std::size_t m = 0; m < nDim; ++m)
	    arrayU[i][j][k][varEng] += arrayU[i][j][k][varDen] * pow(arrayU[i][j][k][varVel[m]],2)/2 ;
	  arrayU[i][j][k][varEng] += arrayU[i][j][k][varPrs] / (gama-1);
	}
  }

  void Fluid::getPtclInfo(std::vector<Particle *> &ptcls) {
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

  void Fluid::calcPtclForce(std::vector<Particle *> &ptcls) {
    // must clear forces each loop, otherwise Fluid::plot prints wrong values;
    // but Fluid::penalize works OK since it uses masks.
    for (std::size_t i = 0; i < gridNx ; ++i)
      for (std::size_t j = 0; j < gridNy; ++j)
	for (std::size_t k = 0; k < gridNz; ++k)
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
      REAL avgDen = 0, avgVel = 0, avgPrs = 0, avgVelGap = 0;
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
	avgVelGap += vfabs(Vec(uxFluid-ux, uyFluid-uy, uzFluid-uz));

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
 	if (i > 0 && i < gridNx-1 && j > 0 && j < gridNy-1 && k > 0 && k < gridNz-1 ) { // do not use (i-1) for std::size_t because (i-1) is postive when i=0
	  arrayPressureForce[i][j][k][0] = -(arrayU[i+1][j][k][varPrs] - arrayU[i-1][j][k][varPrs])/(2*gridDx);
	  arrayPressureForce[i][j][k][1] = -(arrayU[i][j+1][k][varPrs] - arrayU[i][j-1][k][varPrs])/(2*gridDy);
	  arrayPressureForce[i][j][k][2] = -(arrayU[i][j][k+1][varPrs] - arrayU[i][j][k-1][varPrs])/(2*gridDz);
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
      avgVelGap /= fluidGrid.size();

      penalForce *= gridDx*gridDy*gridDz;
      presForce  *= gridDx*gridDy*gridDz;
      (*it)->addForce(penalForce*(Cdi/Cd+1));
      (*it)->addForce(presForce);

      penalMoment *= gridDx*gridDy*gridDz;
      presMoment  *= gridDx*gridDy*gridDz;
      (*it)->addMoment(penalMoment*(Cdi/Cd+1));
      (*it)->addMoment(presMoment);

      for (std::size_t iPrn = 0; iPrn < printPtcls.size(); ++iPrn) {
	if ((*it)->getId() == printPtcls[iPrn]) {
	  char cstr[50];
	  std::fstream pfs;
	  pfs.open (dem::combineStr(cstr, "particle_", printPtcls[iPrn], 7), std::fstream::out | std::fstream::app);
	  if(!pfs) { debugInf << "stream error: Fluid::calcPtclForce" << std::endl; exit(-1); }
	  pfs.setf(std::ios::scientific, std::ios::floatfield);
	  if (iteration == 1) {
	    pfs << std::setw(OWID) << "iteration"
		<< std::setw(OWID) << "accruedTime"
		<< std::setw(OWID) << "penalFx"
		<< std::setw(OWID) << "penalFy"
		<< std::setw(OWID) << "penalFz"
		<< std::setw(OWID) << "pressureFx"
		<< std::setw(OWID) << "pressureFy"
		<< std::setw(OWID) << "pressureFz"
		<< std::setw(OWID) << "viscousCd"
		<< std::setw(OWID) << "pressureCd"
		<< std::setw(OWID) << "totalCd"
		<< std::setw(OWID) << "penalMx"
		<< std::setw(OWID) << "penalMy"
		<< std::setw(OWID) << "penalMz"
		<< std::setw(OWID) << "pressureMx"
		<< std::setw(OWID) << "pressureMy"
		<< std::setw(OWID) << "pressureMz"
		<< std::setw(OWID) << "accelX"
		<< std::setw(OWID) << "accelY"
		<< std::setw(OWID) << "accelZ"
		<< std::setw(OWID) << "velocX"
		<< std::setw(OWID) << "velocY"
		<< std::setw(OWID) << "velocZ"
		<< std::setw(OWID) << "avgDen"
		<< std::setw(OWID) << "avgVel"
		<< std::setw(OWID) << "avgPrs"
		<< std::setw(OWID) << "avgVelGap"
		<< std::endl;
	  }

	  // refF is only for the case of Rankine-Hugoniot Condition to test drag coefficients
	  REAL refF = 0.5*rhoL*uL*uL*dem::Pi*(*it)->getA()*(*it)->getB();
	  pfs << std::setw(OWID) << iteration
	      << std::setw(OWID) << timeAccrued

	      << std::setw(OWID) << penalForce.getX()
	      << std::setw(OWID) << penalForce.getY()
	      << std::setw(OWID) << penalForce.getZ()
	      << std::setw(OWID) << presForce.getX()
	      << std::setw(OWID) << presForce.getY()
	      << std::setw(OWID) << presForce.getZ()

	      << std::setw(OWID) << penalForce.getZ()/refF
	      << std::setw(OWID) << presForce.getZ()/refF
	      << std::setw(OWID) << (penalForce.getZ() + presForce.getZ())/refF

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
	      << std::setw(OWID) << avgVelGap

	      << std::endl ;
	  pfs.close();
	}
      }

    } // end of particle loop  
  }

  void Fluid::checkMomentum(std::vector<Particle *> &ptcls) {
    // A: mass and momentum of the flow field outside of the particles plus momentum of the particles as a function of time
    // B: momentum of the entire flow including inside of the particles plus momentum of the particles.
    Vec momAFluid=0, momBFluid=0, momPtcl=0, momA=0, momB=0;
    for (std::size_t i = 0; i < gridNx ; ++i)
      for (std::size_t j = 0; j < gridNy; ++j)
	for (std::size_t k = 0; k < gridNz; ++k) {
	  momAFluid += Vec(arrayU[i][j][k][varMom[0]], arrayU[i][j][k][varMom[1]], arrayU[i][j][k][varMom[2]]) * (1-arrayU[i][j][k][varMsk]); 
	  momBFluid += Vec(arrayU[i][j][k][varMom[0]], arrayU[i][j][k][varMom[1]], arrayU[i][j][k][varMom[2]]);
	}
    momAFluid *= gridDx*gridDy*gridDz;
    momBFluid *= gridDx*gridDy*gridDz;

    for (std::vector<Particle *>::const_iterator it = ptcls.begin(); it != ptcls.end(); ++it)
      momPtcl += (*it)->getCurrVeloc() * (*it)->getMass();
    momA = momAFluid + momPtcl;
    momB = momBFluid + momPtcl;

    char cstr[50];
    std::fstream pfs;
    pfs.open ("momentum_progress", std::fstream::out | std::fstream::app);
    if(!pfs) { debugInf << "stream error: momentum_progress" << std::endl; exit(-1); }
    pfs.setf(std::ios::scientific, std::ios::floatfield);
    if (iteration == 1) {
      pfs << std::setw(OWID) << "iteration"
	  << std::setw(OWID) << "momAFluidZ"
	  << std::setw(OWID) << "momBFluidZ"
	  << std::setw(OWID) << "momPtclZ"
	  << std::setw(OWID) << "momAZ"
	  << std::setw(OWID) << "momBZ"
	  << std::setw(OWID) << "momDiffZ"
	  << std::endl;
    }

    pfs << std::setw(OWID) << iteration
	<< std::setw(OWID) << momAFluid.getZ()
	<< std::setw(OWID) << momBFluid.getZ()
	<< std::setw(OWID) << momPtcl.getZ()
	<< std::setw(OWID) << momA.getZ()
	<< std::setw(OWID) << momB.getZ()
	<< std::setw(OWID) << momB.getZ()-momA.getZ()
	<< std::endl;

    pfs.close();
  }

  void Fluid::plot(const char *str) const {
    std::ofstream ofs(str);
    if(!ofs) { debugInf << "stream error: Fluid::plot" << std::endl; exit(-1); }
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

    ofs << "ZONE I=" << gridNx -2
	<< ", J=" << gridNy -2
	<< ", K=" << gridNz -2
	<< ", DATAPACKING=POINT"
	<< std::endl;

    for (std::size_t k = 1; k < gridNz - 1; ++k)
      for (std::size_t j = 1; j < gridNy - 1; ++j)
	for (std::size_t i = 1; i < gridNx -1; ++i) {
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

  void Fluid::findPrsVel(REAL UL[], REAL UR[], REAL &p, REAL &u) {
    // purpose: to compute the solution for pressure and
    //          velocity in the Star Region

    REAL dl = UL[varDen];
    REAL dr = UR[varDen];
    REAL ul = UL[varVel[0]];
    REAL ur = UR[varVel[0]];
    REAL pl = UL[varPrs];
    REAL pr = UR[varPrs];
    REAL cl = sqrt(gama*pl/dl);
    REAL cr = sqrt(gama*pr/dr);

    const int maxIter = 20;
    REAL change, fl, fld, fr, frd, pPrev, pstart, du;

    guessPressure(UL, UR, pstart);
    pPrev = pstart;
    du = ur - ul;

    int i = 1;
    for ( ; i <= maxIter; ++i) {
      evalF(fl, fld, pPrev, dl, pl, cl);
      evalF(fr, frd, pPrev, dr, pr, cr);
      p = pPrev - (fl + fr + du)/(fld + frd);
      change = 2.0*abs((p - pPrev)/(p + pPrev));
      if (change <= TOL)
	break;
      if (p < 0.0)
	p = TOL;
      pPrev = p;
    }
    if (i > maxIter)
      debugInf << " NR diverge!" << std::endl;

    // compute velocity in star region
    u = 0.5*(ul + ur + fr - fl);
  }

  void Fluid::guessPressure(REAL UL[], REAL UR[], REAL &pInit) {
    // purpose: to provide a guessed value for pressure
    //          pInit in the Star Region. The choice is made
    //          according to adaptive Riemann solver using
    //          the PVRS, TRRS and TSRS approximate
    //          Riemann solvers. See Sect. 9.5 of Chapt. 9 of Ref. 1

    REAL dl = UL[varDen];
    REAL dr = UR[varDen];
    REAL ul = UL[varVel[0]];
    REAL ur = UR[varVel[0]];
    REAL pl = UL[varPrs];
    REAL pr = UR[varPrs];
    REAL cl = sqrt(gama*pl/dl);
    REAL cr = sqrt(gama*pr/dr);  

    const REAL qUser = 2.0;

    // compute guess pressure from PVRS Riemann solver
    REAL pPV = std::max(0.5*(pl + pr) + 0.5*(ul - ur)*0.25*(dl + dr)*(cl + cr), 0.0);
    REAL pMin = std::min(pl, pr);
    REAL pMax = std::max(pl, pr);
    REAL qMax = pMax/pMin;

    if (qMax <= qUser && (pMin <= pPV && pPV <= pMax)) // select PVRS Riemann solver
      pInit = pPV;     
    else {
      if (pPV < pMin) { // select Two-Rarefaction Riemann solver
	REAL pLR = pow(pl/pr, (gama-1.0)/(2.0*gama));
	REAL uStar = (pLR*ul/cl + ur/cr + 2.0/(gama-1.0)*(pLR - 1.0))/(pLR/cl + 1.0/cr);
	REAL coefL = 1.0 + (gama-1.0)/2.0*(ul - uStar)/cl;
	REAL coefR = 1.0 + (gama-1.0)/2.0*(uStar - ur)/cr;
	pInit = 0.5*(pl*pow(coefL, static_cast<int> (2.0*gama/(gama-1))) + pr*pow(coefR, static_cast<int> (2.0*gama/(gama-1)))); // Equ.(9.36) of Toro.
	//pInit = pow((cl+cr-(gama-1)/2*(ul-ur))/(cl/pow(pl,(gama-1)/2/gama)+cr/pow(pr,(gama-1)/2/gama)), 2.0*gama/(gama-1)); // Equ.(9.32) of Toro.
	pInit = std::max(pInit, TOL);
	//debugInf << " !!!pInit=" << pInit << " pLR= " << pLR << " uStar=" << uStar << " coefL=" << coefL << " coefR=" << coefR << " cl=" << cl << " cr=" << cr << " ul=" << ul << " ur=" << ur;
      } 
      else { // select Two-Shock Riemann solver with PVRS as estimate
	REAL gL = sqrt((2.0/(gama+1.0)/dl)/((gama-1.0)/(gama+1.0)*pl + pPV));
	REAL gR = sqrt((2.0/(gama+1.0)/dr)/((gama-1.0)/(gama+1.0)*pr + pPV));
	pInit = (gL*pl + gR*pr - (ur - ul))/(gL + gR);   
      }
    }
  }

  void Fluid::evalF(REAL &f,
		    REAL &fd,
		    REAL &p,
		    REAL &dk,
		    REAL &pk,
		    REAL &ck) {
    // purpose: to evaluate the pressure functions
    //          fl and fr in exact Riemann solver
    //          and their first derivatives

    if (p <= pk) {
      // rarefaction wave
      f = 2.0/(gama-1.0)*ck*(pow(p/pk, (gama-1.0)/(2.0*gama)) - 1.0);
      fd = (1.0/(dk*ck))*pow(p/pk, -(gama+1.0)/(2.0*gama));
    } else {
      // shock wave
      REAL ak = 2.0/(gama+1.0)/dk;
      REAL bk = (gama-1.0)/(gama+1.0)*pk;
      f = (p - pk)*sqrt(ak/(bk + p));
      fd = (1.0 - 0.5*(p - pk)/(bk + p))*sqrt(ak/(bk + p));
    }
  }

  void Fluid::sampling(REAL UL[], REAL UR[],
		       const REAL pStar,
		       const REAL uStar,
		       const REAL s,
		       REAL &d,
		       REAL &u,
		       REAL &p) {
    // purpose: to sample the solution throughout the wave
    //          pattern. Pressure pStar and velocity uStar in the
    //          star region are known. Sampling is performed
    //          in terms of the 'speed' s = x/t. Sampled
    //          values are d, u, p

    REAL dl = UL[varDen];
    REAL dr = UR[varDen];
    REAL ul = UL[varVel[0]];
    REAL ur = UR[varVel[0]];
    REAL pl = UL[varPrs];
    REAL pr = UR[varPrs];
    REAL cl = sqrt(gama*pl/dl);
    REAL cr = sqrt(gama*pr/dr);  

    REAL c, pStarL, pStarR, sHL, sHR, sL, sR, sTL, sTR;
    //debugInf << " pStar=" << pStar << " uStar=" << uStar << " s=" << s;
    if (s <= uStar) {
      // sampling point lies to the left of the contact discontinuity
      if (pStar <= pl) {
	// left rarefaction
	sHL = ul - cl;
	//debugInf << "  leftrare  sHL=" << sHL;
	if (s <= sHL) {
	  // sampled point is left data state
	  d = dl;
	  u = ul;
	  p = pl;
	} else {
	  sTL = uStar - cl*pow(pStar/pl, (gama-1.0)/(2.0*gama));
	  //debugInf << " sTL=" << sTL;
	  if (s > sTL) {
	    // sampled point is star left state
	    d = dl*pow(pStar/pl, 1.0/gama);
	    u = uStar;
	    p = pStar;
	    //debugInf << " ul=" << ul << " ur=" << ur << " u=" << u << " d=" << d;
	  } else {
	    // sampled point is inside left fan
	    u = 2.0/(gama+1.0)*(cl + (gama-1.0)/2.0*ul + s);
	    c = 2.0/(gama+1.0)*(cl + (gama-1.0)/2.0*(ul - s));
	    d = dl*pow(c/cl, 2.0/(gama-1.0));
	    p = pl*pow(c/cl, 2.0*gama/(gama-1.0));
	    //debugInf << " ul=" << ul << " ur=" << ur << " u=" << u << " d=" << d;
	  }
	}
      } else {
	// left shock
	pStarL = pStar/pl;
	sL = ul - cl*sqrt((gama+1.0)/(2.0*gama)*pStarL + (gama-1.0)/(2.0*gama));
	//debugInf << "  leftshock sL=" << sL;
	if (s <= sL) {
	  // sampled point is left data state
	  d = dl;
	  u = ul;
	  p = pl;
	} else {
	  // sampled point is star left state
	  d = dl*(pStarL + (gama-1.0)/(gama+1.0))/(pStarL*(gama-1.0)/(gama+1.0) + 1.0);
	  u = uStar;
	  p = pStar;
	}
      }
    } else {
      // sampling point lies to the right of the contact discontinuity
      if (pStar > pr) {
	// right shock
	pStarR = pStar/pr;
	sR  = ur + cr*sqrt((gama+1.0)/(2.0*gama)*pStarR + (gama-1.0)/(2.0*gama));
	//debugInf << " rightshock sR=" << sR;
	if (s >= sR) {
	  // sampled point is right data state
	  d = dr;
	  u = ur;
	  p = pr;
	} else {
	  // sampled point is star right state
	  d = dr*(pStarR + (gama-1.0)/(gama+1.0))/(pStarR*(gama-1.0)/(gama+1.0) + 1.0);
	  u = uStar;
	  p = pStar;
	}
      } else {
	// right rarefaction
	sHR = ur + cr;
	//debugInf << " rightrare  sHR=" << sHR;
	if (s >= sHR) {
	  // sampled point is right data state
	  d = dr;
	  u = ur;
	  p = pr;
	} else {
	  sTR = uStar + cr*pow(pStar/pr, (gama-1.0)/(2.0*gama));
	  if (s <= sTR) {
	    // sampled point is star right state
	    d = dr*pow(pStar/pr, 1.0/gama);
	    u = uStar;
	    p = pStar;
	  } else {
	    // sampled point is inside left fan
	    u = 2.0/(gama+1.0)*(-cr + (gama-1.0)/2.0*ur + s);
	    c = 2.0/(gama+1.0)*(cr - (gama-1.0)/2.0*(ur - s));
	    d = dr*pow(c/cr, 2.0/(gama-1.0));
	    p = pr*pow(c/cr, 2.0*gama/(gama-1.0));
	  }
	}
      }
    }
  }

} // name space dem
