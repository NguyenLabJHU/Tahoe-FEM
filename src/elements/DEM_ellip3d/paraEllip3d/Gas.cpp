#include "Gas.h"
#include "const.h"
#include "block.h"
#include <cmath>
#include <string>
#include <algorithm>

namespace dem {

  std::string combineStr(const char *str, std::size_t num, std::size_t width) {
    std::stringstream ss;
    ss << std::setw(width) << std::setfill('0') << std::right << num;
    std::string obj(str);
    obj += ss.str();
    return obj;
  }

  const REAL Gas::Rs = 287.06; // specific gas constant

  void Gas::initPureGasParameter() {
    initSharedParameter();
    gridDz = dem::Parameter::getSingleton().parameter["gridSize"];;
    calcGrid();
    haloGrid = 1;
    printSharedParameter();
  }

  void Gas::initParameter(Gradation &gradation) {
    initSharedParameter();
    ptclGrid = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["ptclGrid"]);
    porosity = dem::Parameter::getSingleton().parameter["porosity"];
    Cdi = dem::Parameter::getSingleton().parameter["Cdi"];
    printPtcls = dem::Parameter::getSingleton().cfdPrintPtcls;

    gridDz = (gradation.getPtclMinRadius(1) * 2) / ptclGrid; // estimate
    calcGrid(); // re-calculate gridDz
    haloGrid = static_cast<std::size_t> (round(gradation.getPtclMaxRadius() / gridDz) ) + 1; // round or ceil? round is more accurate; necessary + 1 for pressure gradient calculation.
  
    debugInf << std::setw(OWID) << "ptclGrid" << std::setw(OWID) << ptclGrid << std::endl;
    debugInf << std::setw(OWID) << "porosity" << std::setw(OWID) << porosity << std::endl;
    debugInf << std::setw(OWID) << "Cdi" << std::setw(OWID) << Cdi << std::endl;
    printSharedParameter();
  }

  void Gas::initSharedParameter() {
    Cd   = dem::Parameter::getSingleton().parameter["Cd"];
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
  }

  void Gas::calcGrid() {
    allGridNz = static_cast<std::size_t> (ceil((z2F - z1F) / gridDz));
    gridDz = (z2F - z1F) / allGridNz; // re-calculate
    gridDx = gridDz;
    gridDy = gridDz;
    allGridNx = static_cast<std::size_t> (ceil((x2F - x1F) / gridDx)); // ceil rather than round to cover particle domain.
    allGridNy = static_cast<std::size_t> (ceil((y2F - y1F) / gridDy));

    allGridNx += 2; // add two boundary cells
    allGridNy += 2;
    allGridNz += 2;
  }

  void Gas::printSharedParameter() {
    if (mpi.mpiRank == 0) {
      debugInf << std::setw(OWID) << "nVar" << std::setw(OWID) << nVar << std::endl;
      debugInf << std::setw(OWID) << "nInteg" << std::setw(OWID) << nInteg << std::endl;

      debugInf << std::setw(OWID) << "Cd" << std::setw(OWID) << Cd << std::endl;
      debugInf << std::setw(OWID) << "Runge-Kutta" << std::setw(OWID) << (int) RK << std::endl;
      debugInf << std::setw(OWID) << "CFL" << std::setw(OWID) << CFL << std::endl;
      debugInf << std::setw(OWID) << "gama" << std::setw(OWID) << gama << std::endl;
      debugInf << std::setw(OWID) << "rhoR" << std::setw(OWID) << rhoR << std::endl;
      debugInf << std::setw(OWID) << "pR" << std::setw(OWID) << pR << std::endl;
      debugInf << std::setw(OWID) << "uR" << std::setw(OWID) << uR << std::endl;

      debugInf << std::setw(OWID) << "gridDx" << std::setw(OWID) << gridDx << std::endl;
      debugInf << std::setw(OWID) << "gridDy" << std::setw(OWID) << gridDy << std::endl;
      debugInf << std::setw(OWID) << "gridDz" << std::setw(OWID) << gridDz << std::endl;
      debugInf << std::setw(OWID) << "allGridNx" << std::setw(OWID) << allGridNx << std::endl;
      debugInf << std::setw(OWID) << "allGridNy" << std::setw(OWID) << allGridNy << std::endl;
      debugInf << std::setw(OWID) << "allGridNz" << std::setw(OWID) << allGridNz << std::endl;

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
      if (aNum != ceil(aNum) && aNum != floor(aNum)) // not a perfect method
	debugInf << "fluid domain in z direction is NOT precisely grided!" << std::endl;
      aNum = (x2F - x1F) / gridDx;
      if (aNum != ceil(aNum) && aNum != floor(aNum))
	debugInf << "fluid domain in x direction is NOT precisely grided!" << std::endl;
      aNum = (y2F - y1F) / gridDy;
      if (aNum != ceil(aNum) && aNum != floor(aNum))
	debugInf << "fluid domain in y direction is NOT precisely grided!" << std::endl;

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
    }
  }

  void Gas::allocArray() {   
    // step 1: determine haloGrid
    haloGridX = haloGrid;
    haloGridY = haloGrid;
    haloGridZ = haloGrid;
    if (mpi.mpiProcX == 1) haloGridX = 0;
    if (mpi.mpiProcY == 1) haloGridY = 0;
    if (mpi.mpiProcZ == 1) haloGridZ = 0;

    // step 2: determine gridN
    gridNx = BLOCK_SIZE(mpi.mpiCoords[0], mpi.mpiProcX, allGridNx) + 2*haloGridX;
    gridNy = BLOCK_SIZE(mpi.mpiCoords[1], mpi.mpiProcY, allGridNy) + 2*haloGridY;
    gridNz = BLOCK_SIZE(mpi.mpiCoords[2], mpi.mpiProcZ, allGridNz) + 2*haloGridZ;
 
    if (mpi.isBdryProcessXMin() || mpi.isBdryProcessXMax())
      gridNx = BLOCK_SIZE(mpi.mpiCoords[0], mpi.mpiProcX, allGridNx) + haloGridX;
    if (mpi.isBdryProcessYMin() || mpi.isBdryProcessYMax())
      gridNy = BLOCK_SIZE(mpi.mpiCoords[1], mpi.mpiProcY, allGridNy) + haloGridY;
    if (mpi.isBdryProcessZMin() || mpi.isBdryProcessZMax())
      gridNz = BLOCK_SIZE(mpi.mpiCoords[2], mpi.mpiProcZ, allGridNz) + haloGridZ;

    /*
    // for ceil/floor method, no longer needed, but keep here for record.
    gridNx = (int) ceil((double) allGridNx / mpi.mpiProcX) + 2*haloGridX;
    gridNy = (int) ceil((double) allGridNy / mpi.mpiProcY) + 2*haloGridY;
    gridNz = (int) ceil((double) allGridNz / mpi.mpiProcZ) + 2*haloGridZ;
 
    if (mpi.isBdryProcessXMin())
      gridNx = (int) ceil((double) allGridNx / mpi.mpiProcX) + haloGridX;
    if (mpi.isBdryProcessYMin())
      gridNy = (int) ceil((double) allGridNy / mpi.mpiProcY) + haloGridY;
    if (mpi.isBdryProcessZMin())
      gridNz = (int) ceil((double) allGridNz / mpi.mpiProcZ) + haloGridZ;

    if (mpi.isBdryProcessXMax()) {
      gridNx = allGridNx - (int) ceil((double) allGridNx / mpi.mpiProcX) * (mpi.mpiProcX-1);
      if (gridNx < haloGridX) {
	std::cout << "Error! mpiRank=" << mpi.mpiRank << ", gridNx=" << gridNx << " smaller than haloGridX=" << haloGridX << ", use smaller mpiProcX!" << std::endl;
	abort();
      }
      gridNx += haloGridX;
    }
    if (mpi.isBdryProcessYMax()) {
      gridNy = allGridNy - (int) ceil((double) allGridNy / mpi.mpiProcY) * (mpi.mpiProcY-1);
      if (gridNy < haloGridY) {
	std::cout << "Error! mpiRank=" << mpi.mpiRank << ", gridNy=" << gridNy << " smaller than haloGridY=" << haloGridY << ", use smaller mpiProcY!" << std::endl;
	abort();
      }
      gridNy += haloGridY;
    }
    if (mpi.isBdryProcessZMax()) {
      gridNz = allGridNz - (int) ceil((double) allGridNz / mpi.mpiProcZ) * (mpi.mpiProcZ-1);
      if (gridNz < haloGridZ) {
	std::cout << "Error! mpiRank=" << mpi.mpiRank << ", gridNz=" << gridNz << " smaller than haloGridZ=" << haloGridZ << ", use smaller mpiProcZ!" << std::endl;
	abort();
      }
      gridNz += haloGridZ;
    }
    */

    if (mpi.mpiRank == 0) {
      debugInf << std::setw(OWID) << "haloGridX=" << std::setw(OWID) << haloGridX << std::endl;
      debugInf << std::setw(OWID) << "haloGridY=" << std::setw(OWID) << haloGridY << std::endl;
      debugInf << std::setw(OWID) << "haloGridZ=" << std::setw(OWID) << haloGridZ << std::endl;
    }

    // step 3: determine bound
    // boundCup: size = gridN, for both serial and parallel.
    boundCup = Bound(0, 0, 0, gridNx - 1, gridNy - 1, gridNz - 1);

    // boundPrn: size = gridN - haloGrid*2, for parallel.
    boundPrn = Bound(haloGridX, haloGridY, haloGridZ, gridNx - haloGridX - 1, gridNy - haloGridY - 1, gridNz - haloGridZ - 1);
    if (mpi.isBdryProcessXMin()) {boundPrn.lowX = 1; boundPrn.sizeX = gridNx - haloGridX - 1;}
    if (mpi.isBdryProcessYMin()) {boundPrn.lowY = 1; boundPrn.sizeY = gridNy - haloGridY - 1;}
    if (mpi.isBdryProcessZMin()) {boundPrn.lowZ = 1; boundPrn.sizeZ = gridNz - haloGridZ - 1;}
    if (mpi.isBdryProcessXMax()) {boundPrn.uppX = gridNx - 2; boundPrn.sizeX = gridNx - haloGridX - 1;}
    if (mpi.isBdryProcessYMax()) {boundPrn.uppY = gridNy - 2; boundPrn.sizeY = gridNy - haloGridY - 1;}
    if (mpi.isBdryProcessZMax()) {boundPrn.uppZ = gridNz - 2; boundPrn.sizeZ = gridNz - haloGridZ - 1;}

    // boundPrn: for serial in a specific direction.
    if (mpi.mpiProcX == 1) {
      boundPrn.lowX = 1; boundPrn.uppX = gridNx - 2; boundPrn.sizeX = gridNx - 2;
    }
    if (mpi.mpiProcY == 1) {
      boundPrn.lowY = 1; boundPrn.uppY = gridNy - 2; boundPrn.sizeY = gridNy - 2;
    }
    if (mpi.mpiProcZ == 1) {
      boundPrn.lowZ = 1; boundPrn.uppZ = gridNz - 2; boundPrn.sizeZ = gridNz - 2;
    }

    // boundGod: size = gridN -1, for both serial and parallel.
    boundGod = Bound(0, 0, 0, gridNx - 2, gridNy - 2, gridNz - 2);

    /*
    std::cout << std::endl << "mpiRank = " << mpi.mpiRank << std::endl;
    IJK globalLow, globalUpp;
    IJK local;
    local = IJK(boundCup.lowX, boundCup.lowY, boundCup.lowZ);
    localIndexToGlobal(local, globalLow);
    local = IJK(boundCup.uppX, boundCup.uppY, boundCup.uppZ);
    localIndexToGlobal(local, globalUpp);
    std::cout << std::setw(OWID) << "boundCupX: size, local[ ], global[ ]" << std::setw(OWID) << boundCup.sizeX << std::setw(OWID) << boundCup.lowX << std::setw(OWID) << boundCup.uppX << std::setw(OWID) << globalLow.i << std::setw(OWID) << globalUpp.i << std::endl;
    std::cout << std::setw(OWID) << "boundCupY: size, local[ ], global[ ]" << std::setw(OWID) << boundCup.sizeY << std::setw(OWID) << boundCup.lowY << std::setw(OWID) << boundCup.uppY << std::setw(OWID) << globalLow.j << std::setw(OWID) << globalUpp.j << std::endl;
    std::cout << std::setw(OWID) << "boundCupZ: size, local[ ], global[ ]" << std::setw(OWID) << boundCup.sizeZ << std::setw(OWID) << boundCup.lowZ << std::setw(OWID) << boundCup.uppZ << std::setw(OWID) << globalLow.k << std::setw(OWID) << globalUpp.k << std::endl;

    local = IJK(boundPrn.lowX, boundPrn.lowY, boundPrn.lowZ);
    localIndexToGlobal(local, globalLow);
    local = IJK(boundPrn.uppX, boundPrn.uppY, boundPrn.uppZ);
    localIndexToGlobal(local, globalUpp);
    std::cout << std::setw(OWID) << "boundPrnX: size, local[ ], global[ ]" << std::setw(OWID) << boundPrn.sizeX << std::setw(OWID) << boundPrn.lowX << std::setw(OWID) << boundPrn.uppX << std::setw(OWID) << globalLow.i << std::setw(OWID) << globalUpp.i << std::endl;
    std::cout << std::setw(OWID) << "boundPrnY: size, local[ ], global[ ]" << std::setw(OWID) << boundPrn.sizeY << std::setw(OWID) << boundPrn.lowY << std::setw(OWID) << boundPrn.uppY << std::setw(OWID) << globalLow.j << std::setw(OWID) << globalUpp.j << std::endl;
    std::cout << std::setw(OWID) << "boundPrnZ: size, local[ ], global[ ]" << std::setw(OWID) << boundPrn.sizeZ << std::setw(OWID) << boundPrn.lowZ << std::setw(OWID) << boundPrn.uppZ << std::setw(OWID) << globalLow.k << std::setw(OWID) << globalUpp.k << std::endl;

    std::cout << std::endl;
    */

    // boundCup.sizeX, boundCup.sizeY, boundCup.sizeZ, nDim
    // resize does zero-initialization
    arrayGridCoord.resize(boundCup.sizeX);
    for (std::size_t i = 0; i < arrayGridCoord.size(); ++i) {
      arrayGridCoord[i].resize(boundCup.sizeY);
      for (std::size_t j = 0; j < arrayGridCoord[i].size(); ++j) {
	arrayGridCoord[i][j].resize(boundCup.sizeZ);
	for (std::size_t k = 0; k < arrayGridCoord[i][j].size(); ++k) 
	  arrayGridCoord[i][j][k].resize(nDim);
      }
    }
    // coordinates
    for (std::size_t i = 0; i < arrayGridCoord.size(); ++i)
      for (std::size_t j = 0; j < arrayGridCoord[i].size(); ++j)
	for (std::size_t k = 0; k < arrayGridCoord[i][j].size(); ++k) {
	  // i, j, k are local, must be mapped to global
	  IJK local(i,j,k);
	  IJK global;
	  localIndexToGlobal(local, global);
	  int iG = global.i;
	  int jG = global.j;
	  int kG = global.k;

	  arrayGridCoord[i][j][k][0] = (x1F - gridDx/2) + iG * gridDx;
	  arrayGridCoord[i][j][k][1] = (y1F - gridDy/2) + jG * gridDy;
	  arrayGridCoord[i][j][k][2] = (z1F - gridDz/2) + kG * gridDz;
	}

    // boundCup.sizeX, boundCup.sizeY, boundCup.sizeZ, nDim
    arrayPenalForce.resize(boundCup.sizeX);
    for (std::size_t i = 0; i < arrayPenalForce.size(); ++i) {
      arrayPenalForce[i].resize(boundCup.sizeY);
      for (std::size_t j = 0; j < arrayPenalForce[i].size(); ++j) {
	arrayPenalForce[i][j].resize(boundCup.sizeZ);
	for (std::size_t k = 0; k < arrayPenalForce[i][j].size(); ++k) 
	  arrayPenalForce[i][j][k].resize(nDim);
      }
    }

    // boundCup.sizeX, boundCup.sizeY, boundCup.sizeZ, nDim
    arrayPressureForce.resize(boundCup.sizeX);
    for (std::size_t i = 0; i < arrayPressureForce.size(); ++i) {
      arrayPressureForce[i].resize(boundCup.sizeY);
      for (std::size_t j = 0; j < arrayPressureForce[i].size(); ++j) {
	arrayPressureForce[i][j].resize(boundCup.sizeZ);
	for (std::size_t k = 0; k < arrayPressureForce[i][j].size(); ++k) 
	  arrayPressureForce[i][j][k].resize(nDim);
      }
    }

    // boundCup.sizeX, boundCup.sizeY, boundCup.sizeZ, nVar
    arrayU.resize(boundCup.sizeX);
    for (std::size_t i = 0; i < arrayU.size(); ++i) {
      arrayU[i].resize(boundCup.sizeY);
      for (std::size_t j = 0; j < arrayU[i].size(); ++j) {
	arrayU[i][j].resize(boundCup.sizeZ);
	for (std::size_t k = 0; k < arrayU[i][j].size(); ++k) 
	  arrayU[i][j][k].resize(nVar);
      }
    }

    // boundCup.sizeX, boundCup.sizeY, boundCup.sizeZ, nVar
    arrayURota.resize(boundCup.sizeX);
    for (std::size_t i = 0; i < arrayURota.size(); ++i) {
      arrayURota[i].resize(boundCup.sizeY);
      for (std::size_t j = 0; j < arrayURota[i].size(); ++j) {
	arrayURota[i][j].resize(boundCup.sizeZ);
	for (std::size_t k = 0; k < arrayURota[i][j].size(); ++k) 
	  arrayURota[i][j][k].resize(nVar);
      }
    }

    // boundCup.sizeX, boundCup.sizeY, boundCup.sizeZ, nInteg
    arrayFlux.resize(boundCup.sizeX);
    for (std::size_t i = 0; i < arrayFlux.size(); ++i) {
      arrayFlux[i].resize(boundCup.sizeY);
      for (std::size_t j = 0; j < arrayFlux[i].size(); ++j) {
	arrayFlux[i][j].resize(boundCup.sizeZ);
	for (std::size_t k = 0; k < arrayFlux[i][j].size(); ++k) 
	  arrayFlux[i][j][k].resize(nInteg);
      }
    }

    // boundGod.sizeX, boundGod.sizeY, boundGod.sizeZ, nInteg, nDim
    arrayGodFlux.resize(boundGod.sizeX);
    for (std::size_t i = 0; i < arrayGodFlux.size(); ++i) {
      arrayGodFlux[i].resize(boundGod.sizeY);
      for (std::size_t j = 0; j < arrayGodFlux[i].size(); ++j) {
	arrayGodFlux[i][j].resize(boundGod.sizeZ);
	for (std::size_t k = 0; k < arrayGodFlux[i][j].size(); ++k) {
	  arrayGodFlux[i][j][k].resize(nInteg);
	  for (std::size_t m = 0; m < arrayGodFlux[i][j][k].size(); ++m)
	    arrayGodFlux[i][j][k][m].resize(nDim);
	}
      }
    }

    if (RK >= 2) {
      // boundCup.sizeX, boundCup.sizeY, boundCup.sizeZ, nVar
      arrayUN.resize(boundCup.sizeX);
      for (std::size_t i = 0; i < arrayUN.size(); ++i) {
	arrayUN[i].resize(boundCup.sizeY);
	for (std::size_t j = 0; j < arrayUN[i].size(); ++j) {
	  arrayUN[i][j].resize(boundCup.sizeZ);
	  for (std::size_t k = 0; k < arrayUN[i][j].size(); ++k) 
	    arrayUN[i][j][k].resize(nVar);
	}
      }
    }

    if (RK == 3) {
      // boundCup.sizeX, boundCup.sizeY, boundCup.sizeZ, nVar
      arrayUStar.resize(boundCup.sizeX);
      for (std::size_t i = 0; i < arrayUStar.size(); ++i) {
	arrayUStar[i].resize(boundCup.sizeY);
	for (std::size_t j = 0; j < arrayUStar[i].size(); ++j) {
	  arrayUStar[i][j].resize(boundCup.sizeZ);
	  for (std::size_t k = 0; k < arrayUStar[i][j].size(); ++k) 
	    arrayUStar[i][j][k].resize(nVar);
	}
      }

      // boundCup.sizeX, boundCup.sizeY, boundCup.sizeZ, nVar
      arrayUStar2.resize(boundCup.sizeX);
      for (std::size_t i = 0; i < arrayUStar2.size(); ++i) {
	arrayUStar2[i].resize(boundCup.sizeY);
	for (std::size_t j = 0; j < arrayUStar2[i].size(); ++j) {
	  arrayUStar2[i][j].resize(boundCup.sizeZ);
	  for (std::size_t k = 0; k < arrayUStar2[i][j].size(); ++k) 
	    arrayUStar2[i][j][k].resize(nVar);
	}
      }
    }

    // boundGod.sizeX, boundGod.sizeY, boundGod.sizeZ, nInteg
    arrayGodFluxTmp.resize(boundGod.sizeX);
    for (std::size_t i = 0; i < arrayGodFluxTmp.size(); ++i) {
      arrayGodFluxTmp[i].resize(boundGod.sizeY);
      for (std::size_t j = 0; j < arrayGodFluxTmp[i].size(); ++j) {
	arrayGodFluxTmp[i][j].resize(boundGod.sizeZ);
	for (std::size_t k = 0; k < arrayGodFluxTmp[i][j].size(); ++k) {
	  arrayGodFluxTmp[i][j][k].resize(nInteg);
	}
      }
    }
    
    // boundCup.sizeX, boundCup.sizeY, boundCup.sizeZ
    arrayH.resize(boundCup.sizeX);
    for (std::size_t i = 0; i < arrayH.size(); ++i) {
      arrayH[i].resize(boundCup.sizeY);
      for (std::size_t j = 0; j < arrayH[i].size(); ++j)
	arrayH[i][j].resize(boundCup.sizeZ);
    }

    // boundCup.sizeX, boundCup.sizeY, boundCup.sizeZ
    arraySoundSpeed.resize(boundCup.sizeX);
    for (std::size_t i = 0; i < arraySoundSpeed.size(); ++i) {
      arraySoundSpeed[i].resize(boundCup.sizeY);
      for (std::size_t j = 0; j < arraySoundSpeed[i].size(); ++j)
	arraySoundSpeed[i][j].resize(boundCup.sizeZ);
    }
  }

  void Gas::coordToGlobalIndex(REAL x, REAL y, REAL z, std::size_t &i, std::size_t &j, std::size_t &k) {
    // important to use round!
    i = static_cast<std::size_t> ( round((x - (x1F - gridDx/2)) / gridDx ));
    j = static_cast<std::size_t> ( round((y - (y1F - gridDy/2)) / gridDy ));
    k = static_cast<std::size_t> ( round((z - (z1F - gridDz/2)) / gridDz ));
  }

  void Gas::coordToGlobalIndex(Vec v, std::size_t &i, std::size_t &j, std::size_t &k) {
    // important to use round!
    i = static_cast<std::size_t> ( round((v.getX() - (x1F - gridDx/2)) / gridDx ));
    j = static_cast<std::size_t> ( round((v.getY() - (y1F - gridDy/2)) / gridDy ));
    k = static_cast<std::size_t> ( round((v.getZ() - (z1F - gridDz/2)) / gridDz ));
  }

  void Gas::coordToGlobalIndex(Vec v, IJK &t) {
    // important to use round!
    t.i = static_cast<std::size_t> ( round((v.getX() - (x1F - gridDx/2)) / gridDx ));
    t.j = static_cast<std::size_t> ( round((v.getY() - (y1F - gridDy/2)) / gridDy ));
    t.k = static_cast<std::size_t> ( round((v.getZ() - (z1F - gridDz/2)) / gridDz ));
  }

  void Gas::localIndexToGlobal(IJK &local, IJK &global) {
    // do no use std::size_t for computing
    int i = local.i;
    int j = local.j;
    int k = local.k;

    int I = BLOCK_LOW(mpi.mpiCoords[0], mpi.mpiProcX, allGridNx) - (int) haloGridX + i;
    int J = BLOCK_LOW(mpi.mpiCoords[1], mpi.mpiProcY, allGridNy) - (int) haloGridY + j;
    int K = BLOCK_LOW(mpi.mpiCoords[2], mpi.mpiProcZ, allGridNz) - (int) haloGridZ + k;

    /*
    // for ceil/floor method, no longer needed, but keep here for record.
    int I = (int) ceil((double) allGridNx / mpi.mpiProcX) * mpi.mpiCoords[0] - (int) haloGridX + i;
    int J = (int) ceil((double) allGridNy / mpi.mpiProcY) * mpi.mpiCoords[1] - (int) haloGridY + j;
    int K = (int) ceil((double) allGridNz / mpi.mpiProcZ) * mpi.mpiCoords[2] - (int) haloGridZ + k;
    */

    // ensure starting from 0 for min boundary process
    if (mpi.isBdryProcessXMin()) I = i;
    if (mpi.isBdryProcessYMin()) J = j;
    if (mpi.isBdryProcessZMin()) K = k;
    // automatically end at (allGrid -1) for max boundary process

    /*
    if (I < 0 || I > allGridNx -1 || J < 0 || J > allGridNy -1 || K < 0 || K > allGridNz -1) 
      std::cout << "Gas::localIndexToGlobal: global IJK out of range, IJK= " << I << " " << J << " " << K << std::endl;
    */

    global.i = I;
    global.j = J;
    global.k = K;
  }

  bool Gas::globalIndexToLocal(IJK &global, IJK &local) {
    // do no use std::size_t for computing
    int i = (int) global.i - BLOCK_LOW(mpi.mpiCoords[0], mpi.mpiProcX, allGridNx) + (int) haloGridX;
    int j = (int) global.j - BLOCK_LOW(mpi.mpiCoords[1], mpi.mpiProcY, allGridNy) + (int) haloGridY;
    int k = (int) global.k - BLOCK_LOW(mpi.mpiCoords[2], mpi.mpiProcZ, allGridNz) + (int) haloGridZ;

    /*
    // for ceil/floor method, no longer needed, but keep here for record. The BLOCK method uses BLOCK_LOW! 
    int segX = (int) ceil((double) allGridNx / mpi.mpiProcX);
    int segY = (int) ceil((double) allGridNy / mpi.mpiProcY);
    int segZ = (int) ceil((double) allGridNz / mpi.mpiProcZ);
    int i = (int) global.i % segX + (int) haloGridX;
    int j = (int) global.j % segY + (int) haloGridY;
    int k = (int) global.k % segZ + (int) haloGridZ;

    // when it nears process boundary, a particle is duplicated to multiple processes by MPI communication,
    // the local indexes of the particle center in different MPI processes are different and must be handled.
    IJK localLowBound(boundPrn.lowX, boundPrn.lowY, boundPrn.lowZ);
    IJK localUppBound(boundPrn.uppX, boundPrn.uppY, boundPrn.uppZ);
    IJK globalLowBound;
    IJK globalUppBound;
    localIndexToGlobal(localLowBound, globalLowBound);
    localIndexToGlobal(localUppBound, globalUppBound);
    //std::cout << "globalIndexToLocal: mpiRank=" << mpi.mpiRank << " boundPrn.lowX=" << boundPrn.lowX << " globalLowBound.i= " << globalLowBound.i << " i=" << i; 

    if (global.i > globalUppBound.i) // particle center is located in the "right" or "upper"process.
      i += segX;
    if (global.i < globalLowBound.i) // particle center is located in the "left" or "lower" process.
      i -= segX;

    if (global.j > globalUppBound.j)
      j += segY;
    if (global.j < globalLowBound.j)
      j -= segY;

    if (global.k > globalUppBound.k)
      k += segZ;
    if (global.k < globalLowBound.k)
      k -= segZ;
    // end of dealing with a particle nearing process boundary.
    */

    // ensure starting from 0 for min boundary process
    if (mpi.isBdryProcessXMin()) i = (int) global.i;
    if (mpi.isBdryProcessYMin()) j = (int) global.j;
    if (mpi.isBdryProcessZMin()) k = (int) global.k;
    // automatically end at (allGrid -1) for max boundary process

    local.i = i;
    local.j = j;
    local.k = k;
    
    if (i < boundCup.lowX || i > boundCup.uppX ||
	j < boundCup.lowY || j > boundCup.uppY ||
	k < boundCup.lowZ || k > boundCup.uppZ) // outside of halo
      return false;

    return true; // inside of halo
  }

  void Gas::initializePureGas() {
    negPrsDen = false;
    RankineHugoniot();
    initialCondition(); 
    soundSpeed(); // for printing Mach number

    if (mpi.mpiRank == 0) {
      debugInf << std::setw(OWID) << "iteration" 
	       << std::setw(OWID) << "demTimeStep"
	       << std::setw(OWID) << "cfdTimeStep" 
	       << std::setw(OWID) << "timeStep" 
	       << std::setw(OWID) << "timeAccrued"

	       << std::setw(OWID) << "cfdCommuTime"
	       << std::setw(OWID) << "cfdRunTime"
	       << std::setw(OWID) << "cfdCommu%";
    }
  }

  void Gas::initialize() {
    negPrsDen = false;
    RankineHugoniot();
    initialCondition(); 
    soundSpeed(); // for printing Mach number

    if (mpi.mpiRank == 0) {
      debugInf << std::setw(OWID) << "iteration" 
	       << std::setw(OWID) << "demTimeStep"
	       << std::setw(OWID) << "cfdTimeStep" 
	       << std::setw(OWID) << "timeStep" 
	       << std::setw(OWID) << "timeAccrued"

	       << std::setw(OWID) << "cfdCommuT"
	       << std::setw(OWID) << "getPtclInfoT"
	       << std::setw(OWID) << "runOneStepT"
	       << std::setw(OWID) << "calcPtclForceT"
	       << std::setw(OWID) << "penalizeT"
	       << std::setw(OWID) << "cfdTotalT"
	       << std::setw(OWID) << "demCommuT"
	       << std::setw(OWID) << "demCompuT"
	       << std::setw(OWID) << "demMigraT"
	       << std::setw(OWID) << "demTotalT"
	       << std::setw(OWID) << "totalT";
    }
  }

  void Gas::runOneStep(std::vector<Particle *> &ptcls) {
    if (RK == 1) {
      inteStep1(ptcls);
    }
    else if (RK == 2) {
      arrayUN = arrayU;
      inteStep1(ptcls);     // arrayU->arrayU* & arrayGodFlux->arrayGodFlux*
      RK2InteStep2(ptcls); 
    }
    else if (RK == 3) {
      arrayUN = arrayU;
      inteStep1(ptcls);     // arrayU->arrayU* & arrayGodFlux->arrayGodFlux*
      arrayUStar = arrayU;  // arrayU*
      RK3InteStep2(ptcls);
      arrayUStar2 = arrayU; // arrayU**
      RK3InteStep3(ptcls);
    }
  }

  void Gas::inteStep1(std::vector<Particle *> &ptcls) { 
    initGhostPoints();
    soundSpeed();
    calcTimeStep();
    enthalpy();
    rotateIJK(ptcls);

    // update conserved variables at the next time step
    std::size_t i, j, k, m;
    int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
#pragma omp parallel for num_threads(ompThreads) private(i, j, k, m) schedule(dynamic)
    for (i = boundCup.lowX + 1; i <= boundCup.uppX - 1; ++i)
      for (j = boundCup.lowY + 1; j <= boundCup.uppY - 1; ++j)
	for (k = boundCup.lowZ + 1; k <= boundCup.uppZ - 1; ++k)
	  for (m = 0; m < nInteg; ++m) {
	    arrayU[i][j][k][m] -= (   timeStep / gridDx * (arrayGodFlux[i][j][k][m][0] - arrayGodFlux[i-1][j][k][m][0])
				    + timeStep / gridDy * (arrayGodFlux[i][j][k][m][1] - arrayGodFlux[i][j-1][k][m][1])
				    + timeStep / gridDz * (arrayGodFlux[i][j][k][m][2] - arrayGodFlux[i][j][k-1][m][2]) );
	  }

    UtoW();
  }
  
  void Gas::RK2InteStep2(std::vector<Particle *> &ptcls) { 
    initGhostPoints();
    soundSpeed();
    enthalpy();
    rotateIJK(ptcls);

    // update conserved variables at the next time step
    std::size_t i, j, k, m;
    int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
#pragma omp parallel for num_threads(ompThreads) private(i, j, k, m) schedule(dynamic)
    for (i = boundCup.lowX + 1; i <= boundCup.uppX - 1; ++i)
      for (j = boundCup.lowY + 1; j <= boundCup.uppY - 1; ++j)
	for (k = boundCup.lowZ + 1; k <= boundCup.uppZ -1; ++k)
	  for (m = 0; m < nInteg; ++m) {
	    arrayU[i][j][k][m] = 0.5*( arrayUN[i][j][k][m] + arrayU[i][j][k][m] -
				       (  timeStep / gridDx * (arrayGodFlux[i][j][k][m][0] - arrayGodFlux[i-1][j][k][m][0]) 
				        + timeStep / gridDy * (arrayGodFlux[i][j][k][m][1] - arrayGodFlux[i][j-1][k][m][1]) 
					+ timeStep / gridDz * (arrayGodFlux[i][j][k][m][2] - arrayGodFlux[i][j][k-1][m][2]) ) );
	  } 

    UtoW();
  }

  void Gas::RK3InteStep2(std::vector<Particle *> &ptcls) { 
    initGhostPoints();
    soundSpeed();
    enthalpy();
    rotateIJK(ptcls);

    // update conserved variables at the next time step
    std::size_t i, j, k, m;
    int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
#pragma omp parallel for num_threads(ompThreads) private(i, j, k, m) schedule(dynamic)
    for (i = boundCup.lowX + 1; i <= boundCup.uppX - 1; ++i)
      for (j = boundCup.lowY + 1; j <= boundCup.uppY - 1; ++j)
	for (k = boundCup.lowZ + 1; k <= boundCup.uppZ - 1; ++k)
	  for (m = 0; m < nInteg; ++m) {
	    arrayU[i][j][k][m] = 0.75*arrayUN[i][j][k][m] + 0.25*arrayUStar[i][j][k][m] - 
	      0.25* (  timeStep / gridDx * (arrayGodFlux[i][j][k][m][0] - arrayGodFlux[i-1][j][k][m][0])
		     + timeStep / gridDy * (arrayGodFlux[i][j][k][m][1] - arrayGodFlux[i][j-1][k][m][1])
		     + timeStep / gridDz * (arrayGodFlux[i][j][k][m][2] - arrayGodFlux[i][j][k-1][m][2]) );
	  }

    UtoW();
  }

  void Gas::RK3InteStep3(std::vector<Particle *> &ptcls) { 
    initGhostPoints();
    soundSpeed();
    enthalpy();
    rotateIJK(ptcls);

    // update conserved variables at the next time step
    std::size_t i, j, k, m;
    int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
#pragma omp parallel for num_threads(ompThreads) private(i, j, k, m) schedule(dynamic)
    for (i = boundCup.lowX + 1; i <= boundCup.uppX - 1; ++i)
      for (j = boundCup.lowY + 1; j <= boundCup.uppY - 1; ++j)
	for (k = boundCup.lowZ + 1; k <= boundCup.uppZ - 1; ++k)
	  for (m = 0; m < nInteg; ++m) {
	    arrayU[i][j][k][m] = ( arrayUN[i][j][k][m] + 2*arrayUStar2[i][j][k][m] - 
                                   2* (  timeStep / gridDx * (arrayGodFlux[i][j][k][m][0] - arrayGodFlux[i-1][j][k][m][0])
				       + timeStep / gridDy * (arrayGodFlux[i][j][k][m][1] - arrayGodFlux[i][j-1][k][m][1])
				       + timeStep / gridDz * (arrayGodFlux[i][j][k][m][2] - arrayGodFlux[i][j][k-1][m][2]) ) )/3.0;
	  }

    UtoW();
  }

  void Gas::rotateIJK(std::vector<Particle *> &ptcls) {
    std::size_t id[3][3] = {{0,1,2},{1,0,2},{2,1,0}};

    // for x, y, z sweeps: 
    // iDim=0 for x sweep
    // iDim=1 for y sweep
    // iDim=2 for z sweep
    for (std::size_t iDim = 0; iDim < nDim; ++iDim) {

      arrayURota = arrayU; // works well for std::vector; for std::valarray they must have same rank and extent;
      
      // switch components
      for (std::size_t i = 0; i < arrayURota.size(); ++i)
	for (std::size_t j = 0; j < arrayURota[i].size(); ++j)
	  for (std::size_t k = 0; k < arrayURota[i][j].size(); ++k) 
	    for (std::size_t jdim = 0; jdim < nDim; ++jdim) {
	      arrayURota[i][j][k][  varMom[jdim]  ] = arrayU[i][j][k][  varMom[id[iDim][jdim]]  ];
	      arrayURota[i][j][k][  varVel[jdim]  ] = arrayU[i][j][k][  varVel[id[iDim][jdim]]  ];
	    }

      flux(iDim, ptcls); // variables defined at cell centers

      int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
      // the value of a private variable is undefined when parallel construct is entered, so do not assign value here.
      std::size_t iGod, jGod, kGod; // defined for arrayGodFlux
      std::size_t IL[3], IR[3];
      REAL UL[9], UR[9], FL[5], FR[5], HL, HR, avgH, avgU, avgV, avgW, avgh, pStar;
      int solver;

      // for local Riemann problem
#pragma omp parallel for num_threads(ompThreads) private(iGod,jGod,kGod,IL,IR,UL,UR,FL,FR,HL,HR,avgH,avgU,avgV,avgW,avgh,pStar,solver) schedule(dynamic)
      for (iGod = boundGod.lowX; iGod <= boundGod.uppX; ++iGod) { // variables defined at cell faces
	for (jGod = boundGod.lowY; jGod <= boundGod.uppY; ++jGod) {
	  for (kGod = boundGod.lowZ; kGod <= boundGod.uppZ; ++kGod) {

	    IL[0]=iGod; IL[1]=jGod; IL[2]=kGod; 
	    IR[0]=iGod; IR[1]=jGod; IR[2]=kGod; 
	    IR[iDim] += 1;    
	    HL = arrayH[IL[0]] [IL[1]] [IL[2]];
	    HR = arrayH[IR[0]] [IR[1]] [IR[2]];

	    for (std::size_t m = 0; m < nVar - 1; ++m) { // NOT nVar, because varMsk is not needed; using nVar(=10) goes out of bound.
	      UL[m] = arrayURota[IL[0]] [IL[1]] [IL[2]] [m];
	      UR[m] = arrayURota[IR[0]] [IR[1]] [IR[2]] [m];
	    }	
	    for (std::size_t m = 0; m < nInteg; ++m) {
	      FL[m] = arrayFlux[IL[0]] [IL[1]] [IL[2]] [m];
	      FR[m] = arrayFlux[IR[0]] [IR[1]] [IR[2]] [m];
	    }    

	    avgH = (sqrt(UL[varDen])*HL + sqrt(UR[varDen])*HR)/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
	    avgU = (sqrt(UL[varDen])*UL[varVel[0]] + sqrt(UR[varDen])*UR[varVel[0]])/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
	    avgV = (sqrt(UL[varDen])*UL[varVel[1]] + sqrt(UR[varDen])*UR[varVel[1]])/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
	    avgW = (sqrt(UL[varDen])*UL[varVel[2]] + sqrt(UR[varDen])*UR[varVel[2]])/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
	    avgh = avgH - 0.5*(avgU*avgU + avgV*avgV + avgW*avgW); // static specific enthalpy
	    guessPressure(UL, UR, pStar);

	    solver = static_cast<int> (dem::Parameter::getSingleton().parameter["solver"]);
	    switch (solver) {
	    case 0: 
	      RoeSolver(UL, UR, FL, FR, HL, HR, iDim, iGod, jGod, kGod);
	      break;
	    case 1:
	      HllcSolver(UL, UR, FL, FR, HL, HR, iDim, iGod, jGod, kGod);
	      break;
	    case 2:
	      HlleSolver(UL, UR, FL, FR, HL, HR, iDim, iGod, jGod, kGod);
	      break;
	    case 3:
	      if (pStar <= UL[varPrs] || pStar <= UR[varPrs] || avgh <= 0 || negPrsDen) // first 2 expressions cover rarefaction and contact waves
		HllcSolver(UL, UR, FL, FR, HL, HR, iDim, iGod, jGod, kGod);
	      else // use Roe solver for shock waves
		RoeSolver(UL, UR, FL, FR, HL, HR, iDim, iGod, jGod, kGod); 
	      break;
	    case 4:
	      if (pStar <= UL[varPrs] || pStar <= UR[varPrs] || avgh <= 0 || negPrsDen) // first 2 expressions cover rarefaction and contact waves
		HlleSolver(UL, UR, FL, FR, HL, HR, iDim, iGod, jGod, kGod);
	      else // use Roe solver for shock waves
		RoeSolver(UL, UR, FL, FR, HL, HR, iDim, iGod, jGod, kGod);
	      break;
	    case 5:
	      if (pStar <= UL[varPrs] || pStar <= UR[varPrs] || avgh <= 0 || negPrsDen) // first 2 expressions cover rarefaction and contact waves
		exactSolver(UL, UR, 0, iDim, iGod, jGod, kGod);
	      else // use Roe solver for shock waves
		RoeSolver(UL, UR, FL, FR, HL, HR, iDim, iGod, jGod, kGod);
	      break;	
	    case -1:
	      // exactSolver itself can provide exact solution to the full domain, if run separately.
	      // herein exactSolver is used at each discritized face to give exact solution to local Riemann problem, for the purpose 
	      // of comparing with other solvers.
	      exactSolver(UL, UR, 0, iDim, iGod, jGod, kGod); // 0 = x/t = s
	      break;
	    case -2:
	      LaxFrieScheme(UL, UR, FL, FR, iDim, iGod, jGod, kGod); // Lax-Friedrichs scheme
	      break;
	    case -3:
	      LaxWendScheme(UL, UR, FL, FR, iDim, iGod, jGod, kGod); // Lax-Wendroff scheme (two-step Richtmyer version)
	      break;
	    }

	  }
	}
      } // end of local Riemann problem

      for (std::size_t i = 0; i < arrayGodFluxTmp.size(); ++i)
	for (std::size_t j = 0; j < arrayGodFluxTmp[i].size(); ++j)
	  for (std::size_t k = 0; k < arrayGodFluxTmp[i][j].size(); ++k)
	    for (std::size_t m = 0; m < nInteg; ++m)
	      arrayGodFluxTmp[i][j][k][m] = arrayGodFlux[i][j][k][m][iDim];

      // switch components back for consistency with u
      for (std::size_t i = 0; i < arrayGodFlux.size(); ++i)
	for (std::size_t j = 0; j < arrayGodFlux[i].size(); ++j)
	  for (std::size_t k = 0; k < arrayGodFlux[i][j].size(); ++k)
	    for (std::size_t m = 0; m < nDim; ++m)
	      arrayGodFlux[i][j][k][varMom[m]][iDim] = arrayGodFluxTmp[i][j][k][ varMom[id[iDim][m]] ];
    } // end of for x, y, z directions

  }

  void Gas::penalize(std::vector<Particle *> &ptcls) {
    // this implementation has higher efficiency than that of checking every fluid cell.
    // for cells that are occupied by particle volumes
    for (std::vector<Particle *>::const_iterator it = ptcls.begin(); it != ptcls.end(); ++it) {
      std::vector< std::vector<REAL> > fluidGrid = (*it)->getFluidGrid();
      for (std::size_t iter = 0; iter < fluidGrid.size(); ++iter) {

	std::size_t i = static_cast<std::size_t> (fluidGrid[iter][0]);
	std::size_t j = static_cast<std::size_t> (fluidGrid[iter][1]);
	std::size_t k = static_cast<std::size_t> (fluidGrid[iter][2]);

	REAL denQuot[3], u0[3];
 	if (i >= boundCup.lowX+1 && i+1 <= boundCup.uppX) // necessary if getPtclInfo does not use +/-1; unnecessary if getPtclInfo uses +/-1, albeit kept here.
	  denQuot[0] = (arrayU[i+1][j][k][varDen] - arrayU[i-1][j][k][varDen]) / (2*gridDx);
	if (j >= boundCup.lowY+1 && j+1 <= boundCup.uppY)
	  denQuot[1] = (arrayU[i][j+1][k][varDen] - arrayU[i][j-1][k][varDen]) / (2*gridDy);
	if (k >= boundCup.lowZ+1 && k+1 <= boundCup.uppZ)
	  denQuot[2] = (arrayU[i][j][k+1][varDen] - arrayU[i][j][k-1][varDen]) / (2*gridDz);
	REAL coordX = arrayGridCoord[i][j][k][0];
	REAL coordY = arrayGridCoord[i][j][k][1];
	REAL coordZ = arrayGridCoord[i][j][k][2];
	Vec dist = Vec(coordX, coordY, coordZ) - (*it)->getCurrPos();
	Vec omgar = (*it)->getCurrOmga() % dist; // w X r = omga % dist, where % is overloaded as cross product
	u0[0] = (*it)->getCurrVeloc().getX() + omgar.getX(); 
	u0[1] = (*it)->getCurrVeloc().getY() + omgar.getY(); 
	u0[2] = (*it)->getCurrVeloc().getZ() + omgar.getZ();

	// 1. momentum penalization
	for (std::size_t m = 0; m < nDim; ++m) {
	  // a. momentum penalization, note Cdi is controlled by porosity while Cd is not (Cd/Cd=1)
	  arrayU[i][j][k][varMom[m]] -= arrayU[i][j][k][varMsk]*arrayPenalForce[i][j][k][m]*timeStep * (Cdi/Cd/porosity + 1);
	  // b. influence of momentum penalization on energy, note Cdi is controlled by porosity while Cd is not (Cd/Cd=1)
	  arrayU[i][j][k][varEng]    -= arrayU[i][j][k][varMsk]*arrayPenalForce[i][j][k][m]*arrayU[i][j][k][varVel[m]]*timeStep * (Cdi/Cd/porosity + 1);
	}

	// 2. mass penalization
	//   a. mass source term is incorporated by changes in flux(), so it is not computed here.

	//   b. influence of mass source term u0[m]*denQuot[m] on energy, it is approximate to d(rho*u0)/dx.
	for (std::size_t m = 0; m < nDim; ++m)
	  arrayU[i][j][k][varEng]  += arrayU[i][j][k][varMsk] * (-0.5*pow(arrayU[i][j][k][varVel[m]],2)) * (u0[m]*denQuot[m]) * timeStep / porosity;

      }
    }
  }
  
  void Gas::initGhostPoints() {
    // non-reflecting BCs
    // -x
    if (mpi.isBdryProcessXMin()) {
      for (std::size_t j = boundCup.lowY; j <= boundCup.uppY; ++j)
	for (std::size_t k = boundCup.lowZ; k <= boundCup.uppZ; ++k)
	  for (std::size_t m = 0; m < nVar; ++m) {
	    arrayU[boundCup.lowX][j][k][m] = arrayU[boundCup.lowX+1][j][k][m]; 
	  }
    }

    // +x
    if (mpi.isBdryProcessXMax()) {
      for (std::size_t j = boundCup.lowY; j <= boundCup.uppY; ++j)
	for (std::size_t k = boundCup.lowZ; k <= boundCup.uppZ; ++k)
	  for (std::size_t m = 0; m < nVar; ++m) {
	    arrayU[boundCup.uppX][j][k][m] = arrayU[boundCup.uppX-1][j][k][m]; 
	  }
    }
    // -y
    if (mpi.isBdryProcessYMin()) {
      for (std::size_t i = boundCup.lowX; i <= boundCup.uppX; ++i)
	for (std::size_t k = boundCup.lowZ; k <= boundCup.uppZ; ++k)
	  for (std::size_t m = 0; m < nVar; ++m) {
	    arrayU[i][boundCup.lowY][k][m] = arrayU[i][boundCup.lowY+1][k][m]; 
	  }
    }
    // +y
    if (mpi.isBdryProcessYMax()) {
      for (std::size_t i = boundCup.lowX; i <= boundCup.uppX; ++i)
	for (std::size_t k = boundCup.lowZ; k <= boundCup.uppZ; ++k)
	  for (std::size_t m = 0; m < nVar; ++m) {
	    arrayU[i][boundCup.uppY][k][m] = arrayU[i][boundCup.uppY-1][k][m]; 
	  }
    }
    // -z
    if (mpi.isBdryProcessZMin()) {
      for (std::size_t i = boundCup.lowX; i <= boundCup.uppX; ++i)
	for (std::size_t j = boundCup.lowY; j <= boundCup.uppY; ++j)
	  for (std::size_t m = 0; m < nVar; ++m) {
	    arrayU[i][j][boundCup.lowZ][m] = arrayU[i][j][boundCup.lowZ+1][m]; 
	  }
    }
    // +z
    if (mpi.isBdryProcessZMax()) {
      for (std::size_t i = boundCup.lowX; i <= boundCup.uppX; ++i)
	for (std::size_t j = boundCup.lowY; j <= boundCup.uppY; ++j)
	  for (std::size_t m = 0; m < nVar; ++m) {
	    arrayU[i][j][boundCup.uppZ][m] = arrayU[i][j][boundCup.uppZ-1][m]; 
	  }
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
      // -x
      if (mpi.isBdryProcessXMin()) {
	for (std::size_t j = boundCup.lowY; j <= boundCup.uppY; ++j)
	  for (std::size_t k = boundCup.lowZ; k <= boundCup.uppZ; ++k)
	    for (std::size_t m = 0; m < 1; ++m) {
	      arrayU[boundCup.lowX][j][k][varMom[m]]    *= (1-2*arrayBC[0]); 
	      arrayU[boundCup.lowX][j][k][varVel[m]]    *= (1-2*arrayBC[0]); 
	    }
      }
      // +x
      if (mpi.isBdryProcessXMax()) {
	for (std::size_t j = boundCup.lowY; j <= boundCup.uppY; ++j)
	  for (std::size_t k = boundCup.lowZ; k <= boundCup.uppZ; ++k)
	    for (std::size_t m = 0; m < 1; ++m) {
	      arrayU[boundCup.uppX][j][k][varMom[m]] *= (1-2*arrayBC[1]); 
	      arrayU[boundCup.uppX][j][k][varVel[m]] *= (1-2*arrayBC[1]); 
	    }
      }
      // -y
      if (mpi.isBdryProcessYMin()) {
	for (std::size_t i = boundCup.lowX; i <= boundCup.uppX; ++i)
	  for (std::size_t k = boundCup.lowZ; k <= boundCup.uppZ; ++k)
	    for (std::size_t m = 1; m < 2; ++m) {
	      arrayU[i][boundCup.lowY][k][varMom[m]]    *= (1-2*arrayBC[2]); 
	      arrayU[i][boundCup.lowY][k][varVel[m]]    *= (1-2*arrayBC[2]); 
	    }
      }
      // +y
      if (mpi.isBdryProcessYMax()) {
	for (std::size_t i = boundCup.lowX; i <= boundCup.uppX; ++i)
	  for (std::size_t k = boundCup.lowZ; k <= boundCup.uppZ; ++k)
	    for (std::size_t m = 1; m < 2; ++m) {
	      arrayU[i][boundCup.uppY][k][varMom[m]] *= (1-2*arrayBC[3]);
	      arrayU[i][boundCup.uppY][k][varVel[m]] *= (1-2*arrayBC[3]);  
	    }
      }
      // -z
      if (mpi.isBdryProcessZMin()) {
	for (std::size_t i = boundCup.lowX; i <= boundCup.uppX; ++i)
	  for (std::size_t j = boundCup.lowY; j <= boundCup.uppY; ++j)
	    for (std::size_t m = 2; m < 3; ++m) {
	      arrayU[i][j][boundCup.lowZ][varMom[m]]    *= (1-2*arrayBC[4]); 
	      arrayU[i][j][boundCup.lowZ][varVel[m]]    *= (1-2*arrayBC[4]); 
	    }
      }
      // +z
      if (mpi.isBdryProcessZMax()) {
	for (std::size_t i = boundCup.lowX; i <= boundCup.uppX; ++i)
	  for (std::size_t j = boundCup.lowY; j <= boundCup.uppY; ++j)
	    for (std::size_t m = 2; m < 3; ++m) {
	      arrayU[i][j][boundCup.uppZ][varMom[m]] *= (1-2*arrayBC[5]); 
	      arrayU[i][j][boundCup.uppZ][varVel[m]] *= (1-2*arrayBC[5]); 
	    }
      }

    }  

  }

  void Gas::calcTimeStep() {
    std::vector<REAL> valX(boundPrn.sizeX * boundPrn.sizeY * boundPrn.sizeZ);
    std::vector<REAL> valY(boundPrn.sizeX * boundPrn.sizeY * boundPrn.sizeZ);
    std::vector<REAL> valZ(boundPrn.sizeX * boundPrn.sizeY * boundPrn.sizeZ);

    std::size_t i, j, k;
    int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
#pragma omp parallel for num_threads(ompThreads) private(i, j, k) schedule(dynamic)
    for (i = boundPrn.lowX; i <= boundPrn.uppX; ++i)
      for (j = boundPrn.lowY; j <= boundPrn.uppY; ++j)
	for (k = boundPrn.lowZ; k <= boundPrn.uppZ; ++k) {
	  std::size_t i1 = i - boundPrn.lowX;
	  std::size_t j1 = j - boundPrn.lowY;
	  std::size_t k1 = k - boundPrn.lowZ;

	  valX[i1 + j1 * boundPrn.sizeX + k1 * boundPrn.sizeX * boundPrn.sizeY] = fabs(arrayU[i][j][k][varVel[0]]) + arraySoundSpeed[i][j][k];
	  valY[i1 + j1 * boundPrn.sizeX + k1 * boundPrn.sizeX * boundPrn.sizeY] = fabs(arrayU[i][j][k][varVel[1]]) + arraySoundSpeed[i][j][k];
	  valZ[i1 + j1 * boundPrn.sizeX + k1 * boundPrn.sizeX * boundPrn.sizeY] = fabs(arrayU[i][j][k][varVel[2]]) + arraySoundSpeed[i][j][k];
	}

    std::vector<REAL> dtMin(3);
    dtMin[0] = gridDx / ( *std::max_element(valX.begin(), valX.end()) );
    dtMin[1] = gridDy / ( *std::max_element(valY.begin(), valY.end()) ); 
    dtMin[2] = gridDz / ( *std::max_element(valZ.begin(), valZ.end()) ); 
    
    REAL demTimeStep = timeStep; // minimum across all processes
    REAL pCFDTimeStep = CFL * (*std::min_element(dtMin.begin(), dtMin.end()));
    REAL cfdTimeStep; // mininum across all processes

    // this guarantees that all processes use the same timeStep
    MPI_Allreduce(&pCFDTimeStep, &cfdTimeStep, 1, MPI_DOUBLE, MPI_MIN, mpi.mpiWorld);
    timeStep = std::min(demTimeStep, cfdTimeStep);
    timeAccrued += timeStep;

    if (mpi.mpiRank == 0) {
      debugInf << std::endl
	       << std::setw(OWID) << iteration 
	       << std::setw(OWID) << demTimeStep 
	       << std::setw(OWID) << cfdTimeStep 
	       << std::setw(OWID) << timeStep
	       << std::setw(OWID) << timeAccrued;
    }
  }

  /*
  // std::valarray also works well.  
  void Gas::calcTimeStep() {
    std::valarray<REAL> valX(boundPrn.sizeX * boundPrn.sizeY * boundPrn.sizeZ);
    std::valarray<REAL> valY(boundPrn.sizeX * boundPrn.sizeY * boundPrn.sizeZ);
    std::valarray<REAL> valZ(boundPrn.sizeX * boundPrn.sizeY * boundPrn.sizeZ);

    std::size_t i, j, k;
    int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
#pragma omp parallel for num_threads(ompThreads) private(i, j, k) schedule(dynamic)
    for (i = boundPrn.lowX; i <= boundPrn.uppX; ++i)
      for (j = boundPrn.lowY; j <= boundPrn.uppY; ++j)
	for (k = boundPrn.lowZ; k <= boundPrn.uppZ; ++k) {
	  std::size_t i1 = i - boundPrn.lowX;
	  std::size_t j1 = j - boundPrn.lowY;
	  std::size_t k1 = k - boundPrn.lowZ;

	  valX[i1 + j1 * boundPrn.sizeX + k1 * boundPrn.sizeX * boundPrn.sizeY] = fabs(arrayU[i][j][k][varVel[0]]) + arraySoundSpeed[i][j][k];
	  valY[i1 + j1 * boundPrn.sizeX + k1 * boundPrn.sizeX * boundPrn.sizeY] = fabs(arrayU[i][j][k][varVel[1]]) + arraySoundSpeed[i][j][k];
	  valZ[i1 + j1 * boundPrn.sizeX + k1 * boundPrn.sizeX * boundPrn.sizeY] = fabs(arrayU[i][j][k][varVel[2]]) + arraySoundSpeed[i][j][k];
	}

    std::valarray<REAL> dtMin(3);
    dtMin[0] = gridDx / valX.max();
    dtMin[1] = gridDy / valY.max();
    dtMin[2] = gridDz / valZ.max();

    REAL ptclTimeStep = timeStep;
    REAL fluidTimeStep = CFL * dtMin.min();
    timeStep = std::min(ptclTimeStep, fluidTimeStep);
    timeAccrued += timeStep;
    debugInf << std::endl
	     << std::setw(OWID) << iteration 
	     << std::setw(OWID) << ptclTimeStep 
	     << std::setw(OWID) << fluidTimeStep 
	     << std::setw(OWID) << timeStep 
	     << std::setw(OWID) << timeAccrued;
  }
  */

  void Gas::soundSpeed() {
    std::size_t i, j, k;
    int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
#pragma omp parallel for num_threads(ompThreads) private(i, j, k) schedule(dynamic)
    for (i = 0; i < arraySoundSpeed.size(); ++i)
      for (j = 0; j < arraySoundSpeed[i].size(); ++j)
	for (k = 0; k < arraySoundSpeed[i][j].size(); ++k)
	  arraySoundSpeed[i][j][k] = sqrt(gama * arrayU[i][j][k][varPrs] / arrayU[i][j][k][varDen]);
  }

  // total specific enthalphy, NOT static specific enthalpy
  void Gas::enthalpy() {
    std::size_t i, j, k;
    int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
#pragma omp parallel for num_threads(ompThreads) private(i, j, k) schedule(dynamic)
    for (i = 0; i < arrayH.size(); ++i)
      for (j = 0; j < arrayH[i].size(); ++j)
	for (k = 0; k < arrayH[i][j].size(); ++k)
	  arrayH[i][j][k] = (arrayU[i][j][k][varEng] + arrayU[i][j][k][varPrs]) / arrayU[i][j][k][varDen];
  }

  void Gas::initialCondition() {
    if (leftType == 1 || leftType == 2) { // normal shock w/ and w/o Rankine-Hugoniot conditions
      for (std::size_t i = 0; i < arrayU.size(); ++i)
	for (std::size_t j = 0; j < arrayU[i].size(); ++j)
	  for (std::size_t k = 0; k < arrayU[i][j].size(); ++k) {
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
      for (std::size_t i = 0; i < arrayU.size(); ++i)
	for (std::size_t j = 0; j < arrayU[i].size(); ++j)
	  for (std::size_t k = 0; k < arrayU[i][j].size(); ++k) {
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
      for (std::size_t i = 0; i < arrayU.size(); ++i)
	for (std::size_t j = 0; j < arrayU[i].size(); ++j)
	  for (std::size_t k = 0; k < arrayU[i][j].size(); ++k) {
	    REAL radius = sqrt(pow(arrayGridCoord[i][j][k][0]-x0L,2) + pow(arrayGridCoord[i][j][k][1]-y0L,2) + pow(arrayGridCoord[i][j][k][2]-z0L,2));
	    if (radius <= r0L) {
	      arrayU[i][j][k][varDen]    = rhoL;
	      arrayU[i][j][k][varPrs]    = pL;
	      arrayU[i][j][k][varVel[0]] = uL*(arrayGridCoord[i][j][k][0]-x0L)/r0L;
	      arrayU[i][j][k][varVel[1]] = uL*(arrayGridCoord[i][j][k][1]-y0L)/r0L;
	      arrayU[i][j][k][varVel[2]] = uL*(arrayGridCoord[i][j][k][2]-z0L)/r0L;

	    } else {
	      arrayU[i][j][k][varDen]    = rhoR;
	      arrayU[i][j][k][varPrs]    = pR;
	      arrayU[i][j][k][varVel[0]] = 0;
	      arrayU[i][j][k][varVel[1]] = 0;
	      arrayU[i][j][k][varVel[2]] = 0;
	    }
	  }
    } else if (leftType == 5) { // normal shock z directions, three initial zones
      for (std::size_t i = 0; i < arrayU.size(); ++i)
	for (std::size_t j = 0; j < arrayU[i].size(); ++j)
	  for (std::size_t k = 0; k < arrayU[i][j].size(); ++k) {
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

  void Gas::RankineHugoniot() { // Rankine-Hugoniot conditions
    if (leftType == 1) {
      shockSpeed = MachShock*sqrt(gama*pR/rhoR);
      rhoL = ( pow(rhoR*(shockSpeed-uR),2)*(1+gama) ) / ( rhoR*pow(shockSpeed-uR,2)*(gama-1) + 2*pR*gama);
      pL   = (pR*(1-gama)+2*rhoR*pow(shockSpeed-uR,2)) / (1+gama);
      uL   = ( rhoR*(shockSpeed-uR)*(2*shockSpeed + uR*(gama-1)) - 2*pR*gama ) / (rhoR * (shockSpeed-uR) * (1+gama));
      if (mpi.mpiRank == 0) {
	debugInf << std::setw(OWID) << "shockSpeed" << std::setw(OWID) << shockSpeed << std::endl;
	debugInf << std::setw(OWID) << "rhoL" << std::setw(OWID) << rhoL << std::endl;
	debugInf << std::setw(OWID) << "pL" << std::setw(OWID) << pL << std::endl;
	debugInf << std::setw(OWID) << "uL" << std::setw(OWID) << uL << std::endl;
      }
    } 

    MachL = uL / sqrt(gama*pL/rhoL);
    if (mpi.mpiRank == 0) {
      debugInf << std::setw(OWID) << "MachL" << std::setw(OWID) << MachL << std::endl << std::endl;
    }
  }

  void Gas::flux(std::size_t iDim, std::vector<Particle *> &ptcls) {
    for (std::size_t i = 0; i < arrayFlux.size(); ++i)
      for (std::size_t j = 0; j < arrayFlux[i].size(); ++j)
	for (std::size_t k = 0; k < arrayFlux[i][j].size(); ++k) {
	  arrayFlux[i][j][k][varDen]    = arrayURota[i][j][k][varDen] * arrayURota[i][j][k][varVel[0]]; // rho*u
	  arrayFlux[i][j][k][varMom[0]] = arrayURota[i][j][k][varDen] * pow(arrayURota[i][j][k][varVel[0]],2) + arrayURota[i][j][k][varPrs]; // rho*u^2 + p
	  arrayFlux[i][j][k][varMom[1]] = arrayURota[i][j][k][varDen] * arrayURota[i][j][k][varVel[0]] * arrayURota[i][j][k][varVel[1]]; // rho*u*v
	  arrayFlux[i][j][k][varMom[2]] = arrayURota[i][j][k][varDen] * arrayURota[i][j][k][varVel[0]] * arrayURota[i][j][k][varVel[2]]; // rho*u*w
	  arrayFlux[i][j][k][varEng]    = arrayURota[i][j][k][varVel[0]] * (arrayURota[i][j][k][varEng] + arrayURota[i][j][k][varPrs]);  // u*(E + p)
	}  

    // for cells that are occupied by particle volumes
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

	// equations are modified in terms of porosity and superficial velocity (Darcy's velocity) for porous material
	// continuity equation modified on porosity
	arrayFlux[i][j][k][varDen]    = 1.0/porosity * arrayURota[i][j][k][varDen] * (arrayURota[i][j][k][varVel[0]] - u0[iDim]) ; // rho*(u-u0)/porosity

	// momentum equations modified on porosity
	arrayFlux[i][j][k][varMom[0]] = arrayURota[i][j][k][varDen] * pow(arrayURota[i][j][k][varVel[0]],2) + arrayURota[i][j][k][varPrs]/porosity; // rho*u^2 + p/porosity
	arrayFlux[i][j][k][varMom[1]] = arrayURota[i][j][k][varDen] * arrayURota[i][j][k][varVel[0]] * arrayURota[i][j][k][varVel[1]]; // (rho*u)*v
	arrayFlux[i][j][k][varMom[2]] = arrayURota[i][j][k][varDen] * arrayURota[i][j][k][varVel[0]] * arrayURota[i][j][k][varVel[2]]; // (rho*u)*w

	// energy equation modified on porosity
	arrayFlux[i][j][k][varEng]    = arrayURota[i][j][k][varVel[0]] * (arrayURota[i][j][k][varEng] + arrayURota[i][j][k][varPrs]/porosity);  // u*(E + p/porosity)
      }
    }
  }

  void Gas::LaxFrieScheme(REAL UL[], REAL UR[], REAL FL[], REAL FR[], std::size_t iDim, std::size_t it, std::size_t jt, std::size_t kt) {
    for (std::size_t ie = 0; ie < nInteg; ++ie)
      arrayGodFlux[it][jt][kt][ie][iDim] = (FL[ie]+FR[ie])/2 + (UL[ie]-UR[ie])/2*gridDz/timeStep; // Equ.(5.77) of Toro
  }

  void Gas::LaxWendScheme(REAL UL[], REAL UR[], REAL FL[], REAL FR[], std::size_t iDim, std::size_t it, std::size_t jt, std::size_t kt) { 
    // Equ.(5.79) of Toro
    REAL r  = (UL[0]+UR[0])/2 + (FL[0]-FR[0])/2*timeStep/gridDz;
    REAL ru = (UL[1]+UR[1])/2 + (FL[1]-FR[1])/2*timeStep/gridDz;
    REAL rv = (UL[2]+UR[2])/2 + (FL[2]-FR[2])/2*timeStep/gridDz;
    REAL rw = (UL[3]+UR[3])/2 + (FL[3]-FR[3])/2*timeStep/gridDz;
    REAL E  = (UL[4]+UR[4])/2 + (FL[4]-FR[4])/2*timeStep/gridDz;
    REAL p  = (E-0.5*r*(ru*ru+rv*rv+rw*rw)/(r*r))*(gama-1);
    arrayGodFlux[it][jt][kt][varDen][iDim]    = ru;
    arrayGodFlux[it][jt][kt][varMom[0]][iDim] = ru*ru/r + p;
    arrayGodFlux[it][jt][kt][varMom[1]][iDim] = ru*rv/r;
    arrayGodFlux[it][jt][kt][varMom[2]][iDim] = ru*rw/r;
    arrayGodFlux[it][jt][kt][varEng][iDim]    = ru/r*(E + p);
  }

  void Gas::exactSolver(REAL UL[], REAL UR[], REAL relaCoord, std::size_t iDim, std::size_t it, std::size_t jt, std::size_t kt) {
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

  void Gas::RoeSolver(REAL UL[], REAL UR[], REAL FL[], REAL FR[], REAL HL, REAL HR, std::size_t iDim, std::size_t it, std::size_t jt, std::size_t kt) {
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

  void Gas::HlleSolver(REAL UL[], REAL UR[], REAL FL[], REAL FR[], REAL HL, REAL HR, std::size_t iDim, std::size_t it, std::size_t jt, std::size_t kt) {
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
    } else if (SR <= 0) {
      for (std::size_t ie = 0; ie < nInteg; ++ie)
	arrayGodFlux[it][jt][kt][ie][iDim] = FR[ie];
    }
    /*
    // modification to HLL
    for (std::size_t ie = 0; ie < nInteg; ++ie)
      arrayGodFlux[it][jt][kt][ie][iDim] -= bPos*bNeg*timeStep/2*delta*eta*avgR[ie];
    */
  }

  void Gas::HllcSolver(REAL UL[], REAL UR[], REAL FL[], REAL FR[], REAL HL, REAL HR, std::size_t iDim, std::size_t it, std::size_t jt, std::size_t kt) {
    REAL avgU = (sqrt(UL[varDen])*UL[varVel[0]] + sqrt(UR[varDen])*UR[varVel[0]])/(sqrt(UL[varDen]) + sqrt(UR[varDen]));
    REAL aL = sqrt(gama*UL[varPrs]/UL[varDen]);
    REAL aR = sqrt(gama*UR[varPrs]/UR[varDen]);

    ///*
    // HLLC part
    // Pressure-based estimate by Toro
    REAL pStar=std::max(0.0, 0.5*(UL[varPrs]+UR[varPrs])-0.5*(UR[varVel[0]]-UL[varVel[0]])*0.25*(UL[varDen]+UR[varDen])*(aL+aR)); // Equ.(10.67) of Toro

    // Pressure-based estimate by TRRS
    //REAL z = 0.5*(gama-1)/gama;
    //REAL pStar = pow((aL+aR-0.5*(gama-1)*(UR[varVel[0]]-UL[varVel[0]]))/(aL/pow(UL[varPrs],z)+aR/pow(UR[varPrs],z)), 1/z); // Equ.(10.63) of Toro

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

  void Gas::UtoW() { // converting conserved variables into primitive
    negPrsDen = false;

    for (std::size_t i = 0; i < arrayU.size() && !negPrsDen; ++i)  // stop if negPrsDen
      for (std::size_t j = 0; j < arrayU[i].size() && !negPrsDen; ++j)  // stop if negPrsDen
	for (std::size_t k = 0; k <  arrayU[i][j].size() && !negPrsDen; ++k) {  // stop if negPrsDen

	  for (std::size_t m = 0; m < nDim; ++m)
	    arrayU[i][j][k][varVel[m]] = arrayU[i][j][k][varMom[m]] / arrayU[i][j][k][varDen];

	  arrayU[i][j][k][varPrs] = 0;
	  for (std::size_t m = 0; m < nDim; ++m) {
	    arrayU[i][j][k][varPrs] += pow(arrayU[i][j][k][varVel[m]],2)/2 ;
	  }
	  arrayU[i][j][k][varPrs] = (arrayU[i][j][k][varEng] - arrayU[i][j][k][varDen]*arrayU[i][j][k][varPrs]) * (gama-1);

	  if (arrayU[i][j][k][varPrs] <= 0) {
	    std::cout << std::setw(OWID) << " UtoW:prs<0";
	    negPrsDen = true;
	  }
	  if (arrayU[i][j][k][varDen] <= 0) {
	    std::cout << std::setw(OWID) << " UtoW:den<0";
	    negPrsDen = true;
	  }

	  /*
	  if (isnan(arrayU[i][j][k][varPrs])) {
	    std::cout << std::setw(OWID) << " UtoW:prs=nan" << std::setw(OWID) << i << std::setw(OWID) << j << std::setw(OWID) << k << std::endl;
	  }
	  if (isnan(arrayU[i][j][k][varDen])) {
	    std::cout << std::setw(OWID) << " UtoW:den=nan" << std::setw(OWID) << i << std::setw(OWID) << j << std::setw(OWID) << k << std::endl;
	  }
	  */
	  
	}
  }

  void Gas::WtoU() { // converting primitive variables into conserved
    // necessary to tranverse the full boundCup domain
    for (std::size_t i = 0; i < arrayU.size(); ++i)
      for (std::size_t j = 0; j < arrayU[i].size(); ++j)
	for (std::size_t k = 0; k < arrayU[i][j].size(); ++k) {
	  for (std::size_t m = 0; m < nDim; ++m)
	    arrayU[i][j][k][varMom[m]] = arrayU[i][j][k][varDen] * arrayU[i][j][k][varVel[m]];

	  arrayU[i][j][k][varEng] = 0;
	  for (std::size_t m = 0; m < nDim; ++m)
	    arrayU[i][j][k][varEng] += arrayU[i][j][k][varDen] * pow(arrayU[i][j][k][varVel[m]],2)/2 ;
	  arrayU[i][j][k][varEng] += arrayU[i][j][k][varPrs] / (gama-1);
	}
  }

  void Gas::getPtclInfo(std::vector<Particle *> &ptcls) {
    // two particles could have overlapping fluid grids in the following algorithm using maxGrid in 
    // all 3 directions, so do not clear masks inside the loops, instead, clear masks here.
    for (std::size_t i = 0; i < arrayU.size(); ++i)
      for (std::size_t j = 0; j < arrayU[i].size(); ++j)
	for (std::size_t k = 0; k < arrayU[i][j].size(); ++k)
	  arrayU[i][j][k][varMsk] = 0;

    for (std::vector<Particle*>::const_iterator it = ptcls.begin(); it != ptcls.end(); ++it)
      (*it)->clearFluidGrid();

    std::size_t maxGrid = (int) haloGrid - 1; // because haloGrid already increases by 1.
    std::size_t ip, jp, kp;

    for (std::vector<Particle*>::iterator it = ptcls.begin(); it != ptcls.end(); ++it) {
      coordToGlobalIndex((*it)->getCurrPos(), ip, jp, kp); 
      IJK global(ip, jp, kp);
      IJK local;
      bool isInHalo = globalIndexToLocal(global, local);
      
      // Note:
      // 1. local.ijk could be infinite (i.e., < 0) if printed out, which implies a particle is duplicated by adjacent processses, but is still located outside of halo zone.
      // 2. In that case, the particle has no DEM-CFD coupling effect, so it is filtered out by isInhalo condition. 
      // 3. In particular, if it is not filtered out, it causes out-of-bound memory error, of course. That is, it must be filtered out.
      // 4. However, the particle has pure DEM effect, which is accounted for by commuParticle(1) and migrateParticle(1).
      // 5. the following the std::cout line is kept for purpose of demonstrating this situation.
      //std::cout << "iter=" << iteration << " getPtclInfo: mpiRank=" << mpi.mpiRank << " ptcl=" << (*it)->getId() << " global=" << ip << " " << jp << " " << kp << " local=" << local.i << " " << local.j << " " << local.k << std::endl; 

      if (isInHalo) {
	// ensure each grid is in the valid range
	// never subtract two numbers of type std::size_t
	int lowX = std::max((int)local.i - (int)maxGrid, (int)boundCup.lowX + 1); // necessary or unnecessary? Using +/- 1 or not determines whether calcPtclForce and penalize need "if" conditions.
	int lowY = std::max((int)local.j - (int)maxGrid, (int)boundCup.lowY + 1);
	int lowZ = std::max((int)local.k - (int)maxGrid, (int)boundCup.lowZ + 1);
	int uppX = std::min((int)local.i + (int)maxGrid, (int)boundCup.uppX - 1);
	int uppY = std::min((int)local.j + (int)maxGrid, (int)boundCup.uppY - 1);
	int uppZ = std::min((int)local.k + (int)maxGrid, (int)boundCup.uppZ - 1);
	//std::cout << "x, y, z local range=" << lowX << " " << uppX << " " << lowY << " " << uppY << " " << lowZ << " " << uppZ << std::endl;

	for (std::size_t i = lowX; i <= uppX ; ++i)
	  for (std::size_t j = lowY; j <= uppY; ++j)
	    for (std::size_t k = lowZ; k <= uppZ; ++k) {
	      REAL coordX = arrayGridCoord[i][j][k][0];
	      REAL coordY = arrayGridCoord[i][j][k][1];
	      REAL coordZ = arrayGridCoord[i][j][k][2];

	      if ( (*it)->surfaceError(Vec(coordX, coordY, coordZ)) <= -EPS ) { // inside particle surface; -EPS is better than 0.
		arrayU[i][j][k][varMsk] = 1; 
		(*it)->recordFluidGrid(i, j, k);
	      }
	    }
	//std::cout << "getPtclInfo: mpiRank=" << mpi.mpiRank << " fluidGrid=" << (*it)->getFluidGrid().size() << " maxGrid= " << maxGrid << std::endl;
      }
    }
  }


  /*
  // a very inefficient algorithm
  void Gas::getPtclInfo(std::vector<Particle *> &ptcls) {
    for (std::vector<Particle*>::const_iterator it = ptcls.begin(); it != ptcls.end(); ++it)
      (*it)->clearFluidGrid();

    // 0 ~ (n-1), including boundaries
    for (std::size_t i = 0; i < arrayGridCoord.size() ; ++i)
      for (std::size_t j = 0; j < arrayGridCoord[i].size(); ++j)
	for (std::size_t k = 0; k < arrayGridCoord[i][j].size(); ++k) {

	  arrayU[i][j][k][varMsk] = 0;
	  REAL coordX = arrayGridCoord[i][j][k][0];
	  REAL coordY = arrayGridCoord[i][j][k][1];
	  REAL coordZ = arrayGridCoord[i][j][k][2];

	  for (std::vector<Particle*>::iterator it = ptcls.begin(); it != ptcls.end(); ++it)
	    if ( (*it)->surfaceError(Vec(coordX, coordY, coordZ)) <= 0 ) { // inside particle surface
	      arrayU[i][j][k][varMsk] = 1; 
	      (*it)->recordFluidGrid(i, j, k);
	      break; // save time
	    }
	}
  }
  */

  void Gas::calcPtclForce(std::vector<Particle *> &ptcls) {
    // must clear forces each loop, otherwise Gas::plot prints wrong values;
    // but Gas::penalize works OK since it uses masks.
    for (std::size_t i = 0; i < arrayPenalForce.size() ; ++i)
      for (std::size_t j = 0; j < arrayPenalForce[i].size(); ++j)
	for (std::size_t k = 0; k < arrayPenalForce[i][j].size(); ++k)
	  for (std::size_t m = 0; m < nDim; ++m) {
	    arrayPenalForce[i][j][k][m] = 0;
	    arrayPressureForce[i][j][k][m] = 0;
	  }

    for (std::vector<Particle *>::const_iterator it = ptcls.begin(); it != ptcls.end(); ++it) {
      //std::cout << "calcPtclForce: mpiRank=" << mpi.mpiRank << " ptcl=" << (*it)->getId() << std::endl; 

      REAL etaBx = 8.0/3.0 * (*it)->getA() / Cd; // local direction x (i.e. a)
      REAL etaBy = 8.0/3.0 * (*it)->getB() / Cd; // local direction y (i.e. b)
      REAL etaBz = 8.0/3.0 * (*it)->getC() / Cd; // local direction z (i.e. c)

      Vec penalForce  = 0, presForce  = 0;
      Vec penalMoment = 0, presMoment = 0;
      REAL avgDen = 0, avgVel = 0, avgPrs = 0, avgVelGap = 0;
      std::vector< std::vector<REAL> > fluidGrid = (*it)->getFluidGrid();
      //std::cout << "iter=" << iteration << " calcPtclInfo: mpiRank=" << mpi.mpiRank << " fluidGrid=" << (*it)->getFluidGrid().size() << std::endl; 

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

	REAL uxGas = arrayU[i][j][k][varVel[0]];
	REAL uyGas = arrayU[i][j][k][varVel[1]];
	REAL uzGas = arrayU[i][j][k][varVel[2]];

	Vec dist = Vec(coordX, coordY, coordZ) - (*it)->getCurrPos();
	Vec omgar = (*it)->getCurrOmga() % dist; // w X r = omga % dist, where % is overloaded as cross product

	REAL ux = (*it)->getCurrVeloc().getX() + omgar.getX(); 
	REAL uy = (*it)->getCurrVeloc().getY() + omgar.getY(); 
	REAL uz = (*it)->getCurrVeloc().getZ() + omgar.getZ();
	avgVelGap += vfabs(Vec(uxGas-ux, uyGas-uy, uzGas-uz));

	// principal axis decomposition
	// Vec globalDelta = Vec(fabs(uxGas - ux)*(uxGas - ux), fabs(uyGas - uy)*(uyGas - uy), fabs(uzGas - uz)*(uzGas - uz));
	Vec velocityGap = Vec(uxGas - ux, uyGas - uy, uzGas - uz);
	Vec globalDelta = vfabs(velocityGap) * velocityGap;

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

	// restrict pressure gradient grids; do not use (i-1) for std::size_t because (i-1) is postive when i=0
 	if (i >= boundCup.lowX+1 && i+1 <= boundCup.uppX) // necessary if getPtclInfo does not use +/-1; unnecessary if getPtclInfo uses +/-1, albeit kept here.
	  arrayPressureForce[i][j][k][0] = -(arrayU[i+1][j][k][varPrs] - arrayU[i-1][j][k][varPrs])/(2*gridDx);
	if (j >= boundCup.lowY+1 && j+1 <= boundCup.uppY)
	  arrayPressureForce[i][j][k][1] = -(arrayU[i][j+1][k][varPrs] - arrayU[i][j-1][k][varPrs])/(2*gridDy);
	if (k >= boundCup.lowZ+1 && k+1 <= boundCup.uppZ)
	  arrayPressureForce[i][j][k][2] = -(arrayU[i][j][k+1][varPrs] - arrayU[i][j][k-1][varPrs])/(2*gridDz);

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

      // print particle info
      if (!(*it)->isReceived()) { // internal particle, not received from adjacent process
	for (std::size_t iPrn = 0; iPrn < printPtcls.size(); ++iPrn) {
	  if ((*it)->getId() == printPtcls[iPrn]) {
	    std::fstream pfs;
	    pfs.open (dem::combineStr("particle_", printPtcls[iPrn], 10).c_str(), std::fstream::out | std::fstream::app);
	    if(!pfs) { std::cout << "stream error: Gas::calcPtclForce" << std::endl; exit(-1); }
	    pfs.setf(std::ios::scientific, std::ios::floatfield);
	    if (iteration == 1) {
	      pfs << std::setw(OWID) << "iteration"
		  << std::setw(OWID) << "accruedTime"
		  << std::setw(OWID) << "x"
		  << std::setw(OWID) << "y"
		  << std::setw(OWID) << "z"
		  << std::setw(OWID) << "penalFx"
		  << std::setw(OWID) << "penalFy"
		  << std::setw(OWID) << "penalFz"
		  << std::setw(OWID) << "pressureFx"
		  << std::setw(OWID) << "pressureFy"
		  << std::setw(OWID) << "pressureFz"
		  << std::setw(OWID) << "internalFx"
		  << std::setw(OWID) << "internalFy"
		  << std::setw(OWID) << "internalFz"
		  << std::setw(OWID) << "viscousCd"
		  << std::setw(OWID) << "pressureCd"
		  << std::setw(OWID) << "internalCd"
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
		//<< std::setw(OWID) << "process"
		  << std::endl;
	    }

	    // refF is only for the case of Rankine-Hugoniot Condition to test drag coefficients
	    REAL refF = 0.5*rhoL*uL*uL*dem::Pi*(*it)->getA()*(*it)->getB();
	    pfs << std::setw(OWID) << iteration
		<< std::setw(OWID) << timeAccrued
		<< std::setw(OWID) << (*it)->getCurrPos().getX()
		<< std::setw(OWID) << (*it)->getCurrPos().getY()
		<< std::setw(OWID) << (*it)->getCurrPos().getZ()

		<< std::setw(OWID) << penalForce.getX()
		<< std::setw(OWID) << penalForce.getY()
		<< std::setw(OWID) << penalForce.getZ()
		<< std::setw(OWID) << presForce.getX()
		<< std::setw(OWID) << presForce.getY()
		<< std::setw(OWID) << presForce.getZ()
		<< std::setw(OWID) << penalForce.getX()*(Cdi/Cd)
		<< std::setw(OWID) << penalForce.getY()*(Cdi/Cd)
		<< std::setw(OWID) << penalForce.getZ()*(Cdi/Cd)

		<< std::setw(OWID) << penalForce.getZ()/refF // viscousCd
		<< std::setw(OWID) << presForce.getZ()/refF  // pressureCd
		<< std::setw(OWID) << (penalForce.getZ()*(Cdi/Cd))/refF // internalCd
		<< std::setw(OWID) << ( penalForce.getZ()*(Cdi/Cd+1) + presForce.getZ() )/refF // totalCd

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
	      //<< std::setw(OWID) << mpi.mpiRank

		<< std::endl ;
	    pfs.close();

	  } // end of if ((*it)->getId() == printPtcls[iPrn])
	} // end of for (std::size_t iPrn = 0; iPrn < printPtcls.size(); ++iPrn)
      } // end of if (!(*it)->isReceived())
      
    } // end of particle loop  
  }

  void Gas::checkMomentum(std::vector<Particle *> &ptcls) {
    // A: mass and momentum of the flow field outside of the particles plus momentum of the particles as a function of time
    // B: momentum of the entire flow including inside of the particles plus momentum of the particles.
    Vec momAGas=0, momBGas=0, momPtcl=0, momA=0, momB=0;
    for (std::size_t i = 0; i < boundCup.sizeX; ++i)
      for (std::size_t j = 0; j < boundCup.sizeY; ++j)
	for (std::size_t k = 0; k < boundCup.sizeZ; ++k) {
	  momAGas += Vec(arrayU[i][j][k][varMom[0]], arrayU[i][j][k][varMom[1]], arrayU[i][j][k][varMom[2]]) * (1-arrayU[i][j][k][varMsk]); 
	  momBGas += Vec(arrayU[i][j][k][varMom[0]], arrayU[i][j][k][varMom[1]], arrayU[i][j][k][varMom[2]]);
	}
    momAGas *= gridDx*gridDy*gridDz;
    momBGas *= gridDx*gridDy*gridDz;

    for (std::vector<Particle *>::const_iterator it = ptcls.begin(); it != ptcls.end(); ++it)
      momPtcl += (*it)->getCurrVeloc() * (*it)->getMass();
    momA = momAGas + momPtcl;
    momB = momBGas + momPtcl;

    char cstr[50];
    std::fstream pfs;
    pfs.open ("momentum_progress", std::fstream::out | std::fstream::app);
    if(!pfs) { debugInf << "stream error: momentum_progress" << std::endl; exit(-1); }
    pfs.setf(std::ios::scientific, std::ios::floatfield);
    if (iteration == 1) {
      pfs << std::setw(OWID) << "iteration"
	  << std::setw(OWID) << "momAGasZ"
	  << std::setw(OWID) << "momBGasZ"
	  << std::setw(OWID) << "momPtclZ"
	  << std::setw(OWID) << "momAZ"
	  << std::setw(OWID) << "momBZ"
	  << std::setw(OWID) << "momDiffZ"
	  << std::endl;
    }

    pfs << std::setw(OWID) << iteration
	<< std::setw(OWID) << momAGas.getZ()
	<< std::setw(OWID) << momBGas.getZ()
	<< std::setw(OWID) << momPtcl.getZ()
	<< std::setw(OWID) << momA.getZ()
	<< std::setw(OWID) << momB.getZ()
	<< std::setw(OWID) << momB.getZ()-momA.getZ()
	<< std::endl;

    pfs.close();
  }

  void Gas::plot(const char *str, std::size_t snap) {

    if (mpi.mpiSize == 1) { // serial computing

      std::ofstream ofs(str);
      if(!ofs) { debugInf << "stream error: Gas::plot" << std::endl; exit(-1); }
      ofs.setf(std::ios::scientific, std::ios::floatfield);
      ofs.precision(OPREC);
    
      ofs << std::setw(OWID) << "VARIABLES = x"
	  << std::setw(OWID) << "y"
	  << std::setw(OWID) << "z"
	  << std::setw(OWID) << "Mach"
	  << std::setw(OWID) << "density"
	  << std::setw(OWID) << "momentumX"
	  << std::setw(OWID) << "momentumY"
	  << std::setw(OWID) << "momentumZ"
	  << std::setw(OWID) << "energy"
	  << std::setw(OWID) << "velocityX"
	  << std::setw(OWID) << "velocityY"
	  << std::setw(OWID) << "velocityZ"
	  << std::setw(OWID) << "pressure"
	  << std::setw(OWID) << "temperature"
	  << std::setw(OWID) << "mask"
	  << std::setw(OWID) << "penalFx"
	  << std::setw(OWID) << "penalFy"
	  << std::setw(OWID) << "penalFz"
	  << std::setw(OWID) << "pressureFx"
	  << std::setw(OWID) << "pressureFy"
	  << std::setw(OWID) << "pressureFz"
	  << std::endl;

      ofs << "ZONE I=" << std::setw(OWID) << allGridNx -2
	  << ", J="    << std::setw(OWID) << allGridNy -2
	  << ", K="    << std::setw(OWID) << allGridNz -2
	  << ", DATAPACKING=POINT, STRANDID=1, SOLUTIONTIME="
	  << std::setw(OWID) << snap
	  << std::endl;

      // write in k,j,i order, only in this order does Tecplot plot correctly (Tecplot requires I vary fastest for IJK-ordered data).
      for (std::size_t k = 1; k < allGridNz - 1; ++k)
	for (std::size_t j = 1; j < allGridNy - 1; ++j)
	  for (std::size_t i = 1; i < allGridNx -1; ++i) {
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

    } // end of serial
    else { // parallel IO

      int boundX = boundPrn.sizeX;
      int boundY = boundPrn.sizeY;
      int boundZ = boundPrn.sizeZ;
      int allBoundX, allBoundY, allBoundZ;
      MPI_Allreduce(&boundX, &allBoundX, 1, MPI_INT, MPI_SUM, mpi.mpiWorld);
      MPI_Allreduce(&boundY, &allBoundY, 1, MPI_INT, MPI_SUM, mpi.mpiWorld);
      MPI_Allreduce(&boundZ, &allBoundZ, 1, MPI_INT, MPI_SUM, mpi.mpiWorld);
      allBoundX /= (mpi.mpiProcY * mpi.mpiProcZ);
      allBoundY /= (mpi.mpiProcZ * mpi.mpiProcX);
      allBoundZ /= (mpi.mpiProcX * mpi.mpiProcY);

      MPI_Status status;
      MPI_File file;
      MPI_File_open(mpi.mpiWorld, const_cast<char *> (str), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
      if(mpi.mpiRank == 0 && !file) { debugInf << "stream error: Gas::plot" << std::endl; MPI_Abort(mpi.mpiWorld, -1);}
      std::stringstream inf;
      inf.setf(std::ios::scientific, std::ios::floatfield);

      int length = 0; // unsigned long long int length = 0; // it does not help because MPI_File_write_ordered only accepts int.

      // for MPI_File_write_ordered, each process can write different amount of data!
      if (mpi.mpiRank == 0) { // root process
	// OWID*24 + std::endl = 361
	inf << std::setw(OWID) << "VARIABLES = x"
	    << std::setw(OWID) << "y"
	    << std::setw(OWID) << "z"
	    << std::setw(OWID) << "iGlobal"
	    << std::setw(OWID) << "jGlobal"
	    << std::setw(OWID) << "kGlobal"
	    << std::setw(OWID) << "Mach"
	    << std::setw(OWID) << "density"
	    << std::setw(OWID) << "momentumX"
	    << std::setw(OWID) << "momentumY"
	    << std::setw(OWID) << "momentumZ"
	    << std::setw(OWID) << "energy"
	    << std::setw(OWID) << "velocityX"
	    << std::setw(OWID) << "velocityY"
	    << std::setw(OWID) << "velocityZ"
	    << std::setw(OWID) << "pressure"
	    << std::setw(OWID) << "temperature"
	    << std::setw(OWID) << "mask"
	    << std::setw(OWID) << "penalFx"
	    << std::setw(OWID) << "penalFy"
	    << std::setw(OWID) << "penalFz"
	    << std::setw(OWID) << "pressureFx"
	    << std::setw(OWID) << "pressureFy"
	    << std::setw(OWID) << "pressureFz"
	    << std::endl;

	std::string vac(239, ' '); // 239 = 361 -122 
	// OWID*4 + 7 + 4 + 4 + 46 + std::endl + 239 = 361
	inf << "ZONE I=" << std::setw(OWID) << allBoundX
	    << ", J="    << std::setw(OWID) << allBoundY
	    << ", K="    << std::setw(OWID) << allBoundZ
	    << ", DATAPACKING=POINT, STRANDID=1, SOLUTIONTIME="
	    << std::setw(OWID) << snap
	    << vac
	    << std::endl;
      
	length += 2; // make two records, each containing 361 chars
      }

      // (OWID*24 + std::endl) * ( boundPrn.sizeX * boundPrn.sizeY * boundPrn.sizeZ )
      for (std::size_t k = boundPrn.lowZ; k <= boundPrn.uppZ; ++k)
	for (std::size_t j = boundPrn.lowY; j <= boundPrn.uppY; ++j)
	  for (std::size_t i = boundPrn.lowX; i <= boundPrn.uppX; ++i) {
	    // i, j, k are local, mapped to global in case.
	    IJK local(i,j,k);
	    IJK global;
	    localIndexToGlobal(local, global);

	    inf << std::setw(OWID) << arrayGridCoord[i][j][k][0]
		<< std::setw(OWID) << arrayGridCoord[i][j][k][1]
		<< std::setw(OWID) << arrayGridCoord[i][j][k][2]
		<< std::setw(OWID) << global.i
		<< std::setw(OWID) << global.j
		<< std::setw(OWID) << global.k
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

      MPI_Datatype record;
      MPI_Type_contiguous(361, MPI_CHAR, &record);
      MPI_Type_commit(&record);

      length += (boundPrn.sizeX * boundPrn.sizeY * boundPrn.sizeZ); // * (OWID*24 + 1);
      //std::cout << "mpiRank=" << mpi.mpiRank << " length=" << length << std::endl;
      MPI_File_write_ordered(file, const_cast<char*> (inf.str().c_str()), length, record, &status); // in this way length is 361 times smaller
      MPI_File_close(&file);

    } // end of parallel IO

  }

  void Gas::findPrsVel(REAL UL[], REAL UR[], REAL &p, REAL &u) {
    REAL dl = UL[varDen];
    REAL dr = UR[varDen];
    REAL ul = UL[varVel[0]];
    REAL ur = UR[varVel[0]];
    REAL pl = UL[varPrs];
    REAL pr = UR[varPrs];
    REAL cl = sqrt(gama*pl/dl);
    REAL cr = sqrt(gama*pr/dr);

    const int maxIter = 20;

    // compute pressure in Star Region
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

  void Gas::guessPressure(REAL UL[], REAL UR[], REAL &pInit) {
    // provide a guessed value for pressure pInit in the Star Region. The choice is made
    // according to adaptive Riemann solver using the PVRS, TRRS and TSRS approximate
    // Riemann solvers. See Sect. 9.5 of Toro

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
    REAL pPV  = std::max(0.5*(pl + pr) + 0.5*(ul - ur)*0.25*(dl + dr)*(cl + cr), 0.0); // Equ.(4.47) of Toro.
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
      else { // select Two-Shock Riemann solver with PVRS as estimate, Equ.(9.42) of Toro.
	REAL gL = sqrt((2.0/(gama+1.0)/dl)/((gama-1.0)/(gama+1.0)*pl + pPV));
	REAL gR = sqrt((2.0/(gama+1.0)/dr)/((gama-1.0)/(gama+1.0)*pr + pPV));
	pInit = (gL*pl + gR*pr - (ur - ul))/(gL + gR);   
      }
    }
  }

  void Gas::evalF(REAL &f, REAL &fd, REAL &p, REAL &dk, REAL &pk, REAL &ck) {
    // evaluate pressure functions fl and fr in exact Riemann solver and their first derivatives
    if (p <= pk) { // rarefaction wave
      f = 2.0/(gama-1.0)*ck*(pow(p/pk, (gama-1.0)/(2.0*gama)) - 1.0);
      fd = (1.0/(dk*ck))*pow(p/pk, -(gama+1.0)/(2.0*gama));
    } else { // shock wave
      REAL ak = 2.0/(gama+1.0)/dk;
      REAL bk = (gama-1.0)/(gama+1.0)*pk;
      f = (p - pk)*sqrt(ak/(bk + p));
      fd = (1.0 - 0.5*(p - pk)/(bk + p))*sqrt(ak/(bk + p));
    }
  }

  void Gas::sampling(REAL UL[], REAL UR[], const REAL pStar, const REAL uStar, const REAL s, REAL &d, REAL &u, REAL &p) {
    // sample the solution throughout the wave pattern. Pressure pStar and velocity uStar in the star region are known. 
    // Sampling is performed in terms of the 'speed' s = x/t. Sampled values are d, u, p

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
    if (s <= uStar) { // sampling point lies to the left of the contact discontinuity
      if (pStar <= pl) { // left rarefaction
	sHL = ul - cl;
	//debugInf << "  leftrare  sHL=" << sHL;
	if (s <= sHL) { // sampled point is left data state
	  d = dl;
	  u = ul;
	  p = pl;
	} else {
	  sTL = uStar - cl*pow(pStar/pl, (gama-1.0)/(2.0*gama));
	  //debugInf << " sTL=" << sTL;
	  if (s > sTL) { // sampled point is star left state
	    d = dl*pow(pStar/pl, 1.0/gama);
	    u = uStar;
	    p = pStar;
	  } else { // sampled point is inside left fan
	    u = 2.0/(gama+1.0)*(cl + (gama-1.0)/2.0*ul + s);
	    c = 2.0/(gama+1.0)*(cl + (gama-1.0)/2.0*(ul - s));
	    d = dl*pow(c/cl, 2.0/(gama-1.0));
	    p = pl*pow(c/cl, 2.0*gama/(gama-1.0));
	  }
	  //debugInf << " ul=" << ul << " ur=" << ur << " u=" << u << " d=" << d;
	}
      } else { // left shock
	pStarL = pStar/pl;
	sL = ul - cl*sqrt((gama+1.0)/(2.0*gama)*pStarL + (gama-1.0)/(2.0*gama));
	//debugInf << "  leftshock sL=" << sL;
	if (s <= sL) { // sampled point is left data state
	  d = dl;
	  u = ul;
	  p = pl;
	} else { // sampled point is star left state
	  d = dl*(pStarL + (gama-1.0)/(gama+1.0))/(pStarL*(gama-1.0)/(gama+1.0) + 1.0);
	  u = uStar;
	  p = pStar;
	}
      }
    } else { // sampling point lies to the right of the contact discontinuity
      if (pStar > pr) { // right shock
	pStarR = pStar/pr;
	sR  = ur + cr*sqrt((gama+1.0)/(2.0*gama)*pStarR + (gama-1.0)/(2.0*gama));
	//debugInf << " rightshock sR=" << sR;
	if (s >= sR) { // sampled point is right data state
	  d = dr;
	  u = ur;
	  p = pr;
	} else { // sampled point is star right state
	  d = dr*(pStarR + (gama-1.0)/(gama+1.0))/(pStarR*(gama-1.0)/(gama+1.0) + 1.0);
	  u = uStar;
	  p = pStar;
	}
      } else { 	// right rarefaction
	sHR = ur + cr;
	//debugInf << " rightrare  sHR=" << sHR;
	if (s >= sHR) { // sampled point is right data state
	  d = dr;
	  u = ur;
	  p = pr;
	} else {
	  sTR = uStar + cr*pow(pStar/pr, (gama-1.0)/(2.0*gama));
	  if (s <= sTR) { // sampled point is star right state
	    d = dr*pow(pStar/pr, 1.0/gama);
	    u = uStar;
	    p = pStar;
	  } else { // sampled point is inside left fan
	    u = 2.0/(gama+1.0)*(-cr + (gama-1.0)/2.0*ur + s);
	    c = 2.0/(gama+1.0)*(cr - (gama-1.0)/2.0*(ur - s));
	    d = dr*pow(c/cr, 2.0/(gama-1.0));
	    p = pr*pow(c/cr, 2.0*gama/(gama-1.0));
	  }
	}
      }
    }
  }

  void Gas::findGasInRectangle(IJK &start, IJK &end, std::vector<GasVar> &foundGas) {
    for (std::size_t i = start.i; i <= end.i; ++i)
      for (std::size_t j = start.j; j <= end.j; ++j)
	for (std::size_t k = start.k; k <= end.k; ++k) {

	  REAL den = arrayU[i][j][k][varDen];
	  REAL mx  = arrayU[i][j][k][varMom[0]];
	  REAL my  = arrayU[i][j][k][varMom[1]];
	  REAL mz  = arrayU[i][j][k][varMom[2]];
	  REAL eng = arrayU[i][j][k][varEng];
	  REAL vx  = arrayU[i][j][k][varVel[0]];
	  REAL vy  = arrayU[i][j][k][varVel[1]];
	  REAL vz  = arrayU[i][j][k][varVel[2]];
	  REAL prs = arrayU[i][j][k][varPrs];
	  Vec  coords(arrayGridCoord[i][j][k][0], arrayGridCoord[i][j][k][1], arrayGridCoord[i][j][k][2]);
	  foundGas.push_back(GasVar(den, Vec(mx,my,mz), eng, Vec(vx,vy,vz), prs, coords));
	}
  }

  void Gas::findPenalInRectangle(IJK &start, IJK &end, std::vector<ValCoord> &foundPenal) {
    for (std::size_t i = start.i; i <= end.i; ++i)
      for (std::size_t j = start.j; j <= end.j; ++j)
	for (std::size_t k = start.k; k <= end.k; ++k) {
	  Vec penal(arrayPenalForce[i][j][k][0], arrayPenalForce[i][j][k][1], arrayPenalForce[i][j][k][2]);
	  Vec coords(arrayGridCoord[i][j][k][0], arrayGridCoord[i][j][k][1], arrayGridCoord[i][j][k][2]);
	  foundPenal.push_back(ValCoord(penal, coords));
	}
  }

  void Gas::commu6() {
    // if a neighbor exists (by findMPINeighbor), communicate with neighboring blocks.
    std::vector<GasVar> gasX1, gasX2;
    std::vector<GasVar> gasY1, gasY2;
    std::vector<GasVar> gasZ1, gasZ2;

    std::vector<GasVar> rGasX1, rGasX2;
    std::vector<GasVar> rGasY1, rGasY2;
    std::vector<GasVar> rGasZ1, rGasZ2;

    boost::mpi::request reqX1[2], reqX2[2];
    boost::mpi::request reqY1[2], reqY2[2];
    boost::mpi::request reqZ1[2], reqZ2[2];

    // 1. communication only occurs when there are multiple partitions in a specific direction;
    //    hereby it is safe to use haloGrid rather than haloGridX, haloGridY, haloGridZ, 
    //    because "if" condition prevents haloGrid to be zero.
    // 2. must use boundPrn for findGasInRectangle, NOT boundCup

    // 6 surfaces
    if (mpi.rankX1 >= 0) { // surface x1
      IJK start(boundPrn.lowX, boundPrn.lowY, boundPrn.lowZ);
      IJK end(boundPrn.lowX+haloGrid-1, boundPrn.uppY, boundPrn.uppZ);
      findGasInRectangle(start, end, gasX1);
      reqX1[0] = mpi.boostWorld.isend(mpi.rankX1, mpi.mpiTag,  gasX1);
      reqX1[1] = mpi.boostWorld.irecv(mpi.rankX1, mpi.mpiTag, rGasX1);
    }
    if (mpi.rankX2 >= 0) { // surface x2
      IJK start(boundPrn.uppX+1-haloGrid, boundPrn.lowY, boundPrn.lowZ);
      IJK end(boundPrn.uppX, boundPrn.uppY, boundPrn.uppZ);
      findGasInRectangle(start, end, gasX2);
      reqX2[0] = mpi.boostWorld.isend(mpi.rankX2, mpi.mpiTag,  gasX2);
      reqX2[1] = mpi.boostWorld.irecv(mpi.rankX2, mpi.mpiTag, rGasX2);
    }
    if (mpi.rankY1 >= 0) {  // surface y1
      IJK start(boundPrn.lowX, boundPrn.lowY, boundPrn.lowZ);
      IJK end(boundPrn.uppX, boundPrn.lowY+haloGrid-1, boundPrn.uppZ);
      findGasInRectangle(start, end, gasY1);
      reqY1[0] = mpi.boostWorld.isend(mpi.rankY1, mpi.mpiTag,  gasY1);
      reqY1[1] = mpi.boostWorld.irecv(mpi.rankY1, mpi.mpiTag, rGasY1);
    }
    if (mpi.rankY2 >= 0) {  // surface y2
      IJK start(boundPrn.lowX, boundPrn.uppY+1-haloGrid, boundPrn.lowZ);
      IJK end(boundPrn.uppX, boundPrn.uppY, boundPrn.uppZ);
      findGasInRectangle(start, end, gasY2);
      reqY2[0] = mpi.boostWorld.isend(mpi.rankY2, mpi.mpiTag,  gasY2);
      reqY2[1] = mpi.boostWorld.irecv(mpi.rankY2, mpi.mpiTag, rGasY2);
    }
    if (mpi.rankZ1 >= 0) {  // surface z1
      IJK start(boundPrn.lowX, boundPrn.lowY, boundPrn.lowZ);
      IJK end(boundPrn.uppX, boundPrn.uppY, boundPrn.lowZ+haloGrid-1);
      findGasInRectangle(start, end, gasZ1);
      reqZ1[0] = mpi.boostWorld.isend(mpi.rankZ1, mpi.mpiTag,  gasZ1);
      reqZ1[1] = mpi.boostWorld.irecv(mpi.rankZ1, mpi.mpiTag, rGasZ1);
    }
    if (mpi.rankZ2 >= 0) {  // surface z2
      IJK start(boundPrn.lowX, boundPrn.lowY, boundPrn.uppZ+1-haloGrid);
      IJK end(boundPrn.uppX, boundPrn.uppY, boundPrn.uppZ);
      findGasInRectangle(start, end, gasZ2);
      reqZ2[0] = mpi.boostWorld.isend(mpi.rankZ2, mpi.mpiTag,  gasZ2);
      reqZ2[1] = mpi.boostWorld.irecv(mpi.rankZ2, mpi.mpiTag, rGasZ2);
    }

    // 6 surfaces
    if (mpi.rankX1 >= 0) boost::mpi::wait_all(reqX1, reqX1 + 2);
    if (mpi.rankX2 >= 0) boost::mpi::wait_all(reqX2, reqX2 + 2);
    if (mpi.rankY1 >= 0) boost::mpi::wait_all(reqY1, reqY1 + 2);
    if (mpi.rankY2 >= 0) boost::mpi::wait_all(reqY2, reqY2 + 2);
    if (mpi.rankZ1 >= 0) boost::mpi::wait_all(reqZ1, reqZ1 + 2);
    if (mpi.rankZ2 >= 0) boost::mpi::wait_all(reqZ2, reqZ2 + 2);

    // merge using std::map
    std::vector<GasVar> recvGasVarVec;
    // 6 surfaces                                                                                          
    if (mpi.rankX1 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX1.begin(), rGasX1.end());            
    if (mpi.rankX2 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX2.begin(), rGasX2.end());            
    if (mpi.rankY1 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasY1.begin(), rGasY1.end());            
    if (mpi.rankY2 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasY2.begin(), rGasY2.end());            
    if (mpi.rankZ1 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasZ1.begin(), rGasZ1.end());            
    if (mpi.rankZ2 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasZ2.begin(), rGasZ2.end());            

    IJK_GasVar_Map recvGasMap;
    for (int i = 0; i < recvGasVarVec.size(); ++i) {
      IJK tri;
      coordToGlobalIndex(recvGasVarVec[i].coords, tri);
      recvGasMap.insert(std::make_pair(tri, recvGasVarVec[i]));
      //recvGasMap[tri] = recvGasVarVec[i]; // equivalent to map::insert.
    }
    
    // construct gas halo with received info
    for (std::size_t i = 0; i < arrayU.size(); ++i)
      for (std::size_t j = 0; j < arrayU[i].size(); ++j)
	for (std::size_t k = 0; k < arrayU[i][j].size(); ++k) {
	  if (!(i >= boundPrn.lowX && i <= boundPrn.uppX &&
		j >= boundPrn.lowY && j <= boundPrn.uppY &&
		k >= boundPrn.lowZ && k <= boundPrn.uppZ )) {

	    // i, j, k are local
	    IJK local(i,j,k);
	    IJK global;
	    localIndexToGlobal(local, global);

	    if (recvGasMap.count(global)) { // important to judge, otherwise it creates a new map elment!
	      arrayU[i][j][k][varDen]    = recvGasMap[global].density;
	      arrayU[i][j][k][varMom[0]] = recvGasMap[global].momentum.getX();
	      arrayU[i][j][k][varMom[1]] = recvGasMap[global].momentum.getY() ;
	      arrayU[i][j][k][varMom[2]] = recvGasMap[global].momentum.getZ() ;
	      arrayU[i][j][k][varEng]    = recvGasMap[global].energy;
	      arrayU[i][j][k][varVel[0]] = recvGasMap[global].velocity.getX();
	      arrayU[i][j][k][varVel[1]] = recvGasMap[global].velocity.getY() ;
	      arrayU[i][j][k][varVel[2]] = recvGasMap[global].velocity.getZ() ;
	      arrayU[i][j][k][varPrs]    = recvGasMap[global].pressure;
	    }
	  } 
	}

  } // end of Gas::commu6

  void Gas::commu26() {
    // if a neighbor exists (by findMPINeighbor), communicate with neighboring blocks.
    std::vector<GasVar> gasX1, gasX2;
    std::vector<GasVar> gasY1, gasY2;
    std::vector<GasVar> gasZ1, gasZ2;
    std::vector<GasVar> gasX1Y1, gasX1Y2, gasX1Z1, gasX1Z2; 
    std::vector<GasVar> gasX2Y1, gasX2Y2, gasX2Z1, gasX2Z2; 
    std::vector<GasVar> gasY1Z1, gasY1Z2, gasY2Z1, gasY2Z2; 
    std::vector<GasVar> gasX1Y1Z1, gasX1Y1Z2, gasX1Y2Z1, gasX1Y2Z2; 
    std::vector<GasVar> gasX2Y1Z1, gasX2Y1Z2, gasX2Y2Z1, gasX2Y2Z2; 

    std::vector<GasVar> rGasX1, rGasX2;
    std::vector<GasVar> rGasY1, rGasY2;
    std::vector<GasVar> rGasZ1, rGasZ2;
    std::vector<GasVar> rGasX1Y1, rGasX1Y2, rGasX1Z1, rGasX1Z2; 
    std::vector<GasVar> rGasX2Y1, rGasX2Y2, rGasX2Z1, rGasX2Z2; 
    std::vector<GasVar> rGasY1Z1, rGasY1Z2, rGasY2Z1, rGasY2Z2; 
    std::vector<GasVar> rGasX1Y1Z1, rGasX1Y1Z2, rGasX1Y2Z1, rGasX1Y2Z2; 
    std::vector<GasVar> rGasX2Y1Z1, rGasX2Y1Z2, rGasX2Y2Z1, rGasX2Y2Z2; 

    boost::mpi::request reqX1[2], reqX2[2];
    boost::mpi::request reqY1[2], reqY2[2];
    boost::mpi::request reqZ1[2], reqZ2[2];
    boost::mpi::request reqX1Y1[2], reqX1Y2[2], reqX1Z1[2], reqX1Z2[2];
    boost::mpi::request reqX2Y1[2], reqX2Y2[2], reqX2Z1[2], reqX2Z2[2];
    boost::mpi::request reqY1Z1[2], reqY1Z2[2], reqY2Z1[2], reqY2Z2[2];
    boost::mpi::request reqX1Y1Z1[2], reqX1Y1Z2[2], reqX1Y2Z1[2], reqX1Y2Z2[2];
    boost::mpi::request reqX2Y1Z1[2], reqX2Y1Z2[2], reqX2Y2Z1[2], reqX2Y2Z2[2];

    // 1. communication only occurs when there are multiple partitions in a specific direction;
    //    hereby it is safe to use haloGrid rather than haloGridX, haloGridY, haloGridZ, 
    //    because "if" condition prevents haloGrid to be zero.
    // 2. must use boundPrn for findGasInRectangle, NOT boundCup

    // 6 surfaces
    if (mpi.rankX1 >= 0) { // surface x1
      IJK start(boundPrn.lowX, boundPrn.lowY, boundPrn.lowZ);
      IJK end(boundPrn.lowX+haloGrid-1, boundPrn.uppY, boundPrn.uppZ);
      findGasInRectangle(start, end, gasX1);
      reqX1[0] = mpi.boostWorld.isend(mpi.rankX1, mpi.mpiTag,  gasX1);
      reqX1[1] = mpi.boostWorld.irecv(mpi.rankX1, mpi.mpiTag, rGasX1);
    }
    if (mpi.rankX2 >= 0) { // surface x2
      IJK start(boundPrn.uppX+1-haloGrid, boundPrn.lowY, boundPrn.lowZ);
      IJK end(boundPrn.uppX, boundPrn.uppY, boundPrn.uppZ);
      findGasInRectangle(start, end, gasX2);
      reqX2[0] = mpi.boostWorld.isend(mpi.rankX2, mpi.mpiTag,  gasX2);
      reqX2[1] = mpi.boostWorld.irecv(mpi.rankX2, mpi.mpiTag, rGasX2);
    }
    if (mpi.rankY1 >= 0) {  // surface y1
      IJK start(boundPrn.lowX, boundPrn.lowY, boundPrn.lowZ);
      IJK end(boundPrn.uppX, boundPrn.lowY+haloGrid-1, boundPrn.uppZ);
      findGasInRectangle(start, end, gasY1);
      reqY1[0] = mpi.boostWorld.isend(mpi.rankY1, mpi.mpiTag,  gasY1);
      reqY1[1] = mpi.boostWorld.irecv(mpi.rankY1, mpi.mpiTag, rGasY1);
    }
    if (mpi.rankY2 >= 0) {  // surface y2
      IJK start(boundPrn.lowX, boundPrn.uppY+1-haloGrid, boundPrn.lowZ);
      IJK end(boundPrn.uppX, boundPrn.uppY, boundPrn.uppZ);
      findGasInRectangle(start, end, gasY2);
      reqY2[0] = mpi.boostWorld.isend(mpi.rankY2, mpi.mpiTag,  gasY2);
      reqY2[1] = mpi.boostWorld.irecv(mpi.rankY2, mpi.mpiTag, rGasY2);
    }
    if (mpi.rankZ1 >= 0) {  // surface z1
      IJK start(boundPrn.lowX, boundPrn.lowY, boundPrn.lowZ);
      IJK end(boundPrn.uppX, boundPrn.uppY, boundPrn.lowZ+haloGrid-1);
      findGasInRectangle(start, end, gasZ1);
      reqZ1[0] = mpi.boostWorld.isend(mpi.rankZ1, mpi.mpiTag,  gasZ1);
      reqZ1[1] = mpi.boostWorld.irecv(mpi.rankZ1, mpi.mpiTag, rGasZ1);
    }
    if (mpi.rankZ2 >= 0) {  // surface z2
      IJK start(boundPrn.lowX, boundPrn.lowY, boundPrn.uppZ+1-haloGrid);
      IJK end(boundPrn.uppX, boundPrn.uppY, boundPrn.uppZ);
      findGasInRectangle(start, end, gasZ2);
      reqZ2[0] = mpi.boostWorld.isend(mpi.rankZ2, mpi.mpiTag,  gasZ2);
      reqZ2[1] = mpi.boostWorld.irecv(mpi.rankZ2, mpi.mpiTag, rGasZ2);
    }
    // 12 edges
    if (mpi.rankX1Y1 >= 0) { // edge x1y1
      IJK start(boundPrn.lowX, boundPrn.lowY, boundPrn.lowZ);
      IJK end(boundPrn.lowX+haloGrid-1, boundPrn.lowY+haloGrid-1, boundPrn.uppZ);
      findGasInRectangle(start, end, gasX1Y1);
      reqX1Y1[0] = mpi.boostWorld.isend(mpi.rankX1Y1, mpi.mpiTag,  gasX1Y1);
      reqX1Y1[1] = mpi.boostWorld.irecv(mpi.rankX1Y1, mpi.mpiTag, rGasX1Y1);
    }
    if (mpi.rankX1Y2 >= 0) { // edge x1y2
      IJK start(boundPrn.lowX, boundPrn.uppY+1-haloGrid, boundPrn.lowZ);
      IJK end(boundPrn.lowX+haloGrid-1, boundPrn.uppY, boundPrn.uppZ);
      findGasInRectangle(start, end, gasX1Y2);
      reqX1Y2[0] = mpi.boostWorld.isend(mpi.rankX1Y2, mpi.mpiTag,  gasX1Y2);
      reqX1Y2[1] = mpi.boostWorld.irecv(mpi.rankX1Y2, mpi.mpiTag, rGasX1Y2);
    }
    if (mpi.rankX1Z1 >= 0) { // edge x1z1
      IJK start(boundPrn.lowX, boundPrn.lowY, boundPrn.lowZ);
      IJK end(boundPrn.lowX+haloGrid-1, boundPrn.uppY, boundPrn.lowZ+haloGrid-1);
      findGasInRectangle(start, end, gasX1Z1);
      reqX1Z1[0] = mpi.boostWorld.isend(mpi.rankX1Z1, mpi.mpiTag,  gasX1Z1);
      reqX1Z1[1] = mpi.boostWorld.irecv(mpi.rankX1Z1, mpi.mpiTag, rGasX1Z1);
    }
    if (mpi.rankX1Z2 >= 0) { // edge x1z2
      IJK start(boundPrn.lowX, boundPrn.lowY, boundPrn.uppZ+1-haloGrid);
      IJK end(boundPrn.lowX+haloGrid-1, boundPrn.uppY, boundPrn.uppZ);
      findGasInRectangle(start, end, gasX1Z2);
      reqX1Z2[0] = mpi.boostWorld.isend(mpi.rankX1Z2, mpi.mpiTag,  gasX1Z2);
      reqX1Z2[1] = mpi.boostWorld.irecv(mpi.rankX1Z2, mpi.mpiTag, rGasX1Z2);
    }
    if (mpi.rankX2Y1 >= 0) { // edge x2y1
      IJK start(boundPrn.uppX+1-haloGrid, boundPrn.lowY, boundPrn.lowZ);
      IJK end(boundPrn.uppX, boundPrn.lowY+haloGrid-1, boundPrn.uppZ);
      findGasInRectangle(start, end, gasX2Y1);
      reqX2Y1[0] = mpi.boostWorld.isend(mpi.rankX2Y1, mpi.mpiTag,  gasX2Y1);
      reqX2Y1[1] = mpi.boostWorld.irecv(mpi.rankX2Y1, mpi.mpiTag, rGasX2Y1);
    }
    if (mpi.rankX2Y2 >= 0) { // edge x2y2
      IJK start(boundPrn.uppX+1-haloGrid, boundPrn.uppX+1-haloGrid, boundPrn.lowZ);
      IJK end(boundPrn.uppX, boundPrn.uppY, boundPrn.uppZ);
      findGasInRectangle(start, end, gasX2Y2);
      reqX2Y2[0] = mpi.boostWorld.isend(mpi.rankX2Y2, mpi.mpiTag,  gasX2Y2);
      reqX2Y2[1] = mpi.boostWorld.irecv(mpi.rankX2Y2, mpi.mpiTag, rGasX2Y2);
    }
    if (mpi.rankX2Z1 >= 0) { // edge x2z1
      IJK start(boundPrn.uppX+1-haloGrid, boundPrn.lowY, boundPrn.lowZ);
      IJK end(boundPrn.uppX, boundPrn.uppY, boundPrn.lowZ+haloGrid-1);
      findGasInRectangle(start, end, gasX2Z1);
      reqX2Z1[0] = mpi.boostWorld.isend(mpi.rankX2Z1, mpi.mpiTag,  gasX2Z1);
      reqX2Z1[1] = mpi.boostWorld.irecv(mpi.rankX2Z1, mpi.mpiTag, rGasX2Z1);
    }
    if (mpi.rankX2Z2 >= 0) { // edge x2z2
      IJK start(boundPrn.uppX+1-haloGrid, boundPrn.lowY, boundPrn.uppZ+1-haloGrid);
      IJK end(boundPrn.uppX, boundPrn.uppY, boundPrn.uppZ);
      findGasInRectangle(start, end, gasX2Z2);
      reqX2Z2[0] = mpi.boostWorld.isend(mpi.rankX2Z2, mpi.mpiTag,  gasX2Z2);
      reqX2Z2[1] = mpi.boostWorld.irecv(mpi.rankX2Z2, mpi.mpiTag, rGasX2Z2);
    }
    if (mpi.rankY1Z1 >= 0) { // edge y1z1
      IJK start(boundPrn.lowX, boundPrn.lowY, boundPrn.lowZ);
      IJK end(boundPrn.uppX, boundPrn.lowY+haloGrid-1,boundPrn.lowZ+ haloGrid-1);
      findGasInRectangle(start, end, gasY1Z1);
      reqY1Z1[0] = mpi.boostWorld.isend(mpi.rankY1Z1, mpi.mpiTag,  gasY1Z1);
      reqY1Z1[1] = mpi.boostWorld.irecv(mpi.rankY1Z1, mpi.mpiTag, rGasY1Z1);
    }
    if (mpi.rankY1Z2 >= 0) { // edge y1z2
      IJK start(boundPrn.lowX, boundPrn.lowY, boundPrn.uppZ+1-haloGrid);
      IJK end(boundPrn.uppX, boundPrn.lowY+haloGrid-1, boundPrn.uppZ);
      findGasInRectangle(start, end, gasY1Z2);
      reqY1Z2[0] = mpi.boostWorld.isend(mpi.rankY1Z2, mpi.mpiTag,  gasY1Z2);
      reqY1Z2[1] = mpi.boostWorld.irecv(mpi.rankY1Z2, mpi.mpiTag, rGasY1Z2);
    }
    if (mpi.rankY2Z1 >= 0) { // edge y2z1
      IJK start(boundPrn.lowX, boundPrn.uppY+1-haloGrid, boundPrn.lowZ);
      IJK end(boundPrn.uppX, boundPrn.uppY, boundPrn.lowZ+haloGrid-1);
      findGasInRectangle(start, end, gasY2Z1);
      reqY2Z1[0] = mpi.boostWorld.isend(mpi.rankY2Z1, mpi.mpiTag,  gasY2Z1);
      reqY2Z1[1] = mpi.boostWorld.irecv(mpi.rankY2Z1, mpi.mpiTag, rGasY2Z1);
    }
    if (mpi.rankY2Z2 >= 0) { // edge y2z2
      IJK start(boundPrn.lowX, boundPrn.uppY+1-haloGrid, boundPrn.uppZ+1-haloGrid);
      IJK end(boundPrn.uppX, boundPrn.uppY, boundPrn.uppZ);
      findGasInRectangle(start, end, gasY2Z2);
      reqY2Z2[0] = mpi.boostWorld.isend(mpi.rankY2Z2, mpi.mpiTag,  gasY2Z2);
      reqY2Z2[1] = mpi.boostWorld.irecv(mpi.rankY2Z2, mpi.mpiTag, rGasY2Z2);
    }
    // 8 vertices
    if (mpi.rankX1Y1Z1 >= 0) { // vertice x1y1z1
      IJK start(boundPrn.lowX, boundPrn.lowY, boundPrn.lowZ);
      IJK end(boundPrn.lowX+haloGrid-1, boundPrn.lowY+haloGrid-1, boundPrn.lowZ+haloGrid-1);
      findGasInRectangle(start, end, gasX1Y1Z1);
      reqX1Y1Z1[0] = mpi.boostWorld.isend(mpi.rankX1Y1Z1, mpi.mpiTag,  gasX1Y1Z1);
      reqX1Y1Z1[1] = mpi.boostWorld.irecv(mpi.rankX1Y1Z1, mpi.mpiTag, rGasX1Y1Z1);
    }
    if (mpi.rankX1Y1Z2 >= 0) { // vertice x1y1z2
      IJK start(boundPrn.lowX, boundPrn.lowY,  boundPrn.uppZ+1-haloGrid);
      IJK end(boundPrn.lowX+haloGrid-1, boundPrn.lowY+haloGrid-1, boundPrn.uppZ);
      findGasInRectangle(start, end, gasX1Y1Z2);
      reqX1Y1Z2[0] = mpi.boostWorld.isend(mpi.rankX1Y1Z2, mpi.mpiTag,  gasX1Y1Z2);
      reqX1Y1Z2[1] = mpi.boostWorld.irecv(mpi.rankX1Y1Z2, mpi.mpiTag, rGasX1Y1Z2);
    }
    if (mpi.rankX1Y2Z1 >= 0) { // vertice x1y2z1
      IJK start(boundPrn.lowX, boundPrn.uppY+1-haloGrid, boundPrn.lowZ);
      IJK end(boundPrn.lowX+haloGrid-1, boundPrn.uppY, boundPrn.lowZ+haloGrid-1);
      findGasInRectangle(start, end, gasX1Y2Z1);
      reqX1Y2Z1[0] = mpi.boostWorld.isend(mpi.rankX1Y2Z1, mpi.mpiTag,  gasX1Y2Z1);
      reqX1Y2Z1[1] = mpi.boostWorld.irecv(mpi.rankX1Y2Z1, mpi.mpiTag, rGasX1Y2Z1);
    }
    if (mpi.rankX1Y2Z2 >= 0) { // vertice x1y2z2
      IJK start(boundPrn.lowX, boundPrn.uppY+1-haloGrid, boundPrn.uppZ+1-haloGrid);
      IJK end(boundPrn.lowX+haloGrid-1, boundPrn.uppY, boundPrn.uppZ);
      findGasInRectangle(start, end, gasX1Y2Z2);
      reqX1Y2Z2[0] = mpi.boostWorld.isend(mpi.rankX1Y2Z2, mpi.mpiTag,  gasX1Y2Z2);
      reqX1Y2Z2[1] = mpi.boostWorld.irecv(mpi.rankX1Y2Z2, mpi.mpiTag, rGasX1Y2Z2);
    }
    if (mpi.rankX2Y1Z1 >= 0) { // vertice x2y1z1
      IJK start(boundPrn.uppX+1-haloGrid, boundPrn.lowY, boundPrn.lowZ);
      IJK end(boundPrn.uppX, boundPrn.lowY+haloGrid-1, boundPrn.lowZ+haloGrid-1);
      findGasInRectangle(start, end, gasX2Y1Z1);
      reqX2Y1Z1[0] = mpi.boostWorld.isend(mpi.rankX2Y1Z1, mpi.mpiTag,  gasX2Y1Z1);
      reqX2Y1Z1[1] = mpi.boostWorld.irecv(mpi.rankX2Y1Z1, mpi.mpiTag, rGasX2Y1Z1);
    }
    if (mpi.rankX2Y1Z2 >= 0) { // vertice x2y1z2
      IJK start(boundPrn.uppX+1-haloGrid, boundPrn.lowY, boundPrn.uppZ+1-haloGrid);
      IJK end(boundPrn.uppX, boundPrn.lowY+haloGrid-1, boundPrn.uppZ);
      findGasInRectangle(start, end, gasX2Y1Z2);
      reqX2Y1Z2[0] = mpi.boostWorld.isend(mpi.rankX2Y1Z2, mpi.mpiTag,  gasX2Y1Z2);
      reqX2Y1Z2[1] = mpi.boostWorld.irecv(mpi.rankX2Y1Z2, mpi.mpiTag, rGasX2Y1Z2);
    }
    if (mpi.rankX2Y2Z1 >= 0) { // vertice x2y2z1
      IJK start(boundPrn.uppX+1-haloGrid, boundPrn.uppY+1-haloGrid, boundPrn.lowZ);
      IJK end(boundPrn.uppX, boundPrn.uppY, boundPrn.lowZ+haloGrid-1);
      findGasInRectangle(start, end, gasX2Y2Z1);
      reqX2Y2Z1[0] = mpi.boostWorld.isend(mpi.rankX2Y2Z1, mpi.mpiTag,  gasX2Y2Z1);
      reqX2Y2Z1[1] = mpi.boostWorld.irecv(mpi.rankX2Y2Z1, mpi.mpiTag, rGasX2Y2Z1);
    }
    if (mpi.rankX2Y2Z2 >= 0) { // vertice x2y2z2
      IJK start(boundPrn.uppX+1-haloGrid, boundPrn.uppY+1-haloGrid, boundPrn.uppZ+1-haloGrid);
      IJK end(boundPrn.uppX, boundPrn.uppY, boundPrn.uppZ);
      findGasInRectangle(start, end, gasX2Y2Z2);
      reqX2Y2Z2[0] = mpi.boostWorld.isend(mpi.rankX2Y2Z2, mpi.mpiTag,  gasX2Y2Z2);
      reqX2Y2Z2[1] = mpi.boostWorld.irecv(mpi.rankX2Y2Z2, mpi.mpiTag, rGasX2Y2Z2);
    }

    // 6 surfaces
    if (mpi.rankX1 >= 0) boost::mpi::wait_all(reqX1, reqX1 + 2);
    if (mpi.rankX2 >= 0) boost::mpi::wait_all(reqX2, reqX2 + 2);
    if (mpi.rankY1 >= 0) boost::mpi::wait_all(reqY1, reqY1 + 2);
    if (mpi.rankY2 >= 0) boost::mpi::wait_all(reqY2, reqY2 + 2);
    if (mpi.rankZ1 >= 0) boost::mpi::wait_all(reqZ1, reqZ1 + 2);
    if (mpi.rankZ2 >= 0) boost::mpi::wait_all(reqZ2, reqZ2 + 2);
    // 12 edges
    if (mpi.rankX1Y1 >= 0) boost::mpi::wait_all(reqX1Y1, reqX1Y1 + 2);
    if (mpi.rankX1Y2 >= 0) boost::mpi::wait_all(reqX1Y2, reqX1Y2 + 2);  
    if (mpi.rankX1Z1 >= 0) boost::mpi::wait_all(reqX1Z1, reqX1Z1 + 2);
    if (mpi.rankX1Z2 >= 0) boost::mpi::wait_all(reqX1Z2, reqX1Z2 + 2);
    if (mpi.rankX2Y1 >= 0) boost::mpi::wait_all(reqX2Y1, reqX2Y1 + 2);
    if (mpi.rankX2Y2 >= 0) boost::mpi::wait_all(reqX2Y2, reqX2Y2 + 2);  
    if (mpi.rankX2Z1 >= 0) boost::mpi::wait_all(reqX2Z1, reqX2Z1 + 2);
    if (mpi.rankX2Z2 >= 0) boost::mpi::wait_all(reqX2Z2, reqX2Z2 + 2); 
    if (mpi.rankY1Z1 >= 0) boost::mpi::wait_all(reqY1Z1, reqY1Z1 + 2);
    if (mpi.rankY1Z2 >= 0) boost::mpi::wait_all(reqY1Z2, reqY1Z2 + 2);
    if (mpi.rankY2Z1 >= 0) boost::mpi::wait_all(reqY2Z1, reqY2Z1 + 2);
    if (mpi.rankY2Z2 >= 0) boost::mpi::wait_all(reqY2Z2, reqY2Z2 + 2); 
    // 8 vertices
    if (mpi.rankX1Y1Z1 >= 0) boost::mpi::wait_all(reqX1Y1Z1, reqX1Y1Z1 + 2);
    if (mpi.rankX1Y1Z2 >= 0) boost::mpi::wait_all(reqX1Y1Z2, reqX1Y1Z2 + 2);
    if (mpi.rankX1Y2Z1 >= 0) boost::mpi::wait_all(reqX1Y2Z1, reqX1Y2Z1 + 2);
    if (mpi.rankX1Y2Z2 >= 0) boost::mpi::wait_all(reqX1Y2Z2, reqX1Y2Z2 + 2);
    if (mpi.rankX2Y1Z1 >= 0) boost::mpi::wait_all(reqX2Y1Z1, reqX2Y1Z1 + 2);
    if (mpi.rankX2Y1Z2 >= 0) boost::mpi::wait_all(reqX2Y1Z2, reqX2Y1Z2 + 2);
    if (mpi.rankX2Y2Z1 >= 0) boost::mpi::wait_all(reqX2Y2Z1, reqX2Y2Z1 + 2);
    if (mpi.rankX2Y2Z2 >= 0) boost::mpi::wait_all(reqX2Y2Z2, reqX2Y2Z2 + 2);

    // merge using std::map
    std::vector<GasVar> recvGasVarVec;
    // 6 surfaces                                                                                          
    if (mpi.rankX1 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX1.begin(), rGasX1.end());            
    if (mpi.rankX2 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX2.begin(), rGasX2.end());            
    if (mpi.rankY1 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasY1.begin(), rGasY1.end());            
    if (mpi.rankY2 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasY2.begin(), rGasY2.end());            
    if (mpi.rankZ1 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasZ1.begin(), rGasZ1.end());            
    if (mpi.rankZ2 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasZ2.begin(), rGasZ2.end());            
    // 12 edges                                                                                            
    if (mpi.rankX1Y1 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX1Y1.begin(), rGasX1Y1.end());      
    if (mpi.rankX1Y2 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX1Y2.begin(), rGasX1Y2.end());      
    if (mpi.rankX1Z1 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX1Z1.begin(), rGasX1Z1.end());      
    if (mpi.rankX1Z2 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX1Z2.begin(), rGasX1Z2.end());      
    if (mpi.rankX2Y1 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX2Y1.begin(), rGasX2Y1.end());      
    if (mpi.rankX2Y2 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX2Y2.begin(), rGasX2Y2.end());      
    if (mpi.rankX2Z1 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX2Z1.begin(), rGasX2Z1.end());      
    if (mpi.rankX2Z2 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX2Z2.begin(), rGasX2Z2.end());      
    if (mpi.rankY1Z1 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasY1Z1.begin(), rGasY1Z1.end());      
    if (mpi.rankY1Z2 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasY1Z2.begin(), rGasY1Z2.end());      
    if (mpi.rankY2Z1 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasY2Z1.begin(), rGasY2Z1.end());      
    if (mpi.rankY2Z2 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasY2Z2.begin(), rGasY2Z2.end());      
    // 8 vertices                                                                                          
    if (mpi.rankX1Y1Z1 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX1Y1Z1.begin(), rGasX1Y1Z1.end());
    if (mpi.rankX1Y1Z2 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX1Y1Z2.begin(), rGasX1Y1Z2.end());
    if (mpi.rankX1Y2Z1 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX1Y2Z1.begin(), rGasX1Y2Z1.end());
    if (mpi.rankX1Y2Z2 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX1Y2Z2.begin(), rGasX1Y2Z2.end());
    if (mpi.rankX2Y1Z1 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX2Y1Z1.begin(), rGasX2Y1Z1.end());
    if (mpi.rankX2Y1Z2 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX2Y1Z2.begin(), rGasX2Y1Z2.end());
    if (mpi.rankX2Y2Z1 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX2Y2Z1.begin(), rGasX2Y2Z1.end());
    if (mpi.rankX2Y2Z2 >= 0) recvGasVarVec.insert(recvGasVarVec.end(), rGasX2Y2Z2.begin(), rGasX2Y2Z2.end());

    IJK_GasVar_Map recvGasMap;
    for (int i = 0; i < recvGasVarVec.size(); ++i) {
      IJK tri;
      coordToGlobalIndex(recvGasVarVec[i].coords, tri);
      recvGasMap[tri] = recvGasVarVec[i];
    }
    
    // construct gas halo with received info
    for (std::size_t i = 0; i < arrayU.size(); ++i)
      for (std::size_t j = 0; j < arrayU[i].size(); ++j)
	for (std::size_t k = 0; k < arrayU[i][j].size(); ++k) {
	  if (!(i >= boundPrn.lowX && i <= boundPrn.uppX &&
		j >= boundPrn.lowY && j <= boundPrn.uppY &&
		k >= boundPrn.lowZ && k <= boundPrn.uppZ )) {

	    // i, j, k are local
	    IJK local(i,j,k);
	    IJK global;
	    localIndexToGlobal(local, global);

	    if (recvGasMap.count(global)) { // important to judge, otherwise it creates a new map elment!
	      arrayU[i][j][k][varDen]    = recvGasMap[global].density;
	      arrayU[i][j][k][varMom[0]] = recvGasMap[global].momentum.getX();
	      arrayU[i][j][k][varMom[1]] = recvGasMap[global].momentum.getY() ;
	      arrayU[i][j][k][varMom[2]] = recvGasMap[global].momentum.getZ() ;
	      arrayU[i][j][k][varEng]    = recvGasMap[global].energy;
	      arrayU[i][j][k][varVel[0]] = recvGasMap[global].velocity.getX();
	      arrayU[i][j][k][varVel[1]] = recvGasMap[global].velocity.getY() ;
	      arrayU[i][j][k][varVel[2]] = recvGasMap[global].velocity.getZ() ;
	      arrayU[i][j][k][varPrs]    = recvGasMap[global].pressure;
	    }
	  } 
	}

  } // end of Gas::commu26

  // This function is useless because inverse communication is not needed in the algorithm.
  // Leave it here for demonstration.
  void Gas::backCommu26() {
    // if a neighbor exists (by findMPINeighbor), communicate with neighboring blocks.

    // in this member function, gas and rGas are used to store arrayPenalForce.
    std::vector<ValCoord> gasX1, gasX2;
    std::vector<ValCoord> gasY1, gasY2;
    std::vector<ValCoord> gasZ1, gasZ2;
    std::vector<ValCoord> gasX1Y1, gasX1Y2, gasX1Z1, gasX1Z2; 
    std::vector<ValCoord> gasX2Y1, gasX2Y2, gasX2Z1, gasX2Z2; 
    std::vector<ValCoord> gasY1Z1, gasY1Z2, gasY2Z1, gasY2Z2; 
    std::vector<ValCoord> gasX1Y1Z1, gasX1Y1Z2, gasX1Y2Z1, gasX1Y2Z2; 
    std::vector<ValCoord> gasX2Y1Z1, gasX2Y1Z2, gasX2Y2Z1, gasX2Y2Z2; 

    std::vector<ValCoord> rGasX1, rGasX2;
    std::vector<ValCoord> rGasY1, rGasY2;
    std::vector<ValCoord> rGasZ1, rGasZ2;
    std::vector<ValCoord> rGasX1Y1, rGasX1Y2, rGasX1Z1, rGasX1Z2; 
    std::vector<ValCoord> rGasX2Y1, rGasX2Y2, rGasX2Z1, rGasX2Z2; 
    std::vector<ValCoord> rGasY1Z1, rGasY1Z2, rGasY2Z1, rGasY2Z2; 
    std::vector<ValCoord> rGasX1Y1Z1, rGasX1Y1Z2, rGasX1Y2Z1, rGasX1Y2Z2; 
    std::vector<ValCoord> rGasX2Y1Z1, rGasX2Y1Z2, rGasX2Y2Z1, rGasX2Y2Z2; 

    boost::mpi::request reqX1[2], reqX2[2];
    boost::mpi::request reqY1[2], reqY2[2];
    boost::mpi::request reqZ1[2], reqZ2[2];
    boost::mpi::request reqX1Y1[2], reqX1Y2[2], reqX1Z1[2], reqX1Z2[2];
    boost::mpi::request reqX2Y1[2], reqX2Y2[2], reqX2Z1[2], reqX2Z2[2];
    boost::mpi::request reqY1Z1[2], reqY1Z2[2], reqY2Z1[2], reqY2Z2[2];
    boost::mpi::request reqX1Y1Z1[2], reqX1Y1Z2[2], reqX1Y2Z1[2], reqX1Y2Z2[2];
    boost::mpi::request reqX2Y1Z1[2], reqX2Y1Z2[2], reqX2Y2Z1[2], reqX2Y2Z2[2];

    // 1. communication only occurs when there are multiple partitions in a specific direction;
    //    hereby it is safe to use haloGrid rather than haloGridX, haloGridY, haloGridZ, 
    //    because "if" condition prevents haloGrid to be zero.
    // 2. must use boundCup for findPenalInRectangle, NOT boundCup

    // 6 surfaces
    if (mpi.rankX1 >= 0) { // surface x1
      IJK start(boundCup.lowX, boundCup.lowY, boundCup.lowZ);
      IJK end(boundCup.lowX+haloGrid-1, boundCup.uppY, boundCup.uppZ);
      findPenalInRectangle(start, end, gasX1);
      reqX1[0] = mpi.boostWorld.isend(mpi.rankX1, mpi.mpiTag,  gasX1);
      reqX1[1] = mpi.boostWorld.irecv(mpi.rankX1, mpi.mpiTag, rGasX1);
    }
    if (mpi.rankX2 >= 0) { // surface x2
      IJK start(boundCup.uppX+1-haloGrid, boundCup.lowY, boundCup.lowZ);
      IJK end(boundCup.uppX, boundCup.uppY, boundCup.uppZ);
      findPenalInRectangle(start, end, gasX2);
      reqX2[0] = mpi.boostWorld.isend(mpi.rankX2, mpi.mpiTag,  gasX2);
      reqX2[1] = mpi.boostWorld.irecv(mpi.rankX2, mpi.mpiTag, rGasX2);
    }
    if (mpi.rankY1 >= 0) {  // surface y1
      IJK start(boundCup.lowX, boundCup.lowY, boundCup.lowZ);
      IJK end(boundCup.uppX, boundCup.lowY+haloGrid-1, boundCup.uppZ);
      findPenalInRectangle(start, end, gasY1);
      reqY1[0] = mpi.boostWorld.isend(mpi.rankY1, mpi.mpiTag,  gasY1);
      reqY1[1] = mpi.boostWorld.irecv(mpi.rankY1, mpi.mpiTag, rGasY1);
    }
    if (mpi.rankY2 >= 0) {  // surface y2
      IJK start(boundCup.lowX, boundCup.uppY+1-haloGrid, boundCup.lowZ);
      IJK end(boundCup.uppX, boundCup.uppY, boundCup.uppZ);
      findPenalInRectangle(start, end, gasY2);
      reqY2[0] = mpi.boostWorld.isend(mpi.rankY2, mpi.mpiTag,  gasY2);
      reqY2[1] = mpi.boostWorld.irecv(mpi.rankY2, mpi.mpiTag, rGasY2);
    }
    if (mpi.rankZ1 >= 0) {  // surface z1
      IJK start(boundCup.lowX, boundCup.lowY, boundCup.lowZ);
      IJK end(boundCup.uppX, boundCup.uppY, boundCup.lowZ+haloGrid-1);
      findPenalInRectangle(start, end, gasZ1);
      reqZ1[0] = mpi.boostWorld.isend(mpi.rankZ1, mpi.mpiTag,  gasZ1);
      reqZ1[1] = mpi.boostWorld.irecv(mpi.rankZ1, mpi.mpiTag, rGasZ1);
    }
    if (mpi.rankZ2 >= 0) {  // surface z2
      IJK start(boundCup.lowX, boundCup.lowY, boundCup.uppZ+1-haloGrid);
      IJK end(boundCup.uppX, boundCup.uppY, boundCup.uppZ);
      findPenalInRectangle(start, end, gasZ2);
      reqZ2[0] = mpi.boostWorld.isend(mpi.rankZ2, mpi.mpiTag,  gasZ2);
      reqZ2[1] = mpi.boostWorld.irecv(mpi.rankZ2, mpi.mpiTag, rGasZ2);
    }
    // 12 edges
    if (mpi.rankX1Y1 >= 0) { // edge x1y1
      IJK start(boundCup.lowX, boundCup.lowY, boundCup.lowZ);
      IJK end(boundCup.lowX+haloGrid-1, boundCup.lowY+haloGrid-1, boundCup.uppZ);
      findPenalInRectangle(start, end, gasX1Y1);
      reqX1Y1[0] = mpi.boostWorld.isend(mpi.rankX1Y1, mpi.mpiTag,  gasX1Y1);
      reqX1Y1[1] = mpi.boostWorld.irecv(mpi.rankX1Y1, mpi.mpiTag, rGasX1Y1);
    }
    if (mpi.rankX1Y2 >= 0) { // edge x1y2
      IJK start(boundCup.lowX, boundCup.uppY+1-haloGrid, boundCup.lowZ);
      IJK end(boundCup.lowX+haloGrid-1, boundCup.uppY, boundCup.uppZ);
      findPenalInRectangle(start, end, gasX1Y2);
      reqX1Y2[0] = mpi.boostWorld.isend(mpi.rankX1Y2, mpi.mpiTag,  gasX1Y2);
      reqX1Y2[1] = mpi.boostWorld.irecv(mpi.rankX1Y2, mpi.mpiTag, rGasX1Y2);
    }
    if (mpi.rankX1Z1 >= 0) { // edge x1z1
      IJK start(boundCup.lowX, boundCup.lowY, boundCup.lowZ);
      IJK end(boundCup.lowX+haloGrid-1, boundCup.uppY, boundCup.lowZ+haloGrid-1);
      findPenalInRectangle(start, end, gasX1Z1);
      reqX1Z1[0] = mpi.boostWorld.isend(mpi.rankX1Z1, mpi.mpiTag,  gasX1Z1);
      reqX1Z1[1] = mpi.boostWorld.irecv(mpi.rankX1Z1, mpi.mpiTag, rGasX1Z1);
    }
    if (mpi.rankX1Z2 >= 0) { // edge x1z2
      IJK start(boundCup.lowX, boundCup.lowY, boundCup.uppZ+1-haloGrid);
      IJK end(boundCup.lowX+haloGrid-1, boundCup.uppY, boundCup.uppZ);
      findPenalInRectangle(start, end, gasX1Z2);
      reqX1Z2[0] = mpi.boostWorld.isend(mpi.rankX1Z2, mpi.mpiTag,  gasX1Z2);
      reqX1Z2[1] = mpi.boostWorld.irecv(mpi.rankX1Z2, mpi.mpiTag, rGasX1Z2);
    }
    if (mpi.rankX2Y1 >= 0) { // edge x2y1
      IJK start(boundCup.uppX+1-haloGrid, boundCup.lowY, boundCup.lowZ);
      IJK end(boundCup.uppX, boundCup.lowY+haloGrid-1, boundCup.uppZ);
      findPenalInRectangle(start, end, gasX2Y1);
      reqX2Y1[0] = mpi.boostWorld.isend(mpi.rankX2Y1, mpi.mpiTag,  gasX2Y1);
      reqX2Y1[1] = mpi.boostWorld.irecv(mpi.rankX2Y1, mpi.mpiTag, rGasX2Y1);
    }
    if (mpi.rankX2Y2 >= 0) { // edge x2y2
      IJK start(boundCup.uppX+1-haloGrid, boundCup.uppX+1-haloGrid, boundCup.lowZ);
      IJK end(boundCup.uppX, boundCup.uppY, boundCup.uppZ);
      findPenalInRectangle(start, end, gasX2Y2);
      reqX2Y2[0] = mpi.boostWorld.isend(mpi.rankX2Y2, mpi.mpiTag,  gasX2Y2);
      reqX2Y2[1] = mpi.boostWorld.irecv(mpi.rankX2Y2, mpi.mpiTag, rGasX2Y2);
    }
    if (mpi.rankX2Z1 >= 0) { // edge x2z1
      IJK start(boundCup.uppX+1-haloGrid, boundCup.lowY, boundCup.lowZ);
      IJK end(boundCup.uppX, boundCup.uppY, boundCup.lowZ+haloGrid-1);
      findPenalInRectangle(start, end, gasX2Z1);
      reqX2Z1[0] = mpi.boostWorld.isend(mpi.rankX2Z1, mpi.mpiTag,  gasX2Z1);
      reqX2Z1[1] = mpi.boostWorld.irecv(mpi.rankX2Z1, mpi.mpiTag, rGasX2Z1);
    }
    if (mpi.rankX2Z2 >= 0) { // edge x2z2
      IJK start(boundCup.uppX+1-haloGrid, boundCup.lowY, boundCup.uppZ+1-haloGrid);
      IJK end(boundCup.uppX, boundCup.uppY, boundCup.uppZ);
      findPenalInRectangle(start, end, gasX2Z2);
      reqX2Z2[0] = mpi.boostWorld.isend(mpi.rankX2Z2, mpi.mpiTag,  gasX2Z2);
      reqX2Z2[1] = mpi.boostWorld.irecv(mpi.rankX2Z2, mpi.mpiTag, rGasX2Z2);
    }
    if (mpi.rankY1Z1 >= 0) { // edge y1z1
      IJK start(boundCup.lowX, boundCup.lowY, boundCup.lowZ);
      IJK end(boundCup.uppX, boundCup.lowY+haloGrid-1,boundCup.lowZ+ haloGrid-1);
      findPenalInRectangle(start, end, gasY1Z1);
      reqY1Z1[0] = mpi.boostWorld.isend(mpi.rankY1Z1, mpi.mpiTag,  gasY1Z1);
      reqY1Z1[1] = mpi.boostWorld.irecv(mpi.rankY1Z1, mpi.mpiTag, rGasY1Z1);
    }
    if (mpi.rankY1Z2 >= 0) { // edge y1z2
      IJK start(boundCup.lowX, boundCup.lowY, boundCup.uppZ+1-haloGrid);
      IJK end(boundCup.uppX, boundCup.lowY+haloGrid-1, boundCup.uppZ);
      findPenalInRectangle(start, end, gasY1Z2);
      reqY1Z2[0] = mpi.boostWorld.isend(mpi.rankY1Z2, mpi.mpiTag,  gasY1Z2);
      reqY1Z2[1] = mpi.boostWorld.irecv(mpi.rankY1Z2, mpi.mpiTag, rGasY1Z2);
    }
    if (mpi.rankY2Z1 >= 0) { // edge y2z1
      IJK start(boundCup.lowX, boundCup.uppY+1-haloGrid, boundCup.lowZ);
      IJK end(boundCup.uppX, boundCup.uppY, boundCup.lowZ+haloGrid-1);
      findPenalInRectangle(start, end, gasY2Z1);
      reqY2Z1[0] = mpi.boostWorld.isend(mpi.rankY2Z1, mpi.mpiTag,  gasY2Z1);
      reqY2Z1[1] = mpi.boostWorld.irecv(mpi.rankY2Z1, mpi.mpiTag, rGasY2Z1);
    }
    if (mpi.rankY2Z2 >= 0) { // edge y2z2
      IJK start(boundCup.lowX, boundCup.uppY+1-haloGrid, boundCup.uppZ+1-haloGrid);
      IJK end(boundCup.uppX, boundCup.uppY, boundCup.uppZ);
      findPenalInRectangle(start, end, gasY2Z2);
      reqY2Z2[0] = mpi.boostWorld.isend(mpi.rankY2Z2, mpi.mpiTag,  gasY2Z2);
      reqY2Z2[1] = mpi.boostWorld.irecv(mpi.rankY2Z2, mpi.mpiTag, rGasY2Z2);
    }
    // 8 vertices
    if (mpi.rankX1Y1Z1 >= 0) { // vertice x1y1z1
      IJK start(boundCup.lowX, boundCup.lowY, boundCup.lowZ);
      IJK end(boundCup.lowX+haloGrid-1, boundCup.lowY+haloGrid-1, boundCup.lowZ+haloGrid-1);
      findPenalInRectangle(start, end, gasX1Y1Z1);
      reqX1Y1Z1[0] = mpi.boostWorld.isend(mpi.rankX1Y1Z1, mpi.mpiTag,  gasX1Y1Z1);
      reqX1Y1Z1[1] = mpi.boostWorld.irecv(mpi.rankX1Y1Z1, mpi.mpiTag, rGasX1Y1Z1);
    }
    if (mpi.rankX1Y1Z2 >= 0) { // vertice x1y1z2
      IJK start(boundCup.lowX, boundCup.lowY,  boundCup.uppZ+1-haloGrid);
      IJK end(boundCup.lowX+haloGrid-1, boundCup.lowY+haloGrid-1, boundCup.uppZ);
      findPenalInRectangle(start, end, gasX1Y1Z2);
      reqX1Y1Z2[0] = mpi.boostWorld.isend(mpi.rankX1Y1Z2, mpi.mpiTag,  gasX1Y1Z2);
      reqX1Y1Z2[1] = mpi.boostWorld.irecv(mpi.rankX1Y1Z2, mpi.mpiTag, rGasX1Y1Z2);
    }
    if (mpi.rankX1Y2Z1 >= 0) { // vertice x1y2z1
      IJK start(boundCup.lowX, boundCup.uppY+1-haloGrid, boundCup.lowZ);
      IJK end(boundCup.lowX+haloGrid-1, boundCup.uppY, boundCup.lowZ+haloGrid-1);
      findPenalInRectangle(start, end, gasX1Y2Z1);
      reqX1Y2Z1[0] = mpi.boostWorld.isend(mpi.rankX1Y2Z1, mpi.mpiTag,  gasX1Y2Z1);
      reqX1Y2Z1[1] = mpi.boostWorld.irecv(mpi.rankX1Y2Z1, mpi.mpiTag, rGasX1Y2Z1);
    }
    if (mpi.rankX1Y2Z2 >= 0) { // vertice x1y2z2
      IJK start(boundCup.lowX, boundCup.uppY+1-haloGrid, boundCup.uppZ+1-haloGrid);
      IJK end(boundCup.lowX+haloGrid-1, boundCup.uppY, boundCup.uppZ);
      findPenalInRectangle(start, end, gasX1Y2Z2);
      reqX1Y2Z2[0] = mpi.boostWorld.isend(mpi.rankX1Y2Z2, mpi.mpiTag,  gasX1Y2Z2);
      reqX1Y2Z2[1] = mpi.boostWorld.irecv(mpi.rankX1Y2Z2, mpi.mpiTag, rGasX1Y2Z2);
    }
    if (mpi.rankX2Y1Z1 >= 0) { // vertice x2y1z1
      IJK start(boundCup.uppX+1-haloGrid, boundCup.lowY, boundCup.lowZ);
      IJK end(boundCup.uppX, boundCup.lowY+haloGrid-1, boundCup.lowZ+haloGrid-1);
      findPenalInRectangle(start, end, gasX2Y1Z1);
      reqX2Y1Z1[0] = mpi.boostWorld.isend(mpi.rankX2Y1Z1, mpi.mpiTag,  gasX2Y1Z1);
      reqX2Y1Z1[1] = mpi.boostWorld.irecv(mpi.rankX2Y1Z1, mpi.mpiTag, rGasX2Y1Z1);
    }
    if (mpi.rankX2Y1Z2 >= 0) { // vertice x2y1z2
      IJK start(boundCup.uppX+1-haloGrid, boundCup.lowY, boundCup.uppZ+1-haloGrid);
      IJK end(boundCup.uppX, boundCup.lowY+haloGrid-1, boundCup.uppZ);
      findPenalInRectangle(start, end, gasX2Y1Z2);
      reqX2Y1Z2[0] = mpi.boostWorld.isend(mpi.rankX2Y1Z2, mpi.mpiTag,  gasX2Y1Z2);
      reqX2Y1Z2[1] = mpi.boostWorld.irecv(mpi.rankX2Y1Z2, mpi.mpiTag, rGasX2Y1Z2);
    }
    if (mpi.rankX2Y2Z1 >= 0) { // vertice x2y2z1
      IJK start(boundCup.uppX+1-haloGrid, boundCup.uppY+1-haloGrid, boundCup.lowZ);
      IJK end(boundCup.uppX, boundCup.uppY, boundCup.lowZ+haloGrid-1);
      findPenalInRectangle(start, end, gasX2Y2Z1);
      reqX2Y2Z1[0] = mpi.boostWorld.isend(mpi.rankX2Y2Z1, mpi.mpiTag,  gasX2Y2Z1);
      reqX2Y2Z1[1] = mpi.boostWorld.irecv(mpi.rankX2Y2Z1, mpi.mpiTag, rGasX2Y2Z1);
    }
    if (mpi.rankX2Y2Z2 >= 0) { // vertice x2y2z2
      IJK start(boundCup.uppX+1-haloGrid, boundCup.uppY+1-haloGrid, boundCup.uppZ+1-haloGrid);
      IJK end(boundCup.uppX, boundCup.uppY, boundCup.uppZ);
      findPenalInRectangle(start, end, gasX2Y2Z2);
      reqX2Y2Z2[0] = mpi.boostWorld.isend(mpi.rankX2Y2Z2, mpi.mpiTag,  gasX2Y2Z2);
      reqX2Y2Z2[1] = mpi.boostWorld.irecv(mpi.rankX2Y2Z2, mpi.mpiTag, rGasX2Y2Z2);
    }

    // 6 surfaces
    if (mpi.rankX1 >= 0) boost::mpi::wait_all(reqX1, reqX1 + 2);
    if (mpi.rankX2 >= 0) boost::mpi::wait_all(reqX2, reqX2 + 2);
    if (mpi.rankY1 >= 0) boost::mpi::wait_all(reqY1, reqY1 + 2);
    if (mpi.rankY2 >= 0) boost::mpi::wait_all(reqY2, reqY2 + 2);
    if (mpi.rankZ1 >= 0) boost::mpi::wait_all(reqZ1, reqZ1 + 2);
    if (mpi.rankZ2 >= 0) boost::mpi::wait_all(reqZ2, reqZ2 + 2);
    // 12 edges
    if (mpi.rankX1Y1 >= 0) boost::mpi::wait_all(reqX1Y1, reqX1Y1 + 2);
    if (mpi.rankX1Y2 >= 0) boost::mpi::wait_all(reqX1Y2, reqX1Y2 + 2);  
    if (mpi.rankX1Z1 >= 0) boost::mpi::wait_all(reqX1Z1, reqX1Z1 + 2);
    if (mpi.rankX1Z2 >= 0) boost::mpi::wait_all(reqX1Z2, reqX1Z2 + 2);
    if (mpi.rankX2Y1 >= 0) boost::mpi::wait_all(reqX2Y1, reqX2Y1 + 2);
    if (mpi.rankX2Y2 >= 0) boost::mpi::wait_all(reqX2Y2, reqX2Y2 + 2);  
    if (mpi.rankX2Z1 >= 0) boost::mpi::wait_all(reqX2Z1, reqX2Z1 + 2);
    if (mpi.rankX2Z2 >= 0) boost::mpi::wait_all(reqX2Z2, reqX2Z2 + 2); 
    if (mpi.rankY1Z1 >= 0) boost::mpi::wait_all(reqY1Z1, reqY1Z1 + 2);
    if (mpi.rankY1Z2 >= 0) boost::mpi::wait_all(reqY1Z2, reqY1Z2 + 2);
    if (mpi.rankY2Z1 >= 0) boost::mpi::wait_all(reqY2Z1, reqY2Z1 + 2);
    if (mpi.rankY2Z2 >= 0) boost::mpi::wait_all(reqY2Z2, reqY2Z2 + 2); 
    // 8 vertices
    if (mpi.rankX1Y1Z1 >= 0) boost::mpi::wait_all(reqX1Y1Z1, reqX1Y1Z1 + 2);
    if (mpi.rankX1Y1Z2 >= 0) boost::mpi::wait_all(reqX1Y1Z2, reqX1Y1Z2 + 2);
    if (mpi.rankX1Y2Z1 >= 0) boost::mpi::wait_all(reqX1Y2Z1, reqX1Y2Z1 + 2);
    if (mpi.rankX1Y2Z2 >= 0) boost::mpi::wait_all(reqX1Y2Z2, reqX1Y2Z2 + 2);
    if (mpi.rankX2Y1Z1 >= 0) boost::mpi::wait_all(reqX2Y1Z1, reqX2Y1Z1 + 2);
    if (mpi.rankX2Y1Z2 >= 0) boost::mpi::wait_all(reqX2Y1Z2, reqX2Y1Z2 + 2);
    if (mpi.rankX2Y2Z1 >= 0) boost::mpi::wait_all(reqX2Y2Z1, reqX2Y2Z1 + 2);
    if (mpi.rankX2Y2Z2 >= 0) boost::mpi::wait_all(reqX2Y2Z2, reqX2Y2Z2 + 2);

    // merge using std::map
    std::vector<ValCoord> recvPenalVec;
    // 6 surfaces                                                                                          
    if (mpi.rankX1 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX1.begin(), rGasX1.end());            
    if (mpi.rankX2 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX2.begin(), rGasX2.end());            
    if (mpi.rankY1 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasY1.begin(), rGasY1.end());            
    if (mpi.rankY2 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasY2.begin(), rGasY2.end());            
    if (mpi.rankZ1 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasZ1.begin(), rGasZ1.end());            
    if (mpi.rankZ2 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasZ2.begin(), rGasZ2.end());            
    // 12 edges                                                                                            
    if (mpi.rankX1Y1 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX1Y1.begin(), rGasX1Y1.end());      
    if (mpi.rankX1Y2 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX1Y2.begin(), rGasX1Y2.end());      
    if (mpi.rankX1Z1 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX1Z1.begin(), rGasX1Z1.end());      
    if (mpi.rankX1Z2 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX1Z2.begin(), rGasX1Z2.end());      
    if (mpi.rankX2Y1 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX2Y1.begin(), rGasX2Y1.end());      
    if (mpi.rankX2Y2 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX2Y2.begin(), rGasX2Y2.end());      
    if (mpi.rankX2Z1 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX2Z1.begin(), rGasX2Z1.end());      
    if (mpi.rankX2Z2 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX2Z2.begin(), rGasX2Z2.end());      
    if (mpi.rankY1Z1 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasY1Z1.begin(), rGasY1Z1.end());      
    if (mpi.rankY1Z2 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasY1Z2.begin(), rGasY1Z2.end());      
    if (mpi.rankY2Z1 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasY2Z1.begin(), rGasY2Z1.end());      
    if (mpi.rankY2Z2 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasY2Z2.begin(), rGasY2Z2.end());      
    // 8 vertices                                                                                          
    if (mpi.rankX1Y1Z1 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX1Y1Z1.begin(), rGasX1Y1Z1.end());
    if (mpi.rankX1Y1Z2 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX1Y1Z2.begin(), rGasX1Y1Z2.end());
    if (mpi.rankX1Y2Z1 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX1Y2Z1.begin(), rGasX1Y2Z1.end());
    if (mpi.rankX1Y2Z2 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX1Y2Z2.begin(), rGasX1Y2Z2.end());
    if (mpi.rankX2Y1Z1 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX2Y1Z1.begin(), rGasX2Y1Z1.end());
    if (mpi.rankX2Y1Z2 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX2Y1Z2.begin(), rGasX2Y1Z2.end());
    if (mpi.rankX2Y2Z1 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX2Y2Z1.begin(), rGasX2Y2Z1.end());
    if (mpi.rankX2Y2Z2 >= 0) recvPenalVec.insert(recvPenalVec.end(), rGasX2Y2Z2.begin(), rGasX2Y2Z2.end());

    IJK_arrayPenalForce_Map recvGasMap;
    for (int i = 0; i < recvPenalVec.size(); ++i) {
      IJK tri;
      coordToGlobalIndex(recvPenalVec[i].coords, tri);
      recvGasMap[tri] = recvPenalVec[i];
    }
    
    // construct penalized gas with received info
    for (std::size_t i = boundPrn.lowX; i <= boundPrn.uppX; ++i)
      for (std::size_t j = boundPrn.lowY; j <= boundPrn.uppY; ++j)
	for (std::size_t k = boundPrn.lowZ; k <= boundPrn.uppZ; ++k) {
	  // i, j, k are local
	  IJK local(i,j,k);
	  IJK global;
	  localIndexToGlobal(local, global);

	  if (recvGasMap.count(global)) { // important to judge, otherwise it creates a new map elment!
	    arrayPenalForce[i][j][k][0] = recvGasMap[global].coords.getX();
	    arrayPenalForce[i][j][k][1] = recvGasMap[global].coords.getY();
	    arrayPenalForce[i][j][k][2] = recvGasMap[global].coords.getZ();
	  }
	}

  } // end of Gas::backCommu26

} // name space dem
