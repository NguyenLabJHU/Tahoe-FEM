///////////////////////////////////////////////////////////////////////////////////////////////////////
//                                   Code: ParaEllip3d-CFD                                           //
//                                 Author: Dr. Beichuan Yan                                          //
//                                  Email: beichuan.yan@colorado.edu                                 //
//                              Institute: University of Colorado Boulder                            //
///////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MPIFRAME_H
#define MPIFRAME_H

#include <vector>
#include <boost/mpi.hpp>

namespace dem {
  
  class MPIFrame {
    
  public:
    boost::mpi::communicator boostWorld;
    MPI_Comm mpiWorld, cartComm;
    int mpiProcX, mpiProcY, mpiProcZ;
    int mpiRank, mpiSize, mpiTag, mpiCoords[3];
    std::vector<std::size_t> bdryProcess;

    int rankX1, rankX2, rankY1, rankY2, rankZ1, rankZ2;
    int rankX1Y1, rankX1Y2, rankX1Z1, rankX1Z2; 
    int rankX2Y1, rankX2Y2, rankX2Z1, rankX2Z2; 
    int rankY1Z1, rankY1Z2, rankY2Z1, rankY2Z2; 
    int rankX1Y1Z1, rankX1Y1Z2, rankX1Y2Z1, rankX1Y2Z2; 
    int rankX2Y1Z1, rankX2Y1Z2, rankX2Y2Z1, rankX2Y2Z2;

  public:
    void setCommunicator(boost::mpi::communicator &comm);
    void findNeighborProcess();
    void setBdryProcess();
    bool isBdryProcess();
    bool isBdryProcessMin();
    bool isBdryProcessMax();
    bool isBdryProcessXMin();
    bool isBdryProcessYMin();
    bool isBdryProcessZMin();
    bool isBdryProcessXMax();
    bool isBdryProcessYMax();
    bool isBdryProcessZMax();
  };

} // namespace dem

#endif
