#include "MPIFrame.h"
#include "const.h"
#include "Parameter.h"

namespace dem {

  void MPIFrame::setCommunicator(boost::mpi::communicator &comm) {
    boostWorld = comm;
    mpiWorld = MPI_Comm(comm);
    mpiProcX = static_cast<int> (dem::Parameter::getSingleton().parameter["mpiProcX"]);
    mpiProcY = static_cast<int> (dem::Parameter::getSingleton().parameter["mpiProcY"]);
    mpiProcZ = static_cast<int> (dem::Parameter::getSingleton().parameter["mpiProcZ"]);
  
    // create Cartesian virtual topology (unavailable in boost.mpi) 
    int ndim = 3;
    int dims[3] = {mpiProcX, mpiProcY, mpiProcZ};
    int periods[3] = {0, 0, 0};
    int reorder = 0; // mpiRank not reordered
    MPI_Cart_create(mpiWorld, ndim, dims, periods, reorder, &cartComm);
    MPI_Comm_rank(cartComm, &mpiRank); 
    MPI_Comm_size(cartComm, &mpiSize);
    MPI_Cart_coords(cartComm, mpiRank, ndim, mpiCoords);
    mpiTag = 0;
    assert(mpiRank == boostWorld.rank());
    //std::cout << "mpiRank=" << mpiRank << " " << mpiCoords[0] << " " << mpiCoords[1] << " " << mpiCoords[2] << std::endl;
  }

  void MPIFrame::setBdryProcess() {
    for (int iRank = 0; iRank < mpiSize; ++iRank) {
      int ndim = 3;
      int coords[3];
      MPI_Cart_coords(cartComm, iRank, ndim, coords);
      if (coords[0] == 0 || coords[0] == mpiProcX - 1 ||
	  coords[1] == 0 || coords[1] == mpiProcY - 1 ||
	  coords[2] == 0 || coords[2] == mpiProcZ - 1)
	bdryProcess.push_back(iRank);
    }
  }

  void MPIFrame::findNeighborProcess() {
    // find neighboring process
    rankX1 = -1; rankX2 = -1; rankY1 = -1; rankY2 = -1; rankZ1 = -1; rankZ2 = -1;
    rankX1Y1 = -1; rankX1Y2 = -1; rankX1Z1 = -1; rankX1Z2 = -1; 
    rankX2Y1 = -1; rankX2Y2 = -1; rankX2Z1 = -1; rankX2Z2 = -1; 
    rankY1Z1 = -1; rankY1Z2 = -1; rankY2Z1 = -1; rankY2Z2 = -1; 
    rankX1Y1Z1 = -1; rankX1Y1Z2 = -1; rankX1Y2Z1 = -1; rankX1Y2Z2 = -1; 
    rankX2Y1Z1 = -1; rankX2Y1Z2 = -1; rankX2Y2Z1 = -1; rankX2Y2Z2 = -1;
    // x1: -x direction
    int neighborCoords[3] = {mpiCoords[0], mpiCoords[1], mpiCoords[2]};
    --neighborCoords[0];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1);
    // x2: +x direction
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2);
    // y1: -y direction
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[1];
    MPI_Cart_rank(cartComm, neighborCoords, &rankY1);
    // y2: +y direction
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[1];
    MPI_Cart_rank(cartComm, neighborCoords, &rankY2);
    // z1: -z direction
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankZ1);
    // z2: +z direction
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankZ2);
    // x1y1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[0]; --neighborCoords[1];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y1);
    // x1y2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[0]; ++neighborCoords[1];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y2);
    // x1z1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[0]; --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1Z1);
    // x1z2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[0]; ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1Z2);
    // x2y1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0]; --neighborCoords[1];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y1);
    // x2y2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0]; ++neighborCoords[1];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y2);
    // x2z1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0]; --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2Z1);
    // x2z2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0]; ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2Z2);
    // y1z1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[1]; --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankY1Z1);
    // y1z2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[1]; ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankY1Z2);
    // y2z1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[1]; --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankY2Z1);
    // y2z2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[1]; ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankY2Z2);
    // x1y1z1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[0]; --neighborCoords[1]; --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y1Z1);
    // x1y1z2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[0]; --neighborCoords[1]; ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y1Z2);
    // x1y2z1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[0]; ++neighborCoords[1]; --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y2Z1);
    // x1y2z2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[0]; ++neighborCoords[1]; ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y2Z2);
    // x2y1z1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0]; --neighborCoords[1]; --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y1Z1);
    // x2y1z2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0]; --neighborCoords[1]; ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y1Z2);
    // x2y2z1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0]; ++neighborCoords[1]; --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y2Z1);
    // x2y2z2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0]; ++neighborCoords[1]; ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y2Z2);
    //std::cout << std::endl << "findMPINeighbor:rankX1=" << rankX1 << std::endl;
  }

  bool MPIFrame::isBdryProcess() {
    return (mpiCoords[0] == 0 || mpiCoords[0] == mpiProcX - 1 ||
	    mpiCoords[1] == 0 || mpiCoords[1] == mpiProcY - 1 ||
	    mpiCoords[2] == 0 || mpiCoords[2] == mpiProcZ - 1);
  }

  bool MPIFrame::isBdryProcessMin() {
    return (mpiCoords[0] == 0 ||
	    mpiCoords[1] == 0 ||
	    mpiCoords[2] == 0);
  }

  bool MPIFrame::isBdryProcessMax() {
    return (mpiCoords[0] == mpiProcX - 1 ||
	    mpiCoords[1] == mpiProcY - 1 ||
	    mpiCoords[2] == mpiProcZ - 1);
  }

  bool MPIFrame::isBdryProcessXMin() {
    return (mpiCoords[0] == 0);
  }

  bool MPIFrame::isBdryProcessYMin() {
    return (mpiCoords[1] == 0);
  }

  bool MPIFrame::isBdryProcessZMin() {
    return (mpiCoords[2] == 0);
  }

  bool MPIFrame::isBdryProcessXMax() {
    return (mpiCoords[0] == mpiProcX - 1);
  }

  bool MPIFrame::isBdryProcessYMax() {
    return (mpiCoords[1] == mpiProcY - 1);
  }

  bool MPIFrame::isBdryProcessZMax() {
    return (mpiCoords[2] == mpiProcZ - 1);
  }

} // namespace dem
