//                          ---------------------
//                         /                    /|
//                        /                    / |
//                       /                    /  |
//                      /         z2 (5)     /   |
//                     /                    /    |height
//                    /                    /     |                    z (sigma3)
//                   /                    /      |                    |
//                  |---------------------       |                    |
//                  |                    | y2(2) |                    |____ y (sigma1)
//                  |                    |       /                   /
//                  |                    |      /                   /
//                  |        x2 (1)      |     /                   x (sigma2) 
//                  |                    |    /length
//                  |                    |   /
//                  |                    |  /
//                  |                    | /
//                  |                    |/
//                  ----------------------
//                         width
//
//    It is preferable to use the description of surface x1, x2, y1, y2, z1, z2,
//    where x1 < x2, y1 < y2, z1 < z2
//
//    sigma1_1 & sigma1_2 refers to side 2 & side 4 respectively,
//    sigma2_1 & sigma2_2 refers to side 1 & side 3 respectively,
//    sigma3_1 & sigma3_2 refers to side 5 & side 6 respectively,
//
//    int mid[2]={1,3};    // boundary 1 and 3
//    int max[2]={2,4};    // boundary 2 and 4
//    int min[2]={5,6};    // boundary 5 and 6
//    min/mid/max does not mean actual magnitude of values, just signs

#include "Assembly.h"
#include "parameter.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <ctime>
#include <cassert>
#include <utility>
#include <sys/time.h>

#ifdef OPENMP
#include <omp.h>
#define OPENMP_IMPL 2
#endif

// OPENMP_IMPL: 
// 0: OpenMP implementation 0, ts partitions, based on linked list
// 1: OpenMP implementation 1, ts partitions, based on vector
// 2: OpenMP implementation 2, no partition, each thread leaps by ts until completed
// 3: OpenMP implementation 3, no partition, each thread leaps by ts until num/2 and handles two particles.
// 4: OpenMP implementation 4, no partition, parallel for, various loop scheduling: (static), (static,1), (dynamic), (dynamic,1)

//#define BINNING
#define TIME_PROFILE

using std::cout;
using std::setw;
using std::endl;
using std::flush;
using std::vector;
using std::pair;

static time_t timeStamp; // for file timestamping
static struct timeval time_w1, time_w2; // for wall-clock time record
static struct timeval time_p1, time_p2; // for internal wall-clock time profiling, can be used on any piece of code
static struct timeval time_r1, time_r2; // for internal wall-clock time profiling for contact resolution only (excluding space search)

namespace dem {

std::ofstream progressinf;

struct timeval timediff(const struct timeval &time1, const struct timeval &time2) {
  struct timeval diff;
  diff.tv_sec = time2.tv_sec - time1.tv_sec; 
  if( ( diff.tv_usec = time2.tv_usec - time1.tv_usec) < 0) {
    diff.tv_usec += 1000000;
    --diff.tv_sec;
  }
  return(diff); 
}


long int timediffmsec(const struct timeval &time1, const struct timeval &time2) {
  struct timeval diff = timediff(time1, time2);
  return(diff.tv_sec * 1000000 + diff.tv_usec); 
}


REAL timediffsec(const struct timeval &time1, const struct timeval &time2) {
  return( (REAL) timediffmsec(time1,time2) / 1.0e+6);
}


const char* combineString(const char *str, int num) {
  std::string obj(str);
  std::stringstream ss;
  ss << std::setw(3) << std::setfill('0') << std::right << num;
  obj += ss.str();
  return obj.c_str();
}


void Assembly::setCommunicator(boost::mpi::communicator &comm) {
  boostWorld = comm;
  mpiWorld = MPI_Comm(comm);
  
  // create Cartesian virtual topology (unavailable in boost.mpi) 
  int ndim = 3;
  int dims[3] = {NPX, NPY, NPZ};
  int periods[3] = {0, 0, 0};
  int reorder = 1;
  MPI_Cart_create(mpiWorld, ndim, dims, periods, reorder, &cartComm);
  MPI_Comm_rank(cartComm, &mpiRank); // mpiRank probably reordered
  MPI_Comm_size(cartComm, &mpiSize);
  MPI_Cart_coords(cartComm, mpiRank, ndim, mpiCoords);
  mpiTag = 0;
  assert(mpiRank == boostWorld.rank());
}


void Assembly::depositIntoContainer() 
{
  int simuType, freeType, totalSteps, snapNum, statInterval;
  REAL trimHeight;

  if (mpiRank == 0) {
    readParameter(simuType, freeType, trimHeight, totalSteps, snapNum, statInterval);
    buildBoundary(5, "dep_boundary");
    generateParticle(freeType, "float_particle"); 
  }
    
  broadcast(boostWorld, totalSteps, 0);
  broadcast(boostWorld, snapNum, 0);
  broadcast(boostWorld, statInterval, 0);
  deposit(totalSteps,
	  snapNum,
	  statInterval,
	  "float_particle",
	  "dep_boundary");

  if (mpiRank == 0) {
    setContainer(Rectangle(allContainer.getMinCorner().getx(),
			   allContainer.getMinCorner().gety(),
			   allContainer.getMinCorner().getz(),
			   allContainer.getMaxCorner().getx(),
			   allContainer.getMaxCorner().gety(),
			   trimHeight));
    buildBoundary("trm_boundary");
    trim(false,
	 "dep_particle_end",
	 "trm_particle");
  }
  
}


void Assembly::
readParameter(int  &simuType,
	      int  &freeType,
	      REAL &trimHeight,
	      int  &totalSteps,
	      int  &snapNum,
	      int  &statInterval) 
{
  std::ifstream ifs;
  ifs.open("parameter");
  if (!ifs) { cout << "stream error!" << endl; exit(-1); }

  std::string str;
  ifs >> str;
  ifs >> simuType;

  ifs >> str;
  REAL x1, y1, z1, x2, y2, z2;
  ifs >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;
  setContainer(Rectangle(x1, y1, z1, x2, y2, z2));

  ifs >> str;
  int sieveNum;
  ifs >> sieveNum;
  std::vector<REAL> percent(sieveNum), size(sieveNum);
  for (int i = 0; i < sieveNum; ++i)
    ifs >> percent[i] >> size[i];

  ifs >> str;
  REAL ratio_ba, ratio_ca;
  ifs >> ratio_ba >> ratio_ca;
  setGradation(Gradation(sieveNum, percent, size, ratio_ba, ratio_ca));

  ifs >> str;
  ifs >> freeType;

  ifs >> str;
  ifs >> trimHeight;

  ifs >> str;
  ifs >> totalSteps;

  ifs >> str;
  ifs >> snapNum;

  ifs >> str;
  ifs >> statInterval;

  ifs >> str;
  ifs >> TIMESTEP;

  ifs >> str;
  ifs >> MASS_SCL;

  ifs >> str;
  ifs >> MNT_SCL;

  ifs >> str;
  ifs >> GRVT_SCL;

  ifs >> str;
  ifs >> DMP_F;

  ifs >> str;
  ifs >> DMP_M;

  ifs.close();

}

// freeType:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
void Assembly::
generateParticle(int freeType,
		 const char *genParticle)
{
  REAL x,y,z;
  Particle* newptcl;
  int particleNum = 0;
  REAL dimx     = allContainer.getDimx();
  REAL dimy     = allContainer.getDimy();
  REAL dimz     = allContainer.getDimz();
  REAL diameter = gradation.getPtclMaxRadius()*2.0;

  REAL offset   = 0;
  REAL edge     = diameter;
  if (gradation.getSize().size() == 1 &&
      gradation.getPtclRatioBA() == 1.0 && 
      gradation.getPtclRatioCA() == 1.0) {
    edge   = diameter*2.0;
    offset = diameter*0.25;
  }
  
  REAL x1 = allContainer.getMinCorner().getx() + edge;
  REAL y1 = allContainer.getMinCorner().gety() + edge;
  REAL z1 = allContainer.getMinCorner().getz() + diameter;
  REAL x2 = allContainer.getMaxCorner().getx() - edge;
  REAL y2 = allContainer.getMaxCorner().gety() - edge;
  REAL z2 = allContainer.getMaxCorner().getz() - diameter;
  REAL x0 = allContainer.getCenter().getx();
  REAL y0 = allContainer.getCenter().gety();
  REAL z0 = allContainer.getCenter().getz();

  if (freeType == 0) {      // just one free particle
    newptcl = new Particle(particleNum+1, 0, Vec(x0,y0,z0), gradation, YOUNG, POISSON);
    allParticleVec.push_back(newptcl);
    particleNum++;
  }
  else if (freeType == 1) { // a horizontal layer of free particles
    for (x = x1; x - x2 < EPS; x += diameter)
      for (y = y1; y - y2 < EPS; y += diameter) {
	newptcl = new Particle(particleNum+1, 0, Vec(x,y,z0), gradation, YOUNG, POISSON);
	allParticleVec.push_back(newptcl);
	particleNum++;
      }
  }
  else if (freeType == 2) { // multiple layers of free particles
    for (z = z1; z - z2 < EPS; z += diameter) {
      for (x = x1 + offset; x - x2 < EPS; x += diameter)
	for (y = y1 + offset; y - y2 < EPS; y += diameter) {
	  newptcl = new Particle(particleNum+1, 0, Vec(x,y,z), gradation, YOUNG, POISSON);
	  allParticleVec.push_back(newptcl);
	  particleNum++;
	}	
      offset *= -1;
    }
  }
  
  printParticle(genParticle); 
}


void Assembly::
trim(bool toRebuild,
     const char* inputParticle,
     const char* trmParticle)
{
  if (toRebuild) readParticle(inputParticle);
  trimHistoryNum = allParticleVec.size();

  Vec  v1 = allContainer.getMinCorner();
  Vec  v2 = allContainer.getMaxCorner();
  Vec  v0 = allContainer.getCenter();
  REAL x1 = v1.getx();
  REAL y1 = v1.gety();
  REAL z1 = v1.getz();
  REAL x2 = v2.getx();
  REAL y2 = v2.gety();
  REAL z2 = v2.getz();
  REAL x0 = v0.getx();
  REAL y0 = v0.gety();
  REAL z0 = v0.getz();
 
  std::vector<Particle*>::iterator itr;
  Vec center;

  for (itr = allParticleVec.begin(); itr != allParticleVec.end(); ){
    center=(*itr)->getCurrPos();
    if(center.getx() <= x1 || center.getx() >= x2 ||
       center.gety() <= y1 || center.gety() >= y2 ||
       center.getz() <= z1 || center.getz() + gradation.getPtclMaxRadius() > z2)
      {
	delete (*itr); // release memory
	itr = allParticleVec.erase(itr); 
      }
    else
      ++itr;
  }
  
  printParticle("trm_particle");
}


void Assembly::
deposit(int totalSteps,  
	int snapNum,
	int statInterval,
	const char *inputParticle,
	const char *inputBoundary) 
{
  if (mpiRank == 0) {
    readBoundary(inputBoundary);
    readParticle(inputParticle); 
  }

  iteration = 0;
  int iterSnap = 0;
  do {
    
    if (mpiRank == 0) {
      if (iteration % (totalSteps/snapNum) == 0) {
	printParticle(combineString("dep_particle_", iterSnap));
	++iterSnap;
      }
      clearForce();
      findParticleOnBoundary();
      boundaryForce();
      //cout << "\n" << flush;
    }

    partiCommuParticle();
    findContact();
    internalForce();
    updateParticle();

    gatherParticle();
    
  } while (++iteration < totalSteps);

  printParticle("dep_particle_end");

}


void Assembly::
findParticleInRectangle(Rectangle &container,
			std::vector<Particle*> &inputParticle,
			std::vector<Particle*> &foundParticle) {
  Vec  v1 = container.getMinCorner();
  Vec  v2 = container.getMaxCorner();
  REAL x1 = v1.getx();
  REAL y1 = v1.gety();
  REAL z1 = v1.getz();
  REAL x2 = v2.getx();
  REAL y2 = v2.gety();
  REAL z2 = v2.getz();
  for (int pt = 0; pt < inputParticle.size(); ++pt) {
    Vec center = inputParticle[pt]->getCurrPos();
    if (center.getx() >= x1 && center.getx() < x2 &&
	center.gety() >= y1 && center.gety() < y2 &&
	center.getz() >= z1 && center.getz() < z2)
      foundParticle.push_back(inputParticle[pt]);
  }
}


REAL Assembly::getPtclMaxZ() const {
  std::vector<Particle*>::const_iterator it = allParticleVec.begin();
  REAL z0 = (*it)->getCurrPos().getz();
  for (; it != allParticleVec.end(); ++it) {
    if ( (*it)->getCurrPos().getz() > z0 )
      z0 = (*it)->getCurrPos().getz();
  }
  return z0;
}


void Assembly::partiCommuParticle() {

  // update allContainer dynamically due to particle motion
  if (mpiRank == 0) {
    setContainer(Rectangle(allContainer.getMinCorner().getx(),
			   allContainer.getMinCorner().gety(),
			   allContainer.getMinCorner().getz(),
			   allContainer.getMaxCorner().getx(),
			   allContainer.getMaxCorner().gety(),
			   getPtclMaxZ() + gradation.getPtclMaxRadius() ));
  }

  // broadcast
  broadcast(boostWorld, allContainer, 0);
  broadcast(boostWorld, gradation, 0);

  // partition particles and send to each process
  if (mpiRank == 0) { // process 0
    Vec v1 = allContainer.getMinCorner();
    Vec v2 = allContainer.getMaxCorner();
    Vec vspan = v2 - v1;

    for (int iRank = mpiSize - 1; iRank >= 0; --iRank) {
      particleVec.clear(); // do not release memory!
      int ndim = 3;
      int coords[3];
      MPI_Cart_coords(cartComm, iRank, ndim, coords);
      Rectangle container(v1.getx() + vspan.getx() / NPX * coords[0],
			  v1.gety() + vspan.gety() / NPY * coords[1],
			  v1.getz() + vspan.getz() / NPZ * coords[2],
			  v1.getx() + vspan.getx() / NPX * (coords[0] + 1),
			  v1.gety() + vspan.gety() / NPY * (coords[1] + 1),
			  v1.getz() + vspan.getz() / NPZ * (coords[2] + 1));
      findParticleInRectangle(container, allParticleVec, particleVec);
      if (iRank != 0)
	boostWorld.send(iRank, mpiTag, particleVec);
    }
  } else { // other processes except 0
    boostWorld.recv(0, mpiTag, particleVec); // need to release memory
  }

  // determine container of each process
  Vec v1 = allContainer.getMinCorner();
  Vec v2 = allContainer.getMaxCorner();
  Vec vspan = v2 - v1;
  container = Rectangle(v1.getx() + vspan.getx() / NPX * mpiCoords[0],
			v1.gety() + vspan.gety() / NPY * mpiCoords[1],
			v1.getz() + vspan.getz() / NPZ * mpiCoords[2],
			v1.getx() + vspan.getx() / NPX * (mpiCoords[0] + 1),
			v1.gety() + vspan.gety() / NPY * (mpiCoords[1] + 1),
			v1.getz() + vspan.getz() / NPZ * (mpiCoords[2] + 1));

  // find neighboring blocks
  int rankX1 = -1, rankX2 = -1, rankY1 = -1, rankY2 = -1, rankZ1 = -1, rankZ2 = -1;
  int rankX1Y1 = -1, rankX1Y2 = -1, rankX1Z1 = -1, rankX1Z2 = -1; 
  int rankX2Y1 = -1, rankX2Y2 = -1, rankX2Z1 = -1, rankX2Z2 = -1; 
  int rankY1Z1 = -1, rankY1Z2 = -1, rankY2Z1 = -1, rankY2Z2 = -1; 
  int rankX1Y1Z1 = -1, rankX1Y1Z2 = -1, rankX1Y2Z1 = -1, rankX1Y2Z2 = -1; 
  int rankX2Y1Z1 = -1, rankX2Y1Z2 = -1, rankX2Y2Z1 = -1, rankX2Y2Z2 = -1;
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

  // if found, communicate with neighboring blocks
  std::vector<Particle*> particleX1, particleX2;
  std::vector<Particle*> particleY1, particleY2;
  std::vector<Particle*> particleZ1, particleZ2;
  std::vector<Particle*> particleX1Y1, particleX1Y2, particleX1Z1, particleX1Z2; 
  std::vector<Particle*> particleX2Y1, particleX2Y2, particleX2Z1, particleX2Z2; 
  std::vector<Particle*> particleY1Z1, particleY1Z2, particleY2Z1, particleY2Z2; 
  std::vector<Particle*> particleX1Y1Z1, particleX1Y1Z2, particleX1Y2Z1, particleX1Y2Z2; 
  std::vector<Particle*> particleX2Y1Z1, particleX2Y1Z2, particleX2Y2Z1, particleX2Y2Z2; 
  boost::mpi::request reqX1[2], reqX2[2];
  boost::mpi::request reqY1[2], reqY2[2];
  boost::mpi::request reqZ1[2], reqZ2[2];
  boost::mpi::request reqX1Y1[2], reqX1Y2[2], reqX1Z1[2], reqX1Z2[2];
  boost::mpi::request reqX2Y1[2], reqX2Y2[2], reqX2Z1[2], reqX2Z2[2];
  boost::mpi::request reqY1Z1[2], reqY1Z2[2], reqY2Z1[2], reqY2Z2[2];
  boost::mpi::request reqX1Y1Z1[2], reqX1Y1Z2[2], reqX1Y2Z1[2], reqX1Y2Z2[2];
  boost::mpi::request reqX2Y1Z1[2], reqX2Y1Z2[2], reqX2Y2Z1[2], reqX2Y2Z2[2];
  v1 = container.getMinCorner(); // redefine v1, v2, vspan in terms of process
  v2 = container.getMaxCorner();   
  vspan = v2 - v1;
  REAL cellSize = gradation.getPtclMaxRadius() * 2;
  // 6 surfaces
  if (rankX1 >= 0) { // surface x1
    Rectangle containerX1(v1.getx(), v1.gety(), v1.getz(), 
			  v1.getx() + cellSize, v2.gety(), v2.getz());
    findParticleInRectangle(containerX1, particleVec, particleX1);
    reqX1[0] = boostWorld.isend(rankX1, mpiTag,  particleX1);
    reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rParticleX1);
  }
  if (rankX2 >= 0) { // surface x2
    Rectangle containerX2(v2.getx() - cellSize, v1.gety(), v1.getz(),
			  v2.getx(), v2.gety(), v2.getz());
    findParticleInRectangle(containerX2, particleVec, particleX2);
    reqX2[0] = boostWorld.isend(rankX2, mpiTag,  particleX2);
    reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rParticleX2);
  }
  if (rankY1 >= 0) {  // surface y1
    Rectangle containerY1(v1.getx(), v1.gety(), v1.getz(), 
			  v2.getx(), v1.gety() + cellSize, v2.getz());
    findParticleInRectangle(containerY1, particleVec, particleY1);
    reqY1[0] = boostWorld.isend(rankY1, mpiTag,  particleY1);
    reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rParticleY1);
  }
  if (rankY2 >= 0) {  // surface y2
    Rectangle containerY2(v1.getx(), v2.gety() - cellSize, v1.getz(),
			  v2.getx(), v2.gety(), v2.getz());
    findParticleInRectangle(containerY2, particleVec, particleY2);
    reqY2[0] = boostWorld.isend(rankY2, mpiTag,  particleY2);
    reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rParticleY2);
  }
  if (rankZ1 >= 0) {  // surface z1
    Rectangle containerZ1(v1.getx(), v1.gety(), v1.getz(),
			  v2.getx(), v2.gety(), v1.getz() + cellSize);
    findParticleInRectangle(containerZ1, particleVec, particleZ1);
    reqZ1[0] = boostWorld.isend(rankZ1, mpiTag,  particleZ1);
    reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rParticleZ1);
  }
  if (rankZ2 >= 0) {  // surface z2
    Rectangle containerZ2(v1.getx(), v1.gety(), v2.getz() - cellSize,
			  v2.getx(), v2.gety(), v2.getz());
    findParticleInRectangle(containerZ2, particleVec, particleZ2);
    reqZ2[0] = boostWorld.isend(rankZ2, mpiTag,  particleZ2);
    reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rParticleZ2);
  }
  // 12 edges
  if (rankX1Y1 >= 0) { // edge x1y1
    Rectangle containerX1Y1(v1.getx(), v1.gety(), v1.getz(),
			    v1.getx() + cellSize, v1.gety() + cellSize, v2.getz());
    findParticleInRectangle(containerX1Y1, particleVec, particleX1Y1);
    reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag,  particleX1Y1);
    reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rParticleX1Y1);
  }
  if (rankX1Y2 >= 0) { // edge x1y2
    Rectangle containerX1Y2(v1.getx(), v2.gety() - cellSize, v1.getz(),
			    v1.getx() + cellSize, v2.gety(), v2.getz());
    findParticleInRectangle(containerX1Y2, particleVec, particleX1Y2);
    reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag,  particleX1Y2);
    reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rParticleX1Y2);
  }
  if (rankX1Z1 >= 0) { // edge x1z1
    Rectangle containerX1Z1(v1.getx(), v1.gety(), v1.getz(),
			    v1.getx() + cellSize, v2.gety(), v1.getz() + cellSize);
    findParticleInRectangle(containerX1Z1, particleVec, particleX1Z1);
    reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag,  particleX1Z1);
    reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rParticleX1Z1);
  }
  if (rankX1Z2 >= 0) { // edge x1z2
    Rectangle containerX1Z2(v1.getx(), v1.gety(), v2.getz() - cellSize,
			    v1.getx() + cellSize, v2.gety(), v2.getz());
    findParticleInRectangle(containerX1Z2, particleVec, particleX1Z2);
    reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag,  particleX1Z2);
    reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rParticleX1Z2);
  }
  if (rankX2Y1 >= 0) { // edge x2y1
    Rectangle containerX2Y1(v2.getx() - cellSize, v1.gety(), v1.getz(),
			    v2.getx(), v1.gety() + cellSize, v2.getz());
    findParticleInRectangle(containerX2Y1, particleVec, particleX2Y1);
    reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag,  particleX2Y1);
    reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rParticleX2Y1);
  }
  if (rankX2Y2 >= 0) { // edge x2y2
    Rectangle containerX2Y2(v2.getx() - cellSize, v2.gety() - cellSize, v1.getz(),
			    v2.getx(), v2.gety(), v2.getz());
    findParticleInRectangle(containerX2Y2, particleVec, particleX2Y2);
    reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag,  particleX2Y2);
    reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rParticleX2Y2);
  }
  if (rankX2Z1 >= 0) { // edge x2z1
    Rectangle containerX2Z1(v2.getx() - cellSize, v1.gety(), v1.getz(),
			    v2.getx(), v2.gety(), v1.getz() + cellSize);
    findParticleInRectangle(containerX2Z1, particleVec, particleX2Z1);
    reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag,  particleX2Z1);
    reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rParticleX2Z1);
  }
  if (rankX2Z2 >= 0) { // edge x2z2
    Rectangle containerX2Z2(v2.getx() - cellSize, v1.gety(), v2.getz() - cellSize,
			    v2.getx(), v2.gety(), v2.getz());
    findParticleInRectangle(containerX2Z2, particleVec, particleX2Z2);
    reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag,  particleX2Z2);
    reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rParticleX2Z2);
  }
  if (rankY1Z1 >= 0) { // edge y1z1
    Rectangle containerY1Z1(v1.getx(), v1.gety(), v1.getz(),
			    v2.getx(), v1.gety() + cellSize, v1.getz() + cellSize);
    findParticleInRectangle(containerY1Z1, particleVec, particleY1Z1);
    reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag,  particleY1Z1);
    reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rParticleY1Z1);
  }
  if (rankY1Z2 >= 0) { // edge y1z2
    Rectangle containerY1Z2(v1.getx(), v1.gety(), v2.getz() - cellSize,
			    v2.getx(), v1.gety() + cellSize, v2.getz());
    findParticleInRectangle(containerY1Z2, particleVec, particleY1Z2);
    reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag,  particleY1Z2);
    reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rParticleY1Z2);
  }
  if (rankY2Z1 >= 0) { // edge y2z1
    Rectangle containerY2Z1(v1.getx(), v2.gety() - cellSize, v1.getz(),
			    v2.getx(), v2.gety(), v1.getz() + cellSize);
    findParticleInRectangle(containerY2Z1, particleVec, particleY2Z1);
    reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag,  particleY2Z1);
    reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rParticleY2Z1);
  }
  if (rankY2Z2 >= 0) { // edge y2z2
    Rectangle containerY2Z2(v1.getx(), v2.gety() - cellSize, v2.getz() - cellSize,
			    v2.getx(), v2.gety(), v2.getz());
    findParticleInRectangle(containerY2Z2, particleVec, particleY2Z2);
    reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag,  particleY2Z2);
    reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rParticleY2Z2);
  }
  // 8 vertices
  if (rankX1Y1Z1 >= 0) { // edge x1y1z1
    Rectangle containerX1Y1Z1(v1.getx(), v1.gety(), v1.getz(),
			      v1.getx() + cellSize, v1.gety() + cellSize, v1.getz() + cellSize);
    findParticleInRectangle(containerX1Y1Z1, particleVec, particleX1Y1Z1);
    reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag,  particleX1Y1Z1);
    reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rParticleX1Y1Z1);
  }
  if (rankX1Y1Z2 >= 0) { // edge x1y1z2
    Rectangle containerX1Y1Z2(v1.getx(), v1.gety(), v2.getz() - cellSize,
			      v1.getx() + cellSize, v1.gety() + cellSize, v2.getz());
    findParticleInRectangle(containerX1Y1Z2, particleVec, particleX1Y1Z2);
    reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag,  particleX1Y1Z2);
    reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rParticleX1Y1Z2);
  }
  if (rankX1Y2Z1 >= 0) { // edge x1y2z1
    Rectangle containerX1Y2Z1(v1.getx(), v2.gety() - cellSize, v1.getz(),
			      v1.getx() + cellSize, v2.gety(), v1.getz() + cellSize);
    findParticleInRectangle(containerX1Y2Z1, particleVec, particleX1Y2Z1);
    reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag,  particleX1Y2Z1);
    reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rParticleX1Y2Z1);
  }
  if (rankX1Y2Z2 >= 0) { // edge x1y2z2
    Rectangle containerX1Y2Z2(v1.getx(), v2.gety() - cellSize, v2.getz() - cellSize,
			      v1.getx() + cellSize, v2.gety() + cellSize, v2.getz());
    findParticleInRectangle(containerX1Y2Z2, particleVec, particleX1Y2Z2);
    reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag,  particleX1Y2Z2);
    reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rParticleX1Y2Z2);
  }
  if (rankX2Y1Z1 >= 0) { // edge x2y1z1
    Rectangle containerX2Y1Z1(v2.getx() - cellSize, v1.gety(), v1.getz(),
			      v2.getx(), v1.gety() + cellSize, v1.getz() + cellSize);
    findParticleInRectangle(containerX2Y1Z1, particleVec, particleX2Y1Z1);
    reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag,  particleX2Y1Z1);
    reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rParticleX2Y1Z1);
  }
  if (rankX2Y1Z2 >= 0) { // edge x2y1z2
    Rectangle containerX2Y1Z2(v2.getx() - cellSize, v1.gety(), v2.getz() - cellSize,
			      v2.getx(), v1.gety() + cellSize, v2.getz());
    findParticleInRectangle(containerX2Y1Z2, particleVec, particleX2Y1Z2);
    reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag,  particleX2Y1Z2);
    reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rParticleX2Y1Z2);
  }
  if (rankX2Y2Z1 >= 0) { // edge x2y2z1
    Rectangle containerX2Y2Z1(v2.getx() - cellSize, v2.gety() - cellSize, v1.getz(),
			      v2.getx(), v2.gety(), v1.getz() + cellSize);
    findParticleInRectangle(containerX2Y2Z1, particleVec, particleX2Y2Z1);
    reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag,  particleX2Y2Z1);
    reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rParticleX2Y2Z1);
  }
  if (rankX2Y2Z2 >= 0) { // edge x2y2z2
    Rectangle containerX2Y2Z2(v2.getx() - cellSize, v2.gety() - cellSize, v2.getz() - cellSize,
			      v2.getx(), v2.gety(), v2.getz());
    findParticleInRectangle(containerX2Y2Z2, particleVec, particleX2Y2Z2);
    reqX2Y2Z2[0] = boostWorld.isend(rankX2Y2Z2, mpiTag,  particleX2Y2Z2);
    reqX2Y2Z2[1] = boostWorld.irecv(rankX2Y2Z2, mpiTag, rParticleX2Y2Z2);
  }
  // 6 surfaces
  if (rankX1 >= 0) boost::mpi::wait_all(reqX1, reqX1 + 2);
  if (rankX2 >= 0) boost::mpi::wait_all(reqX2, reqX2 + 2);
  if (rankY1 >= 0) boost::mpi::wait_all(reqY1, reqY1 + 2);
  if (rankY2 >= 0) boost::mpi::wait_all(reqY2, reqY2 + 2);
  if (rankZ1 >= 0) boost::mpi::wait_all(reqZ1, reqZ1 + 2);
  if (rankZ2 >= 0) boost::mpi::wait_all(reqZ2, reqZ2 + 2);
  // 12 edges
  if (rankX1Y1 >= 0) boost::mpi::wait_all(reqX1Y1, reqX1Y1 + 2);
  if (rankX1Y2 >= 0) boost::mpi::wait_all(reqX1Y2, reqX1Y2 + 2);  
  if (rankX1Z1 >= 0) boost::mpi::wait_all(reqX1Z1, reqX1Z1 + 2);
  if (rankX1Z2 >= 0) boost::mpi::wait_all(reqX1Z2, reqX1Z2 + 2);
  if (rankX2Y1 >= 0) boost::mpi::wait_all(reqX2Y1, reqX2Y1 + 2);
  if (rankX2Y2 >= 0) boost::mpi::wait_all(reqX2Y2, reqX2Y2 + 2);  
  if (rankX2Z1 >= 0) boost::mpi::wait_all(reqX2Z1, reqX2Z1 + 2);
  if (rankX2Z2 >= 0) boost::mpi::wait_all(reqX2Z2, reqX2Z2 + 2); 
  if (rankY1Z1 >= 0) boost::mpi::wait_all(reqY1Z1, reqY1Z1 + 2);
  if (rankY1Z2 >= 0) boost::mpi::wait_all(reqY1Z2, reqY1Z2 + 2);
  if (rankY2Z1 >= 0) boost::mpi::wait_all(reqY2Z1, reqY2Z1 + 2);
  if (rankY2Z2 >= 0) boost::mpi::wait_all(reqY2Z2, reqY2Z2 + 2); 
  // 8 vertices
  if (rankX1Y1Z1 >= 0) boost::mpi::wait_all(reqX1Y1Z1, reqX1Y1Z1 + 2);
  if (rankX1Y1Z2 >= 0) boost::mpi::wait_all(reqX1Y1Z2, reqX1Y1Z2 + 2);
  if (rankX1Y2Z1 >= 0) boost::mpi::wait_all(reqX1Y2Z1, reqX1Y2Z1 + 2);
  if (rankX1Y2Z2 >= 0) boost::mpi::wait_all(reqX1Y2Z2, reqX1Y2Z2 + 2);
  if (rankX2Y1Z1 >= 0) boost::mpi::wait_all(reqX2Y1Z1, reqX2Y1Z1 + 2);
  if (rankX2Y1Z2 >= 0) boost::mpi::wait_all(reqX2Y1Z2, reqX2Y1Z2 + 2);
  if (rankX2Y2Z1 >= 0) boost::mpi::wait_all(reqX2Y2Z1, reqX2Y2Z1 + 2);
  if (rankX2Y2Z2 >= 0) boost::mpi::wait_all(reqX2Y2Z2, reqX2Y2Z2 + 2);  

  // merge: particles inside container (at front) + particles from neighoring blocks (at end)
  recvParticleVec.insert(recvParticleVec.end(), particleVec.begin(), particleVec.end());
  // 6 surfaces
  recvParticleVec.insert(recvParticleVec.end(), rParticleX1.begin(), rParticleX1.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleX2.begin(), rParticleX2.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleY1.begin(), rParticleY1.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleY2.begin(), rParticleY2.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleZ1.begin(), rParticleZ1.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleZ2.begin(), rParticleZ2.end());
  // 12 edges
  recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1.begin(), rParticleX1Y1.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2.begin(), rParticleX1Y2.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z1.begin(), rParticleX1Z1.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z2.begin(), rParticleX1Z2.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1.begin(), rParticleX2Y1.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2.begin(), rParticleX2Y2.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z1.begin(), rParticleX2Z1.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z2.begin(), rParticleX2Z2.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z1.begin(), rParticleY1Z1.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z2.begin(), rParticleY1Z2.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z1.begin(), rParticleY2Z1.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z2.begin(), rParticleY2Z2.end());
  // 8 vertices
  recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z1.begin(), rParticleX1Y1Z1.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z2.begin(), rParticleX1Y1Z2.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z1.begin(), rParticleX1Y2Z1.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z2.begin(), rParticleX1Y2Z2.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z1.begin(), rParticleX2Y1Z1.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z2.begin(), rParticleX2Y1Z2.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z1.begin(), rParticleX2Y2Z1.end());
  recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z2.begin(), rParticleX2Y2Z2.end());

  ///*
  std::cout << "iter=" << iteration << " rank=" << setw(4) << mpiRank 
	    << " particleNum=" << setw(4) << particleVec.size() << " surface="
	    << setw(4) << particleX1.size()  << setw(4) << particleX2.size()
	    << setw(4) << particleY1.size()  << setw(4) << particleY2.size()
	    << setw(4) << particleZ1.size()  << setw(4) << particleZ2.size()  << " recv_surface="
	    << setw(4) << rParticleX1.size() << setw(4) << rParticleX2.size()
	    << setw(4) << rParticleY1.size() << setw(4) << rParticleY2.size()
	    << setw(4) << rParticleZ1.size() << setw(4) << rParticleZ2.size() << " rNum="    
	    << setw(4) << recvParticleVec.size()   
	    << "\n" << std::flush;

  //for (std::vector<Particle*>::const_iterator it = recvParticleVec.begin(); it != recvParticleVec.end();++it)
  //  cout << (*it)->getId() << ' ' << flush;
  //cout << "\n" << flush;
  //*/
}


void Assembly::gatherParticle() {
  // release memory
  // 6 surfaces
  std::vector<Particle*>::iterator it;
  for (it = rParticleX1.begin(); it != rParticleX1.end(); ++it)
    delete (*it);
  for (it = rParticleX2.begin(); it != rParticleX2.end(); ++it)
    delete (*it);
  for (it = rParticleY1.begin(); it != rParticleY1.end(); ++it)
    delete (*it);
  for (it = rParticleY2.begin(); it != rParticleY2.end(); ++it)
    delete (*it);
  for (it = rParticleZ1.begin(); it != rParticleZ1.end(); ++it)
    delete (*it);
  for (it = rParticleZ2.begin(); it != rParticleZ2.end(); ++it)
    delete (*it);
  // 12 edges
  for (it = rParticleX1Y1.begin(); it != rParticleX1Y1.end(); ++it)
    delete (*it);
  for (it = rParticleX1Y2.begin(); it != rParticleX1Y2.end(); ++it)
    delete (*it);
  for (it = rParticleX1Z1.begin(); it != rParticleX1Z1.end(); ++it)
    delete (*it);
  for (it = rParticleX1Z2.begin(); it != rParticleX1Z2.end(); ++it)
    delete (*it);
  for (it = rParticleX2Y1.begin(); it != rParticleX2Y1.end(); ++it)
    delete (*it);
  for (it = rParticleX2Y2.begin(); it != rParticleX2Y2.end(); ++it)
    delete (*it);
  for (it = rParticleX2Z1.begin(); it != rParticleX2Z1.end(); ++it)
    delete (*it);
  for (it = rParticleX2Z2.begin(); it != rParticleX2Z2.end(); ++it)
    delete (*it);
  for (it = rParticleY1Z1.begin(); it != rParticleY1Z1.end(); ++it)
    delete (*it);
  for (it = rParticleY1Z2.begin(); it != rParticleY1Z2.end(); ++it)
    delete (*it);
  for (it = rParticleY2Z1.begin(); it != rParticleY2Z1.end(); ++it)
    delete (*it);
  for (it = rParticleY2Z2.begin(); it != rParticleY2Z2.end(); ++it)
    delete (*it);
  // 8 vertices
  for (it = rParticleX1Y1Z1.begin(); it != rParticleX1Y1Z1.end(); ++it)
    delete (*it);
  for (it = rParticleX1Y1Z2.begin(); it != rParticleX1Y1Z2.end(); ++it)
    delete (*it);
  for (it = rParticleX1Y2Z1.begin(); it != rParticleX1Y2Z1.end(); ++it)
    delete (*it);
  for (it = rParticleX1Y2Z2.begin(); it != rParticleX1Y2Z2.end(); ++it)
    delete (*it);
  for (it = rParticleX2Y1Z1.begin(); it != rParticleX2Y1Z1.end(); ++it)
    delete (*it);
  for (it = rParticleX2Y1Z2.begin(); it != rParticleX2Y1Z2.end(); ++it)
    delete (*it);
  for (it = rParticleX2Y2Z1.begin(); it != rParticleX2Y2Z1.end(); ++it)
    delete (*it);
  for (it = rParticleX2Y2Z2.begin(); it != rParticleX2Y2Z2.end(); ++it)
    delete (*it);
  // 6 surfaces
  rParticleX1.clear();
  rParticleX2.clear();
  rParticleY1.clear();
  rParticleY2.clear();
  rParticleZ1.clear();
  rParticleZ2.clear();
  // 12 edges
  rParticleX1Y1.clear();
  rParticleX1Y2.clear();
  rParticleX1Z1.clear();
  rParticleX1Z2.clear();
  rParticleX2Y1.clear();
  rParticleX2Y2.clear();
  rParticleX2Z1.clear();
  rParticleX2Z2.clear();
  rParticleY1Z1.clear();
  rParticleY1Z2.clear();
  rParticleY2Z1.clear();
  rParticleY2Z2.clear();
  // 8 vertices
  rParticleX1Y1Z1.clear();
  rParticleX1Y1Z2.clear();
  rParticleX1Y2Z1.clear();
  rParticleX1Y2Z2.clear();
  rParticleX2Y1Z1.clear();
  rParticleX2Y1Z2.clear();
  rParticleX2Y2Z1.clear();
  rParticleX2Y2Z2.clear();

  recvParticleVec.clear();

  // update allParticleVec: process 0 collects all updated particles from each process  
  if (mpiRank != 0) { // each process except 0
    boostWorld.send(0, mpiTag, particleVec);
    // release memory
    for (it = particleVec.begin(); it != particleVec.end(); ++it)
      delete (*it);
    particleVec.clear();
  } else { // process 0

    // copy particleVec
    std::vector<Particle*> copyParticleVec;
    for (int i = 0; i < particleVec.size(); ++i)
      copyParticleVec.push_back( new Particle(*particleVec[i]) );

    // clear allParticleVec
    std::vector<Particle*>::iterator it;
    for(it = allParticleVec.begin(); it != allParticleVec.end(); ++it)
      delete (*it);
    allParticleVec.clear();

    // fill allParticleVec with copyparticleVec and received particles
    allParticleVec.insert(allParticleVec.end(), copyParticleVec.begin(), copyParticleVec.end());
    copyParticleVec.clear();  // do not release memory!
    for (int iRank = 1; iRank < mpiSize; ++iRank) {
      recvParticleVec.clear();// do not release memory!
      boostWorld.recv(iRank, mpiTag, recvParticleVec); // new use of recvParticleVec
      allParticleVec.insert(allParticleVec.end(), recvParticleVec.begin(), recvParticleVec.end());
    }
    recvParticleVec.clear();  // do not release memory!

  }
  
}


void Assembly::readParticle(const char *inputParticle){
  std::ifstream ifs(inputParticle);
  if(!ifs) {
    cout << "stream error!" << endl; exit(-1);
  }
  int particleNum;
  ifs >> particleNum;
  std::string str;
  ifs >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str
      >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str
      >> str >> str >> str >> str >> str >> str >> str >> str >> str;
  
  std::vector<Particle*>::iterator it;
  for(it = allParticleVec.begin(); it != allParticleVec.end(); ++it)
    delete (*it);
  allParticleVec.clear();

  int id, type;
  REAL a, b, c, px, py, pz, dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz;
  REAL vx, vy, vz, omx, omy, omz, fx, fy, fz, mx, my, mz;
  for (int i = 0; i < particleNum; ++i){
    ifs >> id >> type >> a >> b >> c >> px >> py >> pz >> dax >> day >> daz >> dbx >> dby >> dbz >> dcx >> dcy >> dcz
	>> vx >> vy >> vz >> omx >> omy >> omz >> fx >> fy >> fz >> mx >> my >> mz;
    Particle* pt= new Particle(id,type,Vec(a,b,c),Vec(px,py,pz),Vec(dax,day,daz),Vec(dbx,dby,dbz),Vec(dcx,dcy,dcz),YOUNG,POISSON);
    
    // optional settings for a particle's initial status
    //pt->setPrevVelocity(Vec(vx,vy,vz));
    //pt->setCurrVelocity(Vec(vx,vy,vz));
    //pt->setPrevOmga(Vec(omx,omy,omz));
    //pt->setCurrOmga(Vec(omx,omy,omz));
    
    //pt->setForce(Vec(fx,fy,fz));  // initial force
    //pt->setMoment(Vec(mx,my,mz)); // initial moment
    
    //pt->setConstForce(Vec(fx,fy,fz));  // constant force, not initial force
    //pt->setConstMoment(Vec(mx,my,mz)); // constant moment, not initial moment

    allParticleVec.push_back(pt);
  }
  
  int sieveNum;
  ifs >> sieveNum;
  std::vector<REAL> percent(sieveNum), size(sieveNum);
  for (int i = 0; i < sieveNum; ++i)
    ifs >> percent[i] >> size[i];
  REAL ratio_ba, ratio_ca;
  ifs >> ratio_ba >> ratio_ca;
  setGradation(Gradation(sieveNum, percent, size, ratio_ba, ratio_ca));
  ifs.close();
}


void Assembly::printParticle(const char* str) const {
  std::ofstream ofs(str);
  if(!ofs) {
    cout << "stream error!" << endl; exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  ofs << setw(OWID) << allParticleVec.size() << endl;
  ofs << setw(OWID) << "id"
      << setw(OWID) << "type"
      << setw(OWID) << "radius_a"
      << setw(OWID) << "radius_b"
      << setw(OWID) << "radius_c"
      << setw(OWID) << "position_x"
      << setw(OWID) << "position_y"
      << setw(OWID) << "position_z"
      << setw(OWID) << "axle_a_x"
      << setw(OWID) << "axle_a_y"
      << setw(OWID) << "axle_a_z"
      << setw(OWID) << "axle_b_x"
      << setw(OWID) << "axle_b_y"
      << setw(OWID) << "axle_b_z"
      << setw(OWID) << "axle_c_x"
      << setw(OWID) << "axle_c_y"
      << setw(OWID) << "axle_c_z"
      << setw(OWID) << "velocity_x"
      << setw(OWID) << "velocity_y"
      << setw(OWID) << "velocity_z"
      << setw(OWID) << "omga_x"
      << setw(OWID) << "omga_y"
      << setw(OWID) << "omga_z"
      << setw(OWID) << "force_x"
      << setw(OWID) << "force_y"
      << setw(OWID) << "force_z"
      << setw(OWID) << "moment_x"
      << setw(OWID) << "moment_y"
      << setw(OWID) << "moment_z"
      << endl;
  
  Vec vObj;
  std::vector<Particle*>::const_iterator  it;
  for (it = allParticleVec.begin(); it != allParticleVec.end();++it)  {
    ofs << setw(OWID) << (*it)->getId()
	<< setw(OWID) << (*it)->getType()
	<< setw(OWID) << (*it)->getA()
	<< setw(OWID) << (*it)->getB()
	<< setw(OWID) << (*it)->getC();
    
    vObj=(*it)->getCurrPos();
    ofs << setw(OWID) << vObj.getx()
	<< setw(OWID) << vObj.gety()
	<< setw(OWID) << vObj.getz();
    
    vObj=(*it)->getCurrDirecA();
    ofs << setw(OWID) << vObj.getx()
	<< setw(OWID) << vObj.gety()
	<< setw(OWID) << vObj.getz();
    
    vObj=(*it)->getCurrDirecB();
    ofs << setw(OWID) << vObj.getx()
	<< setw(OWID) << vObj.gety()
	<< setw(OWID) << vObj.getz();
    
    vObj=(*it)->getCurrDirecC();
    ofs << setw(OWID) << vObj.getx()
	<< setw(OWID) << vObj.gety()
	<< setw(OWID) << vObj.getz();
    
    vObj=(*it)->getCurrVeloc();
    ofs << setw(OWID) << vObj.getx()
	<< setw(OWID) << vObj.gety()
	<< setw(OWID) << vObj.getz();
    
    vObj=(*it)->getCurrOmga();
    ofs << setw(OWID) << vObj.getx()
	<< setw(OWID) << vObj.gety()
	<< setw(OWID) << vObj.getz();
    
    vObj=(*it)->getForce();
    ofs << setw(OWID) << vObj.getx()
	<< setw(OWID) << vObj.gety()
	<< setw(OWID) << vObj.getz();
    
    vObj=(*it)->getMoment();
    ofs << setw(OWID) << vObj.getx()
	<< setw(OWID) << vObj.gety()
	<< setw(OWID) << vObj.getz() << endl;
  }
  
  int sieveNum = gradation.getSieveNum();
  std::vector<REAL> percent = gradation.getPercent();
  std::vector<REAL> size    = gradation.getSize();
  ofs << endl << setw(OWID) << sieveNum << endl;
  for (int i = 0; i < sieveNum; ++i)
    ofs << setw(OWID) << percent[i] << setw(OWID) << size[i] << endl;
  ofs << endl << setw(OWID) << gradation.getPtclRatioBA() << setw(OWID) << gradation.getPtclRatioCA() << endl;

  ofs.close();
}


void Assembly::readBoundary(const char* str) {
  std::ifstream ifs(str);
  if(!ifs) {
    cout << "stream error!" << endl; exit(-1);
  }

  REAL x1, y1, z1, x2, y2, z2;
  ifs >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;
  setContainer(Rectangle(x1, y1, z1, x2, y2, z2)); // update allContainer

  Boundary<Particle>* rbptr;
  int type;
  boundaryVec.clear();
  int boundaryNum;
  ifs >> boundaryNum;
  for(int i = 0; i < boundaryNum; ++i){
    ifs >> type;
    if(type == 1) // plane boundary
      rbptr = new plnBoundary<Particle>(ifs);
    else          // cylindrical boundary
      rbptr = new cylBoundary<Particle>(ifs);
    boundaryVec.push_back(rbptr);
  }

  ifs.close();
}


void Assembly::printBoundary(const char* str) const {
  std::ofstream ofs(str);
  if(!ofs) {cout << "stream error!" << endl; exit(-1);}
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  
  Vec v1 = allContainer.getMinCorner();
  Vec v2 = allContainer.getMaxCorner();
  REAL x1 = v1.getx();
  REAL y1 = v1.gety();
  REAL z1 = v1.getz();
  REAL x2 = v2.getx();
  REAL y2 = v2.gety();
  REAL z2 = v2.getz();
  
  ofs << setw(OWID) << x1 << setw(OWID) << y1 << setw(OWID) << z1
      << setw(OWID) << x2 << setw(OWID) << y2 << setw(OWID) << z2 << endl
      << setw(OWID) << boundaryVec.size() << endl;
  
  std::vector<BOUNDARY*>::const_iterator rt;
  for(rt=boundaryVec.begin();rt!=boundaryVec.end();++rt)
    (*rt)->display(ofs);
  ofs << endl;
  
  ofs.close();
}


void Assembly::plotBoundary(const char *str) const {
  std::ofstream ofs(str);
  if(!ofs) {
    cout << "stream error!" << endl; exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  REAL x1,x2,y1,y2,z1,z2,x0,y0,z0;
  x1 = allContainer.getMinCorner().getx();
  y1 = allContainer.getMinCorner().gety();
  z1 = allContainer.getMinCorner().getz();
  x2 = allContainer.getMaxCorner().getx();
  y2 = allContainer.getMaxCorner().gety();
  z2 = allContainer.getMaxCorner().getz();

  ofs << "ZONE N=8, E=1, DATAPACKING=POINT, ZONETYPE=FEBRICK" << endl;
  ofs << setw(OWID) << x2 << setw(OWID) << y1 << setw(OWID) << z1 << endl;
  ofs << setw(OWID) << x2 << setw(OWID) << y2 << setw(OWID) << z1 << endl;
  ofs << setw(OWID) << x1 << setw(OWID) << y2 << setw(OWID) << z1 << endl;
  ofs << setw(OWID) << x1 << setw(OWID) << y1 << setw(OWID) << z1 << endl;
  ofs << setw(OWID) << x2 << setw(OWID) << y1 << setw(OWID) << z2 << endl;
  ofs << setw(OWID) << x2 << setw(OWID) << y2 << setw(OWID) << z2 << endl;
  ofs << setw(OWID) << x1 << setw(OWID) << y2 << setw(OWID) << z2 << endl;
  ofs << setw(OWID) << x1 << setw(OWID) << y1 << setw(OWID) << z2 << endl;
  ofs << "1 2 3 4 5 6 7 8" << endl;

  ofs.close();
}


void Assembly::findContact() { // serial version, O(n x n), n is the number of particles.
  contactVec.clear();
  possContactNum = 0;
  
#ifdef TIME_PROFILE
  double time_r = 0; // time consumed in contact resolution, i.e., tmpct.isOverlapped()
  gettimeofday(&time_p1,NULL); 
#endif
  
  int num1 = particleVec.size();  // particles inside container
  int num2 = recvParticleVec.size(); // particles inside container (at front) + particles from neighboring blocks (at end)
  for (int i = 0; i < num1 - 1; ++i) {
    Vec u = particleVec[i]->getCurrPos();
    for (int j = i + 1; j < num2; ++j){
      Vec v = recvParticleVec[j]->getCurrPos();
      if (  ( vfabs(v - u) < particleVec[i]->getA() + recvParticleVec[j]->getA())
	    && ( particleVec[i]->getType() !=  1 || recvParticleVec[j]->getType() != 1  )      // not both are fixed particles
	    && ( particleVec[i]->getType() !=  5 || recvParticleVec[j]->getType() != 5  )      // not both are free boundary particles
	    && ( particleVec[i]->getType() != 10 || recvParticleVec[j]->getType() != 10 )  ) { // not both are ghost particles
	Contact<Particle> tmpct(particleVec[i], recvParticleVec[j]); // a local and temparory object
	++possContactNum;
#ifdef TIME_PROFILE
	gettimeofday(&time_r1,NULL); 
#endif
	if(tmpct.isOverlapped())
	  contactVec.push_back(tmpct);    // containers use value semantics, so a "copy" is pushed back.
#ifdef TIME_PROFILE
	gettimeofday(&time_r2,NULL); 
	time_r += timediffsec(time_r1, time_r2);
#endif
      }
    }
  }	
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf << setw(OWID) << "findContact=" << setw(OWID) << timediffsec(time_p1, time_p2) << setw(OWID) << "isOverlapped=" << setw(OWID) << time_r; 
#endif
  
  actualContactNum = contactVec.size();
}


void Assembly::internalForce(){
  REAL avgNormal, avgShear;
  
  int contactNum = contactVec.size();
  if(contactNum == 0){
    avgNormal = 0;
    avgShear = 0;
  }
  else{
    for (std::vector<CONTACT>::iterator it = contactVec.begin(); it != contactVec.end(); ++it)
      it->checkinPreTgt(contactTgtVec); // checkin previous tangential force and displacment    
    
    contactTgtVec.clear(); // contactTgtVec must be cleared before filling in new values.
    
#ifdef TIME_PROFILE
    gettimeofday(&time_p1,NULL); 
#endif 
    for (std::vector<CONTACT>::iterator it = contactVec.begin(); it != contactVec.end(); ++ it){
      it->contactForce();             // cannot be parallelized as it may change a particle's force simultaneously.
      it->checkoutTgt(contactTgtVec); // checkout current tangential force and displacment
      avgNormal += it->getNormalForce();
      avgShear += it->getTgtForce();
    }
    avgNormal /= contactNum;
    avgShear /= contactNum;
    
#ifdef TIME_PROFILE
    gettimeofday(&time_p2,NULL);
    debugInf << setw(OWID) << "internalForce=" << setw(OWID) << timediffsec(time_p1, time_p2) << endl; 
#endif
    
  }
}


void Assembly::updateParticle() {
  for(std::vector<Particle*>::iterator it = particleVec.begin(); it != particleVec.end(); ++it)
    (*it)->update();
}


void Assembly::clearForce() {
  for(std::vector<Particle*>::iterator it = allParticleVec.begin(); it != allParticleVec.end(); ++it)
    (*it)->clearForce();
}


void Assembly::findParticleOnBoundary() {
  for(std::vector<BOUNDARY*>::iterator rt = boundaryVec.begin(); rt != boundaryVec.end(); ++rt)
    (*rt)->findParticleOnBoundary(allParticleVec);
}


void Assembly::boundaryForce() {
  for(std::vector<BOUNDARY*>::iterator rt=boundaryVec.begin();rt!=boundaryVec.end();++rt)
    (*rt)->boundaryForce(boundaryTgtMap);
}


// rule out
/*
void Assembly::plotCavity(const char *str) const {
  std::ofstream ofs(str);
  if(!ofs) {
    cout << "stream error!" << endl; exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  REAL x1,x2,y1,y2,z1,z2,x0,y0,z0;
  x1 = cavity.getMinCorner().getx();
  y1 = cavity.getMinCorner().gety();
  z1 = cavity.getMinCorner().getz();
  x2 = cavity.getMaxCorner().getx();
  y2 = cavity.getMaxCorner().gety();
  z2 = cavity.getMaxCorner().getz();

  ofs << "ZONE N=8, E=1, DATAPACKING=POINT, ZONETYPE=FEBRICK" << endl;
  ofs << setw(OWID) << x2 << setw(OWID) << y1 << setw(OWID) << z1 << endl;
  ofs << setw(OWID) << x2 << setw(OWID) << y2 << setw(OWID) << z1 << endl;
  ofs << setw(OWID) << x1 << setw(OWID) << y2 << setw(OWID) << z1 << endl;
  ofs << setw(OWID) << x1 << setw(OWID) << y1 << setw(OWID) << z1 << endl;
  ofs << setw(OWID) << x2 << setw(OWID) << y1 << setw(OWID) << z2 << endl;
  ofs << setw(OWID) << x2 << setw(OWID) << y2 << setw(OWID) << z2 << endl;
  ofs << setw(OWID) << x1 << setw(OWID) << y2 << setw(OWID) << z2 << endl;
  ofs << setw(OWID) << x1 << setw(OWID) << y1 << setw(OWID) << z2 << endl;
  ofs << "1 2 3 4 5 6 7 8" << endl;

  ofs.close();
}

void Assembly::plotSpring(const char *str) const {
  std::ofstream ofs(str);
  if(!ofs) {
    cout << "stream error!" << endl; exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  int totalMemParticle = 0;
  for (int i = 0; i < memBoundary.size(); ++i) 
    for (int j = 0; j < memBoundary[i].size(); ++j) 
      for (int k = 0; k < memBoundary[i][j].size(); ++k) 
	++totalMemParticle;
  int totalSpring = springVec.size();
  ofs << "ZONE N=" << totalMemParticle << ", E=" << totalSpring << ", DATAPACKING=POINT, ZONETYPE=FELINESEG" << endl;
  Particle *pt = NULL;
  Vec vt;
  for (int i = 0; i < memBoundary.size(); ++i) 
    for (int j = 0; j < memBoundary[i].size(); ++j) 
      for (int k = 0; k < memBoundary[i][j].size(); ++k) {
	pt = memBoundary[i][j][k]; 
	vt = pt->getCurrPos();
	ofs << setw(OWID) << vt.getx() << setw(OWID) << vt.gety() << setw(OWID) << vt.getz() << endl;
      }
  for (int i = 0; i < springVec.size(); ++i) {
    ofs << setw(OWID) << springVec[i]->getParticleId1() - trimHistoryNum  << setw(OWID) << springVec[i]->getParticleId2() - trimHistoryNum << endl;
  }

  ofs.close();
}  

void Assembly::printMemParticle(const char* str) const  {
  std::ofstream ofs(str);
  if(!ofs) {
    cout << "stream error!" << endl; exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  
  int totalMemParticle = 0;
  for (int i = 0; i < memBoundary.size(); ++i) 
    for (int j = 0; j < memBoundary[i].size(); ++j) 
      for (int k = 0; k < memBoundary[i][j].size(); ++k) 
	++totalMemParticle;
  
  ofs << setw(OWID) << totalMemParticle << setw(OWID) << 1 << endl;
  ofs << setw(OWID) << container.getCenter().getx()
      << setw(OWID) << container.getCenter().gety()
      << setw(OWID) << container.getCenter().getz()
      << setw(OWID) << container.getDimx()
      << setw(OWID) << container.getDimy()
      << setw(OWID) << container.getDimz() << endl;
  
  ofs << setw(OWID) << "ID"
      << setw(OWID) << "type"
      << setw(OWID) << "radius_a"
      << setw(OWID) << "radius_b"
      << setw(OWID) << "radius_c"
      << setw(OWID) << "position_x"
      << setw(OWID) << "position_y"
      << setw(OWID) << "position_z"
      << setw(OWID) << "axle_a_x"
      << setw(OWID) << "axle_a_y"
      << setw(OWID) << "axle_a_z"
      << setw(OWID) << "axle_b_x"
      << setw(OWID) << "axle_b_y"
      << setw(OWID) << "axle_b_z"
      << setw(OWID) << "axle_c_x"
      << setw(OWID) << "axle_c_y"
      << setw(OWID) << "axle_c_z"
      << setw(OWID) << "velocity_x"
      << setw(OWID) << "velocity_y"
      << setw(OWID) << "velocity_z"
      << setw(OWID) << "omga_x"
      << setw(OWID) << "omga_y"
      << setw(OWID) << "omga_z"
      << setw(OWID) << "force_x"
      << setw(OWID) << "force_y"
      << setw(OWID) << "force_z"
      << setw(OWID) << "moment_x"
      << setw(OWID) << "moment_y"
      << setw(OWID) << "moment_z"
      << endl;
  
  Particle *it = NULL;
  Vec vObj;
  for (int i = 0; i < memBoundary.size(); ++i) 
    for (int j = 0; j < memBoundary[i].size(); ++j) 
      for (int k = 0; k < memBoundary[i][j].size(); ++k) {
	it = memBoundary[i][j][k];
	ofs << setw(OWID) << it->getId()
	    << setw(OWID) << it->getType()
	    << setw(OWID) << it->getA()
	    << setw(OWID) << it->getB()
	    << setw(OWID) << it->getC();
	
	vObj=it->getCurrPos();
	ofs << setw(OWID) << vObj.getx()
	    << setw(OWID) << vObj.gety()
	    << setw(OWID) << vObj.getz();
	
	vObj=it->getCurrDirecA();
	ofs << setw(OWID) << vObj.getx()
	    << setw(OWID) << vObj.gety()
	    << setw(OWID) << vObj.getz();
	
	vObj=it->getCurrDirecB();
	ofs << setw(OWID) << vObj.getx()
	    << setw(OWID) << vObj.gety()
	    << setw(OWID) << vObj.getz();
	
	vObj=it->getCurrDirecC();
	ofs << setw(OWID) << vObj.getx()
	    << setw(OWID) << vObj.gety()
	    << setw(OWID) << vObj.getz();
	
	vObj=it->getCurrVeloc();
	ofs << setw(OWID) << vObj.getx()
	    << setw(OWID) << vObj.gety()
	    << setw(OWID) << vObj.getz();
	
	vObj=it->getCurrOmga();
	ofs << setw(OWID) << vObj.getx()
	    << setw(OWID) << vObj.gety()
	    << setw(OWID) << vObj.getz();
	
	vObj=it->getForce();
	ofs << setw(OWID) << vObj.getx()
	    << setw(OWID) << vObj.gety()
	    << setw(OWID) << vObj.getz();
	
	vObj=it->getMoment();
	ofs << setw(OWID) << vObj.getx()
	    << setw(OWID) << vObj.gety()
	      << setw(OWID) << vObj.getz() << endl;
      }
  ofs.close();  
}
 
// vector elements are in the order of:
// x1: inner, outer
// x2: inner, outer
// y1: inner, outer
// y2: inner, outer
// z1: inner, outer
// z2: inner, outer
void Assembly::checkMembrane(vector<REAL> &vx ) const {
  vector<Particle*> vec1d;  // 1-dimension
  vector< vector<Particle*>  > vec2d; // 2-dimension
  REAL in, out, tmp;
  REAL x1_in, x1_out, x2_in, x2_out;
  REAL y1_in, y1_out, y2_in, y2_out;
  REAL z1_in, z1_out, z2_in, z2_out;

  // surface x1
  vec2d = memBoundary[0];
  in = vec2d[0][0]->getCurrPos().getx();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPos().getx();
      if (tmp < out) out = tmp;
      if (tmp > in ) in  = tmp;
    }
  vx.push_back(in);
  vx.push_back(out);
  x1_in  = in;
  x1_out = out;

  // surface x2
  vec2d.clear();
  vec2d = memBoundary[1];
  in = vec2d[0][0]->getCurrPos().getx();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPos().getx();
      if (tmp > out) out = tmp;
      if (tmp < in ) in  = tmp;
    }
  vx.push_back(in);
  vx.push_back(out);
  x2_in  = in;
  x2_out = out;

  // surface y1
  vec2d.clear();
  vec2d = memBoundary[2];
  in = vec2d[0][0]->getCurrPos().gety();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPos().gety();
      if (tmp < out) out = tmp;
      if (tmp > in ) in  = tmp;
    }
  vx.push_back(in);
  vx.push_back(out);
  y1_in  = in;
  y1_out = out;

  // surface y2
  vec2d.clear();
  vec2d = memBoundary[3];
  in = vec2d[0][0]->getCurrPos().gety();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPos().gety();
      if (tmp > out) out = tmp;
      if (tmp < in ) in  = tmp;
    }
  vx.push_back(in);
  vx.push_back(out);
  y2_in  = in;
  y2_out = out;
  
  // surface z1
  vec2d.clear();
  vec2d = memBoundary[4];
  in = vec2d[0][0]->getCurrPos().getz();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPos().getz();
      if (tmp < out) out = tmp;
      if (tmp > in ) in  = tmp;
    }
  vx.push_back(in);
  vx.push_back(out);
  z1_in  = in;
  z1_out = out;

  // surface z2
  vec2d.clear();
  vec2d = memBoundary[5];
  in = vec2d[0][0]->getCurrPos().getz();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPos().getz();
      if (tmp > out) out = tmp;
      if (tmp < in ) in  = tmp;
    }
  vx.push_back(in);
  vx.push_back(out);
  z2_in  = in;
  z2_out = out;

}
  
//  1. it is important and helpful to mark a member function as const
//     if it does NOT change member data.
//  2. when a constant member function traverses member data, it can
//     NOT change the data.
//  3. then if it traverses a member data of a list, it should use a
//     const_iterator, otherwise compiler will give errors.
//  4. a const_iterator such as it also guarantees that (*it) will NOT
//     change any data. if (*it) call a modification function, the 
//     compiler will give errors.
void Assembly::printContact(const char* str) const
{
    std::ofstream ofs(str);
    if(!ofs) {
	cout << "stream error!" << endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    ofs << setw(OWID) << actualContactNum << endl;
    ofs << setw(OWID) << "ptcl_1"
        << setw(OWID) << "ptcl_2"
        << setw(OWID) << "point1_x"
        << setw(OWID) << "point1_y"
        << setw(OWID) << "point1_z"
        << setw(OWID) << "point2_x"
        << setw(OWID) << "point2_y"
        << setw(OWID) << "point2_z"
        << setw(OWID) << "radius_1"
        << setw(OWID) << "radius_2"
        << setw(OWID) << "penetration"
        << setw(OWID) << "tangt_dispmt"
        << setw(OWID) << "contact_radius"
        << setw(OWID) << "R0"
        << setw(OWID) << "E0"
        << setw(OWID) << "normal_force"
        << setw(OWID) << "tangt_force"
        << setw(OWID) << "contact_x"
        << setw(OWID) << "contact_y"
        << setw(OWID) << "contact_z"
        << setw(OWID) << "normal_x"
        << setw(OWID) << "normal_y"
        << setw(OWID) << "normal_z"
        << setw(OWID) << "tangt_x"
        << setw(OWID) << "tangt_y"
        << setw(OWID) << "tangt_z"
        << setw(OWID) << "vibra_t_step"
        << setw(OWID) << "impact_t_step"
        << endl;
   std::vector<CONTACT>::const_iterator it;
    for (it=contactVec.begin();it!=contactVec.end();++it)
	ofs << setw(OWID) << it->getP1()->getId()
	    << setw(OWID) << it->getP2()->getId()
	    << setw(OWID) << it->getPoint1().getx()
	    << setw(OWID) << it->getPoint1().gety()
	    << setw(OWID) << it->getPoint1().getz()
	    << setw(OWID) << it->getPoint2().getx()
	    << setw(OWID) << it->getPoint2().gety()
	    << setw(OWID) << it->getPoint2().getz()
	    << setw(OWID) << it->getRadius1()
	    << setw(OWID) << it->getRadius2()
	    << setw(OWID) << it->getPenetration()
	    << setw(OWID) << it->getTgtDisp()
	    << setw(OWID) << it->getContactRadius()
	    << setw(OWID) << it->getR0()
	    << setw(OWID) << it->getE0()
	    << setw(OWID) << it->getNormalForce()
	    << setw(OWID) << it->getTgtForce()
	    << setw(OWID) << ( it->getPoint1().getx()+it->getPoint2().getx() )/2
	    << setw(OWID) << ( it->getPoint1().gety()+it->getPoint2().gety() )/2
	    << setw(OWID) << ( it->getPoint1().getz()+it->getPoint2().getz() )/2
	    << setw(OWID) << it->NormalForceVec().getx()
	    << setw(OWID) << it->NormalForceVec().gety()
	    << setw(OWID) << it->NormalForceVec().getz()
	    << setw(OWID) << it->TgtForceVec().getx()
	    << setw(OWID) << it->TgtForceVec().gety()
	    << setw(OWID) << it->TgtForceVec().getz()
	    << setw(OWID) << it->getVibraTimeStep()
	    << setw(OWID) << it->getImpactTimeStep()
	    << endl;
    ofs.close();
}

//start of def OPENMP 
#ifdef OPENMP	

#if OPENMP_IMPL == 0
// OpenMP implementation 0: ts partitions, each thread handles a partition, max diff = n*n*(1-1/ts)/ts
// implementation is based on linked list, also works for vector but not efficient.
void Assembly::findContact() { 
  contactVec.clear();
  int possContact = 0;
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p1,NULL); 
#endif

  int tid;   // thread id
  int ts;    // number of threads
  int num;   // number of particles
  int tnum;  // number of particles per thread
  int i, j;
  Vec u, v;
  num = particleVec.size();
  ot = particleVec.begin();
  std::vector<Particle*>::iterator ot, it, pt;
  
#pragma omp parallel num_threads(NUM_THREADS) private(tid, ts, tnum, it, pt, i, j, u, v) shared(num) reduction(+: possContact)
  {
    tid = omp_get_thread_num();
    ts  = omp_get_num_threads();
    tnum = num / ts;  // divide itso ts partitions
    rnum = num % ts;  // remainder of the division
    it = ot;          // start particle of each thread
    
    // determine starting point and extend of each partition
    // this algorithm applies to both list and vector
    if (rnum == 0) {
      for (i = 0; i < tid * tnum; ++i)
	++it;         // starting point of each partition
    }
    else {
      if (tid < rnum) {
	tnum += 1;    // tnum changed
	for (i = 0; i < tid * tnum ; ++i)
	  ++it;
      }
      else {
	for (i = 0; i < rnum * (tnum + 1) + (tid - rnum) * tnum; ++ i)
	  ++it;
      }
    }
    
    // explore each partition
    for (j = 0 ; j < tnum; ++j, ++it) { 
      u=(*it)->getCurrPos();
      for (pt = it, ++pt; pt != particleVec.end(); ++pt){
	v=(*pt)->getCurrPos();
	if (   ( vfabs(v-u) < (*it)->getA() + (*pt)->getA())
	       && ( (*it)->getType() !=  1 || (*pt)->getType() != 1  )      // not both are fixed particles
	       && ( (*it)->getType() !=  5 || (*pt)->getType() != 5  )      // not both are free boundary particles
	       && ( (*it)->getType() != 10 || (*pt)->getType() != 10 )  ) { // not both are ghost particles
	  contact<Particle> tmpct(*it, *pt); // a local and temparory object
	  ++possContact;
	  if(tmpct.isOverlapped())
#pragma omp critical
	    contactVec.push_back(tmpct);    // containers use value semantics, so a "copy" is pushed back.
	}
      }
    }
  }
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf << setw(OWID) << "findContact=" << setw(OWID) << timediffsec(time_p1, time_p2); 
#endif
  possContactNum   = possContact;
  actualContactNum = contactVec.size();
} // end of OpenMP implementation 0

#elif OPENMP_IMPL == 1
// OpenMP implementation 1: ts partitions, each thread handles a partition, max diff = n*n*(1-1/ts)/ts
// implementation is based on vector index.
void Assembly::findContact() { 
  contactVec.clear();
  int possContact = 0;
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p1,NULL); 
#endif
  int tid;   // thread id
  int ts;    // number of threads
  int num;   // number of particles
  int start; // start particle index of each thread
  int end;   // last particle index of each thread
  int i, j;
  Vec u, v;
  num = particleVec.size();
  
#pragma omp parallel num_threads(NUM_THREADS) private(tid, ts, start, end, i, j, u, v) shared(num) reduction(+: possContact)
  {
    tid = omp_get_thread_num();
    ts  = omp_get_num_threads();
    start = tid * num / ts;
    end   = (tid + 1) * num / ts - 1;
    
    // explore each partition
    for (i = start; i <= end; ++i) { 
      u = particleVec[i]->getCurrPos();
      for (j = i + 1; j < num; ++j) {
	v = particleVec[j]->getCurrPos();
	if (   ( vfabs(v-u) < particleVec[i]->getA() + particleVec[j]->getA() )
	       && ( particleVec[i]->getType() !=  1 || particleVec[j]->getType() != 1  )      // not both are fixed particles
	       && ( particleVec[i]->getType() !=  5 || particleVec[j]->getType() != 5  )      // not both are free boundary particles
	       && ( particleVec[i]->getType() != 10 || particleVec[j]->getType() != 10 )  ) { // not both are ghost particles
	  contact<Particle> tmpct(particleVec[i], particleVec[j]); // a local and temparory object
	  ++possContact;
	  if(tmpct.isOverlapped())
#pragma omp critical
	    contactVec.push_back(tmpct);    // containers use value semantics, so a "copy" is pushed back.
	}
      }
    }
  }
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf <<  setw(OWID) << "findContact=" << setw(OWID) << timediffsec(time_p1, time_p2); 
#endif
  possContactNum   = possContact;  
  actualContactNum = contactVec.size();
} // end of OpenMP implementation 1

#elif OPENMP_IMPL == 2
// OpenMP implementation 2: no partitions, each thread leaps by ts until completed, max diff = n*(ts-1)/ts
void Assembly::findContact() { 
  contactVec.clear();
  int possContact = 0;
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p1,NULL); 
#endif
  int tid;   // thread id
  int ts;    // number of threads
  int num;   // number of particles
  int i, j;
  Vec u, v;
  num = particleVec.size();
  
#pragma omp parallel num_threads(NUM_THREADS) private(tid, ts, i, j, u, v) shared(num) reduction(+: possContact)
  {
    tid = omp_get_thread_num();
    ts  = omp_get_num_threads();
    
    // explore each partition
    for (i = tid; i < num; i += ts) { 
      u = particleVec[i]->getCurrPos();
      for (j = i + 1; j < num; ++j) {
	v = particleVec[j]->getCurrPos();
	if (   ( vfabs(v-u) < particleVec[i]->getA() + particleVec[j]->getA() )
	       && ( particleVec[i]->getType() !=  1 || particleVec[j]->getType() != 1  )      // not both are fixed particles
	       && ( particleVec[i]->getType() !=  5 || particleVec[j]->getType() != 5  )      // not both are free boundary particles
	       && ( particleVec[i]->getType() != 10 || particleVec[j]->getType() != 10 )  ) { // not both are ghost particles
	  contact<Particle> tmpct(particleVec[i], particleVec[j]); // a local and temparory object
	  ++possContact;
	  if(tmpct.isOverlapped())
#pragma omp critical
	    contactVec.push_back(tmpct);    // containers use value semantics, so a "copy" is pushed back.
	}
      }
    }
  }
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf << setw(OWID) << "findContact=" << setw(OWID) << timediffsec(time_p1, time_p2); 
#endif
  possContactNum   = possContact;  
  actualContactNum = contactVec.size();
} // end of OpenMP implementation 2

#elif OPENMP_IMPL == 3
// OpenMP implementation 3: no partitions, each thread leaps by ts until num/2 and handles two particles, max diff = 0
void Assembly::findContact() { 
  contactVec.clear();
  int possContact = 0;
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p1,NULL); 
#endif
  int tid;   // thread id
  int ts;    // number of threads
  int num;   // number of particles
  int i, j, k;
  Vec u, v;
  num = particleVec.size();
  
#pragma omp parallel num_threads(NUM_THREADS) private(tid, ts, i, j, k, u, v) shared(num) reduction(+: possContact)
  {
    tid = omp_get_thread_num();
    ts  = omp_get_num_threads();
    
    // explore each partition, works whether num is odd or even
    for (i = tid; i <= num / 2; i += ts) {
      int inc = num - 1 - 2*i;
      if (inc == 0) inc = 1; // avoid infinite loop when num is odd
      for (k = i; k <= num - 1 - i; k += inc ) {
	u = particleVec[k]->getCurrPos();
	for (j = k + 1; j < num; ++j) {
	  v = particleVec[j]->getCurrPos();
	  if (   ( vfabs(v-u) < particleVec[k]->getA() + particleVec[j]->getA() )
		 && ( particleVec[k]->getType() !=  1 || particleVec[j]->getType() != 1  )      // not both are fixed particles
		 && ( particleVec[k]->getType() !=  5 || particleVec[j]->getType() != 5  )      // not both are free boundary particles
		 && ( particleVec[k]->getType() != 10 || particleVec[j]->getType() != 10 )  ) { // not both are ghost particles
	    contact<Particle> tmpct(particleVec[k], particleVec[j]); // a local and temparory object
	    ++possContact;
	    if(tmpct.isOverlapped())
#pragma omp critical
	    contactVec.push_back(tmpct);    // containers use value semantics, so a "copy" is pushed back.
	  }
	}
      }
    }
  }
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf << setw(OWID) << "findContact=" << setw(OWID) << timediffsec(time_p1, time_p2); 
#endif
  possContactNum   = possContact;  
  actualContactNum = contactVec.size();
} // end of OpenMP implementation 3

#elif OPENMP_IMPL == 4
// OpenMP implementation 4: no partitions, parallel for, various loop scheduling: (static), (static,1), (dynamic), (dynamic,1)
void Assembly::findContact() { 
  contactVec.clear();
  int possContact = 0;
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p1,NULL); 
#endif

  int num;   // number of particles
  int i, j;
  Vec u, v;
  num = particleVec.size();
  
#pragma omp parallel for num_threads(NUM_THREADS) private(i, j, u, v) shared(num) reduction(+: possContact) schedule(dynamic)
  for (i = 0; i < num - 1; ++i) { 
    u = particleVec[i]->getCurrPos();
    for (j = i + 1; j < num; ++j) {
      v = particleVec[j]->getCurrPos();
      if (   ( vfabs(v-u) < particleVec[i]->getA() + particleVec[j]->getA() )
	     && ( particleVec[i]->getType() !=  1 || particleVec[j]->getType() != 1  )      // not both are fixed particles
	     && ( particleVec[i]->getType() !=  5 || particleVec[j]->getType() != 5  )      // not both are free boundary particles
	     && ( particleVec[i]->getType() != 10 || particleVec[j]->getType() != 10 )  ) { // not both are ghost particles
	contact<Particle> tmpct(particleVec[i], particleVec[j]); // a local and temparory object
	++possContact;
	if(tmpct.isOverlapped())
#pragma omp critical
	  contactVec.push_back(tmpct);    // containers use value semantics, so a "copy" is pushed back.
      }
    }
  }
  
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf << setw(OWID) << "findContact=" << setw(OWID) << timediffsec(time_p1, time_p2); 
#endif
  possContactNum   = possContact;  
  actualContactNum = contactVec.size();
} // end of OpenMP implementation 4

#endif

#else //else of def OPENMP, i.e., serial versions start here:

//start of ndef BINNING
#ifndef BINNING
void Assembly::findContact(){ // serial version, O(n x n), n is the number of particles.
    contactVec.clear();
    possContactNum = 0;

#ifdef TIME_PROFILE
    double time_r = 0; // time consumed in contact resolution, i.e., tmpct.isOverlapped()
    gettimeofday(&time_p1,NULL); 
#endif
    
    int num1 = particleVec.size();  // particles inside container
    int num2 = recvParticleVec.size(); // particles inside container (at front) + particles from neighboring blocks (at end)
    for (int i = 0; i < num1 - 1; ++i) {
      Vec u = particleVec[i]->getCurrPos();
      for (int j = i + 1; j < num2; ++j){
	Vec v = recvParticleVec[j]->getCurrPos();
	if (   ( vfabs(v - u) < particleVec[i]->getA() + recvParticleVec[j]->getA())
	    && ( particleVec[i]->getType() !=  1 || recvParticleVec[j]->getType() != 1  )      // not both are fixed particles
	    && ( particleVec[i]->getType() !=  5 || recvParticleVec[j]->getType() != 5  )      // not both are free boundary particles
	    && ( particleVec[i]->getType() != 10 || recvParticleVec[j]->getType() != 10 )  ) { // not both are ghost particles
	  Contact<Particle> tmpct(particleVec[i], recvParticleVec[j]); // a local and temparory object
	  ++possContactNum;
#ifdef TIME_PROFILE
	  gettimeofday(&time_r1,NULL); 
#endif
	  if(tmpct.isOverlapped())
	    contactVec.push_back(tmpct);    // containers use value semantics, so a "copy" is pushed back.
#ifdef TIME_PROFILE
	  gettimeofday(&time_r2,NULL); 
	  time_r += timediffsec(time_r1, time_r2);
#endif
	}
      }
    }	
    
#ifdef TIME_PROFILE
    gettimeofday(&time_p2,NULL);
    debugInf << setw(OWID) << "findContact=" << setw(OWID) << timediffsec(time_p1, time_p2) << setw(OWID) << "isOverlapped=" << setw(OWID) << time_r; 
#endif
    
    actualContactNum = contactVec.size();
}

//else of ndef BINNING
#else
void Assembly::findContact(){ // serial version, binning methods, cell slightly larger than maximum particle
  contactVec.clear();
  possContactNum = 0;
  
#ifdef TIME_PROFILE
  double time_r = 0;
  gettimeofday(&time_p1,NULL); 
#endif
  REAL maxDiameter = gradation.getPtclMaxRadius() * 2;
  int  nx = floor (container.getDimx() / maxDiameter);
  int  ny = floor (container.getDimy() / maxDiameter);
  int  nz = floor (container.getDimz() *1.5 / maxDiameter);
  REAL dx = container.getDimx() / nx;
  REAL dy = container.getDimx() / ny;
  REAL dz = container.getDimx() *1.5 / nz;
  Vec  minCorner= container.getMinCorner();
  REAL x0 = minCorner.getx();
  REAL y0 = minCorner.gety();
  REAL z0 = minCorner.getz();
  
  // 26 neighbors of each cell
  int neighbor[26][3];
  int count = 0;
  for (int i = -1; i < 2; ++i)
    for (int j = -1; j < 2; ++j)
      for (int k = -1; k < 2; ++k) {
	if (! (i == 0 && j == 0 && k==0 ) ) {
	  neighbor[count][0] = i;
	  neighbor[count][1] = j;
	  neighbor[count][2] = k;
	  ++count;
	}
      }
 
  // 4-dimensional array of cellVec
  typedef pair<bool, vector<Particle*> > cellT;
  vector< vector< vector < cellT > > > cellVec;
  cellVec.resize(nx);
  for (int i = 0; i < cellVec.size(); ++i) {
    cellVec[i].resize(ny);
    for (int j = 0; j < cellVec[i].size(); ++j)
      cellVec[i][j].resize(nz);
  }
  // mark each cell as not searched
  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      for (int k = 0; k < nz; ++k)
	cellVec[i][j][k].first = false; // has not ever been searched

  // find particles in each cell
  Vec center;
  REAL x1, x2, y1, y2, z1, z2;
  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      for (int k = 0; k < nz; ++k) {
	x1 = x0 + dx * i;
	x2 = x0 + dx * (i + 1);
	y1 = y0 + dy * j;
	y2 = y0 + dy * (j + 1);
	z1 = z0 + dz * k;
	z2 = z0 + dz * (k + 1);
	for (int pt = 0; pt < particleVec.size(); ++pt) {
	  center = particleVec[pt]->getCurrPos();
	  if (center.getx() >= x1 && center.getx() < x2 &&
	      center.gety() >= y1 && center.gety() < y2 &&
	      center.getz() >= z1 && center.getz() < z2)
	    cellVec[i][j][k].second.push_back( particleVec[pt] );
	}
      }
  
  // for each cell:
  Particle *it, *pt;
  Vec u, v;
  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      for (int k = 0; k < nz; ++k) {
	// for particles inside the cell	  
	for (int m = 0; m < cellVec[i][j][k].second.size(); ++m) {
	  it = cellVec[i][j][k].second[m];
	  u  = it->getCurrPos();
	  
	  // for particles inside the cell itself   
	  for (int n = m + 1; n < cellVec[i][j][k].second.size(); ++n) {
	    //cout <<  i << " " << j << " " << k << " " << "m n size=" << m << " " << n << " " <<  cellVec[i][j][k].size() << endl;
	    pt = cellVec[i][j][k].second[n];
	    v  = pt->getCurrPos();
	    if ( ( vfabs(u-v) < it->getA() + pt->getA() )  &&
		 ( it->getType() !=  1 || pt->getType() != 1 ) &&   // not both are fixed particles
		 ( it->getType() !=  5 || pt->getType() != 5 ) &&   // not both are free boundary particles
		 ( it->getType() != 10 || pt->getType() != 10)  ) { // not both are ghost particles
	      contact<Particle> tmpct(it, pt); // a local and temparory object
	      ++possContactNum;
#ifdef TIME_PROFILE
		gettimeofday(&time_r1,NULL); 
#endif
	      if(tmpct.isOverlapped())
		contactVec.push_back(tmpct);   // containers use value semantics, so a "copy" is pushed back.
#ifdef TIME_PROFILE
		gettimeofday(&time_r2,NULL); 
		time_r += timediffsec(time_r1, time_r2);
#endif
	    }
	  }
	  
	  // for 26 neighboring cells
	  for (int ncell = 0; ncell < 26; ++ncell ) {
	    int ci = i + neighbor[ncell][0];
	    int cj = j + neighbor[ncell][1];
	    int ck = k + neighbor[ncell][2];
	    if (ci > -1 && ci < nx && cj > -1 && cj < ny && ck > -1 && ck < nz && cellVec[ci][cj][ck].first == false ) {
	      //cout << "i j k m ncell ci cj ck size contacts= " << i << " " << j << " " << k << " " << m  << " " << ncell << " " << ci << " " << cj << " " << ck << " " << cellVec[ci][cj][ck].second.size() << " "  << contactVec.size() << endl;
	      vector<Particle*> vt = cellVec[ci][cj][ck].second;
	      for (int n = 0; n < vt.size(); ++n) {
		pt = vt[n];
		v  = pt->getCurrPos();
		if ( ( vfabs(u-v) < it->getA() + pt->getA() )  &&
		     ( it->getType() !=  1 || pt->getType() != 1 ) &&   // not both are fixed particles
		     ( it->getType() !=  5 || pt->getType() != 5 ) &&   // not both are free boundary particles
		     ( it->getType() != 10 || pt->getType() != 10)  ) { // not both are ghost particles
		  contact<Particle> tmpct(it, pt); // a local and temparory object
		  ++possContactNum;
#ifdef TIME_PROFILE
		gettimeofday(&time_r1,NULL); 
#endif
		  if(tmpct.isOverlapped())
		    contactVec.push_back(tmpct);   // containers use value semantics, so a "copy" is pushed back.
#ifdef TIME_PROFILE
		gettimeofday(&time_r2,NULL); 
		time_r += timediffsec(time_r1, time_r2);
#endif
		  
		}
	      }
	    }
	  }
	}
	cellVec[i][j][k].first = true; // searched, will not be searched again
	
      }
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf << setw(OWID) << "findContact=" << setw(OWID) << timediffsec(time_p1, time_p2) << setw(OWID) << "isOverlapped=" << setw(OWID) << time_r; 
#endif
  
  actualContactNum = contactVec.size();
}

//end of ndef BINNING
#endif

//end of def OPENMP 
#endif

REAL Assembly::getDensity() const{
    REAL dens=0;
    std::vector<Particle*>::const_iterator it;
    for(it=particleVec.begin();it!=particleVec.end();++it)
	dens+=(*it)->getMass();
    return dens/=bulkVolume;
}

REAL Assembly::getAvgPenetration() const{
    int totalcntct = contactVec.size();
    if (totalcntct==0)
	return 0;
    else {
	REAL pene=0;
	for (std::vector<CONTACT>::const_iterator it=contactVec.begin();it!=contactVec.end();++it)
	    pene += it->getPenetration(); 
	return pene/totalcntct;
    }
}

REAL Assembly::getVibraTimeStep() const {
    int totalcntct = contactVec.size();
    if (totalcntct == 0)
	return 0;
    else {
	std::vector<CONTACT>::const_iterator it=contactVec.begin();
        REAL minTimeStep = it->getVibraTimeStep();
	for (++it; it != contactVec.end(); ++it) {
	  REAL val = it->getVibraTimeStep(); 
	  minTimeStep =  val < minTimeStep ? val : minTimeStep;
	}
	return minTimeStep;
    }
}

REAL Assembly::getImpactTimeStep() const {
    int totalcntct = contactVec.size();
    if (totalcntct == 0)
	return 0;
    else {
	std::vector<CONTACT>::const_iterator it=contactVec.begin();
        REAL minTimeStep = it->getImpactTimeStep();
	for (++it; it != contactVec.end(); ++it) {
	  REAL val = it->getImpactTimeStep(); 
	  minTimeStep =  val < minTimeStep ? val : minTimeStep;
	}
	return minTimeStep;
    }
}

REAL Assembly::getAvgVelocity() const {
    REAL avgv=0;
    int count=0;
    std::vector<Particle*>::const_iterator it;
    for(it=particleVec.begin();it!=particleVec.end();++it)
	if ((*it)->getType()==0) {
	    avgv+=vfabs((*it)->getCurrVeloc());
	    count++;
	}
    return avgv/=count;
}

REAL Assembly::getAvgOmga() const {
    REAL avgv=0;
    int count=0;
    std::vector<Particle*>::const_iterator it;
    for(it=particleVec.begin();it!=particleVec.end();++it)
	if ((*it)->getType()==0){
	    avgv+=vfabs((*it)->getCurrOmga());
	    count++;
	}
    return avgv/=count;
}

REAL Assembly::getAvgForce() const {
    REAL avgv=0;
    int count=0;
    std::vector<Particle*>::const_iterator it;
    for(it=particleVec.begin();it!=particleVec.end();++it)
	if ((*it)->getType()==0){
	    avgv+=vfabs((*it)->getForce());
	    count++;
	}
    return avgv/count;
}

REAL Assembly::getAvgMoment() const {
    REAL avgv=0;
    int count=0;
    std::vector<Particle*>::const_iterator it;
    for(it=particleVec.begin();it!=particleVec.end();++it)
	if ((*it)->getType()==0){
	    avgv+=vfabs((*it)->getMoment());
	    count++;
	}
    return avgv/=count;
}

REAL Assembly::getParticleVolume() const {
    REAL avgv=0;
    std::vector<Particle*>::const_iterator it;
    for(it=particleVec.begin();it!=particleVec.end();++it)
	if ((*it)->getType()==0)
	    avgv+=(*it)->getVolume();
    return avgv;
}

Vec Assembly::getTopFreeParticlePosition() const {
    std::vector<Particle*>::const_iterator it,jt,kt;
    it=particleVec.begin();
    while (it!=particleVec.end() && (*it)->getType()!=0)   // find the 1st free particle
	++it;

    if (it==particleVec.end())    // no free particles
	return 0;

    jt=it; 
    kt=it;
    
    // two cases:
    // 1: 1st particle is not free
    // 2: 1st particle is free
    if (++kt!=particleVec.end()){ // case1: more than 2 particles; case 2: more than 1 particle
	for(++it;it!=particleVec.end();++it){
	    if ((*it)->getType()==0)
		if ((*it)->getCurrPos().getz() > (*jt)->getCurrPos().getz())
		    jt=it;
	}
	return (*jt)->getCurrPos();
    }
    else {
	if ((*it)->getType()==0)  // case1: only 2 particles, the 2nd one is free; case2: only 1 particle
	    return (*it)->getCurrPos();
	else
	    return 0;
    }

}



REAL Assembly::ellipPileForce() {
    REAL val=0;
    for(std::vector<Particle*>::iterator it=particleVec.begin();it!=particleVec.end();++it)
	if ((*it)->getType()==3) {
	    val = (*it)->getForce().getz();
	    break;
	}
    return val;
}

Vec Assembly::ellipPileDimn() {
    Vec val;
    for(std::vector<Particle*>::iterator it=particleVec.begin();it!=particleVec.end();++it)
	if ((*it)->getType()==3) {
	    val = Vec((*it)->getA(), (*it)->getB(), (*it)->getC());
	    break;
	}
    return val;
}

REAL Assembly::ellipPileTipZ() {
    REAL val=0;
    for(std::vector<Particle*>::iterator it=particleVec.begin();it!=particleVec.end();++it)
	if ((*it)->getType()==3) {
	    val = (*it)->getCurrPos().getz()-(*it)->getA();
	    break;
	}
    return val;
}

REAL Assembly::ellipPilePeneVol() {
    REAL val=0;
    if (getTopFreeParticlePosition().getz()-ellipPileTipZ()<=0)
	val=0;
    else{
	// low: a signed number as lower limit for volumetric integration
	REAL low=ellipPileTipZ() + ellipPileDimn().getx() - getTopFreeParticlePosition().getz(); 
	REAL lowint=low-pow(low,3)/3.0/pow(ellipPileDimn().getx(),2);
	val = PI * ellipPileDimn().gety() * ellipPileDimn().getz()
	      *(2.0/3*ellipPileDimn().getx()-lowint);
    }
    return val;
}

void Assembly::ellipPileUpdate(){
  for(std::vector<Particle*>::iterator it=particleVec.begin();it!=particleVec.end();++it){
    if ((*it)->getType()==3) {
      (*it)->setCurrVeloc(Vec(0, 0, -PILE_RATE));
      (*it)->setCurrPos( (*it)->getPrevPos() + (*it)->getCurrVeloc() * TIMESTEP);
    }
  }
}

REAL Assembly::getTransEnergy() const{
    REAL engy=0;
    std::vector<Particle*>::const_iterator it;
    for(it=particleVec.begin();it!=particleVec.end();++it){
	if ((*it)->getType()==0)
	    engy+=(*it)->getTransEnergy();
    }
    return engy;
}

REAL Assembly::getRotatEnergy() const{
    REAL engy=0;
    std::vector<Particle*>::const_iterator it;
    for(it=particleVec.begin();it!=particleVec.end();++it){
	if ((*it)->getType()==0)
	    engy+=(*it)->getRotatEnergy();
    }
    return engy;
}

REAL Assembly::getKinetEnergy() const{
    REAL engy=0;
    std::vector<Particle*>::const_iterator it;
    for(it=particleVec.begin();it!=particleVec.end();++it){
	if ((*it)->getType()==0)
	    engy+=(*it)->getKinetEnergy();
    }
    return engy;
}

REAL Assembly::getPotenEnergy(REAL ref) const{
    REAL engy=0;
    std::vector<Particle*>::const_iterator it;
    for(it=particleVec.begin();it!=particleVec.end();++it){
	if ((*it)->getType()==0)
	    engy+=(*it)->getPotenEnergy(ref);
    }
    return engy;
}



void Assembly::springForce() {
  for (vector<Spring*>::iterator it = springVec.begin(); it != springVec.end(); ++it)
    (*it)->applyForce();
}

void Assembly::readCavityBoundary(const char* str) {
  std::ifstream ifs(str);
  if(!ifs) {
    cout << "stream error!" << endl; exit(-1);
  }  

  Boundary<Particle>* rbptr;
  int type;
  cavityBoundaryVec.clear();
  int boundaryNum;
  ifs >> boundaryNum;
  for(int i = 0; i < boundaryNum; i++){
    ifs >> type;
    if(type == 1) // plane boundary
      rbptr = new plnBoundary<Particle>(ifs);
    cavityBoundaryVec.push_back(rbptr);
  }

  ifs.close();
}


void Assembly::printCavityBoundary(const char* str) const {
  std::ofstream ofs(str);
  if(!ofs) {
    cout << "stream error!" << endl; exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  
  ofs << setw(OWID) << cavityBoundaryVec.size() << endl;
  std::vector<BOUNDARY*>::const_iterator rt;
  for(rt = cavityBoundaryVec.begin(); rt != cavityBoundaryVec.end(); ++rt)
    (*rt)->display(ofs);
  ofs << endl;
  
  ofs.close();
}



void Assembly::findParticleOnCavity(){
  std::vector<BOUNDARY*>::iterator rt;
  for(rt = cavityBoundaryVec.begin(); rt != cavityBoundaryVec.end(); ++rt)
    (*rt)->findParticleOnBoundary(allParticleVec);
}



void Assembly::boundaryForce(REAL penetr[],int cntnum[]){
  std::vector<BOUNDARY*>::iterator rt;
  for(rt=boundaryVec.begin();rt!=boundaryVec.end();++rt){	
    (*rt)->boundaryForce(boundaryTgtMap);
    for (int i = 1; i <= 6; ++i) {
      if ((*rt)->getBdryID() == i){
	penetr[i] = (*rt)->getAvgPenetr();
	cntnum[i] = (*rt)->getCntnum();
	break;
      }
    }
  }
}

void Assembly::cavityBoundaryForce(){
  std::vector<BOUNDARY*>::iterator rt;
  for(rt = cavityBoundaryVec.begin(); rt != cavityBoundaryVec.end(); ++rt)
    (*rt)->boundaryForce(boundaryTgtMap);
}

Vec Assembly::getNormalForce(int bdry) const{
    std::vector<BOUNDARY*>::const_iterator it;
    for(it=boundaryVec.begin();it!=boundaryVec.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getNormalForce();
    }
    return 0;
}

Vec Assembly::getShearForce(int bdry) const{
    std::vector<BOUNDARY*>::const_iterator it;
    for(it=boundaryVec.begin();it!=boundaryVec.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getShearForce();
    }
    return 0;
}

REAL Assembly::getAvgNormal(int bdry) const{
    std::vector<BOUNDARY*>::const_iterator it;
    for(it=boundaryVec.begin();it!=boundaryVec.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getAvgNormal();
    }
    return 0;
}

Vec Assembly::getApt(int bdry) const{
    std::vector<BOUNDARY*>::const_iterator it;
    for(it=boundaryVec.begin();it!=boundaryVec.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getApt();
    }
    return 0;
}


Vec Assembly::getDirc(int bdry) const{
    std::vector<BOUNDARY*>::const_iterator it;
    for(it=boundaryVec.begin();it!=boundaryVec.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getDirc();
    }
    return 0;
}

REAL Assembly::getArea(int n) const{
    std::vector<BOUNDARY*>::const_iterator it;
    for(it=boundaryVec.begin();it!=boundaryVec.end();++it){
	if((*it)->getBdryID()==n)
	    return (*it)->area;
    }
    return 0;
}

void Assembly::setArea(int n, REAL a){
    std::vector<BOUNDARY*>::iterator it;
    for(it=boundaryVec.begin();it!=boundaryVec.end();++it){
	if((*it)->getBdryID()==n)
	    (*it)->area=a;
    }
}

REAL Assembly::getAvgPressure() const{
    std::vector<BOUNDARY*>::const_iterator rt;
    REAL avgpres=0;
    for(rt=boundaryVec.begin();rt!=boundaryVec.end();++rt)
	avgpres+=vfabs((*rt)->getNormalForce())/(*rt)->getArea();
    return avgpres/=boundaryVec.size();
}

// only update CoefOfLimits[0] for specified boundaries
void Assembly::updateBoundary(int bn[], UPDATECTL rbctl[], int num){
    for(int i=0;i<num;i++){
	for(std::vector<BOUNDARY*>::iterator rt=boundaryVec.begin();rt!=boundaryVec.end();++rt){
	    if((*rt)->getBdryID()==bn[i]){
		(*rt)->update(rbctl[i]);
		break;
	    }
	}
    }
}

// update CoefOfLimits[1,2,3,4] for all 6 boundaries
void Assembly::updateBoundary6(){
    for(std::vector<BOUNDARY*>::iterator rt=boundaryVec.begin();rt!=boundaryVec.end();++rt){
	if((*rt)->getBdryID()==1 || (*rt)->getBdryID()==3){
	    for(std::vector<BOUNDARY*>::iterator lt=boundaryVec.begin();lt!=boundaryVec.end();++lt){
		if((*lt)->getBdryID()==4)
		    (*rt)->CoefOfLimits[1].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==2)
		    (*rt)->CoefOfLimits[2].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==5)
		    (*rt)->CoefOfLimits[3].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==6)
		    (*rt)->CoefOfLimits[4].apt=(*lt)->CoefOfLimits[0].apt;
	    }
	}
	else if((*rt)->getBdryID()==2 || (*rt)->getBdryID()==4){
	    for(std::vector<BOUNDARY*>::iterator lt=boundaryVec.begin();lt!=boundaryVec.end();++lt){
		if((*lt)->getBdryID()==1)
		    (*rt)->CoefOfLimits[1].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==3)
		    (*rt)->CoefOfLimits[2].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==5)
		    (*rt)->CoefOfLimits[3].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==6)
		    (*rt)->CoefOfLimits[4].apt=(*lt)->CoefOfLimits[0].apt;
	    }

	}
	else if((*rt)->getBdryID()==5 || (*rt)->getBdryID()==6){
	    for(std::vector<BOUNDARY*>::iterator lt=boundaryVec.begin();lt!=boundaryVec.end();++lt){
		if((*lt)->getBdryID()==1)
		    (*rt)->CoefOfLimits[1].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==3)
		    (*rt)->CoefOfLimits[2].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==2)
		    (*rt)->CoefOfLimits[3].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==4)
		    (*rt)->CoefOfLimits[4].apt=(*lt)->CoefOfLimits[0].apt;
	    }

	}
	
    }
}

void Assembly::deposit_repose(int   interval,
			      const char* inibdryfile,
			      const char* ParticleFile, 
			      const char* contactfile,
			      const char* progressfile, 
			      const char* debugfile)
{
  this->container = container;
  
  buildBoundary(5, inibdryfile); // container unchanged
  
  angleOfRepose(interval,           // print interval
		inibdryfile,        // input file, initial boundaries
		ParticleFile,       // output file, resulted particles, including snapNum 
		contactfile,        // output file, resulted contacts, including snapNum 
		progressfile,       // output file, statistical info
		debugfile);         // output file, debug info
}

void Assembly::angleOfRepose(int   interval,
			     const char* inibdryfile,
			     const char* ParticleFile, 
			     const char* contactfile,
			     const char* progressfile, 
			     const char* debugfile)
{
  // pre_1: open streams for output.
  progressinf.open(progressfile); 
  if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
  progressinf.setf(std::ios::scientific, std::ios::floatfield);
  progressinf.precision(OPREC);
  progressinf << setw(OWID) << "iteration"
	      << setw(OWID) << "poss_contact"
	      << setw(OWID) << "actual_contact"
	      << setw(OWID) << "penetration"
	      << setw(OWID) << "avg_normal"
	      << setw(OWID) << "avg_tangt"
	      << setw(OWID) << "avg_velocity"
	      << setw(OWID) << "avg_omga"
	      << setw(OWID) << "avg_force"
	      << setw(OWID) << "avg_moment"
	      << setw(OWID) << "trans_energy"
	      << setw(OWID) << "rotat_energy"
	      << setw(OWID) << "kinet_energy"
	      << setw(OWID) << "poten_energy"
	      << setw(OWID) << "total_energy"
	      << setw(OWID) << "void_ratio"
	      << setw(OWID) << "porosity"
	      << setw(OWID) << "coord_number"
	      << setw(OWID) << "density"
	      << setw(OWID) << "sigma_y1"
	      << setw(OWID) << "sigma_y2"
	      << setw(OWID) << "sigma_x1"
	      << setw(OWID) << "sigma_x2"
	      << setw(OWID) << "sigma_z1"
	      << setw(OWID) << "sigma_z2"
	      << setw(OWID) << "mean_stress"
	      << setw(OWID) << "dimx"
	      << setw(OWID) << "dimy"
	      << setw(OWID) << "dimz"
	      << setw(OWID) << "volume"
	      << setw(OWID) << "epsilon_x"
	      << setw(OWID) << "epsilon_y"
	      << setw(OWID) << "epsilon_z"
	      << setw(OWID) << "epsilon_v"
	      << setw(OWID) << "vibra_t_step"
	      << setw(OWID) << "impact_t_step"
	      << setw(OWID) << "wall_time" << endl;
  
  debugInf.open(debugfile);
  if(!debugInf) { cout << "stream error!" << endl; exit(-1); }
  debugInf.setf(std::ios::scientific, std::ios::floatfield);
  
  // pre_2. create boundaries from existing files.
  readBoundary(inibdryfile);

  // pre_3: define variables used in iterations.
  REAL avgNormal=0;
  REAL avgTangt=0;
  int  stepsnum=0;
  char stepsstr[4];
  bool toSnapshot=false;
  char stepsfp[50];
  REAL void_ratio=0;
  REAL bdry_penetr[7] = {0,0,0,0,0,0,0};
  int  bdry_cntnum[7] = {0,0,0,0,0,0,0};

  REAL maxRadius = gradation.getPtclMaxRadius();
  REAL maxDiameter = maxRadius * 2.0;
  REAL z0 = container.getMinCorner().getz();
  vector<Particle*> lastPtcls;
  Particle *newPtcl = NULL;
  int layers = 1; // how many layers of new particles to generate each time

  iteration = 0; 
  int particleNum = 0;
  REAL zCurr;
  gettimeofday(&time_w1,NULL);
  // iterations starting ...
  do
    {
      // 1. add particle
      if ( particleNum == 0 ) {
	zCurr = z0 + maxRadius;

	for ( int i = 0; i != layers; ++i) {
	  newPtcl = new Particle(particleNum+1, 0, Vec(0, 0, zCurr + maxDiameter * i), gradation, YOUNG, POISSON);
	  particleVec.push_back(newPtcl);
	  ++particleNum;
	  lastPtcls.push_back(newPtcl);
	  
	  newPtcl = new Particle(particleNum+1, 0, Vec(maxDiameter, 0, zCurr + maxDiameter * i), gradation, YOUNG, POISSON);
	  particleVec.push_back(newPtcl);
	  ++particleNum;
	  //lastPtcls.push_back(newPtcl);
	  
	  newPtcl = new Particle(particleNum+1, 0, Vec(-maxDiameter, 0, zCurr + maxDiameter * i), gradation, YOUNG, POISSON);
	  particleVec.push_back(newPtcl);
	  ++particleNum;
	  //lastPtcls.push_back(newPtcl);
	  
	  newPtcl = new Particle(particleNum+1, 0, Vec(0, maxDiameter, zCurr + maxDiameter * i), gradation, YOUNG, POISSON);
	  particleVec.push_back(newPtcl);
	  ++particleNum;
	  //lastPtcls.push_back(newPtcl);
	  
	  newPtcl = new Particle(particleNum+1, 0, Vec(0, -maxDiameter, zCurr + maxDiameter * i), gradation, YOUNG, POISSON);
	  particleVec.push_back(newPtcl);
	  ++particleNum;
	  //lastPtcls.push_back(newPtcl);
	}
	toSnapshot = true;

      }
      else {
	vector<Particle*>::iterator it;
	bool allInContact = false;
	for ( it = lastPtcls.begin(); it != lastPtcls.end(); ++it) {
	  if ( (*it)->isInContact() ) 
	    allInContact = true;
	  else {
	    allInContact = false;
	    break;
	  }
	}

	if ( allInContact ) {

	  lastPtcls.clear(); // do not delete those pointers to release memory; particleVec will do it.
	  zCurr = getPtclMaxZ() + maxDiameter;

	  for ( int i = 0; i != layers; ++i) {
	    newPtcl = new Particle(particleNum+1, 0, Vec(0, 0, zCurr + maxDiameter * i), gradation, YOUNG, POISSON);
	    particleVec.push_back(newPtcl);
	    ++particleNum;
	    lastPtcls.push_back(newPtcl);
	    
	    newPtcl = new Particle(particleNum+1, 0, Vec(maxDiameter, 0, zCurr + maxDiameter * i), gradation, YOUNG, POISSON);
	    particleVec.push_back(newPtcl);
	    ++particleNum;
	    //lastPtcls.push_back(newPtcl);
	    
	    newPtcl = new Particle(particleNum+1, 0, Vec(-maxDiameter, 0, zCurr + maxDiameter * i), gradation, YOUNG, POISSON);
	    particleVec.push_back(newPtcl);
	    ++particleNum;
	    //lastPtcls.push_back(newPtcl);
	    
	    newPtcl = new Particle(particleNum+1, 0, Vec(0, maxDiameter, zCurr + maxDiameter * i), gradation, YOUNG, POISSON);
	    particleVec.push_back(newPtcl);
	    ++particleNum;
	    //lastPtcls.push_back(newPtcl);
	    
	    newPtcl = new Particle(particleNum+1, 0, Vec(0, -maxDiameter, zCurr + maxDiameter * i), gradation, YOUNG, POISSON);
	    particleVec.push_back(newPtcl);
	    ++particleNum;
	    //lastPtcls.push_back(newPtcl);	
	  }
	  toSnapshot = true;
	}
      }
      // 2. create possible boundary particles and contacts between particles.
      findContact();
      findParticleOnBoundary();
      
      // 3. set particle forces/moments as zero before each re-calculation,
      clearForce();	
      
      // 4. calculate contact forces/moments and apply them to particles.
      internalForce(avgNormal, avgTangt);
      
      // 5. calculate boundary forces/moments and apply them to particles.
      boundaryForce(bdry_penetr, bdry_cntnum);
      
      // 6. update particles' velocity/omga/position/orientation based on force/moment.
      updateParticle();
      
      // 7. (1) output particles and contacts information as snapNum.
      if (toSnapshot){
	sprintf(stepsstr, "%03d", stepsnum); 
	strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	printParticle(stepsfp);
	
	sprintf(stepsstr, "%03d", stepsnum); 
	strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	printContact(stepsfp);
	time(&timeStamp);
	g_timeinf << setw(4) << stepsnum << " " << ctime(&timeStamp) << flush;
	++stepsnum;
	toSnapshot = false;
      }
      
      // 8. (2) output stress and strain info.
      if (iteration % interval == 0) {
	gettimeofday(&time_w2,NULL);
	REAL t1=getTransEnergy();
	REAL t2=getRotatEnergy();
	REAL t3=getPotenEnergy(-0.025);
	progressinf << setw(OWID) << iteration
		    << setw(OWID) << getPossContactNum()
		    << setw(OWID) << getActualContactNum()
		    << setw(OWID) << getAvgPenetration()
		    << setw(OWID) << avgNormal
		    << setw(OWID) << avgTangt
		    << setw(OWID) << getAvgVelocity() 
		    << setw(OWID) << getAvgOmga()
		    << setw(OWID) << getAvgForce()   
		    << setw(OWID) << getAvgMoment()
		    << setw(OWID) << t1
		    << setw(OWID) << t2
		    << setw(OWID) << (t1+t2)
		    << setw(OWID) << t3
		    << setw(OWID) << (t1+t2+t3)
		    << setw(OWID) << void_ratio
		    << setw(OWID) << void_ratio/(1+void_ratio)
		    << setw(OWID) << 2.0*(getActualContactNum()
				      +bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
				      +bdry_cntnum[4]+bdry_cntnum[6])/allParticleVec.size()
		    << setw(OWID) << "0"
		    << setw(OWID) << "0"
		    << setw(OWID) << "0"
		    << setw(OWID) << "0"
		    << setw(OWID) << "0"
		    << setw(OWID) << "0"
		    << setw(OWID) << "0"
		    << setw(OWID) << "0"
		    << setw(OWID) << "0"
		    << setw(OWID) << "0"
		    << setw(OWID) << "0"
		    << setw(OWID) << "0"
		    << setw(OWID) << "0"
		    << setw(OWID) << "0"
		    << setw(OWID) << "0"
		    << setw(OWID) << "0"
		    << setw(OWID) << getVibraTimeStep()
		    << setw(OWID) << getImpactTimeStep()
		    << setw(OWID) << timediffsec(time_w1,time_w2)
		        << endl;
      }
      
      // 7. loop break conditions.
      ++iteration;
      
    } while (particleNum < 2000); //( zCurr < container.getMaxCorner().getz() );  //(++iteration < totalSteps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    g_timeinf << setw(4) << "end" << " " << ctime(&timeStamp) << flush;

    // post_2. close streams
    progressinf.close();
    debugInf.close();
}







void Assembly::scale_PtclBdry(int   totalSteps,  
			      int   snapNum,
			      int   interval,
			      REAL dimn,
			      REAL rsize,
			      const char* iniptclfile,   
			      const char* ParticleFile, 
			      const char* contactfile,
			      const char* progressfile, 
			      const char* debugfile)
{
    deposit_p(totalSteps,        // totalSteps
	      snapNum,          // number of snapNum
	      interval,           // print interval
	      dimn,               // dimension of particle-composed-boundary
	      rsize,              // relative container size
	      iniptclfile,        // input file, initial particles
	      ParticleFile,       // output file, resulted particles, including snapNum 
	      contactfile,        // output file, resulted contacts, including snapNum 
	      progressfile,       // output file, statistical info
	      debugfile);         // output file, debug info
}


// collapse a deposited specimen through gravitation
void Assembly::collapse(int   totalSteps,  
			int   snapNum,
			int   interval,
			const char* iniptclfile,
			const char* initboundary,
			const char* ParticleFile,
			const char* contactfile,
			const char* progressfile,
			const char* debugfile)
{
  buildBoundary(1,              // 1-only bottom boundary; 5-no top boundary;6-boxed 6 boundaries
		initboundary);  // output file, containing boundaries info
  
  deposit(totalSteps,        // number of iterations
	  snapNum,          // number of snapNum
	  interval,           // print interval
	  iniptclfile,        // input file, initial particles
	  initboundary,       // input file, boundaries
	  ParticleFile,       // output file, resulted particles, including snapNum 
	  contactfile,        // output file, resulted contacts, including snapNum 
	  progressfile,       // output file, statistical info
	  debugfile);         // output file, debug info
}

  


// make a cavity inside the sample and remove particles in the cavity
void Assembly::trimCavity(bool toRebuild,
			  const char* ParticleFile,
			  const char* cavParticleFile)
{
  if (toRebuild) readParticle(ParticleFile);
  trimHistoryNum = allParticleVec.size();

  REAL x1,x2,y1,y2,z1,z2,x0,y0,z0;
  x1 = cavity.getMinCorner().getx();
  y1 = cavity.getMinCorner().gety();
  z1 = cavity.getMinCorner().getz();
  x2 = cavity.getMaxCorner().getx();
  y2 = cavity.getMaxCorner().gety();
  z2 = cavity.getMaxCorner().getz();
  x0 = cavity.getCenter().getx();
  y0 = cavity.getCenter().gety();
  z0 = cavity.getCenter().getz();
 
  std::vector<Particle*>::iterator itr;
  Vec center;
  REAL delta = gradation.getPtclMaxRadius();

  for (itr = particleVec.begin(); itr != particleVec.end(); ){
    center=(*itr)->getCurrPos();
    if(center.getx() + delta  >= x1 && center.getx() - delta <= x2 &&
       center.gety() + delta  >= y1 && center.gety() - delta <= y2 &&
       center.getz() + delta  >= z1 && center.getz() - delta <= z2 )
      {
	delete (*itr); // release memory
	itr = particleVec.erase(itr); 
      }
    else
      ++itr;
  }
  
  printParticle(cavParticleFile);
}


// expand partcile size by some percentage for particles inside cavity
void Assembly::expandCavityParticles(bool toRebuild,
				     REAL percent,
				     const char* cavityptclfile,
				     const char* ParticleFile,
				     const char* newptclfile)
{
  if (toRebuild) readParticle(ParticleFile);
  trimHistoryNum = allParticleVec.size();

  REAL x1,x2,y1,y2,z1,z2;
  x1 = cavity.getMinCorner().getx();
  y1 = cavity.getMinCorner().gety();
  z1 = cavity.getMinCorner().getz();
  x2 = cavity.getMaxCorner().getx();
  y2 = cavity.getMaxCorner().gety();
  z2 = cavity.getMaxCorner().getz();
 
  std::vector<Particle*>::iterator itr;
  Vec center;

  int cavityPtclNum = 0;
  for (itr = particleVec.begin(); itr != particleVec.end(); ++itr ){
    center=(*itr)->getCurrPos();
    if(center.getx() > x1 && center.getx() < x2 &&
       center.gety() > y1 && center.gety() < y2 &&
       center.getz() > z1 && center.getz() < z2 )
      ++cavityPtclNum;
  }

  printCavityParticle(cavityPtclNum, cavityptclfile);

  for (itr = particleVec.begin(); itr != particleVec.end(); ++itr ){
    center=(*itr)->getCurrPos();
    if(center.getx() > x1 && center.getx() < x2 &&
       center.gety() > y1 && center.gety() < y2 &&
       center.getz() > z1 && center.getz() < z2 )
      (*itr)->expand(percent);
  }

  printParticle(newptclfile);
}


void Assembly::printCavityParticle(int total, const char* str) const {
  std::ofstream ofs(str);
  if(!ofs) {
    cout << "stream error!" << endl; exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  ofs << setw(OWID) << total << setw(OWID) << 1 << endl;
  ofs << setw(OWID) << cavity.getCenter().getx()
      << setw(OWID) << cavity.getCenter().gety()
      << setw(OWID) << cavity.getCenter().getz()
      << setw(OWID) << cavity.getDimx()
      << setw(OWID) << cavity.getDimy()
      << setw(OWID) << cavity.getDimz() << endl;
  
  ofs << setw(OWID) << "ID"
      << setw(OWID) << "type"
      << setw(OWID) << "radius_a"
      << setw(OWID) << "radius_b"
      << setw(OWID) << "radius_c"
      << setw(OWID) << "position_x"
      << setw(OWID) << "position_y"
      << setw(OWID) << "position_z"
      << setw(OWID) << "axle_a_x"
      << setw(OWID) << "axle_a_y"
      << setw(OWID) << "axle_a_z"
      << setw(OWID) << "axle_b_x"
      << setw(OWID) << "axle_b_y"
      << setw(OWID) << "axle_b_z"
      << setw(OWID) << "axle_c_x"
      << setw(OWID) << "axle_c_y"
      << setw(OWID) << "axle_c_z"
      << setw(OWID) << "velocity_x"
      << setw(OWID) << "velocity_y"
      << setw(OWID) << "velocity_z"
      << setw(OWID) << "omga_x"
      << setw(OWID) << "omga_y"
      << setw(OWID) << "omga_z"
      << setw(OWID) << "force_x"
      << setw(OWID) << "force_y"
      << setw(OWID) << "force_z"
      << setw(OWID) << "moment_x"
      << setw(OWID) << "moment_y"
      << setw(OWID) << "moment_z"
      << endl;

  REAL x1,x2,y1,y2,z1,z2;
  x1 = cavity.getMinCorner().getx();
  y1 = cavity.getMinCorner().gety();
  z1 = cavity.getMinCorner().getz();
  x2 = cavity.getMaxCorner().getx();
  y2 = cavity.getMaxCorner().gety();
  z2 = cavity.getMaxCorner().getz();
  
  Vec tmp;
  std::vector<Particle*>::const_iterator  it;
  for (it=particleVec.begin();it!=particleVec.end();++it)  {
    Vec center=(*it)->getCurrPos();
    if(center.getx() > x1 && center.getx() < x2 &&
       center.gety() > y1 && center.gety() < y2 &&
       center.getz() > z1 && center.getz() < z2 ) {

    ofs << setw(OWID) << (*it)->getId()
	<< setw(OWID) << (*it)->getType()
	<< setw(OWID) << (*it)->getA()
	<< setw(OWID) << (*it)->getB()
	<< setw(OWID) << (*it)->getC();
    
    tmp=(*it)->getCurrPos();
    ofs << setw(OWID) << tmp.getx()
	<< setw(OWID) << tmp.gety()
	<< setw(OWID) << tmp.getz();
    
    tmp=(*it)->getCurrDirecA();
    ofs << setw(OWID) << tmp.getx()
	<< setw(OWID) << tmp.gety()
	<< setw(OWID) << tmp.getz();
    
    tmp=(*it)->getCurrDirecB();
    ofs << setw(OWID) << tmp.getx()
	<< setw(OWID) << tmp.gety()
	<< setw(OWID) << tmp.getz();
    
    tmp=(*it)->getCurrDirecC();
    ofs << setw(OWID) << tmp.getx()
	<< setw(OWID) << tmp.gety()
	<< setw(OWID) << tmp.getz();
    
    tmp=(*it)->getCurrVeloc();
    ofs << setw(OWID) << tmp.getx()
	<< setw(OWID) << tmp.gety()
	<< setw(OWID) << tmp.getz();
    
    tmp=(*it)->getCurrOmga();
    ofs << setw(OWID) << tmp.getx()
	<< setw(OWID) << tmp.gety()
	<< setw(OWID) << tmp.getz();
    
    tmp=(*it)->getForce();
    ofs << setw(OWID) << tmp.getx()
	<< setw(OWID) << tmp.gety()
	<< setw(OWID) << tmp.getz();
    
    tmp=(*it)->getMoment();
    ofs << setw(OWID) << tmp.getx()
	<< setw(OWID) << tmp.gety()
	<< setw(OWID) << tmp.getz() << endl;
    }
  }
  
  ofs.close();
}


// bdrymum = 6 by default
// the variable existMaxID is important because cavity and container
// use the same boundaryTgtMap.
void Assembly::buildCavityBoundary(int existMaxId, const char* boundaryFile)
{
  std::ofstream ofs(boundaryFile);
  if(!ofs) { cout << "stream error!" << endl; exit(-1);}

  REAL x1,x2,y1,y2,z1,z2,x0,y0,z0;
  x1 = cavity.getMinCorner().getx();
  y1 = cavity.getMinCorner().gety();
  z1 = cavity.getMinCorner().getz();
  x2 = cavity.getMaxCorner().getx();
  y2 = cavity.getMaxCorner().gety();
  z2 = cavity.getMaxCorner().getz();
  x0 = cavity.getCenter().getx();
  y0 = cavity.getCenter().gety();
  z0 = cavity.getCenter().getz();

  int boundaryNum = 6;

  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs << setw(OWID) << 0
      << setw(OWID) << boundaryNum << endl << endl;

  // boundary 1
  ofs << setw(OWID) << 1 << endl
      << setw(OWID) << existMaxId + 1
      << setw(OWID) << 5
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << -1
      << setw(OWID) << 0
      << setw(OWID) << 0     
      << setw(OWID) << x2
      << setw(OWID) << y0
      << setw(OWID) << z0     
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << -1
      << setw(OWID) << 0
      << setw(OWID) << x0    
      << setw(OWID) << y1
      << setw(OWID) << z0     
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0 
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << x0      
      << setw(OWID) << y2
      << setw(OWID) << z0     
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << 1
      << setw(OWID) << x0     
      << setw(OWID) << y0    
      << setw(OWID) << z2 
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0 
      << setw(OWID) << -1
      << setw(OWID) << x0     
      << setw(OWID) << y0     
      << setw(OWID) << z1
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl << endl
    
    // boundary 2
      << setw(OWID) << 1 << endl
      << setw(OWID) << existMaxId + 2
      << setw(OWID) << 5
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << -1
      << setw(OWID) << 0     
      << setw(OWID) << x0    
      << setw(OWID) << y2
      << setw(OWID) << z0     
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 1
      << setw(OWID) << 0 
      << setw(OWID) << 0
      << setw(OWID) << x2
      << setw(OWID) << y0
      << setw(OWID) << z0     
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << -1
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << x1 
      << setw(OWID) << y0
      << setw(OWID) << z0     
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << 1
      << setw(OWID) << x0     
      << setw(OWID) << y0    
      << setw(OWID) << z2 
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0 
      << setw(OWID) << -1
      << setw(OWID) << x0     
      << setw(OWID) << y0      
      << setw(OWID) << z1
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl << endl
    
    // boundary 3
      << setw(OWID) << 1 << endl
      << setw(OWID) << existMaxId + 3
      << setw(OWID) << 5
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0     
      << setw(OWID) << x1
      << setw(OWID) << y0
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << -1
      << setw(OWID) << 0
      << setw(OWID) << x0     
      << setw(OWID) << y1
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0  
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << x0       
      << setw(OWID) << y2
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << 1
      << setw(OWID) << x0      
      << setw(OWID) << y0     
      << setw(OWID) << z2 
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0 
      << setw(OWID) << -1
      << setw(OWID) << x0      
      << setw(OWID) << y0      
      << setw(OWID) << z1
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl << endl
    
    // boundary 4
      << setw(OWID) << 1 << endl
      << setw(OWID) << existMaxId + 4
      << setw(OWID) << 5
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0 
      << setw(OWID) << 1
      << setw(OWID) << 0     
      << setw(OWID) << x0      
      << setw(OWID) << y1
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 1
      << setw(OWID) << 0 
      << setw(OWID) << 0
      << setw(OWID) << x2
      << setw(OWID) << y0
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << -1 
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << x1 
      << setw(OWID) << y0
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << 1
      << setw(OWID) << x0      
      << setw(OWID) << y0     
      << setw(OWID) << z2 
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0 
      << setw(OWID) << -1
      << setw(OWID) << x0      
      << setw(OWID) << y0      
      << setw(OWID) << z1
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl << endl
    
    // boundary 5
      << setw(OWID) << 1 << endl
      << setw(OWID) << existMaxId + 5
      << setw(OWID) << 5
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0 
      << setw(OWID) << 0
      << setw(OWID) << -1     
      << setw(OWID) << x0      
      << setw(OWID) << y0
      << setw(OWID) << z2 
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 1
      << setw(OWID) << 0 
      << setw(OWID) << 0
      << setw(OWID) << x2
      << setw(OWID) << y0
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << -1 
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << x1 
      << setw(OWID) << y0
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << x0      
      << setw(OWID) << y2
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << -1
      << setw(OWID) << 0
      << setw(OWID) << x0      
      << setw(OWID) << y1
      << setw(OWID) << z0
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl << endl
    
    // boundary 6
      << setw(OWID) << 1 << endl
      << setw(OWID) << existMaxId + 6
      << setw(OWID) << 5
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0 
      << setw(OWID) << 0
      << setw(OWID) << 1    
      << setw(OWID) << x0      
      << setw(OWID) << y0
      << setw(OWID) << z1
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 1
      << setw(OWID) << 0 
      << setw(OWID) << 0
      << setw(OWID) << x2
      << setw(OWID) << y0
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << -1 
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << x1 
      << setw(OWID) << y0
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << x0      
      << setw(OWID) << y2
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << -1
      << setw(OWID) << 0 
      << setw(OWID) << x0      
      << setw(OWID) << y1
      << setw(OWID) << z0
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl << endl; 

  ofs.close();
}


// create boundary particles and springs connecting those boundary particles
void Assembly::createMemParticle(REAL rRadius,
				 bool toRebuild,
				 const char* ParticleFile,
				 const char* allParticle)
{
  if (toRebuild) readParticle(ParticleFile);

  REAL radius = gradation.getMinPtclRadius();
  if (gradation.getSize().size() == 1 &&
      gradation.getPtclRatioBA() == 1.0 && 
      gradation.getPtclRatioCA() == 1.0)
    radius *= rRadius; // determine how tiny the boundary particles are
  REAL diameter = radius*2;
  Vec v1 = allContainer.getMinCorner();
  Vec v2 = allContainer.getMaxCorner();
  Vec v0 = allContainer.getCenter();
  REAL x1 = v1.getx();
  REAL y1 = v1.gety();
  REAL z1 = v1.getz();
  REAL x2 = v2.getx();
  REAL y2 = v2.gety();
  REAL z2 = v2.getz();
  REAL x0 = v0.getx();
  REAL y0 = v0.gety();
  REAL z0 = v0.getz();

  Particle* newptcl = NULL;
  REAL x, y, z;
  
  vector<Particle*> vec1d;  // 1-dimension
  vector< vector<Particle*>  > vec2d; // 2-dimension
  Spring* newSpring = NULL;
  REAL young   = YOUNG;   //memYOUNG;
  REAL poisson = POISSON; //memPOISSON;
  REAL modulus = memYOUNG;

  int memPtclIndex = trimHistoryNum;
  // process in the order of surfaces: x1 x2 y1 y2 z1 z2
  // surface x1
  x = x1 - radius;
  for (z = z1 + radius; z <= z2 - radius + EPS; z += diameter) {
    vec1d.clear();
    for (y = y1 + radius; y <= y2 - radius + EPS; y += diameter) {
      newptcl = new Particle(++memPtclIndex, 5, Vec(x,y,z), radius, young, poisson);
      vec1d.push_back(newptcl);
      particleVec.push_back(newptcl);
    }
    vec2d.push_back(vec1d);
  }
  memBoundary.push_back(vec2d);
  for (int i = 0; i != vec2d.size() ; ++i)
    for (int j = 0; j != vec2d[i].size() ; ++ j) {
      if (j + 1 < vec2d[i].size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], modulus);
	springVec.push_back(newSpring);
      }
      if (i + 1 < vec2d.size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], modulus);
	springVec.push_back(newSpring);
      }
    }

  // surface x2     
  vec2d.clear();
  x = x2 + radius;
  for (z = z1 + radius; z <= z2 - radius + EPS; z += diameter) {
    vec1d.clear();
    for (y = y1 + radius; y <= y2 - radius + EPS; y += diameter) {
      newptcl = new Particle(++memPtclIndex, 5, Vec(x,y,z), radius, young, poisson);
      vec1d.push_back(newptcl);
      particleVec.push_back(newptcl);
    }
    vec2d.push_back(vec1d);
  }
  memBoundary.push_back(vec2d);
  for (int i = 0; i != vec2d.size() ; ++i)
    for (int j = 0; j != vec2d[i].size() ; ++ j) {
      if (j + 1 < vec2d[i].size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], modulus);
	springVec.push_back(newSpring);
      }
      if (i + 1 < vec2d.size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], modulus);
	springVec.push_back(newSpring);
      }
    }
 
  // surface y1
  vec2d.clear();
  y = y1 - radius;
  for (z = z1 + radius; z <= z2 - radius + EPS; z += diameter) {
    vec1d.clear();
    for (x = x1 + radius; x <= x2 - radius + EPS; x += diameter) {
      newptcl = new Particle(++memPtclIndex, 5, Vec(x,y,z), radius, young, poisson);
      vec1d.push_back(newptcl);
      particleVec.push_back(newptcl);
    }
    vec2d.push_back(vec1d);
  }
  memBoundary.push_back(vec2d);
  for (int i = 0; i != vec2d.size() ; ++i)
    for (int j = 0; j != vec2d[i].size() ; ++ j) {
      if (j + 1 < vec2d[i].size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], modulus);
	springVec.push_back(newSpring);
      }
      if (i + 1 < vec2d.size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], modulus);
	springVec.push_back(newSpring);
      }
    }

  // surface y2
  vec2d.clear();
  y = y2 + radius;
  for (z = z1 + radius; z <= z2 - radius + EPS; z += diameter) {
    vec1d.clear();
    for (x = x1 + radius; x <= x2 - radius + EPS; x += diameter) {
      newptcl = new Particle(++memPtclIndex, 5, Vec(x,y,z), radius, young, poisson);
      vec1d.push_back(newptcl);
      particleVec.push_back(newptcl);
    }
    vec2d.push_back(vec1d);
  }
  memBoundary.push_back(vec2d);
  for (int i = 0; i != vec2d.size() ; ++i)
    for (int j = 0; j != vec2d[i].size() ; ++ j) {
      if (j + 1 < vec2d[i].size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], modulus);
	springVec.push_back(newSpring);
      }
      if (i + 1 < vec2d.size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], modulus);
	springVec.push_back(newSpring);
      }
    }

  // surface z1
  vec2d.clear();
  z = z1 - radius;
  for (y = y1 + radius; y <= y2 - radius + EPS; y += diameter) {
    vec1d.clear();
    for (x = x1 + radius; x <= x2 - radius + EPS; x += diameter) {
      newptcl = new Particle(++memPtclIndex, 5, Vec(x,y,z), radius, young, poisson);
      vec1d.push_back(newptcl);
      particleVec.push_back(newptcl);
    }
    vec2d.push_back(vec1d);
  }
  memBoundary.push_back(vec2d);
  for (int i = 0; i != vec2d.size() ; ++i)
    for (int j = 0; j != vec2d[i].size() ; ++ j) {
      if (j + 1 < vec2d[i].size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], modulus);
	springVec.push_back(newSpring);
      }
      if (i + 1 < vec2d.size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], modulus);
	springVec.push_back(newSpring);
      }
    }

  // surface z2
  vec2d.clear();
  z = z2 + radius;
  for (y = y1 + radius; y <= y2 - radius + EPS; y += diameter) {
    vec1d.clear();
    for (x = x1 + radius; x <= x2 - radius + EPS; x += diameter) {
      newptcl = new Particle(++memPtclIndex, 5, Vec(x,y,z), radius, young, poisson);
      vec1d.push_back(newptcl);
      particleVec.push_back(newptcl);
    }
    vec2d.push_back(vec1d);
  }
  memBoundary.push_back(vec2d);
  for (int i = 0; i != vec2d.size() ; ++i)
    for (int j = 0; j != vec2d[i].size() ; ++ j) {
      if (j + 1 < vec2d[i].size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], modulus);
	springVec.push_back(newSpring);
      }
      if (i + 1 < vec2d.size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], modulus);
	springVec.push_back(newSpring);
      }
    }

  // membrane particles at the edges of each surface, for example,
  // x1y1 means particles on surface x1 connecting to particles on surface y1
  vector<Particle *> x1y1;
  vector<Particle *> x1y2;
  vector<Particle *> x1z1;
  vector<Particle *> x1z2;

  vector<Particle *> x2y1;
  vector<Particle *> x2y2;
  vector<Particle *> x2z1;
  vector<Particle *> x2z2;

  vector<Particle *> y1x1;
  vector<Particle *> y1x2;
  vector<Particle *> y1z1;
  vector<Particle *> y1z2;

  vector<Particle *> y2x1;
  vector<Particle *> y2x2;
  vector<Particle *> y2z1;
  vector<Particle *> y2z2;

  vector<Particle *> z1x1;
  vector<Particle *> z1x2;
  vector<Particle *> z1y1;
  vector<Particle *> z1y2;

  vector<Particle *> z2x1;
  vector<Particle *> z2x2;
  vector<Particle *> z2y1;
  vector<Particle *> z2y2;

  // find edge particles for each surface
  // memBoundary[0, 1, 2, 3, 4, 5] correspond to 
  // surface     x1 x2 y1 y2 z1 z2 respectively
  // surface x1
  vec2d.clear();
  vec2d = memBoundary[0];
  x1z1  = vec2d[0];
  x1z2  = vec2d[vec2d.size() - 1];
  for (int i = 0; i < vec2d.size(); ++i) {
    x1y1.push_back(vec2d[i][0]);
    x1y2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
  }
  // surface x2
  vec2d.clear();
  vec2d = memBoundary[1];
  x2z1  = vec2d[0];
  x2z2  = vec2d[vec2d.size() - 1];
  for (int i = 0; i < vec2d.size(); ++i) {
    x2y1.push_back(vec2d[i][0]);
    x2y2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
  }
  // surface y1
  vec2d.clear();
  vec2d = memBoundary[2];
  y1z1  = vec2d[0];
  y1z2  = vec2d[vec2d.size() - 1];
  for (int i = 0; i < vec2d.size(); ++i) {
    y1x1.push_back(vec2d[i][0]);
    y1x2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
  }
  // surface y2
  vec2d.clear();
  vec2d = memBoundary[3];
  y2z1  = vec2d[0];
  y2z2  = vec2d[vec2d.size() - 1];
  for (int i = 0; i < vec2d.size(); ++i) {
    y2x1.push_back(vec2d[i][0]);
    y2x2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
  }
  // surface z1
  vec2d.clear();
  vec2d = memBoundary[4];
  z1y1  = vec2d[0];
  z1y2  = vec2d[vec2d.size() - 1];
  for (int i = 0; i < vec2d.size(); ++i) {
    z1x1.push_back(vec2d[i][0]);
    z1x2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
  }
  // surface z2
  vec2d.clear();
  vec2d = memBoundary[5];
  z2y1  = vec2d[0];
  z2y2  = vec2d[vec2d.size() - 1];
  for (int i = 0; i < vec2d.size(); ++i) {
    z2x1.push_back(vec2d[i][0]);
    z2x2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
  }

  // create springs connecting 12 edges of a cube
  // 4 edges on surface x1
  assert(x1y1.size() == y1x1.size());
  for (int i = 0; i < x1y1.size(); ++i) {
    newSpring = new Spring(*x1y1[i], *y1x1[i], modulus);
    springVec.push_back(newSpring);
  }
  assert(x1y2.size() == y2x1.size());
  for (int i = 0; i < x1y2.size(); ++i) {
    newSpring = new Spring(*x1y2[i], *y2x1[i], modulus);
    springVec.push_back(newSpring);
  }
  assert(x1z1.size() == z1x1.size());
  for (int i = 0; i < x1z1.size(); ++i) {
    newSpring = new Spring(*x1z1[i], *z1x1[i], modulus);
    springVec.push_back(newSpring);
  }
  assert(x1z2.size() == z2x1.size());
  for (int i = 0; i < x1z2.size(); ++i) {
    newSpring = new Spring(*x1z2[i], *z2x1[i], modulus);
    springVec.push_back(newSpring);
  }
  // 4 edges on surface x2  
  assert(x2y1.size() == y1x2.size());
  for (int i = 0; i < x2y1.size(); ++i) {
    newSpring = new Spring(*x2y1[i], *y1x2[i], modulus);
    springVec.push_back(newSpring);
  }
  assert(x2y2.size() == y2x2.size());
  for (int i = 0; i < x2y2.size(); ++i) {
    newSpring = new Spring(*x2y2[i], *y2x2[i], modulus);
    springVec.push_back(newSpring);
  }
  assert(x2z1.size() == z1x2.size());
  for (int i = 0; i < x2z1.size(); ++i) {
    newSpring = new Spring(*x2z1[i], *z1x2[i], modulus);
    springVec.push_back(newSpring);
  }
  assert(x2z2.size() == z2x2.size());
  for (int i = 0; i < x2z2.size(); ++i) {
    newSpring = new Spring(*x2z2[i], *z2x2[i], modulus);
    springVec.push_back(newSpring);
  }
  // 2 edges on surface y1 
  assert(y1z1.size() == z1y1.size());
  for (int i = 0; i < y1z1.size(); ++i) {
    newSpring = new Spring(*y1z1[i], *z1y1[i], modulus);
    springVec.push_back(newSpring);
  }
  assert(y1z2.size() == z2y1.size());
  for (int i = 0; i < y1z2.size(); ++i) {
    newSpring = new Spring(*y1z2[i], *z2y1[i], modulus);
    springVec.push_back(newSpring);
  }
  // 2 edges on surface y2
  assert(y2z1.size() == z1y2.size());
  for (int i = 0; i < y2z1.size(); ++i) {
    newSpring = new Spring(*y2z1[i], *z1y2[i], modulus);
    springVec.push_back(newSpring);
  }
  assert(y2z2.size() == z2y2.size());
  for (int i = 0; i < y2z2.size(); ++i) {
    newSpring = new Spring(*y2z2[i], *z2y2[i], modulus);
    springVec.push_back(newSpring);
  }

  printParticle(allParticle);

}


void Assembly::TrimPtclBdryByHeight(REAL height,
			    const char* iniptclfile,
			    const char* ParticleFile)
{
  readParticle(iniptclfile);
  
  std::vector<Particle*>::iterator itr;
  for (itr = particleVec.begin(); itr != particleVec.end(); ){
    if ( (*itr)->getType() == 1 ) { // 1-fixed
      Vec center=(*itr)->getCurrPos();
      if(center.getz() > height)
	{
	  delete (*itr); // release memory
	  itr = particleVec.erase(itr); 
	}
      else {
	(*itr)->setType(10); // 10-ghost
	++itr;
      }
    }
  }
  
  printParticle(ParticleFile);
}


void Assembly::deposit(int   totalSteps,  
		       int   snapNum,
		       int   interval) {
    // pre_1: open streams for output.
    // particlefile and contactfile are used for snapNum at the end.
    progressinf.open("dep_progress"); 
    if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf.precision(OPREC);
    progressinf << setw(OWID) << "iteration"
	        << setw(OWID) << "poss_contact"
	        << setw(OWID) << "actual_contact"
	        << setw(OWID) << "penetration"
	        << setw(OWID) << "avg_normal"
	        << setw(OWID) << "avg_tangt"
	        << setw(OWID) << "avg_velocity"
	        << setw(OWID) << "avg_omga"
	        << setw(OWID) << "avg_force"
	        << setw(OWID) << "avg_moment"
	        << setw(OWID) << "trans_energy"
	        << setw(OWID) << "rotat_energy"
	        << setw(OWID) << "kinet_energy"
	        << setw(OWID) << "poten_energy"
	        << setw(OWID) << "total_energy"
	        << setw(OWID) << "void_ratio"
	        << setw(OWID) << "porosity"
	        << setw(OWID) << "coord_number"
	        << setw(OWID) << "density"
	        << setw(OWID) << "sigma_y1"
	        << setw(OWID) << "sigma_y2"
	        << setw(OWID) << "sigma_x1"
	        << setw(OWID) << "sigma_x2"
	        << setw(OWID) << "sigma_z1"
	        << setw(OWID) << "sigma_z2"
	        << setw(OWID) << "mean_stress"
	        << setw(OWID) << "dimx"
	        << setw(OWID) << "dimy"
	        << setw(OWID) << "dimz"
	        << setw(OWID) << "volume"
	        << setw(OWID) << "epsilon_x"
	        << setw(OWID) << "epsilon_y"
	        << setw(OWID) << "epsilon_z"
	        << setw(OWID) << "epsilon_v"
	        << setw(OWID) << "vibra_t_step"
	        << setw(OWID) << "impact_t_step"
	        << setw(OWID) << "wall_time" << endl;

    debugInf.open("dep_debug");
    if(!debugInf) { cout << "stream error!" << endl; exit(-1); }
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from existing files.
    readParticle("float_particle"); // create container and particles, velocity and omga are set zero. 
    readBoundary("dep_boundary"); // create boundaries

    // pre_3: define variables used in iterations.
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int  stepsnum=0;
    char stepsstr[4];
    char stepsfp[50];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int  bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    // iterations starting ...
    iteration=0; 
    gettimeofday(&time_w1,NULL);
    do
    {
	// 1. create possible boundary particles and contacts between particles.
        //gettimeofday(&time_w1,NULL);
        findContact();
	//gettimeofday(&time_w2,NULL);
	//cout << setw(OWID) << timediffsec(time_w1,time_w2);
	//gettimeofday(&time_w1,NULL);
        findParticleOnBoundary();

	// 2. set particles' forces/moments as zero before each re-calculation,
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles.
	internalForce(avgNormal, avgTangt);

	// 4. calculate boundary forces/moments and apply them to particles.
	boundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();
	//gettimeofday(&time_w2,NULL);
	//cout << setw(OWID) << timediffsec(time_w1,time_w2) << endl;

	// 6. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx(); bulkVolume=l13*l24*l56;
	void_ratio=bulkVolume/getParticleVolume()-1;

	// 7. (1) output particles and contacts information as snapNum.
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, "dep_particle"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, "dep_contact"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    time(&timeStamp);
	    g_timeinf << setw(4) << stepsnum << " " << ctime(&timeStamp) << flush;
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info.
	if (iteration % interval == 0) {
	    gettimeofday(&time_w2,NULL);
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/allParticleVec.size()
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
	                << setw(OWID) << getVibraTimeStep()
	                << setw(OWID) << getImpactTimeStep()
		        << setw(OWID) << timediffsec(time_w1,time_w2)
		        << endl;

	    
	    //	    debugInf << setw(OWID) << iteration
	    //       << setw(OWID) << bdry_penetr[1]
	    //       << setw(OWID) << bdry_penetr[2]
	    //       << setw(OWID) << bdry_penetr[3]
	    //       << setw(OWID) << bdry_penetr[4]
	    //       << setw(OWID) << bdry_penetr[6]
	    //       << setw(OWID) << bdry_cntnum[1]
	    //       << setw(OWID) << bdry_cntnum[2]
	    //       << setw(OWID) << bdry_cntnum[3]
	    //       << setw(OWID) << bdry_cntnum[4]
	    //       << setw(OWID) << bdry_cntnum[6]
	    //       << endl;

	}

	// 8. loop break conditions.

    } while (++iteration < totalSteps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, "dep_particle"); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, "dep_contact"); strcat(stepsfp, "_end");
    printContact(stepsfp);
    g_timeinf << setw(4) << "end" << " " << ctime(&timeStamp) << flush;

    // post_2. close streams
    progressinf.close();
    debugInf.close();
}


void Assembly::depositAfterCavity(int   totalSteps,  
				  int   snapNum,
				  int   interval,
				  const char* iniptclfile,   
				  const char* inibdryfile,
				  const char* inicavefile,
				  const char* ParticleFile, 
				  const char* contactfile,
				  const char* progressfile, 
				  const char* debugfile)
{
    // pre_1: open streams for output.
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf.precision(OPREC);
    progressinf << setw(OWID) << "iteration"
	        << setw(OWID) << "poss_contact"
	        << setw(OWID) << "actual_contact"
	        << setw(OWID) << "penetration"
	        << setw(OWID) << "avg_normal"
	        << setw(OWID) << "avg_tangt"
	        << setw(OWID) << "avg_velocity"
	        << setw(OWID) << "avg_omga"
	        << setw(OWID) << "avg_force"
	        << setw(OWID) << "avg_moment"
	        << setw(OWID) << "trans_energy"
	        << setw(OWID) << "rotat_energy"
	        << setw(OWID) << "kinet_energy"
	        << setw(OWID) << "poten_energy"
	        << setw(OWID) << "total_energy"
	        << setw(OWID) << "void_ratio"
	        << setw(OWID) << "porosity"
	        << setw(OWID) << "coord_number"
	        << setw(OWID) << "density"
	        << setw(OWID) << "sigma_y1"
	        << setw(OWID) << "sigma_y2"
	        << setw(OWID) << "sigma_x1"
	        << setw(OWID) << "sigma_x2"
	        << setw(OWID) << "sigma_z1"
	        << setw(OWID) << "sigma_z2"
	        << setw(OWID) << "mean_stress"
	        << setw(OWID) << "dimx"
	        << setw(OWID) << "dimy"
	        << setw(OWID) << "dimz"
	        << setw(OWID) << "volume"
	        << setw(OWID) << "epsilon_x"
	        << setw(OWID) << "epsilon_y"
	        << setw(OWID) << "epsilon_z"
	        << setw(OWID) << "epsilon_v"
	        << setw(OWID) << "vibra_t_step"
	        << setw(OWID) << "impact_t_step"
	        << setw(OWID) << "wall_time" << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1); }
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from existing files.
    readParticle(iniptclfile);   // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile); // create boundaries.
    readCavityBoundary(inicavefile); // create cavity boundaries

    // pre_3: define variables used in iterations.
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int  stepsnum=0;
    char stepsstr[4];
    char stepsfp[50];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int  bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    // iterations starting ...
    iteration=0; 
    gettimeofday(&time_w1,NULL);
    do
    {
	// 1. create possible boundary particles and contacts between particles.
        findContact();
        findParticleOnBoundary();
	findParticleOnCavity();

	// 2. set particles' forces/moments as zero before each re-calculation,
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles.
	internalForce(avgNormal, avgTangt);

	// 4. calculate boundary forces/moments and apply them to particles.
	boundaryForce(bdry_penetr, bdry_cntnum);
	cavityBoundaryForce();
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();

	// 6. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx(); bulkVolume=l13*l24*l56;
	void_ratio=bulkVolume/getParticleVolume()-1;

	// 7. (1) output particles and contacts information as snapNum.
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    time(&timeStamp);
	    g_timeinf << setw(4) << stepsnum << " " << ctime(&timeStamp) << flush;
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info.
	if (iteration % interval == 0) {
	    gettimeofday(&time_w2,NULL);
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/allParticleVec.size()
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
	                << setw(OWID) << getVibraTimeStep()
	                << setw(OWID) << getImpactTimeStep()
		        << setw(OWID) << timediffsec(time_w1,time_w2)
		        << endl;

	    
	    //debugInf << setw(OWID) << iteration
	    //       << setw(OWID) << bdry_penetr[1]
	    //       << setw(OWID) << bdry_penetr[2]
	    //       << setw(OWID) << bdry_penetr[3]
	    //       << setw(OWID) << bdry_penetr[4]
	    //       << setw(OWID) << bdry_penetr[6]
	    //       << setw(OWID) << bdry_cntnum[1]
	    //       << setw(OWID) << bdry_cntnum[2]
	    //       << setw(OWID) << bdry_cntnum[3]
	    //       << setw(OWID) << bdry_cntnum[4]
	    //       << setw(OWID) << bdry_cntnum[6]
	    //       << endl;
	}

	// 8. loop break conditions.

    } while (++iteration < totalSteps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    g_timeinf << setw(4) << "end" << " " << ctime(&timeStamp) << flush;

    // post_2. close streams
    progressinf.close();
    debugInf.close();
}


void Assembly::deGravitation(int   totalSteps,  
			     int   snapNum,
			     int   interval,
			     bool  toRebuild,
			     const char* iniptclfile,   
			     const char* ParticleFile, 
			     const char* contactfile,
			     const char* progressfile, 
			     const char* debugfile)
{
  // pre_1: open streams for output.
  progressinf.open(progressfile); 
  if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
  progressinf.setf(std::ios::scientific, std::ios::floatfield);
  progressinf.precision(OPREC);
  progressinf << setw(OWID) << "iteration"
	      << setw(OWID) << "poss_contact"
	      << setw(OWID) << "actual_contact"
	      << setw(OWID) << "penetration"
	      << setw(OWID) << "avg_normal"
	      << setw(OWID) << "avg_tangt"
	      << setw(OWID) << "avg_velocity"
	      << setw(OWID) << "avg_omga"
	      << setw(OWID) << "avg_force"
	      << setw(OWID) << "avg_moment"
	      << setw(OWID) << "trans_energy"
	      << setw(OWID) << "rotat_energy"
	      << setw(OWID) << "kinet_energy"
	     << endl;
  
  debugInf.open(debugfile);
  if(!debugInf) { cout << "stream error!" << endl; exit(-1); }
  debugInf.setf(std::ios::scientific, std::ios::floatfield);
  
  // pre_2. create particles from existing files.
  if (toRebuild) readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
  
  // pre_3. define variables used in iterations
  REAL avgNormal=0;
  REAL avgTangt=0;
  int  stepsnum=0;
  char stepsstr[4];
  char stepsfp[50];

  // iterations starting ...
  iteration=0; 
  do
    {
      // 1. find contacts between particles.
      findContact();
      
      // 2. set particles' forces/moments as zero before each re-calculation,
      clearForce();	
      
      // 3. calculate contact forces/moments and apply them to particles.
      internalForce(avgNormal, avgTangt);
      
      // 4. update particles' velocity/omga/position/orientation based on force/moment.
      updateParticle();
      
      // 5. (1) output particles and contacts information as snapNum.
      if (iteration % (totalSteps/snapNum) == 0){
	sprintf(stepsstr, "%03d", stepsnum); 
	strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	printParticle(stepsfp);
	
	sprintf(stepsstr, "%03d", stepsnum); 
	strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	printContact(stepsfp);
	time(&timeStamp);
	g_timeinf << setw(4) << "end" << " " << ctime(&timeStamp) << flush;
	++stepsnum;
      }
      
      // 5. (2) output stress and strain info.
      if (iteration % interval == 0) {
	progressinf << setw(OWID) << iteration
		    << setw(OWID) << getPossContactNum()
		    << setw(OWID) << getActualContactNum()
		    << setw(OWID) << getAvgPenetration()
		    << setw(OWID) << avgNormal
		    << setw(OWID) << avgTangt
		    << setw(OWID) << getAvgVelocity() 
		    << setw(OWID) << getAvgOmga()
		    << setw(OWID) << getAvgForce()   
		    << setw(OWID) << getAvgMoment()
		    << setw(OWID) << getTransEnergy()
		    << setw(OWID) << getRotatEnergy()
		    << setw(OWID) << getKinetEnergy()
		    << endl;
      }
      
      // 7. loop break conditions.
      if (contactVec.size() == 0) break;
      
    } while (++iteration < totalSteps);
  
  // post_1. store the final snapshot of particles & contacts.
  strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
  printParticle(stepsfp);
  
  strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
  printContact(stepsfp);
  g_timeinf << setw(4) << "end" << " " << ctime(&timeStamp) << flush;
  
  // post_2. close streams
  progressinf.close();
  debugInf.close();
}


// actual deposit function for the case of fixed particle boundaries
void Assembly::deposit_p(int   totalSteps,  
			 int   snapNum,
			 int   interval,
			 REAL dimn,
			 REAL rsize,
			 const char* iniptclfile,   
			 const char* ParticleFile, 
			 const char* contactfile,
			 const char* progressfile, 
			 const char* debugfile)
{
    // pre_1: open streams for output.
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "deposit..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average       translational    rotational       "
	        << "kinetic        potential         total           void            sample       coordination"
	        << "       sample           sample          sample          sample          sample          sample"
	        << "          sample          sample          sample         sample           sample         "
	        << " sample          sample          sample          sample          sample" << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "         omga            force           moment         energy           energy          "
	        << "energy         energy            energy          ratio          porosity         number       "
	        << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
	        << "epsilon_v" << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1); }
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from existing files.
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 

    // pre_3: define variables used in iterations.
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    REAL void_ratio=0;

    // iterations starting ...
    iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles.
	findContact();

	// 2. set particles' forces/moments as zero before each re-calculation,
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles.
	internalForce(avgNormal, avgTangt);

	// 4. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();

	// 5. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - (-dimn/2);
	l24=dimn*rsize;
	l13=dimn*rsize;
	bulkVolume=l13*l24*l56;
	void_ratio=bulkVolume/getParticleVolume()-1;

	// 6. (1) output particles and contacts information as snapNum.
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 6. (2) output statistics info.
	if (iteration % interval == 0) {
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
		        << endl;
	}

	// 7. loop break conditions.


    } while (++iteration < totalSteps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    // post_2. close streams
    progressinf.close();
    debugInf.close();
}


// squeeze paticles inside a container by moving the boundaries
void Assembly::squeeze(int   totalSteps,  
		       int   init_steps,
		       int   snapNum,
		       int   interval,
		       int   flag,
		       const char* iniptclfile,   
		       const char* inibdryfile,
		       const char* ParticleFile, 
		       const char* boundaryfile,
		       const char* contactfile,
		       const char* progressfile, 
		       const char* debugfile)
{
    // pre_1: open streams for output.
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "deposit..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average       translational    rotational       "
	        << "kinetic        potential         total           void            sample       coordination"
	        << "       sample           sample          sample          sample          sample          sample"
	        << "          sample          sample          sample         sample           sample         "
	        << " sample          sample          sample          sample          sample" << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "         omga            force           moment         energy           energy          "
	        << "energy         energy            energy          ratio          porosity         number       "
	        << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
	        << "epsilon_v" << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1); }
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from existing files.
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries.

    // pre_3: define variables used in iterations.
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];

    int         mid[2]={1,3};    // boundary 1 and 3
    UPDATECTL   midctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    // iterations starting ...
    iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles.
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces/moments as zero before each re-calculation,
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles.
	internalForce(avgNormal, avgTangt);

	// 4. calculate boundary forces/moments and apply them to particles.
	boundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();

	// 6. calculate sample void ratio.
	l56=getTopFreeParticlePosition().getz() -getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx(); bulkVolume=l13*l24*l56;
	void_ratio=bulkVolume/getParticleVolume()-1;

	// displacement control
	if (iteration > init_steps) {
	    if (flag==1) // loosen, totally remove the wall
		midctl[0].tran=Vec(TIMESTEP*1.0e+0*flag,0,0);
	    else         // squeeze
		midctl[0].tran=Vec(TIMESTEP*5.0e-3*flag,0,0);
	    updateBoundary(mid,midctl,2);
	}

	// 7. (1) output particles and contacts information as snapNum.
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info.
	if (iteration % interval == 0) {
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    debugInf << setw(OWID) << iteration
		       << setw(OWID) << bdry_penetr[1]
		       << setw(OWID) << bdry_penetr[2]
		       << setw(OWID) << bdry_penetr[3]
		       << setw(OWID) << bdry_penetr[4]
		       << setw(OWID) << bdry_penetr[6]
		       << setw(OWID) << bdry_cntnum[1]
		       << setw(OWID) << bdry_cntnum[2]
		       << setw(OWID) << bdry_cntnum[3]
		       << setw(OWID) << bdry_cntnum[4]
		       << setw(OWID) << bdry_cntnum[6]
		       << endl;

	}

	// 8. loop break conditions.

    } while (++iteration < totalSteps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);

    // post_2. close streams
    progressinf.close();
    debugInf.close();
}


// Isotropically compress floating particles to a specific confining pressure, which is usually a low
// value in order to create an intial status. Force boundaries are used. This process may be not 
// physically true.
void Assembly::isotropic(int   totalSteps,
			 int   snapNum, 
			 int   interval,
			 REAL  sigma,			  
			 const char* iniptclfile,   
			 const char* inibdryfile,
			 const char* ParticleFile, 
			 const char* boundaryfile,
			 const char* contactfile,  
			 const char* progressfile,
			 const char* balancedfile, 
			 const char* debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout << "stream error!" << endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "isotropic..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout << "stream error!" << endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "isotropic..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3: define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         mid[2]={1,3};    // boundary 1 and 3
    int         max[2]={2,4};    // boundary 2 and 4
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   midctl[2];
    UPDATECTL   maxctl[2];
    UPDATECTL   minctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    // iterations start here...
    iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	boundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    bulkVolume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=bulkVolume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=Vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=Vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=Vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=Vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	if (sigma2_1<sigma)
	    midctl[0].tran=Vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=Vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma)
	    midctl[1].tran=Vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=Vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma)
	    maxctl[0].tran=Vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=Vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma)
	    maxctl[1].tran=Vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=Vec(0,-TIMESTEP*RELEASE_RATE,0);
	
	updateBoundary(min,minctl,2);
	updateBoundary(mid,midctl,2);
	updateBoundary(max,maxctl,2);
	updateBoundary6();
	
	// 7. (1) output particles and contacts information
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (iteration % interval == 0 ){
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    debugInf << setw(OWID) << iteration
		       << setw(OWID) << bdry_penetr[1]
		       << setw(OWID) << bdry_penetr[2]
		       << setw(OWID) << bdry_penetr[3]
		       << setw(OWID) << bdry_penetr[4]
		       << setw(OWID) << bdry_penetr[5]
		       << setw(OWID) << bdry_penetr[6]
		       << setw(OWID) << bdry_cntnum[1]
		       << setw(OWID) << bdry_cntnum[2]
		       << setw(OWID) << bdry_cntnum[3]
		       << setw(OWID) << bdry_cntnum[4]
		       << setw(OWID) << bdry_cntnum[5]
		       << setw(OWID) << bdry_cntnum[6]
		       << endl;
	}

	// 8. loop break condition
	if (   fabs(sigma1_1-sigma)/sigma < STRESS_ERROR && fabs(sigma1_2-sigma)/sigma < STRESS_ERROR
	    && fabs(sigma2_1-sigma)/sigma < STRESS_ERROR && fabs(sigma2_2-sigma)/sigma < STRESS_ERROR
	    && fabs(sigma3_1-sigma)/sigma < STRESS_ERROR && fabs(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()    
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()    
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    break;
	}

    } while (++iteration < totalSteps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile);  strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    debugInf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_a. This function
// increases confining pressure step by step to sigma_b, making it possible to find equilibrium 
// state where particle pressure equals confining pressure. Force boundaries are used
void Assembly::isotropic(int   totalSteps,
			 int   snapNum, 
			 int   interval,
			 REAL sigma_a,
			 REAL sigma_b,
			 int   sigma_division,
			 const char* iniptclfile,   
			 const char* inibdryfile,
			 const char* ParticleFile, 
			 const char* boundaryfile,
			 const char* contactfile,  
			 const char* progressfile,
			 const char* balancedfile, 
			 const char* debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout << "stream error!" << endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "isotropic..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout << "stream error!" << endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "isotropic..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3: define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         mid[2]={1,3};    // boundary 1 and 3
    int         max[2]={2,4};    // boundary 2 and 4
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   midctl[2];
    UPDATECTL   maxctl[2];
    UPDATECTL   minctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    REAL sigma=sigma_a;
    REAL sigma_inc=(sigma_b-sigma_a)/sigma_division;

    // iterations start here...
    iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	boundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    bulkVolume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=bulkVolume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=Vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=Vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=Vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=Vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	if (sigma2_1<sigma)
	    midctl[0].tran=Vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=Vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma)
	    midctl[1].tran=Vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=Vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma)
	    maxctl[0].tran=Vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=Vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma)
	    maxctl[1].tran=Vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=Vec(0,-TIMESTEP*RELEASE_RATE,0);
	
	updateBoundary(min,minctl,2);
	updateBoundary(mid,midctl,2);
	updateBoundary(max,maxctl,2);
	updateBoundary6();
	
	// 7. (1) output particles and contacts information
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (iteration % interval == 0 ){
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    debugInf << setw(OWID) << iteration
		       << setw(OWID) << bdry_penetr[1]
		       << setw(OWID) << bdry_penetr[2]
		       << setw(OWID) << bdry_penetr[3]
		       << setw(OWID) << bdry_penetr[4]
		       << setw(OWID) << bdry_penetr[5]
		       << setw(OWID) << bdry_penetr[6]
		       << setw(OWID) << bdry_cntnum[1]
		       << setw(OWID) << bdry_cntnum[2]
		       << setw(OWID) << bdry_cntnum[3]
		       << setw(OWID) << bdry_cntnum[4]
		       << setw(OWID) << bdry_cntnum[5]
		       << setw(OWID) << bdry_cntnum[6]
		       << endl;
	}

	// 8. find the balanced status and increase confining pressure
	if (   fabs(sigma1_1-sigma)/sigma < STRESS_ERROR && fabs(sigma1_2-sigma)/sigma < STRESS_ERROR
	    && fabs(sigma2_1-sigma)/sigma < STRESS_ERROR && fabs(sigma2_2-sigma)/sigma < STRESS_ERROR
	    && fabs(sigma3_1-sigma)/sigma < STRESS_ERROR && fabs(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()    
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    sigma += sigma_inc;
	}

	// 9. loop break condition
	if (   fabs(sigma1_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma1_2-sigma_b)/sigma_b < STRESS_ERROR
	    && fabs(sigma2_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma2_2-sigma_b)/sigma_b < STRESS_ERROR
	    && fabs(sigma3_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma3_2-sigma_b)/sigma_b < STRESS_ERROR ) {
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()    
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    break;
	}
	
    } while (++iteration < totalSteps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    debugInf.close();
}


// loading-unloading-reloading of isotropic compression
// the stress path is defined by sigma_points and sigma_values[]
void Assembly::isotropic(int   totalSteps,  
			 int   snapNum, 
			 int   interval,
			 int   sigma_points,  
			 REAL sigma_values[],  
			 int   sigma_division,	  
			 const char* iniptclfile,  
			 const char* inibdryfile,
			 const char* ParticleFile, 
			 const char* boundaryfile,
			 const char* contactfile,  
			 const char* progressfile,
			 const char* balancedfile, 
			 const char* debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout << "stream error!" << endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "isotropic..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout << "stream error!" << endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "isotropic..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3: define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         mid[2]={1,3};    // boundary 1 and 3
    int         max[2]={2,4};    // boundary 2 and 4
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   midctl[2];
    UPDATECTL   maxctl[2];
    UPDATECTL   minctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    int  i=0;
    REAL sigma=sigma_values[i];
    REAL sigma_inc=(sigma_values[i+1]-sigma_values[i])/sigma_division;
    REAL sigma_b=sigma_values[sigma_points-1];

    // iterations start here...
    iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();
	
	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	boundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    bulkVolume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=bulkVolume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=Vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=Vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=Vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=Vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	if (sigma2_1<sigma)
	    midctl[0].tran=Vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=Vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma)
	    midctl[1].tran=Vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=Vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma)
	    maxctl[0].tran=Vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=Vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma)
	    maxctl[1].tran=Vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=Vec(0,-TIMESTEP*RELEASE_RATE,0);
	
	updateBoundary(min,minctl,2);
	updateBoundary(mid,midctl,2);
	updateBoundary(max,maxctl,2);
	updateBoundary6();
	
	// 7. (1) output particles and contacts information
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (iteration % interval == 0 ){
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    debugInf << setw(OWID) << iteration
		       << setw(OWID) << bdry_penetr[1]
		       << setw(OWID) << bdry_penetr[2]
		       << setw(OWID) << bdry_penetr[3]
		       << setw(OWID) << bdry_penetr[4]
		       << setw(OWID) << bdry_penetr[5]
		       << setw(OWID) << bdry_penetr[6]
		       << setw(OWID) << bdry_cntnum[1]
		       << setw(OWID) << bdry_cntnum[2]
		       << setw(OWID) << bdry_cntnum[3]
		       << setw(OWID) << bdry_cntnum[4]
		       << setw(OWID) << bdry_cntnum[5]
		       << setw(OWID) << bdry_cntnum[6]
		       << endl;
	}

	// 8. find the balanced status and increase confining pressure
	if (   fabs(sigma1_1-sigma)/sigma < STRESS_ERROR && fabs(sigma1_2-sigma)/sigma < STRESS_ERROR
	    && fabs(sigma2_1-sigma)/sigma < STRESS_ERROR && fabs(sigma2_2-sigma)/sigma < STRESS_ERROR
	    && fabs(sigma3_1-sigma)/sigma < STRESS_ERROR && fabs(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()    
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    sigma += sigma_inc;
	    if (sigma==sigma_values[i+1]) {
		i++;
		sigma=sigma_values[i];
		sigma_inc=(sigma_values[i+1]-sigma_values[i])/sigma_division;
	    }

	}

	// 9. loop break condition
	if (   fabs(sigma1_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma1_2-sigma_b)/sigma_b < STRESS_ERROR
	    && fabs(sigma2_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma2_2-sigma_b)/sigma_b < STRESS_ERROR
	    && fabs(sigma3_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma3_2-sigma_b)/sigma_b < STRESS_ERROR ) {
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()    
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    break;
	}
	
    } while (++iteration < totalSteps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    debugInf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_3. This function
// increases vertical pressure step by step to sigma_1, thus making it possible to find out
// balanced status where top & bottom particle pressure equals major principle stress. 
// Side boundaries are fixed, top and bottom plates are force-controlled.
void Assembly::odometer(int   totalSteps,  
			int   snapNum, 
			int   interval,
			REAL sigma_3,     
			REAL sigma_1,    
			int   sigma_division,			  
			const char* iniptclfile,  
			const char* inibdryfile,
			const char* ParticleFile,
			const char* boundaryfile,
			const char* contactfile,  
			const char* progressfile,
			const char* balancedfile, 
			const char* debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout << "stream error!" << endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "odometer..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout << "stream error!" << endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "odometer..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries
 
    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];

    int min[2]={5,6};    // minimum stress acting on boundary 5 and 6
    UPDATECTL minctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    REAL sigma=sigma_3;
    REAL sigma_inc=(sigma_1-sigma_3)/sigma_division;

    // iterations start here...
    iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();
	
	// 2. set particles' forces and moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces and moments
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces
	boundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/displacement based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    bulkVolume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=bulkVolume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=Vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=Vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=Vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=Vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	updateBoundary(min,minctl,2);
	updateBoundary6();

	// 7. (1) output particles and contacts information
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (iteration % interval == 0){
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	}

	// 8. find balanced status of odometer compression
	if (fabs(sigma3_1-sigma)/sigma < STRESS_ERROR && fabs(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()    
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    sigma += sigma_inc;
	}

	// 9. loop break condition
	if (fabs(sigma3_1-sigma_1)/sigma_1 < STRESS_ERROR && fabs(sigma3_2-sigma_1)/sigma_1 < STRESS_ERROR) {
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()    
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    break;
	}
    } while (++iteration < totalSteps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    debugInf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_3. This function
// increases vertical pressure step by step to sigma_1, thus making it possible to find out
// balanced status where top & bottom particle pressure equals major principle stress. 
// Side boundaries are fixed, top and bottom plates are force-controlled. Unloading path is
// applied.
void Assembly::odometer(int   totalSteps,  
			int   snapNum, 
			int   interval,
			int   sigma_points,  
			REAL sigma_values[],  
			int   sigma_division,			  
			const char* iniptclfile,  
			const char* inibdryfile,
			const char* ParticleFile, 
			const char* boundaryfile,
			const char* contactfile,  
			const char* progressfile,
			const char* balancedfile, 
			const char* debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout << "stream error!" << endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "odometer..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout << "stream error!" << endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "odometer..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries
 
    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];

    int min[2]={5,6};    // minimum stress acting on boundary 5 and 6
    UPDATECTL minctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }


    int  i=0;
    REAL sigma=sigma_values[i];
    REAL sigma_inc=(sigma_values[i+1]-sigma_values[i])/sigma_division;
    REAL sigma_b=sigma_values[sigma_points-1];

    // iterations start here...
    iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces and moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces and moments
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces
	boundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/displacement based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    bulkVolume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=bulkVolume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=Vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=Vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=Vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=Vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	updateBoundary(min,minctl,2);
	updateBoundary6();

	// 7. (1) output particles and contacts information
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (iteration % interval == 0){
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	}

	// 8. find balanced status of odometer compression
	if (fabs(sigma3_1-sigma)/sigma < STRESS_ERROR && fabs(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()    
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    sigma += sigma_inc;
	    if (sigma==sigma_values[i+1]) {
		i++;
		sigma=sigma_values[i];
		sigma_inc=(sigma_values[i+1]-sigma_values[i])/sigma_division;
	    }
	}

	// 9. loop break condition
	if (fabs(sigma3_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma3_2-sigma_b)/sigma_b < STRESS_ERROR) {
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()    
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    break;
	}
    } while (++iteration < totalSteps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    debugInf.close();
}


void Assembly::iso_MemBdry(int   totalSteps,  
			   int   snapNum, 
			   int   interval,
			   REAL  sigma3,
			   REAL  rRadius,
			   bool  toRebuild,
			   const char* iniptclfile, 
			   const char* ParticleFile,
			   const char* contactfile, 
			   const char* progressfile,
			   const char* debugfile) 
{
  // pre_1: open streams for output
  // ParticleFile and contactfile are used for snapNum at the end.
  progressinf.open(progressfile);
  if(!progressinf) { cout << "stream error!" << endl; exit(-1);}
  progressinf.setf(std::ios::scientific, std::ios::floatfield);
  progressinf.precision(OPREC);
  progressinf << setw(OWID) << "iteration"
	      << setw(OWID) << "poss_contact"
	      << setw(OWID) << "actual_contact"
	      << setw(OWID) << "penetration"
	      << setw(OWID) << "avg_normal"
	      << setw(OWID) << "avg_tangt"
	      << setw(OWID) << "avg_velocity"
	      << setw(OWID) << "avg_omga"
	      << setw(OWID) << "avg_force"
	      << setw(OWID) << "avg_moment"
	      << setw(OWID) << "trans_energy"
	      << setw(OWID) << "rotat_energy"
	      << setw(OWID) << "kinet_energy"
	      << endl;
  
  debugInf.open(debugfile);
  if(!debugInf) { cout << "stream error!" << endl; exit(-1);}
  debugInf.setf(std::ios::scientific, std::ios::floatfield);
  
  // pre_2. create particles from file and calculate forces caused by hydraulic pressure
  if (toRebuild) readParticle(iniptclfile);
  REAL radius = gradation.getMinPtclRadius();
  if (gradation.getSize().size() == 1 &&
      gradation.getPtclRatioBA() == 1.0 && 
      gradation.getPtclRatioCA() == 1.0)
    radius *= rRadius; // determine how tiny the boundary particles are
  REAL mag  = radius*radius*4*sigma3;
  Vec v1 = allContainer.getMinCorner();
  Vec v2 = allContainer.getMaxCorner();
  REAL x1 = v1.getx();
  REAL y1 = v1.gety();
  REAL z1 = v1.getz();
  REAL x2 = v2.getx();
  REAL y2 = v2.gety();
  REAL z2 = v2.getz();
  std::vector<Particle*>::const_iterator  it;
  Vec pos;
  for (it=particleVec.begin();it!=particleVec.end();++it)
    {
      pos = (*it)->getCurrPos();
      if (pos.getx() < x1)
	(*it)->setConstForce( Vec(mag, 0, 0) );
      else if (pos.getx() > x2)
	(*it)->setConstForce( Vec(-mag, 0, 0) );
      else if (pos.gety() < y1)
	(*it)->setConstForce( Vec(0, mag, 0) );
      else if (pos.gety() > y2)
	(*it)->setConstForce( Vec(0, -mag, 0) );
      else if (pos.getz() < z1)
	(*it)->setConstForce( Vec(0, 0, mag) );
      else if (pos.getz() > z2)
	(*it)->setConstForce( Vec(0, 0, -mag) );
    }

  // pre_3. define variables used in iterations
  REAL avgNormal=0;
  REAL avgTangt=0;
  int  stepsnum=0;
  char stepsstr[4];
  char stepsfp[50];
  
  // iterations start here...
  iteration=0;
  do 
    {
      // 1. find contacts between particles
      findContact();

      // 2. set particles forces/moments to zero before each re-calculation
      clearForce(); // const_force/moment NOT cleared by this call	
      
      // 3. calculate contact forces/moments and apply them to particles
      internalForce(avgNormal, avgTangt);

      // 4. calculate and apply spring forces to boundary particles
      springForce();
      
      // 5. update particles velocity/omga/position/orientation based on force/moment
      updateParticle();
      
      // 6. (1) output particles and contacts information
      if (iteration % (totalSteps/snapNum) == 0){
	sprintf(stepsstr, "%03d", stepsnum); 
	strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	printParticle(stepsfp);

	strcpy(stepsfp,"iso_membrane_"); strcat(stepsfp, stepsstr);
	printMemParticle(stepsfp);
	strcpy(stepsfp,"iso_spring_"); strcat(stepsfp, stepsstr); strcat(stepsfp, ".dat");
	plotSpring(stepsfp);
	
	sprintf(stepsstr, "%03d", stepsnum); 
	strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr); 
	printContact(stepsfp);
	time(&timeStamp);
	g_timeinf << setw(4) << stepsnum << " " << ctime(&timeStamp) << flush;
	++stepsnum;
      }
      
      // 6. (2) output stress and strain info
	if (iteration % interval == 0 ){
	  progressinf << setw(OWID) << iteration
		      << setw(OWID) << getPossContactNum()
		      << setw(OWID) << getActualContactNum()
		      << setw(OWID) << getAvgPenetration()
		      << setw(OWID) << avgNormal
		      << setw(OWID) << avgTangt
		      << setw(OWID) << getAvgVelocity() 
		      << setw(OWID) << getAvgOmga()
		      << setw(OWID) << getAvgForce()   
		      << setw(OWID) << getAvgMoment()
		      << setw(OWID) << getTransEnergy()
		      << setw(OWID) << getRotatEnergy()
		      << setw(OWID) << getKinetEnergy()
		      << endl;
	}
	  
    } while (++iteration < totalSteps);
  
  // post_1. store the final snapshot of particles, contacts and boundaries.
  strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
  printParticle(stepsfp);

  strcpy(stepsfp, "iso_membrane_end");
  printMemParticle(stepsfp);
  strcpy(stepsfp, "iso_spring_end.dat");
  plotSpring(stepsfp);  

  strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
  printContact(stepsfp);
  g_timeinf << setw(4) << "end" << " " << ctime(&timeStamp) << flush;
  
  // post_2. close streams
  progressinf.close();
  debugInf.close();
}


// This function initializes triaxial sample to a certain confining pressure.
void Assembly::triaxialPtclBdryIni(int   totalSteps,  
				   int   snapNum, 
				   int   interval,
				   REAL  sigma,
				   const char* iniptclfile, 
				   const char* inibdryfile,
				   const char* ParticleFile,
				   const char* boundaryfile,
				   const char* contactfile, 
				   const char* progressfile,
				   const char* debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout << "stream error!" << endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "triaxial..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l56= 0;
    REAL sigma3_1, sigma3_2;
    REAL epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   minctl[2];

    // iterations start here...
    iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	boundaryForce();

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	sigma3_1=vfabs(getNormalForce(5))/2.5e-3; sigma3_2=vfabs(getNormalForce(6))/2.5e-3;

	// force control
	if (sigma3_1 < sigma)
	    minctl[0].tran=Vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=Vec(0,0, TIMESTEP*COMPRESS_RATE);

	if (sigma3_2 < sigma)
	    minctl[1].tran=Vec(0,0, TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=Vec(0,0,-TIMESTEP*COMPRESS_RATE);

	updateBoundary(min,minctl,2);
	
	// 7. (1) output particles and contacts information
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	l56=getApt(5).getz()-getApt(6).getz();
	epsilon_h = (H0-l56)/H0;
	if (iteration % interval == 0 ){
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << 0 << setw(OWID) << 0
		        << setw(OWID) << 0 << setw(OWID) << 0
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()
		        << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << l56
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
		        << endl;

	}

	// 9. loop break condition: through displacement control mechanism
	if (   fabs(sigma3_1-sigma)/sigma < STRESS_ERROR && fabs(sigma3_2-sigma)/sigma < STRESS_ERROR )
	       break;
	
    } while (++iteration < totalSteps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    debugInf.close();
}


// This function performs triaxial compression test.
// Displacement boundaries are used in axial direction.
void Assembly::triaxialPtclBdry(int   totalSteps,  
				int   snapNum, 
				int   interval,
				const char* iniptclfile, 
				const char* inibdryfile,
				const char* ParticleFile,
				const char* boundaryfile,
				const char* contactfile, 
				const char* progressfile,
				const char* balancedfile,
				const char* debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout << "stream error!" << endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "triaxial..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout << "stream error!" << endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "triaxial..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l56= 0;
    REAL sigma3_1, sigma3_2;
    REAL epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   minctl[2];

    // iterations start here...
    iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	boundaryForce();

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	sigma3_1=vfabs(getNormalForce(5))/2.5e-3; sigma3_2=vfabs(getNormalForce(6))/2.5e-3;

	// displacement control
	if(iteration < 100001) {
	minctl[0].tran=Vec(0,0,-TIMESTEP*COMPRESS_RATE);
	minctl[1].tran=Vec(0,0, TIMESTEP*COMPRESS_RATE);

	updateBoundary(min,minctl,2);
	}
	// 7. (1) output particles and contacts information
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	l56=getApt(5).getz()-getApt(6).getz();
	epsilon_h = (H0-l56)/H0;
	if (iteration % interval == 0 ){
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << 0 << setw(OWID) << 0
		        << setw(OWID) << 0 << setw(OWID) << 0
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()
		        << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << l56
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
		        << endl;

	}

	// 9. loop break condition: through displacement control mechanism
	
    } while (++iteration < totalSteps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    debugInf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_a. This function
// performs triaxial compression test. Displacement boundaries are used in axial direction.
void Assembly::triaxial(int   totalSteps,  
			int   snapNum, 
			int   interval,
			REAL  sigma_a,	  
			const char* iniptclfile, 
			const char* inibdryfile,
			const char* ParticleFile,
			const char* boundaryfile,
			const char* contactfile, 
			const char* progressfile,
			const char* balancedfile,
			const char* debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout << "stream error!" << endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf.precision(OPREC);
    progressinf << setw(OWID) << "iteration"
	        << setw(OWID) << "possible"
	        << setw(OWID) << "actual"
	        << setw(OWID) << "average"
	        << setw(OWID) << "average"
	        << setw(OWID) << "average"
	        << setw(OWID) << "average"
	        << setw(OWID) << "average"
	        << setw(OWID) << "average"
	        << setw(OWID) << "average"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "void"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "coordinate"
	        << setw(OWID) << "vibra"
	        << setw(OWID) << "impact"
	        << setw(OWID) << "wall-clock" << endl
	        << setw(OWID) << "number"
	        << setw(OWID) << "contacts"
	        << setw(OWID) << "contacts"
	        << setw(OWID) << "penetration"
	        << setw(OWID) << "contact_normal"
	        << setw(OWID) << "contact_tangt"
	        << setw(OWID) << "velocity"
	        << setw(OWID) << "omga"
	        << setw(OWID) << "force"
	        << setw(OWID) << "moment"
	        << setw(OWID) << "density"
	        << setw(OWID) << "sigma1_1"
	        << setw(OWID) << "sigma1_2"
	        << setw(OWID) << "sigma2_1"
	        << setw(OWID) << "sigma2_2"
	        << setw(OWID) << "sigma3_1"
	        << setw(OWID) << "sigma3_2"
	        << setw(OWID) << "mean_stress"
	        << setw(OWID) << "width"
	        << setw(OWID) << "length"
	        << setw(OWID) << "height"
	        << setw(OWID) << "volume"
	        << setw(OWID) << "epsilon_w"
	        << setw(OWID) << "epsilon_l"
	        << setw(OWID) << "epsilon_h"
	        << setw(OWID) << "epsilon_v"
	        << setw(OWID) << "ratio"
	        << setw(OWID) << "porosity"
	        << setw(OWID) << "number"
	        << setw(OWID) << "t_step"
	        << setw(OWID) << "t_step"
	        << setw(OWID) << "time" << endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout << "stream error!" << endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf.precision(OPREC);
    balancedinf << setw(OWID) << "iteration"
	        << setw(OWID) << "possible"
	        << setw(OWID) << "actual"
	        << setw(OWID) << "average"
	        << setw(OWID) << "average"
	        << setw(OWID) << "average"
	        << setw(OWID) << "average"
	        << setw(OWID) << "average"
	        << setw(OWID) << "average"
	        << setw(OWID) << "average"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "void"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "coordinate"
	        << setw(OWID) << "vibra"
	        << setw(OWID) << "impact"
	        << setw(OWID) << "wall-clock" << endl
	        << setw(OWID) << "number"
	        << setw(OWID) << "contacts"
	        << setw(OWID) << "contacts"
	        << setw(OWID) << "penetration"
	        << setw(OWID) << "contact_normal"
	        << setw(OWID) << "contact_tangt"
	        << setw(OWID) << "velocity"
	        << setw(OWID) << "omga"
	        << setw(OWID) << "force"
	        << setw(OWID) << "moment"
	        << setw(OWID) << "density"
	        << setw(OWID) << "sigma1_1"
	        << setw(OWID) << "sigma1_2"
	        << setw(OWID) << "sigma2_1"
	        << setw(OWID) << "sigma2_2"
	        << setw(OWID) << "sigma3_1"
	        << setw(OWID) << "sigma3_2"
	        << setw(OWID) << "mean_stress"
	        << setw(OWID) << "width"
	        << setw(OWID) << "length"
	        << setw(OWID) << "height"
	        << setw(OWID) << "volume"
	        << setw(OWID) << "epsilon_w"
	        << setw(OWID) << "epsilon_l"
	        << setw(OWID) << "epsilon_h"
	        << setw(OWID) << "epsilon_v"
	        << setw(OWID) << "ratio"
	        << setw(OWID) << "porosity"
	        << setw(OWID) << "number"
	        << setw(OWID) << "t_step"
	        << setw(OWID) << "t_step"
	        << setw(OWID) << "time" << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         mid[2]={1,3};    // boundary 1 and 3
    int         max[2]={2,4};    // boundary 2 and 4
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   midctl[2];
    UPDATECTL   maxctl[2];
    UPDATECTL   minctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    // iterations start here...
    iteration=0;
    gettimeofday(&time_w1,NULL);
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	boundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    bulkVolume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=bulkVolume/getParticleVolume()-1;

	// displacement control
	minctl[0].tran=Vec(0,0,-TIMESTEP*COMPRESS_RATE);
	minctl[1].tran=Vec(0,0, TIMESTEP*COMPRESS_RATE);

	// force control
	if (sigma2_1<sigma_a)
	    midctl[0].tran=Vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=Vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma_a)
	    midctl[1].tran=Vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=Vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma_a)
	    maxctl[0].tran=Vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=Vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma_a)
	    maxctl[1].tran=Vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=Vec(0,-TIMESTEP*RELEASE_RATE,0);
	
	updateBoundary(min,minctl,2);
	updateBoundary(mid,midctl,2);
	updateBoundary(max,maxctl,2);
	updateBoundary6();
	
	// 7. (1) output particles and contacts information
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    time(&timeStamp);
	    g_timeinf << setw(4) << stepsnum << " " << ctime(&timeStamp) << flush;
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (iteration % interval == 0 ){
	    gettimeofday(&time_w2,NULL);
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
	                << setw(OWID) << getVibraTimeStep()
	                << setw(OWID) << getImpactTimeStep()
		        << setw(OWID) << timediffsec(time_w1,time_w2)
		        << endl;
	    debugInf << setw(OWID) << iteration
		       << setw(OWID) << getTransEnergy()
		       << setw(OWID) << getRotatEnergy()
		       << setw(OWID) << bdry_penetr[1]
		       << setw(OWID) << bdry_penetr[2]
		       << setw(OWID) << bdry_penetr[3]
		       << setw(OWID) << bdry_penetr[4]
		       << setw(OWID) << bdry_penetr[5]
		       << setw(OWID) << bdry_penetr[6]
		       << setw(OWID) << bdry_cntnum[1]
		       << setw(OWID) << bdry_cntnum[2]
		       << setw(OWID) << bdry_cntnum[3]
		       << setw(OWID) << bdry_cntnum[4]
		       << setw(OWID) << bdry_cntnum[5]
		       << setw(OWID) << bdry_cntnum[6]
		       << endl;
	}

	// Most time it is balanced, so use progressinf instead.
	// 8. find the balanced status and increase confining pressure
	if (   fabs(sigma1_1-sigma_a)/sigma_a < STRESS_ERROR && fabs(sigma1_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabs(sigma2_1-sigma_a)/sigma_a < STRESS_ERROR && fabs(sigma2_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabs(sigma3_1-sigma3_2)/(sigma3_1+sigma3_2)*2<=0.05) {
	    balancedinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()    
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
	                << setw(OWID) << getVibraTimeStep()
	                << setw(OWID) << getImpactTimeStep()
		        << setw(OWID) << timediffsec(time_w1,time_w2)
		        << endl;
	}

	// 9. loop break condition: through displacement control mechanism
	
    } while (++iteration < totalSteps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    g_timeinf << setw(4) << "end" << " " << ctime(&timeStamp) << flush;

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    debugInf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_a. This function
// performs triaxial compression test with unloading. Displacement boundaries are used in 
// axial direction.
void Assembly::triaxial(int   totalSteps,  
			int   unload_step,
			int   snapNum, 
			int   interval,
			REAL sigma_a,	  
			const char* iniptclfile,  
			const char* inibdryfile,
			const char* ParticleFile,
			const char* boundaryfile,
			const char* contactfile,
			const char* progressfile,
			const char* balancedfile,
			const char* debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout << "stream error!" << endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "triaxial..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout << "stream error!" << endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "triaxial..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         mid[2]={1,3};    // boundary 1 and 3
    int         max[2]={2,4};    // boundary 2 and 4
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   midctl[2];
    UPDATECTL   maxctl[2];
    UPDATECTL   minctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    // iterations start here...
    bool reload=false;
    iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();
	
	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	boundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    bulkVolume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=bulkVolume/getParticleVolume()-1;

	// displacement control
	if (iteration <= unload_step){ //loading
	    minctl[0].tran=Vec(0,0,-TIMESTEP*COMPRESS_RATE);
	    minctl[1].tran=Vec(0,0, TIMESTEP*COMPRESS_RATE);
	}
	else { 
	    if (reload==false) { // unloading
		if (fabs(sigma3_1-sigma_a)/sigma_a > STRESS_ERROR && 
		    fabs(sigma3_2-sigma_a)/sigma_a > STRESS_ERROR){
		    minctl[0].tran=Vec(0,0, TIMESTEP*COMPRESS_RATE);
		    minctl[1].tran=Vec(0,0,-TIMESTEP*COMPRESS_RATE);
		}
		else  // reloading
		    reload=true;
	    }
	    else {
		minctl[0].tran=Vec(0,0,-TIMESTEP*COMPRESS_RATE);
		minctl[1].tran=Vec(0,0, TIMESTEP*COMPRESS_RATE);
	    }
	}
	
	// force control
	if (sigma2_1<sigma_a)
	    midctl[0].tran=Vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=Vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma_a)
	    midctl[1].tran=Vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=Vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma_a)
	    maxctl[0].tran=Vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=Vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma_a)
	    maxctl[1].tran=Vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=Vec(0,-TIMESTEP*RELEASE_RATE,0);
	
	updateBoundary(min,minctl,2);
	updateBoundary(mid,midctl,2);
	updateBoundary(max,maxctl,2);
	updateBoundary6();
	
	// 7. (1) output particles and contacts information
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (iteration % interval == 0 ){
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    debugInf << setw(OWID) << iteration
		       << setw(OWID) << bdry_penetr[1]
		       << setw(OWID) << bdry_penetr[2]
		       << setw(OWID) << bdry_penetr[3]
		       << setw(OWID) << bdry_penetr[4]
		       << setw(OWID) << bdry_penetr[5]
		       << setw(OWID) << bdry_penetr[6]
		       << setw(OWID) << bdry_cntnum[1]
		       << setw(OWID) << bdry_cntnum[2]
		       << setw(OWID) << bdry_cntnum[3]
		       << setw(OWID) << bdry_cntnum[4]
		       << setw(OWID) << bdry_cntnum[5]
		       << setw(OWID) << bdry_cntnum[6]
		       << endl;
	}

	// 9. loop break condition: through displacement control mechanism
	
    } while (++iteration < totalSteps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    debugInf.close();
}


// The specimen has been deposited with gravitation within boundaries composed of particles.
// A rectangular pile is then drived into the particles using displacement control.
void Assembly::rectPile_Disp(int   totalSteps,  
			     int   snapNum, 
			     int   interval,
			     const char* iniptclfile,  
			     const char* inibdryfile,
			     const char* ParticleFile, 
			     const char* boundaryfile,
			     const char* contactfile,  
			     const char* progressfile,
			     const char* debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "pile penetrate..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average       translational    rotational       "
	        << "kinetic        potential        total           sample           sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample"
	        << "          sample          sample          sample" << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "         omga            force           moment         energy           energy          "
	        << "energy         energy          energy          density         "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
	        << "epsilon_v" << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);
    debugInf << " iteration    end_bearing     side_friction   total_force" << endl;

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    int    stepsnum=0;
    char   stepsstr[4];
    char   stepsfp[50];
    REAL avgNormal=0;
    REAL avgTangt=0;
    
    int pile[2]={11,12}; // top and down boundaries
    UPDATECTL pilectl[2];

    // iterations start here...
    iteration=0;
    do
    {
      // 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();
	
	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	boundaryForce();

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation

	// displacement control of the pile
	pilectl[0].tran=Vec(0,0,-TIMESTEP*PILE_RATE);
	pilectl[1].tran=Vec(0,0,-TIMESTEP*PILE_RATE);

	updateBoundary(pile, pilectl, 2); 
	updateRectPile();
	if (iteration % interval == 0) {
	    REAL  f7=getShearForce( 7).getz();
	    REAL  f8=getShearForce( 8).getz();
	    REAL  f9=getShearForce( 9).getz();
	    REAL f10=getShearForce(10).getz();
	    REAL  fn=getNormalForce(12).getz();
	    debugInf << setw(OWID) << iteration
		       << setw(OWID) << fn
		       << setw(OWID) << (f7+f8+f9+f10)
		       << setw(OWID) << (fn+f7+f8+f9+f10)
		       << endl;
	}

	// 7. (1) output particles and contacts information
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    printRectPile(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output statistics info.
	if (iteration % interval == 0) {
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3) << endl;
	}

	// 8. loop break condition
	
    } while (++iteration < totalSteps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);
    printRectPile(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    debugInf.close();
}


// The specimen has been deposited with gravitation within boundaries composed of particles.
// An ellipsoidal pile is then drived into the particles using displacement control.
void Assembly::ellipPile_Disp(int   totalSteps,  
			      int   snapNum, 
			      int   interval,
			      REAL dimn,
			      REAL rsize,
			      const char* iniptclfile,
			      const char* ParticleFile, 
			      const char* contactfile,  
			      const char* progressfile,
			      const char* debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "pile penetrate..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average       translational    rotational       "
	        << "kinetic        potential         total           void            sample       coordination"
	        << "       sample           sample          sample          sample          sample          sample"
	        << "          sample          sample          sample         sample           sample         "
	        << " sample          sample          sample          sample          sample" << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "         omga            force           moment         energy           energy          "
	        << "energy         energy            energy          ratio          porosity         number       "
	        << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
	        << "epsilon_v" << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 

    // pre_3. define variables used in iterations
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    REAL void_ratio=0;
    
    // iterations start here...
    iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
        findContact();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 5. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - (-dimn/2);
	l24=dimn*rsize;
	l13=dimn*rsize;
	bulkVolume=l13*l24*l56-ellipPilePeneVol();
	void_ratio=bulkVolume/getParticleVolume()-1;

	// 6. (1) output particles and contacts information
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 6. (2) output statistics info.
	if (iteration % interval == 0) {
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
		        << endl;
	    debugInf << setw(OWID) << iteration
		       << setw(OWID) << getTopFreeParticlePosition().getz()
		       << setw(OWID) << ellipPileTipZ()
		       << setw(OWID) << getTopFreeParticlePosition().getz()-ellipPileTipZ()
		       << setw(OWID) << l13*l24*l56
		       << setw(OWID) << ellipPilePeneVol()
		       << setw(OWID) << bulkVolume
		       << endl;
	}

	// 7. loop break condition
	
    } while (++iteration < totalSteps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    debugInf.close();
}


// The specimen has been deposited with gravitation within rigid boundaries.
// An ellipsoidal penetrator is then impacted into the particles with initial velocity.
void Assembly::ellipPile_Impact(int   totalSteps,  
				int   snapNum, 
				int   interval,
				REAL dimn,
				const char* iniptclfile,
				const char* inibdryfile,
				const char* ParticleFile, 
				const char* contactfile,  
				const char* progressfile,
				const char* debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "penetrator impact..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average       translational    rotational       "
	        << "kinetic        potential         total           void            sample       coordination"
	        << "       sample           sample          sample          sample          sample          sample"
	        << "          sample          sample          sample         sample           sample         "
	        << " sample          sample          sample          sample          sample" << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "         omga            force           moment         energy           energy          "
	        << "energy         energy            energy          ratio          porosity         number       "
	        << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
	        << "epsilon_v" << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles
    readBoundary(inibdryfile);   // create boundaries.

    // pre_3. define variables used in iterations
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }
    
    // iterations start here...
    iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);

	// 4. calculate boundary forces/moments and apply them to particles.
	boundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - (-dimn/2);
	l24=dimn;
	l13=dimn;
	bulkVolume=l13*l24*l56-ellipPilePeneVol();
	void_ratio=bulkVolume/getParticleVolume()-1;

	// 7. (1) output particles and contacts information
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output statistics info.
	if (iteration % interval == 0) {
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    debugInf << setw(OWID) << iteration
		       << setw(OWID) << bdry_penetr[1]
		       << setw(OWID) << bdry_penetr[2]
		       << setw(OWID) << bdry_penetr[3]
		       << setw(OWID) << bdry_penetr[4]
		       << setw(OWID) << bdry_penetr[6]
		       << setw(OWID) << bdry_cntnum[1]
		       << setw(OWID) << bdry_cntnum[2]
		       << setw(OWID) << bdry_cntnum[3]
		       << setw(OWID) << bdry_cntnum[4]
		       << setw(OWID) << bdry_cntnum[6]
		       << endl;

	}

	// 8. loop break condition
	
    } while (++iteration < totalSteps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    debugInf.close();
}


// The specimen has been deposited with gravitation within particle boundaries.
// An ellipsoidal penetrator is then impacted into the particles with initial velocity.
void Assembly::ellipPile_Impact_p(int   totalSteps,  
				  int   snapNum, 
				  int   interval,
				  REAL dimn,
				  const char* iniptclfile,
				  const char* ParticleFile, 
				  const char* contactfile,  
				  const char* progressfile,
				  const char* debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "penetrator impact..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average       translational    rotational       "
	        << "kinetic        potential         total           void            sample       coordination"
	        << "       sample           sample          sample          sample          sample          sample"
	        << "          sample          sample          sample         sample           sample         "
	        << " sample          sample          sample          sample          sample" << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "         omga            force           moment         energy           energy          "
	        << "energy         energy            energy          ratio          porosity         number       "
	        << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
	        << "epsilon_v" << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles

    // pre_3. define variables used in iterations
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    REAL void_ratio=0;
    
    // iterations start here...
    iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 5. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - (-dimn/2);
	l24=dimn;
	l13=dimn;
	bulkVolume=l13*l24*l56-ellipPilePeneVol();
	void_ratio=bulkVolume/getParticleVolume()-1;

	// 6. (1) output particles and contacts information
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 6. (2) output statistics info.
	if (iteration % interval == 0) {
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
		        << endl;
	    debugInf << setw(OWID) << iteration
		       << setw(OWID) << getTopFreeParticlePosition().getz()
		       << setw(OWID) << ellipPileTipZ()
		       << setw(OWID) << getTopFreeParticlePosition().getz()-ellipPileTipZ()
		       << setw(OWID) << l13*l24*l56
		       << setw(OWID) << ellipPilePeneVol()
		       << setw(OWID) << bulkVolume
		       << endl;
	}

	// 7. loop break condition
	
    } while (++iteration < totalSteps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    debugInf.close();
}



// The specimen has been deposited with gravitation within boundaries composed of particles.
// An ellipsoidal pile is then drived into the particles using force control.
// Not recommended.
void Assembly::ellipPile_Force(int   totalSteps,  
			       int   snapNum,
			       int   interval,
			       REAL dimn,
			       REAL force,
			       int   division,
			       const char* iniptclfile,
			       const char* ParticleFile, 
			       const char* contactfile,  
			       const char* progressfile,
			       const char* balancedfile,
			       const char* debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "pile penetrate..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average       translational    rotational       "
	        << "kinetic        potential         total           void            sample       coordination"
	        << "       sample           sample          sample          sample          sample          sample"
	        << "          sample          sample          sample         sample           sample         "
	        << " sample          sample          sample          sample          sample" << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "         omga            force           moment         energy           energy          "
	        << "energy         energy            energy          ratio          porosity         number       "
	        << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
	        << "epsilon_v" << endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout << "stream error!" << endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "pile penetrate..." << endl
	        << "   iteration   apply_force    pile_tip_pos     pile_force" << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 

    // pre_3. define variables used in iterations
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    REAL void_ratio=0;

    REAL zforce_inc=force/division;
    REAL zforce=zforce_inc;

    // iterations start here...
    iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();

	// 5. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - (-dimn/2);
	l24=dimn;
	l13=dimn;
	bulkVolume=l13*l24*l56-ellipPilePeneVol();
	void_ratio=bulkVolume/getParticleVolume()-1;
	
	// 6. update pile external force and position
	if(zforce>ellipPileForce())
	    ellipPileUpdate();

	if(fabs(ellipPileForce()-zforce)/zforce < STRESS_ERROR ){
	    balancedinf << setw(OWID) << iteration
		        << setw(OWID) << zforce
		        << setw(OWID) << getTopFreeParticlePosition().getz()-ellipPileTipZ()
		        << setw(OWID) << ellipPileForce()
		        << endl;
	    zforce += zforce_inc;
	}

	if( iteration % interval == 0){
	    debugInf << setw(OWID) << iteration
		       << setw(OWID) << zforce
		       << setw(OWID) << getTopFreeParticlePosition().getz()-ellipPileTipZ()
		       << setw(OWID) << ellipPileForce()
		       << endl;
	}

	// 7. (1) output particles and contacts information
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output statistics info.
	if (iteration % interval == 0) {
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
		        << endl;
	}

	// 8. loop break condition
	if (fabs((zforce-force)/force)<0.001)
	    break;
	
    } while (++iteration < totalSteps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    debugInf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_a. This function
// performs true triaxial test. Force boundaries are used.
void Assembly::truetriaxial(int   totalSteps,  
			    int   snapNum, 
			    int   interval,
			    REAL sigma_a,     
			    REAL sigma_w,
			    REAL sigma_l,     
			    REAL sigma_h,   
			    int   sigma_division,
			    const char* iniptclfile,  
			    const char* inibdryfile,
			    const char* ParticleFile, 
			    const char* boundaryfile,
			    const char* contactfile,  
			    const char* progressfile,
			    const char* balancedfile, 
			    const char* debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout << "stream error!" << endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "true triaxial..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout << "stream error!" << endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "true triaxial..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;

    debugInf.open(debugfile);
    if(!debugInf) { cout << "stream error!" << endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int mid[2]={1,3};    // boundary 1 and 3
    int max[2]={2,4};    // boundary 2 and 4
    int min[2]={5,6};    // boundary 5 and 6
    UPDATECTL midctl[2];
    UPDATECTL maxctl[2];
    UPDATECTL minctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    REAL sigma_w1=sigma_a;
    REAL sigma_l1=sigma_a;
    REAL sigma_h1=sigma_a;
    REAL sigma_w_inc=(sigma_w-sigma_a)/sigma_division;
    REAL sigma_l_inc=(sigma_l-sigma_a)/sigma_division;
    REAL sigma_h_inc=(sigma_h-sigma_a)/sigma_division;

    // iterations start here...
    iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	boundaryForce();
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    bulkVolume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=bulkVolume/getParticleVolume()-1;

	if (sigma3_1<sigma_h1)
	    minctl[0].tran=Vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=Vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma_h1)
	    minctl[1].tran=Vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=Vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	if (sigma2_1<sigma_l1)
	    midctl[0].tran=Vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=Vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma_l1)
	    midctl[1].tran=Vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=Vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma_w1)
	    maxctl[0].tran=Vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=Vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma_w1)
	    maxctl[1].tran=Vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=Vec(0,-TIMESTEP*RELEASE_RATE,0);
	
	updateBoundary(min,minctl,2);
	updateBoundary(mid,midctl,2);
	updateBoundary(max,maxctl,2);
	updateBoundary6();
	
	// 7. (1) output particles and contacts information
	if (iteration % (totalSteps/snapNum) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (iteration % interval == 0 ){
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()   
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    debugInf << setw(OWID) << iteration
		       << setw(OWID) << bdry_penetr[1]
		       << setw(OWID) << bdry_penetr[2]
		       << setw(OWID) << bdry_penetr[3]
		       << setw(OWID) << bdry_penetr[4]
		       << setw(OWID) << bdry_penetr[5]
		       << setw(OWID) << bdry_penetr[6]
		       << setw(OWID) << bdry_cntnum[1]
		       << setw(OWID) << bdry_cntnum[2]
		       << setw(OWID) << bdry_cntnum[3]
		       << setw(OWID) << bdry_cntnum[4]
		       << setw(OWID) << bdry_cntnum[5]
		       << setw(OWID) << bdry_cntnum[6]
		       << endl;
	}

	// 8. find the balanced status and increase confining pressure
	if (   fabs(sigma1_1-sigma_w1)/sigma_w1 < STRESS_ERROR && fabs(sigma1_2-sigma_w1)/sigma_w1 < STRESS_ERROR
	    && fabs(sigma2_1-sigma_l1)/sigma_l1 < STRESS_ERROR && fabs(sigma2_2-sigma_l1)/sigma_l1 < STRESS_ERROR
	    && fabs(sigma3_1-sigma_h1)/sigma_h1 < STRESS_ERROR && fabs(sigma3_2-sigma_h1)/sigma_h1 < STRESS_ERROR ) {
	    balancedinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()    
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    sigma_w1 += sigma_w_inc;
	    sigma_l1 += sigma_l_inc;
	    sigma_h1 += sigma_h_inc;
	}

	// 9. loop break condition
	if (   fabs(sigma1_1-sigma_w)/sigma_w < STRESS_ERROR && fabs(sigma1_2-sigma_w)/sigma_w < STRESS_ERROR
	    && fabs(sigma2_1-sigma_l)/sigma_l < STRESS_ERROR && fabs(sigma2_2-sigma_l)/sigma_l < STRESS_ERROR
	    && fabs(sigma3_1-sigma_h)/sigma_h < STRESS_ERROR && fabs(sigma3_2-sigma_h)/sigma_h < STRESS_ERROR ) {
	    progressinf << setw(OWID) << iteration
		        << setw(OWID) << getPossContactNum()
		        << setw(OWID) << getActualContactNum()
		        << setw(OWID) << getAvgPenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAvgVelocity() 
		        << setw(OWID) << getAvgOmga()
		        << setw(OWID) << getAvgForce()    
		        << setw(OWID) << getAvgMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAvgPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << bulkVolume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << endl;
	    break;
	}
	
    } while (++iteration < totalSteps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    debugInf.close();
}

*/ 
// rule out

void Assembly::
buildBoundary(int boundaryNum,
	      const char* boundaryFile)
{
  std::ofstream ofs(boundaryFile);
  if(!ofs) { cout << "stream error!" << endl; exit(-1);}
  ofs.setf(std::ios::scientific, std::ios::floatfield);

  Vec  v1 = allContainer.getMinCorner();
  Vec  v2 = allContainer.getMaxCorner();
  Vec  v0 = allContainer.getCenter();
  REAL x1 = v1.getx();
  REAL y1 = v1.gety();
  REAL z1 = v1.getz();
  REAL x2 = v2.getx();
  REAL y2 = v2.gety();
  REAL z2 = v2.getz();
  REAL x0 = v0.getx();
  REAL y0 = v0.gety();
  REAL z0 = v0.getz();

  ofs << setw(OWID) << x1 << setw(OWID) << y1 << setw(OWID) << z1
      << setw(OWID) << x2 << setw(OWID) << y2 << setw(OWID) << z2 << endl
      << setw(OWID) << boundaryNum << endl << endl;
  
  if (boundaryNum == 1){   // only a bottom boundary
    ofs << setw(OWID) << 1 << endl
        << setw(OWID) << 6
        << setw(OWID) << 1 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 0
        << setw(OWID) << -1
        << setw(OWID) << 0
        << setw(OWID) << 0
        << setw(OWID) << z1
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl;
    
  }
  else if (boundaryNum == 5){ // no top boundary
    // boundary 1
    ofs << setw(OWID) << 1 << endl
        << setw(OWID) << 1
        << setw(OWID) << 4 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 0 
        << setw(OWID) << x2
        << setw(OWID) << y0
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << -1
        << setw(OWID) << 0 
        << setw(OWID) << x0
        << setw(OWID) << y1
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << x0
        << setw(OWID) << y2
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 0
        << setw(OWID) << -1
        << setw(OWID) << x0     
        << setw(OWID) << y0
        << setw(OWID) << z1
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl << endl
      
      // boundary 2
        << setw(OWID) << 1 << endl
        << setw(OWID) << 2
        << setw(OWID) << 4 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << x0     
        << setw(OWID) << y2
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << 0 
        << setw(OWID) << x2
        << setw(OWID) << y0      
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << -1
        << setw(OWID) << 0
        << setw(OWID) << 0 
        << setw(OWID) << x1
        << setw(OWID) << y0     
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << 0
        << setw(OWID) << -1
        << setw(OWID) << x0      
        << setw(OWID) << y0     
        << setw(OWID) << z1
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl << endl
      
      // boundary 3
        << setw(OWID) << 1 << endl
        << setw(OWID) << 3
        << setw(OWID) << 4 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << -1
        << setw(OWID) << 0
        << setw(OWID) << 0 
        << setw(OWID) << x1
        << setw(OWID) << y0
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << -1
        << setw(OWID) << 0 
        << setw(OWID) << x0     
        << setw(OWID) << y1
        << setw(OWID) << y0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << x0      
        << setw(OWID) << y2
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << 0
        << setw(OWID) << -1
        << setw(OWID) << x0      
        << setw(OWID) << y0     
        << setw(OWID) << z1
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl << endl
      
      // boundary 4
        << setw(OWID) << 1 << endl
        << setw(OWID) << 4
        << setw(OWID) << 4 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << -1
        << setw(OWID) << 0 
        << setw(OWID) << x0      
        << setw(OWID) << y1
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << 0 
        << setw(OWID) << x2
        << setw(OWID) << y0
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << -1
        << setw(OWID) << 0
        << setw(OWID) << 0 
        << setw(OWID) << x1
        << setw(OWID) << y0
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << 0
        << setw(OWID) << -1
        << setw(OWID) << x0      
        << setw(OWID) << y0     
        << setw(OWID) << z1
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl << endl
      
      // boundary 6
        << setw(OWID) << 1 << endl
        << setw(OWID) << 6
        << setw(OWID) << 5 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 0
        << setw(OWID) << -1
        << setw(OWID) << x0
        << setw(OWID) << y0
        << setw(OWID) << z1
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 0
        << setw(OWID) << x2
        << setw(OWID) << y0
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << -1
        << setw(OWID) << 0
        << setw(OWID) << 0
        << setw(OWID) << x1 
        << setw(OWID) << y0
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << x0      
        << setw(OWID) << y2
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << -1
        << setw(OWID) << 0
        << setw(OWID) << x0      
        << setw(OWID) << y1
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl;
  }
  else if (boundaryNum == 6){ // all 6 boundaries
    // boundary 1
    ofs << setw(OWID) << 1 << endl
        << setw(OWID) << 1
        << setw(OWID) << 5 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 0     
        << setw(OWID) << x2
        << setw(OWID) << y0
        << setw(OWID) << z0     
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << -1
        << setw(OWID) << 0
        << setw(OWID) << x0    
        << setw(OWID) << y1
        << setw(OWID) << z0     
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << x0      
        << setw(OWID) << y2
        << setw(OWID) << z0     
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 0
        << setw(OWID) << 1
        << setw(OWID) << x0     
        << setw(OWID) << y0    
        << setw(OWID) << z2 
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 0 
        << setw(OWID) << -1
        << setw(OWID) << x0     
        << setw(OWID) << y0     
        << setw(OWID) << z1
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl << endl
      
      // boundary 2
        << setw(OWID) << 1 << endl
        << setw(OWID) << 2
        << setw(OWID) << 5 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 1
        << setw(OWID) << 0     
        << setw(OWID) << x0    
        << setw(OWID) << y2
        << setw(OWID) << z0     
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << 0
        << setw(OWID) << x2
        << setw(OWID) << y0
        << setw(OWID) << z0     
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << -1
        << setw(OWID) << 0
        << setw(OWID) << 0
        << setw(OWID) << x1 
        << setw(OWID) << y0
        << setw(OWID) << z0     
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 0
        << setw(OWID) << 1
        << setw(OWID) << x0     
        << setw(OWID) << y0    
        << setw(OWID) << z2 
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 0 
        << setw(OWID) << -1
        << setw(OWID) << x0     
        << setw(OWID) << y0      
        << setw(OWID) << z1
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl << endl
      
      // boundary 3
        << setw(OWID) << 1 << endl
        << setw(OWID) << 3
        << setw(OWID) << 5 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << -1
        << setw(OWID) << 0
        << setw(OWID) << 0     
        << setw(OWID) << x1
        << setw(OWID) << y0
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << -1
        << setw(OWID) << 0
        << setw(OWID) << x0     
        << setw(OWID) << y1
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0  
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << x0       
        << setw(OWID) << y2
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 0
        << setw(OWID) << 1
        << setw(OWID) << x0      
        << setw(OWID) << y0     
        << setw(OWID) << z2 
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 0 
        << setw(OWID) << -1
        << setw(OWID) << x0      
        << setw(OWID) << y0      
        << setw(OWID) << z1
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl << endl
      
      // boundary 4
        << setw(OWID) << 1 << endl
        << setw(OWID) << 4
        << setw(OWID) << 5 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << -1
        << setw(OWID) << 0     
        << setw(OWID) << x0      
        << setw(OWID) << y1
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << 0
        << setw(OWID) << x2
        << setw(OWID) << y0
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << -1 
        << setw(OWID) << 0
        << setw(OWID) << 0
        << setw(OWID) << x1 
        << setw(OWID) << y0
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 0
        << setw(OWID) << 1
        << setw(OWID) << x0      
        << setw(OWID) << y0     
        << setw(OWID) << z2 
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 0 
        << setw(OWID) << -1
        << setw(OWID) << x0      
        << setw(OWID) << y0      
        << setw(OWID) << z1
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl << endl
      
      // boundary 5
        << setw(OWID) << 1 << endl
        << setw(OWID) << 5
        << setw(OWID) << 5 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << 0
        << setw(OWID) << 1     
        << setw(OWID) << x0      
        << setw(OWID) << y0
        << setw(OWID) << z2 
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << 0
        << setw(OWID) << x2
        << setw(OWID) << y0
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << -1 
        << setw(OWID) << 0
        << setw(OWID) << 0
        << setw(OWID) << x1 
        << setw(OWID) << y0
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << x0      
        << setw(OWID) << y2
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << -1
        << setw(OWID) << 0
        << setw(OWID) << x0      
        << setw(OWID) << y1
        << setw(OWID) << z0
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl << endl
      
      // boundary 6
        << setw(OWID) << 1 << endl
        << setw(OWID) << 6
        << setw(OWID) << 5 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << 0
        << setw(OWID) << -1    
        << setw(OWID) << x0      
        << setw(OWID) << y0
        << setw(OWID) << z1
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 1
        << setw(OWID) << 0 
        << setw(OWID) << 0
        << setw(OWID) << x2
        << setw(OWID) << y0
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << -1 
        << setw(OWID) << 0
        << setw(OWID) << 0
        << setw(OWID) << x1 
        << setw(OWID) << y0
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << x0      
        << setw(OWID) << y2
        << setw(OWID) << z0      
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl
      
        << setw(OWID) << 1
        << setw(OWID) << 0
        << setw(OWID) << -1
        << setw(OWID) << 0 
        << setw(OWID) << x0      
        << setw(OWID) << y1
        << setw(OWID) << z0
        << setw(OWID) << 0
        << setw(OWID) << 0 << endl << endl;
  }
  
  ofs.close();
}

// boundaryNum = 6 by default
void Assembly::buildBoundary(const char* boundaryFile)
{
  std::ofstream ofs(boundaryFile);
  if(!ofs) { cout << "stream error!" << endl; exit(-1);}
  ofs.setf(std::ios::scientific, std::ios::floatfield);

  int boundaryNum = 6;

  Vec  v1 = allContainer.getMinCorner();
  Vec  v2 = allContainer.getMaxCorner();
  Vec  v0 = allContainer.getCenter();
  REAL x1 = v1.getx();
  REAL y1 = v1.gety();
  REAL z1 = v1.getz();
  REAL x2 = v2.getx();
  REAL y2 = v2.gety();
  REAL z2 = v2.getz();
  REAL x0 = v0.getx();
  REAL y0 = v0.gety();
  REAL z0 = v0.getz();
  
  ofs << setw(OWID) << x1 << setw(OWID) << y1 << setw(OWID) << z1
      << setw(OWID) << x2 << setw(OWID) << y2 << setw(OWID) << z2 << endl
      << setw(OWID) << boundaryNum << endl << endl;

  // boundary 1
  ofs << setw(OWID) << 1 << endl
      << setw(OWID) << 1
      << setw(OWID) << 5 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0     
      << setw(OWID) << x2
      << setw(OWID) << y0
      << setw(OWID) << z0     
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << -1
      << setw(OWID) << 0
      << setw(OWID) << x0    
      << setw(OWID) << y1
      << setw(OWID) << z0     
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0 
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << x0      
      << setw(OWID) << y2
      << setw(OWID) << z0     
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << 1
      << setw(OWID) << x0     
      << setw(OWID) << y0    
      << setw(OWID) << z2 
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0 
      << setw(OWID) << -1
      << setw(OWID) << x0     
      << setw(OWID) << y0     
      << setw(OWID) << z1
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl << endl
    
    // boundary 2
      << setw(OWID) << 1 << endl
      << setw(OWID) << 2
      << setw(OWID) << 5 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 1
      << setw(OWID) << 0     
      << setw(OWID) << x0    
      << setw(OWID) << y2
      << setw(OWID) << z0     
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 1
      << setw(OWID) << 0 
      << setw(OWID) << 0
      << setw(OWID) << x2
      << setw(OWID) << y0
      << setw(OWID) << z0     
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << -1
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << x1 
      << setw(OWID) << y0
      << setw(OWID) << z0     
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << 1
      << setw(OWID) << x0     
      << setw(OWID) << y0    
      << setw(OWID) << z2 
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0 
      << setw(OWID) << -1
      << setw(OWID) << x0     
      << setw(OWID) << y0      
      << setw(OWID) << z1
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl << endl
    
    // boundary 3
      << setw(OWID) << 1 << endl
      << setw(OWID) << 3
      << setw(OWID) << 5 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << -1
      << setw(OWID) << 0
      << setw(OWID) << 0     
      << setw(OWID) << x1
      << setw(OWID) << y0
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << -1
      << setw(OWID) << 0
      << setw(OWID) << x0     
      << setw(OWID) << y1
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0  
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << x0       
      << setw(OWID) << y2
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << 1
      << setw(OWID) << x0      
      << setw(OWID) << y0     
      << setw(OWID) << z2 
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0 
      << setw(OWID) << -1
      << setw(OWID) << x0      
      << setw(OWID) << y0      
      << setw(OWID) << z1
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl << endl
    
    // boundary 4
      << setw(OWID) << 1 << endl
      << setw(OWID) << 4
      << setw(OWID) << 5 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0 
      << setw(OWID) << -1
      << setw(OWID) << 0     
      << setw(OWID) << x0      
      << setw(OWID) << y1
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 1
      << setw(OWID) << 0 
      << setw(OWID) << 0
      << setw(OWID) << x2
      << setw(OWID) << y0
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << -1 
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << x1 
      << setw(OWID) << y0
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << 1
      << setw(OWID) << x0      
      << setw(OWID) << y0     
      << setw(OWID) << z2 
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0 
      << setw(OWID) << -1
      << setw(OWID) << x0      
      << setw(OWID) << y0      
      << setw(OWID) << z1
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl << endl
    
    // boundary 5
      << setw(OWID) << 1 << endl
      << setw(OWID) << 5
      << setw(OWID) << 5 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0 
      << setw(OWID) << 0
      << setw(OWID) << 1     
      << setw(OWID) << x0      
      << setw(OWID) << y0
      << setw(OWID) << z2 
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 1
      << setw(OWID) << 0 
      << setw(OWID) << 0
      << setw(OWID) << x2
      << setw(OWID) << y0
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << -1 
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << x1 
      << setw(OWID) << y0
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << x0      
      << setw(OWID) << y2
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << -1
      << setw(OWID) << 0
      << setw(OWID) << x0      
      << setw(OWID) << y1
      << setw(OWID) << z0
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl << endl
    
    // boundary 6
      << setw(OWID) << 1 << endl
      << setw(OWID) << 6
      << setw(OWID) << 5 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0 
      << setw(OWID) << 0
      << setw(OWID) << -1    
      << setw(OWID) << x0      
      << setw(OWID) << y0
      << setw(OWID) << z1
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 1
      << setw(OWID) << 0 
      << setw(OWID) << 0
      << setw(OWID) << x2
      << setw(OWID) << y0
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << -1 
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << x1 
      << setw(OWID) << y0
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << x0      
      << setw(OWID) << y2
      << setw(OWID) << z0      
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl
    
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << -1
      << setw(OWID) << 0 
      << setw(OWID) << x0      
      << setw(OWID) << y1
      << setw(OWID) << z0
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl << endl; 

  ofs.close();
}

} // namespace dem ends


/*
// freeType:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
// ht- how many times of size would be the floating height
void Assembly::generateParticle(gradation& grad,
			const char* ParticleFile,
			int freeType,
			REAL ht)
{
    REAL x,y,z;
    Particle* newptcl;
    particleNum = 0;
    REAL est =1.02;
    int grid=9;  
    // grid: dimension of free particle array.
    // 7 - small dimn container
    // 9 - medium dimn container 
    // 11- large dimn container 

    REAL dimn=grad.dimn;
    if (freeType == 0) {      // just one free particle
	newptcl = new Particle(particleNum+1, 0, Vec(dimn/2/40,dimn/2/20,dimn/2), grad, YOUNG, POISSON);
	particleVec.push_back(newptcl);
	particleNum++;
    }
    else if (freeType == 1) { // a horizontal layer of free particles
	z=dimn/2;
	for (x=-dimn/2*(grid-1)/10; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
	    for (y=-dimn/2*(grid-1)/10; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5){
		newptcl = new Particle(particleNum+1, 0, Vec(x,y,z), grad, YOUNG, POISSON);
		particleVec.push_back(newptcl);
		particleNum++;
	    }
    }
    else if (freeType == 2) { // multiple layers of free particles
	REAL offset=0; // 0 for ellipsoids; dimn/2/5/5 for spheres
	if (grad.ratio_ba==1.0 && grad.ratio_ca==1.0)
	    offset = dimn/2/5/5;
	REAL z0 = -dimn/2*9/10 ;// dimn/2;
	for (z=z0; z<z0 + dimn*ht; z+=dimn/2/5) {
	//for (z=-dimn/2*4/5; z<dimn/2 + dimn*ht; z+=dimn/2/10) { // spheres
	    for (x=-dimn/2*(grid-1)/10+offset; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
		for (y=-dimn/2*(grid-1)/10+offset; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5){
		    newptcl = new Particle(particleNum+1, 0, Vec(x,y,z), grad, YOUNG, POISSON);
		    particleVec.push_back(newptcl);
		    particleNum++;
		}	
	    offset *= -1;
	}
    }

    printParticle(ParticleFile);

}
*/
 /*
// create a specimen from discreate particles through floating and then gravitation,
// boundaries are composed of fixed particles.
void Assembly::deposit_PtclBdry(gradation& grad,
				int   freeType,
				REAL rsize,
				int   totalSteps,  
				int   snapNum,
				int   interval,
				const char* iniptclfile,   
				const char* ParticleFile, 
				const char* contactfile,
				const char* progressfile, 
				const char* debugfile)
{
    if (grad.rorc == 1) {
	RORC = grad.rorc;
	container.setCenter(Vec(0,0,0));
	container.setDimx(grad.dimn);
	container.setDimy(grad.dimn);
	container.setDimz(grad.dimn);
	
	generate_p(grad, iniptclfile, freeType, rsize, 4.0);
	deposit_p(totalSteps,        // totalSteps
		  snapNum,          // number of snapNum
		  interval,           // print interval
		  grad.dimn,          // dimension of particle-composed-boundary
		  rsize,              // relative container size
		  iniptclfile,        // input file, initial particles
		  ParticleFile,       // output file, resulted particles, including snapNum 
		  contactfile,        // output file, resulted contacts, including snapNum 
		  progressfile,       // output file, statistical info
		  debugfile);         // output file, debug info
    }
}
 */
  /*
// freeType:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
// ht- how many times of size would be the floating height
void Assembly::generate_p(gradation&  grad,
			 const char* ParticleFile,
			 int freeType,
			 REAL rsize,
			 REAL ht)
{
    REAL x,y,z;
    Particle* newptcl;
    particleNum = 0;
    REAL wall=2.2; // wall - wall height; ht - free particle height
    REAL est =1.02;
    int grid=static_cast<int> (nearbyint(rsize*10)-1);  

    // grid: dimension of free particle array.
    // 7 - small dimn container
    // 9 - medium dimn container 
    // 11- large dimn container 

    REAL dimn=grad.dimn;
    // particle boundary 1
    x=dimn/2*(grid+1)/10;
    for (y=-dimn/2*grid/10; y<dimn/2*grid/10*est; y+=dimn/2/5)
	for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5){
	    newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99, YOUNG, POISSON);
	    particleVec.push_back(newptcl);
	    particleNum++;
	}

    // particle boundary 2
    y=dimn/2*(grid+1)/10;
    for (x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5)
	for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5){
	    newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99, YOUNG, POISSON);
	    particleVec.push_back(newptcl);
	    particleNum++;
	}

    // particle boundary 3
    x=-dimn/2*(grid+1)/10;
    for (y=-dimn/2*grid/10; y<dimn/2*grid/10*est; y+=dimn/2/5)
	for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5){
	    newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99, YOUNG, POISSON);
	    particleVec.push_back(newptcl);
	    particleNum++;
	}

    // particle boundary 4
    y=-dimn/2*(grid+1)/10;
    for (x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5)
	for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5){
	    newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99, YOUNG, POISSON);
	    particleVec.push_back(newptcl);
	    particleNum++;
	}

    // particle boundary 6
    z=-dimn/2;
    for (y=-dimn/2*grid/10; y<dimn/2*grid/10*est; y+=dimn/2/5)
	for( x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5){
	    newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99, YOUNG, POISSON);
	    particleVec.push_back(newptcl);
	    particleNum++;
	}

    if (freeType == 0) {      // just one free particle
	newptcl = new Particle(particleNum+1, 0, Vec(dimn/2/40,dimn/2/20,dimn/2), grad, YOUNG, POISSON);
	particleVec.push_back(newptcl);
	particleNum++;
    }
    else if (freeType == 1) { // a horizontal layer of free particles
	z=dimn/2;
	for (x=-dimn/2*(grid-1)/10; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
	    for (y=-dimn/2*(grid-1)/10; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5){
		newptcl = new Particle(particleNum+1, 0, Vec(x,y,z), grad, YOUNG, POISSON);
		particleVec.push_back(newptcl);
		particleNum++;
	    }
    }
    else if (freeType == 2) { // multiple layers of free particles
	for (z=dimn/2; z<dimn/2 + dimn*ht; z+=dimn/2/5)
	    for (x=-dimn/2*(grid-1)/10; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
		for (y=-dimn/2*(grid-1)/10; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5){
		    newptcl = new Particle(particleNum+1, 0, Vec(x,y,z), grad, YOUNG, POISSON);
		    particleVec.push_back(newptcl);
		    particleNum++;
		}	
    }
    
    printParticle(ParticleFile);
    
}
*/

   /*
void Assembly::boundaryForce(){
  std::vector<BOUNDARY*>::iterator rt;
  for(rt=boundaryVec.begin();rt!=boundaryVec.end();++rt)
    (*rt)->boundaryForce(boundaryTgtMap);

  
  vector<boundarytgt>::iterator it;
  std::vector<BOUNDARY*>::iterator rt;

  for(rt=boundaryVec.begin();rt!=boundaryVec.end();++rt){	
    (*rt)->boundaryForce(boundaryTgtMap);
    for (it=boundaryTgtMap[(*rt)->bdry_id].begin();it!=boundaryTgtMap[(*rt)->bdry_id].end();++it){
      debugInf << setw(OWID) << iteration
		 << setw(OWID) << (*rt)->bdry_id
		 << setw(OWID) << boundaryTgtMap[(*rt)->bdry_id].size()
		 << setw(OWID) << it->TgtForce.getx()
		 << setw(OWID) << it->TgtForce.gety()
		 << setw(OWID) << it->TgtForce.getz()
		 << endl;
      // << setw(OWID) << it->TgtPeak << endl;
    }
  }

  }*/
