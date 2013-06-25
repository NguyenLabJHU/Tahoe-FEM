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
#include "const.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <ctime>
#include <cassert>
#include <utility>
#include <sys/time.h>
#include <omp.h>

//#define BINNING
//#define TIME_PROFILE

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


char *combineString(char *cstr, const char *str, int num, int width) {
  std::string obj(str);
  std::stringstream ss;
  ss << std::setw(width) << std::setfill('0') << std::right << num;
  obj += ss.str();
  return strcpy( cstr, obj.c_str() );
}


void Assembly::setCommunicator(boost::mpi::communicator &comm) {
  boostWorld = comm;
  mpiWorld = MPI_Comm(comm);
  mpiProcX = static_cast<int> (dem::Parameter::getSingleton().parameter["mpiProcX"]);
  mpiProcY = static_cast<int> (dem::Parameter::getSingleton().parameter["mpiProcY"]);
  mpiProcZ = static_cast<int> (dem::Parameter::getSingleton().parameter["mpiProcZ"]);
  
  // create Cartesian virtual topology (unavailable in boost.mpi) 
  int ndim = 3;
  int dims[3] = {mpiProcX, mpiProcY, mpiProcZ};
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
  if (mpiRank == 0) {
    REAL minX = dem::Parameter::getSingleton().parameter["minX"];
    REAL minY = dem::Parameter::getSingleton().parameter["minY"];
    REAL minZ = dem::Parameter::getSingleton().parameter["minZ"];
    REAL maxX = dem::Parameter::getSingleton().parameter["maxX"];
    REAL maxY = dem::Parameter::getSingleton().parameter["maxY"];
    REAL maxZ = dem::Parameter::getSingleton().parameter["maxZ"];
    int particleLayers = dem::Parameter::getSingleton().parameter["particleLayers"];

    setContainer(Rectangle(minX, minY, minZ, maxX, maxY, maxZ));

    buildBoundary(5, "dep_boundary");
    
    int sieveNum = static_cast<int> (dem::Parameter::getSingleton().parameter["sieveNum"]);
    std::vector<REAL> percent(sieveNum), size(sieveNum);
    std::vector<std::pair<REAL, REAL> > &grada = dem::Parameter::getSingleton().gradation;
    assert(grada.size() == sieveNum);
    for (int i = 0; i < sieveNum; ++i) {
      percent[i] = grada[i].first;
      size[i] = grada[i].second;
    }
    REAL ratioBA = dem::Parameter::getSingleton().parameter["ratioBA"];
    REAL ratioCA = dem::Parameter::getSingleton().parameter["ratioCA"];
    setGradation(Gradation(sieveNum, percent, size, ratioBA, ratioCA));

    generateParticle(particleLayers, "float_particle"); 
  }

  deposit(static_cast<int> (dem::Parameter::getSingleton().parameter["totalSteps"]),
	  static_cast<int> (dem::Parameter::getSingleton().parameter["snapNum"]),
	  static_cast<int> (dem::Parameter::getSingleton().parameter["statInterv"]),
	  "dep_boundary",
	  "float_particle");
  
  if (mpiRank == 0) {
    setContainer(Rectangle(allContainer.getMinCorner().getX(),
			   allContainer.getMinCorner().getY(),
			   allContainer.getMinCorner().getZ(),
			   allContainer.getMaxCorner().getX(),
			   allContainer.getMaxCorner().getY(),
			   dem::Parameter::getSingleton().parameter["trimHeight"]));
    buildBoundary("trm_boundary");
    trim(false,
	 "dep_particle_end",
	 "trm_particle");
  }
    
}


void Assembly::expandCavityParticle() 
{
  if (mpiRank == 0) {
    const char *inputParticle = dem::Parameter::getSingleton().datafile["particleFile"].c_str();
    REAL percent = dem::Parameter::getSingleton().parameter["expandPercent"];
    readParticle(inputParticle); 
    
    REAL x1 = dem::Parameter::getSingleton().parameter["cavityMinX"];
    REAL y1 = dem::Parameter::getSingleton().parameter["cavityMinY"];
    REAL z1 = dem::Parameter::getSingleton().parameter["cavityMinZ"];
    REAL x2 = dem::Parameter::getSingleton().parameter["cavityMaxX"];
    REAL y2 = dem::Parameter::getSingleton().parameter["cavityMaxY"];
    REAL z2 = dem::Parameter::getSingleton().parameter["cavityMaxZ"];
    
    std::vector<Particle*> cavityParticleVec;
    std::vector<Particle*>::iterator it;
    Vec center;
    
    for (it = allParticleVec.begin(); it != allParticleVec.end(); ++it ){
      center=(*it)->getCurrPos();
      if(center.getX() > x1 && center.getX() < x2 &&
	 center.getY() > y1 && center.getY() < y2 &&
	 center.getZ() > z1 && center.getZ() < z2 ) {
	(*it)->expand(percent);
	cavityParticleVec.push_back(*it);
      }
    }
    
    printParticle("cav_particle", cavityParticleVec);
    printParticle("exp_particle");
  }
  
  deposit(static_cast<int> (dem::Parameter::getSingleton().parameter["totalSteps"]),
	  static_cast<int> (dem::Parameter::getSingleton().parameter["snapNum"]),
	  static_cast<int> (dem::Parameter::getSingleton().parameter["statInterv"]),
	  dem::Parameter::getSingleton().datafile["boundaryFile"].c_str(),
	  "exp_particle");

}


// particleLayers:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
void Assembly::
generateParticle(int particleLayers,
		 const char *genParticle)
{
  REAL young = dem::Parameter::getSingleton().parameter["young"];
  REAL poisson = dem::Parameter::getSingleton().parameter["poisson"];

  REAL x,y,z;
  Particle* newptcl;
  int particleNum = 0;
  REAL diameter = gradation.getPtclMaxRadius()*2.0;

  REAL offset   = 0;
  REAL edge     = diameter;
  if (gradation.getSize().size() == 1 &&
      gradation.getPtclRatioBA() == 1.0 && 
      gradation.getPtclRatioCA() == 1.0) {
    edge   = diameter*2.0;
    offset = diameter*0.25;
  }
  
  REAL x1 = allContainer.getMinCorner().getX() + edge;
  REAL y1 = allContainer.getMinCorner().getY() + edge;
  REAL z1 = allContainer.getMinCorner().getZ() + diameter;
  REAL x2 = allContainer.getMaxCorner().getX() - edge;
  REAL y2 = allContainer.getMaxCorner().getY() - edge;
  REAL z2 = allContainer.getMaxCorner().getZ() - diameter;
  REAL x0 = allContainer.getCenter().getX();
  REAL y0 = allContainer.getCenter().getY();
  REAL z0 = allContainer.getCenter().getZ();

  if (particleLayers == 0) {      // just one free particle
    newptcl = new Particle(particleNum+1, 0, Vec(x0,y0,z0), gradation, young, poisson);
    allParticleVec.push_back(newptcl);
    particleNum++;
  }
  else if (particleLayers == 1) { // a horizontal layer of free particles
    for (x = x1; x - x2 < EPS; x += diameter)
      for (y = y1; y - y2 < EPS; y += diameter) {
	newptcl = new Particle(particleNum+1, 0, Vec(x,y,z0), gradation, young, poisson);
	allParticleVec.push_back(newptcl);
	particleNum++;
      }
  }
  else if (particleLayers == 2) { // multiple layers of free particles
    for (z = z1; z - z2 < EPS; z += diameter) {
      for (x = x1 + offset; x - x2 < EPS; x += diameter)
	for (y = y1 + offset; y - y2 < EPS; y += diameter) {
	  newptcl = new Particle(particleNum+1, 0, Vec(x,y,z), gradation, young, poisson);
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
     const char *inputParticle,
     const char *trmParticle)
{
  if (toRebuild) readParticle(inputParticle);
  trimHistoryNum = allParticleVec.size();

  Vec  v1 = allContainer.getMinCorner();
  Vec  v2 = allContainer.getMaxCorner();
  REAL x1 = v1.getX();
  REAL y1 = v1.getY();
  REAL z1 = v1.getZ();
  REAL x2 = v2.getX();
  REAL y2 = v2.getY();
  REAL z2 = v2.getZ();
  REAL maxR = gradation.getPtclMaxRadius();
 
  std::vector<Particle*>::iterator itr;
  Vec center;

  for (itr = allParticleVec.begin(); itr != allParticleVec.end(); ) {
    center=(*itr)->getCurrPos();
    if(center.getX() < x1 || center.getX() > x2 ||
       center.getY() < y1 || center.getY() > y2 ||
       center.getZ() < z1 || center.getZ() + maxR > z2)
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
	int statInterv,
	const char *inputBoundary,
	const char *inputParticle) 
{
  if (mpiRank == 0) {
    readBoundary(inputBoundary);
    readParticle(inputParticle); 
  }
  scatterParticle(); // scatter particles only once; also updates grid for the first time

  iteration = 1;
  int iterSnap = static_cast<int> (dem::Parameter::getSingleton().parameter["snapResumeNum"]);
  double time0, time1, time2, commuT, migraT, gatherT, totalT;
  while (iteration <= totalSteps) {
    commuT = migraT = gatherT = totalT = 0;
    time0 = MPI_Wtime();

    time1 = MPI_Wtime();
    commuParticle();
    time2 = MPI_Wtime(); commuT = time2 - time1;

    findContact();
    if (isBdryProcess()) findBdryContact();

    clearContactForce();
    internalForce();
    if (isBdryProcess()) boundaryForce();
    updateParticle();
    updateGridMaxZ();

    if (iteration % (totalSteps / snapNum) == 0) {
      ++iterSnap;
      char cstr[50];

      time1 = MPI_Wtime();
      gatherParticle();
      time2 = MPI_Wtime(); gatherT = time2 - time1;
      if (mpiRank == 0) {
	plotBoundary(combineString(cstr, "dep_bdryplot_", iterSnap, 3));
	plotGrid(combineString(cstr, "dep_gridplot_", iterSnap, 3));
	printParticle(combineString(cstr, "dep_particle_", iterSnap, 3));
	releaseGatheredParticle();
      }

      printContact(combineString(cstr, "dep_contact_", iterSnap, 3));

    }

    releaseRecvParticle(); // late release because printContact refers to received particles
    time1 = MPI_Wtime();
    migrateParticle();
    time2 = MPI_Wtime(); migraT = time2 - time1;
 
    time2 = MPI_Wtime(); totalT = time2 - time0;
    if (mpiRank == 0 && (iteration + 1) % (totalSteps / snapNum) == 0) // ignore gather and print time
      debugInf << "iter=" << std::setw(8) << iteration << std::setprecision(2)
	       << " commu=" << commuT
	       << " gather=" << gatherT
	       << " migra=" << migraT
	       << " total=" << totalT 
	       << " overhead=" << std::fixed << (commuT + gatherT + migraT)/totalT*100 << '%' 
	       << std::scientific << std::setprecision(6) << std::endl;

    ++iteration;
  } 

}


void Assembly::
findParticleInRectangle(const Rectangle &container,
			const std::vector<Particle*> &inputParticle,
			std::vector<Particle*> &foundParticle) {
  Vec  v1 = container.getMinCorner();
  Vec  v2 = container.getMaxCorner();
  REAL x1 = v1.getX();
  REAL y1 = v1.getY();
  REAL z1 = v1.getZ();
  REAL x2 = v2.getX();
  REAL y2 = v2.getY();
  REAL z2 = v2.getZ();
  for (int pt = 0; pt < inputParticle.size(); ++pt) {
    Vec center = inputParticle[pt]->getCurrPos();
    // it is critical to use EPS
    if (center.getX() - x1 >= -EPS && center.getX() - x2 < -EPS &&
	center.getY() - y1 >= -EPS && center.getY() - y2 < -EPS &&
	center.getZ() - z1 >= -EPS && center.getZ() - z2 < -EPS)
      foundParticle.push_back(inputParticle[pt]);
  }
}


void Assembly::removeParticleOutRectangle() {
  Vec  v1 = container.getMinCorner();
  Vec  v2 = container.getMaxCorner();
  REAL x1 = v1.getX();
  REAL y1 = v1.getY();
  REAL z1 = v1.getZ();
  REAL x2 = v2.getX();
  REAL y2 = v2.getY();
  REAL z2 = v2.getZ();

  std::vector<Particle*>::iterator itr;
  Vec center;
  //int flag = 0;

  for (itr = particleVec.begin(); itr != particleVec.end(); ) {
    center=(*itr)->getCurrPos();
    // it is critical to use EPS
    if ( !(center.getX() - x1 >= -EPS && center.getX() - x2 < -EPS &&
	   center.getY() - y1 >= -EPS && center.getY() - y2 < -EPS &&
	   center.getZ() - z1 >= -EPS && center.getZ() - z2 < -EPS) )
      {
	/*
	std::cout << "iter=" << std::setw(8) << iteration << " rank=" << std::setw(2) << mpiRank
		  << " removed=" << std::setw(3) << (*itr)->getId();	
	flag = 1;
	*/
	delete (*itr); // release memory
	itr = particleVec.erase(itr); 
      }
    else
      ++itr;
  }
  /*
  if (flag == 1) {
    std::cout << " now " << particleVec.size() << ": ";
    for (std::vector<Particle*>::const_iterator it = particleVec.begin(); it != particleVec.end(); ++it)
      std::cout << std::setw(3) << (*it)->getId();
    std::cout << std::endl;
  }
  */

}


REAL Assembly::getPtclMaxX(const std::vector<Particle*> &inputParticle) const {
  std::vector<Particle*>::const_iterator it = inputParticle.begin();
  REAL x0 = (*it)->getCurrPos().getX();
  for (; it != inputParticle.end(); ++it) {
    if ( (*it)->getCurrPos().getX() > x0 )
      x0 = (*it)->getCurrPos().getX();
  }
  return x0;
}


REAL Assembly::getPtclMinX(const std::vector<Particle*> &inputParticle) const {
  std::vector<Particle*>::const_iterator it = inputParticle.begin();
  REAL x0 = (*it)->getCurrPos().getX();
  for (; it != inputParticle.end(); ++it) {
    if ( (*it)->getCurrPos().getX() < x0 )
      x0 = (*it)->getCurrPos().getX();
  }
  return x0;
}


REAL Assembly::getPtclMaxY(const std::vector<Particle*> &inputParticle) const {
  std::vector<Particle*>::const_iterator it = inputParticle.begin();
  REAL y0 = (*it)->getCurrPos().getY();
  for (; it != inputParticle.end(); ++it) {
    if ( (*it)->getCurrPos().getY() > y0 )
      y0 = (*it)->getCurrPos().getY();
  }
  return y0;
}


REAL Assembly::getPtclMinY(const std::vector<Particle*> &inputParticle) const {
  std::vector<Particle*>::const_iterator it = inputParticle.begin();
  REAL y0 = (*it)->getCurrPos().getY();
  for (; it != inputParticle.end(); ++it) {
    if ( (*it)->getCurrPos().getY() < y0 )
      y0 = (*it)->getCurrPos().getY();
  }
  return y0;
}


REAL Assembly::getPtclMaxZ(const std::vector<Particle*> &inputParticle) const {
  std::vector<Particle*>::const_iterator it = inputParticle.begin();
  REAL z0 = (*it)->getCurrPos().getZ();
  for (; it != inputParticle.end(); ++it) {
    if ( (*it)->getCurrPos().getZ() > z0 )
      z0 = (*it)->getCurrPos().getZ();
  }
  return z0;
}


REAL Assembly::getPtclMinZ(const std::vector<Particle*> &inputParticle) const {
  std::vector<Particle*>::const_iterator it = inputParticle.begin();
  REAL z0 = (*it)->getCurrPos().getZ();
  for (; it != inputParticle.end(); ++it) {
    if ( (*it)->getCurrPos().getZ() < z0 )
      z0 = (*it)->getCurrPos().getZ();
  }
  return z0;
}


void Assembly::scatterParticle() {
  // partition particles and send to each process
  if (mpiRank == 0) { // process 0
    
    setGrid(Rectangle(grid.getMinCorner().getX(),
		      grid.getMinCorner().getY(),
		      grid.getMinCorner().getZ(),
		      grid.getMaxCorner().getX(),
		      grid.getMaxCorner().getY(),
		      getPtclMaxZ(allParticleVec) + gradation.getPtclMaxRadius() ));
    
    Vec v1 = grid.getMinCorner();
    Vec v2 = grid.getMaxCorner();
    Vec vspan = v2 - v1;

    boost::mpi::request *reqs = new boost::mpi::request [mpiSize - 1];
    std::vector<Particle*> tmpParticleVec;
    for (int iRank = mpiSize - 1; iRank >= 0; --iRank) {
      tmpParticleVec.clear(); // do not release memory!
      int ndim = 3;
      int coords[3];
      MPI_Cart_coords(cartComm, iRank, ndim, coords);
      Rectangle container(v1.getX() + vspan.getX() / mpiProcX * coords[0],
			  v1.getY() + vspan.getY() / mpiProcY * coords[1],
			  v1.getZ() + vspan.getZ() / mpiProcZ * coords[2],
			  v1.getX() + vspan.getX() / mpiProcX * (coords[0] + 1),
			  v1.getY() + vspan.getY() / mpiProcY * (coords[1] + 1),
			  v1.getZ() + vspan.getZ() / mpiProcZ * (coords[2] + 1));
      findParticleInRectangle(container, allParticleVec, tmpParticleVec);
      if (iRank != 0)
	reqs[iRank - 1] = boostWorld.isend(iRank, mpiTag, tmpParticleVec); // non-blocking send
      if (iRank == 0) {
	particleVec.resize(tmpParticleVec.size());
	for (int i = 0; i < particleVec.size(); ++i)
	  particleVec[i] = new Particle(*tmpParticleVec[i]); // default synthesized copy constructor
      } // now particleVec do not share memeory with allParticleVec
    }
    boost::mpi::wait_all(reqs, reqs + mpiSize - 1); // for non-blocking send
    delete [] reqs;

  } else { // other processes except 0
    boostWorld.recv(0, mpiTag, particleVec);
  }

  // content of allParticleVec is not needed any more 
  // but the variable is used for later gathering to print, so clear it.
  if (mpiRank == 0) // process 0
    releaseGatheredParticle();

  // broadcast necessary info
  broadcast(boostWorld, gradation, 0);
  broadcast(boostWorld, boundaryVec, 0);
  broadcast(boostWorld, allContainer, 0);
  broadcast(boostWorld, grid, 0);
}


bool Assembly::isBdryProcess() {
  return (mpiCoords[0] == 0 || mpiCoords[0] == mpiProcX - 1 ||
	  mpiCoords[1] == 0 || mpiCoords[1] == mpiProcY - 1 ||
	  mpiCoords[2] == 0 || mpiCoords[2] == mpiProcZ - 1);
}


void Assembly::commuParticle() 
{
  // determine container of each process
  Vec v1 = grid.getMinCorner();
  Vec v2 = grid.getMaxCorner();
  Vec vspan = v2 - v1;
  container = Rectangle(v1.getX() + vspan.getX() / mpiProcX * mpiCoords[0],
			v1.getY() + vspan.getY() / mpiProcY * mpiCoords[1],
			v1.getZ() + vspan.getZ() / mpiProcZ * mpiCoords[2],
			v1.getX() + vspan.getX() / mpiProcX * (mpiCoords[0] + 1),
			v1.getY() + vspan.getY() / mpiProcY * (mpiCoords[1] + 1),
			v1.getZ() + vspan.getZ() / mpiProcZ * (mpiCoords[2] + 1));

  // find neighboring blocks
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
  v1 = container.getMinCorner(); // redefine v1, v2 in terms of process
  v2 = container.getMaxCorner();   
  //std::cout << "rank=" << mpiRank << ' ' << v1.getX() << ' ' << v1.getY() << ' ' << v1.getZ() << ' '  << v2.getX() << ' ' << v2.getY() << ' ' << v2.getZ() << std::endl;
  REAL cellSize = gradation.getPtclMaxRadius() * 2;
  // 6 surfaces
  if (rankX1 >= 0) { // surface x1
    Rectangle containerX1(v1.getX(), v1.getY(), v1.getZ(), 
			  v1.getX() + cellSize, v2.getY(), v2.getZ());
    findParticleInRectangle(containerX1, particleVec, particleX1);
    reqX1[0] = boostWorld.isend(rankX1, mpiTag,  particleX1);
    reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rParticleX1);
  }
  if (rankX2 >= 0) { // surface x2
    Rectangle containerX2(v2.getX() - cellSize, v1.getY(), v1.getZ(),
			  v2.getX(), v2.getY(), v2.getZ());
    findParticleInRectangle(containerX2, particleVec, particleX2);
    reqX2[0] = boostWorld.isend(rankX2, mpiTag,  particleX2);
    reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rParticleX2);
  }
  if (rankY1 >= 0) {  // surface y1
    Rectangle containerY1(v1.getX(), v1.getY(), v1.getZ(), 
			  v2.getX(), v1.getY() + cellSize, v2.getZ());
    findParticleInRectangle(containerY1, particleVec, particleY1);
    reqY1[0] = boostWorld.isend(rankY1, mpiTag,  particleY1);
    reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rParticleY1);
  }
  if (rankY2 >= 0) {  // surface y2
    Rectangle containerY2(v1.getX(), v2.getY() - cellSize, v1.getZ(),
			  v2.getX(), v2.getY(), v2.getZ());
    findParticleInRectangle(containerY2, particleVec, particleY2);
    reqY2[0] = boostWorld.isend(rankY2, mpiTag,  particleY2);
    reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rParticleY2);
  }
  if (rankZ1 >= 0) {  // surface z1
    Rectangle containerZ1(v1.getX(), v1.getY(), v1.getZ(),
			  v2.getX(), v2.getY(), v1.getZ() + cellSize);
    findParticleInRectangle(containerZ1, particleVec, particleZ1);
    reqZ1[0] = boostWorld.isend(rankZ1, mpiTag,  particleZ1);
    reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rParticleZ1);
  }
  if (rankZ2 >= 0) {  // surface z2
    Rectangle containerZ2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
			  v2.getX(), v2.getY(), v2.getZ());
    findParticleInRectangle(containerZ2, particleVec, particleZ2);
    reqZ2[0] = boostWorld.isend(rankZ2, mpiTag,  particleZ2);
    reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rParticleZ2);
  }
  // 12 edges
  if (rankX1Y1 >= 0) { // edge x1y1
    Rectangle containerX1Y1(v1.getX(), v1.getY(), v1.getZ(),
			    v1.getX() + cellSize, v1.getY() + cellSize, v2.getZ());
    findParticleInRectangle(containerX1Y1, particleVec, particleX1Y1);
    reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag,  particleX1Y1);
    reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rParticleX1Y1);
  }
  if (rankX1Y2 >= 0) { // edge x1y2
    Rectangle containerX1Y2(v1.getX(), v2.getY() - cellSize, v1.getZ(),
			    v1.getX() + cellSize, v2.getY(), v2.getZ());
    findParticleInRectangle(containerX1Y2, particleVec, particleX1Y2);
    reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag,  particleX1Y2);
    reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rParticleX1Y2);
  }
  if (rankX1Z1 >= 0) { // edge x1z1
    Rectangle containerX1Z1(v1.getX(), v1.getY(), v1.getZ(),
			    v1.getX() + cellSize, v2.getY(), v1.getZ() + cellSize);
    findParticleInRectangle(containerX1Z1, particleVec, particleX1Z1);
    reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag,  particleX1Z1);
    reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rParticleX1Z1);
  }
  if (rankX1Z2 >= 0) { // edge x1z2
    Rectangle containerX1Z2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
			    v1.getX() + cellSize, v2.getY(), v2.getZ());
    findParticleInRectangle(containerX1Z2, particleVec, particleX1Z2);
    reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag,  particleX1Z2);
    reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rParticleX1Z2);
  }
  if (rankX2Y1 >= 0) { // edge x2y1
    Rectangle containerX2Y1(v2.getX() - cellSize, v1.getY(), v1.getZ(),
			    v2.getX(), v1.getY() + cellSize, v2.getZ());
    findParticleInRectangle(containerX2Y1, particleVec, particleX2Y1);
    reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag,  particleX2Y1);
    reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rParticleX2Y1);
  }
  if (rankX2Y2 >= 0) { // edge x2y2
    Rectangle containerX2Y2(v2.getX() - cellSize, v2.getY() - cellSize, v1.getZ(),
			    v2.getX(), v2.getY(), v2.getZ());
    findParticleInRectangle(containerX2Y2, particleVec, particleX2Y2);
    reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag,  particleX2Y2);
    reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rParticleX2Y2);
  }
  if (rankX2Z1 >= 0) { // edge x2z1
    Rectangle containerX2Z1(v2.getX() - cellSize, v1.getY(), v1.getZ(),
			    v2.getX(), v2.getY(), v1.getZ() + cellSize);
    findParticleInRectangle(containerX2Z1, particleVec, particleX2Z1);
    reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag,  particleX2Z1);
    reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rParticleX2Z1);
  }
  if (rankX2Z2 >= 0) { // edge x2z2
    Rectangle containerX2Z2(v2.getX() - cellSize, v1.getY(), v2.getZ() - cellSize,
			    v2.getX(), v2.getY(), v2.getZ());
    findParticleInRectangle(containerX2Z2, particleVec, particleX2Z2);
    reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag,  particleX2Z2);
    reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rParticleX2Z2);
  }
  if (rankY1Z1 >= 0) { // edge y1z1
    Rectangle containerY1Z1(v1.getX(), v1.getY(), v1.getZ(),
			    v2.getX(), v1.getY() + cellSize, v1.getZ() + cellSize);
    findParticleInRectangle(containerY1Z1, particleVec, particleY1Z1);
    reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag,  particleY1Z1);
    reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rParticleY1Z1);
  }
  if (rankY1Z2 >= 0) { // edge y1z2
    Rectangle containerY1Z2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
			    v2.getX(), v1.getY() + cellSize, v2.getZ());
    findParticleInRectangle(containerY1Z2, particleVec, particleY1Z2);
    reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag,  particleY1Z2);
    reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rParticleY1Z2);
  }
  if (rankY2Z1 >= 0) { // edge y2z1
    Rectangle containerY2Z1(v1.getX(), v2.getY() - cellSize, v1.getZ(),
			    v2.getX(), v2.getY(), v1.getZ() + cellSize);
    findParticleInRectangle(containerY2Z1, particleVec, particleY2Z1);
    reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag,  particleY2Z1);
    reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rParticleY2Z1);
  }
  if (rankY2Z2 >= 0) { // edge y2z2
    Rectangle containerY2Z2(v1.getX(), v2.getY() - cellSize, v2.getZ() - cellSize,
			    v2.getX(), v2.getY(), v2.getZ());
    findParticleInRectangle(containerY2Z2, particleVec, particleY2Z2);
    reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag,  particleY2Z2);
    reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rParticleY2Z2);
  }
  // 8 vertices
  if (rankX1Y1Z1 >= 0) { // edge x1y1z1
    Rectangle containerX1Y1Z1(v1.getX(), v1.getY(), v1.getZ(),
			      v1.getX() + cellSize, v1.getY() + cellSize, v1.getZ() + cellSize);
    findParticleInRectangle(containerX1Y1Z1, particleVec, particleX1Y1Z1);
    reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag,  particleX1Y1Z1);
    reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rParticleX1Y1Z1);
  }
  if (rankX1Y1Z2 >= 0) { // edge x1y1z2
    Rectangle containerX1Y1Z2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
			      v1.getX() + cellSize, v1.getY() + cellSize, v2.getZ());
    findParticleInRectangle(containerX1Y1Z2, particleVec, particleX1Y1Z2);
    reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag,  particleX1Y1Z2);
    reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rParticleX1Y1Z2);
  }
  if (rankX1Y2Z1 >= 0) { // edge x1y2z1
    Rectangle containerX1Y2Z1(v1.getX(), v2.getY() - cellSize, v1.getZ(),
			      v1.getX() + cellSize, v2.getY(), v1.getZ() + cellSize);
    findParticleInRectangle(containerX1Y2Z1, particleVec, particleX1Y2Z1);
    reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag,  particleX1Y2Z1);
    reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rParticleX1Y2Z1);
  }
  if (rankX1Y2Z2 >= 0) { // edge x1y2z2
    Rectangle containerX1Y2Z2(v1.getX(), v2.getY() - cellSize, v2.getZ() - cellSize,
			      v1.getX() + cellSize, v2.getY() + cellSize, v2.getZ());
    findParticleInRectangle(containerX1Y2Z2, particleVec, particleX1Y2Z2);
    reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag,  particleX1Y2Z2);
    reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rParticleX1Y2Z2);
  }
  if (rankX2Y1Z1 >= 0) { // edge x2y1z1
    Rectangle containerX2Y1Z1(v2.getX() - cellSize, v1.getY(), v1.getZ(),
			      v2.getX(), v1.getY() + cellSize, v1.getZ() + cellSize);
    findParticleInRectangle(containerX2Y1Z1, particleVec, particleX2Y1Z1);
    reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag,  particleX2Y1Z1);
    reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rParticleX2Y1Z1);
  }
  if (rankX2Y1Z2 >= 0) { // edge x2y1z2
    Rectangle containerX2Y1Z2(v2.getX() - cellSize, v1.getY(), v2.getZ() - cellSize,
			      v2.getX(), v1.getY() + cellSize, v2.getZ());
    findParticleInRectangle(containerX2Y1Z2, particleVec, particleX2Y1Z2);
    reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag,  particleX2Y1Z2);
    reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rParticleX2Y1Z2);
  }
  if (rankX2Y2Z1 >= 0) { // edge x2y2z1
    Rectangle containerX2Y2Z1(v2.getX() - cellSize, v2.getY() - cellSize, v1.getZ(),
			      v2.getX(), v2.getY(), v1.getZ() + cellSize);
    findParticleInRectangle(containerX2Y2Z1, particleVec, particleX2Y2Z1);
    reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag,  particleX2Y2Z1);
    reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rParticleX2Y2Z1);
  }
  if (rankX2Y2Z2 >= 0) { // edge x2y2z2
    Rectangle containerX2Y2Z2(v2.getX() - cellSize, v2.getY() - cellSize, v2.getZ() - cellSize,
			      v2.getX(), v2.getY(), v2.getZ());
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
  recvParticleVec.clear();
  // 6 surfaces
  if (rankX1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1.begin(), rParticleX1.end());
  if (rankX2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2.begin(), rParticleX2.end());
  if (rankY1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1.begin(), rParticleY1.end());
  if (rankY2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2.begin(), rParticleY2.end());
  if (rankZ1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleZ1.begin(), rParticleZ1.end());
  if (rankZ2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleZ2.begin(), rParticleZ2.end());
  // 12 edges
  if (rankX1Y1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1.begin(), rParticleX1Y1.end());
  if (rankX1Y2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2.begin(), rParticleX1Y2.end());
  if (rankX1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z1.begin(), rParticleX1Z1.end());
  if (rankX1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z2.begin(), rParticleX1Z2.end());
  if (rankX2Y1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1.begin(), rParticleX2Y1.end());
  if (rankX2Y2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2.begin(), rParticleX2Y2.end());
  if (rankX2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z1.begin(), rParticleX2Z1.end());
  if (rankX2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z2.begin(), rParticleX2Z2.end());
  if (rankY1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z1.begin(), rParticleY1Z1.end());
  if (rankY1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z2.begin(), rParticleY1Z2.end());
  if (rankY2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z1.begin(), rParticleY2Z1.end());
  if (rankY2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z2.begin(), rParticleY2Z2.end());
  // 8 vertices
  if (rankX1Y1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z1.begin(), rParticleX1Y1Z1.end());
  if (rankX1Y1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z2.begin(), rParticleX1Y1Z2.end());
  if (rankX1Y2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z1.begin(), rParticleX1Y2Z1.end());
  if (rankX1Y2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z2.begin(), rParticleX1Y2Z2.end());
  if (rankX2Y1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z1.begin(), rParticleX2Y1Z1.end());
  if (rankX2Y1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z2.begin(), rParticleX2Y1Z2.end());
  if (rankX2Y2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z1.begin(), rParticleX2Y2Z1.end());
  if (rankX2Y2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z2.begin(), rParticleX2Y2Z2.end());

  mergeParticleVec.clear();
  mergeParticleVec = particleVec; // duplicate pointers, pointing to the same memory
  mergeParticleVec.insert(mergeParticleVec.end(), recvParticleVec.begin(), recvParticleVec.end());

  /*
  std::vector<Particle*> testParticleVec;
  testParticleVec.insert(testParticleVec.end(), rParticleX1.begin(), rParticleX1.end());
  testParticleVec.insert(testParticleVec.end(), rParticleX2.begin(), rParticleX2.end());
  testParticleVec.insert(testParticleVec.end(), rParticleY1.begin(), rParticleY1.end());
  testParticleVec.insert(testParticleVec.end(), rParticleY2.begin(), rParticleY2.end());
  testParticleVec.insert(testParticleVec.end(), rParticleZ1.begin(), rParticleZ1.end());
  testParticleVec.insert(testParticleVec.end(), rParticleZ2.begin(), rParticleZ2.end());
  std::cout << "iter=" << std::setw(4) << iteration << " rank=" << std::setw(4) << mpiRank 
	    << " ptclNum=" << std::setw(4) << particleVec.size() 
	    << " surface="
	    << std::setw(4) << particleX1.size()  << std::setw(4) << particleX2.size()
	    << std::setw(4) << particleY1.size()  << std::setw(4) << particleY2.size()
	    << std::setw(4) << particleZ1.size()  << std::setw(4) << particleZ2.size()  
	    << " recv="
	    << std::setw(4) << rParticleX1.size() << std::setw(4) << rParticleX2.size()
	    << std::setw(4) << rParticleY1.size() << std::setw(4) << rParticleY2.size()
	    << std::setw(4) << rParticleZ1.size() << std::setw(4) << rParticleZ2.size() 
	    << " rNum="    
	    << std::setw(4) << recvParticleVec.size() << ": ";   

  for (std::vector<Particle*>::const_iterator it = testParticleVec.begin(); it != testParticleVec.end();++it)
    std::cout << (*it)->getId() << ' ';
  std::cout << std::endl;
  testParticleVec.clear();
  */
}


void Assembly::releaseRecvParticle() {
  // release memory of received particles
  for (std::vector<Particle*>::iterator it = recvParticleVec.begin(); it != recvParticleVec.end(); ++it)
    delete (*it);
  recvParticleVec.clear();
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
}


void Assembly::updateContainerMinX() {
  REAL pMinX = getPtclMinX(particleVec);
  REAL minX = 0;
  MPI_Allreduce(&pMinX, &minX, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);

  setContainer(Rectangle(minX - gradation.getPtclMaxRadius(),
			 allContainer.getMinCorner().getY(),
			 allContainer.getMinCorner().getZ(),
			 allContainer.getMaxCorner().getX(),
			 allContainer.getMaxCorner().getY(),
			 allContainer.getMaxCorner().getZ() ));
}


void Assembly::updateContainerMaxX() {
  REAL pMaxX = getPtclMaxX(particleVec);
  REAL maxX = 0;
  MPI_Allreduce(&pMaxX, &maxX, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);

  setContainer(Rectangle(allContainer.getMinCorner().getX(),
			 allContainer.getMinCorner().getY(),
			 allContainer.getMinCorner().getZ(),
			 maxX + gradation.getPtclMaxRadius(),
			 allContainer.getMaxCorner().getY(),
			 allContainer.getMaxCorner().getZ() ));
}


void Assembly::updateContainerMinY() {
  REAL pMinY = getPtclMinY(particleVec);
  REAL minY = 0;
  MPI_Allreduce(&pMinY, &minY, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);

  setContainer(Rectangle(allContainer.getMinCorner().getX(),
			 minY - gradation.getPtclMaxRadius(),
			 allContainer.getMinCorner().getZ(),
			 allContainer.getMaxCorner().getX(),
			 allContainer.getMaxCorner().getY(),
			 allContainer.getMaxCorner().getZ() ));
}


void Assembly::updateContainerMaxY() {
  REAL pMaxY = getPtclMaxY(particleVec);
  REAL maxY = 0;
  MPI_Allreduce(&pMaxY, &maxY, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);

  setContainer(Rectangle(allContainer.getMinCorner().getX(),
			 allContainer.getMinCorner().getY(),
			 allContainer.getMinCorner().getZ(),
			 allContainer.getMaxCorner().getX(),
			 maxY + gradation.getPtclMaxRadius(),
			 allContainer.getMaxCorner().getZ() ));
}


void Assembly::updateContainerMinZ() {
  REAL pMinZ = getPtclMinZ(particleVec);
  REAL minZ = 0;
  MPI_Allreduce(&pMinZ, &minZ, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);

  setContainer(Rectangle(allContainer.getMinCorner().getX(),
			 allContainer.getMinCorner().getY(),
			 minZ - gradation.getPtclMaxRadius(),
			 allContainer.getMaxCorner().getX(),
			 allContainer.getMaxCorner().getY(),
			 allContainer.getMaxCorner().getZ() ));
}


void Assembly::updateContainerMaxZ() {
  // update allContainer adaptively due to particle motion
  REAL pMaxZ = getPtclMaxZ(particleVec);
  REAL maxZ = 0;
  MPI_Allreduce(&pMaxZ, &maxZ, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);

  // no need to broadcast allContainer as it is updated in each process
  setContainer(Rectangle(allContainer.getMinCorner().getX(),
			 allContainer.getMinCorner().getY(),
			 allContainer.getMinCorner().getZ(),
			 allContainer.getMaxCorner().getX(),
			 allContainer.getMaxCorner().getY(),
			 maxZ + gradation.getPtclMaxRadius() ));
  /*
  std::vector<int> procGroup;
  for (int iRank = 0; iRank < mpiSize; ++iRank) {
    int ndim = 3;
    int coords[3];
    MPI_Cart_coords(cartComm, iRank, ndim, coords);
    if (coords[2] == mpiProcZ - 1) // process at top level
      procGroup.push_back(iRank);
  }

  int size = procGroup.size();
  int *ranks = new int[size];
  for (int i = 0; i < size; ++i)
    ranks[i] = procGroup[i];
  
  MPI_Group groupWorld, groupWorker;
  MPI_Comm commWorker;
  MPI_Comm_group(mpiWorld, &groupWorld);
  MPI_Group_incl(groupWorld, size, ranks, &groupWorker);
  MPI_Comm_create(mpiWorld, groupWorker, &commWorker);

  REAL pMaxZ = getPtclMaxZ(particleVec);
  REAL maxZ = 0;
  MPI_Allreduce(&pMaxZ, &maxZ, 1, MPI_DOUBLE, MPI_MAX, commWorker);

  std::vector<int>::iterator it;
  it = find(procGroup.begin(), procGroup.end(), mpiRank);
  if (it != procGroup.end())
  std::cout << "iter=" << std::setw(8) << iteration << " rank=" << std::setw(2) << mpiRank
	    << " pNum=" << particleVec.size() << " topProc=" << size
	    << " maxZ=" << std::scientific << std::setprecision(10) << pMaxZ << std::fixed << std::endl;

  MPI_Group_free(&groupWorld);
  MPI_Group_free(&groupWorker);
  MPI_Comm_free(&commWorker);
  delete [] ranks;
  */
}


void Assembly::updateGridMaxZ() {
  // update compute grids adaptively due to particle motion
  REAL pMaxZ = getPtclMaxZ(particleVec);
  REAL maxZ = 0;
  MPI_Allreduce(&pMaxZ, &maxZ, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);

  // no need to broadcast grid as it is updated in each process
  setGrid(Rectangle(grid.getMinCorner().getX(),
		    grid.getMinCorner().getY(),
		    grid.getMinCorner().getZ(),
		    grid.getMaxCorner().getX(),
		    grid.getMaxCorner().getY(),
		    maxZ + gradation.getPtclMaxRadius() ));
}


void Assembly::migrateParticle() 
{
  Vec vspan = grid.getMaxCorner() - grid.getMinCorner();
  double segX = vspan.getX() / mpiProcX;
  double segY = vspan.getY() / mpiProcY;
  double segZ = vspan.getZ() / mpiProcZ;
  Vec v1 = container.getMinCorner(); // v1, v2 in terms of process
  Vec v2 = container.getMaxCorner();  

  // if a neighbor exists, transfer particles crossing the boundary in between.
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

  // 6 surfaces
  if (rankX1 >= 0) { // surface x1
    Rectangle containerX1(v1.getX() - segX, v1.getY(), v1.getZ(), 
			  v1.getX(), v2.getY(), v2.getZ());
    findParticleInRectangle(containerX1, particleVec, particleX1);
    reqX1[0] = boostWorld.isend(rankX1, mpiTag,  particleX1);
    reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rParticleX1);
  }
  if (rankX2 >= 0) { // surface x2
    Rectangle containerX2(v2.getX(), v1.getY(), v1.getZ(),
			  v2.getX() + segX, v2.getY(), v2.getZ());
    findParticleInRectangle(containerX2, particleVec, particleX2);
    reqX2[0] = boostWorld.isend(rankX2, mpiTag,  particleX2);
    reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rParticleX2);
  }
  if (rankY1 >= 0) {  // surface y1
    Rectangle containerY1(v1.getX(), v1.getY() - segY, v1.getZ(), 
			  v2.getX(), v1.getY(), v2.getZ());
    findParticleInRectangle(containerY1, particleVec, particleY1);
    reqY1[0] = boostWorld.isend(rankY1, mpiTag,  particleY1);
    reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rParticleY1);
  }
  if (rankY2 >= 0) {  // surface y2
    Rectangle containerY2(v1.getX(), v2.getY(), v1.getZ(),
			  v2.getX(), v2.getY() + segY, v2.getZ());
    findParticleInRectangle(containerY2, particleVec, particleY2);
    reqY2[0] = boostWorld.isend(rankY2, mpiTag,  particleY2);
    reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rParticleY2);
  }
  if (rankZ1 >= 0) {  // surface z1
    Rectangle containerZ1(v1.getX(), v1.getY(), v1.getZ() - segZ,
			  v2.getX(), v2.getY(), v1.getZ());
    findParticleInRectangle(containerZ1, particleVec, particleZ1);
    reqZ1[0] = boostWorld.isend(rankZ1, mpiTag,  particleZ1);
    reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rParticleZ1);
  }
  if (rankZ2 >= 0) {  // surface z2
    Rectangle containerZ2(v1.getX(), v1.getY(), v2.getZ(),
			  v2.getX(), v2.getY(), v2.getZ() + segZ);
    findParticleInRectangle(containerZ2, particleVec, particleZ2);
    reqZ2[0] = boostWorld.isend(rankZ2, mpiTag,  particleZ2);
    reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rParticleZ2);
  }
  // 12 edges
  if (rankX1Y1 >= 0) { // edge x1y1
    Rectangle containerX1Y1(v1.getX() - segX, v1.getY() - segY, v1.getZ(),
			    v1.getX(), v1.getY(), v2.getZ());
    findParticleInRectangle(containerX1Y1, particleVec, particleX1Y1);
    reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag,  particleX1Y1);
    reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rParticleX1Y1);
  }
  if (rankX1Y2 >= 0) { // edge x1y2
    Rectangle containerX1Y2(v1.getX() - segX, v2.getY(), v1.getZ(),
			    v1.getX(), v2.getY() + segY, v2.getZ());
    findParticleInRectangle(containerX1Y2, particleVec, particleX1Y2);
    reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag,  particleX1Y2);
    reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rParticleX1Y2);
  }
  if (rankX1Z1 >= 0) { // edge x1z1
    Rectangle containerX1Z1(v1.getX() - segX, v1.getY(), v1.getZ() -segZ,
			    v1.getX(), v2.getY(), v1.getZ());
    findParticleInRectangle(containerX1Z1, particleVec, particleX1Z1);
    reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag,  particleX1Z1);
    reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rParticleX1Z1);
  }
  if (rankX1Z2 >= 0) { // edge x1z2
    Rectangle containerX1Z2(v1.getX() - segX, v1.getY(), v2.getZ(),
			    v1.getX(), v2.getY(), v2.getZ() + segZ);
    findParticleInRectangle(containerX1Z2, particleVec, particleX1Z2);
    reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag,  particleX1Z2);
    reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rParticleX1Z2);
  }
  if (rankX2Y1 >= 0) { // edge x2y1
    Rectangle containerX2Y1(v2.getX(), v1.getY() - segY, v1.getZ(),
			    v2.getX() + segX, v1.getY(), v2.getZ());
    findParticleInRectangle(containerX2Y1, particleVec, particleX2Y1);
    reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag,  particleX2Y1);
    reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rParticleX2Y1);
  }
  if (rankX2Y2 >= 0) { // edge x2y2
    Rectangle containerX2Y2(v2.getX(), v2.getY(), v1.getZ(),
			    v2.getX() + segX, v2.getY() + segY, v2.getZ());
    findParticleInRectangle(containerX2Y2, particleVec, particleX2Y2);
    reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag,  particleX2Y2);
    reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rParticleX2Y2);
  }
  if (rankX2Z1 >= 0) { // edge x2z1
    Rectangle containerX2Z1(v2.getX(), v1.getY(), v1.getZ() - segZ,
			    v2.getX() + segX, v2.getY(), v1.getZ());
    findParticleInRectangle(containerX2Z1, particleVec, particleX2Z1);
    reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag,  particleX2Z1);
    reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rParticleX2Z1);
  }
  if (rankX2Z2 >= 0) { // edge x2z2
    Rectangle containerX2Z2(v2.getX(), v1.getY(), v2.getZ(),
			    v2.getX() + segX, v2.getY(), v2.getZ() + segZ);
    findParticleInRectangle(containerX2Z2, particleVec, particleX2Z2);
    reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag,  particleX2Z2);
    reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rParticleX2Z2);
  }
  if (rankY1Z1 >= 0) { // edge y1z1
    Rectangle containerY1Z1(v1.getX(), v1.getY() - segY, v1.getZ() - segZ,
			    v2.getX(), v1.getY(), v1.getZ());
    findParticleInRectangle(containerY1Z1, particleVec, particleY1Z1);
    reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag,  particleY1Z1);
    reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rParticleY1Z1);
  }
  if (rankY1Z2 >= 0) { // edge y1z2
    Rectangle containerY1Z2(v1.getX(), v1.getY() - segY, v2.getZ(),
			    v2.getX(), v1.getY(), v2.getZ() + segZ);
    findParticleInRectangle(containerY1Z2, particleVec, particleY1Z2);
    reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag,  particleY1Z2);
    reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rParticleY1Z2);
  }
  if (rankY2Z1 >= 0) { // edge y2z1
    Rectangle containerY2Z1(v1.getX(), v2.getY(), v1.getZ() - segZ,
			    v2.getX(), v2.getY() + segY, v1.getZ());
    findParticleInRectangle(containerY2Z1, particleVec, particleY2Z1);
    reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag,  particleY2Z1);
    reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rParticleY2Z1);
  }
  if (rankY2Z2 >= 0) { // edge y2z2
    Rectangle containerY2Z2(v1.getX(), v2.getY(), v2.getZ(),
			    v2.getX(), v2.getY() + segY, v2.getZ() + segZ);
    findParticleInRectangle(containerY2Z2, particleVec, particleY2Z2);
    reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag,  particleY2Z2);
    reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rParticleY2Z2);
  }
  // 8 vertices
  if (rankX1Y1Z1 >= 0) { // edge x1y1z1
    Rectangle containerX1Y1Z1(v1.getX() - segX, v1.getY() - segY, v1.getZ() - segZ,
			      v1.getX(), v1.getY(), v1.getZ());
    findParticleInRectangle(containerX1Y1Z1, particleVec, particleX1Y1Z1);
    reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag,  particleX1Y1Z1);
    reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rParticleX1Y1Z1);
  }
  if (rankX1Y1Z2 >= 0) { // edge x1y1z2
    Rectangle containerX1Y1Z2(v1.getX() - segX, v1.getY() - segY, v2.getZ(),
			      v1.getX(), v1.getY(), v2.getZ() + segZ);
    findParticleInRectangle(containerX1Y1Z2, particleVec, particleX1Y1Z2);
    reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag,  particleX1Y1Z2);
    reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rParticleX1Y1Z2);
  }
  if (rankX1Y2Z1 >= 0) { // edge x1y2z1
    Rectangle containerX1Y2Z1(v1.getX() - segX, v2.getY(), v1.getZ() - segZ,
			      v1.getX(), v2.getY() + segY, v1.getZ());
    findParticleInRectangle(containerX1Y2Z1, particleVec, particleX1Y2Z1);
    reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag,  particleX1Y2Z1);
    reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rParticleX1Y2Z1);
  }
  if (rankX1Y2Z2 >= 0) { // edge x1y2z2
    Rectangle containerX1Y2Z2(v1.getX() - segX, v2.getY(), v2.getZ(),
			      v1.getX(), v2.getY() + segY, v2.getZ() + segZ);
    findParticleInRectangle(containerX1Y2Z2, particleVec, particleX1Y2Z2);
    reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag,  particleX1Y2Z2);
    reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rParticleX1Y2Z2);
  }
  if (rankX2Y1Z1 >= 0) { // edge x2y1z1
    Rectangle containerX2Y1Z1(v2.getX(), v1.getY() - segY, v1.getZ() - segZ,
			      v2.getX() + segX, v1.getY(), v1.getZ());
    findParticleInRectangle(containerX2Y1Z1, particleVec, particleX2Y1Z1);
    reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag,  particleX2Y1Z1);
    reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rParticleX2Y1Z1);
  }
  if (rankX2Y1Z2 >= 0) { // edge x2y1z2
    Rectangle containerX2Y1Z2(v2.getX(), v1.getY() - segY, v2.getZ(),
			      v2.getX() + segX, v1.getY(), v2.getZ() + segZ);
    findParticleInRectangle(containerX2Y1Z2, particleVec, particleX2Y1Z2);
    reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag,  particleX2Y1Z2);
    reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rParticleX2Y1Z2);
  }
  if (rankX2Y2Z1 >= 0) { // edge x2y2z1
    Rectangle containerX2Y2Z1(v2.getX(), v2.getY(), v1.getZ() - segZ,
			      v2.getX() + segX, v2.getY() + segY, v1.getZ());
    findParticleInRectangle(containerX2Y2Z1, particleVec, particleX2Y2Z1);
    reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag,  particleX2Y2Z1);
    reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rParticleX2Y2Z1);
  }
  if (rankX2Y2Z2 >= 0) { // edge x2y2z2
    Rectangle containerX2Y2Z2(v2.getX(), v2.getY(), v2.getZ(),
			      v2.getX() + segX, v2.getY() + segY, v2.getZ() + segZ);
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

  // delete outgoing particles
  removeParticleOutRectangle();

  // add incoming particles
  recvParticleVec.clear(); // new use of recvParticleVec
  // 6 surfaces
  if (rankX1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1.begin(), rParticleX1.end());
  if (rankX2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2.begin(), rParticleX2.end());
  if (rankY1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1.begin(), rParticleY1.end());
  if (rankY2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2.begin(), rParticleY2.end());
  if (rankZ1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleZ1.begin(), rParticleZ1.end());
  if (rankZ2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleZ2.begin(), rParticleZ2.end());
  // 12 edges
  if (rankX1Y1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1.begin(), rParticleX1Y1.end());
  if (rankX1Y2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2.begin(), rParticleX1Y2.end());
  if (rankX1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z1.begin(), rParticleX1Z1.end());
  if (rankX1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z2.begin(), rParticleX1Z2.end());
  if (rankX2Y1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1.begin(), rParticleX2Y1.end());
  if (rankX2Y2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2.begin(), rParticleX2Y2.end());
  if (rankX2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z1.begin(), rParticleX2Z1.end());
  if (rankX2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z2.begin(), rParticleX2Z2.end());
  if (rankY1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z1.begin(), rParticleY1Z1.end());
  if (rankY1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z2.begin(), rParticleY1Z2.end());
  if (rankY2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z1.begin(), rParticleY2Z1.end());
  if (rankY2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z2.begin(), rParticleY2Z2.end());
  // 8 vertices
  if (rankX1Y1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z1.begin(), rParticleX1Y1Z1.end());
  if (rankX1Y1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z2.begin(), rParticleX1Y1Z2.end());
  if (rankX1Y2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z1.begin(), rParticleX1Y2Z1.end());
  if (rankX1Y2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z2.begin(), rParticleX1Y2Z2.end());
  if (rankX2Y1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z1.begin(), rParticleX2Y1Z1.end());
  if (rankX2Y1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z2.begin(), rParticleX2Y1Z2.end());
  if (rankX2Y2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z1.begin(), rParticleX2Y2Z1.end());
  if (rankX2Y2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z2.begin(), rParticleX2Y2Z2.end());

  particleVec.insert(particleVec.end(), recvParticleVec.begin(), recvParticleVec.end());

  /*
  if (recvParticleVec.size() > 0) {    
    std::cout << "iter=" << std::setw(8) << iteration << " rank=" << std::setw(2) << mpiRank 
	     << "   added=";
    for (std::vector<Particle*>::const_iterator it = recvParticleVec.begin(); it != recvParticleVec.end(); ++it)
      std::cout << std::setw(3) << (*it)->getId();
    std::cout << " now " << particleVec.size() << ": ";
    for (std::vector<Particle*>::const_iterator it = particleVec.begin(); it != particleVec.end(); ++it)
      std::cout << std::setw(3) << (*it)->getId();
    std::cout << std::endl;
  }
  */

  // do not release memory of received particles because they are part of and managed by particleVec
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
}


void Assembly::gatherParticle() {

  // update allParticleVec: process 0 collects all updated particles from each other process  
  if (mpiRank != 0) {// each process except 0
    boostWorld.send(0, mpiTag, particleVec);
  }
  else { // process 0
    // allParticleVec is cleared first time in scatterParticle, and at the end of each gathering.

    // duplicate particleVec so that it is not destroyed by allParticleVec in next iteration,
    // otherwise it causes memory error.
    std::vector<Particle*> dupParticleVec(particleVec.size());
    for (int i = 0; i < dupParticleVec.size(); ++i)
      dupParticleVec[i] = new Particle(*particleVec[i]);

    // fill allParticleVec with dupParticleVec and received particles
    allParticleVec.insert(allParticleVec.end(), dupParticleVec.begin(), dupParticleVec.end());

    std::vector<Particle*> tmpParticleVec;
    long gatherRam = 0;
    for (int iRank = 1; iRank < mpiSize; ++iRank) {

      tmpParticleVec.clear();// do not destroy particles!
      boostWorld.recv(iRank, mpiTag, tmpParticleVec);
      allParticleVec.insert(allParticleVec.end(), tmpParticleVec.begin(), tmpParticleVec.end());
      gatherRam += tmpParticleVec.size();

    }
    //std::cout << "gather: particleNum = " << gatherRam <<  " particleRam = " << gatherRam * sizeof(Particle) << std::endl;
  }
  
}

void Assembly::releaseGatheredParticle() {
  // clear allParticleVec, avoid long time memory footprint.
  for (std::vector<Particle*>::iterator it = allParticleVec.begin(); it != allParticleVec.end(); ++it)
    delete (*it);
  allParticleVec.clear();
  std::vector<Particle*>().swap(allParticleVec); // actual memory release
}


void Assembly::readParticle(const char *inputParticle) {

  REAL young = dem::Parameter::getSingleton().parameter["young"];
  REAL poisson = dem::Parameter::getSingleton().parameter["poisson"];

  std::ifstream ifs(inputParticle);
  if(!ifs) { std::cout << "stream error: readParticle" << std::endl; exit(-1); }
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
    Particle* pt= new Particle(id, type, Vec(a,b,c), Vec(px,py,pz), Vec(dax,day,daz), Vec(dbx,dby,dbz), Vec(dcx,dcy,dcz), young, poisson);
    
    // optional settings for a particle's initial status
    if ( (static_cast<int> (dem::Parameter::getSingleton().parameter["toInitParticle"])) == 1 ) {
      pt->setPrevVeloc(Vec(vx,vy,vz));
      pt->setCurrVeloc(Vec(vx,vy,vz));
      pt->setPrevOmga(Vec(omx,omy,omz));
      pt->setCurrOmga(Vec(omx,omy,omz));
      pt->setForce(Vec(fx,fy,fz));  // initial force
      pt->setMoment(Vec(mx,my,mz)); // initial moment
    }
    
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


void Assembly::printParticle(const char *str) const {
  std::ofstream ofs(str);
  if(!ofs) { std::cout << "stream error: printParticle" << std::endl; exit(-1); }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  ofs << std::setw(OWID) << allParticleVec.size() << std::endl;
  ofs << std::setw(OWID) << "id"
      << std::setw(OWID) << "type"
      << std::setw(OWID) << "radius_a"
      << std::setw(OWID) << "radius_b"
      << std::setw(OWID) << "radius_c"
      << std::setw(OWID) << "position_x"
      << std::setw(OWID) << "position_y"
      << std::setw(OWID) << "position_z"
      << std::setw(OWID) << "axle_a_x"
      << std::setw(OWID) << "axle_a_y"
      << std::setw(OWID) << "axle_a_z"
      << std::setw(OWID) << "axle_b_x"
      << std::setw(OWID) << "axle_b_y"
      << std::setw(OWID) << "axle_b_z"
      << std::setw(OWID) << "axle_c_x"
      << std::setw(OWID) << "axle_c_y"
      << std::setw(OWID) << "axle_c_z"
      << std::setw(OWID) << "velocity_x"
      << std::setw(OWID) << "velocity_y"
      << std::setw(OWID) << "velocity_z"
      << std::setw(OWID) << "omga_x"
      << std::setw(OWID) << "omga_y"
      << std::setw(OWID) << "omga_z"
      << std::setw(OWID) << "force_x"
      << std::setw(OWID) << "force_y"
      << std::setw(OWID) << "force_z"
      << std::setw(OWID) << "moment_x"
      << std::setw(OWID) << "moment_y"
      << std::setw(OWID) << "moment_z"
      << std::endl;
  
  Vec vObj;
  std::vector<Particle*>::const_iterator  it;
  for (it = allParticleVec.begin(); it != allParticleVec.end();++it)  {
    ofs << std::setw(OWID) << (*it)->getId()
	<< std::setw(OWID) << (*it)->getType()
	<< std::setw(OWID) << (*it)->getA()
	<< std::setw(OWID) << (*it)->getB()
	<< std::setw(OWID) << (*it)->getC();
    
    vObj=(*it)->getCurrPos();
    ofs << std::setw(OWID) << vObj.getX()
	<< std::setw(OWID) << vObj.getY()
	<< std::setw(OWID) << vObj.getZ();
    
    vObj=(*it)->getCurrDirecA();
    ofs << std::setw(OWID) << vObj.getX()
	<< std::setw(OWID) << vObj.getY()
	<< std::setw(OWID) << vObj.getZ();
    
    vObj=(*it)->getCurrDirecB();
    ofs << std::setw(OWID) << vObj.getX()
	<< std::setw(OWID) << vObj.getY()
	<< std::setw(OWID) << vObj.getZ();
    
    vObj=(*it)->getCurrDirecC();
    ofs << std::setw(OWID) << vObj.getX()
	<< std::setw(OWID) << vObj.getY()
	<< std::setw(OWID) << vObj.getZ();
    
    vObj=(*it)->getCurrVeloc();
    ofs << std::setw(OWID) << vObj.getX()
	<< std::setw(OWID) << vObj.getY()
	<< std::setw(OWID) << vObj.getZ();
    
    vObj=(*it)->getCurrOmga();
    ofs << std::setw(OWID) << vObj.getX()
	<< std::setw(OWID) << vObj.getY()
	<< std::setw(OWID) << vObj.getZ();
    
    vObj=(*it)->getForce();
    ofs << std::setw(OWID) << vObj.getX()
	<< std::setw(OWID) << vObj.getY()
	<< std::setw(OWID) << vObj.getZ();
    
    vObj=(*it)->getMoment();
    ofs << std::setw(OWID) << vObj.getX()
	<< std::setw(OWID) << vObj.getY()
	<< std::setw(OWID) << vObj.getZ() << std::endl;
  }
  
  int sieveNum = gradation.getSieveNum();
  std::vector<REAL> percent = gradation.getPercent();
  std::vector<REAL> size    = gradation.getSize();
  ofs << std::endl << std::setw(OWID) << sieveNum << std::endl;
  for (int i = 0; i < sieveNum; ++i)
    ofs << std::setw(OWID) << percent[i] << std::setw(OWID) << size[i] << std::endl;
  ofs << std::endl << std::setw(OWID) << gradation.getPtclRatioBA() << std::setw(OWID) << gradation.getPtclRatioCA() << std::endl;

  ofs.close();
}


void Assembly::printParticle(const char *str, std::vector<Particle*>  &particleVec) const {
  std::ofstream ofs(str);
  if(!ofs) { std::cout << "stream error: printParticle" << std::endl; exit(-1); }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  ofs << std::setw(OWID) << particleVec.size() << std::endl;
  ofs << std::setw(OWID) << "id"
      << std::setw(OWID) << "type"
      << std::setw(OWID) << "radius_a"
      << std::setw(OWID) << "radius_b"
      << std::setw(OWID) << "radius_c"
      << std::setw(OWID) << "position_x"
      << std::setw(OWID) << "position_y"
      << std::setw(OWID) << "position_z"
      << std::setw(OWID) << "axle_a_x"
      << std::setw(OWID) << "axle_a_y"
      << std::setw(OWID) << "axle_a_z"
      << std::setw(OWID) << "axle_b_x"
      << std::setw(OWID) << "axle_b_y"
      << std::setw(OWID) << "axle_b_z"
      << std::setw(OWID) << "axle_c_x"
      << std::setw(OWID) << "axle_c_y"
      << std::setw(OWID) << "axle_c_z"
      << std::setw(OWID) << "velocity_x"
      << std::setw(OWID) << "velocity_y"
      << std::setw(OWID) << "velocity_z"
      << std::setw(OWID) << "omga_x"
      << std::setw(OWID) << "omga_y"
      << std::setw(OWID) << "omga_z"
      << std::setw(OWID) << "force_x"
      << std::setw(OWID) << "force_y"
      << std::setw(OWID) << "force_z"
      << std::setw(OWID) << "moment_x"
      << std::setw(OWID) << "moment_y"
      << std::setw(OWID) << "moment_z"
      << std::endl;
  
  Vec vObj;
  std::vector<Particle*>::const_iterator  it;
  for (it = particleVec.begin(); it != particleVec.end();++it)  {
    ofs << std::setw(OWID) << (*it)->getId()
	<< std::setw(OWID) << (*it)->getType()
	<< std::setw(OWID) << (*it)->getA()
	<< std::setw(OWID) << (*it)->getB()
	<< std::setw(OWID) << (*it)->getC();
    
    vObj=(*it)->getCurrPos();
    ofs << std::setw(OWID) << vObj.getX()
	<< std::setw(OWID) << vObj.getY()
	<< std::setw(OWID) << vObj.getZ();
    
    vObj=(*it)->getCurrDirecA();
    ofs << std::setw(OWID) << vObj.getX()
	<< std::setw(OWID) << vObj.getY()
	<< std::setw(OWID) << vObj.getZ();
    
    vObj=(*it)->getCurrDirecB();
    ofs << std::setw(OWID) << vObj.getX()
	<< std::setw(OWID) << vObj.getY()
	<< std::setw(OWID) << vObj.getZ();
    
    vObj=(*it)->getCurrDirecC();
    ofs << std::setw(OWID) << vObj.getX()
	<< std::setw(OWID) << vObj.getY()
	<< std::setw(OWID) << vObj.getZ();
    
    vObj=(*it)->getCurrVeloc();
    ofs << std::setw(OWID) << vObj.getX()
	<< std::setw(OWID) << vObj.getY()
	<< std::setw(OWID) << vObj.getZ();
    
    vObj=(*it)->getCurrOmga();
    ofs << std::setw(OWID) << vObj.getX()
	<< std::setw(OWID) << vObj.getY()
	<< std::setw(OWID) << vObj.getZ();
    
    vObj=(*it)->getForce();
    ofs << std::setw(OWID) << vObj.getX()
	<< std::setw(OWID) << vObj.getY()
	<< std::setw(OWID) << vObj.getZ();
    
    vObj=(*it)->getMoment();
    ofs << std::setw(OWID) << vObj.getX()
	<< std::setw(OWID) << vObj.getY()
	<< std::setw(OWID) << vObj.getZ() << std::endl;
  }
  
  ofs.close();
}


void Assembly::readBoundary(const char *str) {
  std::ifstream ifs(str);
  if(!ifs) { std::cout << "stream error: readBoundary" << std::endl; exit(-1); }

  REAL x1, y1, z1, x2, y2, z2;
  ifs >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;
  setContainer(Rectangle(x1, y1, z1, x2, y2, z2));
  // compute grid assumed to be the same as container, change in scatterParticle() if necessary.
  setGrid(Rectangle(x1, y1, z1, x2, y2, z2)); 

  Boundary *rbptr;
  int type;
  boundaryVec.clear();
  int boundaryNum;
  ifs >> boundaryNum;
  for(int i = 0; i < boundaryNum; ++i){
    ifs >> type;
    if(type == 1) // plane boundary
      rbptr = new plnBoundary(ifs);
    else          // cylindrical boundary
      rbptr = new cylBoundary(ifs);
    boundaryVec.push_back(rbptr);
  }

  ifs.close();
}


void Assembly::printBoundary(const char *str) const {
  std::ofstream ofs(str);
  if(!ofs) { std::cout << "stream error: printBoundary" << std::endl; exit(-1); }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  
  Vec v1 = allContainer.getMinCorner();
  Vec v2 = allContainer.getMaxCorner();
  REAL x1 = v1.getX();
  REAL y1 = v1.getY();
  REAL z1 = v1.getZ();
  REAL x2 = v2.getX();
  REAL y2 = v2.getY();
  REAL z2 = v2.getZ();
  
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z1
      << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z2 << std::endl
      << std::setw(OWID) << boundaryVec.size() << std::endl;
  
  std::vector<Boundary*>::const_iterator rt;
  for(rt = boundaryVec.begin(); rt != boundaryVec.end(); ++rt)
    (*rt)->display(ofs);
  ofs << std::endl;
  
  ofs.close();
}


void Assembly::plotBoundary(const char *str) const {
  std::ofstream ofs(str);
  if(!ofs) { std::cout << "stream error: plotBoundary" << std::endl; exit(-1); }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  REAL x1, y1, z1, x2, y2, z2;
  x1 = allContainer.getMinCorner().getX();
  y1 = allContainer.getMinCorner().getY();
  z1 = allContainer.getMinCorner().getZ();
  x2 = allContainer.getMaxCorner().getX();
  y2 = allContainer.getMaxCorner().getY();
  z2 = allContainer.getMaxCorner().getZ();

  ofs << "ZONE N=8, E=1, DATAPACKING=POINT, ZONETYPE=FEBRICK" << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y1 << std::setw(OWID) << z1 << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z1 << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y2 << std::setw(OWID) << z1 << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z1 << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y1 << std::setw(OWID) << z2 << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z2 << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y2 << std::setw(OWID) << z2 << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z2 << std::endl;
  ofs << "1 2 3 4 5 6 7 8" << std::endl;

  ofs.close();
}


void Assembly::plotGrid(const char *str) const {
  std::ofstream ofs(str);
  if(!ofs) { std::cout << "stream error: plotGrid" << std::endl; exit(-1); }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  Vec v1 = grid.getMinCorner();
  Vec v2 = grid.getMaxCorner();
  Vec vspan = v2 - v1;

  ofs << "ZONE N=" << (mpiProcX + 1) * (mpiProcY + 1) * (mpiProcZ + 1)
      << ", E=" << mpiProcX * mpiProcY * mpiProcZ << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" << std::endl;

  std::vector<Vec> coords((mpiProcX + 1) * (mpiProcY + 1) * (mpiProcZ + 1));
  int index = 0;
  for (int i = 0; i < mpiProcX + 1; ++i)
    for (int j = 0; j < mpiProcY + 1; ++j)
      for (int k = 0; k < mpiProcZ + 1; ++k)
	coords[index++] = Vec(v1.getX() + vspan.getX() / mpiProcX * i,
			      v1.getY() + vspan.getY() / mpiProcY * j,
			      v1.getZ() + vspan.getZ() / mpiProcZ * k);

  for (int i = 0; i < (mpiProcX + 1) * (mpiProcY + 1) * (mpiProcZ + 1); ++i)
    ofs << std::setw(OWID) << coords[i].getX() 
	<< std::setw(OWID) << coords[i].getY() 
	<< std::setw(OWID) << coords[i].getZ() << std::endl;

  for (int iRank = 0; iRank < mpiSize; ++iRank) {
      int coords[3];
      MPI_Cart_coords(cartComm, iRank, 3, coords);

      int id4 = 1 + coords[0]*(mpiProcZ+1)*(mpiProcY+1) + coords[1]*(mpiProcZ+1) + coords[2];
      int id1 = 1 + (coords[0]+1)*(mpiProcZ+1)*(mpiProcY+1) + coords[1]*(mpiProcZ+1) + coords[2];
      int id3 = 1 + coords[0]*(mpiProcZ+1)*(mpiProcY+1) + (coords[1]+1)*(mpiProcZ+1) + coords[2];
      int id2 = 1 + (coords[0]+1)*(mpiProcZ+1)*(mpiProcY+1) + (coords[1]+1)*(mpiProcZ+1) + coords[2];

      int id8 = 1 + coords[0]*(mpiProcZ+1)*(mpiProcY+1) + coords[1]*(mpiProcZ+1) + (coords[2]+1);
      int id5 = 1 + (coords[0]+1)*(mpiProcZ+1)*(mpiProcY+1) + coords[1]*(mpiProcZ+1) + (coords[2]+1);
      int id7 = 1 + coords[0]*(mpiProcZ+1)*(mpiProcY+1) + (coords[1]+1)*(mpiProcZ+1) + (coords[2]+1);
      int id6 = 1 + (coords[0]+1)*(mpiProcZ+1)*(mpiProcY+1) + (coords[1]+1)*(mpiProcZ+1) + (coords[2]+1);

      ofs << std::setw(8) << id1 << std::setw(8) << id2 << std::setw(8) << id3 << std::setw(8) << id4 
	  << std::setw(8) << id5 << std::setw(8) << id6 << std::setw(8) << id7 << std::setw(8) << id8 << std::endl;
  }

  ofs.close();
}


void Assembly::findContact() { // various implementations

  int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];

  if (ompThreads == 1) { // non-openmp single-thread version, time complexity bigO(n x n), n is the number of particles.  
    contactVec.clear();
    possContactNum = 0;
    
#ifdef TIME_PROFILE
    double time_r = 0; // time consumed in contact resolution, i.e., tmpContact.isOverlapped()
    gettimeofday(&time_p1, NULL); 
#endif
    
    int num1 = particleVec.size();      // particles inside container
    int num2 = mergeParticleVec.size(); // paticles inside container (at front) + particles from neighboring blocks (at end)
    for (int i = 0; i < num1; ++i) {    // NOT (num1 - 1), in parallel situation where one particle could contact received particles!
      Vec u = particleVec[i]->getCurrPos();
      for (int j = i + 1; j < num2; ++j){
	Vec v = mergeParticleVec[j]->getCurrPos();
	if ( ( vfabs(v - u) < particleVec[i]->getA() + mergeParticleVec[j]->getA())
	     && ( particleVec[i]->getType() !=  1 || mergeParticleVec[j]->getType() != 1  )      // not both are fixed particles
	     && ( particleVec[i]->getType() !=  5 || mergeParticleVec[j]->getType() != 5  )      // not both are free boundary particles
	     && ( particleVec[i]->getType() != 10 || mergeParticleVec[j]->getType() != 10 )  ) { // not both are ghost particles
	  Contact tmpContact(particleVec[i], mergeParticleVec[j]); // a local and temparory object
	  ++possContactNum;
#ifdef TIME_PROFILE
	  gettimeofday(&time_r1, NULL); 
#endif
	  if(tmpContact.isOverlapped())
	    contactVec.push_back(tmpContact);    // containers use value semantics, so a "copy" is pushed back.
#ifdef TIME_PROFILE
	  gettimeofday(&time_r2, NULL); 
	  time_r += timediffsec(time_r1, time_r2);
#endif
	}
      }
    }	
    
#ifdef TIME_PROFILE
    gettimeofday(&time_p2, NULL);
    debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2) << std::setw(OWID) << "isOverlapped=" << std::setw(OWID) << time_r; 
#endif
    
    actualContactNum = contactVec.size();
  }

  else if (ompThreads > 1) { // openmp implementation: various loop scheduling - (static), (static,1), (dynamic), (dynamic,1)
    contactVec.clear();
    int possContact = 0;
    
#ifdef TIME_PROFILE
    gettimeofday(&time_p1, NULL); 
#endif
    
    int i, j;
    Vec u, v;
    int num1 = particleVec.size();
    int num2 = mergeParticleVec.size();  
    int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
    
#pragma omp parallel for num_threads(ompThreads) private(i, j, u, v) shared(num1, num2) reduction(+: possContact) schedule(dynamic)
    for (i = 0; i < num1; ++i) { 
      u = particleVec[i]->getCurrPos();
      for (j = i + 1; j < num2; ++j) {
	v = mergeParticleVec[j]->getCurrPos();
	if ( ( vfabs(v - u) < particleVec[i]->getA() + mergeParticleVec[j]->getA() )
	     && ( particleVec[i]->getType() !=  1 || mergeParticleVec[j]->getType() != 1  )      // not both are fixed particles
	     && ( particleVec[i]->getType() !=  5 || mergeParticleVec[j]->getType() != 5  )      // not both are free boundary particles
	     && ( particleVec[i]->getType() != 10 || mergeParticleVec[j]->getType() != 10 )  ) { // not both are ghost particles
	  Contact tmpContact(particleVec[i], mergeParticleVec[j]); // a local and temparory object
	  ++possContact;
	  if(tmpContact.isOverlapped())
#pragma omp critical
	    contactVec.push_back(tmpContact);    // containers use value semantics, so a "copy" is pushed back.
	}
      }
    }
    
#ifdef TIME_PROFILE
    gettimeofday(&time_p2, NULL);
    debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2); 
#endif
    possContactNum   = possContact;  
    actualContactNum = contactVec.size();
  } // end of openmp implementation

}


void Assembly::internalForce(){
  REAL avgNormal, avgShear;
  
  int contactNum = contactVec.size();
  if(contactNum == 0){
    avgNormal = 0;
    avgShear = 0;
  }
  else{
    for (std::vector<Contact>::iterator it = contactVec.begin(); it != contactVec.end(); ++it)
      it->checkinPrevTgt(contactTgtVec); // checkin previous tangential force and displacment    
    
#ifdef TIME_PROFILE
    gettimeofday(&time_p1,NULL); 
#endif 

    contactTgtVec.clear(); // contactTgtVec must be cleared before filling in new values.
    for (std::vector<Contact>::iterator it = contactVec.begin(); it != contactVec.end(); ++ it){
      it->contactForce();             // cannot be parallelized as it may change a particle's force simultaneously.
      it->checkoutTgt(contactTgtVec); // checkout current tangential force and displacment
      avgNormal += it->getNormalForce();
      avgShear += it->getTgtForce();
    }
    avgNormal /= contactNum;
    avgShear /= contactNum;
    
#ifdef TIME_PROFILE
    gettimeofday(&time_p2,NULL);
    debugInf << std::setw(OWID) << "internalForce=" << std::setw(OWID) << timediffsec(time_p1, time_p2) << std::endl; 
#endif
    
  }
}


void Assembly::updateParticle() {
  for(std::vector<Particle*>::iterator it = particleVec.begin(); it != particleVec.end(); ++it)
    (*it)->update();
}


void Assembly::clearContactForce() {
  for(std::vector<Particle*>::iterator it = particleVec.begin(); it != particleVec.end(); ++it)
    (*it)->clearContactForce();
}


void Assembly::findBdryContact() {
  for(std::vector<Boundary*>::iterator rt = boundaryVec.begin(); rt != boundaryVec.end(); ++rt)
    (*rt)->findBdryContact(particleVec);
}


void Assembly::boundaryForce() {
  for(std::vector<Boundary*>::iterator rt = boundaryVec.begin(); rt != boundaryVec.end(); ++rt)
    (*rt)->boundaryForce(boundaryTgtMap);
}


void Assembly::printContact(char *str) const
{
  char csuf[10];
  combineString(csuf, ".p", mpiRank, 5);
  strcat(str, csuf);

  std::ofstream ofs(str);
  if(!ofs) { std::cout << "stream error: printContact" << std::endl; exit(-1); }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  
  ofs << std::setw(OWID) << contactVec.size() << std::endl;
  ofs << std::setw(OWID) << "ptcl_1"
      << std::setw(OWID) << "ptcl_2"
      << std::setw(OWID) << "point1_x"
      << std::setw(OWID) << "point1_y"
      << std::setw(OWID) << "point1_z"
      << std::setw(OWID) << "point2_x"
      << std::setw(OWID) << "point2_y"
      << std::setw(OWID) << "point2_z"
      << std::setw(OWID) << "radius_1"
      << std::setw(OWID) << "radius_2"
      << std::setw(OWID) << "penetration"
      << std::setw(OWID) << "tangt_disp"
      << std::setw(OWID) << "contact_radius"
      << std::setw(OWID) << "R0"
      << std::setw(OWID) << "E0"
      << std::setw(OWID) << "normal_force"
      << std::setw(OWID) << "tangt_force"
      << std::setw(OWID) << "contact_x"
      << std::setw(OWID) << "contact_y"
      << std::setw(OWID) << "contact_z"
      << std::setw(OWID) << "normal_x"
      << std::setw(OWID) << "normal_y"
      << std::setw(OWID) << "normal_z"
      << std::setw(OWID) << "tangt_x"
      << std::setw(OWID) << "tangt_y"
      << std::setw(OWID) << "tangt_z"
      << std::setw(OWID) << "vibra_t_step"
      << std::setw(OWID) << "impact_t_step"
      << std::endl;
  
  std::vector<Contact>::const_iterator it;
  for (it = contactVec.begin(); it != contactVec.end(); ++it)
    ofs << std::setw(OWID) << it->getP1()->getId()
	<< std::setw(OWID) << it->getP2()->getId()
	<< std::setw(OWID) << it->getPoint1().getX()
	<< std::setw(OWID) << it->getPoint1().getY()
	<< std::setw(OWID) << it->getPoint1().getZ()
	<< std::setw(OWID) << it->getPoint2().getX()
	<< std::setw(OWID) << it->getPoint2().getY()
	<< std::setw(OWID) << it->getPoint2().getZ()
	<< std::setw(OWID) << it->getRadius1()
	<< std::setw(OWID) << it->getRadius2()
	<< std::setw(OWID) << it->getPenetration()
	<< std::setw(OWID) << it->getTgtDisp()
	<< std::setw(OWID) << it->getContactRadius()
	<< std::setw(OWID) << it->getR0()
	<< std::setw(OWID) << it->getE0()
	<< std::setw(OWID) << it->getNormalForce()
	<< std::setw(OWID) << it->getTgtForce()
	<< std::setw(OWID) << ( it->getPoint1().getX() + it->getPoint2().getX() )/2
	<< std::setw(OWID) << ( it->getPoint1().getY() + it->getPoint2().getY() )/2
	<< std::setw(OWID) << ( it->getPoint1().getZ() + it->getPoint2().getZ() )/2
	<< std::setw(OWID) << it->normalForceVec().getX()
	<< std::setw(OWID) << it->normalForceVec().getY()
	<< std::setw(OWID) << it->normalForceVec().getZ()
	<< std::setw(OWID) << it->tgtForceVec().getX()
	<< std::setw(OWID) << it->tgtForceVec().getY()
	<< std::setw(OWID) << it->tgtForceVec().getZ()
	<< std::setw(OWID) << it->getVibraTimeStep()
	<< std::setw(OWID) << it->getImpactTimeStep()
	<< std::endl;
  ofs.close();
}


// rule out
/*
void Assembly::plotCavity(const char *str) const {
  std::ofstream ofs(str);
  if(!ofs) { std::cout << "stream error: plotCavity" << std::endl; exit(-1); }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  REAL x1,x2,y1,y2,z1,z2,x0,y0,z0;
  x1 = cavity.getMinCorner().getX();
  y1 = cavity.getMinCorner().getY();
  z1 = cavity.getMinCorner().getZ();
  x2 = cavity.getMaxCorner().getX();
  y2 = cavity.getMaxCorner().getY();
  z2 = cavity.getMaxCorner().getZ();

  ofs << "ZONE N=8, E=1, DATAPACKING=POINT, ZONETYPE=FEBRICK" << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y1 << std::setw(OWID) << z1 << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z1 << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y2 << std::setw(OWID) << z1 << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z1 << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y1 << std::setw(OWID) << z2 << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z2 << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y2 << std::setw(OWID) << z2 << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z2 << std::endl;
  ofs << "1 2 3 4 5 6 7 8" << std::endl;

  ofs.close();
}

void Assembly::plotSpring(const char *str) const {
  std::ofstream ofs(str);
  if(!ofs) { std::cout << "stream error: plotSpring" << std::endl; exit(-1); }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  int totalMemParticle = 0;
  for (int i = 0; i < memBoundary.size(); ++i) 
    for (int j = 0; j < memBoundary[i].size(); ++j) 
      for (int k = 0; k < memBoundary[i][j].size(); ++k) 
	++totalMemParticle;
  int totalSpring = springVec.size();
  ofs << "ZONE N=" << totalMemParticle << ", E=" << totalSpring << ", DATAPACKING=POINT, ZONETYPE=FELINESEG" << std::endl;
  Particle *pt = NULL;
  Vec vt;
  for (int i = 0; i < memBoundary.size(); ++i) 
    for (int j = 0; j < memBoundary[i].size(); ++j) 
      for (int k = 0; k < memBoundary[i][j].size(); ++k) {
	pt = memBoundary[i][j][k]; 
	vt = pt->getCurrPos();
	ofs << std::setw(OWID) << vt.getX() << std::setw(OWID) << vt.getY() << std::setw(OWID) << vt.getZ() << std::endl;
      }
  for (int i = 0; i < springVec.size(); ++i) {
    ofs << std::setw(OWID) << springVec[i]->getParticleId1() - trimHistoryNum  << std::setw(OWID) << springVec[i]->getParticleId2() - trimHistoryNum << std::endl;
  }

  ofs.close();
}  

void Assembly::printMemParticle(const char *str) const  {
  std::ofstream ofs(str);
  if(!ofs) { std::cout << "stream error: printMemParticle" << std::endl; exit(-1); }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  
  int totalMemParticle = 0;
  for (int i = 0; i < memBoundary.size(); ++i) 
    for (int j = 0; j < memBoundary[i].size(); ++j) 
      for (int k = 0; k < memBoundary[i][j].size(); ++k) 
	++totalMemParticle;
  
  ofs << std::setw(OWID) << totalMemParticle << std::setw(OWID) << 1 << std::endl;
  ofs << std::setw(OWID) << container.getCenter().getX()
      << std::setw(OWID) << container.getCenter().getY()
      << std::setw(OWID) << container.getCenter().getZ()
      << std::setw(OWID) << container.getDimx()
      << std::setw(OWID) << container.getDimy()
      << std::setw(OWID) << container.getDimz() << std::endl;
  
  ofs << std::setw(OWID) << "ID"
      << std::setw(OWID) << "type"
      << std::setw(OWID) << "radius_a"
      << std::setw(OWID) << "radius_b"
      << std::setw(OWID) << "radius_c"
      << std::setw(OWID) << "position_x"
      << std::setw(OWID) << "position_y"
      << std::setw(OWID) << "position_z"
      << std::setw(OWID) << "axle_a_x"
      << std::setw(OWID) << "axle_a_y"
      << std::setw(OWID) << "axle_a_z"
      << std::setw(OWID) << "axle_b_x"
      << std::setw(OWID) << "axle_b_y"
      << std::setw(OWID) << "axle_b_z"
      << std::setw(OWID) << "axle_c_x"
      << std::setw(OWID) << "axle_c_y"
      << std::setw(OWID) << "axle_c_z"
      << std::setw(OWID) << "velocity_x"
      << std::setw(OWID) << "velocity_y"
      << std::setw(OWID) << "velocity_z"
      << std::setw(OWID) << "omga_x"
      << std::setw(OWID) << "omga_y"
      << std::setw(OWID) << "omga_z"
      << std::setw(OWID) << "force_x"
      << std::setw(OWID) << "force_y"
      << std::setw(OWID) << "force_z"
      << std::setw(OWID) << "moment_x"
      << std::setw(OWID) << "moment_y"
      << std::setw(OWID) << "moment_z"
      << std::endl;
  
  Particle *it = NULL;
  Vec vObj;
  for (int i = 0; i < memBoundary.size(); ++i) 
    for (int j = 0; j < memBoundary[i].size(); ++j) 
      for (int k = 0; k < memBoundary[i][j].size(); ++k) {
	it = memBoundary[i][j][k];
	ofs << std::setw(OWID) << it->getId()
	    << std::setw(OWID) << it->getType()
	    << std::setw(OWID) << it->getA()
	    << std::setw(OWID) << it->getB()
	    << std::setw(OWID) << it->getC();
	
	vObj=it->getCurrPos();
	ofs << std::setw(OWID) << vObj.getX()
	    << std::setw(OWID) << vObj.getY()
	    << std::setw(OWID) << vObj.getZ();
	
	vObj=it->getCurrDirecA();
	ofs << std::setw(OWID) << vObj.getX()
	    << std::setw(OWID) << vObj.getY()
	    << std::setw(OWID) << vObj.getZ();
	
	vObj=it->getCurrDirecB();
	ofs << std::setw(OWID) << vObj.getX()
	    << std::setw(OWID) << vObj.getY()
	    << std::setw(OWID) << vObj.getZ();
	
	vObj=it->getCurrDirecC();
	ofs << std::setw(OWID) << vObj.getX()
	    << std::setw(OWID) << vObj.getY()
	    << std::setw(OWID) << vObj.getZ();
	
	vObj=it->getCurrVeloc();
	ofs << std::setw(OWID) << vObj.getX()
	    << std::setw(OWID) << vObj.getY()
	    << std::setw(OWID) << vObj.getZ();
	
	vObj=it->getCurrOmga();
	ofs << std::setw(OWID) << vObj.getX()
	    << std::setw(OWID) << vObj.getY()
	    << std::setw(OWID) << vObj.getZ();
	
	vObj=it->getForce();
	ofs << std::setw(OWID) << vObj.getX()
	    << std::setw(OWID) << vObj.getY()
	    << std::setw(OWID) << vObj.getZ();
	
	vObj=it->getMoment();
	ofs << std::setw(OWID) << vObj.getX()
	    << std::setw(OWID) << vObj.getY()
	      << std::setw(OWID) << vObj.getZ() << std::endl;
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
  std::vector<Particle*> vec1d;  // 1-dimension
  std::vector< std::vector<Particle*>  > vec2d; // 2-dimension
  REAL in, out, tmp;
  REAL x1_in, x1_out, x2_in, x2_out;
  REAL y1_in, y1_out, y2_in, y2_out;
  REAL z1_in, z1_out, z2_in, z2_out;

  // surface x1
  vec2d = memBoundary[0];
  in = vec2d[0][0]->getCurrPos().getX();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPos().getX();
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
  in = vec2d[0][0]->getCurrPos().getX();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPos().getX();
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
  in = vec2d[0][0]->getCurrPos().getY();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPos().getY();
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
  in = vec2d[0][0]->getCurrPos().getY();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPos().getY();
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
  in = vec2d[0][0]->getCurrPos().getZ();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPos().getZ();
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
  in = vec2d[0][0]->getCurrPos().getZ();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPos().getZ();
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

// OPENMP_IMPL: 
// 0: implementation 0, ts partitions, based on linked list
// 1: implementation 1, ts partitions, based on vector
// 2: implementation 2, no partition, each thread leaps by ts until completed
// 3: implementation 3, no partition, each thread leaps by ts until num/2 and handles two particles.
// 4: implementation 4, no partition, parallel for, various loop scheduling: (static), (static,1), (dynamic), (dynamic,1)

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
  
#pragma omp parallel num_threads(nThreads) private(tid, ts, tnum, it, pt, i, j, u, v) shared(num) reduction(+: possContact)
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
	  contact<Particle> tmpContact(*it, *pt); // a local and temparory object
	  ++possContact;
	  if(tmpContact.isOverlapped())
#pragma omp critical
	    contactVec.push_back(tmpContact);    // containers use value semantics, so a "copy" is pushed back.
	}
      }
    }
  }
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2); 
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
  
#pragma omp parallel num_threads(nThreads) private(tid, ts, start, end, i, j, u, v) shared(num) reduction(+: possContact)
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
	  contact<Particle> tmpContact(particleVec[i], particleVec[j]); // a local and temparory object
	  ++possContact;
	  if(tmpContact.isOverlapped())
#pragma omp critical
	    contactVec.push_back(tmpContact);    // containers use value semantics, so a "copy" is pushed back.
	}
      }
    }
  }
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf <<  std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2); 
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
  
#pragma omp parallel num_threads(nThreads) private(tid, ts, i, j, u, v) shared(num) reduction(+: possContact)
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
	  contact<Particle> tmpContact(particleVec[i], particleVec[j]); // a local and temparory object
	  ++possContact;
	  if(tmpContact.isOverlapped())
#pragma omp critical
	    contactVec.push_back(tmpContact);    // containers use value semantics, so a "copy" is pushed back.
	}
      }
    }
  }
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2); 
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
  
#pragma omp parallel num_threads(nThreads) private(tid, ts, i, j, k, u, v) shared(num) reduction(+: possContact)
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
	    contact<Particle> tmpContact(particleVec[k], particleVec[j]); // a local and temparory object
	    ++possContact;
	    if(tmpContact.isOverlapped())
#pragma omp critical
	    contactVec.push_back(tmpContact);    // containers use value semantics, so a "copy" is pushed back.
	  }
	}
      }
    }
  }
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2); 
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
  
#pragma omp parallel for num_threads(nThreads) private(i, j, u, v) shared(num) reduction(+: possContact) schedule(dynamic)
  for (i = 0; i < num - 1; ++i) { 
    u = particleVec[i]->getCurrPos();
    for (j = i + 1; j < num; ++j) {
      v = particleVec[j]->getCurrPos();
      if (   ( vfabs(v-u) < particleVec[i]->getA() + particleVec[j]->getA() )
	     && ( particleVec[i]->getType() !=  1 || particleVec[j]->getType() != 1  )      // not both are fixed particles
	     && ( particleVec[i]->getType() !=  5 || particleVec[j]->getType() != 5  )      // not both are free boundary particles
	     && ( particleVec[i]->getType() != 10 || particleVec[j]->getType() != 10 )  ) { // not both are ghost particles
	contact<Particle> tmpContact(particleVec[i], particleVec[j]); // a local and temparory object
	++possContact;
	if(tmpContact.isOverlapped())
#pragma omp critical
	  contactVec.push_back(tmpContact);    // containers use value semantics, so a "copy" is pushed back.
      }
    }
  }
  
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2); 
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
    double time_r = 0; // time consumed in contact resolution, i.e., tmpContact.isOverlapped()
    gettimeofday(&time_p1,NULL); 
#endif
    
    int num1 = particleVec.size();  // particles inside container
    int num2 = mergeParticleVec.size(); // particles inside container (at front) + particles from neighboring blocks (at end)
    for (int i = 0; i < num1 - 1; ++i) {
      Vec u = particleVec[i]->getCurrPos();
      for (int j = i + 1; j < num2; ++j){
	Vec v = mergeParticleVec[j]->getCurrPos();
	if (   ( vfabs(v - u) < particleVec[i]->getA() + mergeParticleVec[j]->getA())
	    && ( particleVec[i]->getType() !=  1 || mergeParticleVec[j]->getType() != 1  )      // not both are fixed particles
	    && ( particleVec[i]->getType() !=  5 || mergeParticleVec[j]->getType() != 5  )      // not both are free boundary particles
	    && ( particleVec[i]->getType() != 10 || mergeParticleVec[j]->getType() != 10 )  ) { // not both are ghost particles
	  Contact tmpContact(particleVec[i], mergeParticleVec[j]); // a local and temparory object
	  ++possContactNum;
#ifdef TIME_PROFILE
	  gettimeofday(&time_r1,NULL); 
#endif
	  if(tmpContact.isOverlapped())
	    contactVec.push_back(tmpContact);    // containers use value semantics, so a "copy" is pushed back.
#ifdef TIME_PROFILE
	  gettimeofday(&time_r2,NULL); 
	  time_r += timediffsec(time_r1, time_r2);
#endif
	}
      }
    }	
    
#ifdef TIME_PROFILE
    gettimeofday(&time_p2,NULL);
    debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2) << std::setw(OWID) << "isOverlapped=" << std::setw(OWID) << time_r; 
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
  REAL x0 = minCorner.getX();
  REAL y0 = minCorner.getY();
  REAL z0 = minCorner.getZ();
  
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
  typedef std::pair<bool, std::vector<Particle*> > cellT;
  std::vector< std::vector< std::vector < cellT > > > cellVec;
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
	  if (center.getX() >= x1 && center.getX() < x2 &&
	      center.getY() >= y1 && center.getY() < y2 &&
	      center.getZ() >= z1 && center.getZ() < z2)
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
	    //std::cout <<  i << " " << j << " " << k << " " << "m n size=" << m << " " << n << " " <<  cellVec[i][j][k].size() << std::endl;
	    pt = cellVec[i][j][k].second[n];
	    v  = pt->getCurrPos();
	    if ( ( vfabs(u-v) < it->getA() + pt->getA() )  &&
		 ( it->getType() !=  1 || pt->getType() != 1 ) &&   // not both are fixed particles
		 ( it->getType() !=  5 || pt->getType() != 5 ) &&   // not both are free boundary particles
		 ( it->getType() != 10 || pt->getType() != 10)  ) { // not both are ghost particles
	      contact<Particle> tmpContact(it, pt); // a local and temparory object
	      ++possContactNum;
#ifdef TIME_PROFILE
		gettimeofday(&time_r1,NULL); 
#endif
	      if(tmpContact.isOverlapped())
		contactVec.push_back(tmpContact);   // containers use value semantics, so a "copy" is pushed back.
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
	      //std::cout << "i j k m ncell ci cj ck size contacts= " << i << " " << j << " " << k << " " << m  << " " << ncell << " " << ci << " " << cj << " " << ck << " " << cellVec[ci][cj][ck].second.size() << " "  << contactVec.size() << std::endl;
	      std::vector<Particle*> vt = cellVec[ci][cj][ck].second;
	      for (int n = 0; n < vt.size(); ++n) {
		pt = vt[n];
		v  = pt->getCurrPos();
		if ( ( vfabs(u-v) < it->getA() + pt->getA() )  &&
		     ( it->getType() !=  1 || pt->getType() != 1 ) &&   // not both are fixed particles
		     ( it->getType() !=  5 || pt->getType() != 5 ) &&   // not both are free boundary particles
		     ( it->getType() != 10 || pt->getType() != 10)  ) { // not both are ghost particles
		  contact<Particle> tmpContact(it, pt); // a local and temparory object
		  ++possContactNum;
#ifdef TIME_PROFILE
		gettimeofday(&time_r1,NULL); 
#endif
		  if(tmpContact.isOverlapped())
		    contactVec.push_back(tmpContact);   // containers use value semantics, so a "copy" is pushed back.
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
  debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2) << std::setw(OWID) << "isOverlapped=" << std::setw(OWID) << time_r; 
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
	for (std::vector<Contact>::const_iterator it=contactVec.begin();it!=contactVec.end();++it)
	    pene += it->getPenetration(); 
	return pene/totalcntct;
    }
}

REAL Assembly::getVibraTimeStep() const {
    int totalcntct = contactVec.size();
    if (totalcntct == 0)
	return 0;
    else {
	std::vector<Contact>::const_iterator it=contactVec.begin();
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
	std::vector<Contact>::const_iterator it=contactVec.begin();
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
		if ((*it)->getCurrPos().getZ() > (*jt)->getCurrPos().getZ())
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
	    val = (*it)->getForce().getZ();
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
	    val = (*it)->getCurrPos().getZ()-(*it)->getA();
	    break;
	}
    return val;
}

REAL Assembly::ellipPilePeneVol() {
    REAL val=0;
    if (getTopFreeParticlePosition().getZ()-ellipPileTipZ() <= 0)
	val=0;
    else{
	// low: a signed number as lower limit for volumetric integration
	REAL low=ellipPileTipZ() + ellipPileDimn().getX() - getTopFreeParticlePosition().getZ(); 
	REAL lowint=low-pow(low,3)/3.0/pow(ellipPileDimn().getX(),2);
	val = Pi * ellipPileDimn().getY() * ellipPileDimn().getZ()
	      *(2.0/3*ellipPileDimn().getX()-lowint);
    }
    return val;
}

void Assembly::ellipPileUpdate(){
  for(std::vector<Particle*>::iterator it=particleVec.begin();it!=particleVec.end();++it){
    if ((*it)->getType()==3) {
      (*it)->setCurrVeloc(Vec(0, 0, -pileRate));
      (*it)->setCurrPos( (*it)->getPrevPos() + (*it)->getCurrVeloc() * timeStep);
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

void Assembly::readCavityBoundary(const char *str) {
  std::ifstream ifs(str);
  if(!ifs) { std::cout << "stream error: readCavityBoundary" << std::endl; exit(-1); }  

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


void Assembly::printCavityBoundary(const char *str) const {
  std::ofstream ofs(str);
  if(!ofs) { std::cout << "stream error: printCavityBoundary" << std::endl; exit(-1); }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  
  ofs << std::setw(OWID) << cavityBoundaryVec.size() << std::endl;
  std::vector<Boundary*>::const_iterator rt;
  for(rt = cavityBoundaryVec.begin(); rt != cavityBoundaryVec.end(); ++rt)
    (*rt)->display(ofs);
  ofs << std::endl;
  
  ofs.close();
}



void Assembly::findCavityContact(){
  std::vector<Boundary*>::iterator rt;
  for(rt = cavityBoundaryVec.begin(); rt != cavityBoundaryVec.end(); ++rt)
    (*rt)->findBdryContact(allParticleVec);
}



void Assembly::boundaryForce(REAL penetr[],int cntnum[]){
  std::vector<Boundary*>::iterator rt;
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
  std::vector<Boundary*>::iterator rt;
  for(rt = cavityBoundaryVec.begin(); rt != cavityBoundaryVec.end(); ++rt)
    (*rt)->boundaryForce(boundaryTgtMap);
}

Vec Assembly::getNormalForce(int bdry) const{
    std::vector<Boundary*>::const_iterator it;
    for(it=boundaryVec.begin();it!=boundaryVec.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getNormalForce();
    }
    return 0;
}

Vec Assembly::getShearForce(int bdry) const{
    std::vector<Boundary*>::const_iterator it;
    for(it=boundaryVec.begin();it!=boundaryVec.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getShearForce();
    }
    return 0;
}

REAL Assembly::getAvgNormal(int bdry) const{
    std::vector<Boundary*>::const_iterator it;
    for(it=boundaryVec.begin();it!=boundaryVec.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getAvgNormal();
    }
    return 0;
}

Vec Assembly::getApt(int bdry) const{
    std::vector<Boundary*>::const_iterator it;
    for(it=boundaryVec.begin();it!=boundaryVec.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getApt();
    }
    return 0;
}


Vec Assembly::getDirc(int bdry) const{
    std::vector<Boundary*>::const_iterator it;
    for(it=boundaryVec.begin();it!=boundaryVec.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getDirc();
    }
    return 0;
}

REAL Assembly::getArea(int n) const{
    std::vector<Boundary*>::const_iterator it;
    for(it=boundaryVec.begin();it!=boundaryVec.end();++it){
	if((*it)->getBdryID()==n)
	    return (*it)->area;
    }
    return 0;
}

void Assembly::setArea(int n, REAL a){
    std::vector<Boundary*>::iterator it;
    for(it=boundaryVec.begin();it!=boundaryVec.end();++it){
	if((*it)->getBdryID()==n)
	    (*it)->area=a;
    }
}

REAL Assembly::getAvgPressure() const{
    std::vector<Boundary*>::const_iterator rt;
    REAL avgpres=0;
    for(rt=boundaryVec.begin();rt!=boundaryVec.end();++rt)
	avgpres+=vfabs((*rt)->getNormalForce())/(*rt)->getArea();
    return avgpres/=boundaryVec.size();
}

// only update CoefOfLimits[0] for specified boundaries
void Assembly::updateBoundary(int bn[], UPDATECTL rbctl[], int num){
    for(int i=0;i<num;i++){
	for(std::vector<Boundary*>::iterator rt=boundaryVec.begin();rt!=boundaryVec.end();++rt){
	    if((*rt)->getBdryID()==bn[i]){
		(*rt)->update(rbctl[i]);
		break;
	    }
	}
    }
}

// update CoefOfLimits[1,2,3,4] for all 6 boundaries
void Assembly::updateBoundary6(){
    for(std::vector<Boundary*>::iterator rt=boundaryVec.begin();rt!=boundaryVec.end();++rt){
	if((*rt)->getBdryID()==1 || (*rt)->getBdryID()==3){
	    for(std::vector<Boundary*>::iterator lt=boundaryVec.begin();lt!=boundaryVec.end();++lt){
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
	    for(std::vector<Boundary*>::iterator lt=boundaryVec.begin();lt!=boundaryVec.end();++lt){
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
	    for(std::vector<Boundary*>::iterator lt=boundaryVec.begin();lt!=boundaryVec.end();++lt){
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
			      const char *inibdryfile,
			      const char *ParticleFile, 
			      const char *contactfile,
			      const char *progressfile, 
			      const char *debugfile)
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
			     const char *inibdryfile,
			     const char *ParticleFile, 
			     const char *contactfile,
			     const char *progressfile, 
			     const char *debugfile)
{
  // pre_1: open streams for output.
  progressinf.open(progressfile); 
  if(!progressinf) { std::cout << "stream error: angleOfRepose" << std::endl; exit(-1); }
  progressinf.setf(std::ios::scientific, std::ios::floatfield);
  progressinf.precision(OPREC);
  progressinf << std::setw(OWID) << "iteration"
	      << std::setw(OWID) << "poss_contact"
	      << std::setw(OWID) << "actual_contact"
	      << std::setw(OWID) << "penetration"
	      << std::setw(OWID) << "avg_normal"
	      << std::setw(OWID) << "avg_tangt"
	      << std::setw(OWID) << "avg_velocity"
	      << std::setw(OWID) << "avg_omga"
	      << std::setw(OWID) << "avg_force"
	      << std::setw(OWID) << "avg_moment"
	      << std::setw(OWID) << "trans_energy"
	      << std::setw(OWID) << "rotat_energy"
	      << std::setw(OWID) << "kinet_energy"
	      << std::setw(OWID) << "poten_energy"
	      << std::setw(OWID) << "total_energy"
	      << std::setw(OWID) << "void_ratio"
	      << std::setw(OWID) << "porosity"
	      << std::setw(OWID) << "coord_number"
	      << std::setw(OWID) << "density"
	      << std::setw(OWID) << "sigma_y1"
	      << std::setw(OWID) << "sigma_y2"
	      << std::setw(OWID) << "sigma_x1"
	      << std::setw(OWID) << "sigma_x2"
	      << std::setw(OWID) << "sigma_z1"
	      << std::setw(OWID) << "sigma_z2"
	      << std::setw(OWID) << "mean_stress"
	      << std::setw(OWID) << "dimx"
	      << std::setw(OWID) << "dimy"
	      << std::setw(OWID) << "dimz"
	      << std::setw(OWID) << "volume"
	      << std::setw(OWID) << "epsilon_x"
	      << std::setw(OWID) << "epsilon_y"
	      << std::setw(OWID) << "epsilon_z"
	      << std::setw(OWID) << "epsilon_v"
	      << std::setw(OWID) << "vibra_t_step"
	      << std::setw(OWID) << "impact_t_step"
	      << std::setw(OWID) << "wall_time" << std::endl;
  
  debugInf.open(debugfile);
  if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1); }
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
  REAL z0 = container.getMinCorner().getZ();
  std::vector<Particle*> lastPtcls;
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
	  newPtcl = new Particle(particleNum+1, 0, Vec(0, 0, zCurr + maxDiameter * i), gradation, young, poisson);
	  particleVec.push_back(newPtcl);
	  ++particleNum;
	  lastPtcls.push_back(newPtcl);
	  
	  newPtcl = new Particle(particleNum+1, 0, Vec(maxDiameter, 0, zCurr + maxDiameter * i), gradation, young, poisson);
	  particleVec.push_back(newPtcl);
	  ++particleNum;
	  //lastPtcls.push_back(newPtcl);
	  
	  newPtcl = new Particle(particleNum+1, 0, Vec(-maxDiameter, 0, zCurr + maxDiameter * i), gradation, young, poisson);
	  particleVec.push_back(newPtcl);
	  ++particleNum;
	  //lastPtcls.push_back(newPtcl);
	  
	  newPtcl = new Particle(particleNum+1, 0, Vec(0, maxDiameter, zCurr + maxDiameter * i), gradation, young, poisson);
	  particleVec.push_back(newPtcl);
	  ++particleNum;
	  //lastPtcls.push_back(newPtcl);
	  
	  newPtcl = new Particle(particleNum+1, 0, Vec(0, -maxDiameter, zCurr + maxDiameter * i), gradation, young, poisson);
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
	  zCurr = getPtclMaxZ(allParticleVec) + maxDiameter;

	  for ( int i = 0; i != layers; ++i) {
	    newPtcl = new Particle(particleNum+1, 0, Vec(0, 0, zCurr + maxDiameter * i), gradation, young, poisson);
	    particleVec.push_back(newPtcl);
	    ++particleNum;
	    lastPtcls.push_back(newPtcl);
	    
	    newPtcl = new Particle(particleNum+1, 0, Vec(maxDiameter, 0, zCurr + maxDiameter * i), gradation, young, poisson);
	    particleVec.push_back(newPtcl);
	    ++particleNum;
	    //lastPtcls.push_back(newPtcl);
	    
	    newPtcl = new Particle(particleNum+1, 0, Vec(-maxDiameter, 0, zCurr + maxDiameter * i), gradation, young, poisson);
	    particleVec.push_back(newPtcl);
	    ++particleNum;
	    //lastPtcls.push_back(newPtcl);
	    
	    newPtcl = new Particle(particleNum+1, 0, Vec(0, maxDiameter, zCurr + maxDiameter * i), gradation, young, poisson);
	    particleVec.push_back(newPtcl);
	    ++particleNum;
	    //lastPtcls.push_back(newPtcl);
	    
	    newPtcl = new Particle(particleNum+1, 0, Vec(0, -maxDiameter, zCurr + maxDiameter * i), gradation, young, poisson);
	    particleVec.push_back(newPtcl);
	    ++particleNum;
	    //lastPtcls.push_back(newPtcl);	
	  }
	  toSnapshot = true;
	}
      }
      // 2. create possible boundary particles and contacts between particles.
      findContact();
      findBdryContact();
      
      // 3. set particle forces/moments as zero before each re-calculation,
      clearContactForce();	
      
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
	g_timeinf << std::setw(4) << stepsnum << " " << ctime(&timeStamp) << std::flush;
	++stepsnum;
	toSnapshot = false;
      }
      
      // 8. (2) output stress and strain info.
      if (iteration % interval == 0) {
	gettimeofday(&time_w2,NULL);
	REAL t1=getTransEnergy();
	REAL t2=getRotatEnergy();
	REAL t3=getPotenEnergy(-0.025);
	progressinf << std::setw(OWID) << iteration
		    << std::setw(OWID) << getPossContactNum()
		    << std::setw(OWID) << getActualContactNum()
		    << std::setw(OWID) << getAvgPenetration()
		    << std::setw(OWID) << avgNormal
		    << std::setw(OWID) << avgTangt
		    << std::setw(OWID) << getAvgVelocity() 
		    << std::setw(OWID) << getAvgOmga()
		    << std::setw(OWID) << getAvgForce()   
		    << std::setw(OWID) << getAvgMoment()
		    << std::setw(OWID) << t1
		    << std::setw(OWID) << t2
		    << std::setw(OWID) << (t1+t2)
		    << std::setw(OWID) << t3
		    << std::setw(OWID) << (t1+t2+t3)
		    << std::setw(OWID) << void_ratio
		    << std::setw(OWID) << void_ratio/(1+void_ratio)
		    << std::setw(OWID) << 2.0*(getActualContactNum()
				      +bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
				      +bdry_cntnum[4]+bdry_cntnum[6])/allParticleVec.size()
		    << std::setw(OWID) << "0"
		    << std::setw(OWID) << "0"
		    << std::setw(OWID) << "0"
		    << std::setw(OWID) << "0"
		    << std::setw(OWID) << "0"
		    << std::setw(OWID) << "0"
		    << std::setw(OWID) << "0"
		    << std::setw(OWID) << "0"
		    << std::setw(OWID) << "0"
		    << std::setw(OWID) << "0"
		    << std::setw(OWID) << "0"
		    << std::setw(OWID) << "0"
		    << std::setw(OWID) << "0"
		    << std::setw(OWID) << "0"
		    << std::setw(OWID) << "0"
		    << std::setw(OWID) << "0"
		    << std::setw(OWID) << getVibraTimeStep()
		    << std::setw(OWID) << getImpactTimeStep()
		    << std::setw(OWID) << timediffsec(time_w1,time_w2)
		        << std::endl;
      }
      
      // 7. loop break conditions.
      ++iteration;
      
    } while (particleNum < 2000); //( zCurr < container.getMaxCorner().getZ() );  //(++iteration < totalSteps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    g_timeinf << std::setw(4) << "end" << " " << ctime(&timeStamp) << std::flush;

    // post_2. close streams
    progressinf.close();
    debugInf.close();
}







void Assembly::scale_PtclBdry(int   totalSteps,  
			      int   snapNum,
			      int   interval,
			      REAL dimn,
			      REAL rsize,
			      const char *iniptclfile,   
			      const char *ParticleFile, 
			      const char *contactfile,
			      const char *progressfile, 
			      const char *debugfile)
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
			const char *iniptclfile,
			const char *initboundary,
			const char *ParticleFile,
			const char *contactfile,
			const char *progressfile,
			const char *debugfile)
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
			  const char *ParticleFile,
			  const char *cavParticleFile)
{
  if (toRebuild) readParticle(ParticleFile);
  trimHistoryNum = allParticleVec.size();

  REAL x1,x2,y1,y2,z1,z2,x0,y0,z0;
  x1 = cavity.getMinCorner().getX();
  y1 = cavity.getMinCorner().getY();
  z1 = cavity.getMinCorner().getZ();
  x2 = cavity.getMaxCorner().getX();
  y2 = cavity.getMaxCorner().getY();
  z2 = cavity.getMaxCorner().getZ();
  x0 = cavity.getCenter().getX();
  y0 = cavity.getCenter().getY();
  z0 = cavity.getCenter().getZ();
 
  std::vector<Particle*>::iterator itr;
  Vec center;
  REAL delta = gradation.getPtclMaxRadius();

  for (itr = particleVec.begin(); itr != particleVec.end(); ){
    center=(*itr)->getCurrPos();
    if(center.getX() + delta  >= x1 && center.getX() - delta <= x2 &&
       center.getY() + delta  >= y1 && center.getY() - delta <= y2 &&
       center.getZ() + delta  >= z1 && center.getZ() - delta <= z2 )
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
void Assembly::expandCavityParticle(bool toRebuild,
				    REAL percent,
				    const char *cavityptclfile,
				    const char *ParticleFile,
				    const char *newptclfile)
{
  if (toRebuild) readParticle(ParticleFile);
  trimHistoryNum = allParticleVec.size();

  REAL x1,x2,y1,y2,z1,z2;
  x1 = cavity.getMinCorner().getX();
  y1 = cavity.getMinCorner().getY();
  z1 = cavity.getMinCorner().getZ();
  x2 = cavity.getMaxCorner().getX();
  y2 = cavity.getMaxCorner().getY();
  z2 = cavity.getMaxCorner().getZ();
 
  std::vector<Particle*>::iterator itr;
  Vec center;

  int cavityPtclNum = 0;
  for (itr = particleVec.begin(); itr != particleVec.end(); ++itr ){
    center=(*itr)->getCurrPos();
    if(center.getX() > x1 && center.getX() < x2 &&
       center.getY() > y1 && center.getY() < y2 &&
       center.getZ() > z1 && center.getZ() < z2 )
      ++cavityPtclNum;
  }

  printCavityParticle(cavityPtclNum, cavityptclfile);

  for (itr = particleVec.begin(); itr != particleVec.end(); ++itr ){
    center=(*itr)->getCurrPos();
    if(center.getX() > x1 && center.getX() < x2 &&
       center.getY() > y1 && center.getY() < y2 &&
       center.getZ() > z1 && center.getZ() < z2 )
      (*itr)->expand(percent);
  }

  printParticle(newptclfile);
}


void Assembly::printCavityParticle(int total, const char *str) const {
  std::ofstream ofs(str);
  if(!ofs) { std::cout << "stream error: printCavityParticle" << std::endl; exit(-1); }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  ofs << std::setw(OWID) << total << std::setw(OWID) << 1 << std::endl;
  ofs << std::setw(OWID) << cavity.getCenter().getX()
      << std::setw(OWID) << cavity.getCenter().getY()
      << std::setw(OWID) << cavity.getCenter().getZ()
      << std::setw(OWID) << cavity.getDimx()
      << std::setw(OWID) << cavity.getDimy()
      << std::setw(OWID) << cavity.getDimz() << std::endl;
  
  ofs << std::setw(OWID) << "ID"
      << std::setw(OWID) << "type"
      << std::setw(OWID) << "radius_a"
      << std::setw(OWID) << "radius_b"
      << std::setw(OWID) << "radius_c"
      << std::setw(OWID) << "position_x"
      << std::setw(OWID) << "position_y"
      << std::setw(OWID) << "position_z"
      << std::setw(OWID) << "axle_a_x"
      << std::setw(OWID) << "axle_a_y"
      << std::setw(OWID) << "axle_a_z"
      << std::setw(OWID) << "axle_b_x"
      << std::setw(OWID) << "axle_b_y"
      << std::setw(OWID) << "axle_b_z"
      << std::setw(OWID) << "axle_c_x"
      << std::setw(OWID) << "axle_c_y"
      << std::setw(OWID) << "axle_c_z"
      << std::setw(OWID) << "velocity_x"
      << std::setw(OWID) << "velocity_y"
      << std::setw(OWID) << "velocity_z"
      << std::setw(OWID) << "omga_x"
      << std::setw(OWID) << "omga_y"
      << std::setw(OWID) << "omga_z"
      << std::setw(OWID) << "force_x"
      << std::setw(OWID) << "force_y"
      << std::setw(OWID) << "force_z"
      << std::setw(OWID) << "moment_x"
      << std::setw(OWID) << "moment_y"
      << std::setw(OWID) << "moment_z"
      << std::endl;

  REAL x1,x2,y1,y2,z1,z2;
  x1 = cavity.getMinCorner().getX();
  y1 = cavity.getMinCorner().getY();
  z1 = cavity.getMinCorner().getZ();
  x2 = cavity.getMaxCorner().getX();
  y2 = cavity.getMaxCorner().getY();
  z2 = cavity.getMaxCorner().getZ();
  
  Vec tmp;
  std::vector<Particle*>::const_iterator  it;
  for (it=particleVec.begin();it!=particleVec.end();++it)  {
    Vec center=(*it)->getCurrPos();
    if(center.getX() > x1 && center.getX() < x2 &&
       center.getY() > y1 && center.getY() < y2 &&
       center.getZ() > z1 && center.getZ() < z2 ) {

    ofs << std::setw(OWID) << (*it)->getId()
	<< std::setw(OWID) << (*it)->getType()
	<< std::setw(OWID) << (*it)->getA()
	<< std::setw(OWID) << (*it)->getB()
	<< std::setw(OWID) << (*it)->getC();
    
    tmp=(*it)->getCurrPos();
    ofs << std::setw(OWID) << tmp.getX()
	<< std::setw(OWID) << tmp.getY()
	<< std::setw(OWID) << tmp.getZ();
    
    tmp=(*it)->getCurrDirecA();
    ofs << std::setw(OWID) << tmp.getX()
	<< std::setw(OWID) << tmp.getY()
	<< std::setw(OWID) << tmp.getZ();
    
    tmp=(*it)->getCurrDirecB();
    ofs << std::setw(OWID) << tmp.getX()
	<< std::setw(OWID) << tmp.getY()
	<< std::setw(OWID) << tmp.getZ();
    
    tmp=(*it)->getCurrDirecC();
    ofs << std::setw(OWID) << tmp.getX()
	<< std::setw(OWID) << tmp.getY()
	<< std::setw(OWID) << tmp.getZ();
    
    tmp=(*it)->getCurrVeloc();
    ofs << std::setw(OWID) << tmp.getX()
	<< std::setw(OWID) << tmp.getY()
	<< std::setw(OWID) << tmp.getZ();
    
    tmp=(*it)->getCurrOmga();
    ofs << std::setw(OWID) << tmp.getX()
	<< std::setw(OWID) << tmp.getY()
	<< std::setw(OWID) << tmp.getZ();
    
    tmp=(*it)->getForce();
    ofs << std::setw(OWID) << tmp.getX()
	<< std::setw(OWID) << tmp.getY()
	<< std::setw(OWID) << tmp.getZ();
    
    tmp=(*it)->getMoment();
    ofs << std::setw(OWID) << tmp.getX()
	<< std::setw(OWID) << tmp.getY()
	<< std::setw(OWID) << tmp.getZ() << std::endl;
    }
  }
  
  ofs.close();
}


// bdrymum = 6 by default
// the variable existMaxID is important because cavity and container
// use the same boundaryTgtMap.
void Assembly::buildCavityBoundary(int existMaxId, const char *boundaryFile)
{
  std::ofstream ofs(boundaryFile);
  if(!ofs) { std::cout << "stream error: buildCavityBoundary" << std::endl; exit(-1); }

  REAL x1,x2,y1,y2,z1,z2,x0,y0,z0;
  x1 = cavity.getMinCorner().getX();
  y1 = cavity.getMinCorner().getY();
  z1 = cavity.getMinCorner().getZ();
  x2 = cavity.getMaxCorner().getX();
  y2 = cavity.getMaxCorner().getY();
  z2 = cavity.getMaxCorner().getZ();
  x0 = cavity.getCenter().getX();
  y0 = cavity.getCenter().getY();
  z0 = cavity.getCenter().getZ();

  int boundaryNum = 6;

  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs << std::setw(OWID) << 0
      << std::setw(OWID) << boundaryNum << std::endl << std::endl;

  // boundary 1
  ofs << std::setw(OWID) << 1 << std::endl
      << std::setw(OWID) << existMaxId + 1
      << std::setw(OWID) << 5
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << -1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0     
      << std::setw(OWID) << x2
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0     
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << -1
      << std::setw(OWID) << 0
      << std::setw(OWID) << x0    
      << std::setw(OWID) << y1
      << std::setw(OWID) << z0     
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y2
      << std::setw(OWID) << z0     
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << 1
      << std::setw(OWID) << x0     
      << std::setw(OWID) << y0    
      << std::setw(OWID) << z2 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 
      << std::setw(OWID) << -1
      << std::setw(OWID) << x0     
      << std::setw(OWID) << y0     
      << std::setw(OWID) << z1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl << std::endl
    
    // boundary 2
      << std::setw(OWID) << 1 << std::endl
      << std::setw(OWID) << existMaxId + 2
      << std::setw(OWID) << 5
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << -1
      << std::setw(OWID) << 0     
      << std::setw(OWID) << x0    
      << std::setw(OWID) << y2
      << std::setw(OWID) << z0     
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 0
      << std::setw(OWID) << x2
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0     
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << -1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << x1 
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0     
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << 1
      << std::setw(OWID) << x0     
      << std::setw(OWID) << y0    
      << std::setw(OWID) << z2 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 
      << std::setw(OWID) << -1
      << std::setw(OWID) << x0     
      << std::setw(OWID) << y0      
      << std::setw(OWID) << z1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl << std::endl
    
    // boundary 3
      << std::setw(OWID) << 1 << std::endl
      << std::setw(OWID) << existMaxId + 3
      << std::setw(OWID) << 5
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0     
      << std::setw(OWID) << x1
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << -1
      << std::setw(OWID) << 0
      << std::setw(OWID) << x0     
      << std::setw(OWID) << y1
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0  
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << x0       
      << std::setw(OWID) << y2
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << 1
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y0     
      << std::setw(OWID) << z2 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 
      << std::setw(OWID) << -1
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y0      
      << std::setw(OWID) << z1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl << std::endl
    
    // boundary 4
      << std::setw(OWID) << 1 << std::endl
      << std::setw(OWID) << existMaxId + 4
      << std::setw(OWID) << 5
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0     
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y1
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 0
      << std::setw(OWID) << x2
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << -1 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << x1 
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << 1
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y0     
      << std::setw(OWID) << z2 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 
      << std::setw(OWID) << -1
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y0      
      << std::setw(OWID) << z1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl << std::endl
    
    // boundary 5
      << std::setw(OWID) << 1 << std::endl
      << std::setw(OWID) << existMaxId + 5
      << std::setw(OWID) << 5
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 0
      << std::setw(OWID) << -1     
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y0
      << std::setw(OWID) << z2 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 0
      << std::setw(OWID) << x2
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << -1 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << x1 
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y2
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << -1
      << std::setw(OWID) << 0
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y1
      << std::setw(OWID) << z0
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl << std::endl
    
    // boundary 6
      << std::setw(OWID) << 1 << std::endl
      << std::setw(OWID) << existMaxId + 6
      << std::setw(OWID) << 5
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 1    
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y0
      << std::setw(OWID) << z1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 0
      << std::setw(OWID) << x2
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << -1 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << x1 
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y2
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << -1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y1
      << std::setw(OWID) << z0
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl << std::endl; 

  ofs.close();
}


// create boundary particles and springs connecting those boundary particles
void Assembly::createMemParticle(REAL rRadius,
				 bool toRebuild,
				 const char *ParticleFile,
				 const char *allParticle)
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
  REAL x1 = v1.getX();
  REAL y1 = v1.getY();
  REAL z1 = v1.getZ();
  REAL x2 = v2.getX();
  REAL y2 = v2.getY();
  REAL z2 = v2.getZ();
  REAL x0 = v0.getX();
  REAL y0 = v0.getY();
  REAL z0 = v0.getZ();

  Particle* newptcl = NULL;
  REAL x, y, z;
  
  std::vector<Particle*> vec1d;  // 1-dimension
  std::vector< std::vector<Particle*>  > vec2d; // 2-dimension
  Spring* newSpring = NULL;
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
	newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], memYoung);
	springVec.push_back(newSpring);
      }
      if (i + 1 < vec2d.size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], memYoung);
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
	newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], memYoung);
	springVec.push_back(newSpring);
      }
      if (i + 1 < vec2d.size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], memYoung);
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
	newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], memYoung);
	springVec.push_back(newSpring);
      }
      if (i + 1 < vec2d.size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], memYoung);
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
	newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], memYoung);
	springVec.push_back(newSpring);
      }
      if (i + 1 < vec2d.size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], memYoung);
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
	newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], memYoung);
	springVec.push_back(newSpring);
      }
      if (i + 1 < vec2d.size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], memYoung);
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
	newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], memYoung);
	springVec.push_back(newSpring);
      }
      if (i + 1 < vec2d.size() ) {
	newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], memYoung);
	springVec.push_back(newSpring);
      }
    }

  // membrane particles at the edges of each surface, for example,
  // x1y1 means particles on surface x1 connecting to particles on surface y1
  std::vector<Particle *> x1y1;
  std::vector<Particle *> x1y2;
  std::vector<Particle *> x1z1;
  std::vector<Particle *> x1z2;

  std::vector<Particle *> x2y1;
  std::vector<Particle *> x2y2;
  std::vector<Particle *> x2z1;
  std::vector<Particle *> x2z2;

  std::vector<Particle *> y1x1;
  std::vector<Particle *> y1x2;
  std::vector<Particle *> y1z1;
  std::vector<Particle *> y1z2;

  std::vector<Particle *> y2x1;
  std::vector<Particle *> y2x2;
  std::vector<Particle *> y2z1;
  std::vector<Particle *> y2z2;

  std::vector<Particle *> z1x1;
  std::vector<Particle *> z1x2;
  std::vector<Particle *> z1y1;
  std::vector<Particle *> z1y2;

  std::vector<Particle *> z2x1;
  std::vector<Particle *> z2x2;
  std::vector<Particle *> z2y1;
  std::vector<Particle *> z2y2;

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
    newSpring = new Spring(*x1y1[i], *y1x1[i], memYoung);
    springVec.push_back(newSpring);
  }
  assert(x1y2.size() == y2x1.size());
  for (int i = 0; i < x1y2.size(); ++i) {
    newSpring = new Spring(*x1y2[i], *y2x1[i], memYoung);
    springVec.push_back(newSpring);
  }
  assert(x1z1.size() == z1x1.size());
  for (int i = 0; i < x1z1.size(); ++i) {
    newSpring = new Spring(*x1z1[i], *z1x1[i], memYoung);
    springVec.push_back(newSpring);
  }
  assert(x1z2.size() == z2x1.size());
  for (int i = 0; i < x1z2.size(); ++i) {
    newSpring = new Spring(*x1z2[i], *z2x1[i], memYoung);
    springVec.push_back(newSpring);
  }
  // 4 edges on surface x2  
  assert(x2y1.size() == y1x2.size());
  for (int i = 0; i < x2y1.size(); ++i) {
    newSpring = new Spring(*x2y1[i], *y1x2[i], memYoung);
    springVec.push_back(newSpring);
  }
  assert(x2y2.size() == y2x2.size());
  for (int i = 0; i < x2y2.size(); ++i) {
    newSpring = new Spring(*x2y2[i], *y2x2[i], memYoung);
    springVec.push_back(newSpring);
  }
  assert(x2z1.size() == z1x2.size());
  for (int i = 0; i < x2z1.size(); ++i) {
    newSpring = new Spring(*x2z1[i], *z1x2[i], memYoung);
    springVec.push_back(newSpring);
  }
  assert(x2z2.size() == z2x2.size());
  for (int i = 0; i < x2z2.size(); ++i) {
    newSpring = new Spring(*x2z2[i], *z2x2[i], memYoung);
    springVec.push_back(newSpring);
  }
  // 2 edges on surface y1 
  assert(y1z1.size() == z1y1.size());
  for (int i = 0; i < y1z1.size(); ++i) {
    newSpring = new Spring(*y1z1[i], *z1y1[i], memYoung);
    springVec.push_back(newSpring);
  }
  assert(y1z2.size() == z2y1.size());
  for (int i = 0; i < y1z2.size(); ++i) {
    newSpring = new Spring(*y1z2[i], *z2y1[i], memYoung);
    springVec.push_back(newSpring);
  }
  // 2 edges on surface y2
  assert(y2z1.size() == z1y2.size());
  for (int i = 0; i < y2z1.size(); ++i) {
    newSpring = new Spring(*y2z1[i], *z1y2[i], memYoung);
    springVec.push_back(newSpring);
  }
  assert(y2z2.size() == z2y2.size());
  for (int i = 0; i < y2z2.size(); ++i) {
    newSpring = new Spring(*y2z2[i], *z2y2[i], memYoung);
    springVec.push_back(newSpring);
  }

  printParticle(allParticle);

}


void Assembly::TrimPtclBdryByHeight(REAL height,
			    const char *iniptclfile,
			    const char *ParticleFile)
{
  readParticle(iniptclfile);
  
  std::vector<Particle*>::iterator itr;
  for (itr = particleVec.begin(); itr != particleVec.end(); ){
    if ( (*itr)->getType() == 1 ) { // 1-fixed
      Vec center=(*itr)->getCurrPos();
      if(center.getZ() > height)
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
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf.precision(OPREC);
    progressinf << std::setw(OWID) << "iteration"
	        << std::setw(OWID) << "poss_contact"
	        << std::setw(OWID) << "actual_contact"
	        << std::setw(OWID) << "penetration"
	        << std::setw(OWID) << "avg_normal"
	        << std::setw(OWID) << "avg_tangt"
	        << std::setw(OWID) << "avg_velocity"
	        << std::setw(OWID) << "avg_omga"
	        << std::setw(OWID) << "avg_force"
	        << std::setw(OWID) << "avg_moment"
	        << std::setw(OWID) << "trans_energy"
	        << std::setw(OWID) << "rotat_energy"
	        << std::setw(OWID) << "kinet_energy"
	        << std::setw(OWID) << "poten_energy"
	        << std::setw(OWID) << "total_energy"
	        << std::setw(OWID) << "void_ratio"
	        << std::setw(OWID) << "porosity"
	        << std::setw(OWID) << "coord_number"
	        << std::setw(OWID) << "density"
	        << std::setw(OWID) << "sigma_y1"
	        << std::setw(OWID) << "sigma_y2"
	        << std::setw(OWID) << "sigma_x1"
	        << std::setw(OWID) << "sigma_x2"
	        << std::setw(OWID) << "sigma_z1"
	        << std::setw(OWID) << "sigma_z2"
	        << std::setw(OWID) << "mean_stress"
	        << std::setw(OWID) << "dimx"
	        << std::setw(OWID) << "dimy"
	        << std::setw(OWID) << "dimz"
	        << std::setw(OWID) << "volume"
	        << std::setw(OWID) << "epsilon_x"
	        << std::setw(OWID) << "epsilon_y"
	        << std::setw(OWID) << "epsilon_z"
	        << std::setw(OWID) << "epsilon_v"
	        << std::setw(OWID) << "vibra_t_step"
	        << std::setw(OWID) << "impact_t_step"
	        << std::setw(OWID) << "wall_time" << std::endl;

    debugInf.open("dep_debug");
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1); }
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
	//std::cout << std::setw(OWID) << timediffsec(time_w1,time_w2);
	//gettimeofday(&time_w1,NULL);
        findBdryContact();

	// 2. set particles' forces/moments as zero before each re-calculation,
	clearContactForce();	

	// 3. calculate contact forces/moments and apply them to particles.
	internalForce(avgNormal, avgTangt);

	// 4. calculate boundary forces/moments and apply them to particles.
	boundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();
	//gettimeofday(&time_w2,NULL);
	//std::cout << std::setw(OWID) << timediffsec(time_w1,time_w2) << std::endl;

	// 6. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getZ() - getApt(6).getZ();
	l24=getApt(2).getY()-getApt(4).getY();
	l13=getApt(1).getX()-getApt(3).getX(); bulkVolume=l13*l24*l56;
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
	    g_timeinf << std::setw(4) << stepsnum << " " << ctime(&timeStamp) << std::flush;
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info.
	if (iteration % interval == 0) {
	    gettimeofday(&time_w2,NULL);
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << t1
		        << std::setw(OWID) << t2
		        << std::setw(OWID) << (t1+t2)
		        << std::setw(OWID) << t3
		        << std::setw(OWID) << (t1+t2+t3)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/allParticleVec.size()
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
	                << std::setw(OWID) << getVibraTimeStep()
	                << std::setw(OWID) << getImpactTimeStep()
		        << std::setw(OWID) << timediffsec(time_w1,time_w2)
		        << std::endl;

	    
	    //	    debugInf << std::setw(OWID) << iteration
	    //       << std::setw(OWID) << bdry_penetr[1]
	    //       << std::setw(OWID) << bdry_penetr[2]
	    //       << std::setw(OWID) << bdry_penetr[3]
	    //       << std::setw(OWID) << bdry_penetr[4]
	    //       << std::setw(OWID) << bdry_penetr[6]
	    //       << std::setw(OWID) << bdry_cntnum[1]
	    //       << std::setw(OWID) << bdry_cntnum[2]
	    //       << std::setw(OWID) << bdry_cntnum[3]
	    //       << std::setw(OWID) << bdry_cntnum[4]
	    //       << std::setw(OWID) << bdry_cntnum[6]
	    //       << std::endl;

	}

	// 8. loop break conditions.

    } while (++iteration < totalSteps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, "dep_particle"); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, "dep_contact"); strcat(stepsfp, "_end");
    printContact(stepsfp);
    g_timeinf << std::setw(4) << "end" << " " << ctime(&timeStamp) << std::flush;

    // post_2. close streams
    progressinf.close();
    debugInf.close();
}


void Assembly::depositAfterCavity(int   totalSteps,  
				  int   snapNum,
				  int   interval,
				  const char *iniptclfile,   
				  const char *inibdryfile,
				  const char *inicavefile,
				  const char *ParticleFile, 
				  const char *contactfile,
				  const char *progressfile, 
				  const char *debugfile)
{
    // pre_1: open streams for output.
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf.precision(OPREC);
    progressinf << std::setw(OWID) << "iteration"
	        << std::setw(OWID) << "poss_contact"
	        << std::setw(OWID) << "actual_contact"
	        << std::setw(OWID) << "penetration"
	        << std::setw(OWID) << "avg_normal"
	        << std::setw(OWID) << "avg_tangt"
	        << std::setw(OWID) << "avg_velocity"
	        << std::setw(OWID) << "avg_omga"
	        << std::setw(OWID) << "avg_force"
	        << std::setw(OWID) << "avg_moment"
	        << std::setw(OWID) << "trans_energy"
	        << std::setw(OWID) << "rotat_energy"
	        << std::setw(OWID) << "kinet_energy"
	        << std::setw(OWID) << "poten_energy"
	        << std::setw(OWID) << "total_energy"
	        << std::setw(OWID) << "void_ratio"
	        << std::setw(OWID) << "porosity"
	        << std::setw(OWID) << "coord_number"
	        << std::setw(OWID) << "density"
	        << std::setw(OWID) << "sigma_y1"
	        << std::setw(OWID) << "sigma_y2"
	        << std::setw(OWID) << "sigma_x1"
	        << std::setw(OWID) << "sigma_x2"
	        << std::setw(OWID) << "sigma_z1"
	        << std::setw(OWID) << "sigma_z2"
	        << std::setw(OWID) << "mean_stress"
	        << std::setw(OWID) << "dimx"
	        << std::setw(OWID) << "dimy"
	        << std::setw(OWID) << "dimz"
	        << std::setw(OWID) << "volume"
	        << std::setw(OWID) << "epsilon_x"
	        << std::setw(OWID) << "epsilon_y"
	        << std::setw(OWID) << "epsilon_z"
	        << std::setw(OWID) << "epsilon_v"
	        << std::setw(OWID) << "vibra_t_step"
	        << std::setw(OWID) << "impact_t_step"
	        << std::setw(OWID) << "wall_time" << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1); }
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
        findBdryContact();
	findCavityContact();

	// 2. set particles' forces/moments as zero before each re-calculation,
	clearContactForce();	

	// 3. calculate contact forces/moments and apply them to particles.
	internalForce(avgNormal, avgTangt);

	// 4. calculate boundary forces/moments and apply them to particles.
	boundaryForce(bdry_penetr, bdry_cntnum);
	cavityBoundaryForce();
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();

	// 6. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getZ() - getApt(6).getZ();
	l24=getApt(2).getY()-getApt(4).getY();
	l13=getApt(1).getX()-getApt(3).getX(); bulkVolume=l13*l24*l56;
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
	    g_timeinf << std::setw(4) << stepsnum << " " << ctime(&timeStamp) << std::flush;
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info.
	if (iteration % interval == 0) {
	    gettimeofday(&time_w2,NULL);
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << t1
		        << std::setw(OWID) << t2
		        << std::setw(OWID) << (t1+t2)
		        << std::setw(OWID) << t3
		        << std::setw(OWID) << (t1+t2+t3)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/allParticleVec.size()
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
		        << std::setw(OWID) << "0"
	                << std::setw(OWID) << getVibraTimeStep()
	                << std::setw(OWID) << getImpactTimeStep()
		        << std::setw(OWID) << timediffsec(time_w1,time_w2)
		        << std::endl;

	    
	    //debugInf << std::setw(OWID) << iteration
	    //       << std::setw(OWID) << bdry_penetr[1]
	    //       << std::setw(OWID) << bdry_penetr[2]
	    //       << std::setw(OWID) << bdry_penetr[3]
	    //       << std::setw(OWID) << bdry_penetr[4]
	    //       << std::setw(OWID) << bdry_penetr[6]
	    //       << std::setw(OWID) << bdry_cntnum[1]
	    //       << std::setw(OWID) << bdry_cntnum[2]
	    //       << std::setw(OWID) << bdry_cntnum[3]
	    //       << std::setw(OWID) << bdry_cntnum[4]
	    //       << std::setw(OWID) << bdry_cntnum[6]
	    //       << std::endl;
	}

	// 8. loop break conditions.

    } while (++iteration < totalSteps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    g_timeinf << std::setw(4) << "end" << " " << ctime(&timeStamp) << std::flush;

    // post_2. close streams
    progressinf.close();
    debugInf.close();
}


void Assembly::deGravitation(int   totalSteps,  
			     int   snapNum,
			     int   interval,
			     bool  toRebuild,
			     const char *iniptclfile,   
			     const char *ParticleFile, 
			     const char *contactfile,
			     const char *progressfile, 
			     const char *debugfile)
{
  // pre_1: open streams for output.
  progressinf.open(progressfile); 
  if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1); }
  progressinf.setf(std::ios::scientific, std::ios::floatfield);
  progressinf.precision(OPREC);
  progressinf << std::setw(OWID) << "iteration"
	      << std::setw(OWID) << "poss_contact"
	      << std::setw(OWID) << "actual_contact"
	      << std::setw(OWID) << "penetration"
	      << std::setw(OWID) << "avg_normal"
	      << std::setw(OWID) << "avg_tangt"
	      << std::setw(OWID) << "avg_velocity"
	      << std::setw(OWID) << "avg_omga"
	      << std::setw(OWID) << "avg_force"
	      << std::setw(OWID) << "avg_moment"
	      << std::setw(OWID) << "trans_energy"
	      << std::setw(OWID) << "rotat_energy"
	      << std::setw(OWID) << "kinet_energy"
	     << std::endl;
  
  debugInf.open(debugfile);
  if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1); }
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
      clearContactForce();	
      
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
	g_timeinf << std::setw(4) << "end" << " " << ctime(&timeStamp) << std::flush;
	++stepsnum;
      }
      
      // 5. (2) output stress and strain info.
      if (iteration % interval == 0) {
	progressinf << std::setw(OWID) << iteration
		    << std::setw(OWID) << getPossContactNum()
		    << std::setw(OWID) << getActualContactNum()
		    << std::setw(OWID) << getAvgPenetration()
		    << std::setw(OWID) << avgNormal
		    << std::setw(OWID) << avgTangt
		    << std::setw(OWID) << getAvgVelocity() 
		    << std::setw(OWID) << getAvgOmga()
		    << std::setw(OWID) << getAvgForce()   
		    << std::setw(OWID) << getAvgMoment()
		    << std::setw(OWID) << getTransEnergy()
		    << std::setw(OWID) << getRotatEnergy()
		    << std::setw(OWID) << getKinetEnergy()
		    << std::endl;
      }
      
      // 7. loop break conditions.
      if (contactVec.size() == 0) break;
      
    } while (++iteration < totalSteps);
  
  // post_1. store the final snapshot of particles & contacts.
  strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
  printParticle(stepsfp);
  
  strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
  printContact(stepsfp);
  g_timeinf << std::setw(4) << "end" << " " << ctime(&timeStamp) << std::flush;
  
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
			 const char *iniptclfile,   
			 const char *ParticleFile, 
			 const char *contactfile,
			 const char *progressfile, 
			 const char *debugfile)
{
    // pre_1: open streams for output.
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "deposit..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average       translational    rotational       "
	        << "kinetic        potential         total           void            sample       coordination"
	        << "       sample           sample          sample          sample          sample          sample"
	        << "          sample          sample          sample         sample           sample         "
	        << " sample          sample          sample          sample          sample" << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "         omga            force           moment         energy           energy          "
	        << "energy         energy            energy          ratio          porosity         number       "
	        << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
	        << "epsilon_v" << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1); }
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
	clearContactForce();	

	// 3. calculate contact forces/moments and apply them to particles.
	internalForce(avgNormal, avgTangt);

	// 4. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();

	// 5. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getZ() - (-dimn/2);
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
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << t1
		        << std::setw(OWID) << t2
		        << std::setw(OWID) << (t1+t2)
		        << std::setw(OWID) << t3
		        << std::setw(OWID) << (t1+t2+t3)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
		        << std::endl;
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
		       const char *iniptclfile,   
		       const char *inibdryfile,
		       const char *ParticleFile, 
		       const char *boundaryfile,
		       const char *contactfile,
		       const char *progressfile, 
		       const char *debugfile)
{
    // pre_1: open streams for output.
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "deposit..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average       translational    rotational       "
	        << "kinetic        potential         total           void            sample       coordination"
	        << "       sample           sample          sample          sample          sample          sample"
	        << "          sample          sample          sample         sample           sample         "
	        << " sample          sample          sample          sample          sample" << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "         omga            force           moment         energy           energy          "
	        << "energy         energy            energy          ratio          porosity         number       "
	        << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
	        << "epsilon_v" << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1); }
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
	findBdryContact();

	// 2. set particles' forces/moments as zero before each re-calculation,
	clearContactForce();	

	// 3. calculate contact forces/moments and apply them to particles.
	internalForce(avgNormal, avgTangt);

	// 4. calculate boundary forces/moments and apply them to particles.
	boundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();

	// 6. calculate sample void ratio.
	l56=getTopFreeParticlePosition().getZ() -getApt(6).getZ();
	l24=getApt(2).getY()-getApt(4).getY();
	l13=getApt(1).getX()-getApt(3).getX(); bulkVolume=l13*l24*l56;
	void_ratio=bulkVolume/getParticleVolume()-1;

	// displacement control
	if (iteration > init_steps) {
	    if (flag==1) // loosen, totally remove the wall
		midctl[0].tran=Vec(timeStep*1.0e+0*flag,0,0);
	    else         // squeeze
		midctl[0].tran=Vec(timeStep*5.0e-3*flag,0,0);
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
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << t1
		        << std::setw(OWID) << t2
		        << std::setw(OWID) << (t1+t2)
		        << std::setw(OWID) << t3
		        << std::setw(OWID) << (t1+t2+t3)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
	    debugInf << std::setw(OWID) << iteration
		       << std::setw(OWID) << bdry_penetr[1]
		       << std::setw(OWID) << bdry_penetr[2]
		       << std::setw(OWID) << bdry_penetr[3]
		       << std::setw(OWID) << bdry_penetr[4]
		       << std::setw(OWID) << bdry_penetr[6]
		       << std::setw(OWID) << bdry_cntnum[1]
		       << std::setw(OWID) << bdry_cntnum[2]
		       << std::setw(OWID) << bdry_cntnum[3]
		       << std::setw(OWID) << bdry_cntnum[4]
		       << std::setw(OWID) << bdry_cntnum[6]
		       << std::endl;

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
			 const char *iniptclfile,   
			 const char *inibdryfile,
			 const char *ParticleFile, 
			 const char *boundaryfile,
			 const char *contactfile,  
			 const char *progressfile,
			 const char *balancedfile, 
			 const char *debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "isotropic..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "isotropic..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3: define variables used in iterations
    REAL W0 = getApt(2).getY()-getApt(4).getY();
    REAL L0 = getApt(1).getX()-getApt(3).getX();
    REAL H0 = getApt(5).getZ()-getApt(6).getZ();
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
	findBdryContact();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearContactForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	boundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getZ()-getApt(6).getZ();
	l24=getApt(2).getY()-getApt(4).getY();
	l13=getApt(1).getX()-getApt(3).getX();    bulkVolume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=bulkVolume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=Vec(0,0,-timeStep*boundaryRate);
	else
	    minctl[0].tran=Vec(0,0,timeStep*boundaryRate);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=Vec(0,0,timeStep*boundaryRate);
	else
	    minctl[1].tran=Vec(0,0,-timeStep*boundaryRate);
	
	if (sigma2_1<sigma)
	    midctl[0].tran=Vec(-timeStep*boundaryRate,0,0);
	else
	    midctl[0].tran=Vec(timeStep*boundaryRate,0,0);
	
	if (sigma2_2<sigma)
	    midctl[1].tran=Vec(timeStep*boundaryRate,0,0);
	else
	    midctl[1].tran=Vec(-timeStep*boundaryRate,0,0);
	
	if (sigma1_1<sigma)
	    maxctl[0].tran=Vec(0,-timeStep*boundaryRate,0);
	else
	    maxctl[0].tran=Vec(0,timeStep*boundaryRate,0);
	
	if (sigma1_2<sigma)
	    maxctl[1].tran=Vec(0,timeStep*boundaryRate,0);
	else
	    maxctl[1].tran=Vec(0,-timeStep*boundaryRate,0);
	
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
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
	    debugInf << std::setw(OWID) << iteration
		       << std::setw(OWID) << bdry_penetr[1]
		       << std::setw(OWID) << bdry_penetr[2]
		       << std::setw(OWID) << bdry_penetr[3]
		       << std::setw(OWID) << bdry_penetr[4]
		       << std::setw(OWID) << bdry_penetr[5]
		       << std::setw(OWID) << bdry_penetr[6]
		       << std::setw(OWID) << bdry_cntnum[1]
		       << std::setw(OWID) << bdry_cntnum[2]
		       << std::setw(OWID) << bdry_cntnum[3]
		       << std::setw(OWID) << bdry_cntnum[4]
		       << std::setw(OWID) << bdry_cntnum[5]
		       << std::setw(OWID) << bdry_cntnum[6]
		       << std::endl;
	}

	// 8. loop break condition
	if (   fabs(sigma1_1-sigma)/sigma < boundaryStressTol && fabs(sigma1_2-sigma)/sigma < boundaryStressTol
	    && fabs(sigma2_1-sigma)/sigma < boundaryStressTol && fabs(sigma2_2-sigma)/sigma < boundaryStressTol
	    && fabs(sigma3_1-sigma)/sigma < boundaryStressTol && fabs(sigma3_2-sigma)/sigma < boundaryStressTol ) {
	    balancedinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()    
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()  // just the mean stress p
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()    
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()  // just the mean stress p
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
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
			 const char *iniptclfile,   
			 const char *inibdryfile,
			 const char *ParticleFile, 
			 const char *boundaryfile,
			 const char *contactfile,  
			 const char *progressfile,
			 const char *balancedfile, 
			 const char *debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "isotropic..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "isotropic..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3: define variables used in iterations
    REAL W0 = getApt(2).getY()-getApt(4).getY();
    REAL L0 = getApt(1).getX()-getApt(3).getX();
    REAL H0 = getApt(5).getZ()-getApt(6).getZ();
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
	findBdryContact();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearContactForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	boundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getZ()-getApt(6).getZ();
	l24=getApt(2).getY()-getApt(4).getY();
	l13=getApt(1).getX()-getApt(3).getX();    bulkVolume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=bulkVolume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=Vec(0,0,-timeStep*boundaryRate);
	else
	    minctl[0].tran=Vec(0,0,timeStep*boundaryRate);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=Vec(0,0,timeStep*boundaryRate);
	else
	    minctl[1].tran=Vec(0,0,-timeStep*boundaryRate);
	
	if (sigma2_1<sigma)
	    midctl[0].tran=Vec(-timeStep*boundaryRate,0,0);
	else
	    midctl[0].tran=Vec(timeStep*boundaryRate,0,0);
	
	if (sigma2_2<sigma)
	    midctl[1].tran=Vec(timeStep*boundaryRate,0,0);
	else
	    midctl[1].tran=Vec(-timeStep*boundaryRate,0,0);
	
	if (sigma1_1<sigma)
	    maxctl[0].tran=Vec(0,-timeStep*boundaryRate,0);
	else
	    maxctl[0].tran=Vec(0,timeStep*boundaryRate,0);
	
	if (sigma1_2<sigma)
	    maxctl[1].tran=Vec(0,timeStep*boundaryRate,0);
	else
	    maxctl[1].tran=Vec(0,-timeStep*boundaryRate,0);
	
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
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
	    debugInf << std::setw(OWID) << iteration
		       << std::setw(OWID) << bdry_penetr[1]
		       << std::setw(OWID) << bdry_penetr[2]
		       << std::setw(OWID) << bdry_penetr[3]
		       << std::setw(OWID) << bdry_penetr[4]
		       << std::setw(OWID) << bdry_penetr[5]
		       << std::setw(OWID) << bdry_penetr[6]
		       << std::setw(OWID) << bdry_cntnum[1]
		       << std::setw(OWID) << bdry_cntnum[2]
		       << std::setw(OWID) << bdry_cntnum[3]
		       << std::setw(OWID) << bdry_cntnum[4]
		       << std::setw(OWID) << bdry_cntnum[5]
		       << std::setw(OWID) << bdry_cntnum[6]
		       << std::endl;
	}

	// 8. find the balanced status and increase confining pressure
	if (   fabs(sigma1_1-sigma)/sigma < boundaryStressTol && fabs(sigma1_2-sigma)/sigma < boundaryStressTol
	    && fabs(sigma2_1-sigma)/sigma < boundaryStressTol && fabs(sigma2_2-sigma)/sigma < boundaryStressTol
	    && fabs(sigma3_1-sigma)/sigma < boundaryStressTol && fabs(sigma3_2-sigma)/sigma < boundaryStressTol ) {
	    balancedinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()    
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()  // just the mean stress p
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
	    sigma += sigma_inc;
	}

	// 9. loop break condition
	if (   fabs(sigma1_1-sigma_b)/sigma_b < boundaryStressTol && fabs(sigma1_2-sigma_b)/sigma_b < boundaryStressTol
	    && fabs(sigma2_1-sigma_b)/sigma_b < boundaryStressTol && fabs(sigma2_2-sigma_b)/sigma_b < boundaryStressTol
	    && fabs(sigma3_1-sigma_b)/sigma_b < boundaryStressTol && fabs(sigma3_2-sigma_b)/sigma_b < boundaryStressTol ) {
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()    
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()  // just the mean stress p
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
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
			 const char *iniptclfile,  
			 const char *inibdryfile,
			 const char *ParticleFile, 
			 const char *boundaryfile,
			 const char *contactfile,  
			 const char *progressfile,
			 const char *balancedfile, 
			 const char *debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "isotropic..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "isotropic..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3: define variables used in iterations
    REAL W0 = getApt(2).getY()-getApt(4).getY();
    REAL L0 = getApt(1).getX()-getApt(3).getX();
    REAL H0 = getApt(5).getZ()-getApt(6).getZ();
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
	findBdryContact();
	
	// 2. set particles' forces/moments as zero before each re-calculation
	clearContactForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	boundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getZ()-getApt(6).getZ();
	l24=getApt(2).getY()-getApt(4).getY();
	l13=getApt(1).getX()-getApt(3).getX();    bulkVolume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=bulkVolume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=Vec(0,0,-timeStep*boundaryRate);
	else
	    minctl[0].tran=Vec(0,0,timeStep*boundaryRate);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=Vec(0,0,timeStep*boundaryRate);
	else
	    minctl[1].tran=Vec(0,0,-timeStep*boundaryRate);
	
	if (sigma2_1<sigma)
	    midctl[0].tran=Vec(-timeStep*boundaryRate,0,0);
	else
	    midctl[0].tran=Vec(timeStep*boundaryRate,0,0);
	
	if (sigma2_2<sigma)
	    midctl[1].tran=Vec(timeStep*boundaryRate,0,0);
	else
	    midctl[1].tran=Vec(-timeStep*boundaryRate,0,0);
	
	if (sigma1_1<sigma)
	    maxctl[0].tran=Vec(0,-timeStep*boundaryRate,0);
	else
	    maxctl[0].tran=Vec(0,timeStep*boundaryRate,0);
	
	if (sigma1_2<sigma)
	    maxctl[1].tran=Vec(0,timeStep*boundaryRate,0);
	else
	    maxctl[1].tran=Vec(0,-timeStep*boundaryRate,0);
	
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
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
	    debugInf << std::setw(OWID) << iteration
		       << std::setw(OWID) << bdry_penetr[1]
		       << std::setw(OWID) << bdry_penetr[2]
		       << std::setw(OWID) << bdry_penetr[3]
		       << std::setw(OWID) << bdry_penetr[4]
		       << std::setw(OWID) << bdry_penetr[5]
		       << std::setw(OWID) << bdry_penetr[6]
		       << std::setw(OWID) << bdry_cntnum[1]
		       << std::setw(OWID) << bdry_cntnum[2]
		       << std::setw(OWID) << bdry_cntnum[3]
		       << std::setw(OWID) << bdry_cntnum[4]
		       << std::setw(OWID) << bdry_cntnum[5]
		       << std::setw(OWID) << bdry_cntnum[6]
		       << std::endl;
	}

	// 8. find the balanced status and increase confining pressure
	if (   fabs(sigma1_1-sigma)/sigma < boundaryStressTol && fabs(sigma1_2-sigma)/sigma < boundaryStressTol
	    && fabs(sigma2_1-sigma)/sigma < boundaryStressTol && fabs(sigma2_2-sigma)/sigma < boundaryStressTol
	    && fabs(sigma3_1-sigma)/sigma < boundaryStressTol && fabs(sigma3_2-sigma)/sigma < boundaryStressTol ) {
	    balancedinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()    
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()  // just the mean stress p
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
	    sigma += sigma_inc;
	    if (sigma==sigma_values[i+1]) {
		i++;
		sigma=sigma_values[i];
		sigma_inc=(sigma_values[i+1]-sigma_values[i])/sigma_division;
	    }

	}

	// 9. loop break condition
	if (   fabs(sigma1_1-sigma_b)/sigma_b < boundaryStressTol && fabs(sigma1_2-sigma_b)/sigma_b < boundaryStressTol
	    && fabs(sigma2_1-sigma_b)/sigma_b < boundaryStressTol && fabs(sigma2_2-sigma_b)/sigma_b < boundaryStressTol
	    && fabs(sigma3_1-sigma_b)/sigma_b < boundaryStressTol && fabs(sigma3_2-sigma_b)/sigma_b < boundaryStressTol ) {
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()    
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()  // just the mean stress p
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
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
			const char *iniptclfile,  
			const char *inibdryfile,
			const char *ParticleFile,
			const char *boundaryfile,
			const char *contactfile,  
			const char *progressfile,
			const char *balancedfile, 
			const char *debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "odometer..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "odometer..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries
 
    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).getY()-getApt(4).getY();
    REAL L0 = getApt(1).getX()-getApt(3).getX();
    REAL H0 = getApt(5).getZ()-getApt(6).getZ();
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
	findBdryContact();
	
	// 2. set particles' forces and moments as zero before each re-calculation
	clearContactForce();	

	// 3. calculate contact forces and moments
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces
	boundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/displacement based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getZ()-getApt(6).getZ();
	l24=getApt(2).getY()-getApt(4).getY();
	l13=getApt(1).getX()-getApt(3).getX();    bulkVolume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=bulkVolume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=Vec(0,0,-timeStep*boundaryRate);
	else
	    minctl[0].tran=Vec(0,0,timeStep*boundaryRate);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=Vec(0,0,timeStep*boundaryRate);
	else
	    minctl[1].tran=Vec(0,0,-timeStep*boundaryRate);
	
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
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
	}

	// 8. find balanced status of odometer compression
	if (fabs(sigma3_1-sigma)/sigma < boundaryStressTol && fabs(sigma3_2-sigma)/sigma < boundaryStressTol ) {
	    balancedinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()    
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()  // just the mean stress p
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
	    sigma += sigma_inc;
	}

	// 9. loop break condition
	if (fabs(sigma3_1-sigma_1)/sigma_1 < boundaryStressTol && fabs(sigma3_2-sigma_1)/sigma_1 < boundaryStressTol) {
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()    
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()  // just the mean stress p
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
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
			const char *iniptclfile,  
			const char *inibdryfile,
			const char *ParticleFile, 
			const char *boundaryfile,
			const char *contactfile,  
			const char *progressfile,
			const char *balancedfile, 
			const char *debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "odometer..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "odometer..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries
 
    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).getY()-getApt(4).getY();
    REAL L0 = getApt(1).getX()-getApt(3).getX();
    REAL H0 = getApt(5).getZ()-getApt(6).getZ();
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
	findBdryContact();

	// 2. set particles' forces and moments as zero before each re-calculation
	clearContactForce();	

	// 3. calculate contact forces and moments
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces
	boundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/displacement based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getZ()-getApt(6).getZ();
	l24=getApt(2).getY()-getApt(4).getY();
	l13=getApt(1).getX()-getApt(3).getX();    bulkVolume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=bulkVolume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=Vec(0,0,-timeStep*boundaryRate);
	else
	    minctl[0].tran=Vec(0,0,timeStep*boundaryRate);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=Vec(0,0,timeStep*boundaryRate);
	else
	    minctl[1].tran=Vec(0,0,-timeStep*boundaryRate);
	
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
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
	}

	// 8. find balanced status of odometer compression
	if (fabs(sigma3_1-sigma)/sigma < boundaryStressTol && fabs(sigma3_2-sigma)/sigma < boundaryStressTol ) {
	    balancedinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()    
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()  // just the mean stress p
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
	    sigma += sigma_inc;
	    if (sigma==sigma_values[i+1]) {
		i++;
		sigma=sigma_values[i];
		sigma_inc=(sigma_values[i+1]-sigma_values[i])/sigma_division;
	    }
	}

	// 9. loop break condition
	if (fabs(sigma3_1-sigma_b)/sigma_b < boundaryStressTol && fabs(sigma3_2-sigma_b)/sigma_b < boundaryStressTol) {
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()    
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()  // just the mean stress p
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
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
			   const char *iniptclfile, 
			   const char *ParticleFile,
			   const char *contactfile, 
			   const char *progressfile,
			   const char *debugfile) 
{
  // pre_1: open streams for output
  // ParticleFile and contactfile are used for snapNum at the end.
  progressinf.open(progressfile);
  if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1);}
  progressinf.setf(std::ios::scientific, std::ios::floatfield);
  progressinf.precision(OPREC);
  progressinf << std::setw(OWID) << "iteration"
	      << std::setw(OWID) << "poss_contact"
	      << std::setw(OWID) << "actual_contact"
	      << std::setw(OWID) << "penetration"
	      << std::setw(OWID) << "avg_normal"
	      << std::setw(OWID) << "avg_tangt"
	      << std::setw(OWID) << "avg_velocity"
	      << std::setw(OWID) << "avg_omga"
	      << std::setw(OWID) << "avg_force"
	      << std::setw(OWID) << "avg_moment"
	      << std::setw(OWID) << "trans_energy"
	      << std::setw(OWID) << "rotat_energy"
	      << std::setw(OWID) << "kinet_energy"
	      << std::endl;
  
  debugInf.open(debugfile);
  if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
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
  REAL x1 = v1.getX();
  REAL y1 = v1.getY();
  REAL z1 = v1.getZ();
  REAL x2 = v2.getX();
  REAL y2 = v2.getY();
  REAL z2 = v2.getZ();
  std::vector<Particle*>::const_iterator  it;
  Vec pos;
  for (it=particleVec.begin();it!=particleVec.end();++it)
    {
      pos = (*it)->getCurrPos();
      if (pos.getX() < x1)
	(*it)->setConstForce( Vec(mag, 0, 0) );
      else if (pos.getX() > x2)
	(*it)->setConstForce( Vec(-mag, 0, 0) );
      else if (pos.getY() < y1)
	(*it)->setConstForce( Vec(0, mag, 0) );
      else if (pos.getY() > y2)
	(*it)->setConstForce( Vec(0, -mag, 0) );
      else if (pos.getZ() < z1)
	(*it)->setConstForce( Vec(0, 0, mag) );
      else if (pos.getZ() > z2)
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
      clearContactForce(); // const_force/moment NOT cleared by this call	
      
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
	g_timeinf << std::setw(4) << stepsnum << " " << ctime(&timeStamp) << std::flush;
	++stepsnum;
      }
      
      // 6. (2) output stress and strain info
	if (iteration % interval == 0 ){
	  progressinf << std::setw(OWID) << iteration
		      << std::setw(OWID) << getPossContactNum()
		      << std::setw(OWID) << getActualContactNum()
		      << std::setw(OWID) << getAvgPenetration()
		      << std::setw(OWID) << avgNormal
		      << std::setw(OWID) << avgTangt
		      << std::setw(OWID) << getAvgVelocity() 
		      << std::setw(OWID) << getAvgOmga()
		      << std::setw(OWID) << getAvgForce()   
		      << std::setw(OWID) << getAvgMoment()
		      << std::setw(OWID) << getTransEnergy()
		      << std::setw(OWID) << getRotatEnergy()
		      << std::setw(OWID) << getKinetEnergy()
		      << std::endl;
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
  g_timeinf << std::setw(4) << "end" << " " << ctime(&timeStamp) << std::flush;
  
  // post_2. close streams
  progressinf.close();
  debugInf.close();
}


// This function initializes triaxial sample to a certain confining pressure.
void Assembly::triaxialPtclBdryIni(int   totalSteps,  
				   int   snapNum, 
				   int   interval,
				   REAL  sigma,
				   const char *iniptclfile, 
				   const char *inibdryfile,
				   const char *ParticleFile,
				   const char *boundaryfile,
				   const char *contactfile, 
				   const char *progressfile,
				   const char *debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "triaxial..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL H0 = getApt(5).getZ()-getApt(6).getZ();
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
	findBdryContact();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearContactForce();	

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
	    minctl[0].tran=Vec(0,0,-timeStep*boundaryRate);
	else
	    minctl[0].tran=Vec(0,0, timeStep*boundaryRate);

	if (sigma3_2 < sigma)
	    minctl[1].tran=Vec(0,0, timeStep*boundaryRate);
	else
	    minctl[1].tran=Vec(0,0,-timeStep*boundaryRate);

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
	l56=getApt(5).getZ()-getApt(6).getZ();
	epsilon_h = (H0-l56)/H0;
	if (iteration % interval == 0 ){
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << 0 << std::setw(OWID) << 0
		        << std::setw(OWID) << 0 << std::setw(OWID) << 0
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()
		        << std::setw(OWID) << 0 << std::setw(OWID) << 0 << std::setw(OWID) << l56
		        << std::setw(OWID) << 0
		        << std::setw(OWID) << 0
		        << std::setw(OWID) << 0
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << 0
		        << std::setw(OWID) << 0
		        << std::setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
		        << std::endl;

	}

	// 9. loop break condition: through displacement control mechanism
	if (   fabs(sigma3_1-sigma)/sigma < boundaryStressTol && fabs(sigma3_2-sigma)/sigma < boundaryStressTol )
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
				const char *iniptclfile, 
				const char *inibdryfile,
				const char *ParticleFile,
				const char *boundaryfile,
				const char *contactfile, 
				const char *progressfile,
				const char *balancedfile,
				const char *debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "triaxial..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "triaxial..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL H0 = getApt(5).getZ()-getApt(6).getZ();
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
	findBdryContact();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearContactForce();	

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
	minctl[0].tran=Vec(0,0,-timeStep*boundaryRate);
	minctl[1].tran=Vec(0,0, timeStep*boundaryRate);

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
	l56=getApt(5).getZ()-getApt(6).getZ();
	epsilon_h = (H0-l56)/H0;
	if (iteration % interval == 0 ){
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << 0 << std::setw(OWID) << 0
		        << std::setw(OWID) << 0 << std::setw(OWID) << 0
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()
		        << std::setw(OWID) << 0 << std::setw(OWID) << 0 << std::setw(OWID) << l56
		        << std::setw(OWID) << 0
		        << std::setw(OWID) << 0
		        << std::setw(OWID) << 0
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << 0
		        << std::setw(OWID) << 0
		        << std::setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
		        << std::endl;

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
			const char *iniptclfile, 
			const char *inibdryfile,
			const char *ParticleFile,
			const char *boundaryfile,
			const char *contactfile, 
			const char *progressfile,
			const char *balancedfile,
			const char *debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf.precision(OPREC);
    progressinf << std::setw(OWID) << "iteration"
	        << std::setw(OWID) << "possible"
	        << std::setw(OWID) << "actual"
	        << std::setw(OWID) << "average"
	        << std::setw(OWID) << "average"
	        << std::setw(OWID) << "average"
	        << std::setw(OWID) << "average"
	        << std::setw(OWID) << "average"
	        << std::setw(OWID) << "average"
	        << std::setw(OWID) << "average"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "void"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "coordinate"
	        << std::setw(OWID) << "vibra"
	        << std::setw(OWID) << "impact"
	        << std::setw(OWID) << "wall-clock" << std::endl
	        << std::setw(OWID) << "number"
	        << std::setw(OWID) << "contacts"
	        << std::setw(OWID) << "contacts"
	        << std::setw(OWID) << "penetration"
	        << std::setw(OWID) << "contact_normal"
	        << std::setw(OWID) << "contact_tangt"
	        << std::setw(OWID) << "velocity"
	        << std::setw(OWID) << "omga"
	        << std::setw(OWID) << "force"
	        << std::setw(OWID) << "moment"
	        << std::setw(OWID) << "density"
	        << std::setw(OWID) << "sigma1_1"
	        << std::setw(OWID) << "sigma1_2"
	        << std::setw(OWID) << "sigma2_1"
	        << std::setw(OWID) << "sigma2_2"
	        << std::setw(OWID) << "sigma3_1"
	        << std::setw(OWID) << "sigma3_2"
	        << std::setw(OWID) << "mean_stress"
	        << std::setw(OWID) << "width"
	        << std::setw(OWID) << "length"
	        << std::setw(OWID) << "height"
	        << std::setw(OWID) << "volume"
	        << std::setw(OWID) << "epsilon_w"
	        << std::setw(OWID) << "epsilon_l"
	        << std::setw(OWID) << "epsilon_h"
	        << std::setw(OWID) << "epsilon_v"
	        << std::setw(OWID) << "ratio"
	        << std::setw(OWID) << "porosity"
	        << std::setw(OWID) << "number"
	        << std::setw(OWID) << "t_step"
	        << std::setw(OWID) << "t_step"
	        << std::setw(OWID) << "time" << std::endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf.precision(OPREC);
    balancedinf << std::setw(OWID) << "iteration"
	        << std::setw(OWID) << "possible"
	        << std::setw(OWID) << "actual"
	        << std::setw(OWID) << "average"
	        << std::setw(OWID) << "average"
	        << std::setw(OWID) << "average"
	        << std::setw(OWID) << "average"
	        << std::setw(OWID) << "average"
	        << std::setw(OWID) << "average"
	        << std::setw(OWID) << "average"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "void"
	        << std::setw(OWID) << "sample"
	        << std::setw(OWID) << "coordinate"
	        << std::setw(OWID) << "vibra"
	        << std::setw(OWID) << "impact"
	        << std::setw(OWID) << "wall-clock" << std::endl
	        << std::setw(OWID) << "number"
	        << std::setw(OWID) << "contacts"
	        << std::setw(OWID) << "contacts"
	        << std::setw(OWID) << "penetration"
	        << std::setw(OWID) << "contact_normal"
	        << std::setw(OWID) << "contact_tangt"
	        << std::setw(OWID) << "velocity"
	        << std::setw(OWID) << "omga"
	        << std::setw(OWID) << "force"
	        << std::setw(OWID) << "moment"
	        << std::setw(OWID) << "density"
	        << std::setw(OWID) << "sigma1_1"
	        << std::setw(OWID) << "sigma1_2"
	        << std::setw(OWID) << "sigma2_1"
	        << std::setw(OWID) << "sigma2_2"
	        << std::setw(OWID) << "sigma3_1"
	        << std::setw(OWID) << "sigma3_2"
	        << std::setw(OWID) << "mean_stress"
	        << std::setw(OWID) << "width"
	        << std::setw(OWID) << "length"
	        << std::setw(OWID) << "height"
	        << std::setw(OWID) << "volume"
	        << std::setw(OWID) << "epsilon_w"
	        << std::setw(OWID) << "epsilon_l"
	        << std::setw(OWID) << "epsilon_h"
	        << std::setw(OWID) << "epsilon_v"
	        << std::setw(OWID) << "ratio"
	        << std::setw(OWID) << "porosity"
	        << std::setw(OWID) << "number"
	        << std::setw(OWID) << "t_step"
	        << std::setw(OWID) << "t_step"
	        << std::setw(OWID) << "time" << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).getY()-getApt(4).getY();
    REAL L0 = getApt(1).getX()-getApt(3).getX();
    REAL H0 = getApt(5).getZ()-getApt(6).getZ();
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
	findBdryContact();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearContactForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	boundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getZ()-getApt(6).getZ();
	l24=getApt(2).getY()-getApt(4).getY();
	l13=getApt(1).getX()-getApt(3).getX();    bulkVolume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=bulkVolume/getParticleVolume()-1;

	// displacement control
	minctl[0].tran=Vec(0,0,-timeStep*boundaryRate);
	minctl[1].tran=Vec(0,0, timeStep*boundaryRate);

	// force control
	if (sigma2_1<sigma_a)
	    midctl[0].tran=Vec(-timeStep*boundaryRate,0,0);
	else
	    midctl[0].tran=Vec(timeStep*boundaryRate,0,0);
	
	if (sigma2_2<sigma_a)
	    midctl[1].tran=Vec(timeStep*boundaryRate,0,0);
	else
	    midctl[1].tran=Vec(-timeStep*boundaryRate,0,0);
	
	if (sigma1_1<sigma_a)
	    maxctl[0].tran=Vec(0,-timeStep*boundaryRate,0);
	else
	    maxctl[0].tran=Vec(0,timeStep*boundaryRate,0);
	
	if (sigma1_2<sigma_a)
	    maxctl[1].tran=Vec(0,timeStep*boundaryRate,0);
	else
	    maxctl[1].tran=Vec(0,-timeStep*boundaryRate,0);
	
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
	    g_timeinf << std::setw(4) << stepsnum << " " << ctime(&timeStamp) << std::flush;
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (iteration % interval == 0 ){
	    gettimeofday(&time_w2,NULL);
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
	                << std::setw(OWID) << getVibraTimeStep()
	                << std::setw(OWID) << getImpactTimeStep()
		        << std::setw(OWID) << timediffsec(time_w1,time_w2)
		        << std::endl;
	    debugInf << std::setw(OWID) << iteration
		       << std::setw(OWID) << getTransEnergy()
		       << std::setw(OWID) << getRotatEnergy()
		       << std::setw(OWID) << bdry_penetr[1]
		       << std::setw(OWID) << bdry_penetr[2]
		       << std::setw(OWID) << bdry_penetr[3]
		       << std::setw(OWID) << bdry_penetr[4]
		       << std::setw(OWID) << bdry_penetr[5]
		       << std::setw(OWID) << bdry_penetr[6]
		       << std::setw(OWID) << bdry_cntnum[1]
		       << std::setw(OWID) << bdry_cntnum[2]
		       << std::setw(OWID) << bdry_cntnum[3]
		       << std::setw(OWID) << bdry_cntnum[4]
		       << std::setw(OWID) << bdry_cntnum[5]
		       << std::setw(OWID) << bdry_cntnum[6]
		       << std::endl;
	}

	// Most time it is balanced, so use progressinf instead.
	// 8. find the balanced status and increase confining pressure
	if (   fabs(sigma1_1-sigma_a)/sigma_a < boundaryStressTol && fabs(sigma1_2-sigma_a)/sigma_a < boundaryStressTol
	    && fabs(sigma2_1-sigma_a)/sigma_a < boundaryStressTol && fabs(sigma2_2-sigma_a)/sigma_a < boundaryStressTol
	    && fabs(sigma3_1-sigma3_2)/(sigma3_1+sigma3_2)*2<=0.05) {
	    balancedinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()    
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()  // just the mean stress p
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
	                << std::setw(OWID) << getVibraTimeStep()
	                << std::setw(OWID) << getImpactTimeStep()
		        << std::setw(OWID) << timediffsec(time_w1,time_w2)
		        << std::endl;
	}

	// 9. loop break condition: through displacement control mechanism
	
    } while (++iteration < totalSteps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    g_timeinf << std::setw(4) << "end" << " " << ctime(&timeStamp) << std::flush;

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
			const char *iniptclfile,  
			const char *inibdryfile,
			const char *ParticleFile,
			const char *boundaryfile,
			const char *contactfile,
			const char *progressfile,
			const char *balancedfile,
			const char *debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "triaxial..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "triaxial..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).getY()-getApt(4).getY();
    REAL L0 = getApt(1).getX()-getApt(3).getX();
    REAL H0 = getApt(5).getZ()-getApt(6).getZ();
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
	findBdryContact();
	
	// 2. set particles' forces/moments as zero before each re-calculation
	clearContactForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	boundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getZ()-getApt(6).getZ();
	l24=getApt(2).getY()-getApt(4).getY();
	l13=getApt(1).getX()-getApt(3).getX();    bulkVolume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=bulkVolume/getParticleVolume()-1;

	// displacement control
	if (iteration <= unload_step){ //loading
	    minctl[0].tran=Vec(0,0,-timeStep*boundaryRate);
	    minctl[1].tran=Vec(0,0, timeStep*boundaryRate);
	}
	else { 
	    if (reload==false) { // unloading
		if (fabs(sigma3_1-sigma_a)/sigma_a > boundaryStressTol && 
		    fabs(sigma3_2-sigma_a)/sigma_a > boundaryStressTol){
		    minctl[0].tran=Vec(0,0, timeStep*boundaryRate);
		    minctl[1].tran=Vec(0,0,-timeStep*boundaryRate);
		}
		else  // reloading
		    reload=true;
	    }
	    else {
		minctl[0].tran=Vec(0,0,-timeStep*boundaryRate);
		minctl[1].tran=Vec(0,0, timeStep*boundaryRate);
	    }
	}
	
	// force control
	if (sigma2_1<sigma_a)
	    midctl[0].tran=Vec(-timeStep*boundaryRate,0,0);
	else
	    midctl[0].tran=Vec(timeStep*boundaryRate,0,0);
	
	if (sigma2_2<sigma_a)
	    midctl[1].tran=Vec(timeStep*boundaryRate,0,0);
	else
	    midctl[1].tran=Vec(-timeStep*boundaryRate,0,0);
	
	if (sigma1_1<sigma_a)
	    maxctl[0].tran=Vec(0,-timeStep*boundaryRate,0);
	else
	    maxctl[0].tran=Vec(0,timeStep*boundaryRate,0);
	
	if (sigma1_2<sigma_a)
	    maxctl[1].tran=Vec(0,timeStep*boundaryRate,0);
	else
	    maxctl[1].tran=Vec(0,-timeStep*boundaryRate,0);
	
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
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
	    debugInf << std::setw(OWID) << iteration
		       << std::setw(OWID) << bdry_penetr[1]
		       << std::setw(OWID) << bdry_penetr[2]
		       << std::setw(OWID) << bdry_penetr[3]
		       << std::setw(OWID) << bdry_penetr[4]
		       << std::setw(OWID) << bdry_penetr[5]
		       << std::setw(OWID) << bdry_penetr[6]
		       << std::setw(OWID) << bdry_cntnum[1]
		       << std::setw(OWID) << bdry_cntnum[2]
		       << std::setw(OWID) << bdry_cntnum[3]
		       << std::setw(OWID) << bdry_cntnum[4]
		       << std::setw(OWID) << bdry_cntnum[5]
		       << std::setw(OWID) << bdry_cntnum[6]
		       << std::endl;
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
			     const char *iniptclfile,  
			     const char *inibdryfile,
			     const char *ParticleFile, 
			     const char *boundaryfile,
			     const char *contactfile,  
			     const char *progressfile,
			     const char *debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "pile penetrate..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average       translational    rotational       "
	        << "kinetic        potential        total           sample           sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample"
	        << "          sample          sample          sample" << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "         omga            force           moment         energy           energy          "
	        << "energy         energy          energy          density         "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
	        << "epsilon_v" << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);
    debugInf << " iteration    end_bearing     side_friction   total_force" << std::endl;

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
	findBdryContact();
	
	// 2. set particles' forces/moments as zero before each re-calculation
	clearContactForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	boundaryForce();

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation

	// displacement control of the pile
	pilectl[0].tran=Vec(0,0,-timeStep*pileRate);
	pilectl[1].tran=Vec(0,0,-timeStep*pileRate);

	updateBoundary(pile, pilectl, 2); 
	updateRectPile();
	if (iteration % interval == 0) {
	    REAL  f7=getShearForce( 7).getZ();
	    REAL  f8=getShearForce( 8).getZ();
	    REAL  f9=getShearForce( 9).getZ();
	    REAL f10=getShearForce(10).getZ();
	    REAL  fn=getNormalForce(12).getZ();
	    debugInf << std::setw(OWID) << iteration
		       << std::setw(OWID) << fn
		       << std::setw(OWID) << (f7+f8+f9+f10)
		       << std::setw(OWID) << (fn+f7+f8+f9+f10)
		       << std::endl;
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
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << t1
		        << std::setw(OWID) << t2
		        << std::setw(OWID) << (t1+t2)
		        << std::setw(OWID) << t3
		        << std::setw(OWID) << (t1+t2+t3) << std::endl;
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
			      const char *iniptclfile,
			      const char *ParticleFile, 
			      const char *contactfile,  
			      const char *progressfile,
			      const char *debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "pile penetrate..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average       translational    rotational       "
	        << "kinetic        potential         total           void            sample       coordination"
	        << "       sample           sample          sample          sample          sample          sample"
	        << "          sample          sample          sample         sample           sample         "
	        << " sample          sample          sample          sample          sample" << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "         omga            force           moment         energy           energy          "
	        << "energy         energy            energy          ratio          porosity         number       "
	        << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
	        << "epsilon_v" << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
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
	clearContactForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 5. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getZ() - (-dimn/2);
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
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << t1
		        << std::setw(OWID) << t2
		        << std::setw(OWID) << (t1+t2)
		        << std::setw(OWID) << t3
		        << std::setw(OWID) << (t1+t2+t3)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
		        << std::endl;
	    debugInf << std::setw(OWID) << iteration
		       << std::setw(OWID) << getTopFreeParticlePosition().getZ()
		       << std::setw(OWID) << ellipPileTipZ()
		       << std::setw(OWID) << getTopFreeParticlePosition().getZ()-ellipPileTipZ()
		       << std::setw(OWID) << l13*l24*l56
		       << std::setw(OWID) << ellipPilePeneVol()
		       << std::setw(OWID) << bulkVolume
		       << std::endl;
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
				const char *iniptclfile,
				const char *inibdryfile,
				const char *ParticleFile, 
				const char *contactfile,  
				const char *progressfile,
				const char *debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "penetrator impact..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average       translational    rotational       "
	        << "kinetic        potential         total           void            sample       coordination"
	        << "       sample           sample          sample          sample          sample          sample"
	        << "          sample          sample          sample         sample           sample         "
	        << " sample          sample          sample          sample          sample" << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "         omga            force           moment         energy           energy          "
	        << "energy         energy            energy          ratio          porosity         number       "
	        << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
	        << "epsilon_v" << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
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
	clearContactForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);

	// 4. calculate boundary forces/moments and apply them to particles.
	boundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getZ() - (-dimn/2);
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
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << t1
		        << std::setw(OWID) << t2
		        << std::setw(OWID) << (t1+t2)
		        << std::setw(OWID) << t3
		        << std::setw(OWID) << (t1+t2+t3)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
	    debugInf << std::setw(OWID) << iteration
		       << std::setw(OWID) << bdry_penetr[1]
		       << std::setw(OWID) << bdry_penetr[2]
		       << std::setw(OWID) << bdry_penetr[3]
		       << std::setw(OWID) << bdry_penetr[4]
		       << std::setw(OWID) << bdry_penetr[6]
		       << std::setw(OWID) << bdry_cntnum[1]
		       << std::setw(OWID) << bdry_cntnum[2]
		       << std::setw(OWID) << bdry_cntnum[3]
		       << std::setw(OWID) << bdry_cntnum[4]
		       << std::setw(OWID) << bdry_cntnum[6]
		       << std::endl;

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
				  const char *iniptclfile,
				  const char *ParticleFile, 
				  const char *contactfile,  
				  const char *progressfile,
				  const char *debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "penetrator impact..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average       translational    rotational       "
	        << "kinetic        potential         total           void            sample       coordination"
	        << "       sample           sample          sample          sample          sample          sample"
	        << "          sample          sample          sample         sample           sample         "
	        << " sample          sample          sample          sample          sample" << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "         omga            force           moment         energy           energy          "
	        << "energy         energy            energy          ratio          porosity         number       "
	        << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
	        << "epsilon_v" << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
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
	clearContactForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 5. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getZ() - (-dimn/2);
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
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << t1
		        << std::setw(OWID) << t2
		        << std::setw(OWID) << (t1+t2)
		        << std::setw(OWID) << t3
		        << std::setw(OWID) << (t1+t2+t3)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
		        << std::endl;
	    debugInf << std::setw(OWID) << iteration
		       << std::setw(OWID) << getTopFreeParticlePosition().getZ()
		       << std::setw(OWID) << ellipPileTipZ()
		       << std::setw(OWID) << getTopFreeParticlePosition().getZ()-ellipPileTipZ()
		       << std::setw(OWID) << l13*l24*l56
		       << std::setw(OWID) << ellipPilePeneVol()
		       << std::setw(OWID) << bulkVolume
		       << std::endl;
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
			       const char *iniptclfile,
			       const char *ParticleFile, 
			       const char *contactfile,  
			       const char *progressfile,
			       const char *balancedfile,
			       const char *debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "pile penetrate..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average       translational    rotational       "
	        << "kinetic        potential         total           void            sample       coordination"
	        << "       sample           sample          sample          sample          sample          sample"
	        << "          sample          sample          sample         sample           sample         "
	        << " sample          sample          sample          sample          sample" << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "         omga            force           moment         energy           energy          "
	        << "energy         energy            energy          ratio          porosity         number       "
	        << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
	        << "epsilon_v" << std::endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "pile penetrate..." << std::endl
	        << "   iteration   apply_force    pile_tip_pos     pile_force" << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
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
	clearContactForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();

	// 5. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getZ() - (-dimn/2);
	l24=dimn;
	l13=dimn;
	bulkVolume=l13*l24*l56-ellipPilePeneVol();
	void_ratio=bulkVolume/getParticleVolume()-1;
	
	// 6. update pile external force and position
	if(zforce>ellipPileForce())
	    ellipPileUpdate();

	if(fabs(ellipPileForce()-zforce)/zforce < boundaryStressTol ){
	    balancedinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << zforce
		        << std::setw(OWID) << getTopFreeParticlePosition().getZ()-ellipPileTipZ()
		        << std::setw(OWID) << ellipPileForce()
		        << std::endl;
	    zforce += zforce_inc;
	}

	if( iteration % interval == 0){
	    debugInf << std::setw(OWID) << iteration
		       << std::setw(OWID) << zforce
		       << std::setw(OWID) << getTopFreeParticlePosition().getZ()-ellipPileTipZ()
		       << std::setw(OWID) << ellipPileForce()
		       << std::endl;
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
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << t1
		        << std::setw(OWID) << t2
		        << std::setw(OWID) << (t1+t2)
		        << std::setw(OWID) << t3
		        << std::setw(OWID) << (t1+t2+t3)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
		        << std::endl;
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
			    const char *iniptclfile,  
			    const char *inibdryfile,
			    const char *ParticleFile, 
			    const char *boundaryfile,
			    const char *contactfile,  
			    const char *progressfile,
			    const char *balancedfile, 
			    const char *debugfile) 
{
    // pre_1: open streams for output
    // ParticleFile and contactfile are used for snapNum at the end.
    progressinf.open(progressfile);
    if(!progressinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "true triaxial..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { std::cout << "stream error!" << std::endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf << "true triaxial..." << std::endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << std::endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << std::endl;

    debugInf.open(debugfile);
    if(!debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
    debugInf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).getY()-getApt(4).getY();
    REAL L0 = getApt(1).getX()-getApt(3).getX();
    REAL H0 = getApt(5).getZ()-getApt(6).getZ();
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
	findBdryContact();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearContactForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	boundaryForce();
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getZ()-getApt(6).getZ();
	l24=getApt(2).getY()-getApt(4).getY();
	l13=getApt(1).getX()-getApt(3).getX();    bulkVolume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=bulkVolume/getParticleVolume()-1;

	if (sigma3_1<sigma_h1)
	    minctl[0].tran=Vec(0,0,-timeStep*boundaryRate);
	else
	    minctl[0].tran=Vec(0,0,timeStep*boundaryRate);
	
	if (sigma3_2<sigma_h1)
	    minctl[1].tran=Vec(0,0,timeStep*boundaryRate);
	else
	    minctl[1].tran=Vec(0,0,-timeStep*boundaryRate);
	
	if (sigma2_1<sigma_l1)
	    midctl[0].tran=Vec(-timeStep*boundaryRate,0,0);
	else
	    midctl[0].tran=Vec(timeStep*boundaryRate,0,0);
	
	if (sigma2_2<sigma_l1)
	    midctl[1].tran=Vec(timeStep*boundaryRate,0,0);
	else
	    midctl[1].tran=Vec(-timeStep*boundaryRate,0,0);
	
	if (sigma1_1<sigma_w1)
	    maxctl[0].tran=Vec(0,-timeStep*boundaryRate,0);
	else
	    maxctl[0].tran=Vec(0,timeStep*boundaryRate,0);
	
	if (sigma1_2<sigma_w1)
	    maxctl[1].tran=Vec(0,timeStep*boundaryRate,0);
	else
	    maxctl[1].tran=Vec(0,-timeStep*boundaryRate,0);
	
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
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()   
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
	    debugInf << std::setw(OWID) << iteration
		       << std::setw(OWID) << bdry_penetr[1]
		       << std::setw(OWID) << bdry_penetr[2]
		       << std::setw(OWID) << bdry_penetr[3]
		       << std::setw(OWID) << bdry_penetr[4]
		       << std::setw(OWID) << bdry_penetr[5]
		       << std::setw(OWID) << bdry_penetr[6]
		       << std::setw(OWID) << bdry_cntnum[1]
		       << std::setw(OWID) << bdry_cntnum[2]
		       << std::setw(OWID) << bdry_cntnum[3]
		       << std::setw(OWID) << bdry_cntnum[4]
		       << std::setw(OWID) << bdry_cntnum[5]
		       << std::setw(OWID) << bdry_cntnum[6]
		       << std::endl;
	}

	// 8. find the balanced status and increase confining pressure
	if (   fabs(sigma1_1-sigma_w1)/sigma_w1 < boundaryStressTol && fabs(sigma1_2-sigma_w1)/sigma_w1 < boundaryStressTol
	    && fabs(sigma2_1-sigma_l1)/sigma_l1 < boundaryStressTol && fabs(sigma2_2-sigma_l1)/sigma_l1 < boundaryStressTol
	    && fabs(sigma3_1-sigma_h1)/sigma_h1 < boundaryStressTol && fabs(sigma3_2-sigma_h1)/sigma_h1 < boundaryStressTol ) {
	    balancedinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()    
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()  // just the mean stress p
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
	    sigma_w1 += sigma_w_inc;
	    sigma_l1 += sigma_l_inc;
	    sigma_h1 += sigma_h_inc;
	}

	// 9. loop break condition
	if (   fabs(sigma1_1-sigma_w)/sigma_w < boundaryStressTol && fabs(sigma1_2-sigma_w)/sigma_w < boundaryStressTol
	    && fabs(sigma2_1-sigma_l)/sigma_l < boundaryStressTol && fabs(sigma2_2-sigma_l)/sigma_l < boundaryStressTol
	    && fabs(sigma3_1-sigma_h)/sigma_h < boundaryStressTol && fabs(sigma3_2-sigma_h)/sigma_h < boundaryStressTol ) {
	    progressinf << std::setw(OWID) << iteration
		        << std::setw(OWID) << getPossContactNum()
		        << std::setw(OWID) << getActualContactNum()
		        << std::setw(OWID) << getAvgPenetration()
		        << std::setw(OWID) << avgNormal
		        << std::setw(OWID) << avgTangt
		        << std::setw(OWID) << getAvgVelocity() 
		        << std::setw(OWID) << getAvgOmga()
		        << std::setw(OWID) << getAvgForce()    
		        << std::setw(OWID) << getAvgMoment()
		        << std::setw(OWID) << getDensity()
		        << std::setw(OWID) << sigma1_1 << std::setw(OWID) << sigma1_2
		        << std::setw(OWID) << sigma2_1 << std::setw(OWID) << sigma2_2
		        << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
		        << std::setw(OWID) << getAvgPressure()  // just the mean stress p
		        << std::setw(OWID) << l24 << std::setw(OWID) << l13 << std::setw(OWID) << l56
		        << std::setw(OWID) << bulkVolume
		        << std::setw(OWID) << epsilon_w
		        << std::setw(OWID) << epsilon_l
		        << std::setw(OWID) << epsilon_h
		        << std::setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << std::setw(OWID) << void_ratio
		        << std::setw(OWID) << void_ratio/(1+void_ratio)
		        << std::setw(OWID) << 2.0*(getActualContactNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/allParticleVec.size()
		        << std::endl;
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
	      const char *boundaryFile)
{
  std::ofstream ofs(boundaryFile);
  if(!ofs) { std::cout << "stream error!" << std::endl; exit(-1);}
  ofs.setf(std::ios::scientific, std::ios::floatfield);

  Vec  v1 = allContainer.getMinCorner();
  Vec  v2 = allContainer.getMaxCorner();
  Vec  v0 = allContainer.getCenter();
  REAL x1 = v1.getX();
  REAL y1 = v1.getY();
  REAL z1 = v1.getZ();
  REAL x2 = v2.getX();
  REAL y2 = v2.getY();
  REAL z2 = v2.getZ();
  REAL x0 = v0.getX();
  REAL y0 = v0.getY();
  REAL z0 = v0.getZ();

  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z1
      << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z2 << std::endl
      << std::setw(OWID) << boundaryNum << std::endl << std::endl;
  
  if (boundaryNum == 1){   // only a bottom boundary
    ofs << std::setw(OWID) << 1 << std::endl
        << std::setw(OWID) << 6
        << std::setw(OWID) << 1 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0
        << std::setw(OWID) << -1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0
        << std::setw(OWID) << z1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl;
    
  }
  else if (boundaryNum == 5){ // no top boundary
    // boundary 1
    ofs << std::setw(OWID) << 1 << std::endl
        << std::setw(OWID) << 1
        << std::setw(OWID) << 4 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 
        << std::setw(OWID) << x2
        << std::setw(OWID) << y0
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << -1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << x0
        << std::setw(OWID) << y1
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << x0
        << std::setw(OWID) << y2
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0
        << std::setw(OWID) << -1
        << std::setw(OWID) << x0     
        << std::setw(OWID) << y0
        << std::setw(OWID) << z1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl << std::endl
      
      // boundary 2
        << std::setw(OWID) << 1 << std::endl
        << std::setw(OWID) << 2
        << std::setw(OWID) << 4 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << x0     
        << std::setw(OWID) << y2
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << 0 
        << std::setw(OWID) << x2
        << std::setw(OWID) << y0      
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << -1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 
        << std::setw(OWID) << x1
        << std::setw(OWID) << y0     
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << 0
        << std::setw(OWID) << -1
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y0     
        << std::setw(OWID) << z1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl << std::endl
      
      // boundary 3
        << std::setw(OWID) << 1 << std::endl
        << std::setw(OWID) << 3
        << std::setw(OWID) << 4 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << -1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 
        << std::setw(OWID) << x1
        << std::setw(OWID) << y0
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << -1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << x0     
        << std::setw(OWID) << y1
        << std::setw(OWID) << y0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y2
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << 0
        << std::setw(OWID) << -1
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y0     
        << std::setw(OWID) << z1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl << std::endl
      
      // boundary 4
        << std::setw(OWID) << 1 << std::endl
        << std::setw(OWID) << 4
        << std::setw(OWID) << 4 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << -1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y1
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << 0 
        << std::setw(OWID) << x2
        << std::setw(OWID) << y0
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << -1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 
        << std::setw(OWID) << x1
        << std::setw(OWID) << y0
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << 0
        << std::setw(OWID) << -1
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y0     
        << std::setw(OWID) << z1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl << std::endl
      
      // boundary 6
        << std::setw(OWID) << 1 << std::endl
        << std::setw(OWID) << 6
        << std::setw(OWID) << 5 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0
        << std::setw(OWID) << -1
        << std::setw(OWID) << x0
        << std::setw(OWID) << y0
        << std::setw(OWID) << z1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0
        << std::setw(OWID) << x2
        << std::setw(OWID) << y0
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << -1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0
        << std::setw(OWID) << x1 
        << std::setw(OWID) << y0
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y2
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << -1
        << std::setw(OWID) << 0
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y1
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl;
  }
  else if (boundaryNum == 6){ // all 6 boundaries
    // boundary 1
    ofs << std::setw(OWID) << 1 << std::endl
        << std::setw(OWID) << 1
        << std::setw(OWID) << 5 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0     
        << std::setw(OWID) << x2
        << std::setw(OWID) << y0
        << std::setw(OWID) << z0     
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << -1
        << std::setw(OWID) << 0
        << std::setw(OWID) << x0    
        << std::setw(OWID) << y1
        << std::setw(OWID) << z0     
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y2
        << std::setw(OWID) << z0     
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0
        << std::setw(OWID) << 1
        << std::setw(OWID) << x0     
        << std::setw(OWID) << y0    
        << std::setw(OWID) << z2 
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 
        << std::setw(OWID) << -1
        << std::setw(OWID) << x0     
        << std::setw(OWID) << y0     
        << std::setw(OWID) << z1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl << std::endl
      
      // boundary 2
        << std::setw(OWID) << 1 << std::endl
        << std::setw(OWID) << 2
        << std::setw(OWID) << 5 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0     
        << std::setw(OWID) << x0    
        << std::setw(OWID) << y2
        << std::setw(OWID) << z0     
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << 0
        << std::setw(OWID) << x2
        << std::setw(OWID) << y0
        << std::setw(OWID) << z0     
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << -1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0
        << std::setw(OWID) << x1 
        << std::setw(OWID) << y0
        << std::setw(OWID) << z0     
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0
        << std::setw(OWID) << 1
        << std::setw(OWID) << x0     
        << std::setw(OWID) << y0    
        << std::setw(OWID) << z2 
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 
        << std::setw(OWID) << -1
        << std::setw(OWID) << x0     
        << std::setw(OWID) << y0      
        << std::setw(OWID) << z1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl << std::endl
      
      // boundary 3
        << std::setw(OWID) << 1 << std::endl
        << std::setw(OWID) << 3
        << std::setw(OWID) << 5 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << -1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0     
        << std::setw(OWID) << x1
        << std::setw(OWID) << y0
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << -1
        << std::setw(OWID) << 0
        << std::setw(OWID) << x0     
        << std::setw(OWID) << y1
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0  
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << x0       
        << std::setw(OWID) << y2
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0
        << std::setw(OWID) << 1
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y0     
        << std::setw(OWID) << z2 
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 
        << std::setw(OWID) << -1
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y0      
        << std::setw(OWID) << z1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl << std::endl
      
      // boundary 4
        << std::setw(OWID) << 1 << std::endl
        << std::setw(OWID) << 4
        << std::setw(OWID) << 5 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << -1
        << std::setw(OWID) << 0     
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y1
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << 0
        << std::setw(OWID) << x2
        << std::setw(OWID) << y0
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << -1 
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0
        << std::setw(OWID) << x1 
        << std::setw(OWID) << y0
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0
        << std::setw(OWID) << 1
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y0     
        << std::setw(OWID) << z2 
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 
        << std::setw(OWID) << -1
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y0      
        << std::setw(OWID) << z1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl << std::endl
      
      // boundary 5
        << std::setw(OWID) << 1 << std::endl
        << std::setw(OWID) << 5
        << std::setw(OWID) << 5 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << 0
        << std::setw(OWID) << 1     
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y0
        << std::setw(OWID) << z2 
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << 0
        << std::setw(OWID) << x2
        << std::setw(OWID) << y0
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << -1 
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0
        << std::setw(OWID) << x1 
        << std::setw(OWID) << y0
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y2
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << -1
        << std::setw(OWID) << 0
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y1
        << std::setw(OWID) << z0
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl << std::endl
      
      // boundary 6
        << std::setw(OWID) << 1 << std::endl
        << std::setw(OWID) << 6
        << std::setw(OWID) << 5 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << 0
        << std::setw(OWID) << -1    
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y0
        << std::setw(OWID) << z1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << 0
        << std::setw(OWID) << x2
        << std::setw(OWID) << y0
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << -1 
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0
        << std::setw(OWID) << x1 
        << std::setw(OWID) << y0
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y2
        << std::setw(OWID) << z0      
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl
      
        << std::setw(OWID) << 1
        << std::setw(OWID) << 0
        << std::setw(OWID) << -1
        << std::setw(OWID) << 0 
        << std::setw(OWID) << x0      
        << std::setw(OWID) << y1
        << std::setw(OWID) << z0
        << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::endl << std::endl;
  }
  
  ofs.close();
}

// boundaryNum = 6 by default
void Assembly::buildBoundary(const char *boundaryFile)
{
  std::ofstream ofs(boundaryFile);
  if(!ofs) { std::cout << "stream error!" << std::endl; exit(-1);}
  ofs.setf(std::ios::scientific, std::ios::floatfield);

  int boundaryNum = 6;

  Vec  v1 = allContainer.getMinCorner();
  Vec  v2 = allContainer.getMaxCorner();
  Vec  v0 = allContainer.getCenter();
  REAL x1 = v1.getX();
  REAL y1 = v1.getY();
  REAL z1 = v1.getZ();
  REAL x2 = v2.getX();
  REAL y2 = v2.getY();
  REAL z2 = v2.getZ();
  REAL x0 = v0.getX();
  REAL y0 = v0.getY();
  REAL z0 = v0.getZ();
  
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z1
      << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z2 << std::endl
      << std::setw(OWID) << boundaryNum << std::endl << std::endl;

  // boundary 1
  ofs << std::setw(OWID) << 1 << std::endl
      << std::setw(OWID) << 1
      << std::setw(OWID) << 5 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0     
      << std::setw(OWID) << x2
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0     
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << -1
      << std::setw(OWID) << 0
      << std::setw(OWID) << x0    
      << std::setw(OWID) << y1
      << std::setw(OWID) << z0     
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y2
      << std::setw(OWID) << z0     
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << 1
      << std::setw(OWID) << x0     
      << std::setw(OWID) << y0    
      << std::setw(OWID) << z2 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 
      << std::setw(OWID) << -1
      << std::setw(OWID) << x0     
      << std::setw(OWID) << y0     
      << std::setw(OWID) << z1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl << std::endl
    
    // boundary 2
      << std::setw(OWID) << 1 << std::endl
      << std::setw(OWID) << 2
      << std::setw(OWID) << 5 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0     
      << std::setw(OWID) << x0    
      << std::setw(OWID) << y2
      << std::setw(OWID) << z0     
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 0
      << std::setw(OWID) << x2
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0     
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << -1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << x1 
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0     
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << 1
      << std::setw(OWID) << x0     
      << std::setw(OWID) << y0    
      << std::setw(OWID) << z2 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 
      << std::setw(OWID) << -1
      << std::setw(OWID) << x0     
      << std::setw(OWID) << y0      
      << std::setw(OWID) << z1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl << std::endl
    
    // boundary 3
      << std::setw(OWID) << 1 << std::endl
      << std::setw(OWID) << 3
      << std::setw(OWID) << 5 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << -1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0     
      << std::setw(OWID) << x1
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << -1
      << std::setw(OWID) << 0
      << std::setw(OWID) << x0     
      << std::setw(OWID) << y1
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0  
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << x0       
      << std::setw(OWID) << y2
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << 1
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y0     
      << std::setw(OWID) << z2 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 
      << std::setw(OWID) << -1
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y0      
      << std::setw(OWID) << z1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl << std::endl
    
    // boundary 4
      << std::setw(OWID) << 1 << std::endl
      << std::setw(OWID) << 4
      << std::setw(OWID) << 5 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << -1
      << std::setw(OWID) << 0     
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y1
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 0
      << std::setw(OWID) << x2
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << -1 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << x1 
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << 1
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y0     
      << std::setw(OWID) << z2 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 
      << std::setw(OWID) << -1
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y0      
      << std::setw(OWID) << z1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl << std::endl
    
    // boundary 5
      << std::setw(OWID) << 1 << std::endl
      << std::setw(OWID) << 5
      << std::setw(OWID) << 5 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 1     
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y0
      << std::setw(OWID) << z2 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 0
      << std::setw(OWID) << x2
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << -1 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << x1 
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y2
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << -1
      << std::setw(OWID) << 0
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y1
      << std::setw(OWID) << z0
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl << std::endl
    
    // boundary 6
      << std::setw(OWID) << 1 << std::endl
      << std::setw(OWID) << 6
      << std::setw(OWID) << 5 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 0
      << std::setw(OWID) << -1    
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y0
      << std::setw(OWID) << z1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 0
      << std::setw(OWID) << x2
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << -1 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << x1 
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y2
      << std::setw(OWID) << z0      
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl
    
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << -1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y1
      << std::setw(OWID) << z0
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 << std::endl << std::endl; 

  ofs.close();
}

} // namespace dem ends


/*
// particleLayers:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
// ht- how many times of size would be the floating height
void Assembly::generateParticle(gradation& grad,
			const char *ParticleFile,
			int particleLayers,
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
    if (particleLayers == 0) {      // just one free particle
	newptcl = new Particle(particleNum+1, 0, Vec(dimn/2/40,dimn/2/20,dimn/2), grad, young, poisson);
	particleVec.push_back(newptcl);
	particleNum++;
    }
    else if (particleLayers == 1) { // a horizontal layer of free particles
	z=dimn/2;
	for (x=-dimn/2*(grid-1)/10; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
	    for (y=-dimn/2*(grid-1)/10; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5){
		newptcl = new Particle(particleNum+1, 0, Vec(x,y,z), grad, young, poisson);
		particleVec.push_back(newptcl);
		particleNum++;
	    }
    }
    else if (particleLayers == 2) { // multiple layers of free particles
	REAL offset=0; // 0 for ellipsoids; dimn/2/5/5 for spheres
	if (grad.ratio_ba==1.0 && grad.ratio_ca==1.0)
	    offset = dimn/2/5/5;
	REAL z0 = -dimn/2*9/10 ;// dimn/2;
	for (z=z0; z<z0 + dimn*ht; z+=dimn/2/5) {
	//for (z=-dimn/2*4/5; z<dimn/2 + dimn*ht; z+=dimn/2/10) { // spheres
	    for (x=-dimn/2*(grid-1)/10+offset; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
		for (y=-dimn/2*(grid-1)/10+offset; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5){
		    newptcl = new Particle(particleNum+1, 0, Vec(x,y,z), grad, young, poisson);
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
				int   particleLayers,
				REAL rsize,
				int   totalSteps,  
				int   snapNum,
				int   interval,
				const char *iniptclfile,   
				const char *ParticleFile, 
				const char *contactfile,
				const char *progressfile, 
				const char *debugfile)
{
    if (grad.rorc == 1) {
	RORC = grad.rorc;
	container.setCenter(Vec(0,0,0));
	container.setDimx(grad.dimn);
	container.setDimy(grad.dimn);
	container.setDimz(grad.dimn);
	
	generate_p(grad, iniptclfile, particleLayers, rsize, 4.0);
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
// particleLayers:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
// ht- how many times of size would be the floating height
void Assembly::generate_p(gradation&  grad,
			 const char *ParticleFile,
			 int particleLayers,
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
	    newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99, young, poisson);
	    particleVec.push_back(newptcl);
	    particleNum++;
	}

    // particle boundary 2
    y=dimn/2*(grid+1)/10;
    for (x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5)
	for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5){
	    newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99, young, poisson);
	    particleVec.push_back(newptcl);
	    particleNum++;
	}

    // particle boundary 3
    x=-dimn/2*(grid+1)/10;
    for (y=-dimn/2*grid/10; y<dimn/2*grid/10*est; y+=dimn/2/5)
	for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5){
	    newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99, young, poisson);
	    particleVec.push_back(newptcl);
	    particleNum++;
	}

    // particle boundary 4
    y=-dimn/2*(grid+1)/10;
    for (x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5)
	for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5){
	    newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99, young, poisson);
	    particleVec.push_back(newptcl);
	    particleNum++;
	}

    // particle boundary 6
    z=-dimn/2;
    for (y=-dimn/2*grid/10; y<dimn/2*grid/10*est; y+=dimn/2/5)
	for( x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5){
	    newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99, young, poisson);
	    particleVec.push_back(newptcl);
	    particleNum++;
	}

    if (particleLayers == 0) {      // just one free particle
	newptcl = new Particle(particleNum+1, 0, Vec(dimn/2/40,dimn/2/20,dimn/2), grad, young, poisson);
	particleVec.push_back(newptcl);
	particleNum++;
    }
    else if (particleLayers == 1) { // a horizontal layer of free particles
	z=dimn/2;
	for (x=-dimn/2*(grid-1)/10; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
	    for (y=-dimn/2*(grid-1)/10; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5){
		newptcl = new Particle(particleNum+1, 0, Vec(x,y,z), grad, young, poisson);
		particleVec.push_back(newptcl);
		particleNum++;
	    }
    }
    else if (particleLayers == 2) { // multiple layers of free particles
	for (z=dimn/2; z<dimn/2 + dimn*ht; z+=dimn/2/5)
	    for (x=-dimn/2*(grid-1)/10; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
		for (y=-dimn/2*(grid-1)/10; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5){
		    newptcl = new Particle(particleNum+1, 0, Vec(x,y,z), grad, young, poisson);
		    particleVec.push_back(newptcl);
		    particleNum++;
		}	
    }
    
    printParticle(ParticleFile);
    
}
*/

   /*
void Assembly::boundaryForce(){
  std::vector<Boundary*>::iterator rt;
  for(rt=boundaryVec.begin();rt!=boundaryVec.end();++rt)
    (*rt)->boundaryForce(boundaryTgtMap);

  
  std::vector<boundarytgt>::iterator it;
  std::vector<Boundary*>::iterator rt;

  for(rt=boundaryVec.begin();rt!=boundaryVec.end();++rt){	
    (*rt)->boundaryForce(boundaryTgtMap);
    for (it=boundaryTgtMap[(*rt)->bdry_id].begin();it!=boundaryTgtMap[(*rt)->bdry_id].end();++it){
      debugInf << std::setw(OWID) << iteration
		 << std::setw(OWID) << (*rt)->bdry_id
		 << std::setw(OWID) << boundaryTgtMap[(*rt)->bdry_id].size()
		 << std::setw(OWID) << it->TgtForce.getX()
		 << std::setw(OWID) << it->TgtForce.getY()
		 << std::setw(OWID) << it->TgtForce.getZ()
		 << std::endl;
      // << std::setw(OWID) << it->TgtPeak << std::endl;
    }
  }

  }*/
