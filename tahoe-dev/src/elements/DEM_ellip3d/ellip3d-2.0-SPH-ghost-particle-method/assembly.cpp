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

#include "assembly.h"
#include "parameter.h"
#include "timefunc.h"
#include "ran.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <cassert>
#include <utility>
#include <cmath>
#include <stdlib.h>
#include <math.h>

#ifdef OPENMP
#include <omp.h>
#define OPENMP_IMPL 2
#endif

// id is the id of thread, n is problem size, p is the number of threads

#define BLOCK_LOW(id, p, n)    ((id)*(n)/(p))

#define BLOCK_HIGH(id, p, n)   (BLOCK_LOW((id)+1, p, n)-1)

#define BLOCK_SIZE(id, p, n)   (BLOCK_LOW((id)+1, p, n)-BLOCK_LOW(id, p, n))

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
  
void assembly::printParticle(const char* str) const {
  std::ofstream ofs(str);
  if(!ofs) {
    cout << "stream error!" << endl; exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  ofs << setw(OWID) << TotalNum << setw(OWID) << 1 << endl;
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
//      << setw(OWID) << "initial_x"	// for testing granular strain
//      << setw(OWID) << "initial_y"
//      << setw(OWID) << "initial_z"
      << endl;
  
  vec tmp;
  std::vector<particle*>::const_iterator  it;
  for (it=ParticleVec.begin();it!=ParticleVec.end();++it)  {
    ofs << setw(OWID) << (*it)->getID()
	<< setw(OWID) << (*it)->getType()
	<< setw(OWID) << (*it)->getA()
	<< setw(OWID) << (*it)->getB()
	<< setw(OWID) << (*it)->getC();
    
    tmp=(*it)->getCurrPosition();
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
    
    tmp=(*it)->getCurrVelocity();
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
    
//    tmp=(*it)->getInitPosition();
//    ofs << setw(OWID) << tmp.getx()
//	<< setw(OWID) << tmp.gety()
//	<< setw(OWID) << tmp.getz() << endl;
  }
  
  ofs.close();
}


void assembly::plotBoundary(const char *str) const {
  std::ofstream ofs(str);
  if(!ofs) {
    cout << "stream error!" << endl; exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  REAL x1,x2,y1,y2,z1,z2,x0,y0,z0;
  x1 = container.getMinCorner().getx();
  y1 = container.getMinCorner().gety();
  z1 = container.getMinCorner().getz();
  x2 = container.getMaxCorner().getx();
  y2 = container.getMaxCorner().gety();
  z2 = container.getMaxCorner().getz();

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


void assembly::plotCavity(const char *str) const {
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


void assembly::plotSpring(const char *str) const {
  std::ofstream ofs(str);
  if(!ofs) {
    cout << "stream error!" << endl; exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  int totalMemParticle = 0;
  for (int i = 0; i < MemBoundary.size(); ++i) 
    for (int j = 0; j < MemBoundary[i].size(); ++j) 
      for (int k = 0; k < MemBoundary[i][j].size(); ++k) 
	++totalMemParticle;
  int totalSpring = SpringVec.size();
  ofs << "ZONE N=" << totalMemParticle << ", E=" << totalSpring << ", DATAPACKING=POINT, ZONETYPE=FELINESEG" << endl;
  particle *pt = NULL;
  vec vt;
  for (int i = 0; i < MemBoundary.size(); ++i) 
    for (int j = 0; j < MemBoundary[i].size(); ++j) 
      for (int k = 0; k < MemBoundary[i][j].size(); ++k) {
	pt = MemBoundary[i][j][k]; 
	vt = pt->getCurrPosition();
	ofs << setw(OWID) << vt.getx() << setw(OWID) << vt.gety() << setw(OWID) << vt.getz() << endl;
      }
  for (int i = 0; i < SpringVec.size(); ++i) {
    ofs << setw(OWID) << SpringVec[i]->getParticleId1() - trimHistoryNum  << setw(OWID) << SpringVec[i]->getParticleId2() - trimHistoryNum << endl;
  }

  ofs.close();
}  

void assembly::printMemParticle(const char* str) const  {
  std::ofstream ofs(str);
  if(!ofs) {
    cout << "stream error!" << endl; exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  
  int totalMemParticle = 0;
  for (int i = 0; i < MemBoundary.size(); ++i) 
    for (int j = 0; j < MemBoundary[i].size(); ++j) 
      for (int k = 0; k < MemBoundary[i][j].size(); ++k) 
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
  
  particle *it = NULL;
  vec tmp;
  for (int i = 0; i < MemBoundary.size(); ++i) 
    for (int j = 0; j < MemBoundary[i].size(); ++j) 
      for (int k = 0; k < MemBoundary[i][j].size(); ++k) {
	it = MemBoundary[i][j][k];
	ofs << setw(OWID) << it->getID()
	    << setw(OWID) << it->getType()
	    << setw(OWID) << it->getA()
	    << setw(OWID) << it->getB()
	    << setw(OWID) << it->getC();
	
	tmp=it->getCurrPosition();
	ofs << setw(OWID) << tmp.getx()
	    << setw(OWID) << tmp.gety()
	    << setw(OWID) << tmp.getz();
	
	tmp=it->getCurrDirecA();
	ofs << setw(OWID) << tmp.getx()
	    << setw(OWID) << tmp.gety()
	    << setw(OWID) << tmp.getz();
	
	tmp=it->getCurrDirecB();
	ofs << setw(OWID) << tmp.getx()
	    << setw(OWID) << tmp.gety()
	    << setw(OWID) << tmp.getz();
	
	tmp=it->getCurrDirecC();
	ofs << setw(OWID) << tmp.getx()
	    << setw(OWID) << tmp.gety()
	    << setw(OWID) << tmp.getz();
	
	tmp=it->getCurrVelocity();
	ofs << setw(OWID) << tmp.getx()
	    << setw(OWID) << tmp.gety()
	    << setw(OWID) << tmp.getz();
	
	tmp=it->getCurrOmga();
	ofs << setw(OWID) << tmp.getx()
	    << setw(OWID) << tmp.gety()
	    << setw(OWID) << tmp.getz();
	
	tmp=it->getForce();
	ofs << setw(OWID) << tmp.getx()
	    << setw(OWID) << tmp.gety()
	    << setw(OWID) << tmp.getz();
	
	tmp=it->getMoment();
	ofs << setw(OWID) << tmp.getx()
	    << setw(OWID) << tmp.gety()
	      << setw(OWID) << tmp.getz() << endl;
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
void assembly::checkMembrane(vector<REAL> &vx ) const {
  vector<particle*> vec1d;  // 1-dimension
  vector< vector<particle*>  > vec2d; // 2-dimension
  REAL in, out, tmp;
  REAL x1_in, x1_out, x2_in, x2_out;
  REAL y1_in, y1_out, y2_in, y2_out;
  REAL z1_in, z1_out, z2_in, z2_out;

  // surface x1
  vec2d = MemBoundary[0];
  in = vec2d[0][0]->getCurrPosition().getx();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPosition().getx();
      if (tmp < out) out = tmp;
      if (tmp > in ) in  = tmp;
    }
  vx.push_back(in);
  vx.push_back(out);
  x1_in  = in;
  x1_out = out;

  // surface x2
  vec2d.clear();
  vec2d = MemBoundary[1];
  in = vec2d[0][0]->getCurrPosition().getx();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPosition().getx();
      if (tmp > out) out = tmp;
      if (tmp < in ) in  = tmp;
    }
  vx.push_back(in);
  vx.push_back(out);
  x2_in  = in;
  x2_out = out;

  // surface y1
  vec2d.clear();
  vec2d = MemBoundary[2];
  in = vec2d[0][0]->getCurrPosition().gety();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPosition().gety();
      if (tmp < out) out = tmp;
      if (tmp > in ) in  = tmp;
    }
  vx.push_back(in);
  vx.push_back(out);
  y1_in  = in;
  y1_out = out;

  // surface y2
  vec2d.clear();
  vec2d = MemBoundary[3];
  in = vec2d[0][0]->getCurrPosition().gety();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPosition().gety();
      if (tmp > out) out = tmp;
      if (tmp < in ) in  = tmp;
    }
  vx.push_back(in);
  vx.push_back(out);
  y2_in  = in;
  y2_out = out;
  
  // surface z1
  vec2d.clear();
  vec2d = MemBoundary[4];
  in = vec2d[0][0]->getCurrPosition().getz();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPosition().getz();
      if (tmp < out) out = tmp;
      if (tmp > in ) in  = tmp;
    }
  vx.push_back(in);
  vx.push_back(out);
  z1_in  = in;
  z1_out = out;

  // surface z2
  vec2d.clear();
  vec2d = MemBoundary[5];
  in = vec2d[0][0]->getCurrPosition().getz();
  out= in;
  for (int i = 0; i < vec2d.size(); ++i)
    for (int j = 0; j < vec2d[i].size(); ++j) {
      tmp = vec2d[i][j]->getCurrPosition().getz();
      if (tmp > out) out = tmp;
      if (tmp < in ) in  = tmp;
    }
  vx.push_back(in);
  vx.push_back(out);
  z2_in  = in;
  z2_out = out;

}
  
void assembly::printRectPile(const char* str)
{
    std::ofstream ofs(str, std::ios_base::app);
    if(!ofs) {
	cout << "stream error!" << endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);

    ofs << setw(OWID) << 8 << setw(OWID) << 6 << endl;
    vec pos[8];
    for(std::vector<RGDBDRY*>::iterator rt=RBVec.begin();rt!=RBVec.end();++rt){
	if((*rt)->getBdryID()==7){
	    pos[0]=vec((*rt)->CoefOfLimits[0].apt.getx(),
		       (*rt)->CoefOfLimits[1].apt.gety(),
		       (*rt)->CoefOfLimits[4].apt.getz());
	    pos[1]=vec((*rt)->CoefOfLimits[0].apt.getx(),
		       (*rt)->CoefOfLimits[2].apt.gety(),
		       (*rt)->CoefOfLimits[4].apt.getz());
	    pos[5]=vec((*rt)->CoefOfLimits[0].apt.getx(),
		       (*rt)->CoefOfLimits[2].apt.gety(),
		       (*rt)->CoefOfLimits[3].apt.getz());
	    pos[4]=vec((*rt)->CoefOfLimits[0].apt.getx(),
		       (*rt)->CoefOfLimits[1].apt.gety(),
		       (*rt)->CoefOfLimits[3].apt.getz());
	}
	else if((*rt)->getBdryID()==9) {
	    pos[2]=vec((*rt)->CoefOfLimits[0].apt.getx(),
		       (*rt)->CoefOfLimits[2].apt.gety(),
		       (*rt)->CoefOfLimits[4].apt.getz());
	    pos[3]=vec((*rt)->CoefOfLimits[0].apt.getx(),
		       (*rt)->CoefOfLimits[1].apt.gety(),
		       (*rt)->CoefOfLimits[4].apt.getz());
	    pos[7]=vec((*rt)->CoefOfLimits[0].apt.getx(),
		       (*rt)->CoefOfLimits[1].apt.gety(),
		       (*rt)->CoefOfLimits[3].apt.getz());
	    pos[6]=vec((*rt)->CoefOfLimits[0].apt.getx(),
		       (*rt)->CoefOfLimits[2].apt.gety(),
		       (*rt)->CoefOfLimits[3].apt.getz());
	}
    }

    for (int i=0;i<8;i++)
	ofs << setw(OWID) << pos[i].getx() << setw(OWID) << pos[i].gety() << setw(OWID) << pos[i].getz() << endl;

    ofs << setw(OWID) << 1 << setw(OWID) << 2 << setw(OWID) << 6 << setw(OWID) << 5 << endl
        << setw(OWID) << 2 << setw(OWID) << 3 << setw(OWID) << 7 << setw(OWID) << 6 << endl
        << setw(OWID) << 3 << setw(OWID) << 4 << setw(OWID) << 8 << setw(OWID) << 7 << endl
        << setw(OWID) << 4 << setw(OWID) << 1 << setw(OWID) << 5 << setw(OWID) << 8 << endl
        << setw(OWID) << 1 << setw(OWID) << 4 << setw(OWID) << 3 << setw(OWID) << 2 << endl
        << setw(OWID) << 5 << setw(OWID) << 6 << setw(OWID) << 7 << setw(OWID) << 8 << endl;

    ofs.close();
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
void assembly::printContact(const char* str) const
{
    std::ofstream ofs(str);
    if(!ofs) {
	cout << "stream error!" << endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    ofs << setw(OWID) << ActualCntctNum << endl;
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
    for (it=ContactVec.begin();it!=ContactVec.end();++it)
	ofs << setw(OWID) << it->getP1()->getID()
	    << setw(OWID) << it->getP2()->getID()
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

	
void assembly::readSample(const char* str){
    std::ifstream ifs(str);
    if(!ifs) {
	cout << "stream error!" << endl; exit(-1);
    }
    int RORC;
    ifs >> TotalNum >> RORC;

    REAL cx,cy,cz,dx,dy,dz;

    ifs >> cx >> cy >> cz >> dx >> dy >> dz;
    container.set(dx,dy,dz,vec(cx,cy,cz));
    Volume = container.getDimx() * container.getDimy() * container.getDimz();
    
    char s[20];
    ifs>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
       >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s;

    ParticleVec.clear();
    int ID, type;
    REAL a, b, c, px,py,pz,dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz;
    REAL vx,vy,vz,omx,omy,omz,fx,fy,fz,mx,my,mz;
    for (int i=0;i<TotalNum;i++){
	ifs>>ID>>type>>a>>b>>c>>px>>py>>pz>>dax>>day>>daz>>dbx>>dby>>dbz>>dcx>>dcy>>dcz
	   >>vx>>vy>>vz>>omx>>omy>>omz>>fx>>fy>>fz>>mx>>my>>mz;
	particle* pt= new particle(ID,type,vec(a,b,c),vec(px,py,pz),vec(dax,day,daz),vec(dbx,dby,dbz),vec(dcx,dcy,dcz),YOUNG,POISSON);

//      optional settings for a particle's initial status
//	pt->setPrevVelocity(vec(vx,vy,vz));
//	pt->setCurrVelocity(vec(vx,vy,vz));
//	pt->setPrevOmga(vec(omx,omy,omz));
//	pt->setCurrOmga(vec(omx,omy,omz));
//	pt->setConstForce(vec(fx,fy,fz));  // constant force, not initial force
//	pt->setConstMoment(vec(mx,my,mz)); // constant moment, not initial moment

	ParticleVec.push_back(pt);
    }
    ifs.close();
}

// read tessellation information from Qhull output file for granular strain calculation
// written on Feb 18, 2013
void assembly::readTesse(const char* str){
	std::ifstream ifs(str);
	std::cout << "Read tessellation begin!" << std::endl;
    	if(!ifs) {
		std::cout << "stream error!" << std::endl; exit(-1);
   	}
	int m, n, i, j;
	int totalNum;	// total number of cells, ie pyramids
	cell tempCell;
	edge tempEdge;
	ifs>>totalNum;
	edge_map.clear();
	for(int it=0; it!=totalNum; it++){
		ifs>>m>>n>>i>>j;	// read nodes in each line
		m = m+1;
		n = n+1;
		i = i+1;
		j = j+1;	//the ID from Qhull is starting from 0

		tempCell.setNodes(m,n,i,j);
		// edge mn
		tempEdge.setNodes(m,n);
		edge_map[tempEdge].push_back(tempCell);
		// edge mi
		tempEdge.setNodes(m,i);
		edge_map[tempEdge].push_back(tempCell);
		// edge mj
		tempEdge.setNodes(m,j);
		edge_map[tempEdge].push_back(tempCell);
		// edge ni
		tempEdge.setNodes(n,i);
		edge_map[tempEdge].push_back(tempCell);
		// edge nj
		tempEdge.setNodes(n,j);
		edge_map[tempEdge].push_back(tempCell);
		// edge ij
		tempEdge.setNodes(i,j);
		edge_map[tempEdge].push_back(tempCell);
	}
}

// read cell information into std::vector<cell> cellVec from Qhull output file, for finite granular strain calculation
// written on March 26, 2013
void assembly::readTesse_finite(const char* str){
	std::ifstream ifs(str);
	std::cout << "Read cell information begin!" << std::endl;
	if(!ifs) {
		std::cout << "stream error!" << std::endl; exit(-1);
	}
	int m, n, i, j;
	int totalNum;

	ifs>>totalNum;
	cellVec.clear();
	for(int it=0; it!=totalNum; it++){
		ifs>>m>>n>>i>>j;	// read nodes in each line
		m = m+1;
		n = n+1;
		i = i+1;
		j = j+1;	// the ID from Qhull is starting from 0
		cell* ptCell = new cell();
		ptCell->setNodes(m,n,i,j);
		cellVec.push_back(ptCell);
	}
	setNumberingOrder();	// this is important
}

//start of def OPENMP 
#ifdef OPENMP	

#if OPENMP_IMPL == 0
// OpenMP implementation 0: ts partitions, each thread handles a partition, max diff = n*n*(1-1/ts)/ts
// implementation is based on linked list, also works for vector but not efficient.
void assembly::findContact() { 
  ContactVec.clear();
  int possContact = 0;
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p1,NULL); 
#endif

  int tid;   // thread id
  int ts;    // number of threads
  int num;   // number of particles
  int tnum;  // number of particles per thread
  int i, j;
  vec u, v;
  num = ParticleVec.size();
  ot = ParticleVec.begin();
  std::vector<particle*>::iterator ot, it, pt;
  
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
      u=(*it)->getCurrPosition();
      for (pt = it, ++pt; pt != ParticleVec.end(); ++pt){
	v=(*pt)->getCurrPosition();
	if (   ( vfabs(v-u) < (*it)->getA() + (*pt)->getA())
	       && ( (*it)->getType() !=  1 || (*pt)->getType() != 1  )      // not both are fixed particles
	       && ( (*it)->getType() !=  5 || (*pt)->getType() != 5  )      // not both are free boundary particles
	       && ( (*it)->getType() != 10 || (*pt)->getType() != 10 )  ) { // not both are ghost particles
	  contact<particle> tmpct(*it, *pt); // a local and temparory object
	  ++possContact;
	  if(tmpct.isOverlapped())
#pragma omp critical
	    ContactVec.push_back(tmpct);    // containers use value semantics, so a "copy" is pushed back.
	}
      }
    }
  }
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  g_debuginf << setw(OWID) << "findContact=" << setw(OWID) << timediffsec(time_p1, time_p2); 
#endif
  PossCntctNum   = possContact;
  ActualCntctNum = ContactVec.size();
} // end of OpenMP implementation 0

#elif OPENMP_IMPL == 1
// OpenMP implementation 1: ts partitions, each thread handles a partition, max diff = n*n*(1-1/ts)/ts
// implementation is based on vector index.
void assembly::findContact() { 
  ContactVec.clear();
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
  vec u, v;
  num = ParticleVec.size();
  
#pragma omp parallel num_threads(NUM_THREADS) private(tid, ts, start, end, i, j, u, v) shared(num) reduction(+: possContact)
  {
    tid = omp_get_thread_num();
    ts  = omp_get_num_threads();
    start = tid * num / ts;
    end   = (tid + 1) * num / ts - 1;
    
    // explore each partition
    for (i = start; i <= end; ++i) { 
      u = ParticleVec[i]->getCurrPosition();
      for (j = i + 1; j < num; ++j) {
	v = ParticleVec[j]->getCurrPosition();
	if (   ( vfabs(v-u) < ParticleVec[i]->getA() + ParticleVec[j]->getA() )
	       && ( ParticleVec[i]->getType() !=  1 || ParticleVec[j]->getType() != 1  )      // not both are fixed particles
	       && ( ParticleVec[i]->getType() !=  5 || ParticleVec[j]->getType() != 5  )      // not both are free boundary particles
	       && ( ParticleVec[i]->getType() != 10 || ParticleVec[j]->getType() != 10 )  ) { // not both are ghost particles
	  contact<particle> tmpct(ParticleVec[i], ParticleVec[j]); // a local and temparory object
	  ++possContact;
	  if(tmpct.isOverlapped())
#pragma omp critical
	    ContactVec.push_back(tmpct);    // containers use value semantics, so a "copy" is pushed back.
	}
      }
    }
  }
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  g_debuginf <<  setw(OWID) << "findContact=" << setw(OWID) << timediffsec(time_p1, time_p2); 
#endif
  PossCntctNum   = possContact;  
  ActualCntctNum = ContactVec.size();
} // end of OpenMP implementation 1

#elif OPENMP_IMPL == 2
// OpenMP implementation 2: no partitions, each thread leaps by ts until completed, max diff = n*(ts-1)/ts
void assembly::findContact() { 
  ContactVec.clear();
  int possContact = 0;
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p1,NULL); 
#endif
  int tid;   // thread id
  int ts;    // number of threads
  int num;   // number of particles
  int i, j;
  vec u, v;
  num = ParticleVec.size();
  
#pragma omp parallel num_threads(NUM_THREADS) private(tid, ts, i, j, u, v) shared(num) reduction(+: possContact)
  {
    tid = omp_get_thread_num();
    ts  = omp_get_num_threads();
    
    // explore each partition
    for (i = tid; i < num; i += ts) { 
      u = ParticleVec[i]->getCurrPosition();
      for (j = i + 1; j < num; ++j) {
	v = ParticleVec[j]->getCurrPosition();
	if (   ( vfabs(v-u) < ParticleVec[i]->getA() + ParticleVec[j]->getA() )
	       && ( ParticleVec[i]->getType() !=  1 || ParticleVec[j]->getType() != 1  )      // not both are fixed particles
	       && ( ParticleVec[i]->getType() !=  5 || ParticleVec[j]->getType() != 5  )      // not both are free boundary particles
	       && ( ParticleVec[i]->getType() != 10 || ParticleVec[j]->getType() != 10 )  ) { // not both are ghost particles
	  contact<particle> tmpct(ParticleVec[i], ParticleVec[j]); // a local and temparory object
	  ++possContact;
	  if(tmpct.isOverlapped())
#pragma omp critical
	    ContactVec.push_back(tmpct);    // containers use value semantics, so a "copy" is pushed back.
	}
      }
    }
  }
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  g_debuginf << setw(OWID) << "findContact=" << setw(OWID) << timediffsec(time_p1, time_p2); 
#endif
  PossCntctNum   = possContact;  
  ActualCntctNum = ContactVec.size();
} // end of OpenMP implementation 2

#elif OPENMP_IMPL == 3
// OpenMP implementation 3: no partitions, each thread leaps by ts until num/2 and handles two particles, max diff = 0
void assembly::findContact() { 
  ContactVec.clear();
  int possContact = 0;
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p1,NULL); 
#endif
  int tid;   // thread id
  int ts;    // number of threads
  int num;   // number of particles
  int i, j, k;
  vec u, v;
  num = ParticleVec.size();
  
#pragma omp parallel num_threads(NUM_THREADS) private(tid, ts, i, j, k, u, v) shared(num) reduction(+: possContact)
  {
    tid = omp_get_thread_num();
    ts  = omp_get_num_threads();
    
    // explore each partition, works whether num is odd or even
    for (i = tid; i <= num / 2; i += ts) {
      int inc = num - 1 - 2*i;
      if (inc == 0) inc = 1; // avoid infinite loop when num is odd
      for (k = i; k <= num - 1 - i; k += inc ) {
	u = ParticleVec[k]->getCurrPosition();
	for (j = k + 1; j < num; ++j) {
	  v = ParticleVec[j]->getCurrPosition();
	  if (   ( vfabs(v-u) < ParticleVec[k]->getA() + ParticleVec[j]->getA() )
		 && ( ParticleVec[k]->getType() !=  1 || ParticleVec[j]->getType() != 1  )      // not both are fixed particles
		 && ( ParticleVec[k]->getType() !=  5 || ParticleVec[j]->getType() != 5  )      // not both are free boundary particles
		 && ( ParticleVec[k]->getType() != 10 || ParticleVec[j]->getType() != 10 )  ) { // not both are ghost particles
	    contact<particle> tmpct(ParticleVec[k], ParticleVec[j]); // a local and temparory object
	    ++possContact;
	    if(tmpct.isOverlapped())
#pragma omp critical
	    ContactVec.push_back(tmpct);    // containers use value semantics, so a "copy" is pushed back.
	  }
	}
      }
    }
  }
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  g_debuginf << setw(OWID) << "findContact=" << setw(OWID) << timediffsec(time_p1, time_p2); 
#endif
  PossCntctNum   = possContact;  
  ActualCntctNum = ContactVec.size();
} // end of OpenMP implementation 3

#elif OPENMP_IMPL == 4
// OpenMP implementation 4: no partitions, parallel for, various loop scheduling: (static), (static,1), (dynamic), (dynamic,1)
void assembly::findContact() { 
  ContactVec.clear();
  int possContact = 0;
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p1,NULL); 
#endif

  int num;   // number of particles
  int i, j;
  vec u, v;
  num = ParticleVec.size();
  
#pragma omp parallel for num_threads(NUM_THREADS) private(i, j, u, v) shared(num) reduction(+: possContact) schedule(dynamic)
  for (i = 0; i < num - 1; ++i) { 
    u = ParticleVec[i]->getCurrPosition();
    for (j = i + 1; j < num; ++j) {
      v = ParticleVec[j]->getCurrPosition();
      if (   ( vfabs(v-u) < ParticleVec[i]->getA() + ParticleVec[j]->getA() )
	     && ( ParticleVec[i]->getType() !=  1 || ParticleVec[j]->getType() != 1  )      // not both are fixed particles
	     && ( ParticleVec[i]->getType() !=  5 || ParticleVec[j]->getType() != 5  )      // not both are free boundary particles
	     && ( ParticleVec[i]->getType() != 10 || ParticleVec[j]->getType() != 10 )  ) { // not both are ghost particles
	contact<particle> tmpct(ParticleVec[i], ParticleVec[j]); // a local and temparory object
	++possContact;
	if(tmpct.isOverlapped())
#pragma omp critical
	  ContactVec.push_back(tmpct);    // containers use value semantics, so a "copy" is pushed back.
      }
    }
  }
  
  
#ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  g_debuginf << setw(OWID) << "findContact=" << setw(OWID) << timediffsec(time_p1, time_p2); 
#endif
  PossCntctNum   = possContact;  
  ActualCntctNum = ContactVec.size();
} // end of OpenMP implementation 4

#endif

#else //else of def OPENMP, i.e., serial versions start here:

//start of ndef BINNING
#ifndef BINNING
void assembly::findContact(){ // serial version, O(n x n), n is the number of particles.
    ContactVec.clear();
    PossCntctNum = 0;

#ifdef TIME_PROFILE
    double time_r = 0; // time consumed in contact resolution, i.e., tmpct.isOverlapped()
    gettimeofday(&time_p1,NULL); 
#endif
    std::vector<particle*>::iterator it, pt;
    vec u,v;
    for (it=ParticleVec.begin();it!=ParticleVec.end();++it){
	u=(*it)->getCurrPosition();
	for (pt=it,++pt;pt!=ParticleVec.end();++pt){
	    v=(*pt)->getCurrPosition();
	    if (   ( vfabs(v-u) < (*it)->getA() + (*pt)->getA())
		&& ( (*it)->getType() !=  1 || (*pt)->getType() != 1  )      // not both are fixed particles
		&& ( (*it)->getType() !=  5 || (*pt)->getType() != 5  )      // not both are free boundary particles
		&& ( (*it)->getType() != 10 || (*pt)->getType() != 10 )  ) { // not both are ghost particles
		contact<particle> tmpct(*it, *pt); // a local and temparory object
		++PossCntctNum;
#ifdef TIME_PROFILE
		gettimeofday(&time_r1,NULL); 
#endif
		if(tmpct.isOverlapped())
		  ContactVec.push_back(tmpct);    // containers use value semantics, so a "copy" is pushed back.
#ifdef TIME_PROFILE
		gettimeofday(&time_r2,NULL); 
		time_r += timediffsec(time_r1, time_r2);
#endif
	    }
	}
    }	

#ifdef TIME_PROFILE
    gettimeofday(&time_p2,NULL);
    g_debuginf << setw(OWID) << "findContact=" << setw(OWID) << timediffsec(time_p1, time_p2) << setw(OWID) << "isOverlapped=" << setw(OWID) << time_r; 
#endif

    ActualCntctNum = ContactVec.size();
}

//else of ndef BINNING
#else
void assembly::findContact(){ // serial version, binning methods, cell slightly larger than maximum particle
  ContactVec.clear();
  PossCntctNum = 0;
  
#ifdef TIME_PROFILE
  double time_r = 0;
  gettimeofday(&time_p1,NULL); 
#endif
  REAL maxDiameter = gradInfo.getMaxPtclRadius() * 2;
  int  nx = floor (container.getDimx() / maxDiameter);
  int  ny = floor (container.getDimy() / maxDiameter);
  int  nz = floor (container.getDimz() *1.5 / maxDiameter);
  REAL dx = container.getDimx() / nx;
  REAL dy = container.getDimx() / ny;
  REAL dz = container.getDimx() *1.5 / nz;
  vec  minCorner= container.getMinCorner();
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
  typedef pair<bool, vector<particle*> > cellT;
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
  vec center;
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
	for (int pt = 0; pt < ParticleVec.size(); ++pt) {
	  center = ParticleVec[pt]->getCurrPosition();
	  if (center.getx() >= x1 && center.getx() < x2 &&
	      center.gety() >= y1 && center.gety() < y2 &&
	      center.getz() >= z1 && center.getz() < z2)
	    cellVec[i][j][k].second.push_back( ParticleVec[pt] );
	}
      }
  
  // for each cell:
  particle *it, *pt;
  vec u, v;
  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      for (int k = 0; k < nz; ++k) {
	// for particles inside the cell	  
	for (int m = 0; m < cellVec[i][j][k].second.size(); ++m) {
	  it = cellVec[i][j][k].second[m];
	  u  = it->getCurrPosition();
	  
	  // for particles inside the cell itself   
	  for (int n = m + 1; n < cellVec[i][j][k].second.size(); ++n) {
	    //cout <<  i << " " << j << " " << k << " " << "m n size=" << m << " " << n << " " <<  cellVec[i][j][k].size() << endl;
	    pt = cellVec[i][j][k].second[n];
	    v  = pt->getCurrPosition();
	    if ( ( vfabs(u-v) < it->getA() + pt->getA() )  &&
		 ( it->getType() !=  1 || pt->getType() != 1 ) &&   // not both are fixed particles
		 ( it->getType() !=  5 || pt->getType() != 5 ) &&   // not both are free boundary particles
		 ( it->getType() != 10 || pt->getType() != 10)  ) { // not both are ghost particles
	      contact<particle> tmpct(it, pt); // a local and temparory object
	      ++PossCntctNum;
#ifdef TIME_PROFILE
		gettimeofday(&time_r1,NULL); 
#endif
	      if(tmpct.isOverlapped())
		ContactVec.push_back(tmpct);   // containers use value semantics, so a "copy" is pushed back.
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
	      //cout << "i j k m ncell ci cj ck size contacts= " << i << " " << j << " " << k << " " << m  << " " << ncell << " " << ci << " " << cj << " " << ck << " " << cellVec[ci][cj][ck].second.size() << " "  << ContactVec.size() << endl;
	      vector<particle*> vt = cellVec[ci][cj][ck].second;
	      for (int n = 0; n < vt.size(); ++n) {
		pt = vt[n];
		v  = pt->getCurrPosition();
		if ( ( vfabs(u-v) < it->getA() + pt->getA() )  &&
		     ( it->getType() !=  1 || pt->getType() != 1 ) &&   // not both are fixed particles
		     ( it->getType() !=  5 || pt->getType() != 5 ) &&   // not both are free boundary particles
		     ( it->getType() != 10 || pt->getType() != 10)  ) { // not both are ghost particles
		  contact<particle> tmpct(it, pt); // a local and temparory object
		  ++PossCntctNum;
#ifdef TIME_PROFILE
		gettimeofday(&time_r1,NULL); 
#endif
		  if(tmpct.isOverlapped())
		    ContactVec.push_back(tmpct);   // containers use value semantics, so a "copy" is pushed back.
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
  g_debuginf << setw(OWID) << "findContact=" << setw(OWID) << timediffsec(time_p1, time_p2) << setw(OWID) << "isOverlapped=" << setw(OWID) << time_r; 
#endif
  
  ActualCntctNum = ContactVec.size();
}

//end of ndef BINNING
#endif

//end of def OPENMP 
#endif

REAL assembly::getDensity() const{
    REAL dens=0;
    std::vector<particle*>::const_iterator it;
    for(it=ParticleVec.begin();it!=ParticleVec.end();++it)
	dens+=(*it)->getMass();
    return dens/=Volume;
}


REAL assembly::getAveragePenetration() const{
    int totalcntct = ContactVec.size();
    if (totalcntct==0)
	return 0;
    else {
	REAL pene=0;
	for (std::vector<CONTACT>::const_iterator it=ContactVec.begin();it!=ContactVec.end();++it)
	    pene += it->getPenetration(); 
	return pene/totalcntct;
    }
}


REAL assembly::getVibraTimeStep() const {
    int totalcntct = ContactVec.size();
    if (totalcntct == 0)
	return 0;
    else {
	std::vector<CONTACT>::const_iterator it=ContactVec.begin();
        REAL minTimeStep = it->getVibraTimeStep();
	for (++it; it != ContactVec.end(); ++it) {
	  REAL val = it->getVibraTimeStep(); 
	  minTimeStep =  val < minTimeStep ? val : minTimeStep;
	}
	return minTimeStep;
    }
}


REAL assembly::getImpactTimeStep() const {
    int totalcntct = ContactVec.size();
    if (totalcntct == 0)
	return 0;
    else {
	std::vector<CONTACT>::const_iterator it=ContactVec.begin();
        REAL minTimeStep = it->getImpactTimeStep();
	for (++it; it != ContactVec.end(); ++it) {
	  REAL val = it->getImpactTimeStep(); 
	  minTimeStep =  val < minTimeStep ? val : minTimeStep;
	}
	return minTimeStep;
    }
}
 

REAL assembly::getAverageVelocity() const{
    REAL avgv=0;
    int count=0;
    std::vector<particle*>::const_iterator it;
    for(it=ParticleVec.begin();it!=ParticleVec.end();++it)
	if ((*it)->getType()==0) {
	    avgv+=vfabs((*it)->getCurrVelocity());
	    count++;
	}
    return avgv/=count;
}


REAL assembly::getAverageOmga() const{
    REAL avgv=0;
    int count=0;
    std::vector<particle*>::const_iterator it;
    for(it=ParticleVec.begin();it!=ParticleVec.end();++it)
	if ((*it)->getType()==0){
	    avgv+=vfabs((*it)->getCurrOmga());
	    count++;
	}
    return avgv/=count;
}


REAL assembly::getAverageForce() const{
    REAL avgv=0;
    int count=0;
    std::vector<particle*>::const_iterator it;
    for(it=ParticleVec.begin();it!=ParticleVec.end();++it)
	if ((*it)->getType()==0){
	    avgv+=vfabs((*it)->getForce());
	    count++;
	}
    return avgv/count;
}


REAL assembly::getAverageMoment() const{
    REAL avgv=0;
    int count=0;
    std::vector<particle*>::const_iterator it;
    for(it=ParticleVec.begin();it!=ParticleVec.end();++it)
	if ((*it)->getType()==0){
	    avgv+=vfabs((*it)->getMoment());
	    count++;
	}
    return avgv/=count;
}


REAL assembly::getParticleVolume() const{
    REAL avgv=0;
    std::vector<particle*>::const_iterator it;
    for(it=ParticleVec.begin();it!=ParticleVec.end();++it)
	if ((*it)->getType()==0)
	    avgv+=(*it)->getVolume();
    return avgv;
}


vec assembly::getTopFreeParticlePosition() const{
    std::vector<particle*>::const_iterator it,jt,kt;
    it=ParticleVec.begin();
    while (it!=ParticleVec.end() && (*it)->getType()!=0)   // find the 1st free particle
	++it;

    if (it==ParticleVec.end())    // no free particles
	return 0;

    jt=it; 
    kt=it;
    
    // two cases:
    // 1: 1st particle is not free
    // 2: 1st particle is free
    if (++kt!=ParticleVec.end()){ // case1: more than 2 particles; case 2: more than 1 particle
	for(++it;it!=ParticleVec.end();++it){
	    if ((*it)->getType()==0)
		if ((*it)->getCurrPosition().getz() > (*jt)->getCurrPosition().getz())
		    jt=it;
	}
	return (*jt)->getCurrPosition();
    }
    else {
	if ((*it)->getType()==0)  // case1: only 2 particles, the 2nd one is free; case2: only 1 particle
	    return (*it)->getCurrPosition();
	else
	    return 0;
    }

}

REAL assembly::getMaxCenterHeight() const {
  std::vector<particle*>::const_iterator it = ParticleVec.begin();
  REAL z0 = (*it)->getCurrPosition().getz();
  for (; it != ParticleVec.end(); ++it) {
    if ( (*it)->getCurrPosition().getz() > z0 )
      z0 = (*it)->getCurrPosition().getz();
  }
  return z0;
}


REAL assembly::ellipPileForce() {
    REAL val=0;
    for(std::vector<particle*>::iterator it=ParticleVec.begin();it!=ParticleVec.end();++it)
	if ((*it)->getType()==3) {
	    val = (*it)->getForce().getz();
	    break;
	}
    return val;
}


vec assembly::ellipPileDimn() {
    vec val;
    for(std::vector<particle*>::iterator it=ParticleVec.begin();it!=ParticleVec.end();++it)
	if ((*it)->getType()==3) {
	    val = vec((*it)->getA(), (*it)->getB(), (*it)->getC());
	    break;
	}
    return val;
}


REAL assembly::ellipPileTipZ() {
    REAL val=0;
    for(std::vector<particle*>::iterator it=ParticleVec.begin();it!=ParticleVec.end();++it)
	if ((*it)->getType()==3) {
	    val = (*it)->getCurrPosition().getz()-(*it)->getA();
	    break;
	}
    return val;
}


REAL assembly::ellipPilePeneVol() {
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


void assembly::ellipPileUpdate(){
    for(std::vector<particle*>::iterator it=ParticleVec.begin();it!=ParticleVec.end();++it){
	if ((*it)->getType()==3) {
	    (*it)->curr_velocity.setx(0);	
	    (*it)->curr_velocity.sety(0);
	    (*it)->curr_velocity.setz(-PILE_RATE);
	    (*it)->curr_position = (*it)->prev_position + (*it)->curr_velocity*TIMESTEP;
	}
    }
}


REAL assembly::getTransEnergy() const{
    REAL engy=0;
    std::vector<particle*>::const_iterator it;
    for(it=ParticleVec.begin();it!=ParticleVec.end();++it){
	if ((*it)->getType()==0)
	    engy+=(*it)->getTransEnergy();
    }
    return engy;
}


REAL assembly::getRotatEnergy() const{
    REAL engy=0;
    std::vector<particle*>::const_iterator it;
    for(it=ParticleVec.begin();it!=ParticleVec.end();++it){
	if ((*it)->getType()==0)
	    engy+=(*it)->getRotatEnergy();
    }
    return engy;
}


REAL assembly::getKinetEnergy() const{
    REAL engy=0;
    std::vector<particle*>::const_iterator it;
    for(it=ParticleVec.begin();it!=ParticleVec.end();++it){
	if ((*it)->getType()==0)
	    engy+=(*it)->getKinetEnergy();
    }
    return engy;
}


REAL assembly::getPotenEnergy(REAL ref) const{
    REAL engy=0;
    std::vector<particle*>::const_iterator it;
    for(it=ParticleVec.begin();it!=ParticleVec.end();++it){
	if ((*it)->getType()==0)
	    engy+=(*it)->getPotenEnergy(ref);
    }
    return engy;
}


void assembly::clearForce(){
    for(std::vector<particle*>::iterator it=ParticleVec.begin();it!=ParticleVec.end();++it){
	(*it)->clearForce();
    }
}


void assembly::flexiBoundaryForceZero(){
    for(std::vector<particle*>::iterator it=ParticleVec.begin();it!=ParticleVec.end();++it){
	(*it)->flb_force=0;
	(*it)->flb_moment=0;
    }
}


void assembly::initFBForce(){
    for(std::vector<particle*>::iterator it=ParticleVec.begin();it!=ParticleVec.end();++it){
	(*it)->force+=(*it)->flb_force;
	(*it)->moment+=(*it)->flb_moment;
    }
}

void assembly::springForce() {
  vector<spring*> :: iterator it;
  for (it = SpringVec.begin(); it != SpringVec.end(); ++it) {
    (*it)->applyForce();
  }
}

void assembly::internalForce(REAL& avgnm, REAL& avgsh){
    avgnm=0;
    avgsh=0;

    int totalcntct = ContactVec.size();
    if(totalcntct==0){
	avgnm = 0;
	avgsh = 0;
    }
    else{
	std::vector<CONTACT>::iterator it;
	for (it=ContactVec.begin();it!=ContactVec.end();++it)
	    it->checkinPreTgt(CntTgtVec); // checkin previous tangential force and displacment    
	
	CntTgtVec.clear(); // CntTgtVec must be cleared before filling in new values.

#ifdef TIME_PROFILE
	gettimeofday(&time_p1,NULL); 
#endif 
	for (it=ContactVec.begin();it!=ContactVec.end();++it){
            bool exceed = false;
	    it->contactForce(exceed);     // cannot be parallelized as it may change a particle's force simultaneously.
	    it->checkoutTgt(CntTgtVec);   // checkout current tangential force and displacment
	    avgnm += it->getNormalForce();
	    avgsh += it->getTgtForce();
#ifndef NDEBUG
	    if (exceed) {
	      char stepsstr[7];
	      char stepsfp[50];
	      sprintf(stepsstr, "%06d", g_iteration);
	      strcpy(stepsfp,"particle_");
	      strcat(stepsfp, stepsstr);
	      printParticle(stepsfp);
	    }
#endif
	}
	avgnm /= totalcntct;
	avgsh /= totalcntct;

#ifdef TIME_PROFILE
	gettimeofday(&time_p2,NULL);
	g_debuginf << setw(OWID) << "internalForce=" << setw(OWID) << timediffsec(time_p1, time_p2) << endl; 
#endif

    }
}


void assembly::updateParticle(){
    for(std::vector<particle*>::iterator it=ParticleVec.begin();it!=ParticleVec.end();++it){
	(*it)->update();
    }
}


void assembly::readBoundary(const char* str){
    std::ifstream ifs(str);
    if(!ifs) {
	cout << "stream error!" << endl; exit(-1);
    }
    ifs >> BdryType;
    if(BdryType==0){      // rigid boundaries
	readRigidBoundary(ifs);
    }
    else if(BdryType==1){ // flexible boundaries
	readFlexiBoundary(ifs);
    }
    ifs.close();
}


void assembly::readRigidBoundary(std::ifstream &ifs){
    rgd_bdry<particle>* rbptr;
    int type;
    RBVec.clear();
    ifs>>RgdBdryNum;
    
    for(int i=0;i<RgdBdryNum;i++){
	ifs>>type;
	if(type==1) // plane boundary
	    rbptr=new plnrgd_bdry<particle>(ifs);
	else        // cylindrical boundary
	    rbptr=new cylrgd_bdry<particle>(ifs);
	RBVec.push_back(rbptr);
    }
}


void assembly::readFlexiBoundary(std::ifstream &ifs){
    flb_bdry<particle>* fbptr;
    int type;
    FBVec.clear();
    ifs>>FlbBdryNum;
    
    for(int i=0;i<FlbBdryNum;i++){
	ifs>>type;
	if(type==1) // plane boundary
	    fbptr=new plnflb_bdry<particle>(ifs);
	else        // cylindrical boundary
	    fbptr=new cylflb_bdry<particle>(ifs);
	FBVec.push_back(fbptr);
    }
}


void assembly::readCavityBoundary(const char* str){
    std::ifstream ifs(str);
    if(!ifs) {
	cout << "stream error!" << endl; exit(-1);
    }
    ifs >> BdryType;
    if(BdryType == 0){      // rigid boundaries
      rgd_bdry<particle>* rbptr;
      int type;
      CavityRBVec.clear();
      ifs>>RgdBdryNum;
      
      for(int i=0;i<RgdBdryNum;i++){
	ifs>>type;
	if(type==1) // plane boundary
	  rbptr=new plnrgd_bdry<particle>(ifs);
	CavityRBVec.push_back(rbptr);
      }
    }

    ifs.close();
}


void assembly::printBoundary(const char* str) const
{
    std::ofstream ofs(str);
    if(!ofs) {
	cout << "stream error!" << endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);

    ofs << setw(OWID) << BdryType
        << setw(OWID) << RgdBdryNum << endl;
    
    std::vector<RGDBDRY*>::const_iterator rt;
    for(rt=RBVec.begin();rt!=RBVec.end();++rt)
	(*rt)->disp(ofs);
    ofs << endl;

    ofs.close();
}


void assembly::printCavityBoundary(const char* str) const
{
    std::ofstream ofs(str);
    if(!ofs) {
	cout << "stream error!" << endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);

    ofs << setw(OWID) << BdryType
        << setw(OWID) << RgdBdryNum << endl;
    
    std::vector<RGDBDRY*>::const_iterator rt;
    for(rt = CavityRBVec.begin(); rt != CavityRBVec.end(); ++rt)
	(*rt)->disp(ofs);
    ofs << endl;

    ofs.close();
}


void assembly::findParticleOnBoundary(){
    std::vector<RGDBDRY*>::iterator rt;
    std::vector<FLBBDRY*>::iterator ft;
    for(rt=RBVec.begin();rt!=RBVec.end();++rt)
	(*rt)->findParticleOnBoundary(ParticleVec);
    for(ft=FBVec.begin();ft!=FBVec.end();++ft)
	(*ft)->findParticleOnBoundary(ParticleVec);
}


void assembly::findParticleOnCavity(){
    std::vector<RGDBDRY*>::iterator rt;
    for(rt = CavityRBVec.begin(); rt != CavityRBVec.end(); ++rt)
	(*rt)->findParticleOnBoundary(ParticleVec);
}

void assembly::findParticleOnLine(){
    std::vector<FLBBDRY*>::iterator ft;
    for(ft=FBVec.begin();ft!=FBVec.end();++ft)
	(*ft)->findParticleOnLine();
}


void assembly::createFlbNet(){
    std::vector<FLBBDRY*>::iterator ft;
    for(ft=FBVec.begin();ft!=FBVec.end();++ft)
	(*ft)->createFlbNet();
}


void assembly::rigidBoundaryForce(){
  std::vector<RGDBDRY*>::iterator rt;
  for(rt=RBVec.begin();rt!=RBVec.end();++rt)
    (*rt)->rigidBF(BdryTgtMap);

  /*
  vector<boundarytgt>::iterator it;
  std::vector<RGDBDRY*>::iterator rt;

  for(rt=RBVec.begin();rt!=RBVec.end();++rt){	
    (*rt)->rigidBF(BdryTgtMap);
    for (it=BdryTgtMap[(*rt)->bdry_id].begin();it!=BdryTgtMap[(*rt)->bdry_id].end();++it){
      g_debuginf << setw(OWID) << g_iteration
		 << setw(OWID) << (*rt)->bdry_id
		 << setw(OWID) << BdryTgtMap[(*rt)->bdry_id].size()
		 << setw(OWID) << it->TgtForce.getx()
		 << setw(OWID) << it->TgtForce.gety()
		 << setw(OWID) << it->TgtForce.getz()
		 << endl;
      // << setw(OWID) << it->TgtPeak << endl;
    }
  }
  */
}


void assembly::rigidBoundaryForce(REAL penetr[],int cntnum[]){
  std::vector<RGDBDRY*>::iterator rt;
  for(rt=RBVec.begin();rt!=RBVec.end();++rt){	
    (*rt)->rigidBF(BdryTgtMap);
    for (int i = 1; i <= 6; ++i) {
      if ((*rt)->getBdryID() == i){
	penetr[i] = (*rt)->getAvgPenetr();
	cntnum[i] = (*rt)->getCntnum();
	break;
      }
    }
  }
}


void assembly::cavityBoundaryForce(){
  std::vector<RGDBDRY*>::iterator rt;
  for(rt = CavityRBVec.begin(); rt != CavityRBVec.end(); ++rt)
    (*rt)->rigidBF(BdryTgtMap);
}


void assembly::flexiBoundaryForce(){
    std::vector<FLBBDRY*>::iterator ft;
    for(ft=FBVec.begin();ft!=FBVec.end();++ft)
	(*ft)->flxbBF();
}


vec assembly::getNormalForce(int bdry) const{
    std::vector<RGDBDRY*>::const_iterator it;
    for(it=RBVec.begin();it!=RBVec.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getNormalForce();
    }
    return 0;
}

// calculate fabric  tensor Fij=1/Nc*sum(ni*nj), written on March 8, 2013. the unit of fabric tensor is 1 after normalization
matrix assembly::getFabric() const{
	matrix fabric(3,3);
	matrix n(3,1);	// contact vector
	matrix P1center(3,1), P2center(3,1);
	int Nc = 0;	// contact number
	REAL n_magnitude;
	// initialize fabric tensor
	for(int ir=0; ir!=3; ir++)
		for(int ic=0; ic!=3; ic++)
			fabric(ir+1,ic+1) = 0;

	// sum the total particle contact force
 	std::vector<CONTACT>::const_iterator it;

	for (it=ContactVec.begin();it!=ContactVec.end();++it){
		P1center(1,1) = it->getP1()->getCurrPosition().getx();
		P1center(2,1) = it->getP1()->getCurrPosition().gety();
		P1center(3,1) = it->getP1()->getCurrPosition().getz();

		P2center(1,1) = it->getP2()->getCurrPosition().getx();
		P2center(2,1) = it->getP2()->getCurrPosition().gety();
		P2center(3,1) = it->getP2()->getCurrPosition().getz();

		n = P2center-P1center;
		n_magnitude = sqrt(n(1,1)*n(1,1)+n(2,1)*n(2,1)+n(3,1)*n(3,1));
		n = n/n_magnitude;	// normalize vector n
		Nc++;		
		fabric += n*n.getTrans();

	}
std::cout << "Nc: " << Nc << " AtualConctNum: " << ActualCntctNum << std::endl;	// used to debug
	fabric = fabric/Nc;

	return fabric;
}

// reset initial position for each particle
void assembly::resetStartPosition() {
	vec curr_p;
	std::vector<particle*>::const_iterator it;
	for(it=ParticleVec.begin();it!=ParticleVec.end();++it){
		curr_p = (*it)->getCurrPosition();
		(*it)->setStartPosition(curr_p);
	}
}

// create input file for Qhull
void assembly::createInputForQhull() const{
	std::ofstream ofs("input_for_Qhull");
	if(!ofs){
		cout << "Stream error when create input for Qhull!" << endl;
		exit (-1);
	}
	ofs.setf(std::ios::scientific, std::ios::floatfield);
  	ofs.precision(OPREC);
	ofs << setw(OWID) << "3" << endl
	    << setw(OWID) << TotalNum << endl;
	vec tmp;
  	std::vector<particle*>::const_iterator  it;
  	for (it=ParticleVec.begin();it!=ParticleVec.end();++it)  {
    		tmp=(*it)->getCurrPosition();
    		ofs << setw(OWID) << tmp.getx()
		    << setw(OWID) << tmp.gety()
		    << setw(OWID) << tmp.getz() << endl;
	}
/*
	REAL x_max = 5.0e-3;
	REAL y_max = 5.0e-3;
	REAL z_max = 3.0e-3;
	REAL num_zero = 0;
	ofs << setw(OWID) << 0
	    << setw(OWID) << y_max
    	    << setw(OWID) << z_max << endl;
	ofs << setw(OWID) << x_max
	    << setw(OWID) << y_max
    	    << setw(OWID) << z_max << endl;
	ofs << setw(OWID) << x_max
	    << setw(OWID) << 0
    	    << setw(OWID) << z_max << endl;
*/
    	ofs.close();
}

// call Qhull to tessellate and ouput "tess_info"
void assembly::callQhull() const{
	std::cout << "Checking if processor is available..." << std::endl;
	if(system(NULL)) std::cout << "Ok!" << std::endl;
	else exit(EXIT_FAILURE);
	system("./qdelaunay Qt i < input_for_Qhull TO tess_info");	// call the external command qdelaunay
}

// calculate granular stress, written on Feb 13, 2013
matrix assembly::getGranularStress() const{
	matrix totalForce(3,3);
	matrix sigma(3,3);
	matrix lc(3,1), Fc(1,3);	// sigma(i,j) = lixFj
	matrix P1center(3,1), P2center(3,1);

	// initialize totalForce
	for(int ir=0; ir!=3; ir++)
		for(int ic=0; ic!=3; ic++)
			totalForce(ir+1,ic+1) = 0;

	// below is to calculate the granular stress

	// sum the total particle contact force
 	std::vector<CONTACT>::const_iterator it;

	for (it=ContactVec.begin();it!=ContactVec.end();++it){
		P1center(1,1) = it->getP1()->getCurrPosition().getx();
		P1center(2,1) = it->getP1()->getCurrPosition().gety();
		P1center(3,1) = it->getP1()->getCurrPosition().getz();

		P2center(1,1) = it->getP2()->getCurrPosition().getx();
		P2center(2,1) = it->getP2()->getCurrPosition().gety();
		P2center(3,1) = it->getP2()->getCurrPosition().getz();

		lc = P2center-P1center;
		Fc(1,1) = (it->NormalForceVec().getx())+(it->TgtForceVec().getx());
		Fc(1,2) = (it->NormalForceVec().gety())+(it->TgtForceVec().gety());
		Fc(1,3) = (it->NormalForceVec().getz())+(it->TgtForceVec().getz());

		totalForce += lc*Fc;

	}
	sigma = totalForce/Volume;
	return sigma;

// above is to calculate granular stress

}

// granular strain claculation
// written on Feb 18, 2013
matrix assembly::getum(int ID, int type) const{	
// type = 1 means calculate um from starting position, for Bagi's strain and Eulerian strain; type = 0 means calculate um from
// initial position, for lagrangian finite strain.
	if(type != 0 && type != 1){std::cout << "Error: type should be 0 or 1 calling getum()!" << std::endl; exit(-1);}
	std::vector<particle*>::const_iterator it;
	it = ParticleVec.begin();
	for(int i=0; i!=ID-1; i++)
		it++;	// go to the ID particle
	matrix init_xyz(3,1);	// initial center
	if(type == 0){
		init_xyz(1,1) = ((*it)->getInitPosition()).getx();
		init_xyz(2,1) = ((*it)->getInitPosition()).gety();
		init_xyz(3,1) = ((*it)->getInitPosition()).getz();
	}
	else if (type == 1){
		init_xyz(1,1) = ((*it)->getStartPosition()).getx();
		init_xyz(2,1) = ((*it)->getStartPosition()).gety();
		init_xyz(3,1) = ((*it)->getStartPosition()).getz();
	}
	matrix curr_xyz(3,1);	// current center
	curr_xyz(1,1) = ((*it)->getCurrPosition()).getx();
	curr_xyz(2,1) = ((*it)->getCurrPosition()).gety();
	curr_xyz(3,1) = ((*it)->getCurrPosition()).getz();
	matrix u(3,1);	// displacement vector
	u = curr_xyz-init_xyz;
	return u;
}

matrix assembly::getb(cell pyra, int ID_m) const{
	if(ID_m!=pyra.getm() && ID_m!=pyra.getn() && ID_m!=pyra.geti() && ID_m!=pyra.getj()){
		std::cout << "Error when calculate b vector: the node is not in this cell!" << std::endl;
		exit(-1);
	}
	std::vector<int> node_ID;	// the other three nodes in this cell
	std::vector<int>::size_type it_ID = 0;
	if(pyra.getm() != ID_m)
		node_ID.push_back(pyra.getm());
	if(pyra.getn() != ID_m)
		node_ID.push_back(pyra.getn());
	if(pyra.geti() != ID_m)
		node_ID.push_back(pyra.geti());
	if(pyra.getj() != ID_m)
		node_ID.push_back(pyra.getj());
	
	std::vector<particle*>::const_iterator it;
	REAL x1, y1, z1;	// the current position of first node
	it = ParticleVec.begin();
	for(int i=0; i!=node_ID[it_ID]-1; i++)
		it++;	// go to the ID particle
	x1 = ((*it)->getInitPosition()).getx();
	y1 = ((*it)->getInitPosition()).gety();
	z1 = ((*it)->getInitPosition()).getz();
	it_ID++;	// go to the next node
	REAL x2, y2, z2;	// the current position of second node
	it = ParticleVec.begin();
	for(int i=0; i!=node_ID[it_ID]-1; i++)
		it++;	// go to the ID particle
	x2 = ((*it)->getInitPosition()).getx();
	y2 = ((*it)->getInitPosition()).gety();
	z2 = ((*it)->getInitPosition()).getz();
	it_ID++;	// go to the next node
	REAL x3, y3, z3;	// the current position of third node
	it = ParticleVec.begin();
	for(int i=0; i!=node_ID[it_ID]-1; i++)
		it++;	// go to the ID particle
	x3 = ((*it)->getInitPosition()).getx();
	y3 = ((*it)->getInitPosition()).gety();
	z3 = ((*it)->getInitPosition()).getz();

	// calculate area
	REAL area_tri;
	area_tri = 0.5*sqrt( pow((y2-y1)*(z3-z1)-(z2-z1)*(y3-y1),2)
		   +pow((z2-z1)*(x3-x1)-(x2-x1)*(z3-z1),2)+pow((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1),2));
//	area_tri = fabs(area_tri);	// to make sure it is positve since from the formula it is possible that the area is negative
::
	// calculate normal vector
	REAL nx, ny, nz;
	nx = (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1);
	ny = (z2-z1)*(x3-x1)-(x2-x1)*(z3-z1);
	nz = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);
	// to check normal vector is outward or not
	REAL xm, ym, zm;	// the current position of node ID_m
	it = ParticleVec.begin();
	for(int i=0; i!=ID_m-1; i++)
		it++;	// go to the ID particle
	xm = ((*it)->getInitPosition()).getx();
	ym = ((*it)->getInitPosition()).gety();
	zm = ((*it)->getInitPosition()).getz();
	REAL vec_mx, vec_my, vec_mz;
	vec_mx = xm-x3;	// points from node 3 to node m
	vec_my = ym-y3;
	vec_mz = zm-z3;
	if( nx*vec_mx+ny*vec_my+nz*vec_mz > 0) {	// means the n vector is inward
		nx = -nx;
		ny = -ny;
		nz = -nz;
	}
	REAL mag_n = sqrt(nx*nx+ny*ny+nz*nz);	// magnitude of normal vector
	REAL bx, by, bz;	// three components of b vector
	bx = nx/mag_n*area_tri;
	by = ny/mag_n*area_tri;
	bz = nz/mag_n*area_tri;
	matrix b_vec(3,1);
	b_vec(1,1) = bx;
	b_vec(2,1) = by;
	b_vec(3,1) = bz;
	return b_vec;
}

matrix assembly::getdmn(edge mn, std::vector<cell> cellsContain) const{
	int node_m, node_n;	// nodes of edge mn
	node_m = mn.getm();
	node_n = mn.getn();
	matrix dmn(3,1);
	dmn(1,1) = 0;
	dmn(2,1) = 0;
	dmn(3,1) = 0;
	for(std::vector<cell>::const_iterator it=cellsContain.begin(); it!=cellsContain.end(); it++){
		dmn += getb((*it), node_n)-getb((*it), node_m);
	}
	dmn = dmn/12;
	return dmn;
}

matrix assembly::getGranularStrain() {
	int node_m, node_n;
	matrix deltau, dmn;
	edge temp;
	matrix granular_e(3,3), granular_strain(3,3);
	// initialize
	for(int i=0; i!=3; i++){
		for(int j=0; j!=3; j++){
			granular_e(i+1,j+1) = 0;
			granular_strain(i+1,j+1) = 0;
		}
	}
	average_dudx_Bagi = granular_e;	// used to test quadratic terms, April 22, 2013
//std::cout << "point 1!" << std::endl;
	for(std::map<edge, std::vector<cell> >::const_iterator map_it=edge_map.begin(); map_it!=edge_map.end(); map_it++){
		// get the edge
		temp = map_it->first;
		dmn.clear();
		dmn = getdmn(temp, (map_it->second));
		// get the node m,n
		node_m = temp.getm();
		node_n = temp.getn();	// m < n
		deltau.clear();
		deltau = getum(node_m,1)-getum(node_n,1);
		// 1 means calculate um from starting position, for Bagi's strain and Eulerian strain
		granular_e += deltau*dmn.getTrans();
	}
	granular_e = granular_e/Volume;
	average_dudx_Bagi = granular_e;	// used to test quadratic terms
	// get granular strain
	granular_strain = (granular_e+granular_e.getTrans())/2;
	return granular_strain;
}

// above is for granular strain calculation

// finite granular strain calculation, March 26, 2013
void assembly::setNumberingOrder() {
	int nm, nn, ni, nj;	// 4 nodes of a tetrahedron
	int n1, n2 , n3, n4;	// the final odering node 
	std::vector<particle*>::const_iterator it;
//	REAL x1, y1, z1;	// the current position of node 1
//	REAL x2, y2, z2;	// the current position of node 2
//	REAL x3, y3, z3;	// the current position of node 3
//	REAL x4, y4, z4;	// the current position of node 4
	REAL cell_volume;	// area of triangle, volume of tet
	for(std::vector<cell*>::const_iterator iter=cellVec.begin(); iter!=cellVec.end(); iter++){
		nm = (*iter)->getm();
		nn = (*iter)->getn();
		ni = (*iter)->geti();
		nj = (*iter)->getj();
		// 1, choose m, n as node 1, 2
		n1 = nm;
		n2 = nn;
		n3 = ni;
		n4 = nj;
		(*iter)->setNodes(n1,n2,n3,n4);
		cell_volume = getCellVolume(**iter);

		if(cell_volume < 0){	// swap n2 and n3
			n2 = ni;
			n3 = nn;
			(*iter)->setNodes(n1,n2,n3,n4);
			cell_volume = getCellVolume(**iter);
		}
		
//		// get coordinates for node 1, 2
//		it = ParticleVec.begin();
//		for(int i=0; i!=n1-1; i++)
//			it++;	// go to the n1 particle
//		x1 = ((*it)->getCurrPosition()).getx();
//		y1 = ((*it)->getCurrPosition()).gety();
//		z1 = ((*it)->getCurrPosition()).getz();
//
//		it = ParticleVec.begin();
//		for(int i=0; i!=n2-1; i++)
//			it++;	// go to the n2 particle
//		x2 = ((*it)->getCurrPosition()).getx();
//		y2 = ((*it)->getCurrPosition()).gety();
//		z2 = ((*it)->getCurrPosition()).getz();
//		
//		// 2, choose i as node 3 and test if Area n1,n2,n3 is positive
//		n3 = ni;
//		n4 = nj;
//		// get coordinates for node 3
//		it = ParticleVec.begin();
//		for(int i=0; i!=n3-1; i++)
//			it++;	// go to the n3 particle
//		x3 = ((*it)->getCurrPosition()).getx();
//		y3 = ((*it)->getCurrPosition()).gety();
//		z3 = ((*it)->getCurrPosition()).getz();
//		
//		// calculate area
//		area = 0.5*( x1*y2+x3*y1+x2*y3-x3*y2-x1*y3-x2*y1);
//		if(area < 0){	// means the previous three nodes are not numbered counter-clockwise
//			n3 = nj;
//			n4 = ni;
//			// get coordinates for node 3
//			it = ParticleVec.begin();
//			for(int i=0; i!=n3-1; i++)
//				it++;	// go to the n3 particle
//			x3 = ((*it)->getCurrPosition()).getx();
//			y3 = ((*it)->getCurrPosition()).gety();
//			z3 = ((*it)->getCurrPosition()).getz();
//		}
//		// get coordinates for node 4
//		it = ParticleVec.begin();
//		for(int i=0; i!=n4-1; i++)
//			it++;	// go to the n3 particle
//		x4 = ((*it)->getCurrPosition()).getx();
//		y4 = ((*it)->getCurrPosition()).gety();
//		z4 = ((*it)->getCurrPosition()).getz();
//		// test if volume is positive
//		cell_volume = x1*y2*z3+x2*y3*z4+x3*y4*z1+x4*y1*z2
//			     -x1*y4*z3-x2*y1*z4-x3*y2*z1-x4*y3*z2;
		if(cell_volume < 0){	// means n1,n2,n3,n4 are not numbered counter-clockwise
			std::cout << "Error: nodes are not numbered counter-clockwise!" << std::endl;
		}

		(*iter)->setInitialCellVolume(cell_volume);
		matrix bigB(4,3);
		bigB = getBigB(**iter);
		(*iter)->setInitialBigB(bigB);
	}
}

matrix assembly::getBigB(const cell& tempCell) const{
	int n1,n2,n3,n4;
	REAL x1, y1, z1;	// the current position of node 1
	REAL x2, y2, z2;	// the current position of node 2
	REAL x3, y3, z3;	// the current position of node 3
	REAL x4, y4, z4;	// the current position of node 4
	// get ID for 4 nodes
	n1 = tempCell.getm();
	n2 = tempCell.getn();
	n3 = tempCell.geti();
	n4 = tempCell.getj();

	// get coordinates for the 4 nodes
	std::vector<particle*>::const_iterator it;
	
	it = ParticleVec.begin();
	for(int i=0; i!=n1-1; i++)
		it++;	// go to the n1 particle
	x1 = ((*it)->getCurrPosition()).getx();
	y1 = ((*it)->getCurrPosition()).gety();
	z1 = ((*it)->getCurrPosition()).getz();
	
	it = ParticleVec.begin();
	for(int i=0; i!=n2-1; i++)
		it++;	// go to the n2 particle
	x2 = ((*it)->getCurrPosition()).getx();
	y2 = ((*it)->getCurrPosition()).gety();
	z2 = ((*it)->getCurrPosition()).getz();

	it = ParticleVec.begin();
	for(int i=0; i!=n3-1; i++)
		it++;	// go to the n3 particle
	x3 = ((*it)->getCurrPosition()).getx();
	y3 = ((*it)->getCurrPosition()).gety();
	z3 = ((*it)->getCurrPosition()).getz();

	it = ParticleVec.begin();
	for(int i=0; i!=n4-1; i++)
		it++;	// go to the n4 particle
	x4 = ((*it)->getCurrPosition()).getx();
	y4 = ((*it)->getCurrPosition()).gety();
	z4 = ((*it)->getCurrPosition()).getz();

	// get intermediate variables
	REAL a1, a2, a3, a4;
	REAL b1, b2, b3, b4;
	REAL c1, c2, c3, c4;
	a1 = y2*(z4-z3)-y3*(z4-z2)+y4*(z3-z2);
	a2 = -y1*(z4-z3)+y3*(z4-z1)-y4*(z3-z1);
	a3 = y1*(z4-z2)-y2*(z4-z1)+y4*(z2-z1);
	a4 = -y1*(z3-z2)+y2*(z3-z1)-y3*(z2-z1);

	b1 = -x2*(z4-z3)+x3*(z4-z2)-x4*(z3-z2);
	b2 = x1*(z4-z3)-x3*(z4-z1)+x4*(z3-z1);
	b3 = -x1*(z4-z2)+x2*(z4-z1)-x4*(z2-z1);
	b4 = x1*(z3-z2)-x2*(z3-z1)+x3*(z2-z1);

	c1 = x2*(y4-y3)-x3*(y4-y2)+x4*(y3-y2);
	c2 = -x1*(y4-y3)+x3*(y4-y1)-x4*(y3-y1);
	c3 = x1*(y4-y2)-x2*(y4-y1)+x4*(y2-y1);
	c4 = -x1*(y3-y2)+x2*(y3-y1)-x3*(y2-y1);

	// assemble bigB
	matrix bigB(4,3);
	bigB(1,1) = a1;
	bigB(2,1) = a2;
	bigB(3,1) = a3;
	bigB(4,1) = a4;

	bigB(1,2) = b1;
	bigB(2,2) = b2;
	bigB(3,2) = b3;
	bigB(4,2) = b4;
	
	bigB(1,3) = c1;
	bigB(2,3) = c2;
	bigB(3,3) = c3;
	bigB(4,3) = c4;
	
	return bigB;
}

// dudx for F=dudx+I
matrix assembly::getdudx(const cell& tempCell) const{
	int n1,n2,n3,n4;

	// get ID for 4 nodes
	n1 = tempCell.getm();
	n2 = tempCell.getn();
	n3 = tempCell.geti();
	n4 = tempCell.getj();
	// get displacement for 4 nodes
	matrix u1(3,1), u2(3,1), u3(3,1), u4(3,1);
	u1 = getum(n1,0);
	u2 = getum(n2,0);
	u3 = getum(n3,0);	// 0 means calculate um from initial position, for lagrangian finite strain
	u4 = getum(n4,0);

	// get volume of this cell
	REAL cell_volume;
	cell_volume = tempCell.getInitialCellVolume();
	if(cell_volume < 0){
		std::cout << "Error: volume is negative in getdudx()!" << std::endl;
	}

	matrix bigU(3,4), bigB(4,3);
	// assemble bigU
	bigU(1,1) = u1(1,1);
	bigU(2,1) = u1(2,1);
	bigU(3,1) = u1(3,1);

	bigU(1,2) = u2(1,1);
	bigU(2,2) = u2(2,1);
	bigU(3,2) = u2(3,1);

	bigU(1,3) = u3(1,1);
	bigU(2,3) = u3(2,1);
	bigU(3,3) = u3(3,1);

	bigU(1,4) = u4(1,1);
	bigU(2,4) = u4(2,1);
	bigU(3,4) = u4(3,1);

// test assembly of bigU
//std::cout << "u1: " << std::endl;
//std::cout << u1(1,1) << " " << u1(2,1) << " " << u1(3,1) << std::endl;
//std::cout << "u2: " << std::endl;
//std::cout << u2(1,1) << " " << u2(2,1) << " " << u2(3,1) << std::endl;
//std::cout << "u3: " << std::endl;
//std::cout << u3(1,1) << " " << u3(2,1) << " " << u3(3,1) << std::endl;
//
//std::cout << "bigU: " << std::endl;
//std::cout << bigU(1,1) << " " << bigU(1,2) << " " << bigU(1,3) << std::endl; 
//std::cout << bigU(2,1) << " " << bigU(2,2) << " " << bigU(2,3) << std::endl; 
//std::cout << bigU(3,1) << " " << bigU(3,2) << " " << bigU(3,3) << std::endl; 
//std::cout << bigU(4,1) << " " << bigU(4,2) << " " << bigU(4,3) << std::endl; 
//

	bigB = tempCell.getInitialBigB();
	matrix result;
	result = 1/(6.0*cell_volume)*bigU*bigB;
	return result;
}

// dudx for F=dudx+I with respect to current configuration
matrix assembly::getdudx_curr(const cell& tempCell) const{
	int n1,n2,n3,n4;

	// get ID for 4 nodes
	n1 = tempCell.getm();
	n2 = tempCell.getn();
	n3 = tempCell.geti();
	n4 = tempCell.getj();
	// get displacement for 4 nodes
	matrix u1(3,1), u2(3,1), u3(3,1), u4(3,1);
	u1 = getum(n1,1);
	u2 = getum(n2,1);
	u3 = getum(n3,1);
	u4 = getum(n4,1);	// 1 means calculate um from starting position, for Eulerian strain and Bagi's strain

	// get volume of this cell
	REAL cell_volume;
	cell_volume = getCellVolume(tempCell);
	if(cell_volume < 0){
		std::cout << "Error: volume is negative in getdudx_curr()!" << std::endl;
	}

	matrix bigU(3,4), bigB(4,3);
	// assemble bigU
	bigU(1,1) = u1(1,1);
	bigU(2,1) = u1(2,1);
	bigU(3,1) = u1(3,1);

	bigU(1,2) = u2(1,1);
	bigU(2,2) = u2(2,1);
	bigU(3,2) = u2(3,1);

	bigU(1,3) = u3(1,1);
	bigU(2,3) = u3(2,1);
	bigU(3,3) = u3(3,1);

	bigU(1,4) = u4(1,1);
	bigU(2,4) = u4(2,1);
	bigU(3,4) = u4(3,1);

// test assembly of bigU
//std::cout << "u1: " << std::endl;
//std::cout << u1(1,1) << " " << u1(2,1) << " " << u1(3,1) << std::endl;
//std::cout << "u2: " << std::endl;
//std::cout << u2(1,1) << " " << u2(2,1) << " " << u2(3,1) << std::endl;
//std::cout << "u3: " << std::endl;
//std::cout << u3(1,1) << " " << u3(2,1) << " " << u3(3,1) << std::endl;
//
//std::cout << "bigU: " << std::endl;
//std::cout << bigU(1,1) << " " << bigU(1,2) << " " << bigU(1,3) << std::endl; 
//std::cout << bigU(2,1) << " " << bigU(2,2) << " " << bigU(2,3) << std::endl; 
//std::cout << bigU(3,1) << " " << bigU(3,2) << " " << bigU(3,3) << std::endl; 
//std::cout << bigU(4,1) << " " << bigU(4,2) << " " << bigU(4,3) << std::endl; 
//

	bigB = getBigB(tempCell);
	matrix result;
	result = 1/(6.0*cell_volume)*bigU*bigB;
	return result;
}

// get volume of cell
REAL assembly::getCellVolume(const cell& tempCell) const{
	int n1,n2,n3,n4;
	REAL x1, y1, z1;	// the current position of node 1
	REAL x2, y2, z2;	// the current position of node 2
	REAL x3, y3, z3;	// the current position of node 3
	REAL x4, y4, z4;	// the current position of node 4
	// get ID for 4 nodes
	n1 = tempCell.getm();
	n2 = tempCell.getn();
	n3 = tempCell.geti();
	n4 = tempCell.getj();
	
	// get coordinates for the 4 nodes
	std::vector<particle*>::const_iterator it;
	
	it = ParticleVec.begin();
	for(int i=0; i!=n1-1; i++)
		it++;	// go to the n1 particle
	x1 = ((*it)->getCurrPosition()).getx();
	y1 = ((*it)->getCurrPosition()).gety();
	z1 = ((*it)->getCurrPosition()).getz();
	
	it = ParticleVec.begin();
	for(int i=0; i!=n2-1; i++)
		it++;	// go to the n2 particle
	x2 = ((*it)->getCurrPosition()).getx();
	y2 = ((*it)->getCurrPosition()).gety();
	z2 = ((*it)->getCurrPosition()).getz();

	it = ParticleVec.begin();
	for(int i=0; i!=n3-1; i++)
		it++;	// go to the n3 particle
	x3 = ((*it)->getCurrPosition()).getx();
	y3 = ((*it)->getCurrPosition()).gety();
	z3 = ((*it)->getCurrPosition()).getz();

	it = ParticleVec.begin();
	for(int i=0; i!=n4-1; i++)
		it++;	// go to the n4 particle
	x4 = ((*it)->getCurrPosition()).getx();
	y4 = ((*it)->getCurrPosition()).gety();
	z4 = ((*it)->getCurrPosition()).getz();

	// get volume of this cell
	REAL cell_volume;
	cell_volume = (x2*y3*z4-x2*y4*z3-x3*y2*z4+x3*y4*z2+x4*y2*z3-x4*y3*z2
		      -x1*y3*z4+x1*y4*z3+x3*y1*z4-x3*y4*z1-x4*y1*z3+x4*y3*z1
		      +x1*y2*z4-x1*y4*z2-x2*y1*z4+x2*y4*z1+x4*y1*z2-x4*y2*z1
		      -x1*y2*z3+x1*y3*z2+x2*y1*z3-x2*y3*z1-x3*y1*z2+x3*y2*z1)/6.0;
	return cell_volume;
}

// calculate spatial velocity gradient tensor for each tet, June 24, 2013
matrix assembly::getdvdx_curr(const cell& tempCell) const{
	int n1,n2,n3,n4;

	// get ID for 4 nodes
	n1 = tempCell.getm();
	n2 = tempCell.getn();
	n3 = tempCell.geti();
	n4 = tempCell.getj();
	// get velocities for 4 nodes
	matrix v1(3,1), v2(3,1), v3(3,1), v4(3,1);
	std::vector<particle*>::const_iterator it;
	// go to particle n1
	it = ParticleVec.begin();
	for(int i=0; i!=n1-1; i++)
		it++;	// go to the n1 particle
	v1(1,1) = ((*it)->getCurrVelocity()).getx();
	v1(2,1) = ((*it)->getCurrVelocity()).gety();
	v1(3,1) = ((*it)->getCurrVelocity()).getz();
	// go to particle n2
	it = ParticleVec.begin();
	for(int i=0; i!=n2-1; i++)
		it++;	// go to the n2 particle
	v2(1,1) = ((*it)->getCurrVelocity()).getx();
	v2(2,1) = ((*it)->getCurrVelocity()).gety();
	v2(3,1) = ((*it)->getCurrVelocity()).getz();
	// go to particle n3
	it = ParticleVec.begin();
	for(int i=0; i!=n3-1; i++)
		it++;	// go to the n3 particle
	v3(1,1) = ((*it)->getCurrVelocity()).getx();
	v3(2,1) = ((*it)->getCurrVelocity()).gety();
	v3(3,1) = ((*it)->getCurrVelocity()).getz();
	// go to particle n4
	it = ParticleVec.begin();
	for(int i=0; i!=n4-1; i++)
		it++;	// go to the n4 particle
	v4(1,1) = ((*it)->getCurrVelocity()).getx();
	v4(2,1) = ((*it)->getCurrVelocity()).gety();
	v4(3,1) = ((*it)->getCurrVelocity()).getz();

	// get volume of this cell
	REAL cell_volume;
	cell_volume = getCellVolume(tempCell);
	if(cell_volume < 0){
		std::cout << "Error: volume is negative in getdudx_curr()!" << std::endl;
	}

	matrix bigV(3,4), bigB(4,3);
	// assemble bigU
	bigV(1,1) = v1(1,1);
	bigV(2,1) = v1(2,1);
	bigV(3,1) = v1(3,1);

	bigV(1,2) = v2(1,1);
	bigV(2,2) = v2(2,1);
	bigV(3,2) = v2(3,1);

	bigV(1,3) = v3(1,1);
	bigV(2,3) = v3(2,1);
	bigV(3,3) = v3(3,1);

	bigV(1,4) = v4(1,1);
	bigV(2,4) = v4(2,1);
	bigV(3,4) = v4(3,1);

	bigB = getBigB(tempCell);
	matrix result;
	result = 1/(6.0*cell_volume)*bigV*bigB;
	return result;
}


// caclulate finite granular strain
matrix assembly::getFiniteStrain(){
	matrix finite_strain(3,3);
	matrix dudx_trans(3,3);
	REAL init_totalVolume = 0;	// the summed initial volume of all cells
	// initialize finite_strain(3,3)
	for(int ir=1; ir!=4;ir++){
		for(int ic=1; ic!=4; ic++){
			finite_strain(ir,ic) = 0;
			dudx_trans(ir,ic) = 0;
		}
	}
	average_dudx_Lagrangian = finite_strain;	// used to test quadratic terms, April 22, 2013
	matrix temp_dudx(3,3);
	REAL temp_volume;
	for(std::vector<cell*>::const_iterator iter=cellVec_init.begin(); iter!=cellVec_init.end(); iter++){
		temp_volume = (*iter)->getInitialCellVolume();
		temp_dudx = getdudx(**iter);
//		finite_strain += temp_volume*(temp_dudx+temp_dudx.getTrans()+temp_dudx.getTrans()*temp_dudx); 
		average_dudx_Lagrangian += temp_volume*temp_dudx;	// used to test quadratic terms
		init_totalVolume += temp_volume;
	}

//	finite_strain = finite_strain/(2.0*initVolume);
	average_dudx_Lagrangian = average_dudx_Lagrangian/init_totalVolume;	// used to test quadratic terms
	dudx_trans = average_dudx_Lagrangian.getTrans();
	finite_strain = (average_dudx_Lagrangian+dudx_trans+dudx_trans*average_dudx_Lagrangian)/2;
	return finite_strain;
}

// calculate 1/(2V)*sum(dudx'*dudx*vL), for mixed finite strain
matrix assembly::getHigherStrain() const{
	matrix higher_strain(3,3);
	// initialize finite_strain(3,3)
	for(int ir=1; ir!=4;ir++){
		for(int ic=1; ic!=4; ic++){
			higher_strain(ir,ic) = 0;
		}
	}
	matrix temp_dudx(3,3);
	REAL temp_volume;
	for(std::vector<cell*>::const_iterator iter=cellVec_init.begin(); iter!=cellVec_init.end(); iter++){
		temp_volume = (*iter)->getInitialCellVolume();
		temp_dudx = getdudx(**iter);
		higher_strain += temp_volume*(temp_dudx.getTrans()*temp_dudx); 
	}
	higher_strain = higher_strain/(2.0*initVolume);
	return higher_strain;
}

// caclulate finite granular strain
matrix assembly::getEulerianStrain(){
	matrix eulerian_strain(3,3);
	matrix dudx_trans(3,3);
	REAL curr_totalVolume = 0;	// the current summed volume of all cells
	// initialize finite_strain(3,3)
	for(int ir=1; ir!=4;ir++){
		for(int ic=1; ic!=4; ic++){
			eulerian_strain(ir,ic) = 0;
			dudx_trans(ir,ic) = 0;
		}
	}
	average_dudx_Eulerian = eulerian_strain;	// used to test quadratic terms
	matrix temp_dudx(3,3);
	REAL temp_volume;
	for(std::vector<cell*>::const_iterator iter=cellVec.begin(); iter!=cellVec.end(); iter++){
		temp_volume = getCellVolume(**iter);
		temp_dudx = getdudx_curr(**iter);

//		eulerian_strain += temp_volume*(temp_dudx+temp_dudx.getTrans()-temp_dudx.getTrans()*temp_dudx); 
		average_dudx_Eulerian += temp_volume*temp_dudx;	// used to test quadratic terms
		curr_totalVolume += temp_volume;
	}
//	eulerian_strain = eulerian_strain/(2.0*Volume);
	average_dudx_Eulerian = average_dudx_Eulerian/curr_totalVolume;	// used to test quadratic terms
	dudx_trans = average_dudx_Eulerian.getTrans();
	eulerian_strain = (average_dudx_Eulerian+dudx_trans-dudx_trans*average_dudx_Eulerian)/2;
	return eulerian_strain;
}

REAL assembly::getCurrCellVolume() const{

	REAL curr_totalVolume = 0;
	REAL temp_volume;
	for(std::vector<cell*>::const_iterator iter=cellVec.begin(); iter!=cellVec.end(); iter++){
		temp_volume = getCellVolume(**iter);
		curr_totalVolume += temp_volume;
	}

	return curr_totalVolume;

} // getCurrCellVolume

// caclulate average spatial velocity gradient tensor, June 24, 2013
matrix assembly::getAverage_dvdx() const{
	matrix average_dvdx(3,3);
	REAL curr_totalVolume = 0;	// the current summed volume of all cells
	// initialize average_dvdx(3,3)
	for(int ir=1; ir!=4;ir++){
		for(int ic=1; ic!=4; ic++){
			average_dvdx(ir,ic) = 0;
		}
	}
	matrix temp_dvdx(3,3);
	REAL temp_volume;
	for(std::vector<cell*>::const_iterator iter=cellVec.begin(); iter!=cellVec.end(); iter++){
		temp_volume = getCellVolume(**iter);
		temp_dvdx = getdvdx_curr(**iter);

		average_dvdx += temp_volume*temp_dvdx;	
		curr_totalVolume += temp_volume;
	}
	average_dvdx = average_dvdx/curr_totalVolume;	
	return average_dvdx;
}


// caclulate finite granular strain
matrix assembly::getEuler_HOT() const{
	matrix eulerian_strain(3,3);
	// initialize finite_strain(3,3)
	for(int ir=1; ir!=4;ir++){
		for(int ic=1; ic!=4; ic++){
			eulerian_strain(ir,ic) = 0;
		}
	}
	matrix temp_dudx(3,3);
	REAL temp_volume;
	for(std::vector<cell*>::const_iterator iter=cellVec.begin(); iter!=cellVec.end(); iter++){
		temp_volume = getCellVolume(**iter);
		temp_dudx = getdudx_curr(**iter);

		eulerian_strain += temp_volume*(-1*temp_dudx.getTrans()*temp_dudx); 
	}
	eulerian_strain = eulerian_strain/(2.0*Volume);
	return eulerian_strain;
}

// above is for finite granular strain calculation



vec assembly::getShearForce(int bdry) const{
    std::vector<RGDBDRY*>::const_iterator it;
    for(it=RBVec.begin();it!=RBVec.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getShearForce();
    }
    return 0;
}


REAL assembly::getAvgNormal(int bdry) const{
    std::vector<RGDBDRY*>::const_iterator it;
    for(it=RBVec.begin();it!=RBVec.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getAvgNormal();
    }
    return 0;
}


vec assembly::getApt(int bdry) const{
    std::vector<RGDBDRY*>::const_iterator it;
    for(it=RBVec.begin();it!=RBVec.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getApt();
    }
    return 0;
}


vec assembly::getDirc(int bdry) const{
    std::vector<RGDBDRY*>::const_iterator it;
    for(it=RBVec.begin();it!=RBVec.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getDirc();
    }
    return 0;
}


REAL assembly::getArea(int n) const{
    std::vector<RGDBDRY*>::const_iterator it;
    for(it=RBVec.begin();it!=RBVec.end();++it){
	if((*it)->getBdryID()==n)
	    return (*it)->area;
    }
    return 0;
}


void assembly::setArea(int n, REAL a){
    std::vector<RGDBDRY*>::iterator it;
    for(it=RBVec.begin();it!=RBVec.end();++it){
	if((*it)->getBdryID()==n)
	    (*it)->area=a;
    }
}


REAL assembly::getAverageRigidPressure() const{
    std::vector<RGDBDRY*>::const_iterator rt;
    REAL avgpres=0;
    for(rt=RBVec.begin();rt!=RBVec.end();++rt)
	avgpres+=vfabs((*rt)->getNormalForce())/(*rt)->getArea();
    return avgpres/=RgdBdryNum;
}


// only update CoefOfLimits[0] for specified boundaries
void assembly::updateRB(int bn[], UPDATECTL rbctl[], int num){
    for(int i=0;i<num;i++){
	for(std::vector<RGDBDRY*>::iterator rt=RBVec.begin();rt!=RBVec.end();++rt){
	    if((*rt)->getBdryID()==bn[i]){
		(*rt)->update(rbctl[i]);
		break;
	    }
	}
    }
}


// update CoefOfLimits[1,2,3,4] for all 6 boundaries
void assembly::updateRB6(){
    for(std::vector<RGDBDRY*>::iterator rt=RBVec.begin();rt!=RBVec.end();++rt){
	if((*rt)->getBdryID()==1 || (*rt)->getBdryID()==3){
	    for(std::vector<RGDBDRY*>::iterator lt=RBVec.begin();lt!=RBVec.end();++lt){
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
	    for(std::vector<RGDBDRY*>::iterator lt=RBVec.begin();lt!=RBVec.end();++lt){
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
	    for(std::vector<RGDBDRY*>::iterator lt=RBVec.begin();lt!=RBVec.end();++lt){
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


// upgrade CoefOfLimits[1,2,3,4] for rectangular pile
void assembly::updateRectPile(){
    for(std::vector<RGDBDRY*>::iterator rt=RBVec.begin();rt!=RBVec.end();++rt){
	if((*rt)->getBdryID()==7 || (*rt)->getBdryID()==9 ){
	    for(std::vector<RGDBDRY*>::iterator lt=RBVec.begin();lt!=RBVec.end();++lt){
		if((*lt)->getBdryID()==10)
		    (*rt)->CoefOfLimits[1].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==8)
		    (*rt)->CoefOfLimits[2].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==11)
		    (*rt)->CoefOfLimits[3].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==12)
		    (*rt)->CoefOfLimits[4].apt=(*lt)->CoefOfLimits[0].apt;
	    }
	}
	else if((*rt)->getBdryID()==8 || (*rt)->getBdryID()==10){
	    for(std::vector<RGDBDRY*>::iterator lt=RBVec.begin();lt!=RBVec.end();++lt){
		if((*lt)->getBdryID()==7)
		    (*rt)->CoefOfLimits[1].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==9)
		    (*rt)->CoefOfLimits[2].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==11)
		    (*rt)->CoefOfLimits[3].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==12)
		    (*rt)->CoefOfLimits[4].apt=(*lt)->CoefOfLimits[0].apt;
	    }
	}
    }
}


void assembly::updateFB(int bn[], UPDATECTL fbctl[], int num){
    std::vector<FLBBDRY*>::iterator ft;
    int i,k=1;
    for(i=0;i<num;i++){
	for(ft=FBVec.begin();ft!=FBVec.end();++ft){
	    if(k++==bn[i]){
		UPDATECTL ctl[2];
		ctl[0]=fbctl[bn[i]*2-2];
		ctl[1]=fbctl[bn[i]*2-1];
		(*ft)->update(ctl,2);
		break;
	    }
	}
    }
}


// create a specimen from discreate particles through floating and then gravitation,
// file cre_particle contains the final particle information,
// file cre_boundary contains the final boundary information.
void assembly::deposit_RgdBdry(int   freetype,
			       int   total_steps,  
			       int   snapshots,
			       int   interval,
			       REAL  rFloHeight,
//			       REAL sampleVolumeRatio,   // Modified by Boning Zhang
			       const char* iniptclfile,   
			       const char* inibdryfile,
			       const char* particlefile, 
			       const char* contactfile,
			       const char* progressfile, 
			       const char* trmparticle,
			       const char* trmboundary,
			       const char* debugfile)
{
  generate(iniptclfile, freetype, rFloHeight); 
//  generate(iniptclfile, freetype, sampleVolumeRatio);   // Modified by Boning Zhang
  
  buildBoundary(5, inibdryfile); // container unchanged
  
  deposit(total_steps,        // total_steps
	  snapshots,          // number of snapshots
	  interval,           // print interval
	  iniptclfile,        // input file, initial particles
	  inibdryfile,        // input file, initial boundaries
	  particlefile,       // output file, resulted particles, including snapshots 
	  contactfile,        // output file, resulted contacts, including snapshots 
	  progressfile,       // output file, statistical info
	  debugfile);         // output file, debug info
  
  buildBoundary(trmboundary);   // output file, containing boundaries info
  
  trim(false,                 // recreate from input file or not
       "dep_particle_end",    // input file, particles to be trimmed
       trmparticle);          // output file, trimmed particles
}


// read the iniptclfile, then reset the particles' orientations by random 
// rotation method from Professor Brannon. July 27, 2014
void assembly::deposit_Compare(int   freetype,
			       int   total_steps,  
			       int   snapshots,
			       int   interval,
			       REAL  rFloHeight,
			       const char* iniptclfile,   
			       const char* inibdryfile,
			       const char* particlefile, 
			       const char* contactfile,
			       const char* progressfile, 
			       const char* trmparticle,
			       const char* trmboundary,
			       const char* debugfile)
{
  readSample(iniptclfile);

  resetPctlOrientations();

  const char* iniptclfileReset = "ini_particle_reset";

  printParticle(iniptclfileReset);

  buildBoundary(5, inibdryfile); // container unchanged
  
  deposit(total_steps,        // total_steps
	  snapshots,          // number of snapshots
	  interval,           // print interval
	  iniptclfileReset,   // input file, initial particles
	  inibdryfile,        // input file, initial boundaries
	  particlefile,       // output file, resulted particles, including snapshots 
	  contactfile,        // output file, resulted contacts, including snapshots 
	  progressfile,       // output file, statistical info
	  debugfile);         // output file, debug info
  
  buildBoundary(trmboundary);   // output file, containing boundaries info
  
  trim(false,                 // recreate from input file or not
       "dep_particle_end",    // input file, particles to be trimmed
       trmparticle);          // output file, trimmed particles
}


// this random orientation method is from Professor Brannon, July 27, 2014
// http://csmbrannon.net/2013/11/25/visualization-of-rotation/
void assembly::resetPctlOrientations(){

    for(std::vector<particle*>::iterator ip = ParticleVec.begin(); ip!=ParticleVec.end(); ip++){

    	REAL r1 = ran(&idum);
    	REAL r2 = ran(&idum);
    	REAL r3 = ran(&idum);

    	REAL theta = 2*PI*r1;
    	REAL phi = 2*PI*r2;
    	REAL rho = sqrt(r3);

    	matrix n(3,1);
    	n(1,1) = rho*cos(phi);
    	n(2,1) = rho*sin(phi);
    	n(3,1) = sqrt(1.0-rho*rho);

    	matrix identity3x3(3,3);
    	identity3x3(1,1) = 1; identity3x3(2,2) = 1; identity3x3(3,3) = 1;
	
    	matrix thetaMatrix(3,3);
    	thetaMatrix(1,1) = cos(theta); thetaMatrix(1,2) = sin(theta); thetaMatrix(1,3) = 0;
   	thetaMatrix(2,1) =-sin(theta); thetaMatrix(2,2) = cos(theta); thetaMatrix(2,3) = 0;
    	thetaMatrix(3,1) = 0;	   thetaMatrix(3,2) = 0;	  thetaMatrix(3,3) = 1;

    	matrix R = (2*n*n.getTrans()-identity3x3)*thetaMatrix;

//std::cout << std::endl;
//std::cout << "R.getTrans(): " << R.getTrans().print() << std::endl;
//std::cout << "R.getInvs(): " << R.getInvs().print() << std::endl;
//std::cout << "det(R): " << R(1,1)*R(2,2)*R(3,3)-R(1,1)*R(2,3)*R(3,2)-R(1,2)*R(2,1)*R(3,3)+R(1,2)*R(2,3)*R(3,1)+R(1,3)*R(2,1)*R(3,2)-R(1,3)*R(2,2)*R(3,1) << std::endl;

	// assume a set of principle directions are along the positive x, y and z directions
	// then rotate it by a random rotation matrix so that we will have random orientations for ellipsoids
	matrix x_mat(1,3);
	x_mat(1,1) = 1; x_mat(1,2) = 0; x_mat(1,3) = 0;	
	matrix y_mat(1,3);
	y_mat(1,1) = 0; y_mat(1,2) = 1; y_mat(1,3) = 0;
	matrix z_mat(1,3);
	z_mat(1,1) = 0; z_mat(1,2) = 0; z_mat(1,3) = 1;

	x_mat = x_mat*R;
	y_mat = y_mat*R;
	z_mat = z_mat*R;

	(*ip)->setCurrDirecA(vec( acos(x_mat(1,1)), acos(x_mat(1,2)), acos(x_mat(1,3))) );
	(*ip)->setCurrDirecB(vec( acos(y_mat(1,1)), acos(y_mat(1,2)), acos(y_mat(1,3))) );
	(*ip)->setCurrDirecC(vec( acos(z_mat(1,1)), acos(z_mat(1,2)), acos(z_mat(1,3))) );


    }
}



void assembly::deposit_repose(int   interval,
			      const char* inibdryfile,
			      const char* particlefile, 
			      const char* contactfile,
			      const char* progressfile, 
			      const char* debugfile)
{
  this->container = container;
  
  buildBoundary(5, inibdryfile); // container unchanged
  
  angleOfRepose(interval,           // print interval
		inibdryfile,        // input file, initial boundaries
		particlefile,       // output file, resulted particles, including snapshots 
		contactfile,        // output file, resulted contacts, including snapshots 
		progressfile,       // output file, statistical info
		debugfile);         // output file, debug info
}

void assembly::angleOfRepose(int   interval,
			     const char* inibdryfile,
			     const char* particlefile, 
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
  
  g_debuginf.open(debugfile);
  if(!g_debuginf) { cout << "stream error!" << endl; exit(-1); }
  g_debuginf.setf(std::ios::scientific, std::ios::floatfield);
  
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

  REAL maxRadius = gradInfo.getMaxPtclRadius();
  REAL maxDiameter = maxRadius * 2.0;
  REAL z0 = container.getMinCorner().getz();
  vector<particle*> lastPtcls;
  particle *newPtcl = NULL;
  int layers = 1; // how many layers of new particles to generate each time?

  g_iteration = 0; 
  TotalNum = 0;
  REAL zCurr;
  gettimeofday(&time_w1,NULL);
  // iterations starting ...
  do
    {
      // 1. add particle
      if ( TotalNum == 0 ) {
	zCurr = z0 + maxRadius;

	for ( int i = 0; i != layers; ++i) {
	  newPtcl = new particle(TotalNum+1, 0, vec(0, 0, zCurr + maxDiameter * i), gradInfo, YOUNG, POISSON);
	  ParticleVec.push_back(newPtcl);
	  ++TotalNum;
	  lastPtcls.push_back(newPtcl);
	  
	  newPtcl = new particle(TotalNum+1, 0, vec(maxDiameter, 0, zCurr + maxDiameter * i), gradInfo, YOUNG, POISSON);
	  ParticleVec.push_back(newPtcl);
	  ++TotalNum;
	  //lastPtcls.push_back(newPtcl);
	  
	  newPtcl = new particle(TotalNum+1, 0, vec(-maxDiameter, 0, zCurr + maxDiameter * i), gradInfo, YOUNG, POISSON);
	  ParticleVec.push_back(newPtcl);
	  ++TotalNum;
	  //lastPtcls.push_back(newPtcl);
	  
	  newPtcl = new particle(TotalNum+1, 0, vec(0, maxDiameter, zCurr + maxDiameter * i), gradInfo, YOUNG, POISSON);
	  ParticleVec.push_back(newPtcl);
	  ++TotalNum;
	  //lastPtcls.push_back(newPtcl);
	  
	  newPtcl = new particle(TotalNum+1, 0, vec(0, -maxDiameter, zCurr + maxDiameter * i), gradInfo, YOUNG, POISSON);
	  ParticleVec.push_back(newPtcl);
	  ++TotalNum;
	  //lastPtcls.push_back(newPtcl);
	}
	toSnapshot = true;

      }
      else {
	vector<particle*>::iterator it;
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

	  lastPtcls.clear(); // do not delete those pointers to release memory; ParticleVec will do it.
	  zCurr = getMaxCenterHeight() + maxDiameter;

	  for ( int i = 0; i != layers; ++i) {
	    newPtcl = new particle(TotalNum+1, 0, vec(0, 0, zCurr + maxDiameter * i), gradInfo, YOUNG, POISSON);
	    ParticleVec.push_back(newPtcl);
	    ++TotalNum;
	    lastPtcls.push_back(newPtcl);
	    
	    newPtcl = new particle(TotalNum+1, 0, vec(maxDiameter, 0, zCurr + maxDiameter * i), gradInfo, YOUNG, POISSON);
	    ParticleVec.push_back(newPtcl);
	    ++TotalNum;
	    //lastPtcls.push_back(newPtcl);
	    
	    newPtcl = new particle(TotalNum+1, 0, vec(-maxDiameter, 0, zCurr + maxDiameter * i), gradInfo, YOUNG, POISSON);
	    ParticleVec.push_back(newPtcl);
	    ++TotalNum;
	    //lastPtcls.push_back(newPtcl);
	    
	    newPtcl = new particle(TotalNum+1, 0, vec(0, maxDiameter, zCurr + maxDiameter * i), gradInfo, YOUNG, POISSON);
	    ParticleVec.push_back(newPtcl);
	    ++TotalNum;
	    //lastPtcls.push_back(newPtcl);
	    
	    newPtcl = new particle(TotalNum+1, 0, vec(0, -maxDiameter, zCurr + maxDiameter * i), gradInfo, YOUNG, POISSON);
	    ParticleVec.push_back(newPtcl);
	    ++TotalNum;
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
      rigidBoundaryForce(bdry_penetr, bdry_cntnum);
      
      // 6. update particles' velocity/omga/position/orientation based on force/moment.
      updateParticle();
      
      // 7. (1) output particles and contacts information as snapshots.
      if (toSnapshot){
	sprintf(stepsstr, "%03d", stepsnum); 
	strcpy(stepsfp, particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
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
      if (g_iteration % interval == 0) {
	gettimeofday(&time_w2,NULL);
	REAL t1=getTransEnergy();
	REAL t2=getRotatEnergy();
	REAL t3=getPotenEnergy(-0.025);
	progressinf << setw(OWID) << g_iteration
		    << setw(OWID) << getPossCntctNum()
		    << setw(OWID) << getActualCntctNum()
		    << setw(OWID) << getAveragePenetration()
		    << setw(OWID) << avgNormal
		    << setw(OWID) << avgTangt
		    << setw(OWID) << getAverageVelocity() 
		    << setw(OWID) << getAverageOmga()
		    << setw(OWID) << getAverageForce()   
		    << setw(OWID) << getAverageMoment()
		    << setw(OWID) << t1
		    << setw(OWID) << t2
		    << setw(OWID) << (t1+t2)
		    << setw(OWID) << t3
		    << setw(OWID) << (t1+t2+t3)
		    << setw(OWID) << void_ratio
		    << setw(OWID) << void_ratio/(1+void_ratio)
		    << setw(OWID) << 2.0*(getActualCntctNum()
				      +bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
				      +bdry_cntnum[4]+bdry_cntnum[6])/TotalNum
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
      ++g_iteration;
      
    } while (TotalNum < 2000); //( zCurr < container.getMaxCorner().getz() );  //(++g_iteration < total_steps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    g_timeinf << setw(4) << "end" << " " << ctime(&timeStamp) << flush;

    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


// freetype:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
// ht- how many times of size would be the floating height

void assembly::generate(const char* particlefile,
			int freetype,
			REAL rFloHeight)
{
  REAL x,y,z;
  particle* newptcl;
  TotalNum = 0;
  REAL dimx     = container.getDimx();
  REAL dimy     = container.getDimy();
  REAL dimz     = container.getDimz();
  REAL diameter = gradInfo.getMaxPtclRadius()*2.0;

  REAL offset   = 0;
  REAL edge     = diameter;
  if (gradInfo.getSize().size() == 1 &&
      gradInfo.getPtclRatioBA() == 1.0 && 
      gradInfo.getPtclRatioCA() == 1.0) {
    edge   = diameter*2.0;
    offset = diameter*0.25;
  }
  
  REAL x1 = container.getMinCorner().getx() + edge;
  REAL y1 = container.getMinCorner().gety() + edge;
  REAL x2 = container.getMaxCorner().getx() - edge;
  REAL y2 = container.getMaxCorner().gety() - edge;
  REAL z1 = container.getMinCorner().getz() + diameter;

  REAL x0 = container.getCenter().getx();
  REAL y0 = container.getCenter().gety();
  REAL z0 = container.getCenter().getz();

  if (freetype == 0) {      // just one free particle
    newptcl = new particle(TotalNum+1, 0, vec(x0,y0,z0), gradInfo, YOUNG, POISSON);
    ParticleVec.push_back(newptcl);
    TotalNum++;
  }
  else if (freetype == 1) { // a horizontal layer of free particles
    z = container.getMaxCorner().getz();
    for (x = x1; x - x2 < EPS; x += diameter)
      for (y = y1; y - y2 < EPS; y += diameter) {
	newptcl = new particle(TotalNum+1, 0, vec(x,y,z), gradInfo, YOUNG, POISSON);
	ParticleVec.push_back(newptcl);
	TotalNum++;
      }
  }
  else if (freetype == 2) { // multiple layers of free particles
    for (z = z1; z < z1 + dimz*rFloHeight; z += diameter) {
      for (x = x1 + offset; x - x2 < EPS; x += diameter)
// this change is for generating a 2D sample
//	for (y = y1 + offset; y - y2 < EPS; y += diameter) {
	  newptcl = new particle(TotalNum+1, 0, vec(x,y0,z), gradInfo, YOUNG, POISSON);
	  ParticleVec.push_back(newptcl);
	  TotalNum++;
//	}	
      offset *= -1;
    }
  }
  
  printParticle(particlefile);
  
}

/*
// The above comment out is made by Boning Zhang, so as below new generate function. the purpose of this new function is to generate 
//enough particles to fill up our container.

void assembly::generate(const char* particlefile,
			int freetype,
			REAL sampleVolumeRatio)
{
  REAL x,y,z;
  REAL rFloHeight;
  particle* newptcl;
  TotalNum = 0;
  REAL targetTotalNum;
  REAL volume100Ptcl;
  int num10Height;
  int iptcl;
  REAL dimx     = container.getDimx();
  REAL dimy     = container.getDimy();
  REAL dimz     = container.getDimz();
  REAL containerVolume = container.getVolume();
  REAL diameter = gradInfo.getMaxPtclRadius()*2.0;

  REAL offset   = 0;
  REAL edge     = diameter;
  if (gradInfo.getSize().size() == 1 &&
      gradInfo.getPtclRatioBA() == 1.0 && 
      gradInfo.getPtclRatioCA() == 1.0) {
    edge   = diameter*2.0;
    offset = diameter*0.25;
  }
  
  REAL x1 = container.getMinCorner().getx() + edge;
  REAL y1 = container.getMinCorner().gety() + edge;
  REAL x2 = container.getMaxCorner().getx() - edge;
  REAL y2 = container.getMaxCorner().gety() - edge;
  REAL z1 = container.getMinCorner().getz() + diameter;

  REAL x0 = container.getCenter().getx();
  REAL y0 = container.getCenter().gety();
  REAL z0 = container.getCenter().getz();

  // get the number of particles if relative height is 10
  num10Height = 0;
  for (z = z1; z < z1 + dimz*10; z += diameter) {
      for (x = x1 + offset; x - x2 < EPS; x += diameter)
	for (y = y1 + offset; y - y2 < EPS; y += diameter) {
	  num10Height++;
	}	
  //    offset *= -1;
    }

  // get the total volume of 100 particles that are generated based on the gradinfo
  volume100Ptcl = 0;
  for (iptcl = 1; iptcl <= 100; iptcl++) {
	dem::particle newptcl(iptcl, 0, vec(0,0,0), gradInfo, YOUNG, POISSON);
	volume100Ptcl += newptcl.getVolume();
//	delete newptcl;
  }

  targetTotalNum = sampleVolumeRatio*containerVolume/volume100Ptcl*100*1.2;  // here 1.2 is safty factor
  rFloHeight	   = targetTotalNum/num10Height*10;
  
  if (freetype == 0) {      // just one free particle
    newptcl = new particle(TotalNum+1, 0, vec(x0,y0,z0), gradInfo, YOUNG, POISSON);
    ParticleVec.push_back(newptcl);
    TotalNum++;
  }
  else if (freetype == 1) { // a horizontal layer of free particles
    z = container.getMaxCorner().getz();
    for (x = x1; x - x2 < EPS; x += diameter)
      for (y = y1; y - y2 < EPS; y += diameter) {
	newptcl = new particle(TotalNum+1, 0, vec(x,y,z), gradInfo, YOUNG, POISSON);
	ParticleVec.push_back(newptcl);
	TotalNum++;
      }
  }
  else if (freetype == 2) { // multiple layers of free particles
    for (z = z1; z < z1 + dimz*rFloHeight; z += diameter) {
      for (x = x1 + offset; x - x2 < EPS; x += diameter)
	for (y = y1 + offset; y - y2 < EPS; y += diameter) {
	  newptcl = new particle(TotalNum+1, 0, vec(x,y,z), gradInfo, YOUNG, POISSON);
	  ParticleVec.push_back(newptcl);
	  TotalNum++;
	}	
      offset *= -1;
    }
  }
  
  printParticle(particlefile);
  
}
*/


/*
// freetype:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
// ht- how many times of size would be the floating height
void assembly::generate(gradation& grad,
			const char* particlefile,
			int freetype,
			REAL ht)
{
    REAL x,y,z;
    particle* newptcl;
    TotalNum = 0;
    REAL est =1.02;
    int grid=9;  
    // grid: dimension of free particle array.
    // 7 - small dimn container
    // 9 - medium dimn container 
    // 11- large dimn container 

    REAL dimn=grad.dimn;
    if (freetype == 0) {      // just one free particle
	newptcl = new particle(TotalNum+1, 0, vec(dimn/2/40,dimn/2/20,dimn/2), grad, YOUNG, POISSON);
	ParticleVec.push_back(newptcl);
	TotalNum++;
    }
    else if (freetype == 1) { // a horizontal layer of free particles
	z=dimn/2;
	for (x=-dimn/2*(grid-1)/10; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
	    for (y=-dimn/2*(grid-1)/10; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5){
		newptcl = new particle(TotalNum+1, 0, vec(x,y,z), grad, YOUNG, POISSON);
		ParticleVec.push_back(newptcl);
		TotalNum++;
	    }
    }
    else if (freetype == 2) { // multiple layers of free particles
	REAL offset=0; // 0 for ellipsoids; dimn/2/5/5 for spheres
	if (grad.ratio_ba==1.0 && grad.ratio_ca==1.0)
	    offset = dimn/2/5/5;
	REAL z0 = -dimn/2*9/10 ;// dimn/2;
	for (z=z0; z<z0 + dimn*ht; z+=dimn/2/5) {
	//for (z=-dimn/2*4/5; z<dimn/2 + dimn*ht; z+=dimn/2/10) { // spheres
	    for (x=-dimn/2*(grid-1)/10+offset; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
		for (y=-dimn/2*(grid-1)/10+offset; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5){
		    newptcl = new particle(TotalNum+1, 0, vec(x,y,z), grad, YOUNG, POISSON);
		    ParticleVec.push_back(newptcl);
		    TotalNum++;
		}	
	    offset *= -1;
	}
    }

    printParticle(particlefile);

}
*/
 /*
// create a specimen from discreate particles through floating and then gravitation,
// boundaries are composed of fixed particles.
void assembly::deposit_PtclBdry(gradation& grad,
				int   freetype,
				REAL rsize,
				int   total_steps,  
				int   snapshots,
				int   interval,
				const char* iniptclfile,   
				const char* particlefile, 
				const char* contactfile,
				const char* progressfile, 
				const char* debugfile)
{
    if (grad.rorc == 1) {
	RORC = grad.rorc;
	container.setCenter(vec(0,0,0));
	container.setDimx(grad.dimn);
	container.setDimy(grad.dimn);
	container.setDimz(grad.dimn);
	
	generate_p(grad, iniptclfile, freetype, rsize, 4.0);
	deposit_p(total_steps,        // total_steps
		  snapshots,          // number of snapshots
		  interval,           // print interval
		  grad.dimn,          // dimension of particle-composed-boundary
		  rsize,              // relative container size
		  iniptclfile,        // input file, initial particles
		  particlefile,       // output file, resulted particles, including snapshots 
		  contactfile,        // output file, resulted contacts, including snapshots 
		  progressfile,       // output file, statistical info
		  debugfile);         // output file, debug info
    }
}
 */
  /*
// freetype:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
// ht- how many times of size would be the floating height
void assembly::generate_p(gradation&  grad,
			 const char* particlefile,
			 int freetype,
			 REAL rsize,
			 REAL ht)
{
    REAL x,y,z;
    particle* newptcl;
    TotalNum = 0;
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
	    newptcl = new particle(TotalNum+1, 1, vec(x,y,z), grad.ptclsize[0]*0.99, YOUNG, POISSON);
	    ParticleVec.push_back(newptcl);
	    TotalNum++;
	}

    // particle boundary 2
    y=dimn/2*(grid+1)/10;
    for (x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5)
	for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5){
	    newptcl = new particle(TotalNum+1, 1, vec(x,y,z), grad.ptclsize[0]*0.99, YOUNG, POISSON);
	    ParticleVec.push_back(newptcl);
	    TotalNum++;
	}

    // particle boundary 3
    x=-dimn/2*(grid+1)/10;
    for (y=-dimn/2*grid/10; y<dimn/2*grid/10*est; y+=dimn/2/5)
	for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5){
	    newptcl = new particle(TotalNum+1, 1, vec(x,y,z), grad.ptclsize[0]*0.99, YOUNG, POISSON);
	    ParticleVec.push_back(newptcl);
	    TotalNum++;
	}

    // particle boundary 4
    y=-dimn/2*(grid+1)/10;
    for (x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5)
	for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5){
	    newptcl = new particle(TotalNum+1, 1, vec(x,y,z), grad.ptclsize[0]*0.99, YOUNG, POISSON);
	    ParticleVec.push_back(newptcl);
	    TotalNum++;
	}

    // particle boundary 6
    z=-dimn/2;
    for (y=-dimn/2*grid/10; y<dimn/2*grid/10*est; y+=dimn/2/5)
	for( x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5){
	    newptcl = new particle(TotalNum+1, 1, vec(x,y,z), grad.ptclsize[0]*0.99, YOUNG, POISSON);
	    ParticleVec.push_back(newptcl);
	    TotalNum++;
	}

    if (freetype == 0) {      // just one free particle
	newptcl = new particle(TotalNum+1, 0, vec(dimn/2/40,dimn/2/20,dimn/2), grad, YOUNG, POISSON);
	ParticleVec.push_back(newptcl);
	TotalNum++;
    }
    else if (freetype == 1) { // a horizontal layer of free particles
	z=dimn/2;
	for (x=-dimn/2*(grid-1)/10; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
	    for (y=-dimn/2*(grid-1)/10; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5){
		newptcl = new particle(TotalNum+1, 0, vec(x,y,z), grad, YOUNG, POISSON);
		ParticleVec.push_back(newptcl);
		TotalNum++;
	    }
    }
    else if (freetype == 2) { // multiple layers of free particles
	for (z=dimn/2; z<dimn/2 + dimn*ht; z+=dimn/2/5)
	    for (x=-dimn/2*(grid-1)/10; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
		for (y=-dimn/2*(grid-1)/10; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5){
		    newptcl = new particle(TotalNum+1, 0, vec(x,y,z), grad, YOUNG, POISSON);
		    ParticleVec.push_back(newptcl);
		    TotalNum++;
		}	
    }
    
    printParticle(particlefile);
    
}
*/

void assembly::scale_PtclBdry(int   total_steps,  
			      int   snapshots,
			      int   interval,
			      REAL dimn,
			      REAL rsize,
			      const char* iniptclfile,   
			      const char* particlefile, 
			      const char* contactfile,
			      const char* progressfile, 
			      const char* debugfile)
{
    deposit_p(total_steps,        // total_steps
	      snapshots,          // number of snapshots
	      interval,           // print interval
	      dimn,               // dimension of particle-composed-boundary
	      rsize,              // relative container size
	      iniptclfile,        // input file, initial particles
	      particlefile,       // output file, resulted particles, including snapshots 
	      contactfile,        // output file, resulted contacts, including snapshots 
	      progressfile,       // output file, statistical info
	      debugfile);         // output file, debug info
}


// collapse a deposited specimen through gravitation
void assembly::collapse(int   total_steps,  
			int   snapshots,
			int   interval,
			const char* iniptclfile,
			const char* initboundary,
			const char* particlefile,
			const char* contactfile,
			const char* progressfile,
			const char* debugfile)
{
  buildBoundary(1,              // 1-only bottom boundary; 5-no top boundary;6-boxed 6 boundaries
	      initboundary);  // output file, containing boundaries info
  
  deposit(total_steps,        // number of iterations
	  snapshots,          // number of snapshots
	  interval,           // print interval
	  iniptclfile,        // input file, initial particles
	  initboundary,       // input file, boundaries
	  particlefile,       // output file, resulted particles, including snapshots 
	  contactfile,        // output file, resulted contacts, including snapshots 
	  progressfile,       // output file, statistical info
	  debugfile);         // output file, debug info
}

  
void assembly::buildBoundary(int bdrynum,
			   const char* boundaryfile)
{
  std::ofstream ofs(boundaryfile);
  if(!ofs) { cout << "stream error!" << endl; exit(-1);}
  REAL x1,x2,y1,y2,z1,z2,x0,y0,z0;
  x1 = container.getMinCorner().getx();
  y1 = container.getMinCorner().gety();
  z1 = container.getMinCorner().getz();
  x2 = container.getMaxCorner().getx();
  y2 = container.getMaxCorner().gety();
  z2 = container.getMaxCorner().getz();
  x0 = container.getCenter().getx();
  y0 = container.getCenter().gety();
  z0 = container.getCenter().getz();

  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs << setw(OWID) << 0
      << setw(OWID) << bdrynum << endl << endl;
  
  if (bdrynum == 1){   // only a bottom boundary
    ofs << setw(OWID) << 1 << endl
        << setw(OWID) << 6
        << setw(OWID) << 1
        << setw(OWID) << 0 << endl
      
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
  else if (bdrynum == 5){ // no top boundary
    // boundary 1
    ofs << setw(OWID) << 1 << endl
        << setw(OWID) << 1
        << setw(OWID) << 4
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
        << setw(OWID) << 4
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
        << setw(OWID) << 4
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
        << setw(OWID) << 4
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
        << setw(OWID) << 5
        << setw(OWID) << 0 << endl
      
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
  else if (bdrynum == 6){ // all 6 boundaries
    // boundary 1
    ofs << setw(OWID) << 1 << endl
        << setw(OWID) << 1
        << setw(OWID) << 5
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
        << setw(OWID) << 5
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
        << setw(OWID) << 5
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
        << setw(OWID) << 5
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
        << setw(OWID) << 5
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
        << setw(OWID) << 5
        << setw(OWID) << 0 << endl
      
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

// bdrymum = 6 by default
void assembly::buildBoundary(const char* boundaryfile)
{
  std::ofstream ofs(boundaryfile);
  if(!ofs) { cout << "stream error!" << endl; exit(-1);}

  REAL x1,x2,y1,y2,z1,z2,x0,y0,z0;
  x1 = container.getMinCorner().getx();
  y1 = container.getMinCorner().gety();
  z1 = container.getMinCorner().getz();
  x2 = container.getMaxCorner().getx();
  y2 = container.getMaxCorner().gety();
  z2 = container.getMaxCorner().getz();
  x0 = container.getCenter().getx();
  y0 = container.getCenter().gety();
  z0 = container.getCenter().getz();

  int bdrynum = 6;

  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs << setw(OWID) << 0
      << setw(OWID) << bdrynum << endl << endl;

  // boundary 1
  ofs << setw(OWID) << 1 << endl
      << setw(OWID) << 1
      << setw(OWID) << 5
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
      << setw(OWID) << 5
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
      << setw(OWID) << 5
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
      << setw(OWID) << 5
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
      << setw(OWID) << 5
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
      << setw(OWID) << 5
      << setw(OWID) << 0 << endl
    
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

void assembly::trim(bool toRebuild,
		    const char* particlefile,
		    const char* trmparticle)
{
  if (toRebuild) readSample(particlefile);
  trimHistoryNum = TotalNum;

  REAL x1,x2,y1,y2,z1,z2,x0,y0,z0;
  x1 = container.getMinCorner().getx();
  y1 = container.getMinCorner().gety();
  z1 = container.getMinCorner().getz();
  x2 = container.getMaxCorner().getx();
  y2 = container.getMaxCorner().gety();
  z2 = container.getMaxCorner().getz();
  x0 = container.getCenter().getx();
  y0 = container.getCenter().gety();
  z0 = container.getCenter().getz();
 
  std::vector<particle*>::iterator itr;
  vec center;
  REAL mass = 0;

  for (itr = ParticleVec.begin(); itr != ParticleVec.end(); ){
    center=(*itr)->getCurrPosition();
    if(center.getx() <= x1 || center.getx() >= x2 ||
       center.gety() <= y1 || center.gety() >= y2 ||
       center.getz() <= z1 || center.getz() + gradInfo.getMaxPtclRadius() > z2)
      {
	delete (*itr); // release memory
	itr = ParticleVec.erase(itr); 
      }
    else
      ++itr;
  }
  
  for(itr=ParticleVec.begin();itr!=ParticleVec.end();++itr)
    mass += (*itr)->getMass();
  
  Volume = container.getDimx() * container.getDimy() * container.getDimz();
  BulkDensity = mass/Volume;
  
  TotalNum = ParticleVec.size();
  printParticle(trmparticle);
}


// make a cavity inside the sample and remove particles in the cavity
void assembly::trimCavity(bool toRebuild,
			  const char* particlefile,
			  const char* cavparticlefile)
{
  if (toRebuild) readSample(particlefile);
  trimHistoryNum = TotalNum;

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
 
  std::vector<particle*>::iterator itr;
  vec center;
  REAL mass = 0;
  REAL delta = gradInfo.getMaxPtclRadius();

  for (itr = ParticleVec.begin(); itr != ParticleVec.end(); ){
    center=(*itr)->getCurrPosition();
    if(center.getx() + delta  >= x1 && center.getx() - delta <= x2 &&
       center.gety() + delta  >= y1 && center.gety() - delta <= y2 &&
       center.getz() + delta  >= z1 && center.getz() - delta <= z2 )
      {
	delete (*itr); // release memory
	itr = ParticleVec.erase(itr); 
      }
    else
      ++itr;
  }
  
  for(itr=ParticleVec.begin();itr!=ParticleVec.end();++itr)
    mass += (*itr)->getMass();
  
  Volume = container.getDimx() * container.getDimy() * container.getDimz()
          -cavity.getDimx() * cavity.getDimy() * cavity.getDimz();
  BulkDensity = mass/Volume;
  
  TotalNum = ParticleVec.size();
  printParticle(cavparticlefile);
}


// expand partcile size by some percentage for particles inside cavity
void assembly::expandCavityParticles(bool toRebuild,
				     REAL percent,
				     const char* cavityptclfile,
				     const char* particlefile,
				     const char* newptclfile)
{
  if (toRebuild) readSample(particlefile);
  trimHistoryNum = TotalNum;

  REAL x1,x2,y1,y2,z1,z2;
  x1 = cavity.getMinCorner().getx();
  y1 = cavity.getMinCorner().gety();
  z1 = cavity.getMinCorner().getz();
  x2 = cavity.getMaxCorner().getx();
  y2 = cavity.getMaxCorner().gety();
  z2 = cavity.getMaxCorner().getz();
 
  std::vector<particle*>::iterator itr;
  vec center;

  int cavityPtclNum = 0;
  for (itr = ParticleVec.begin(); itr != ParticleVec.end(); ++itr ){
    center=(*itr)->getCurrPosition();
    if(center.getx() > x1 && center.getx() < x2 &&
       center.gety() > y1 && center.gety() < y2 &&
       center.getz() > z1 && center.getz() < z2 )
      ++cavityPtclNum;
  }

  printCavityParticle(cavityPtclNum, cavityptclfile);

  for (itr = ParticleVec.begin(); itr != ParticleVec.end(); ++itr ){
    center=(*itr)->getCurrPosition();
    if(center.getx() > x1 && center.getx() < x2 &&
       center.gety() > y1 && center.gety() < y2 &&
       center.getz() > z1 && center.getz() < z2 )
      (*itr)->expand(percent);
  }

  printParticle(newptclfile);
}


void assembly::printCavityParticle(int total, const char* str) const {
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
  
  vec tmp;
  std::vector<particle*>::const_iterator  it;
  for (it=ParticleVec.begin();it!=ParticleVec.end();++it)  {
    vec center=(*it)->getCurrPosition();
    if(center.getx() > x1 && center.getx() < x2 &&
       center.gety() > y1 && center.gety() < y2 &&
       center.getz() > z1 && center.getz() < z2 ) {

    ofs << setw(OWID) << (*it)->getID()
	<< setw(OWID) << (*it)->getType()
	<< setw(OWID) << (*it)->getA()
	<< setw(OWID) << (*it)->getB()
	<< setw(OWID) << (*it)->getC();
    
    tmp=(*it)->getCurrPosition();
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
    
    tmp=(*it)->getCurrVelocity();
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
// use the same BdryTgtMap.
void assembly::buildCavityBoundary(int existMaxId, const char* boundaryfile)
{
  std::ofstream ofs(boundaryfile);
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

  int bdrynum = 6;

  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs << setw(OWID) << 0
      << setw(OWID) << bdrynum << endl << endl;

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
void assembly::createMemParticle(REAL rRadius,
				 bool toRebuild,
				 const char* particlefile,
				 const char* allparticle)
{
  if (toRebuild) readSample(particlefile);

  REAL radius = gradInfo.getMinPtclRadius();
  if (gradInfo.getSize().size() == 1 &&
      gradInfo.getPtclRatioBA() == 1.0 && 
      gradInfo.getPtclRatioCA() == 1.0)
    radius *= rRadius; // determine how tiny the boundary particles are
  REAL diameter = radius*2;
  REAL x1,x2,y1,y2,z1,z2,x0,y0,z0;
  x1 = container.getMinCorner().getx();
  y1 = container.getMinCorner().gety();
  z1 = container.getMinCorner().getz();
  x2 = container.getMaxCorner().getx();
  y2 = container.getMaxCorner().gety();
  z2 = container.getMaxCorner().getz();
  x0 = container.getCenter().getx();
  y0 = container.getCenter().gety();
  z0 = container.getCenter().getz();

  particle* newptcl = NULL;
  REAL x, y, z;
  
  vector<particle*> vec1d;  // 1-dimension
  vector< vector<particle*>  > vec2d; // 2-dimension
  spring* newspring = NULL;
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
      newptcl = new particle(++memPtclIndex, 5, vec(x,y,z), radius, young, poisson);
      vec1d.push_back(newptcl);
      ParticleVec.push_back(newptcl);
    }
    vec2d.push_back(vec1d);
  }
  MemBoundary.push_back(vec2d);
  for (int i = 0; i != vec2d.size() ; ++i)
    for (int j = 0; j != vec2d[i].size() ; ++ j) {
      if (j + 1 < vec2d[i].size() ) {
	newspring = new spring(*vec2d[i][j], *vec2d[i][j+1], modulus);
	SpringVec.push_back(newspring);
      }
      if (i + 1 < vec2d.size() ) {
	newspring = new spring(*vec2d[i][j], *vec2d[i+1][j], modulus);
	SpringVec.push_back(newspring);
      }
    }

  // surface x2     
  vec2d.clear();
  x = x2 + radius;
  for (z = z1 + radius; z <= z2 - radius + EPS; z += diameter) {
    vec1d.clear();
    for (y = y1 + radius; y <= y2 - radius + EPS; y += diameter) {
      newptcl = new particle(++memPtclIndex, 5, vec(x,y,z), radius, young, poisson);
      vec1d.push_back(newptcl);
      ParticleVec.push_back(newptcl);
    }
    vec2d.push_back(vec1d);
  }
  MemBoundary.push_back(vec2d);
  for (int i = 0; i != vec2d.size() ; ++i)
    for (int j = 0; j != vec2d[i].size() ; ++ j) {
      if (j + 1 < vec2d[i].size() ) {
	newspring = new spring(*vec2d[i][j], *vec2d[i][j+1], modulus);
	SpringVec.push_back(newspring);
      }
      if (i + 1 < vec2d.size() ) {
	newspring = new spring(*vec2d[i][j], *vec2d[i+1][j], modulus);
	SpringVec.push_back(newspring);
      }
    }
 
  // surface y1
  vec2d.clear();
  y = y1 - radius;
  for (z = z1 + radius; z <= z2 - radius + EPS; z += diameter) {
    vec1d.clear();
    for (x = x1 + radius; x <= x2 - radius + EPS; x += diameter) {
      newptcl = new particle(++memPtclIndex, 5, vec(x,y,z), radius, young, poisson);
      vec1d.push_back(newptcl);
      ParticleVec.push_back(newptcl);
    }
    vec2d.push_back(vec1d);
  }
  MemBoundary.push_back(vec2d);
  for (int i = 0; i != vec2d.size() ; ++i)
    for (int j = 0; j != vec2d[i].size() ; ++ j) {
      if (j + 1 < vec2d[i].size() ) {
	newspring = new spring(*vec2d[i][j], *vec2d[i][j+1], modulus);
	SpringVec.push_back(newspring);
      }
      if (i + 1 < vec2d.size() ) {
	newspring = new spring(*vec2d[i][j], *vec2d[i+1][j], modulus);
	SpringVec.push_back(newspring);
      }
    }

  // surface y2
  vec2d.clear();
  y = y2 + radius;
  for (z = z1 + radius; z <= z2 - radius + EPS; z += diameter) {
    vec1d.clear();
    for (x = x1 + radius; x <= x2 - radius + EPS; x += diameter) {
      newptcl = new particle(++memPtclIndex, 5, vec(x,y,z), radius, young, poisson);
      vec1d.push_back(newptcl);
      ParticleVec.push_back(newptcl);
    }
    vec2d.push_back(vec1d);
  }
  MemBoundary.push_back(vec2d);
  for (int i = 0; i != vec2d.size() ; ++i)
    for (int j = 0; j != vec2d[i].size() ; ++ j) {
      if (j + 1 < vec2d[i].size() ) {
	newspring = new spring(*vec2d[i][j], *vec2d[i][j+1], modulus);
	SpringVec.push_back(newspring);
      }
      if (i + 1 < vec2d.size() ) {
	newspring = new spring(*vec2d[i][j], *vec2d[i+1][j], modulus);
	SpringVec.push_back(newspring);
      }
    }

  // surface z1
  vec2d.clear();
  z = z1 - radius;
  for (y = y1 + radius; y <= y2 - radius + EPS; y += diameter) {
    vec1d.clear();
    for (x = x1 + radius; x <= x2 - radius + EPS; x += diameter) {
      newptcl = new particle(++memPtclIndex, 5, vec(x,y,z), radius, young, poisson);
      vec1d.push_back(newptcl);
      ParticleVec.push_back(newptcl);
    }
    vec2d.push_back(vec1d);
  }
  MemBoundary.push_back(vec2d);
  for (int i = 0; i != vec2d.size() ; ++i)
    for (int j = 0; j != vec2d[i].size() ; ++ j) {
      if (j + 1 < vec2d[i].size() ) {
	newspring = new spring(*vec2d[i][j], *vec2d[i][j+1], modulus);
	SpringVec.push_back(newspring);
      }
      if (i + 1 < vec2d.size() ) {
	newspring = new spring(*vec2d[i][j], *vec2d[i+1][j], modulus);
	SpringVec.push_back(newspring);
      }
    }

  // surface z2
  vec2d.clear();
  z = z2 + radius;
  for (y = y1 + radius; y <= y2 - radius + EPS; y += diameter) {
    vec1d.clear();
    for (x = x1 + radius; x <= x2 - radius + EPS; x += diameter) {
      newptcl = new particle(++memPtclIndex, 5, vec(x,y,z), radius, young, poisson);
      vec1d.push_back(newptcl);
      ParticleVec.push_back(newptcl);
    }
    vec2d.push_back(vec1d);
  }
  MemBoundary.push_back(vec2d);
  for (int i = 0; i != vec2d.size() ; ++i)
    for (int j = 0; j != vec2d[i].size() ; ++ j) {
      if (j + 1 < vec2d[i].size() ) {
	newspring = new spring(*vec2d[i][j], *vec2d[i][j+1], modulus);
	SpringVec.push_back(newspring);
      }
      if (i + 1 < vec2d.size() ) {
	newspring = new spring(*vec2d[i][j], *vec2d[i+1][j], modulus);
	SpringVec.push_back(newspring);
      }
    }

  // membrane particles at the edges of each surface, for example,
  // x1y1 means particles on surface x1 connecting to particles on surface y1
  vector<particle *> x1y1;
  vector<particle *> x1y2;
  vector<particle *> x1z1;
  vector<particle *> x1z2;

  vector<particle *> x2y1;
  vector<particle *> x2y2;
  vector<particle *> x2z1;
  vector<particle *> x2z2;

  vector<particle *> y1x1;
  vector<particle *> y1x2;
  vector<particle *> y1z1;
  vector<particle *> y1z2;

  vector<particle *> y2x1;
  vector<particle *> y2x2;
  vector<particle *> y2z1;
  vector<particle *> y2z2;

  vector<particle *> z1x1;
  vector<particle *> z1x2;
  vector<particle *> z1y1;
  vector<particle *> z1y2;

  vector<particle *> z2x1;
  vector<particle *> z2x2;
  vector<particle *> z2y1;
  vector<particle *> z2y2;

  // find edge particles for each surface
  // MemBoundary[0, 1, 2, 3, 4, 5] correspond to 
  // surface     x1 x2 y1 y2 z1 z2 respectively
  // surface x1
  vec2d.clear();
  vec2d = MemBoundary[0];
  x1z1  = vec2d[0];
  x1z2  = vec2d[vec2d.size() - 1];
  for (int i = 0; i < vec2d.size(); ++i) {
    x1y1.push_back(vec2d[i][0]);
    x1y2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
  }
  // surface x2
  vec2d.clear();
  vec2d = MemBoundary[1];
  x2z1  = vec2d[0];
  x2z2  = vec2d[vec2d.size() - 1];
  for (int i = 0; i < vec2d.size(); ++i) {
    x2y1.push_back(vec2d[i][0]);
    x2y2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
  }
  // surface y1
  vec2d.clear();
  vec2d = MemBoundary[2];
  y1z1  = vec2d[0];
  y1z2  = vec2d[vec2d.size() - 1];
  for (int i = 0; i < vec2d.size(); ++i) {
    y1x1.push_back(vec2d[i][0]);
    y1x2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
  }
  // surface y2
  vec2d.clear();
  vec2d = MemBoundary[3];
  y2z1  = vec2d[0];
  y2z2  = vec2d[vec2d.size() - 1];
  for (int i = 0; i < vec2d.size(); ++i) {
    y2x1.push_back(vec2d[i][0]);
    y2x2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
  }
  // surface z1
  vec2d.clear();
  vec2d = MemBoundary[4];
  z1y1  = vec2d[0];
  z1y2  = vec2d[vec2d.size() - 1];
  for (int i = 0; i < vec2d.size(); ++i) {
    z1x1.push_back(vec2d[i][0]);
    z1x2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
  }
  // surface z2
  vec2d.clear();
  vec2d = MemBoundary[5];
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
    newspring = new spring(*x1y1[i], *y1x1[i], modulus);
    SpringVec.push_back(newspring);
  }
  assert(x1y2.size() == y2x1.size());
  for (int i = 0; i < x1y2.size(); ++i) {
    newspring = new spring(*x1y2[i], *y2x1[i], modulus);
    SpringVec.push_back(newspring);
  }
  assert(x1z1.size() == z1x1.size());
  for (int i = 0; i < x1z1.size(); ++i) {
    newspring = new spring(*x1z1[i], *z1x1[i], modulus);
    SpringVec.push_back(newspring);
  }
  assert(x1z2.size() == z2x1.size());
  for (int i = 0; i < x1z2.size(); ++i) {
    newspring = new spring(*x1z2[i], *z2x1[i], modulus);
    SpringVec.push_back(newspring);
  }
  // 4 edges on surface x2  
  assert(x2y1.size() == y1x2.size());
  for (int i = 0; i < x2y1.size(); ++i) {
    newspring = new spring(*x2y1[i], *y1x2[i], modulus);
    SpringVec.push_back(newspring);
  }
  assert(x2y2.size() == y2x2.size());
  for (int i = 0; i < x2y2.size(); ++i) {
    newspring = new spring(*x2y2[i], *y2x2[i], modulus);
    SpringVec.push_back(newspring);
  }
  assert(x2z1.size() == z1x2.size());
  for (int i = 0; i < x2z1.size(); ++i) {
    newspring = new spring(*x2z1[i], *z1x2[i], modulus);
    SpringVec.push_back(newspring);
  }
  assert(x2z2.size() == z2x2.size());
  for (int i = 0; i < x2z2.size(); ++i) {
    newspring = new spring(*x2z2[i], *z2x2[i], modulus);
    SpringVec.push_back(newspring);
  }
  // 2 edges on surface y1 
  assert(y1z1.size() == z1y1.size());
  for (int i = 0; i < y1z1.size(); ++i) {
    newspring = new spring(*y1z1[i], *z1y1[i], modulus);
    SpringVec.push_back(newspring);
  }
  assert(y1z2.size() == z2y1.size());
  for (int i = 0; i < y1z2.size(); ++i) {
    newspring = new spring(*y1z2[i], *z2y1[i], modulus);
    SpringVec.push_back(newspring);
  }
  // 2 edges on surface y2
  assert(y2z1.size() == z1y2.size());
  for (int i = 0; i < y2z1.size(); ++i) {
    newspring = new spring(*y2z1[i], *z1y2[i], modulus);
    SpringVec.push_back(newspring);
  }
  assert(y2z2.size() == z2y2.size());
  for (int i = 0; i < y2z2.size(); ++i) {
    newspring = new spring(*y2z2[i], *z2y2[i], modulus);
    SpringVec.push_back(newspring);
  }

  TotalNum = ParticleVec.size();
  printParticle(allparticle);

}


void assembly::TrimPtclBdryByHeight(REAL height,
			    const char* iniptclfile,
			    const char* particlefile)
{
  readSample(iniptclfile);
  
  std::vector<particle*>::iterator itr;
  for (itr = ParticleVec.begin(); itr != ParticleVec.end(); ){
    if ( (*itr)->getType() == 1 ) { // 1-fixed
      vec center=(*itr)->getCurrPosition();
      if(center.getz() > height)
	{
	  delete (*itr); // release memory
	  itr = ParticleVec.erase(itr); 
	}
      else {
	(*itr)->setType(10); // 10-ghost
	++itr;
      }
    }
  }
  
  TotalNum = ParticleVec.size();  
  printParticle(particlefile);
}


// deposit floating particles into a container through applying gravity,
// the container can be as simple as a bottom plate
void assembly::deposit(int   total_steps,  
		       int   snapshots,
		       int   interval,
		       const char* iniptclfile,   
		       const char* inibdryfile,
		       const char* particlefile, 
		       const char* contactfile,
		       const char* progressfile, 
		       const char* debugfile)
{
    // pre_1: open streams for output.
    // particlefile and contactfile are used for snapshots at the end.
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
//		<< setw(OWID) << "sigma_11"
//		<< setw(OWID) << "sigma_12"
//		<< setw(OWID) << "sigma_13"
//		<< setw(OWID) << "sigma_21"
//		<< setw(OWID) << "sigma_22"
//		<< setw(OWID) << "sigma_23"
//		<< setw(OWID) << "sigma_31"
//		<< setw(OWID) << "sigma_32"
//		<< setw(OWID) << "sigma_33"
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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1); }
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from existing files.
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries.

    // pre_3: define variables used in iterations.
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int  stepsnum=0;
    char stepsstr[4];
    char stepsfp[50];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
//    matrix granularStress;	// granular stress
    int  bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    // iterations starting ...
    g_iteration=0; 
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
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	
///////////////////////////////////////////////
// this change is for generating a 2D sample
	    // should fix the dem particle in y, i.e. force_j = 0; moment_i=moment_k = 0
	    for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
		(*pt)->fixY();
	    }
//////////////////////////////////////////////

	// 5. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();
	//gettimeofday(&time_w2,NULL);
	//cout << setw(OWID) << timediffsec(time_w1,time_w2) << endl;

	// 6. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx(); Volume=l13*l24*l56;
	
	// granular stress calculatioin
//	granularStress.clear();
//	granularStress = getGranularStress();
	void_ratio=Volume/getParticleVolume()-1;

	// 7. (1) output particles and contacts information as snapshots.
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    time(&timeStamp);
	    g_timeinf << setw(4) << stepsnum << " " << ctime(&timeStamp) << flush;
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info.
	if (g_iteration % interval == 0) {
	    gettimeofday(&time_w2,NULL);
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/TotalNum
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
//			<< setw(OWID) << granularStress(1,1)
//			<< setw(OWID) << granularStress(1,2)
//			<< setw(OWID) << granularStress(1,3)
//			<< setw(OWID) << granularStress(2,1)
//			<< setw(OWID) << granularStress(2,2)
//			<< setw(OWID) << granularStress(2,3)
//			<< setw(OWID) << granularStress(3,1)
//			<< setw(OWID) << granularStress(3,2)
//			<< setw(OWID) << granularStress(3,3)
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

	    /*
	    g_debuginf << setw(OWID) << g_iteration
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
	    */

	}

	// 8. loop break conditions.

    } while (++g_iteration < total_steps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    g_timeinf << setw(4) << "end" << " " << ctime(&timeStamp) << flush;

    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


void assembly::depositAfterCavity(int   total_steps,  
				  int   snapshots,
				  int   interval,
				  const char* iniptclfile,   
				  const char* inibdryfile,
				  const char* inicavefile,
				  const char* particlefile, 
				  const char* contactfile,
				  const char* progressfile, 
				  const char* debugfile)
{
    // pre_1: open streams for output.
    // particlefile and contactfile are used for snapshots at the end.
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
//		<< setw(OWID) << "sigma_11"
//		<< setw(OWID) << "sigma_12"
//		<< setw(OWID) << "sigma_13"
//		<< setw(OWID) << "sigma_21"
//		<< setw(OWID) << "sigma_22"
//		<< setw(OWID) << "sigma_23"
//		<< setw(OWID) << "sigma_31"
//		<< setw(OWID) << "sigma_32"
//		<< setw(OWID) << "sigma_33"
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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1); }
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from existing files.
    readSample(iniptclfile);   // create container and particles, velocity and omga are set zero. 
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
//    matrix granularStress;	// granular stress
    int  bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    // iterations starting ...
    g_iteration=0; 
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
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	cavityBoundaryForce();
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();

	// 6. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx(); Volume=l13*l24*l56;

	// granular stress calculation
//	granularStress.clear();
//	granularStress = getGranularStress();
	void_ratio=Volume/getParticleVolume()-1;

	// 7. (1) output particles and contacts information as snapshots.
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    time(&timeStamp);
	    g_timeinf << setw(4) << stepsnum << " " << ctime(&timeStamp) << flush;
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info.
	if (g_iteration % interval == 0) {
	    gettimeofday(&time_w2,NULL);
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/TotalNum
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
//			<< setw(OWID) << granularStress(1,1)
//			<< setw(OWID) << granularStress(1,2)
//			<< setw(OWID) << granularStress(1,3)
//			<< setw(OWID) << granularStress(2,1)
//			<< setw(OWID) << granularStress(2,2)
//			<< setw(OWID) << granularStress(2,3)
//			<< setw(OWID) << granularStress(3,1)
//			<< setw(OWID) << granularStress(3,2)
//			<< setw(OWID) << granularStress(3,3)
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

	    /*
	    g_debuginf << setw(OWID) << g_iteration
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
	    */

	}

	// 8. loop break conditions.

    } while (++g_iteration < total_steps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    g_timeinf << setw(4) << "end" << " " << ctime(&timeStamp) << flush;

    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


void assembly::deGravitation(int   total_steps,  
			     int   snapshots,
			     int   interval,
			     bool  toRebuild,
			     const char* iniptclfile,   
			     const char* particlefile, 
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
  
  g_debuginf.open(debugfile);
  if(!g_debuginf) { cout << "stream error!" << endl; exit(-1); }
  g_debuginf.setf(std::ios::scientific, std::ios::floatfield);
  
  // pre_2. create particles from existing files.
  if (toRebuild) readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
  
  // pre_3. define variables used in iterations
  REAL avgNormal=0;
  REAL avgTangt=0;
  int  stepsnum=0;
  char stepsstr[4];
  char stepsfp[50];

  // iterations starting ...
  g_iteration=0; 
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
      
      // 5. (1) output particles and contacts information as snapshots.
      if (g_iteration % (total_steps/snapshots) == 0){
	sprintf(stepsstr, "%03d", stepsnum); 
	strcpy(stepsfp, particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	printParticle(stepsfp);
	
	sprintf(stepsstr, "%03d", stepsnum); 
	strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	printContact(stepsfp);
	time(&timeStamp);
	g_timeinf << setw(4) << "end" << " " << ctime(&timeStamp) << flush;
	++stepsnum;
      }
      
      // 5. (2) output stress and strain info.
      if (g_iteration % interval == 0) {
	progressinf << setw(OWID) << g_iteration
		    << setw(OWID) << getPossCntctNum()
		    << setw(OWID) << getActualCntctNum()
		    << setw(OWID) << getAveragePenetration()
		    << setw(OWID) << avgNormal
		    << setw(OWID) << avgTangt
		    << setw(OWID) << getAverageVelocity() 
		    << setw(OWID) << getAverageOmga()
		    << setw(OWID) << getAverageForce()   
		    << setw(OWID) << getAverageMoment()
		    << setw(OWID) << getTransEnergy()
		    << setw(OWID) << getRotatEnergy()
		    << setw(OWID) << getKinetEnergy()
		    << endl;
      }
      
      // 7. loop break conditions.
      if (ContactVec.size() == 0) break;
      
    } while (++g_iteration < total_steps);
  
  // post_1. store the final snapshot of particles & contacts.
  strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
  printParticle(stepsfp);
  
  strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
  printContact(stepsfp);
  g_timeinf << setw(4) << "end" << " " << ctime(&timeStamp) << flush;
  
  // post_2. close streams
  progressinf.close();
  g_debuginf.close();
}


// actual deposit function for the case of fixed particle boundaries
void assembly::deposit_p(int   total_steps,  
			 int   snapshots,
			 int   interval,
			 REAL dimn,
			 REAL rsize,
			 const char* iniptclfile,   
			 const char* particlefile, 
			 const char* contactfile,
			 const char* progressfile, 
			 const char* debugfile)
{
    // pre_1: open streams for output.
    // particlefile and contactfile are used for snapshots at the end.
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
	        << "sigma3_1        sigma3_2        "
//		<< "sigma_21        sigma_22        sigma_23        sigma_31        sigma_32        sigma_33        "
		<< "p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
	        << "epsilon_v" << endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1); }
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from existing files.
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 

    // pre_3: define variables used in iterations.
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    REAL void_ratio=0;
//    matrix granularStress;	// granular stress
    // iterations starting ...
    g_iteration=0;
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
	Volume=l13*l24*l56;
	void_ratio=Volume/getParticleVolume()-1;

	// granular stress calculation
//	granularStress.clear();
//	granularStress = getGranularStress();

	// 6. (1) output particles and contacts information as snapshots.
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 6. (2) output statistics info.
	if (g_iteration % interval == 0) {
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*getActualCntctNum()/TotalNum
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
		        << setw(OWID) << "0"
//			<< setw(OWID) << granularStress(1,1)
//			<< setw(OWID) << granularStress(1,2)
//			<< setw(OWID) << granularStress(1,3)
//			<< setw(OWID) << granularStress(2,1)
//			<< setw(OWID) << granularStress(2,2)
//			<< setw(OWID) << granularStress(2,3)
//			<< setw(OWID) << granularStress(3,1)
//			<< setw(OWID) << granularStress(3,2)
//			<< setw(OWID) << granularStress(3,3)
		        << endl;
	}

	// 7. loop break conditions.


    } while (++g_iteration < total_steps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


// squeeze paticles inside a container by moving the boundaries
void assembly::squeeze(int   total_steps,  
		       int   init_steps,
		       int   snapshots,
		       int   interval,
		       int   flag,
		       const char* iniptclfile,   
		       const char* inibdryfile,
		       const char* particlefile, 
		       const char* boundaryfile,
		       const char* contactfile,
		       const char* progressfile, 
		       const char* debugfile)
{
    // pre_1: open streams for output.
    // particlefile and contactfile are used for snapshots at the end.
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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1); }
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from existing files.
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
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
    g_iteration=0;
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
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();

	// 6. calculate sample void ratio.
	l56=getTopFreeParticlePosition().getz() -getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx(); Volume=l13*l24*l56;
	void_ratio=Volume/getParticleVolume()-1;

	// displacement control
	if (g_iteration > init_steps) {
	    if (flag==1) // loosen, totally remove the wall
		midctl[0].tran=vec(TIMESTEP*1.0e+0*flag,0,0);
	    else         // squeeze
		midctl[0].tran=vec(TIMESTEP*5.0e-3*flag,0,0);
	    updateRB(mid,midctl,2);
	}

	// 7. (1) output particles and contacts information as snapshots.
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info.
	if (g_iteration % interval == 0) {
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/TotalNum
		        << endl;
	    g_debuginf << setw(OWID) << g_iteration
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

    } while (++g_iteration < total_steps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);

    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


// Isotropically compress floating particles to a specific confining pressure, which is usually a low
// value in order to create an intial status. Force boundaries are used. This process may be not 
// physically true.
void assembly::isotropic(int   total_steps,
			 int   snapshots, 
			 int   interval,
			 REAL  sigma,			  
			 const char* iniptclfile,   
			 const char* inibdryfile,
			 const char* particlefile, 
			 const char* boundaryfile,
			 const char* contactfile,  
			 const char* progressfile,
			 const char* balancedfile, 
			 const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout << "stream error!" << endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
/*    progressinf << "isotropic..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
		<< "sample          sample         sample           " // for new 3 sigma's
		<< "sample          sample         sample           sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          "
		<< "sample          sample         sample           sample          sample          sample     "
	        << "     sample          sample          sample          "
		<< "void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
		<< "sigma1_1        sigma1_2        sigma2_1        "
	        << "sigma2_2        sigma3_1        sigma3_2        "
	        << "sigma_11        sigma_12        sigma_13        sigma_21        "
	        << "sigma_22        sigma_23        sigma_31        sigma_32        sigma_33        "   
		<< "p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v       "
		<< "epsilon_11        epsilon_12        epsilon_13        epsilon_21        "
	        << "epsilon_22        epsilon_23        epsilon_31        epsilon_32        epsilon_33        "
	        << "        ratio          porosity         number"
	        << endl;
*/
    progressinf << setw(OWID) << "isotropic..." << endl
		<< setw(OWID) << "iteration"
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
		<< setw(OWID) << "sample" // for new 3 sigma's
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
		<< setw(OWID) << "sample"
	        << setw(OWID) << "void"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "coordinate" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag"
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" << endl
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
	        << setw(OWID) << "sigma_11"
	        << setw(OWID) << "sigma_12"
	        << setw(OWID) << "sigma_13"
	        << setw(OWID) << "sigma_21"
	        << setw(OWID) << "sigma_22"
	        << setw(OWID) << "sigma_23"
	        << setw(OWID) << "sigma_31"
	        << setw(OWID) << "sigma_32"
	        << setw(OWID) << "sigma_33"
	        << setw(OWID) << "mean_stress"
	        << setw(OWID) << "width"
	        << setw(OWID) << "length"
	        << setw(OWID) << "height"
	        << setw(OWID) << "volume"
	        << setw(OWID) << "epsilon_w"
	        << setw(OWID) << "epsilon_l"
	        << setw(OWID) << "epsilon_h"
	        << setw(OWID) << "epsilon_v"
	        << setw(OWID) << "epsilon_11"
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "ratio"
	        << setw(OWID) << "porosity"
	        << setw(OWID) << "number" 
	        << setw(OWID) << "epsilon_11"	// for finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "epsilon_11"	// for mixed finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "epsilon_11"	// for eulerian strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33" 
	        << setw(OWID) << "epsilon_11"	// for mixed eulerian strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33" << endl;


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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries
//    readTesse("tess_info");	// create edge set based on tessellation file from qhull
    // pre_3: define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    initVolume = W0*L0*H0;
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    matrix granularStress; // granular stress
    matrix granularStrain; // granular strain
    matrix previousStrain(3,3);
    matrix previousEuler(3,3);
    matrix finiteStrain(3,3);	// finite granular strain
    matrix lag_HOT(3,3);	// mixed finite granular strain
    matrix eulerianStrain(3,3);	// eulerian granular strain
    matrix euler_HOT(3,3);
    matrix previousEuler_HOT(3,3);

    // initialize previousStrain(3,3)
    for(int i_ps=0; i_ps!=2; i_ps++){
	for(int j_ps=0; j_ps!=2; j_ps++){
		previousStrain(i_ps+1,j_ps+1) = 0;
		previousEuler(i_ps+1,j_ps+1) = 0;
		previousEuler_HOT(i_ps+1,j_ps+1) = 0;
	}
    }
    
    // first tessellation
    createInputForQhull();
    callQhull();
    readTesse("tess_info");
    readTesse_finite("tess_info");  
    cellVec_init = cellVec;	// keep cellVec_init for lagrangian strain unchanged  

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
    g_iteration=0;
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
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);

	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	
	// calculate granular stress
	granularStress.clear();
	granularStress = getGranularStress();
//	// calculate granular strain
//	granularStrain.clear();
//	granularStrain = getGranularStrain();

	void_ratio=Volume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	if (sigma2_1<sigma)
	    midctl[0].tran=vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma)
	    midctl[1].tran=vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma)
	    maxctl[0].tran=vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma)
	    maxctl[1].tran=vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=vec(0,-TIMESTEP*RELEASE_RATE,0);

/*
// below is adjustment of the velocities of boundary walls is based on new stress calculation

	// along the z direction	
	if (fabs(sigma_33)<sigma){  // means we can continue to compress
		minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);  // minus is because this is for the upper wall, it needs to go down
		minctl[1].tran=vec(0,0,TIMESTEP*COMPRESS_RATE);
	}
	else{  // means we need to release
		minctl[0].tran=vec(0,0,TIMESTEP*RELEASE_RATE);
		minctl[1].tran=vec(0,0,-TIMESTEP*RELEASE_RATE);
	}
	// along the x direction
	if (fabs(sigma_11)<sigma){
		midctl[0].tran=vec(-TIMESTEP*COMPRESS_RATE,0,0);
		midctl[1].tran=vec(TIMESTEP*COMPRESS_RATE,0,0);	
	}
	else{
		midctl[0].tran=vec(TIMESTEP*RELEASE_RATE,0,0);
		midctl[1].tran=vec(-TIMESTEP*RELEASE_RATE,0,0);
	}
	// along the y direction
	if (fabs(sigma_22)<sigma){
		maxctl[0].tran=vec(0,-TIMESTEP*COMPRESS_RATE,0);
		maxctl[1].tran=vec(0,TIMESTEP*COMPRESS_RATE,0);
	}
	else{
		maxctl[0].tran=vec(0,TIMESTEP*RELEASE_RATE,0);
		maxctl[1].tran=vec(0,-TIMESTEP*RELEASE_RATE,0);
	}

// above is modification about the adjustment of velocities
*/
	updateRB(min,minctl,2);
	updateRB(mid,midctl,2);
	updateRB(max,maxctl,2);
	updateRB6();
	
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0 ){
	    // calculate granular strain and Euler strain
	    granularStrain.clear();
	    eulerianStrain.clear();
	    euler_HOT.clear();
	    granularStrain = getGranularStrain()+previousStrain;
	    eulerianStrain = getEulerianStrain()+previousEuler;
	    euler_HOT = getEuler_HOT()+previousEuler_HOT;
	    if(eulerianStrain(1,1)-previousEuler(1,1)>strainThreshold ||
	       eulerianStrain(2,2)-previousEuler(2,2)>strainThreshold ||
	       eulerianStrain(3,3)-previousEuler(3,3)>strainThreshold ){	// need to tessellate again
			resetStartPosition();	// reset initial position
			// tessellate again
			createInputForQhull();
			callQhull();
			readTesse("tess_info");
			readTesse_finite("tess_info");
			// update previousStrain
			previousStrain = granularStrain;
			previousEuler = eulerianStrain;
			previousEuler_HOT = euler_HOT;
	    }
	    // calculate finite granular strain
	    finiteStrain.clear();
	    finiteStrain = getFiniteStrain();
	    // calculate mixed finite granular strain
	    lag_HOT.clear();
	    lag_HOT = getHigherStrain();

	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
			<< setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2 << setw(OWID) << sigma2_1
			<< setw(OWID) << sigma2_2 << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << granularStress(1,1) << setw(OWID) << granularStress(1,2) << setw(OWID) << granularStress(1,3)
		        << setw(OWID) << granularStress(2,1) << setw(OWID) << granularStress(2,2) << setw(OWID) << granularStress(2,3)
		        << setw(OWID) << granularStress(3,1) << setw(OWID) << granularStress(3,2) << setw(OWID) << granularStress(3,3)
		        << setw(OWID) << getAverageRigidPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
			<< setw(OWID) << granularStrain(1,1) << setw(OWID) << granularStrain(1,2) << setw(OWID) << granularStrain(1,3)
		        << setw(OWID) << granularStrain(2,1) << setw(OWID) << granularStrain(2,2) << setw(OWID) << granularStrain(2,3)
		        << setw(OWID) << granularStrain(3,1) << setw(OWID) << granularStrain(3,2) << setw(OWID) << granularStrain(3,3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
			<< setw(OWID) << finiteStrain(1,1) << setw(OWID) << finiteStrain(1,2) << setw(OWID) << finiteStrain(1,3)
		        << setw(OWID) << finiteStrain(2,1) << setw(OWID) << finiteStrain(2,2) << setw(OWID) << finiteStrain(2,3)
		        << setw(OWID) << finiteStrain(3,1) << setw(OWID) << finiteStrain(3,2) << setw(OWID) << finiteStrain(3,3)
			<< setw(OWID) << lag_HOT(1,1) << setw(OWID) << lag_HOT(1,2) << setw(OWID) << lag_HOT(1,3)
		        << setw(OWID) << lag_HOT(2,1) << setw(OWID) << lag_HOT(2,2) << setw(OWID) << lag_HOT(2,3)
		        << setw(OWID) << lag_HOT(3,1) << setw(OWID) << lag_HOT(3,2) << setw(OWID) << lag_HOT(3,3)
			<< setw(OWID) << eulerianStrain(1,1) << setw(OWID) << eulerianStrain(1,2) << setw(OWID) << eulerianStrain(1,3)
		        << setw(OWID) << eulerianStrain(2,1) << setw(OWID) << eulerianStrain(2,2) << setw(OWID) << eulerianStrain(2,3)
		        << setw(OWID) << eulerianStrain(3,1) << setw(OWID) << eulerianStrain(3,2) << setw(OWID) << eulerianStrain(3,3)
			<< setw(OWID) << euler_HOT(1,1) << setw(OWID) << euler_HOT(1,2) << setw(OWID) << euler_HOT(1,3)
		        << setw(OWID) << euler_HOT(2,1) << setw(OWID) << euler_HOT(2,2) << setw(OWID) << euler_HOT(2,3)
		        << setw(OWID) << euler_HOT(3,1) << setw(OWID) << euler_HOT(3,2) << setw(OWID) << euler_HOT(3,3)
		        << endl;
	    g_debuginf << setw(OWID) << g_iteration
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
/*
// below is based on new stress calculation
	if (	fabs(fabs(granularStress(1,1))-sigma)/sigma < STRESS_ERROR
	    &&  fabs(fabs(granularStress(2,2))-sigma)/sigma < STRESS_ERROR
	    &&  fabs(fabs(granularStress(3,3))-sigma)/sigma < STRESS_ERROR ) {
// above is based on new stress calculation
*/
/*	    balancedinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()    
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		        << endl;*/
	    // calculate granular strain
	    granularStrain.clear();
	    granularStrain = getGranularStrain()+previousStrain;
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()    
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2 << setw(OWID) << sigma2_1
			<< setw(OWID) << sigma2_2 << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << granularStress(1,1) << setw(OWID) << granularStress(1,2) << setw(OWID) << granularStress(1,3)
		        << setw(OWID) << granularStress(2,1) << setw(OWID) << granularStress(2,2) << setw(OWID) << granularStress(2,3)
		        << setw(OWID) << granularStress(3,1) << setw(OWID) << granularStress(3,2) << setw(OWID) << granularStress(3,3)
		        << setw(OWID) << getAverageRigidPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
			<< setw(OWID) << granularStrain(1,1) << setw(OWID) << granularStrain(1,2) << setw(OWID) << granularStrain(1,3)
		        << setw(OWID) << granularStrain(2,1) << setw(OWID) << granularStrain(2,2) << setw(OWID) << granularStrain(2,3)
		        << setw(OWID) << granularStrain(3,1) << setw(OWID) << granularStrain(3,2) << setw(OWID) << granularStrain(3,3)
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
			<< setw(OWID) << finiteStrain(1,1) << setw(OWID) << finiteStrain(1,2) << setw(OWID) << finiteStrain(1,3)
		        << setw(OWID) << finiteStrain(2,1) << setw(OWID) << finiteStrain(2,2) << setw(OWID) << finiteStrain(2,3)
		        << setw(OWID) << finiteStrain(3,1) << setw(OWID) << finiteStrain(3,2) << setw(OWID) << finiteStrain(3,3)
			<< setw(OWID) << lag_HOT(1,1) << setw(OWID) << lag_HOT(1,2) << setw(OWID) << lag_HOT(1,3)
		        << setw(OWID) << lag_HOT(2,1) << setw(OWID) << lag_HOT(2,2) << setw(OWID) << lag_HOT(2,3)
		        << setw(OWID) << lag_HOT(3,1) << setw(OWID) << lag_HOT(3,2) << setw(OWID) << lag_HOT(3,3)
			<< setw(OWID) << eulerianStrain(1,1) << setw(OWID) << eulerianStrain(1,2) << setw(OWID) << eulerianStrain(1,3)
		        << setw(OWID) << eulerianStrain(2,1) << setw(OWID) << eulerianStrain(2,2) << setw(OWID) << eulerianStrain(2,3)
		        << setw(OWID) << eulerianStrain(3,1) << setw(OWID) << eulerianStrain(3,2) << setw(OWID) << eulerianStrain(3,3)
			<< setw(OWID) << euler_HOT(1,1) << setw(OWID) << euler_HOT(1,2) << setw(OWID) << euler_HOT(1,3)
		        << setw(OWID) << euler_HOT(2,1) << setw(OWID) << euler_HOT(2,2) << setw(OWID) << euler_HOT(2,3)
		        << setw(OWID) << euler_HOT(3,1) << setw(OWID) << euler_HOT(3,2) << setw(OWID) << euler_HOT(3,3)
		        << endl;
	    break;
	}

    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile);  strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_a. This function
// increases confining pressure step by step to sigma_b, making it possible to find equilibrium 
// state where particle pressure equals confining pressure. Force boundaries are used
void assembly::isotropic(int   total_steps,
			 int   snapshots, 
			 int   interval,
			 REAL sigma_a,
			 REAL sigma_b,
			 int   sigma_division,
			 const char* iniptclfile,   
			 const char* inibdryfile,
			 const char* particlefile, 
			 const char* boundaryfile,
			 const char* contactfile,  
			 const char* progressfile,
			 const char* balancedfile, 
			 const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout << "stream error!" << endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
/*   
 progressinf << "isotropic..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
		<< "     sample          sample         sample           sample          sample         sample"
		<< "          sample         sample           sample           "           
	        << "sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           "
		<< "sigma_11        sigma_12        sigma_13        sigma_21        sigma_22        sigma_23           "
		<< "sigma_31        sigma_32        sigma_33           "
		<< "p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	        << "        ratio          porosity         number"
	        << endl;
*/
    progressinf << setw(OWID) << "isotropic..." << endl
		<< setw(OWID) << "iteration"
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
		<< setw(OWID) << "sample" // for new 3 sigma's
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
		<< setw(OWID) << "sample"
	        << setw(OWID) << "void"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "coordinate" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag"
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "Logarithmic" 	// for logarithmic strain
		<< setw(OWID) << "Logarithmic"
		<< setw(OWID) << "Logarithmic"
		<< setw(OWID) << "Logarithmic" 
		<< setw(OWID) << "Bagi" 	// for average_dudx_Bagi, used to test quadratic terms, April 22, 2013
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Lagrangian" 	// for average_dudx_Lagrangian, used to test quadratic terms
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian"
		<< setw(OWID) << "Eulerian" 	// for average_dudx_Eulerian, used to test quadratic terms
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" << endl
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
	        << setw(OWID) << "sigma_11"
	        << setw(OWID) << "sigma_12"
	        << setw(OWID) << "sigma_13"
	        << setw(OWID) << "sigma_21"
	        << setw(OWID) << "sigma_22"
	        << setw(OWID) << "sigma_23"
	        << setw(OWID) << "sigma_31"
	        << setw(OWID) << "sigma_32"
	        << setw(OWID) << "sigma_33"
	        << setw(OWID) << "mean_stress"
	        << setw(OWID) << "width"
	        << setw(OWID) << "length"
	        << setw(OWID) << "height"
	        << setw(OWID) << "volume"
	        << setw(OWID) << "epsilon_w"
	        << setw(OWID) << "epsilon_l"
	        << setw(OWID) << "epsilon_h"
	        << setw(OWID) << "epsilon_v"
	        << setw(OWID) << "epsilon_11"
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "ratio"
	        << setw(OWID) << "porosity"
	        << setw(OWID) << "number" 
	        << setw(OWID) << "epsilon_11"	// for finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "epsilon_11"	// for mixed finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "epsilon_11"	// for eulerian finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33" 
	        << setw(OWID) << "epsilon_11"	// for mixed eulerian finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33" 
		<< setw(OWID) << "epsilon_w"	// for logranthmic strain
	        << setw(OWID) << "epsilon_l"
	        << setw(OWID) << "epsilon_h"
	        << setw(OWID) << "epsilon_v" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Bagi, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Lagrangian, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Eulerian, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" << endl;


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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    std::vector<REAL> temp_vec;	// used to test quadratic terms, April 22, 2013
    temp_vec.push_back(0);
    temp_vec.push_back(0);
    temp_vec.push_back(0);
    for(int i=0; i!=3; i++){
	average_dudx_Bagi.appendRow(temp_vec);
	average_dudx_Eulerian.appendRow(temp_vec);
	average_dudx_Lagrangian.appendRow(temp_vec);
    }

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3: define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    initVolume = W0*L0*H0;
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    matrix granularStress;	// granular stress
    matrix granularStrain; // granular strain
    matrix previousStrain(3,3);
    matrix previousEuler(3,3);
    matrix finiteStrain(3,3);	// finite granular strain
    matrix lag_HOT(3,3);	// mixed finite granular strain
    matrix eulerianStrain(3,3);	// eulerian granular strain
    matrix euler_HOT(3,3);
    matrix previousEuler_HOT(3,3);
    matrix previousEuler_dudx(3,3);	// used to test quadratic terms, April 22, 2013
    matrix previousBagi_dudx(3,3);	

    // initialize previousStrain(3,3)
    for(int i_ps=0; i_ps!=2; i_ps++){
	for(int j_ps=0; j_ps!=2; j_ps++){
		previousStrain(i_ps+1,j_ps+1) = 0;
		previousEuler(i_ps+1,j_ps+1) = 0;
		previousEuler_HOT(i_ps+1,j_ps+1) = 0;
		previousEuler_dudx(i_ps+1,j_ps+1) = 0;	// used to test quadratic terms
		previousBagi_dudx(i_ps+1,j_ps+1) = 0;
	}
    }
    
    // first tessellation
    createInputForQhull();
    callQhull();
    readTesse("tess_info");
    readTesse_finite("tess_info");    
    cellVec_init = cellVec;	// keep cellVec_init for lagrangian strain unchanged 

    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL epsilon_w_log, epsilon_l_log, epsilon_h_log;	// lograthmic strain
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
    g_iteration=0;
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
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	// granular stress calculation
	granularStress.clear();
	granularStress = getGranularStress();
	void_ratio=Volume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	if (sigma2_1<sigma)
	    midctl[0].tran=vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma)
	    midctl[1].tran=vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma)
	    maxctl[0].tran=vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma)
	    maxctl[1].tran=vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=vec(0,-TIMESTEP*RELEASE_RATE,0);
	
	updateRB(min,minctl,2);
	updateRB(mid,midctl,2);
	updateRB(max,maxctl,2);
	updateRB6();
	
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	epsilon_w_log = log(l24/W0); epsilon_l_log = log(l13/L0); epsilon_h_log = log(l56/H0);
	if (g_iteration % interval == 0 ){
	    // calculate granular strain and Euler strain
	    average_dudx_Bagi.clear();	// used to test quadratic terms
	    average_dudx_Eulerian.clear();

	    granularStrain.clear();
	    eulerianStrain.clear();
	    euler_HOT.clear();

	    granularStrain = getGranularStrain()+previousStrain;
	    eulerianStrain = getEulerianStrain()+previousEuler;  
	    euler_HOT = getEuler_HOT()+previousEuler_HOT;

	    average_dudx_Bagi = average_dudx_Bagi+previousBagi_dudx;	// used to test quadratic terms
	    average_dudx_Eulerian = average_dudx_Eulerian+previousEuler_dudx;

	    if(eulerianStrain(1,1)-previousEuler(1,1)>strainThreshold ||
	       eulerianStrain(2,2)-previousEuler(2,2)>strainThreshold ||
	       eulerianStrain(3,3)-previousEuler(3,3)>strainThreshold  ){	// need to tessellate again
			resetStartPosition();	// reset initial position
			// tessellate again
			createInputForQhull();
			callQhull();
			readTesse("tess_info");
			readTesse_finite("tess_info");
			// update previousStrain
			previousStrain = granularStrain;
			previousEuler = eulerianStrain;
			previousEuler_HOT = euler_HOT;
			previousBagi_dudx = average_dudx_Bagi;	// used to test quadratic terms
			previousEuler_dudx = average_dudx_Eulerian;	
	    }

	    // calculate finite granular strain
	    finiteStrain.clear();
	    average_dudx_Lagrangian.clear();	// used to test quadratic terms
	    finiteStrain = getFiniteStrain();
	    // calculate mixed finite granular strain
	    lag_HOT.clear();
	    lag_HOT = getHigherStrain();

	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
			<< setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2 << setw(OWID) << sigma2_1
			<< setw(OWID) << sigma2_2 << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << granularStress(1,1) << setw(OWID) << granularStress(1,2) << setw(OWID) << granularStress(1,3)
		        << setw(OWID) << granularStress(2,1) << setw(OWID) << granularStress(2,2) << setw(OWID) << granularStress(2,3)
		        << setw(OWID) << granularStress(3,1) << setw(OWID) << granularStress(3,2) << setw(OWID) << granularStress(3,3)
		        << setw(OWID) << getAverageRigidPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
			<< setw(OWID) << granularStrain(1,1) << setw(OWID) << granularStrain(1,2) << setw(OWID) << granularStrain(1,3)
		        << setw(OWID) << granularStrain(2,1) << setw(OWID) << granularStrain(2,2) << setw(OWID) << granularStrain(2,3)
		        << setw(OWID) << granularStrain(3,1) << setw(OWID) << granularStrain(3,2) << setw(OWID) << granularStrain(3,3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
			<< setw(OWID) << finiteStrain(1,1) << setw(OWID) << finiteStrain(1,2) << setw(OWID) << finiteStrain(1,3)
		        << setw(OWID) << finiteStrain(2,1) << setw(OWID) << finiteStrain(2,2) << setw(OWID) << finiteStrain(2,3)
		        << setw(OWID) << finiteStrain(3,1) << setw(OWID) << finiteStrain(3,2) << setw(OWID) << finiteStrain(3,3)
			<< setw(OWID) << lag_HOT(1,1) << setw(OWID) << lag_HOT(1,2) << setw(OWID) << lag_HOT(1,3)
		        << setw(OWID) << lag_HOT(2,1) << setw(OWID) << lag_HOT(2,2) << setw(OWID) << lag_HOT(2,3)
		        << setw(OWID) << lag_HOT(3,1) << setw(OWID) << lag_HOT(3,2) << setw(OWID) << lag_HOT(3,3)
			<< setw(OWID) << eulerianStrain(1,1) << setw(OWID) << eulerianStrain(1,2) << setw(OWID) << eulerianStrain(1,3)
		        << setw(OWID) << eulerianStrain(2,1) << setw(OWID) << eulerianStrain(2,2) << setw(OWID) << eulerianStrain(2,3)
		        << setw(OWID) << eulerianStrain(3,1) << setw(OWID) << eulerianStrain(3,2) << setw(OWID) << eulerianStrain(3,3)
			<< setw(OWID) << euler_HOT(1,1) << setw(OWID) << euler_HOT(1,2) << setw(OWID) << euler_HOT(1,3)
		        << setw(OWID) << euler_HOT(2,1) << setw(OWID) << euler_HOT(2,2) << setw(OWID) << euler_HOT(2,3)
		        << setw(OWID) << euler_HOT(3,1) << setw(OWID) << euler_HOT(3,2) << setw(OWID) << euler_HOT(3,3)
		        << setw(OWID) << epsilon_w_log
		        << setw(OWID) << epsilon_l_log
		        << setw(OWID) << epsilon_h_log
		        << setw(OWID) << (epsilon_w_log+epsilon_l_log+epsilon_h_log)
			<< setw(OWID) << average_dudx_Bagi(1,1) << setw(OWID) << average_dudx_Bagi(1,2) << setw(OWID) << average_dudx_Bagi(1,3)
		        << setw(OWID) << average_dudx_Bagi(2,1) << setw(OWID) << average_dudx_Bagi(2,2) << setw(OWID) << average_dudx_Bagi(2,3)
		        << setw(OWID) << average_dudx_Bagi(3,1) << setw(OWID) << average_dudx_Bagi(3,2) << setw(OWID) << average_dudx_Bagi(3,3)
		        
			<< setw(OWID) << average_dudx_Lagrangian(1,1) << setw(OWID) << average_dudx_Lagrangian(1,2) << setw(OWID) << average_dudx_Lagrangian(1,3)
		        << setw(OWID) << average_dudx_Lagrangian(2,1) << setw(OWID) << average_dudx_Lagrangian(2,2) << setw(OWID) << average_dudx_Lagrangian(2,3)
		        << setw(OWID) << average_dudx_Lagrangian(3,1) << setw(OWID) << average_dudx_Lagrangian(3,2) << setw(OWID) << average_dudx_Lagrangian(3,3)
			<< setw(OWID) << average_dudx_Eulerian(1,1) << setw(OWID) << average_dudx_Eulerian(1,2) << setw(OWID) << average_dudx_Eulerian(1,3)
		        << setw(OWID) << average_dudx_Eulerian(2,1) << setw(OWID) << average_dudx_Eulerian(2,2) << setw(OWID) << average_dudx_Eulerian(2,3)
		        << setw(OWID) << average_dudx_Eulerian(3,1) << setw(OWID) << average_dudx_Eulerian(3,2) << setw(OWID) << average_dudx_Eulerian(3,3) 
			<< endl;

	    g_debuginf << setw(OWID) << g_iteration
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
/*
	// below is controlled by granular stress
	if (	fabs(fabs(granularStress(1,1))-sigma)/sigma < STRESS_ERROR
	    &&  fabs(fabs(granularStress(2,2))-sigma)/sigma < STRESS_ERROR
	    &&  fabs(fabs(granularStress(3,3))-sigma)/sigma < STRESS_ERROR ) {
*/
	    balancedinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()    
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		        << endl;
	    sigma += sigma_inc;
	}

	// 9. loop break condition
	if (   fabs(sigma1_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma1_2-sigma_b)/sigma_b < STRESS_ERROR
	    && fabs(sigma2_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma2_2-sigma_b)/sigma_b < STRESS_ERROR
	    && fabs(sigma3_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma3_2-sigma_b)/sigma_b < STRESS_ERROR ) {
/*
	// below is controlled by granular stress
	if (	fabs(fabs(granularStress(1,1))-sigma_b)/sigma_b < STRESS_ERROR
	    &&  fabs(fabs(granularStress(2,2))-sigma_b)/sigma_b < STRESS_ERROR
	    &&  fabs(fabs(granularStress(3,3))-sigma_b)/sigma_b < STRESS_ERROR ) {
*/
	    // calculate granular strain
	    granularStrain.clear();
	    granularStrain = getGranularStrain()+previousStrain;
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()    
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2 << setw(OWID) << sigma2_1
			<< setw(OWID) << sigma2_2 << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << granularStress(1,1) << setw(OWID) << granularStress(1,2) << setw(OWID) << granularStress(1,3)
		        << setw(OWID) << granularStress(2,1) << setw(OWID) << granularStress(2,2) << setw(OWID) << granularStress(2,3)
		        << setw(OWID) << granularStress(3,1) << setw(OWID) << granularStress(3,2) << setw(OWID) << granularStress(3,3)
		        << setw(OWID) << getAverageRigidPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
			<< setw(OWID) << granularStrain(1,1) << setw(OWID) << granularStrain(1,2) << setw(OWID) << granularStrain(1,3)
		        << setw(OWID) << granularStrain(2,1) << setw(OWID) << granularStrain(2,2) << setw(OWID) << granularStrain(2,3)
		        << setw(OWID) << granularStrain(3,1) << setw(OWID) << granularStrain(3,2) << setw(OWID) << granularStrain(3,3)
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
			<< setw(OWID) << finiteStrain(1,1) << setw(OWID) << finiteStrain(1,2) << setw(OWID) << finiteStrain(1,3)
		        << setw(OWID) << finiteStrain(2,1) << setw(OWID) << finiteStrain(2,2) << setw(OWID) << finiteStrain(2,3)
		        << setw(OWID) << finiteStrain(3,1) << setw(OWID) << finiteStrain(3,2) << setw(OWID) << finiteStrain(3,3)
			<< setw(OWID) << lag_HOT(1,1) << setw(OWID) << lag_HOT(1,2) << setw(OWID) << lag_HOT(1,3)
		        << setw(OWID) << lag_HOT(2,1) << setw(OWID) << lag_HOT(2,2) << setw(OWID) << lag_HOT(2,3)
		        << setw(OWID) << lag_HOT(3,1) << setw(OWID) << lag_HOT(3,2) << setw(OWID) << lag_HOT(3,3)
			<< setw(OWID) << eulerianStrain(1,1) << setw(OWID) << eulerianStrain(1,2) << setw(OWID) << eulerianStrain(1,3)
		        << setw(OWID) << eulerianStrain(2,1) << setw(OWID) << eulerianStrain(2,2) << setw(OWID) << eulerianStrain(2,3)
		        << setw(OWID) << eulerianStrain(3,1) << setw(OWID) << eulerianStrain(3,2) << setw(OWID) << eulerianStrain(3,3)
			<< setw(OWID) << euler_HOT(1,1) << setw(OWID) << euler_HOT(1,2) << setw(OWID) << euler_HOT(1,3)
		        << setw(OWID) << euler_HOT(2,1) << setw(OWID) << euler_HOT(2,2) << setw(OWID) << euler_HOT(2,3)
		        << setw(OWID) << euler_HOT(3,1) << setw(OWID) << euler_HOT(3,2) << setw(OWID) << euler_HOT(3,3)
		        << setw(OWID) << epsilon_w_log
		        << setw(OWID) << epsilon_l_log
		        << setw(OWID) << epsilon_h_log
		        << setw(OWID) << (epsilon_w_log+epsilon_l_log+epsilon_h_log)
			<< setw(OWID) << average_dudx_Bagi(1,1) << setw(OWID) << average_dudx_Bagi(1,2) << setw(OWID) << average_dudx_Bagi(1,3)
		        << setw(OWID) << average_dudx_Bagi(2,1) << setw(OWID) << average_dudx_Bagi(2,2) << setw(OWID) << average_dudx_Bagi(2,3)
		        << setw(OWID) << average_dudx_Bagi(3,1) << setw(OWID) << average_dudx_Bagi(3,2) << setw(OWID) << average_dudx_Bagi(3,3)
		        
			<< setw(OWID) << average_dudx_Lagrangian(1,1) << setw(OWID) << average_dudx_Lagrangian(1,2) << setw(OWID) << average_dudx_Lagrangian(1,3)
		        << setw(OWID) << average_dudx_Lagrangian(2,1) << setw(OWID) << average_dudx_Lagrangian(2,2) << setw(OWID) << average_dudx_Lagrangian(2,3)
		        << setw(OWID) << average_dudx_Lagrangian(3,1) << setw(OWID) << average_dudx_Lagrangian(3,2) << setw(OWID) << average_dudx_Lagrangian(3,3)
			<< setw(OWID) << average_dudx_Eulerian(1,1) << setw(OWID) << average_dudx_Eulerian(1,2) << setw(OWID) << average_dudx_Eulerian(1,3)
		        << setw(OWID) << average_dudx_Eulerian(2,1) << setw(OWID) << average_dudx_Eulerian(2,2) << setw(OWID) << average_dudx_Eulerian(2,3)
		        << setw(OWID) << average_dudx_Eulerian(3,1) << setw(OWID) << average_dudx_Eulerian(3,2) << setw(OWID) << average_dudx_Eulerian(3,3)
		        << endl; 
	    break;
	}
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


// loading-unloading-reloading of isotropic compression
// the stress path is defined by sigma_points and sigma_values[]
void assembly::isotropic(int   total_steps,  
			 int   snapshots, 
			 int   interval,
			 int   sigma_points,  
			 REAL sigma_values[],  
			 int   sigma_division,	  
			 const char* iniptclfile,  
			 const char* inibdryfile,
			 const char* particlefile, 
			 const char* boundaryfile,
			 const char* contactfile,  
			 const char* progressfile,
			 const char* balancedfile, 
			 const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3: define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    initVolume = W0*L0*H0;
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    matrix granularStress;	// granular stress
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
    g_iteration=0;
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
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	// granular stress calculation
	granularStress.clear();
	granularStress = getGranularStress();
	void_ratio=Volume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	if (sigma2_1<sigma)
	    midctl[0].tran=vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma)
	    midctl[1].tran=vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma)
	    maxctl[0].tran=vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma)
	    maxctl[1].tran=vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=vec(0,-TIMESTEP*RELEASE_RATE,0);
	
	updateRB(min,minctl,2);
	updateRB(mid,midctl,2);
	updateRB(max,maxctl,2);
	updateRB6();
	
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0 ){
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
			<< setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2 << setw(OWID) << sigma2_1
			<< setw(OWID) << sigma2_2 << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << granularStress(1,1) << setw(OWID) << granularStress(1,2) << setw(OWID) << granularStress(1,3)
		        << setw(OWID) << granularStress(2,1) << setw(OWID) << granularStress(2,2) << setw(OWID) << granularStress(2,3)
		        << setw(OWID) << granularStress(3,1) << setw(OWID) << granularStress(3,2) << setw(OWID) << granularStress(3,3)
		        << setw(OWID) << getAverageRigidPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		        << endl;
	    g_debuginf << setw(OWID) << g_iteration
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
/*
	// below is granular stress control
	if (	fabs(fabs(granularStress(1,1))-sigma)/sigma < STRESS_ERROR
	    &&  fabs(fabs(granularStress(2,2))-sigma)/sigma < STRESS_ERROR
	    &&  fabs(fabs(granularStress(3,3))-sigma)/sigma < STRESS_ERROR ) {
*/
	    balancedinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()    
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
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
/*
	// below is granular stress control
	if (	fabs(fabs(granularStress(1,1))-sigma_b)/sigma_b < STRESS_ERROR
	    &&  fabs(fabs(granularStress(2,2))-sigma_b)/sigma_b < STRESS_ERROR
	    &&  fabs(fabs(granularStress(3,3))-sigma_b)/sigma_b < STRESS_ERROR ) {
*/
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()    
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
			<< setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2 << setw(OWID) << sigma2_1
			<< setw(OWID) << sigma2_2 << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << granularStress(1,1) << setw(OWID) << granularStress(1,2) << setw(OWID) << granularStress(1,3)
		        << setw(OWID) << granularStress(2,1) << setw(OWID) << granularStress(2,2) << setw(OWID) << granularStress(2,3)
		        << setw(OWID) << granularStress(3,1) << setw(OWID) << granularStress(3,2) << setw(OWID) << granularStress(3,3)
		        << setw(OWID) << getAverageRigidPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		        << endl;
	    break;
	}
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_3. This function
// increases vertical pressure step by step to sigma_1, thus making it possible to find out
// balanced status where top & bottom particle pressure equals major principle stress. 
// Side boundaries are fixed, top and bottom plates are force-controlled.
void assembly::odometer(int   total_steps,  
			int   snapshots, 
			int   interval,
			REAL sigma_3,     
			REAL sigma_1,    
			int   sigma_division,			  
			const char* iniptclfile,  
			const char* inibdryfile,
			const char* particlefile,
			const char* boundaryfile,
			const char* contactfile,  
			const char* progressfile,
			const char* balancedfile, 
			const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries
 
    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    initVolume = W0*L0*H0;
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
    g_iteration=0;
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
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/displacement based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=Volume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	updateRB(min,minctl,2);
	updateRB6();

	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0){
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		        << endl;
	}

	// 8. find balanced status of odometer compression
	if (fabs(sigma3_1-sigma)/sigma < STRESS_ERROR && fabs(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()    
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		        << endl;
	    sigma += sigma_inc;
	}

	// 9. loop break condition
	if (fabs(sigma3_1-sigma_1)/sigma_1 < STRESS_ERROR && fabs(sigma3_2-sigma_1)/sigma_1 < STRESS_ERROR) {
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()    
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		        << endl;
	    break;
	}
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_3. This function
// increases vertical pressure step by step to sigma_1, thus making it possible to find out
// balanced status where top & bottom particle pressure equals major principle stress. 
// Side boundaries are fixed, top and bottom plates are force-controlled. Unloading path is
// applied.
void assembly::odometer(int   total_steps,  
			int   snapshots, 
			int   interval,
			int   sigma_points,  
			REAL sigma_values[],  
			int   sigma_division,			  
			const char* iniptclfile,  
			const char* inibdryfile,
			const char* particlefile, 
			const char* boundaryfile,
			const char* contactfile,  
			const char* progressfile,
			const char* balancedfile, 
			const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries
 
    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    initVolume = W0*L0*H0;
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
    g_iteration=0;
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
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/displacement based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=Volume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	updateRB(min,minctl,2);
	updateRB6();

	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0){
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		        << endl;
	}

	// 8. find balanced status of odometer compression
	if (fabs(sigma3_1-sigma)/sigma < STRESS_ERROR && fabs(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()    
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
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
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()    
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		        << endl;
	    break;
	}
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


void assembly::unconfined(int   total_steps,  
			  int   snapshots,	
			  int   interval,
			  const char* iniptclfile,  
			  const char* inibdryfile,
			  const char* particlefile,
			  const char* contactfile,  
			  const char* progressfile,
			  const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.  
    progressinf.open(progressfile);
    if(!progressinf) { cout << "stream error!" << endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "unconfined..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample"
	        << "          sample          sample          sample" << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           p             width          length           "
	        << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
	        << "epsilon_v" << endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries
 
    // pre_3. define variables used in iterations
    REAL sigma3_1, sigma3_2;
    int    stepsnum=0;
    char   stepsstr[4];
    char   stepsfp[50];
    REAL avgNormal=0;
    REAL avgTangt=0;
    int    min[2]={5,6};    //  boundary 5 and 6
    UPDATECTL minctl[2];

    // iterations start here...
    g_iteration=0;
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
	rigidBoundaryForce();

	// 5. update particles' velocity/omga/displacement based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	sigma3_1=vfabs(getNormalForce(5))/getArea(5); sigma3_2=vfabs(getNormalForce(6))/getArea(6);
	minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	minctl[1].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);
	updateRB(min,minctl,2);

	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output progress info.
	if (g_iteration % interval == 0)
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << 0
		        << setw(OWID) << 0 << setw(OWID) << 0
		        << setw(OWID) << 0 << setw(OWID) << 0
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << endl;
/*
	// 8. loop break condition
	if (getAverageForce() < 1.0) {
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()    
		        << setw(OWID) << getAverageMoment() << endl;
	    break;
	}
*/
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


void assembly::iso_MemBdry(int   total_steps,  
			   int   snapshots, 
			   int   interval,
			   REAL  sigma3,
			   REAL  rRadius,
			   bool  toRebuild,
			   const char* iniptclfile, 
			   const char* particlefile,
			   const char* contactfile, 
			   const char* progressfile,
			   const char* debugfile) 
{
  // pre_1: open streams for output
  // particlefile and contactfile are used for snapshots at the end.
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
  
  g_debuginf.open(debugfile);
  if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
  g_debuginf.setf(std::ios::scientific, std::ios::floatfield);
  
  // pre_2. create particles from file and calculate forces caused by hydraulic pressure
  if (toRebuild) readSample(iniptclfile);
  REAL radius = gradInfo.getMinPtclRadius();
  if (gradInfo.getSize().size() == 1 &&
      gradInfo.getPtclRatioBA() == 1.0 && 
      gradInfo.getPtclRatioCA() == 1.0)
    radius *= rRadius; // determine how tiny the boundary particles are
  REAL mag  = radius*radius*4*sigma3;
  REAL x1 = container.getMinCorner().getx();
  REAL y1 = container.getMinCorner().gety();
  REAL z1 = container.getMinCorner().getz();
  REAL x2 = container.getMaxCorner().getx();
  REAL y2 = container.getMaxCorner().gety();
  REAL z2 = container.getMaxCorner().getz();
  std::vector<particle*>::const_iterator  it;
  vec pos;
  for (it=ParticleVec.begin();it!=ParticleVec.end();++it)
    {
      pos = (*it)->getCurrPosition();
      if (pos.getx() < x1)
	(*it)->setConstForce( vec(mag, 0, 0) );
      else if (pos.getx() > x2)
	(*it)->setConstForce( vec(-mag, 0, 0) );
      else if (pos.gety() < y1)
	(*it)->setConstForce( vec(0, mag, 0) );
      else if (pos.gety() > y2)
	(*it)->setConstForce( vec(0, -mag, 0) );
      else if (pos.getz() < z1)
	(*it)->setConstForce( vec(0, 0, mag) );
      else if (pos.getz() > z2)
	(*it)->setConstForce( vec(0, 0, -mag) );
    }

  // pre_3. define variables used in iterations
  REAL avgNormal=0;
  REAL avgTangt=0;
  int  stepsnum=0;
  char stepsstr[4];
  char stepsfp[50];
  
  // iterations start here...
  g_iteration=0;
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
      if (g_iteration % (total_steps/snapshots) == 0){
	sprintf(stepsstr, "%03d", stepsnum); 
	strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
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
	if (g_iteration % interval == 0 ){
	  progressinf << setw(OWID) << g_iteration
		      << setw(OWID) << getPossCntctNum()
		      << setw(OWID) << getActualCntctNum()
		      << setw(OWID) << getAveragePenetration()
		      << setw(OWID) << avgNormal
		      << setw(OWID) << avgTangt
		      << setw(OWID) << getAverageVelocity() 
		      << setw(OWID) << getAverageOmga()
		      << setw(OWID) << getAverageForce()   
		      << setw(OWID) << getAverageMoment()
		      << setw(OWID) << getTransEnergy()
		      << setw(OWID) << getRotatEnergy()
		      << setw(OWID) << getKinetEnergy()
		      << endl;
	}
	  
    } while (++g_iteration < total_steps);
  
  // post_1. store the final snapshot of particles, contacts and boundaries.
  strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
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
  g_debuginf.close();
}


// This function initializes triaxial sample to a certain confining pressure.
void assembly::triaxialPtclBdryIni(int   total_steps,  
				   int   snapshots, 
				   int   interval,
				   REAL  sigma,
				   const char* iniptclfile, 
				   const char* inibdryfile,
				   const char* particlefile,
				   const char* boundaryfile,
				   const char* contactfile, 
				   const char* progressfile,
				   const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
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
    g_iteration=0;
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
	rigidBoundaryForce();

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	sigma3_1=vfabs(getNormalForce(5))/2.5e-3; sigma3_2=vfabs(getNormalForce(6))/2.5e-3;

	// force control
	if (sigma3_1 < sigma)
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);

	if (sigma3_2 < sigma)
	    minctl[1].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);

	updateRB(min,minctl,2);
	
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	l56=getApt(5).getz()-getApt(6).getz();
	epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0 ){
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << 0 << setw(OWID) << 0
		        << setw(OWID) << 0 << setw(OWID) << 0
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()
		        << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << l56
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 2.0*getActualCntctNum()/TotalNum
		        << endl;

	}

	// 9. loop break condition: through displacement control mechanism
	if (   fabs(sigma3_1-sigma)/sigma < STRESS_ERROR && fabs(sigma3_2-sigma)/sigma < STRESS_ERROR )
	       break;
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


// This function performs triaxial compression test.
// Displacement boundaries are used in axial direction.
void assembly::triaxialPtclBdry(int   total_steps,  
				int   snapshots, 
				int   interval,
				const char* iniptclfile, 
				const char* inibdryfile,
				const char* particlefile,
				const char* boundaryfile,
				const char* contactfile, 
				const char* progressfile,
				const char* balancedfile,
				const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
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
    g_iteration=0;
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
	rigidBoundaryForce();

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	sigma3_1=vfabs(getNormalForce(5))/2.5e-3; sigma3_2=vfabs(getNormalForce(6))/2.5e-3;

	// displacement control
	if(g_iteration < 100001) {
	minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	minctl[1].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);

	updateRB(min,minctl,2);
	}
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	l56=getApt(5).getz()-getApt(6).getz();
	epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0 ){
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << 0 << setw(OWID) << 0
		        << setw(OWID) << 0 << setw(OWID) << 0
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()
		        << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << l56
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 2.0*getActualCntctNum()/TotalNum
		        << endl;

	}

/* Most time it is balanced, so use progressinf instead.
	// 8. find the balanced status and increase confining pressure
	if (   fabs(sigma1_1-sigma_a)/sigma_a < STRESS_ERROR && fabs(sigma1_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabs(sigma2_1-sigma_a)/sigma_a < STRESS_ERROR && fabs(sigma2_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabs(sigma3_1-sigma3_2)/(sigma3_1+sigma3_2)*2<=0.05) {
	    balancedinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()    
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h) << endl;
	}
*/
	// 9. loop break condition: through displacement control mechanism
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_a. This function
// performs triaxial compression test. Displacement boundaries are used in axial direction.
void assembly::triaxial(int   total_steps,  
			int   snapshots, 
			int   interval,
			REAL  sigma_a,	  
			const char* iniptclfile, 
			const char* inibdryfile,
			const char* particlefile,
			const char* boundaryfile,
			const char* contactfile, 
			const char* progressfile,
			const char* balancedfile,
			const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
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
		<< setw(OWID) << "sample" // for new 3 sigma's
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
		<< setw(OWID) << "sample"
	        << setw(OWID) << "void"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "coordinate"
	        << setw(OWID) << "vibra"
	        << setw(OWID) << "impact"
	        << setw(OWID) << "wall-clock" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag"
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "Logarithmic" 	// for logarithmic strain
		<< setw(OWID) << "Logarithmic"
		<< setw(OWID) << "Logarithmic"
		<< setw(OWID) << "Logarithmic" 
		<< setw(OWID) << "Bagi" 	// for average_dudx_Bagi, used to test quadratic terms, April 22, 2013
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Lagrangian" 	// for average_dudx_Lagrangian, used to test quadratic terms
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian"
		<< setw(OWID) << "Eulerian" 	// for average_dudx_Eulerian, used to test quadratic terms
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "spatial" 	// for spatial velocity gradient tensor, deformation rate tensor
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "from rate" 	// for strain based on spatial deformation rate tensor, July 8, 2013
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" << endl
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
	        << setw(OWID) << "sigma_11"
	        << setw(OWID) << "sigma_12"
	        << setw(OWID) << "sigma_13"
	        << setw(OWID) << "sigma_21"
	        << setw(OWID) << "sigma_22"
	        << setw(OWID) << "sigma_23"
	        << setw(OWID) << "sigma_31"
	        << setw(OWID) << "sigma_32"
	        << setw(OWID) << "sigma_33"
	        << setw(OWID) << "mean_stress"
	        << setw(OWID) << "width"
	        << setw(OWID) << "length"
	        << setw(OWID) << "height"
	        << setw(OWID) << "volume"
	        << setw(OWID) << "epsilon_w"
	        << setw(OWID) << "epsilon_l"
	        << setw(OWID) << "epsilon_h"
	        << setw(OWID) << "epsilon_v"
	        << setw(OWID) << "epsilon_11"
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "ratio"
	        << setw(OWID) << "porosity"
	        << setw(OWID) << "number"
	        << setw(OWID) << "t_step"
	        << setw(OWID) << "t_step"
	        << setw(OWID) << "time" 
	        << setw(OWID) << "epsilon_11"	// for finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "epsilon_11"	// for mixed finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "epsilon_11"	// for eulerian finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33" 
	        << setw(OWID) << "epsilon_11"	// for mixed eulerian finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33" 
		<< setw(OWID) << "epsilon_w"	// for logranthmic strain
	        << setw(OWID) << "epsilon_l"
	        << setw(OWID) << "epsilon_h"
	        << setw(OWID) << "epsilon_v" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Bagi, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Lagrangian, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Eulerian, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" 
	        << setw(OWID) << "dvdx_11"	// for spatial velocity gradient tensor, deformation rate tensor
	        << setw(OWID) << "dvdx_12"
	        << setw(OWID) << "dvdx_13"
	        << setw(OWID) << "dvdx_21"
	        << setw(OWID) << "dvdx_22"
	        << setw(OWID) << "dvdx_23"
	        << setw(OWID) << "dvdx_31"
	        << setw(OWID) << "dvdx_32"
	        << setw(OWID) << "dvdx_33"
	        << setw(OWID) << "epsilon_11"	// for strain based on spatial deformation rate tensor
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33" << endl;

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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    std::vector<REAL> temp_vec;	// used to test quadratic terms, April 22, 2013
    temp_vec.push_back(0);
    temp_vec.push_back(0);
    temp_vec.push_back(0);
    for(int i=0; i!=3; i++){
	average_dudx_Bagi.appendRow(temp_vec);
	average_dudx_Eulerian.appendRow(temp_vec);
	average_dudx_Lagrangian.appendRow(temp_vec);
    }

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    initVolume = W0*L0*H0;
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    matrix granularStress;
    matrix granularStrain;	// granular strain
    matrix previousStrain(3,3);
    matrix previousEuler(3,3);
    matrix finiteStrain(3,3);	// finite granular strain
    matrix lag_HOT(3,3);	// mixed finite granular strain
    matrix eulerianStrain(3,3);	// eulerian granular strain
    matrix euler_HOT(3,3);
    matrix previousEuler_HOT(3,3);
    matrix previousEuler_dudx(3,3);	// used to test quadratic terms, April 22, 2013
    matrix previousBagi_dudx(3,3);
    matrix spatial_dvdx(3,3);	// average spatial velocity gradient tensor

    matrix curr_rate(3,3);	// current spatial deformation rate tensor, July 8, 2013
    matrix curr_dvdx(3,3);	// current spatial velocity gradient tensor
    matrix prev_strain_rate(3,3);
    matrix curr_strain_rate(3,3);	// current strain based on deformation rate tensor


    // initialize previousStrain(3,3)
    for(int i_ps=0; i_ps!=2; i_ps++){
	for(int j_ps=0; j_ps!=2; j_ps++){
		previousStrain(i_ps+1,j_ps+1) = 0;
		previousEuler(i_ps+1,j_ps+1) = 0;
		previousEuler_HOT(i_ps+1,j_ps+1) = 0;
		previousEuler_dudx(i_ps+1,j_ps+1) = 0;	// used to test quadratic terms
		previousBagi_dudx(i_ps+1,j_ps+1) = 0;
		curr_strain_rate(i_ps+1,j_ps+1) = 0;	// calculate strain by rate, July 8, 2013
	}
    }
    
    // first tessellation
    createInputForQhull();
    callQhull();
    readTesse("tess_info");
    readTesse_finite("tess_info");
    cellVec_init = cellVec;	// keep cellVec_init for lagrangian strain unchanged 

    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL epsilon_w_log, epsilon_l_log, epsilon_h_log;	// lograthmic strain
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
    g_iteration=0;
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
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;

//	granularStress.clear();
//	granularStress = getGranularStress();	// calculate granular stress

	void_ratio=Volume/getParticleVolume()-1;

	// displacement control
	minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	minctl[1].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);

	// force control
	if (sigma2_1<sigma_a)
	    midctl[0].tran=vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma_a)
	    midctl[1].tran=vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma_a)
	    maxctl[0].tran=vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma_a)
	    maxctl[1].tran=vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=vec(0,-TIMESTEP*RELEASE_RATE,0);

/*
// below is adjustment of the velocities of boundary walls is based on new stress calculatio

	// along the x direction
	if (fabs(-sigma_11)<sigma_a){
		midctl[0].tran=vec(-TIMESTEP*COMPRESS_RATE,0,0);
		midctl[1].tran=vec(TIMESTEP*COMPRESS_RATE,0,0);	
	}
	else{
		midctl[0].tran=vec(TIMESTEP*RELEASE_RATE,0,0);
		midctl[1].tran=vec(-TIMESTEP*RELEASE_RATE,0,0);
	}
	// along the y direction
	if (fabs(-sigma_22)<sigma_a){
		maxctl[0].tran=vec(0,-TIMESTEP*COMPRESS_RATE,0);
		maxctl[1].tran=vec(0,TIMESTEP*COMPRESS_RATE,0);
	}
	else{
		maxctl[0].tran=vec(0,TIMESTEP*RELEASE_RATE,0);
		maxctl[1].tran=vec(0,-TIMESTEP*RELEASE_RATE,0);
	}

// above is modification about the adjustment of velocities
*/
	
	updateRB(min,minctl,2);
	updateRB(mid,midctl,2);
	updateRB(max,maxctl,2);
	updateRB6();
	
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
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
	epsilon_w_log = log(l24/W0); epsilon_l_log = log(l13/L0); epsilon_h_log = log(l56/H0);
	if (g_iteration % interval == 0 ){

//	    if(eulerianStrain(1,1)-previousEuler(1,1)>strainThreshold ||
//	       eulerianStrain(2,2)-previousEuler(2,2)>strainThreshold ||
//	       eulerianStrain(3,3)-previousEuler(3,3)>strainThreshold ){	// need to tessellate again
//			resetStartPosition();	// reset initial position
			// tessellate again
			createInputForQhull();
			callQhull();
			readTesse("tess_info");
			readTesse_finite("tess_info");
//	    }

	// calculate strain by rate
	prev_strain_rate.clear();
	prev_strain_rate = curr_strain_rate;	// set previous strain
	
	// calculate spatial velocity gradient tensor
	spatial_dvdx.clear();
	spatial_dvdx = getAverage_dvdx();
	curr_dvdx.clear();
	curr_dvdx = spatial_dvdx;
	
	curr_rate = 0.5*(curr_dvdx+curr_dvdx.getTrans());	// current spatial deformation rate

	curr_strain_rate = prev_strain_rate+(curr_rate-curr_dvdx.getTrans()*prev_strain_rate-prev_strain_rate*curr_dvdx)*TIMESTEP*interval;


	    // calculate granular strain and Euler strain
	    average_dudx_Bagi.clear();	// used to test quadratic terms
	    average_dudx_Eulerian.clear();

	    granularStrain.clear();
	    eulerianStrain.clear();
	    euler_HOT.clear();
//	    granularStrain = getGranularStrain()+previousStrain;
//	    eulerianStrain = getEulerianStrain()+previousEuler;
//	    euler_HOT = getEuler_HOT()+previousEuler_HOT;
//
//	    average_dudx_Bagi = average_dudx_Bagi+previousBagi_dudx;	// used to test quadratic terms
//	    average_dudx_Eulerian = average_dudx_Eulerian+previousEuler_dudx;

	    granularStrain = getGranularStrain();
	    eulerianStrain = getEulerianStrain();
	    euler_HOT = getEuler_HOT();

	    average_dudx_Bagi = average_dudx_Bagi;	// used to test quadratic terms
	    average_dudx_Eulerian = average_dudx_Eulerian;

	    // update previousStrain
	    previousStrain = granularStrain;
	    previousEuler = eulerianStrain;
	    previousEuler_HOT = euler_HOT;
	    previousBagi_dudx = average_dudx_Bagi;	// used to test quadratic terms
	    previousEuler_dudx = average_dudx_Eulerian;

	    // calculate finite granular strain
	    finiteStrain.clear();
	    average_dudx_Lagrangian.clear();	// used to test quadratic terms
	    finiteStrain = getFiniteStrain();
	    // calculate mixed finite granular strain
	    lag_HOT.clear();
	    lag_HOT = getHigherStrain();


	    granularStress.clear();
	    granularStress = getGranularStress();	// calculate granular stress

	    gettimeofday(&time_w2,NULL);
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2 << setw(OWID) << sigma2_1
			<< setw(OWID) << sigma2_2 << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << granularStress(1,1) << setw(OWID) << granularStress(1,2) << setw(OWID) << granularStress(1,3)
		        << setw(OWID) << granularStress(2,1) << setw(OWID) << granularStress(2,2) << setw(OWID) << granularStress(2,3)
		        << setw(OWID) << granularStress(3,1) << setw(OWID) << granularStress(3,2) << setw(OWID) << granularStress(3,3)
		        << setw(OWID) << getAverageRigidPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
			<< setw(OWID) << granularStrain(1,1) << setw(OWID) << granularStrain(1,2) << setw(OWID) << granularStrain(1,3)
		        << setw(OWID) << granularStrain(2,1) << setw(OWID) << granularStrain(2,2) << setw(OWID) << granularStrain(2,3)
		        << setw(OWID) << granularStrain(3,1) << setw(OWID) << granularStrain(3,2) << setw(OWID) << granularStrain(3,3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
	                << setw(OWID) << getVibraTimeStep()
	                << setw(OWID) << getImpactTimeStep()
		        << setw(OWID) << timediffsec(time_w1,time_w2)
			<< setw(OWID) << finiteStrain(1,1) << setw(OWID) << finiteStrain(1,2) << setw(OWID) << finiteStrain(1,3)
		        << setw(OWID) << finiteStrain(2,1) << setw(OWID) << finiteStrain(2,2) << setw(OWID) << finiteStrain(2,3)
		        << setw(OWID) << finiteStrain(3,1) << setw(OWID) << finiteStrain(3,2) << setw(OWID) << finiteStrain(3,3)
			<< setw(OWID) << lag_HOT(1,1) << setw(OWID) << lag_HOT(1,2) << setw(OWID) << lag_HOT(1,3)
		        << setw(OWID) << lag_HOT(2,1) << setw(OWID) << lag_HOT(2,2) << setw(OWID) << lag_HOT(2,3)
		        << setw(OWID) << lag_HOT(3,1) << setw(OWID) << lag_HOT(3,2) << setw(OWID) << lag_HOT(3,3)
			<< setw(OWID) << eulerianStrain(1,1) << setw(OWID) << eulerianStrain(1,2) << setw(OWID) << eulerianStrain(1,3)
		        << setw(OWID) << eulerianStrain(2,1) << setw(OWID) << eulerianStrain(2,2) << setw(OWID) << eulerianStrain(2,3)
		        << setw(OWID) << eulerianStrain(3,1) << setw(OWID) << eulerianStrain(3,2) << setw(OWID) << eulerianStrain(3,3)
			<< setw(OWID) << euler_HOT(1,1) << setw(OWID) << euler_HOT(1,2) << setw(OWID) << euler_HOT(1,3)
		        << setw(OWID) << euler_HOT(2,1) << setw(OWID) << euler_HOT(2,2) << setw(OWID) << euler_HOT(2,3)
		        << setw(OWID) << euler_HOT(3,1) << setw(OWID) << euler_HOT(3,2) << setw(OWID) << euler_HOT(3,3)
		        << setw(OWID) << epsilon_w_log
		        << setw(OWID) << epsilon_l_log
		        << setw(OWID) << epsilon_h_log
		        << setw(OWID) << (epsilon_w_log+epsilon_l_log+epsilon_h_log)
			<< setw(OWID) << average_dudx_Bagi(1,1) << setw(OWID) << average_dudx_Bagi(1,2) << setw(OWID) << average_dudx_Bagi(1,3)
		        << setw(OWID) << average_dudx_Bagi(2,1) << setw(OWID) << average_dudx_Bagi(2,2) << setw(OWID) << average_dudx_Bagi(2,3)
		        << setw(OWID) << average_dudx_Bagi(3,1) << setw(OWID) << average_dudx_Bagi(3,2) << setw(OWID) << average_dudx_Bagi(3,3)
		        
			<< setw(OWID) << average_dudx_Lagrangian(1,1) << setw(OWID) << average_dudx_Lagrangian(1,2) << setw(OWID) << average_dudx_Lagrangian(1,3)
		        << setw(OWID) << average_dudx_Lagrangian(2,1) << setw(OWID) << average_dudx_Lagrangian(2,2) << setw(OWID) << average_dudx_Lagrangian(2,3)
		        << setw(OWID) << average_dudx_Lagrangian(3,1) << setw(OWID) << average_dudx_Lagrangian(3,2) << setw(OWID) << average_dudx_Lagrangian(3,3)
			<< setw(OWID) << average_dudx_Eulerian(1,1) << setw(OWID) << average_dudx_Eulerian(1,2) << setw(OWID) << average_dudx_Eulerian(1,3)
		        << setw(OWID) << average_dudx_Eulerian(2,1) << setw(OWID) << average_dudx_Eulerian(2,2) << setw(OWID) << average_dudx_Eulerian(2,3)
		        << setw(OWID) << average_dudx_Eulerian(3,1) << setw(OWID) << average_dudx_Eulerian(3,2) << setw(OWID) << average_dudx_Eulerian(3,3)
			<< setw(OWID) << spatial_dvdx(1,1) << setw(OWID) << spatial_dvdx(1,2) << setw(OWID) << spatial_dvdx(1,3)
		        << setw(OWID) << spatial_dvdx(2,1) << setw(OWID) << spatial_dvdx(2,2) << setw(OWID) << spatial_dvdx(2,3)
		        << setw(OWID) << spatial_dvdx(3,1) << setw(OWID) << spatial_dvdx(3,2) << setw(OWID) << spatial_dvdx(3,3)
			<< setw(OWID) << curr_strain_rate(1,1) << setw(OWID) << curr_strain_rate(1,2) << setw(OWID) << curr_strain_rate(1,3)
		        << setw(OWID) << curr_strain_rate(2,1) << setw(OWID) << curr_strain_rate(2,2) << setw(OWID) << curr_strain_rate(2,3)
		        << setw(OWID) << curr_strain_rate(3,1) << setw(OWID) << curr_strain_rate(3,2) << setw(OWID) << curr_strain_rate(3,3)
		        << endl;
//	    g_debuginf << setw(OWID) << g_iteration
//		       << setw(OWID) << getTransEnergy()
//		       << setw(OWID) << getRotatEnergy()
//		       << setw(OWID) << bdry_penetr[1]
//		       << setw(OWID) << bdry_penetr[2]
//		       << setw(OWID) << bdry_penetr[3]
//		       << setw(OWID) << bdry_penetr[4]
//		       << setw(OWID) << bdry_penetr[5]
//		       << setw(OWID) << bdry_penetr[6]
//		       << setw(OWID) << bdry_cntnum[1]
//		       << setw(OWID) << bdry_cntnum[2]
//		       << setw(OWID) << bdry_cntnum[3]
//		       << setw(OWID) << bdry_cntnum[4]
//		       << setw(OWID) << bdry_cntnum[5]
//		       << setw(OWID) << bdry_cntnum[6]
//		       << endl;
	}

	// Most time it is balanced, so use progressinf instead.
	// 8. find the balanced status and increase confining pressure
/*	if (   fabs(sigma1_1-sigma_a)/sigma_a < STRESS_ERROR && fabs(sigma1_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabs(sigma2_1-sigma_a)/sigma_a < STRESS_ERROR && fabs(sigma2_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabs(sigma3_1-sigma3_2)/(sigma3_1+sigma3_2)*2<=0.05) {
*/
// below is based on new stress calculation
	if (	fabs(fabs(granularStress(1,1))-sigma_a)/sigma_a < STRESS_ERROR
	    &&  fabs(fabs(granularStress(2,2))-sigma_a)/sigma_a < STRESS_ERROR
	    &&  fabs(fabs(granularStress(3,3))-sigma_a)/sigma_a < STRESS_ERROR ) {
// above is based on new stress calculation
	    balancedinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()    
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
	                << setw(OWID) << getVibraTimeStep()
	                << setw(OWID) << getImpactTimeStep()
		        << setw(OWID) << timediffsec(time_w1,time_w2)
		        << endl;
	}

	// 9. loop break condition: through displacement control mechanism
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    g_timeinf << setw(4) << "end" << " " << ctime(&timeStamp) << flush;

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_a. This function
// performs triaxial compression test with unloading. Displacement boundaries are used in 
// axial direction.
void assembly::triaxial(int   total_steps,  
			int   unload_step,
			int   snapshots, 
			int   interval,
			REAL sigma_a,	  
			const char* iniptclfile,  
			const char* inibdryfile,
			const char* particlefile,
			const char* boundaryfile,
			const char* contactfile,
			const char* progressfile,
			const char* balancedfile,
			const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout << "stream error!" << endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf << "triaxial..." << endl
	        << "     iteration possible  actual      average	    average         average         average"
	        << "         average         average         average        sample            sample     "
	        << "     sample          sample          sample          sample          sample          "
	        << "sample          sample         sample           sample          sample          sample     "
		<< "sample          sample         sample           sample          sample          sample          "
		<< "sample          sample         sample           "
	        << "     sample          sample          sample          void            sample        coordinate"
	        << endl
	        << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	        << "          omga            force           moment        density          "
	        << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	        << "sigma3_1        sigma3_2           "
		<< "sigma_11        sigma_12        sigma_13        "
		<< "sigma_21        sigma_22        sigma_23        "
		<< "sigma_31        sigma_32        sigma_33        "
		<< "p             width          length           "
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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    initVolume = W0*L0*H0;
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    matrix granularStress;	// granular stress
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
    g_iteration=0;
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
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	// granular stress calculation
	granularStress.clear();
	granularStress = getGranularStress();
	void_ratio=Volume/getParticleVolume()-1;

	// displacement control
	if (g_iteration <= unload_step){ //loading
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	    minctl[1].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);
	}
	else { 
	    if (reload==false) { // unloading
		if (fabs(sigma3_1-sigma_a)/sigma_a > STRESS_ERROR && 
		    fabs(sigma3_2-sigma_a)/sigma_a > STRESS_ERROR){
		    minctl[0].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);
		    minctl[1].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
		}
		else  // reloading
		    reload=true;
	    }
	    else {
		minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
		minctl[1].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);
	    }
	}
	
	// force control
	if (sigma2_1<sigma_a)
	    midctl[0].tran=vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma_a)
	    midctl[1].tran=vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma_a)
	    maxctl[0].tran=vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma_a)
	    maxctl[1].tran=vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=vec(0,-TIMESTEP*RELEASE_RATE,0);
	
	updateRB(min,minctl,2);
	updateRB(mid,midctl,2);
	updateRB(max,maxctl,2);
	updateRB6();
	
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0 ){
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
			<< setw(OWID) << granularStress(1,1) << setw(OWID) << granularStress(1,2) << setw(OWID) << granularStress(1,3)
		        << setw(OWID) << granularStress(2,1) << setw(OWID) << granularStress(2,2) << setw(OWID) << granularStress(2,3)
		        << setw(OWID) << granularStress(3,1) << setw(OWID) << granularStress(3,2) << setw(OWID) << granularStress(3,3)
		        << setw(OWID) << getAverageRigidPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		        << endl;
	    g_debuginf << setw(OWID) << g_iteration
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

/*
	// 8. find the balanced status and increase confining pressure
	if (   fabs(sigma1_1-sigma_a)/sigma_a < STRESS_ERROR && fabs(sigma1_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabs(sigma2_1-sigma_a)/sigma_a < STRESS_ERROR && fabs(sigma2_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabs(sigma3_1-sigma3_2)/(sigma3_1+sigma3_2)*2<=0.05) {
	    balancedinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()    
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h) << endl;
	}
*/
	// 9. loop break condition: through displacement control mechanism
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


// The specimen has been deposited with gravitation within boundaries composed of particles.
// A rectangular pile is then drived into the particles using displacement control.
void assembly::rectPile_Disp(int   total_steps,  
			     int   snapshots, 
			     int   interval,
			     const char* iniptclfile,  
			     const char* inibdryfile,
			     const char* particlefile, 
			     const char* boundaryfile,
			     const char* contactfile,  
			     const char* progressfile,
			     const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);
    g_debuginf << " iteration    end_bearing     side_friction   total_force" << endl;

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
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
    g_iteration=0;
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
	rigidBoundaryForce();

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation

	// displacement control of the pile
	pilectl[0].tran=vec(0,0,-TIMESTEP*PILE_RATE);
	pilectl[1].tran=vec(0,0,-TIMESTEP*PILE_RATE);

	updateRB(pile, pilectl, 2); 
	updateRectPile();
	if (g_iteration % interval == 0) {
	    REAL  f7=getShearForce( 7).getz();
	    REAL  f8=getShearForce( 8).getz();
	    REAL  f9=getShearForce( 9).getz();
	    REAL f10=getShearForce(10).getz();
	    REAL  fn=getNormalForce(12).getz();
	    g_debuginf << setw(OWID) << g_iteration
		       << setw(OWID) << fn
		       << setw(OWID) << (f7+f8+f9+f10)
		       << setw(OWID) << (fn+f7+f8+f9+f10)
		       << endl;
	}

	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    printRectPile(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output statistics info.
	if (g_iteration % interval == 0) {
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3) << endl;
	}

	// 8. loop break condition
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);
    printRectPile(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


// The specimen has been deposited with gravitation within boundaries composed of particles.
// An ellipsoidal pile is then drived into the particles using displacement control.
void assembly::ellipPile_Disp(int   total_steps,  
			      int   snapshots, 
			      int   interval,
			      REAL dimn,
			      REAL rsize,
			      const char* iniptclfile,
			      const char* particlefile, 
			      const char* contactfile,  
			      const char* progressfile,
			      const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
		<< setw(OWID) << "sample" // for new 3 sigma's
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
		<< setw(OWID) << "sample"
	        << setw(OWID) << "void"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "coordinate"
	        << setw(OWID) << "vibra"
	        << setw(OWID) << "impact"
	        << setw(OWID) << "wall-clock" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag"
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "Logarithmic" 	// for logarithmic strain
		<< setw(OWID) << "Logarithmic"
		<< setw(OWID) << "Logarithmic"
		<< setw(OWID) << "Logarithmic" 
		<< setw(OWID) << "Bagi" 	// for average_dudx_Bagi, used to test quadratic terms, April 22, 2013
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Lagrangian" 	// for average_dudx_Lagrangian, used to test quadratic terms
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian"
		<< setw(OWID) << "Eulerian" 	// for average_dudx_Eulerian, used to test quadratic terms
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" << endl
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
	        << setw(OWID) << "sigma_11"
	        << setw(OWID) << "sigma_12"
	        << setw(OWID) << "sigma_13"
	        << setw(OWID) << "sigma_21"
	        << setw(OWID) << "sigma_22"
	        << setw(OWID) << "sigma_23"
	        << setw(OWID) << "sigma_31"
	        << setw(OWID) << "sigma_32"
	        << setw(OWID) << "sigma_33"
	        << setw(OWID) << "mean_stress"
	        << setw(OWID) << "width"
	        << setw(OWID) << "length"
	        << setw(OWID) << "height"
	        << setw(OWID) << "volume"
	        << setw(OWID) << "epsilon_w"
	        << setw(OWID) << "epsilon_l"
	        << setw(OWID) << "epsilon_h"
	        << setw(OWID) << "epsilon_v"
	        << setw(OWID) << "epsilon_11"
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "ratio"
	        << setw(OWID) << "porosity"
	        << setw(OWID) << "number"
	        << setw(OWID) << "t_step"
	        << setw(OWID) << "t_step"
	        << setw(OWID) << "time" 
	        << setw(OWID) << "epsilon_11"	// for finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "epsilon_11"	// for mixed finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "epsilon_11"	// for eulerian finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33" 
	        << setw(OWID) << "epsilon_11"	// for mixed eulerian finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33" 
		<< setw(OWID) << "epsilon_w"	// for logranthmic strain
	        << setw(OWID) << "epsilon_l"
	        << setw(OWID) << "epsilon_h"
	        << setw(OWID) << "epsilon_v" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Bagi, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Lagrangian, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Eulerian, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" << endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 

    // pre_3. define variables used in iterations
    REAL l13, l24, l56;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    matrix granularStress;
    matrix granularStrain;	// granular strain
    matrix previousStrain(3,3);
    matrix previousEuler(3,3);
    matrix finiteStrain(3,3);	// finite granular strain
    matrix lag_HOT(3,3);	// mixed finite granular strain
    matrix eulerianStrain(3,3);	// eulerian granular strain
    matrix euler_HOT(3,3);
    matrix previousEuler_HOT(3,3);
    matrix previousEuler_dudx(3,3);	// used to test quadratic terms, April 22, 2013
    matrix previousBagi_dudx(3,3);

    // initialize previousStrain(3,3)
    for(int i_ps=0; i_ps!=2; i_ps++){
	for(int j_ps=0; j_ps!=2; j_ps++){
		previousStrain(i_ps+1,j_ps+1) = 0;
		previousEuler(i_ps+1,j_ps+1) = 0;
		previousEuler_HOT(i_ps+1,j_ps+1) = 0;
		previousEuler_dudx(i_ps+1,j_ps+1) = 0;	// used to test quadratic terms
		previousBagi_dudx(i_ps+1,j_ps+1) = 0;
	}
    }
    
    // first tessellation
    createInputForQhull();
    callQhull();
    readTesse("tess_info");
    readTesse_finite("tess_info");
    cellVec_init = cellVec;	// keep cellVec_init for lagrangian strain unchanged 

    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    REAL void_ratio=0;
    
    // iterations start here...
    g_iteration=0;
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
	Volume=l13*l24*l56-ellipPilePeneVol();
	void_ratio=Volume/getParticleVolume()-1;

	granularStress.clear();
	granularStress = getGranularStress();	// calculate granular stress

	// 6. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 6. (2) output statistics info.
	if (g_iteration % interval == 0) {

	    // calculate granular strain and Euler strain
	    average_dudx_Bagi.clear();	// used to test quadratic terms
	    average_dudx_Eulerian.clear();

	    granularStrain.clear();
	    eulerianStrain.clear();
	    euler_HOT.clear();
	    granularStrain = getGranularStrain()+previousStrain;
	    eulerianStrain = getEulerianStrain()+previousEuler;
	    euler_HOT = getEuler_HOT()+previousEuler_HOT;

	    average_dudx_Bagi = average_dudx_Bagi+previousBagi_dudx;	// used to test quadratic terms
	    average_dudx_Eulerian = average_dudx_Eulerian+previousEuler_dudx;

	    if(eulerianStrain(1,1)-previousEuler(1,1)>strainThreshold ||
	       eulerianStrain(2,2)-previousEuler(2,2)>strainThreshold ||
	       eulerianStrain(3,3)-previousEuler(3,3)>strainThreshold ){	// need to tessellate again
			resetStartPosition();	// reset initial position
			// tessellate again
			createInputForQhull();
			callQhull();
			readTesse("tess_info");
			readTesse_finite("tess_info");
			// update previousStrain
			previousStrain = granularStrain;
			previousEuler = eulerianStrain;
			previousEuler_HOT = euler_HOT;
			previousBagi_dudx = average_dudx_Bagi;	// used to test quadratic terms
			previousEuler_dudx = average_dudx_Eulerian;
	    }
	    // calculate finite granular strain
	    finiteStrain.clear();
	    average_dudx_Lagrangian.clear();	// used to test quadratic terms
	    finiteStrain = getFiniteStrain();
	    // calculate mixed finite granular strain
	    lag_HOT.clear();
	    lag_HOT = getHigherStrain();

	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);

	    gettimeofday(&time_w2,NULL);
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << 0
		        << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 0
			<< setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 0
		        << setw(OWID) << granularStress(1,1) << setw(OWID) << granularStress(1,2) << setw(OWID) << granularStress(1,3)
		        << setw(OWID) << granularStress(2,1) << setw(OWID) << granularStress(2,2) << setw(OWID) << granularStress(2,3)
		        << setw(OWID) << granularStress(3,1) << setw(OWID) << granularStress(3,2) << setw(OWID) << granularStress(3,3)
		        << setw(OWID) << getAverageRigidPressure()
		        << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 0
			<< setw(OWID) << granularStrain(1,1) << setw(OWID) << granularStrain(1,2) << setw(OWID) << granularStrain(1,3)
		        << setw(OWID) << granularStrain(2,1) << setw(OWID) << granularStrain(2,2) << setw(OWID) << granularStrain(2,3)
		        << setw(OWID) << granularStrain(3,1) << setw(OWID) << granularStrain(3,2) << setw(OWID) << granularStrain(3,3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 0
	                << setw(OWID) << getVibraTimeStep()
	                << setw(OWID) << getImpactTimeStep()
		        << setw(OWID) << timediffsec(time_w1,time_w2)
			<< setw(OWID) << finiteStrain(1,1) << setw(OWID) << finiteStrain(1,2) << setw(OWID) << finiteStrain(1,3)
		        << setw(OWID) << finiteStrain(2,1) << setw(OWID) << finiteStrain(2,2) << setw(OWID) << finiteStrain(2,3)
		        << setw(OWID) << finiteStrain(3,1) << setw(OWID) << finiteStrain(3,2) << setw(OWID) << finiteStrain(3,3)
			<< setw(OWID) << lag_HOT(1,1) << setw(OWID) << lag_HOT(1,2) << setw(OWID) << lag_HOT(1,3)
		        << setw(OWID) << lag_HOT(2,1) << setw(OWID) << lag_HOT(2,2) << setw(OWID) << lag_HOT(2,3)
		        << setw(OWID) << lag_HOT(3,1) << setw(OWID) << lag_HOT(3,2) << setw(OWID) << lag_HOT(3,3)
			<< setw(OWID) << eulerianStrain(1,1) << setw(OWID) << eulerianStrain(1,2) << setw(OWID) << eulerianStrain(1,3)
		        << setw(OWID) << eulerianStrain(2,1) << setw(OWID) << eulerianStrain(2,2) << setw(OWID) << eulerianStrain(2,3)
		        << setw(OWID) << eulerianStrain(3,1) << setw(OWID) << eulerianStrain(3,2) << setw(OWID) << eulerianStrain(3,3)
			<< setw(OWID) << euler_HOT(1,1) << setw(OWID) << euler_HOT(1,2) << setw(OWID) << euler_HOT(1,3)
		        << setw(OWID) << euler_HOT(2,1) << setw(OWID) << euler_HOT(2,2) << setw(OWID) << euler_HOT(2,3)
		        << setw(OWID) << euler_HOT(3,1) << setw(OWID) << euler_HOT(3,2) << setw(OWID) << euler_HOT(3,3)
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 0
			<< setw(OWID) << average_dudx_Bagi(1,1) << setw(OWID) << average_dudx_Bagi(1,2) << setw(OWID) << average_dudx_Bagi(1,3)
		        << setw(OWID) << average_dudx_Bagi(2,1) << setw(OWID) << average_dudx_Bagi(2,2) << setw(OWID) << average_dudx_Bagi(2,3)
		        << setw(OWID) << average_dudx_Bagi(3,1) << setw(OWID) << average_dudx_Bagi(3,2) << setw(OWID) << average_dudx_Bagi(3,3)
		        
			<< setw(OWID) << average_dudx_Lagrangian(1,1) << setw(OWID) << average_dudx_Lagrangian(1,2) << setw(OWID) << average_dudx_Lagrangian(1,3)
		        << setw(OWID) << average_dudx_Lagrangian(2,1) << setw(OWID) << average_dudx_Lagrangian(2,2) << setw(OWID) << average_dudx_Lagrangian(2,3)
		        << setw(OWID) << average_dudx_Lagrangian(3,1) << setw(OWID) << average_dudx_Lagrangian(3,2) << setw(OWID) << average_dudx_Lagrangian(3,3)
			<< setw(OWID) << average_dudx_Eulerian(1,1) << setw(OWID) << average_dudx_Eulerian(1,2) << setw(OWID) << average_dudx_Eulerian(1,3)
		        << setw(OWID) << average_dudx_Eulerian(2,1) << setw(OWID) << average_dudx_Eulerian(2,2) << setw(OWID) << average_dudx_Eulerian(2,3)
		        << setw(OWID) << average_dudx_Eulerian(3,1) << setw(OWID) << average_dudx_Eulerian(3,2) << setw(OWID) << average_dudx_Eulerian(3,3) << endl;
	    g_debuginf << setw(OWID) << g_iteration
		       << setw(OWID) << getTopFreeParticlePosition().getz()
		       << setw(OWID) << ellipPileTipZ()
		       << setw(OWID) << getTopFreeParticlePosition().getz()-ellipPileTipZ()
		       << setw(OWID) << l13*l24*l56
		       << setw(OWID) << ellipPilePeneVol()
		       << setw(OWID) << Volume
		       << endl;
	}

	// 7. loop break condition
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


// The specimen has been deposited with gravitation within rigid boundaries.
// An ellipsoidal penetrator is then impacted into the particles with initial velocity.
void assembly::ellipPile_Impact(int   total_steps,  
				int   snapshots, 
				int   interval,
				REAL dimn,
				const char* iniptclfile,
				const char* inibdryfile,
				const char* particlefile, 
				const char* contactfile,  
				const char* progressfile,
				const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles
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
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);

	// 4. calculate boundary forces/moments and apply them to particles.
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - (-dimn/2);
	l24=dimn;
	l13=dimn;
	Volume=l13*l24*l56-ellipPilePeneVol();
	void_ratio=Volume/getParticleVolume()-1;

	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output statistics info.
	if (g_iteration % interval == 0) {
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/TotalNum
		        << endl;
	    g_debuginf << setw(OWID) << g_iteration
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
	    /*
	    g_debuginf << setw(OWID) << g_iteration
		       << setw(OWID) << getTopFreeParticlePosition().getz()
		       << setw(OWID) << ellipPileTipZ()
		       << setw(OWID) << getTopFreeParticlePosition().getz()-ellipPileTipZ()
		       << setw(OWID) << l13*l24*l56
		       << setw(OWID) << ellipPilePeneVol()
		       << setw(OWID) << Volume
		       << endl;
	    */
	}

	// 8. loop break condition
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


// The specimen has been deposited with gravitation within particle boundaries.
// An ellipsoidal penetrator is then impacted into the particles with initial velocity.
void assembly::ellipPile_Impact_p(int   total_steps,  
				  int   snapshots, 
				  int   interval,
				  REAL dimn,
				  const char* iniptclfile,
				  const char* particlefile, 
				  const char* contactfile,  
				  const char* progressfile,
				  const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles

    // pre_3. define variables used in iterations
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    REAL void_ratio=0;
    
    // iterations start here...
    g_iteration=0;
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
	Volume=l13*l24*l56-ellipPilePeneVol();
	void_ratio=Volume/getParticleVolume()-1;

	// 6. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 6. (2) output statistics info.
	if (g_iteration % interval == 0) {
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*getActualCntctNum()/TotalNum
		        << endl;
	    g_debuginf << setw(OWID) << g_iteration
		       << setw(OWID) << getTopFreeParticlePosition().getz()
		       << setw(OWID) << ellipPileTipZ()
		       << setw(OWID) << getTopFreeParticlePosition().getz()-ellipPileTipZ()
		       << setw(OWID) << l13*l24*l56
		       << setw(OWID) << ellipPilePeneVol()
		       << setw(OWID) << Volume
		       << endl;
	}

	// 7. loop break condition
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}



// The specimen has been deposited with gravitation within boundaries composed of particles.
// An ellipsoidal pile is then drived into the particles using force control.
// Not recommended.
void assembly::ellipPile_Force(int   total_steps,  
			       int   snapshots,
			       int   interval,
			       REAL dimn,
			       REAL force,
			       int   division,
			       const char* iniptclfile,
			       const char* particlefile, 
			       const char* contactfile,  
			       const char* progressfile,
			       const char* balancedfile,
			       const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 

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
    g_iteration=0;
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
	Volume=l13*l24*l56-ellipPilePeneVol();
	void_ratio=Volume/getParticleVolume()-1;
	
	// 6. update pile external force and position
	if(zforce>ellipPileForce())
	    ellipPileUpdate();

	if(fabs(ellipPileForce()-zforce)/zforce < STRESS_ERROR ){
	    balancedinf << setw(OWID) << g_iteration
		        << setw(OWID) << zforce
		        << setw(OWID) << getTopFreeParticlePosition().getz()-ellipPileTipZ()
		        << setw(OWID) << ellipPileForce()
		        << endl;
	    zforce += zforce_inc;
	}

	if( g_iteration % interval == 0){
	    g_debuginf << setw(OWID) << g_iteration
		       << setw(OWID) << zforce
		       << setw(OWID) << getTopFreeParticlePosition().getz()-ellipPileTipZ()
		       << setw(OWID) << ellipPileForce()
		       << endl;
	}

	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output statistics info.
	if (g_iteration % interval == 0) {
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << t1
		        << setw(OWID) << t2
		        << setw(OWID) << (t1+t2)
		        << setw(OWID) << t3
		        << setw(OWID) << (t1+t2+t3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*getActualCntctNum()/TotalNum
		        << endl;
	}

	// 8. loop break condition
	if (fabs((zforce-force)/force)<0.001)
	    break;
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_a. This function
// performs true triaxial test. Force boundaries are used.
void assembly::truetriaxial(int   total_steps,  
			    int   snapshots, 
			    int   interval,
			    REAL sigma_a,     
			    REAL sigma_w,
			    REAL sigma_l,     
			    REAL sigma_h,   
			    int   sigma_division,
			    const char* iniptclfile,  
			    const char* inibdryfile,
			    const char* particlefile, 
			    const char* boundaryfile,
			    const char* contactfile,  
			    const char* progressfile,
			    const char* balancedfile, 
			    const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
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
    g_iteration=0;
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
	rigidBoundaryForce();
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=Volume/getParticleVolume()-1;

	if (sigma3_1<sigma_h1)
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma_h1)
	    minctl[1].tran=vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	if (sigma2_1<sigma_l1)
	    midctl[0].tran=vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma_l1)
	    midctl[1].tran=vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma_w1)
	    maxctl[0].tran=vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma_w1)
	    maxctl[1].tran=vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=vec(0,-TIMESTEP*RELEASE_RATE,0);
	
	updateRB(min,minctl,2);
	updateRB(mid,midctl,2);
	updateRB(max,maxctl,2);
	updateRB6();
	
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0 ){
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		        << endl;
	    g_debuginf << setw(OWID) << g_iteration
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
	    balancedinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()    
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		        << endl;
	    sigma_w1 += sigma_w_inc;
	    sigma_l1 += sigma_l_inc;
	    sigma_h1 += sigma_h_inc;
	}

	// 9. loop break condition
	if (   fabs(sigma1_1-sigma_w)/sigma_w < STRESS_ERROR && fabs(sigma1_2-sigma_w)/sigma_w < STRESS_ERROR
	    && fabs(sigma2_1-sigma_l)/sigma_l < STRESS_ERROR && fabs(sigma2_2-sigma_l)/sigma_l < STRESS_ERROR
	    && fabs(sigma3_1-sigma_h)/sigma_h < STRESS_ERROR && fabs(sigma3_2-sigma_h)/sigma_h < STRESS_ERROR ) {
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()    
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2
		        << setw(OWID) << sigma2_1 << setw(OWID) << sigma2_2
		        << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << getAverageRigidPressure()  // just the mean stress p
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		        << endl;
	    break;
	}
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}

//} // namespace dem ends

/* 
void assembly::dircShear(REAL rate, REAL roterate,REAL stress,const char* iniptclfile,
						 const char* boundaryfile, const char* responsefile, const char* resultfile,
						 const char* trackfile){
	readSample(iniptclfile);//create particles 
	readBoundary(boundaryfile);//create rigid boundaries

	FILE* fprslt=fopen(responsefile,"w");

	FILE* fp;
	fp=fopen(trackfile,"w");

	clearForce();

	int upanddown[2]={5,6};
	UPDATECTL updownctl[2];
	updownctl[0].expnd=1;
	updownctl[0].fixpt=0;
	updownctl[0].rote=0;
	updownctl[0].tran=vec(0,0,rate*TIMESTEP);
	updownctl[1]=updownctl[0];

	int load[2]={2,4};
	UPDATECTL loadctl[2];
	loadctl[0].expnd=1;
	loadctl[0].fixpt=getApt(2);
	loadctl[0].rote=vec(roterate*TIMESTEP,0,0);
	loadctl[0].tran=0;
	loadctl[1]=loadctl[0];
	loadctl[1].fixpt=getApt(4);

	std::vector<RGDBDRY*>::iterator rt;
	fprintf(fprslt,"bdry_1_norm_x  bdry_1_norm_y  bdry_1_norm_z  bdry_1_shar_x  bdry_1_shar_y  bdry_1_shar_z  \
bdry_2_norm_x  bdry_2_norm_y  bdry_2_norm_z  bdry_2_shar_x  bdry_2_shar_y  bdry_2_shar_z  \
bdry_3_norm_x  bdry_3_norm_y  bdry_3_norm_z  bdry_3_shar_x  bdry_3_shar_y  bdry_3_shar_z  \
bdry_4_norm_x  bdry_4_norm_y  bdry_4_norm_z  bdry_4_shar_x  bdry_4_shar_y  bdry_4_shar_z  \
bdry_5_norm_x  bdry_5_norm_y  bdry_5_norm_z  bdry_5_shar_x  bdry_5_shar_y  bdry_5_shar_z  \
bdry_6_norm_x  bdry_6_norm_y  bdry_6_norm_z  bdry_6_shar_x  bdry_6_shar_y  bdry_6_shar_z\n");

	REAL avgsigma;
	REAL l13, l24, l56, min_area, mid_area, max_area, lead;
	REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
	vec tmpnorm, tmpshar;
	REAL av=0;
	REAL ao=0;
	REAL af=0;
	REAL am=0;
	REAL avgNormal=0;
	REAL avgTangt=0;

	progressinf << "DircShearing..." << endl
	          << "iter_num   "
                  << "init_contact  "
                  << "contact  "
	          << "normal force      "
                  << "velocity        "
	          << "omga            "
	          << "force           "
	          << "moment" << endl;

	g_iteration=0;
	do{
                cout << "DircShearing..." << g_iteration << endl;
		progressinf << setw(OWID) << g_iteration;

		findParticleOnBoundary();
		findContact();

		internalForce(avgNormal, avgTangt);
		rigidBoundaryForce();
		//track(fp,5);
		progressinf << setw(OWID) << getAverageVelocity()
		          << setw(OWID) << getAverageOmga()
		          << setw(OWID) << getAverageForce()
			  << setw(OWID) << getAverageMoment() << endl;
		contactUpdate();
		updateParticle();

		l56=getApt(5).getz()-getApt(6).getz();
		l24=getApt(2).gety()-getApt(4).gety();
		l13=getApt(1).getx()-getApt(3).getx();
		min_area=l13*l24;
		mid_area=l56*l24;
		lead=fabs(normalize(getDirc(2))%vec(0,1,0));
		max_area=l56*l13;
		setArea(5,min_area);
		setArea(6,min_area);
		setArea(1,mid_area);
		setArea(3,mid_area);
		setArea(2,max_area);
		setArea(4,max_area);
		avgsigma=getAverageRigidPressure();
		printf("avgsigma=%15.3lf\n",avgsigma);
		sigma1_1=fabs(getNormalForce(2))/max_area;
		sigma1_2=fabs(getNormalForce(4))/max_area;
		sigma2_1=fabs(getNormalForce(1))/mid_area;
		sigma2_2=fabs(getNormalForce(3))/mid_area;
		sigma3_1=fabs(getNormalForce(5))/min_area;
		sigma3_2=fabs(getNormalForce(6))/min_area;
		if(sigma3_1<stress)
			updownctl[0].tran=vec(0,0,-rate*TIMESTEP);
		else
			updownctl[0].tran=vec(0,0,rate*TIMESTEP);
		if(sigma3_2<stress)
			updownctl[1].tran=vec(0,0,rate*TIMESTEP);
		else
			updownctl[1].tran=vec(0,0,-rate*TIMESTEP);
		updateRB(upanddown,updownctl,2);
		if (1){
			updateRB(load,loadctl,2);
			fprintf(fprslt,"%15.6lf",lead);
			for(rt=RBVec.begin();rt!=RBVec.end();++rt){
				tmpnorm=(*rt)->getNormalForce()/(*rt)->getArea()/1000;
				tmpshar=(*rt)->getShearForce()/(*rt)->getArea()/1000;
				fprintf(fprslt,"%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf",
					tmpnorm.x,tmpnorm.y,tmpnorm.z,
					tmpshar.x,tmpshar.y,tmpshar.z);
			}
			fprintf(fprslt,"\n");
		}
	}while(++g_iteration<10000);
	fclose(fp);
	fclose(fprslt);
	printParticle(resultfile);
}
*/

/* 
void assembly::soft_tric(REAL _sigma3,REAL _b,const char* iniptclfile,
						   const char* boundaryfile,const char* responsefile,
						   const char* resultfile,const char* trackfile){
	readSample(iniptclfile); //create particles 
	readBoundary(boundaryfile);

	FILE* fprslt=fopen(responsefile,"w");
	FILE* fp=fopen(trackfile,"w");
	
	clearForce();

	int pre_it=0;
	int pre_snap=0;
	int snapnum=0;
	char snapfile[80];

	std::vector<RGDBDRY*>::iterator rt;

	int max[2]={1,2};//maximum stress acting on boundary 5 and 6
	UPDATECTL maxctl[2];
	REAL loading_rate=0.01;

	REAL avgsigma;
	REAL af, av, am, ao, adr, pre_af;
	vec disp, tmp;
	av=ao=af=am=adr=pre_af=0;

	progressinf << "Soft_tric..." << endl
	          << "iter_num   "
                  << "init_contact  "
                  << "contact  "
	          << "normal force      "
                  << "velocity        "
	          << "omga            "
	          << "force           "
	          << "moment          "
		  << "friction" << endl;

	g_iteration=0;
	do{
                cout << "Soft_tric..." << g_iteration << endl;
		progressinf << setw(OWID) << g_iteration;

		findParticleOnBoundary();
		findParticleOnLine();
		createFlbNet();
		flexiBoundaryForceZero();
		flexiBoundaryForce();
		findContact();

		initFBForce();
		internalForce();
		rigidBoundaryForce();
		//track(fp,5);
		progressinf << setw(OWID) << (av=getAverageVelocity())
		          << setw(OWID) << (ao=getAverageOmga())
		          << setw(OWID) << (af=getAverageForce())
		          << setw(OWID) << (am=getAverageMoment())
			  << setw(OWID) << (adr=avgDgrFric());
		contactUpdate();
		
		avgsigma=getAverageRigidPressure();
		maxctl[0].tran=TIMESTEP*vec(0,0,-loading_rate);
		maxctl[1].tran=TIMESTEP*vec(0,0,loading_rate);
		if(af<0.03&&g_iteration-pre_it>=20||g_iteration-pre_it>=500){
		//if(g_iteration-pre_it>=50){
			pre_it=g_iteration;
		        updateRB(max,maxctl,2);
			if(g_iteration-pre_snap>=5000){
			    snapnum++;
			    sprintf(snapfile,"%s%d","snap",snapnum);
			    pre_snap=g_iteration;
			    snapshot(snapfile);
			}
			fprintf(fprslt,"%10d%6d%10.3lf%15.6lf%15.6lf%15.6lf%15.6lf",g_iteration,ActualCntctNum,adr,af,am,av,ao);
			for(rt=RBVec.begin();rt!=RBVec.end();++rt){
				disp=(*rt)->getApt();
				tmp=(*rt)->getNormalForce()/1000/(*rt)->getArea();
				fprintf(fprslt,"%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf",
					disp.getx(),disp.gety(),disp.getz(),tmp.getx(),tmp.gety(),tmp.getz());
			}
			fprintf(fprslt,"\n");
		}
	updateParticle();
	pre_af=af;
	}while(++g_iteration<1000000);
	fclose(fp);
	fclose(fprslt);
	printParticle(resultfile);
}//end of soft_tric
*/

/* 
void assembly::shallowFoundation(const char* iniptclfile, const char* boundaryfile,const char* responsefile, 
	const char* resultfile, const char* trackfile)
{
	readSample(iniptclfile);//create particles 
	readBoundary(boundaryfile);

	FILE* fprslt=fopen(responsefile,"w");

	FILE* fp;
	fp=fopen(trackfile,"w");

	int pre_it=0;
	int pre_snap=0;
	int snapnum=0;
	char snapfile[80];

	clearForce();

	std::vector<RGDBDRY*>::iterator rt;

//	int mid[2]={2,4};//intermediate stress acting on boundary 2 and 4
//	UPDATECTL midctl[2];
	int max[2]={5,6};//maximum stress acting on boundary 5 and 6
	UPDATECTL maxctl[2];
//	int min[2]={1,3};//minimum stress acting on boundary 1 and 3
//	UPDATECTL minctl[2];
	REAL loading_rate=0.01;

	REAL avgsigma;
	REAL af, av, am, ao, adr, pre_af;
	int nbdry;
	vec disp, tmp, zbdry_velocity_0;
	av=ao=af=am=adr=pre_af=0;

	progressinf << "Shallow Foundation..." << endl
	          << "iter_num   "
                  << "init_contact  "
                  << "contact  "
	          << "normal force      "
                  << "velocity        "
	          << "omga            "
	          << "force           "
	          << "moment" << endl;

	g_iteration=0;
	do{
                cout << "Shallow Foundation..." << g_iteration << endl;
		progressinf << setw(OWID) << g_iteration;

		findParticleOnBoundary();
		findParticleOnLine();
		createFlbNet();
		flexiBoundaryForceZero();
		flexiBoundaryForce();
		findContact();
		
		initFBForce();

		internalForce();
		rigidBoundaryForce();
		//track(fp,5);
		progressinf << setw(OWID) << (av=getAverageVelocity())
		          << setw(OWID) << (ao=getAverageOmga())
		          << setw(OWID) << (af=getAverageForce())
		          << setw(OWID) << (am=getAverageMoment())
			  << setw(OWID) << (adr=avgDgrFric());

		contactUpdate();
		
		avgsigma=getAverageRigidPressure();
		zbdry_velocity_0=vec(0,0,-loading_rate);
		maxctl[0].tran=TIMESTEP*zbdry_velocity_0;
		if(1){
			if(af<0.05&&g_iteration-pre_it>=20
				  ||g_iteration-pre_it>=500){
			pre_it=g_iteration;
		        updateRB(&max[0],&maxctl[0],1);
			nbdry=1;

			if(g_iteration-pre_snap>=20000){
			    snapnum++;
			    sprintf(snapfile,"%s%d","snap",snapnum);
			    pre_snap=g_iteration;
			    snapshot(snapfile);
			}
			fprintf(fprslt,"%10d%6d%10.3lf%15.6lf%15.6lf%15.6lf%15.6lf",g_iteration,ActualCntctNum,adr,af,am,av,ao);
			for(rt=RBVec.begin();rt!=RBVec.end();++rt,++nbdry){
				disp=(*rt)->getApt();
				tmp=(*rt)->getNormalForce()/1000;
				tmp/=getArea(nbdry);
				fprintf(fprslt,"%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf",
					disp.getx(),disp.gety(),disp.getz(),tmp.getx(),tmp.gety(),tmp.getz());
			}
			fprintf(fprslt,"\n");
			}
		}
		updateParticle();
		pre_af=af;
	}while(++g_iteration<1000000);
	fclose(fp);
	fclose(fprslt);
	printParticle(resultfile);
}
*/

/* 
void assembly::simpleShear(REAL _sigma3,REAL _b,
			const char* iniptclfile,const char* boundaryfile,
			const char* responsefile,const char* resultfile, const char* trackfile)
{
	readSample(iniptclfile);//create particles 
	readBoundary(boundaryfile);
	FILE* fprslt=fopen(responsefile,"w");

	FILE* fp;
	fp=fopen(trackfile,"w");

	clearForce();

	int pre_it=0;
	int pre_snap=0;
	int snapnum=0;
	char snapfile[80];

	std::vector<RGDBDRY*>::iterator rt;

	int mid[2]={2,4};//intermediate stress acting on boundary 2 and 4
	UPDATECTL midctl[2];
	int max[2]={5,6};//maximum stress acting on boundary 5 and 6
	UPDATECTL maxctl[2];
	int min[2]={1,3};//minimum stress acting on boundary 1 and 3
	UPDATECTL minctl[2];
//	REAL loading_rate=0.01;
	REAL increment=0.0001;
	REAL angular_velocity=0.1;
	vec increment_velocity_x(increment,0,0);
	vec increment_velocity_y(0,increment,0);
	vec increment_velocity_z(0,0,increment);
	vec xbdry_velocity_0,xbdry_velocity_1;
	vec ybdry_velocity_0,ybdry_velocity_1;
	vec zbdry_velocity_0,zbdry_velocity_1;

	REAL avgsigma;
	REAL af, av, am, ao, adr, pre_af;
	REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2, sigma2;
	REAL ita1_1, ita1_2, ita2_1, ita2_2, ita3_1, ita3_2;
	REAL ar;
	int nbdry;
	vec disp, angl, nm, sh;
	av=ao=af=am=adr=pre_af=0;

	progressinf << "SimpleShearing..." << endl
	          << "iter_num   "
                  << "init_contact  "
                  << "contact  "
	          << "normal force      "
                  << "velocity        "
	          << "omga            "
	          << "force           "
	          << "moment          "
		  << "friction" << endl;

	g_iteration=0;
	do{
                cout << "SimpleShearinging..." << g_iteration << endl;
		progressinf << setw(OWID) << g_iteration;

		findParticleOnBoundary();
		findContact();

		internalForce();
		rigidBoundaryForce();
		flexiBoundaryForce();
		//track(fp,5);
		progressinf << setw(OWID) << (av=getAverageVelocity())
		          << setw(OWID) << (ao=getAverageOmga())
		          << setw(OWID) << (af=getAverageForce())
		          << setw(OWID) << (am=getAverageMoment())
			  << setw(OWID) << (adr=avgDgrFric());
		contactUpdate();
		
		avgsigma=getAverageRigidPressure();
		sigma1_1=fabs(getNormalForce(5))/getArea(5);
		ita1_1=fabs(getShearForce(5))/getArea(5);
		sigma1_2=fabs(getNormalForce(6))/getArea(6);
		ita1_2=fabs(getShearForce(6))/getArea(6);
		sigma2_1=fabs(getNormalForce(2))/getArea(2);
		ita2_1=fabs(getShearForce(2))/getArea(2);
		sigma2_2=fabs(getNormalForce(4))/getArea(4);
		ita2_2=fabs(getShearForce(4))/getArea(4);
		sigma3_1=fabs(getNormalForce(1))/getArea(1);
		ita3_1=fabs(getShearForce(1))/getArea(1);
		sigma3_2=fabs(getNormalForce(3))/getArea(3);
		ita3_2=fabs(getShearForce(3))/getArea(1);
		if(sigma3_1<_sigma3)
			xbdry_velocity_0=-increment_velocity_x;
		else
			xbdry_velocity_0=increment_velocity_x;
		if(sigma3_2<_sigma3)
			xbdry_velocity_1=increment_velocity_x;
		else
			xbdry_velocity_1=-increment_velocity_x;
		sigma2=_sigma3;//_b*sigma1+(1-_b)*_sigma3;
		if(sigma2_1<sigma2)
			ybdry_velocity_0=-increment_velocity_y;
		else
			ybdry_velocity_0=increment_velocity_y;
		if(sigma2_2<sigma2)
			ybdry_velocity_1=increment_velocity_y;
	 	else
			ybdry_velocity_1=-increment_velocity_y;
		if(sigma1_1<_sigma3)
			zbdry_velocity_0=-increment_velocity_z;
		else
			zbdry_velocity_0=increment_velocity_z;
		if(sigma1_2<_sigma3)
			zbdry_velocity_1=increment_velocity_z;
		else
			zbdry_velocity_1=-increment_velocity_z;
		minctl[0].tran=TIMESTEP*xbdry_velocity_0;
		minctl[0].fixpt=getApt(1);
		minctl[0].rote=TIMESTEP*vec(0,0,angular_velocity);
		minctl[1].fixpt=getApt(3);
		minctl[1].tran=TIMESTEP*xbdry_velocity_1;
		minctl[1].rote=TIMESTEP*vec(0,0,angular_velocity);
		midctl[0].fixpt=getApt(2);
		midctl[0].tran=TIMESTEP*ybdry_velocity_0;
		midctl[0].rote=TIMESTEP*vec(0,0,-angular_velocity);
		midctl[1].fixpt=getApt(4);
		midctl[1].tran=TIMESTEP*ybdry_velocity_1;
		midctl[1].rote=TIMESTEP*vec(0,0,-angular_velocity);
		maxctl[0].tran=TIMESTEP*zbdry_velocity_0;
		maxctl[1].tran=TIMESTEP*zbdry_velocity_1;
		//if(af<0.01){
		if(1){
			UPDATECTL tmpctl;

			tmpctl.tran=minctl[0].tran;
			updateRB(&min[0],&tmpctl,1);

			tmpctl.tran=minctl[1].tran;
			updateRB(&min[1],&tmpctl,1);

			tmpctl.tran=midctl[0].tran;
			updateRB(&mid[0],&tmpctl,1);

			tmpctl.tran=midctl[1].tran;
			updateRB(&mid[1],&tmpctl,1);

			tmpctl.tran=maxctl[0].tran;
			updateRB(&max[0],&tmpctl,1);

			tmpctl.tran=maxctl[1].tran;
			updateRB(&max[1],&tmpctl,1);
			
			if(af<0.02&&fabs(sigma3_1-_sigma3)<0.02*_sigma3
				  &&fabs(sigma3_2-_sigma3)<0.02*_sigma3
				  &&fabs(sigma2_1-sigma2)<0.02*_sigma3
				  &&fabs(sigma2_2-sigma2)<0.02*_sigma3
				  &&fabs(sigma1_1-_sigma3)<0.02*_sigma3
				  &&fabs(sigma1_2-_sigma3)<0.02*_sigma3
				  &&g_iteration-pre_it>=20
				  ||g_iteration-pre_it>=500){
			pre_it=g_iteration;
		        updateRB(min,minctl,2);
			updateRB(mid,midctl,2);
			nbdry=1;

			if(g_iteration-pre_snap>=20000){
			    snapnum++;
			    sprintf(snapfile,"%s%d","snap",snapnum);
			    pre_snap=g_iteration;
			    snapshot(snapfile);
			}
			fprintf(fprslt,"%10d%6d%10.3lf%15.6lf%15.6lf%15.6lf%15.6lf",g_iteration,ActualCntctNum,adr,af,am,av,ao);
			for(rt=RBVec.begin();rt!=RBVec.end();++rt,++nbdry){
				ar=(*rt)->getArea();
				disp=(*rt)->getApt();
				angl=(*rt)->getDirc();
				nm=(*rt)->getNormalForce()/1000/ar;
				sh=(*rt)->getShearForce()/1000/ar;
				fprintf(fprslt,"%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf", 
disp.x,disp.y,disp.z,angl.x,angl.y,angl.z,nm.x,nm.y,nm.z,sh.x,sh.y,sh.z);
			}
			fprintf(fprslt,"\n");
			}
		}
	updateParticle();
	pre_af=af;
	}while(++g_iteration<300000);
	fclose(fp);
	fclose(fprslt);
	printParticle(resultfile);
}
*/

/* 
void assembly::earthPressure(REAL pressure,bool IsPassive, 
				const char* iniptclfile, const char* boundaryfile,
				const char* responsefile, const char* resultfile,
				const char* trackfile)
{
	readSample(iniptclfile);//create particles 
	readBoundary(boundaryfile);

	FILE* fprslt=fopen(responsefile,"w");

	FILE* fp;
	fp=fopen(trackfile,"w");

	int pre_it=0;
	int pre_snap=0;
	int snapnum=0;
	char snapfile[80];
	clearForce();

	std::vector<RGDBDRY*>::iterator rt;

	int wall[1]={1};
	UPDATECTL wallctl[1];
//	REAL loading_rate=0.001;

	REAL avgsigma;
	REAL af, av, am, ao, adr, pre_af;
	int nbdry;
	vec disp, tmp;
	av=ao=af=am=adr=pre_af=0;

	progressinf << "EarthPressure..." << endl
	          << "iter_num   "
                  << "init_contact  "
                  << "contact  "
	          << "normal force      "
                  << "velocity        "
	          << "omga            "
	          << "force           "
	          << "moment          "
		  << "friction" << endl;

	g_iteration=0;
	do{
	        cout << "EarthPressure..." << g_iteration << endl;
		progressinf << setw(OWID) << g_iteration;

		findParticleOnBoundary();
		findParticleOnLine();
		createFlbNet();
		flexiBoundaryForceZero();
		flexiBoundaryForce();
		findContact();
		
                initFBForce();
		internalForce();
		rigidBoundaryForce();
		//track(fp,5);

		progressinf << setw(OWID) << (av=getAverageVelocity())
		          << setw(OWID) << (ao=getAverageOmga())
		          << setw(OWID) << (af=getAverageForce())
		          << setw(OWID) << (am=getAverageMoment())
			  << setw(OWID) << (adr=avgDgrFric());
		contactUpdate();
		
		avgsigma=getAverageRigidPressure();
		wallctl[0].tran=0;
		wallctl[0].fixpt=getApt(1);
		if(IsPassive)
			wallctl[0].rote=TIMESTEP*vec(0,-1.0,0);
		else
			wallctl[0].rote=TIMESTEP*vec(0,1.0,0);
		if(1){
			if(af<0.05 &&g_iteration-pre_it>=20
				  ||g_iteration-pre_it>=500){
			pre_it=g_iteration;
		        updateRB(wall,wallctl,1);
			nbdry=1;

			if(g_iteration-pre_snap>=20000){
			    snapnum++;
			    sprintf(snapfile,"%s%d","snap",snapnum);
			    pre_snap=g_iteration;
			    snapshot(snapfile);
			}
			fprintf(fprslt,"%10d%6d%10.3lf%15.6lf%15.6lf%15.6lf%15.6lf",g_iteration,ActualCntctNum,adr,af,am,av,ao);
			for(rt=RBVec.begin();rt!=RBVec.end();++rt,++nbdry){
				if(nbdry==1)
					disp=(*rt)->getDirc();
				else
					disp=(*rt)->getApt();
				tmp=(*rt)->getNormalForce()/1000/getArea(nbdry);
				fprintf(fprslt,"%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf",
					disp.getx(),disp.gety(),disp.getz(),tmp.getx(),tmp.gety(),tmp.getz());
			}
			fprintf(fprslt,"\n");
			}
		}
	updateParticle();
	pre_af=af;
	}while(++g_iteration<1000000);
	fclose(fp);
	fclose(fprslt);
	printParticle(resultfile);
}
*/

// calculate granular stress for assembly from input files, such as the assembly after deposition.
// writtend on Feb 16, 2013
void assembly::calculateStress(REAL   relativeHeight,	// ratio of height of each subdomain to boundary height (assume they are equal)
			       REAL   relativeWidth,	// ratio of width of each subdomain to boundary width
			       REAL   relativeLength,	// ratio of length of each subdomain to boundary length
			       std::vector<REAL> &relBaseH,	// ratio of the base height of each subdomain above the assembly boundary base
			       const char* inptclfile,	// particle file, can use trm_particle_end for deposited assembly
			       const char* inbdryfile,
			       const char* resultfile)
{
	// pre_1: open streams for output
	std::ofstream resultstress;
	resultstress.open(resultfile);
	if(!resultstress) { cout << "stream error!" << endl; exit(-1);}
	resultstress.setf(std::ios::scientific, std::ios::floatfield);
	resultstress << setw(OWID) << "total"
		     << setw(OWID) << ""
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"	// volume
		     << setw(OWID) << "contact point"
		     << setw(OWID) << "formula"
		     << endl
		     << setw(OWID) << "Height"
		     << setw(OWID) << "Depth"
		     << setw(OWID) << "sigma_11"
		     << setw(OWID) << "sigma_12"
		     << setw(OWID) << "sigma_13"
		     << setw(OWID) << "sigma_21"
		     << setw(OWID) << "sigma_22"
		     << setw(OWID) << "sigma_23"
		     << setw(OWID) << "sigma_31"
		     << setw(OWID) << "sigma_32"
		     << setw(OWID) << "sigma_33"
		     << setw(OWID) << "volume"
		     << setw(OWID) << "number"
		     << setw(OWID) << "stress"
		     << endl;
	
	// pre_2. create particle and boundaries from files
	readSample(inptclfile);	// create container and particles, velocity and omga are set zero. 
	readBoundary(inbdryfile);	// create boundaries

	// pre_3. define variables used in iterations
	REAL bdryW = getApt(2).gety()-getApt(4).gety();
    	REAL bdryL = getApt(1).getx()-getApt(3).getx();
	REAL bdryH = getApt(5).getz()-getApt(6).getz();	// height of boundary box
	bdryH = 0.13;	// after trmByHeight();
	REAL subHeight = bdryH*relativeHeight;	// the actual height of subdomain
	REAL subWidth = bdryW*relativeWidth;
	REAL subLength = bdryL*relativeLength;
	REAL subArea = subWidth*subLength;

	REAL baseCoordH = getApt(6).getz();	// the coordinate of base of boundary box
	REAL baseCoordL = getApt(3).getx()+(bdryL-subLength)*0.5;
	REAL baseCoordW = getApt(4).gety()+(bdryW-subWidth)*0.5;

	REAL subVolume = subHeight*subWidth*subLength;
	vector<REAL> subBaseCoord;	// actual height above the boundary base of base of each subdomain
	vector<REAL> baseH;
	for(vector<REAL>::const_iterator it=relBaseH.begin(); it!=relBaseH.end(); it++){
		subBaseCoord.push_back(bdryH*(*it)+baseCoordH);
		baseH.push_back(bdryH*(*it));
	}

	// not used
	REAL avgNormal=0;
    	REAL avgTangt=0;
	REAL bdry_penetr[7];
        int         bdry_cntnum[7];
        for (int i=0;i<7;++i){
		bdry_penetr[i]=0; bdry_cntnum[i]=0;
        }

	// calculate granular stress
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces/moments as zero 
	clearForce();
	
	// 3. calcualte contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
//	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	// calculate granular stress
	vector<matrix> granularStress;	// granular stress for each subdomain
	vector<matrix> totalForce;
	vector<int> subPctlNum;	// number of particle of each subdomain
	vector<REAL> subMass;	// mass of particles above each subdomain bottom
	matrix lc(3,1), Fc(1,3);	// sigma(i,j) = lixFj
	matrix P1center(3,1), P2center(3,1);
	
	// generate the same number of elements in granularStress and totalForce as relBaseH
	matrix temp_mat(3,3);
	// initialize totalForce
	for(int ir=0; ir!=3; ir++)
		for(int ic=0; ic!=3; ic++)
			temp_mat(ir+1,ic+1) = 0;
	for(vector<REAL>::const_iterator it=relBaseH.begin(); it!=relBaseH.end(); it++){
		granularStress.push_back(temp_mat);
		totalForce.push_back(temp_mat);	// every total Force is initialized to be zero
		subPctlNum.push_back(0);
		subMass.push_back(0);
	}

	// get density of each subdomain
	std::vector<particle*>::const_iterator  it_p;
	REAL x_p, y_p, z_p;
	vector<REAL>::size_type it_mass = 0;
	bool is_pctl_in;
  	for (it_p=ParticleVec.begin();it_p!=ParticleVec.end();++it_p)  {	
		x_p = (*it_p)->getCurrPosition().getx();
		y_p = (*it_p)->getCurrPosition().gety();
		z_p = (*it_p)->getCurrPosition().getz();
		it_mass = 0;
		for(std::vector<REAL>::const_iterator it=subBaseCoord.begin(); it!=subBaseCoord.end(); it++){
			is_pctl_in = z_p>(*it)+subHeight/2 && z_p<(baseCoordH+bdryH)
				   && y_p>baseCoordW && y_p<(baseCoordW+subWidth)
				   && x_p>baseCoordL && x_p<(baseCoordL+subLength);	// l->x, w->y, h->z
			if(is_pctl_in)
				subMass[it_mass] += (*it_p)->getMass();
			it_mass++;
		}
	}
	
	// below is to calculate the granular stress

	// sum the total particle contact force
	vector<matrix>::size_type it_mat = 0;
	vector<int>::size_type it_subnum = 0;

 	std::vector<CONTACT>::const_iterator it;
	bool is_con1_in, is_con2_in;
	REAL x_con1, y_con1, z_con1, x_con2, y_con2, z_con2;	// contact coordinates
	for (it=ContactVec.begin();it!=ContactVec.end();++it){
		P1center(1,1) = it->getP1()->getCurrPosition().getx();
		P1center(2,1) = it->getP1()->getCurrPosition().gety();
		P1center(3,1) = it->getP1()->getCurrPosition().getz();

		P2center(1,1) = it->getP2()->getCurrPosition().getx();
		P2center(2,1) = it->getP2()->getCurrPosition().gety();
		P2center(3,1) = it->getP2()->getCurrPosition().getz();

		lc = P2center-P1center;
		Fc(1,1) = (it->NormalForceVec().getx())+(it->TgtForceVec().getx());
		Fc(1,2) = (it->NormalForceVec().gety())+(it->TgtForceVec().gety());
		Fc(1,3) = (it->NormalForceVec().getz())+(it->TgtForceVec().getz());
		// get contact points coordinates
		x_con1 = it->getPoint1().getx();		
		y_con1 = it->getPoint1().gety();	
		z_con1 = it->getPoint1().getz();	
		x_con2 = it->getPoint2().getx();		
		y_con2 = it->getPoint2().gety();	
		z_con2 = it->getPoint2().getz();	
		// judge that which subdomain this pair of particles belong to
		it_subnum = 0;
		it_mat = 0;	// start from the first subdomain
		for(std::vector<REAL>::const_iterator it_r=subBaseCoord.begin(); it_r!=subBaseCoord.end(); it_r++){
			is_con1_in = z_con1>(*it_r) && z_con1<((*it_r)+subHeight)
				   && y_con1>baseCoordW && y_con1<(baseCoordW+subWidth)
				   && x_con1>baseCoordL && x_con1<(baseCoordL+subLength);	// l->x, w->y, h->z
			is_con2_in = z_con2>(*it_r) && z_con2<((*it_r)+subHeight)
				   && y_con2>baseCoordW && y_con2<(baseCoordW+subWidth)
				   && x_con2>baseCoordL && x_con2<(baseCoordL+subLength);	// l->x, w->y, h->z;
			if( is_con1_in && is_con2_in ){	// tow contact points are both in the subdomain
				totalForce[it_mat] += lc*Fc;
				subPctlNum[it_subnum] += 2;
			}
			else if( is_con1_in || is_con2_in ){// have and only have one point in the subdomain
				totalForce[it_mat] += lc*Fc/2;
				subPctlNum[it_subnum] += 1;
			}
			it_mat++;	// go to the next subdomain to make the same judgement
			it_subnum++;
		}// notation: our subdomains can be overlaped with each other
	}
	it_mat = 0;
	for(it_mat=0; it_mat!=totalForce.size(); it_mat++){
		granularStress[it_mat] = totalForce[it_mat]/subVolume;
	}
	// above is to calculate granular stress

	// output stress
	it_subnum = 0;
	it_mass = 0;
	for(it_mat=0; it_mat!=granularStress.size(); it_mat++, it_subnum++, it_mass++){
		resultstress << setw(OWID) << bdryH
			     << setw(OWID) << bdryH-baseH[it_mat]-subHeight/2
			     << setw(OWID) << granularStress[it_mat](1,1)
			     << setw(OWID) << granularStress[it_mat](1,2)
			     << setw(OWID) << granularStress[it_mat](1,3)
			     << setw(OWID) << granularStress[it_mat](2,1)
			     << setw(OWID) << granularStress[it_mat](2,2)
			     << setw(OWID) << granularStress[it_mat](2,3)
			     << setw(OWID) << granularStress[it_mat](3,1)
			     << setw(OWID) << granularStress[it_mat](3,2)
			     << setw(OWID) << granularStress[it_mat](3,3)
			     << setw(OWID) << subVolume
			     << setw(OWID) << subPctlNum[it_subnum]
			     << setw(OWID) << subMass[it_mass]/subArea
			     << endl;
	}	
}

void assembly::trimByHeight(REAL Height,
			const char* inptclefile,	// initial particle file
			const char* inbdryfile,		// initial boundary file
			const char* resultptcle,
			const char* resultbdry){	// result particle file

	readSample(inptclefile);	// create container and particles, velocity and omga are set zero. 
	readBoundary(inbdryfile);	// create boundaries

	// pre_3. define variables used in iterations
	REAL bdryW = getApt(2).gety()-getApt(4).gety();
    	REAL bdryL = getApt(1).getx()-getApt(3).getx();
	REAL bdryH = getApt(5).getz()-getApt(6).getz();	// height of boundary box

	trimHistoryNum = TotalNum;

	REAL mass_removed = 0;	// the mass of removed particles, July 2, 2013. test the residual stress in stress distribution
 	REAL z1, z2;
 	z1 = getApt(6).getz();
	z2 = z1+Height;
 
  	std::vector<particle*>::iterator itr;
  	vec center;
  	REAL mass = 0;
	int num_iter = 0;
  	for (itr = ParticleVec.begin(); itr != ParticleVec.end(); ){
    		center=(*itr)->getCurrPosition();
    		if(center.getz()+(*itr)->getA() > z2)
      		{
			mass_removed += (*itr)->getMass();
			delete (*itr); // release memory
			itr = ParticleVec.erase(itr); 
	      	}
    		else
      			++itr;
		num_iter++;
//		std::cout << num_iter << std::endl;
 	 }
  
std::cout << "mass of removed particles: " << endl << mass_removed << endl;

  	for(itr=ParticleVec.begin();itr!=ParticleVec.end();++itr)
    		mass += (*itr)->getMass();
  
  	Volume = bdryW*bdryL*Height;
  	BulkDensity = mass/Volume;
  
  	TotalNum = ParticleVec.size();
  	printParticle(resultptcle);

	// write new boundary
  std::ofstream ofs(resultbdry);
  if(!ofs) { cout << "stream error!" << endl; exit(-1);}

  REAL x1,x2,y1,y2,x0,y0,z0;
  x1 = getApt(3).getx();
  y1 = getApt(4).gety();
  x2 = getApt(1).getx();
  y2 = getApt(2).gety();
  x0 = (x1+x2)/2;
  y0 = (y1+y2)/2;
  z0 = (z1+z2)/2;

  int bdrynum = 6;

  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs << setw(OWID) << 0
      << setw(OWID) << bdrynum << endl << endl;

  // boundary 1
  ofs << setw(OWID) << 1 << endl
      << setw(OWID) << 1
      << setw(OWID) << 5
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
      << setw(OWID) << 5
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
      << setw(OWID) << 5
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
      << setw(OWID) << 5
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
      << setw(OWID) << 5
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
      << setw(OWID) << 5
      << setw(OWID) << 0 << endl
    
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


void assembly::calculateMiddleStress(int subNum,	// the number of subdomains that are used to calculate formula stress
			       const char* inptclfile,	// particle file, can use trm_particle_end for deposited assembly
			       const char* inbdryfile,
			       const char* resultfile)
{
	// pre_1: open streams for output
	std::ofstream resultstress;
	resultstress.open(resultfile);
	if(!resultstress) { cout << "stream error!" << endl; exit(-1);}
	resultstress.setf(std::ios::scientific, std::ios::floatfield);
	resultstress << setw(OWID) << "total"
		     << setw(OWID) << ""
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "sample"	// volume
		     << setw(OWID) << "subdomain"
		     << setw(OWID) << "formula"
		     << setw(OWID) << "sample"
		     << endl
		     << setw(OWID) << "Height"
		     << setw(OWID) << "Depth"
		     << setw(OWID) << "sigma_11"
		     << setw(OWID) << "sigma_12"
		     << setw(OWID) << "sigma_13"
		     << setw(OWID) << "sigma_21"
		     << setw(OWID) << "sigma_22"
		     << setw(OWID) << "sigma_23"
		     << setw(OWID) << "sigma_31"
		     << setw(OWID) << "sigma_32"
		     << setw(OWID) << "sigma_33"
		     << setw(OWID) << "volume"
		     << setw(OWID) << "volume"
		     << setw(OWID) << "stress"
		     << setw(OWID) << "density"
		     << endl;
	
	// pre_2. create particle and boundaries from files
	readSample(inptclfile);	// create container and particles, velocity and omga are set zero. 
	readBoundary(inbdryfile);	// create boundaries

	// pre_3. define variables used in iterations
	REAL bdryW = getApt(2).gety()-getApt(4).gety();
    	REAL bdryL = getApt(1).getx()-getApt(3).getx();
	REAL bdryH = getApt(5).getz()-getApt(6).getz();	// height of boundary box
	bdryH = 0.003;	// after trmByHeight();

	REAL subArea = bdryL*bdryW;
	REAL subHeight = bdryH/subNum;

	REAL baseCoordH = getApt(6).getz();	// the coordinate of base of boundary box
	REAL baseCoordL = getApt(3).getx();
	REAL baseCoordW = getApt(4).gety();

	// not used
	REAL avgNormal=0;
    	REAL avgTangt=0;
	REAL bdry_penetr[7];
        int         bdry_cntnum[7];
        for (int i=0;i<7;++i){
		bdry_penetr[i]=0; bdry_cntnum[i]=0;
        }

	// calculate granular stress
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces/moments as zero 
	clearForce();
	
	// 3. calcualte contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
//	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	// calculate granular stress
	matrix granularStress;	// granular stress for each subdomain
	
	vector<REAL> subMass;	// mass of particles above each subdomain bottom
	REAL totalMass = 0;
	vector<REAL> subVolume;
	std::vector<REAL>::size_type it_mass = 0;
	for(int i=0; i!=subNum; i++){
		subMass.push_back(0);
		subVolume.push_back((bdryH-i*subHeight)*subArea);
	}
	// get density of each subdomain
	std::vector<particle*>::const_iterator  it_p;
	REAL x_p, y_p, z_p;

	bool is_pctl_in;
  	for (it_p=ParticleVec.begin();it_p!=ParticleVec.end();++it_p)  {	
		x_p = (*it_p)->getCurrPosition().getx();
		y_p = (*it_p)->getCurrPosition().gety();
		z_p = (*it_p)->getCurrPosition().getz();
		it_mass = 0;
		for(int i=0; i!=subNum; i++){
			is_pctl_in = z_p>baseCoordH+i*subHeight;	// l->x, w->y, h->z
			if(is_pctl_in)
				subMass[it_mass] += (*it_p)->getMass();
			it_mass++;
		}
		totalMass += (*it_p)->getMass();
	}
	
	// below is to calculate the granular stress

	// sum the total particle contact force
	granularStress = getGranularStress();

	// output stress
	it_mass = 0;
	for(int i=0; i!=subNum; i++, it_mass++){
		resultstress << setw(OWID) << bdryH
			     << setw(OWID) << bdryH-i*subHeight
			     << setw(OWID) << granularStress(1,1)
			     << setw(OWID) << granularStress(1,2)
			     << setw(OWID) << granularStress(1,3)
			     << setw(OWID) << granularStress(2,1)
			     << setw(OWID) << granularStress(2,2)
			     << setw(OWID) << granularStress(2,3)
			     << setw(OWID) << granularStress(3,1)
			     << setw(OWID) << granularStress(3,2)
			     << setw(OWID) << granularStress(3,3)
			     << setw(OWID) << bdryH*bdryL*bdryW
			     << setw(OWID) << subVolume[it_mass]
			     << setw(OWID) << subMass[it_mass]/subArea
			     << setw(OWID) << totalMass/(bdryH*bdryL*bdryW)
			     << endl;
	}	
}

void assembly::calculateFabric(const char* inptclfile,
			       const REAL boxLength,	// the height, width and length of box. used for porosity
			       const REAL boxWidth,	// the unit is mm
			       const REAL boxHeight,
			       const char* resultfile){
	// open streams for output
	std::ofstream resultfabric;
	resultfabric.open(resultfile);
	if(!resultfabric){cout << "Stream error when open fabric output file!" << endl; exit(-1);}
	resultfabric.setf(std::ios::scientific, std::ios::floatfield);
	resultfabric << setw(OWID) << "fabric"
		     << setw(OWID) << "fabric"
		     << setw(OWID) << "fabric"
		     << setw(OWID) << "fabric"
		     << setw(OWID) << "fabric"
		     << setw(OWID) << "fabric"
		     << setw(OWID) << "fabric"
		     << setw(OWID) << "fabric"
		     << setw(OWID) << "fabric"	// fabric tensor
		     << setw(OWID) << "first invariant"
		     << setw(OWID) << "second invariant"
		     << setw(OWID) << "third invariant"
		     << setw(OWID) << "sample"
		     << setw(OWID) << "anisotropy"
		     << setw(OWID) << "assembly"
		     << endl
		     << setw(OWID) << "F_11"
		     << setw(OWID) << "F_12"
		     << setw(OWID) << "F_13"
		     << setw(OWID) << "F_21"
		     << setw(OWID) << "F_22"
		     << setw(OWID) << "F_23"
		     << setw(OWID) << "F_31"
		     << setw(OWID) << "F_32"
		     << setw(OWID) << "F_33"
		     << setw(OWID) << "I1(F)"
		     << setw(OWID) << "sqrt[I2(d)]"
		     << setw(OWID) << "I3(F)" 
		     << setw(OWID) << "porosity" 
		     << setw(OWID) << "a"
		     << setw(OWID) << "density" << endl;

	// create particles
	readSample(inptclfile); // create container and particles, velocity and omga are set zero. 

	// define variables
	matrix fabricTensor(3,3);
	matrix Dtensor(3,3);	// D=F-I1(F)/3
	REAL F_first;	// the first invariant of fabric tensor
	REAL D_second;	// the secodn invariant of D tensor
	REAL F_third;	// the third invariant of fabric tensor
	REAL tao;	// tao = sqrt(I2(D))
	REAL boxVolume = boxHeight*boxWidth*boxLength*pow(10,-9);
std::cout << "boxVolume: " << boxVolume << std::endl;
	REAL sample_voidratio = boxVolume/getParticleVolume()-1;
std::cout << "void ratio: " << sample_voidratio << std::endl;
	REAL sample_porosity = sample_voidratio/(1+sample_voidratio);
std::cout << "porosity: " << sample_porosity << std::endl;

	// find contacts between particles
	findContact();
	
	fabricTensor = getFabric();	// get fabric tensor

	// calculate the first invariant of fabric tensor
	F_first = trace(fabricTensor);
	// calculate the second invariant of D and tao
	Dtensor = fabricTensor-F_first/3.0;
	D_second = fabs(trace(Dtensor)*trace(Dtensor)-trace(Dtensor*Dtensor))*0.5;	
	tao = sqrt(2*D_second);
	// calculate the third invariant of fabric tensor
	matrix fabric_3 = fabricTensor*fabricTensor*fabricTensor;
	matrix fabric_2 = fabricTensor*fabricTensor;
	F_third = (pow(trace(fabricTensor),3)+2*trace(fabric_3)-3*trace(fabricTensor)*trace(fabric_2))/6.0;
	
	// calculate a
	REAL a1, a2, a3, a;
	a1 = (5-15*fabricTensor(1,1))/(3-5*fabricTensor(1,1));
	a2 = (5-15*fabricTensor(2,2))/(3-5*fabricTensor(2,2));
	a3 = (15*fabricTensor(3,3)-5)/(5*fabricTensor(3,3)+1);
	a = (a1+a2+a3)/3.0;
std::cout << "a: " << a << std::endl;

    	REAL dens=0;
    	std::vector<particle*>::const_iterator it;
    	for(it=ParticleVec.begin();it!=ParticleVec.end();++it)
	    dens+=(*it)->getMass();

	dens = dens/boxVolume;

	// output information about fabric
	resultfabric << setw(OWID) << fabricTensor(1,1)
		     << setw(OWID) << fabricTensor(1,2)
		     << setw(OWID) << fabricTensor(1,3)
		     << setw(OWID) << fabricTensor(2,1)
		     << setw(OWID) << fabricTensor(2,2)
		     << setw(OWID) << fabricTensor(2,3)
		     << setw(OWID) << fabricTensor(3,1)
		     << setw(OWID) << fabricTensor(3,2)
		     << setw(OWID) << fabricTensor(3,3)
		     << setw(OWID) << F_first
		     << setw(OWID) << tao
		     << setw(OWID) << F_third
		     << setw(OWID) << sample_porosity
		     << setw(OWID) << a
		     << setw(OWID) << dens
		     << endl;
}


// reset coordinates of particles which are contacting with boundaries. And then call Qhull to tessellate. This function is written for Nonlinear FEM term project -- a 3D tet poromechanics Matlab code.
  void assembly::resetCoord_FEM(const char* inptclfile,
		      const char* inbdryfile){

	// pre_2. create particles and boundaries from files
	readSample(inptclfile); // create container and particles, velocity and omga are set zero. 
        readBoundary(inbdryfile);   // create boundaries
	
/*
	// find particle on boundary and then reset their coordinates
	findParticleOnBoundary();

	std::vector<particle*> temp_particle;
	std::vector<particle*>::const_iterator it_all;	// this used to loop over ParticleVec
	std::vector<particle*>::const_iterator  it_p;
	int temp_ID;	// particle ID
	int bdry_ID;	// boundary ID
	REAL x1,y1,z1;
	vec new_coord;	// new coordinate for boundary particle

	std::vector<RGDBDRY*>::iterator rt;    
	plnrgd_bdry<particle>* rbptr;
    	for(rt=RBVec.begin();rt!=RBVec.end();++rt){
		rbptr = static_cast<plnrgd_bdry<particle>* >(rt);	// convert from base class to derived class
		temp_particle = (*rt)->PBVec;
		bdry_ID = (*rt)->getBdryID();	// get boundary ID
		for (it_p=temp_particle.begin();it_p!=temp_particle.end();++it_p)  {	
			temp_ID = (*it_p)->getID();
			it_all = ParticleVec.begin();
			for(int i=0; i!=temp_ID-1; i++)
				it_all++;	// go to the present boundary particle
			x1 = ((*it_all)->getCurrPosition()).getx();
			y1 = ((*it_all)->getCurrPosition()).gety();
			z1 = ((*it_all)->getCurrPosition()).getz();
			// judge with which boundary surface this particle is contacting
			if(bdry_ID == 1)	// x: 1 3
				x1 = getApt(1).getx();	// set x coordinate of this particle as boundary coordinate

			if(bdry_ID == 2)	// y: 2 4
				y1 = getApt(2).gety();	// set y coordinate of this particle as boundary coordinate

			if(bdry_ID == 3)	// x: 1 3
				x1 = getApt(3).getx();	// set x coordinate of this particle as boundary coordinate

			if(bdry_ID == 4)	// y: 2 4
				y1 = getApt(4).gety();	// set y coordinate of this particle as boundary coordinate

			if(bdry_ID == 5)	// z: 5 6
				z1 = getApt(5).getz();	// set z coordinate of this particle as boundary coordinate

			if(bdry_ID == 6)	// z:5 6
				z1 = getApt(6).getz();	// set z coordinate of this particle as boundary coordinate

			// reset coordinates
			new_coord.setx(x1);
			new_coord.sety(y1);
			new_coord.setz(z1);
			(*it_all)->setCurrPosition(new_coord);
		}
	}
	
*/    

	// another way to find particles on boundary and reset their coordinates
	REAL x_min, x_max, y_min, y_max, z_min, z_max;	// the coordinates of boundary surfaces
	REAL x_p, y_p, z_p;	// coordinates of particle
	
	x_max = getApt(1).getx();
	x_min = getApt(3).getx();

	y_max = getApt(2).gety();
	y_min = getApt(4).gety();

	z_max = getApt(5).getz();
	z_min = getApt(6).getz();
	REAL bdry_H = 0.003;	// after trmByHeight
	z_max = z_min+bdry_H;	

	// print boundary information in the screen
	std::cout << "x_max: " << x_max << std::endl;	
	std::cout << "x_min: " << x_min << std::endl;	
	std::cout << "y_max: " << y_max << std::endl;	
	std::cout << "y_min: " << y_min << std::endl;	
	std::cout << "z_max: " << z_max << std::endl;	
	std::cout << "z_min: " << z_min << std::endl;	

	vec new_coord;	// new coordinate for boundary particle
	REAL max_r;	// max radius of particle
	std::vector<particle*>::const_iterator  it_p;
  	for (it_p=ParticleVec.begin();it_p!=ParticleVec.end();++it_p)  {	
		x_p = (*it_p)->getCurrPosition().getx();
		y_p = (*it_p)->getCurrPosition().gety();
		z_p = (*it_p)->getCurrPosition().getz();
		max_r = (*it_p)->getA();
	
		if(x_max-x_p < max_r )	// means particle is near to max x surface
			x_p = x_max;
		if(x_p-x_min < max_r)	// means particle is near to min x surface
			x_p = x_min;
	
		if(y_max-y_p < max_r )	// means particle is near to max y surface
			y_p = y_max;
		if(y_p-y_min < max_r)	// means particle is near to min y surface
			y_p = y_min;
	
		if(z_max-z_p < max_r )	// means particle is near to max z surface
			z_p = z_max;
		if(z_p-z_min < max_r)	// means particle is near to min z surface
			z_p = z_min;

		// reset coordinates
		new_coord.setx(x_p);
		new_coord.sety(y_p);
		new_coord.setz(z_p);
		(*it_p)->setCurrPosition(new_coord);
	}
	

	// first tessellation
	createInputForQhull();
    	callQhull();

	readTesse_finite("tess_info");

/*    
	std::ifstream ifs("tess_info");
	std::cout << "Read cell information begin!" << std::endl;
	if(!ifs) {
		std::cout << "stream error!" << std::endl; exit(-1);
	}
*/
	// write mesh file for abaqus
	std::ofstream ofs("mesh_for_Abaqus");
	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(OPREC);
	if(!ofs) {
		std::cout << "stream error!" << std::endl; exit(-1);
	}

	int m, n, i, j;
	int totalNum;

/*	ifs>>totalNum;
	ofs << totalNum << std::endl;
	for(int it=0; it!=totalNum; it++){
		ifs>>m>>n>>i>>j;	// read nodes in each line
		m = m+1;
		n = n+1;
		i = i+1;
		j = j+1;	// the ID from Qhull is starting from 0

		ofs << setw(OWID) << it+1 << ", "
		    << setw(OWID) << m << ", " 
		    << setw(OWID) << n << ", "
		    << setw(OWID) << i << ", "
		    << setw(OWID) << j << std::endl;
	}
	ifs.close();
*/
	int it_tet = 0;
	for(std::vector<cell*>::const_iterator iter=cellVec.begin(); iter!=cellVec.end(); iter++, it_tet++){
		m = (*iter)->getm();
		n = (*iter)->getn();
		i = (*iter)->geti();
		j = (*iter)->getj();

		ofs << setw(OWID) << it_tet+1 << ", "
		    << setw(OWID) << m << ", " 
		    << setw(OWID) << n << ", "
		    << setw(OWID) << i << ", "
		    << setw(OWID) << j << std::endl;

	
	}

	ofs.close();

	// write node information
	std::ofstream node_info("node_info_Abaqus");
  	if(!node_info) {
   		 cout << "stream error!" << endl; exit(-1);
  	}
 	node_info.setf(std::ios::scientific, std::ios::floatfield);
  	node_info.precision(OPREC);
  
  	vec tmp;
  	std::vector<particle*>::const_iterator  it;
	int ptclID = 1;
  	for (it=ParticleVec.begin();it!=ParticleVec.end();++it, ptclID++)  {
    
    		tmp=(*it)->getCurrPosition();
    		node_info << setw(OWID) << ptclID << ", "
		    	  << setw(OWID) << tmp.getx() << ", "
		    	  << setw(OWID) << tmp.gety() << ", "
		    	  << setw(OWID) << tmp.getz() << std::endl;
    
  	}

	REAL zero_num = 0;
  	node_info << setw(OWID) << ptclID+1 << ", "
		    	  << setw(OWID) << zero_num << ", "
		    	  << setw(OWID) << y_max << ", "
		    	  << setw(OWID) << z_max << std::endl;
	node_info << setw(OWID) << ptclID+2 << ", "
		    	  << setw(OWID) << x_max << ", "
		    	  << setw(OWID) << y_max << ", "
		    	  << setw(OWID) << z_max << std::endl;
	node_info << setw(OWID) << ptclID+3 << ", "
		    	  << setw(OWID) << x_max << ", "
		    	  << setw(OWID) << zero_num << ", "
		    	  << setw(OWID) << z_max << std::endl;
  	node_info.close();

}


// rotate the boundary walls along different axis. Test finite granular strain. May 21, 2013
void assembly::rotateXYZ(REAL angle,	// angle that you want to roate
	       		 vec axis_flag,	// represent axis you want to roate: 0 means not to this axis, 1 means along this axis. If axis_num is 1 1 1, then it means rotate along x,y,z
			 int   snapshots, 
			 int   interval,			  
			 const char* iniptclfile,   
			 const char* inibdryfile,
			 const char* particlefile, 
			 const char* boundaryfile,
			 const char* contactfile,  
			 const char* progressfile,
			 const char* balancedfile, 
			 const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
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
		<< setw(OWID) << "sample" // for new 3 sigma's
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
		<< setw(OWID) << "sample"
	        << setw(OWID) << "void"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "coordinate"
	        << setw(OWID) << "vibra"
	        << setw(OWID) << "impact"
	        << setw(OWID) << "wall-clock" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag"
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "Logarithmic" 	// for logarithmic strain
		<< setw(OWID) << "Logarithmic"
		<< setw(OWID) << "Logarithmic"
		<< setw(OWID) << "Logarithmic" 
		<< setw(OWID) << "Bagi" 	// for average_dudx_Bagi, used to test quadratic terms, April 22, 2013
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Lagrangian" 	// for average_dudx_Lagrangian, used to test quadratic terms
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian"
		<< setw(OWID) << "Eulerian" 	// for average_dudx_Eulerian, used to test quadratic terms
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" << endl
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
	        << setw(OWID) << "sigma_11"
	        << setw(OWID) << "sigma_12"
	        << setw(OWID) << "sigma_13"
	        << setw(OWID) << "sigma_21"
	        << setw(OWID) << "sigma_22"
	        << setw(OWID) << "sigma_23"
	        << setw(OWID) << "sigma_31"
	        << setw(OWID) << "sigma_32"
	        << setw(OWID) << "sigma_33"
	        << setw(OWID) << "mean_stress"
	        << setw(OWID) << "width"
	        << setw(OWID) << "length"
	        << setw(OWID) << "height"
	        << setw(OWID) << "volume"
	        << setw(OWID) << "epsilon_w"
	        << setw(OWID) << "epsilon_l"
	        << setw(OWID) << "epsilon_h"
	        << setw(OWID) << "epsilon_v"
	        << setw(OWID) << "epsilon_11"
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "ratio"
	        << setw(OWID) << "porosity"
	        << setw(OWID) << "number"
	        << setw(OWID) << "t_step"
	        << setw(OWID) << "t_step"
	        << setw(OWID) << "time" 
	        << setw(OWID) << "epsilon_11"	// for finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "epsilon_11"	// for mixed finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "epsilon_11"	// for eulerian finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33" 
	        << setw(OWID) << "epsilon_11"	// for mixed eulerian finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33" 
		<< setw(OWID) << "epsilon_w"	// for logranthmic strain
	        << setw(OWID) << "epsilon_l"
	        << setw(OWID) << "epsilon_h"
	        << setw(OWID) << "epsilon_v" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Bagi, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Lagrangian, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Eulerian, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" << endl;

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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    std::vector<REAL> temp_vec;	// used to test quadratic terms, April 22, 2013
    temp_vec.push_back(0);
    temp_vec.push_back(0);
    temp_vec.push_back(0);
    for(int i=0; i!=3; i++){
	average_dudx_Bagi.appendRow(temp_vec);
	average_dudx_Eulerian.appendRow(temp_vec);
	average_dudx_Lagrangian.appendRow(temp_vec);
    }

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    initVolume = W0*L0*H0;
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    matrix granularStress;
    matrix granularStrain;	// granular strain
    matrix previousStrain(3,3);
    matrix previousEuler(3,3);
    matrix finiteStrain(3,3);	// finite granular strain
    matrix lag_HOT(3,3);	// mixed finite granular strain
    matrix eulerianStrain(3,3);	// eulerian granular strain
    matrix euler_HOT(3,3);
    matrix previousEuler_HOT(3,3);
    matrix previousEuler_dudx(3,3);	// used to test quadratic terms, April 22, 2013
    matrix previousBagi_dudx(3,3);

    int total_steps;	// total steps will be determined by TIMESTEP and angle that you want to rotate
    total_steps = angle/TIMESTEP;
    std::cout << "total steps: " << total_steps << std::endl;


    // initialize previousStrain(3,3)
    for(int i_ps=0; i_ps!=2; i_ps++){
	for(int j_ps=0; j_ps!=2; j_ps++){
		previousStrain(i_ps+1,j_ps+1) = 0;
		previousEuler(i_ps+1,j_ps+1) = 0;
		previousEuler_HOT(i_ps+1,j_ps+1) = 0;
		previousEuler_dudx(i_ps+1,j_ps+1) = 0;	// used to test quadratic terms
		previousBagi_dudx(i_ps+1,j_ps+1) = 0;
	}
    }
    
    // first tessellation
    createInputForQhull();
    callQhull();
    readTesse("tess_info");
    readTesse_finite("tess_info");
    cellVec_init = cellVec;	// keep cellVec_init for lagrangian strain unchanged 

    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL epsilon_w_log, epsilon_l_log, epsilon_h_log;	// lograthmic strain
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
    g_iteration=0;
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
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;

	granularStress.clear();
	granularStress = getGranularStress();	// calculate granular stress

	void_ratio=Volume/getParticleVolume()-1;

	// control the rotate angle at every time step
	// initialize
	minctl[0].rote = TIMESTEP*vec(0,0,0);
	minctl[1].rote = TIMESTEP*vec(0,0,0);

	midctl[0].rote = TIMESTEP*vec(0,0,0);
	midctl[1].rote = TIMESTEP*vec(0,0,0);	

	maxctl[0].rote = TIMESTEP*vec(0,0,0);
	maxctl[1].rote = TIMESTEP*vec(0,0,0);	

	if(axis_flag.getx()){	// if x is 1, means rotate along x
		minctl[0].rote += TIMESTEP*vec(1.0,0,0);
		minctl[1].rote += TIMESTEP*vec(1.0,0,0);

		midctl[0].rote += TIMESTEP*vec(1.0,0,0);
		midctl[1].rote += TIMESTEP*vec(1.0,0,0);	

		maxctl[0].rote += TIMESTEP*vec(1.0,0,0);
		maxctl[1].rote += TIMESTEP*vec(1.0,0,0);	
	}

	if(axis_flag.gety()){	// if y is 1, means rotate along y
		minctl[0].rote += TIMESTEP*vec(0,1.0,0);
		minctl[1].rote += TIMESTEP*vec(0,1.0,0);

		midctl[0].rote += TIMESTEP*vec(0,1.0,0);
		midctl[1].rote += TIMESTEP*vec(0,1.0,0);	

		maxctl[0].rote += TIMESTEP*vec(0,1.0,0);
		maxctl[1].rote += TIMESTEP*vec(0,1.0,0);	
	}

	if(axis_flag.getz()){	// if x=z is 1, means rotate along z
		minctl[0].rote += TIMESTEP*vec(0,0,1.0);
		minctl[1].rote += TIMESTEP*vec(0,0,1.0);

		midctl[0].rote += TIMESTEP*vec(0,0,1.0);
		midctl[1].rote += TIMESTEP*vec(0,0,1.0);	

		maxctl[0].rote += TIMESTEP*vec(0,0,1.0);
		maxctl[1].rote += TIMESTEP*vec(0,0,1.0);	
	}

	updateRB(min,minctl,2);
	updateRB(mid,midctl,2);
	updateRB(max,maxctl,2);
	updateRB6();
	
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
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
	epsilon_w_log = log(l24/W0); epsilon_l_log = log(l13/L0); epsilon_h_log = log(l56/H0);
	if (g_iteration % interval == 0 ){
	    // calculate granular strain and Euler strain
	    average_dudx_Bagi.clear();	// used to test quadratic terms
	    average_dudx_Eulerian.clear();

	    granularStrain.clear();
	    eulerianStrain.clear();
	    euler_HOT.clear();
	    granularStrain = getGranularStrain()+previousStrain;
	    eulerianStrain = getEulerianStrain()+previousEuler;
	    euler_HOT = getEuler_HOT()+previousEuler_HOT;

	    average_dudx_Bagi = average_dudx_Bagi+previousBagi_dudx;	// used to test quadratic terms
	    average_dudx_Eulerian = average_dudx_Eulerian+previousEuler_dudx;

	    if(eulerianStrain(1,1)-previousEuler(1,1)>strainThreshold ||
	       eulerianStrain(2,2)-previousEuler(2,2)>strainThreshold ||
	       eulerianStrain(3,3)-previousEuler(3,3)>strainThreshold ){	// need to tessellate again
			resetStartPosition();	// reset initial position
			// tessellate again
			createInputForQhull();
			callQhull();
			readTesse("tess_info");
			readTesse_finite("tess_info");
			// update previousStrain
			previousStrain = granularStrain;
			previousEuler = eulerianStrain;
			previousEuler_HOT = euler_HOT;
			previousBagi_dudx = average_dudx_Bagi;	// used to test quadratic terms
			previousEuler_dudx = average_dudx_Eulerian;
	    }
	    // calculate finite granular strain
	    finiteStrain.clear();
	    average_dudx_Lagrangian.clear();	// used to test quadratic terms
	    finiteStrain = getFiniteStrain();
	    // calculate mixed finite granular strain
	    lag_HOT.clear();
	    lag_HOT = getHigherStrain();

	    gettimeofday(&time_w2,NULL);
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2 << setw(OWID) << sigma2_1
			<< setw(OWID) << sigma2_2 << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << granularStress(1,1) << setw(OWID) << granularStress(1,2) << setw(OWID) << granularStress(1,3)
		        << setw(OWID) << granularStress(2,1) << setw(OWID) << granularStress(2,2) << setw(OWID) << granularStress(2,3)
		        << setw(OWID) << granularStress(3,1) << setw(OWID) << granularStress(3,2) << setw(OWID) << granularStress(3,3)
		        << setw(OWID) << getAverageRigidPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
			<< setw(OWID) << granularStrain(1,1) << setw(OWID) << granularStrain(1,2) << setw(OWID) << granularStrain(1,3)
		        << setw(OWID) << granularStrain(2,1) << setw(OWID) << granularStrain(2,2) << setw(OWID) << granularStrain(2,3)
		        << setw(OWID) << granularStrain(3,1) << setw(OWID) << granularStrain(3,2) << setw(OWID) << granularStrain(3,3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
	                << setw(OWID) << getVibraTimeStep()
	                << setw(OWID) << getImpactTimeStep()
		        << setw(OWID) << timediffsec(time_w1,time_w2)
			<< setw(OWID) << finiteStrain(1,1) << setw(OWID) << finiteStrain(1,2) << setw(OWID) << finiteStrain(1,3)
		        << setw(OWID) << finiteStrain(2,1) << setw(OWID) << finiteStrain(2,2) << setw(OWID) << finiteStrain(2,3)
		        << setw(OWID) << finiteStrain(3,1) << setw(OWID) << finiteStrain(3,2) << setw(OWID) << finiteStrain(3,3)
			<< setw(OWID) << lag_HOT(1,1) << setw(OWID) << lag_HOT(1,2) << setw(OWID) << lag_HOT(1,3)
		        << setw(OWID) << lag_HOT(2,1) << setw(OWID) << lag_HOT(2,2) << setw(OWID) << lag_HOT(2,3)
		        << setw(OWID) << lag_HOT(3,1) << setw(OWID) << lag_HOT(3,2) << setw(OWID) << lag_HOT(3,3)
			<< setw(OWID) << eulerianStrain(1,1) << setw(OWID) << eulerianStrain(1,2) << setw(OWID) << eulerianStrain(1,3)
		        << setw(OWID) << eulerianStrain(2,1) << setw(OWID) << eulerianStrain(2,2) << setw(OWID) << eulerianStrain(2,3)
		        << setw(OWID) << eulerianStrain(3,1) << setw(OWID) << eulerianStrain(3,2) << setw(OWID) << eulerianStrain(3,3)
			<< setw(OWID) << euler_HOT(1,1) << setw(OWID) << euler_HOT(1,2) << setw(OWID) << euler_HOT(1,3)
		        << setw(OWID) << euler_HOT(2,1) << setw(OWID) << euler_HOT(2,2) << setw(OWID) << euler_HOT(2,3)
		        << setw(OWID) << euler_HOT(3,1) << setw(OWID) << euler_HOT(3,2) << setw(OWID) << euler_HOT(3,3)
		        << setw(OWID) << epsilon_w_log
		        << setw(OWID) << epsilon_l_log
		        << setw(OWID) << epsilon_h_log
		        << setw(OWID) << (epsilon_w_log+epsilon_l_log+epsilon_h_log)
			<< setw(OWID) << average_dudx_Bagi(1,1) << setw(OWID) << average_dudx_Bagi(1,2) << setw(OWID) << average_dudx_Bagi(1,3)
		        << setw(OWID) << average_dudx_Bagi(2,1) << setw(OWID) << average_dudx_Bagi(2,2) << setw(OWID) << average_dudx_Bagi(2,3)
		        << setw(OWID) << average_dudx_Bagi(3,1) << setw(OWID) << average_dudx_Bagi(3,2) << setw(OWID) << average_dudx_Bagi(3,3)
		        
			<< setw(OWID) << average_dudx_Lagrangian(1,1) << setw(OWID) << average_dudx_Lagrangian(1,2) << setw(OWID) << average_dudx_Lagrangian(1,3)
		        << setw(OWID) << average_dudx_Lagrangian(2,1) << setw(OWID) << average_dudx_Lagrangian(2,2) << setw(OWID) << average_dudx_Lagrangian(2,3)
		        << setw(OWID) << average_dudx_Lagrangian(3,1) << setw(OWID) << average_dudx_Lagrangian(3,2) << setw(OWID) << average_dudx_Lagrangian(3,3)
			<< setw(OWID) << average_dudx_Eulerian(1,1) << setw(OWID) << average_dudx_Eulerian(1,2) << setw(OWID) << average_dudx_Eulerian(1,3)
		        << setw(OWID) << average_dudx_Eulerian(2,1) << setw(OWID) << average_dudx_Eulerian(2,2) << setw(OWID) << average_dudx_Eulerian(2,3)
		        << setw(OWID) << average_dudx_Eulerian(3,1) << setw(OWID) << average_dudx_Eulerian(3,2) << setw(OWID) << average_dudx_Eulerian(3,3)
		        << endl;
	    g_debuginf << setw(OWID) << g_iteration
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
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    g_timeinf << setw(4) << "end" << " " << ctime(&timeStamp) << flush;

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}




// rotate the each particle with the same rotation tensor directly. Test finite granular strain. June 10, 2013
void assembly::rotateXYZ_dir(REAL angle,	// angle that you want to roate
	       		 vec axis_flag,	// represent axis you want to roate: 0 means not to this axis, 1 means along this axis. If axis_num is 1 1 1, then it means rotate along x,y,z
			 int   snapshots, 
			 int   interval,			  
			 const char* iniptclfile,   
			 const char* inibdryfile,
			 const char* particlefile, 
			 const char* boundaryfile,
			 const char* contactfile,  
			 const char* progressfile,
			 const char* balancedfile, 
			 const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
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
		<< setw(OWID) << "sample" // for new 3 sigma's
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
		<< setw(OWID) << "sample"
	        << setw(OWID) << "void"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "coordinate"
	        << setw(OWID) << "vibra"
	        << setw(OWID) << "impact"
	        << setw(OWID) << "wall-clock" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag"
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "Logarithmic" 	// for logarithmic strain
		<< setw(OWID) << "Logarithmic"
		<< setw(OWID) << "Logarithmic"
		<< setw(OWID) << "Logarithmic" 
		<< setw(OWID) << "Bagi" 	// for average_dudx_Bagi, used to test quadratic terms, April 22, 2013
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Lagrangian" 	// for average_dudx_Lagrangian, used to test quadratic terms
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian"
		<< setw(OWID) << "Eulerian" 	// for average_dudx_Eulerian, used to test quadratic terms
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "spatial" 	// for spatial velocity gradient tensor, deformation rate tensor
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial"
		<< setw(OWID) << "from rate" 	// for strain based on spatial deformation rate tensor, July 8, 2013
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate"  << endl
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
	        << setw(OWID) << "sigma_11"
	        << setw(OWID) << "sigma_12"
	        << setw(OWID) << "sigma_13"
	        << setw(OWID) << "sigma_21"
	        << setw(OWID) << "sigma_22"
	        << setw(OWID) << "sigma_23"
	        << setw(OWID) << "sigma_31"
	        << setw(OWID) << "sigma_32"
	        << setw(OWID) << "sigma_33"
	        << setw(OWID) << "mean_stress"
	        << setw(OWID) << "width"
	        << setw(OWID) << "length"
	        << setw(OWID) << "height"
	        << setw(OWID) << "volume"
	        << setw(OWID) << "epsilon_w"
	        << setw(OWID) << "epsilon_l"
	        << setw(OWID) << "epsilon_h"
	        << setw(OWID) << "epsilon_v"
	        << setw(OWID) << "epsilon_11"
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "ratio"
	        << setw(OWID) << "porosity"
	        << setw(OWID) << "number"
	        << setw(OWID) << "t_step"
	        << setw(OWID) << "t_step"
	        << setw(OWID) << "time" 
	        << setw(OWID) << "epsilon_11"	// for finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "epsilon_11"	// for mixed finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "epsilon_11"	// for eulerian finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33" 
	        << setw(OWID) << "epsilon_11"	// for mixed eulerian finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33" 
		<< setw(OWID) << "epsilon_w"	// for logranthmic strain
	        << setw(OWID) << "epsilon_l"
	        << setw(OWID) << "epsilon_h"
	        << setw(OWID) << "epsilon_v" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Bagi, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Lagrangian, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Eulerian, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" 
		<< setw(OWID) << "dvdx_11"	// for spatial velocity gradient tensor, deformation rate tensor
	        << setw(OWID) << "dvdx_12"
	        << setw(OWID) << "dvdx_13"
	        << setw(OWID) << "dvdx_21"
	        << setw(OWID) << "dvdx_22"
	        << setw(OWID) << "dvdx_23"
	        << setw(OWID) << "dvdx_31"
	        << setw(OWID) << "dvdx_32"
	        << setw(OWID) << "dvdx_33" 
	        << setw(OWID) << "epsilon_11"	// for strain based on spatial deformation rate tensor
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33" << endl;

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

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    std::vector<REAL> temp_vec;	// used to test quadratic terms, April 22, 2013
    temp_vec.push_back(0);
    temp_vec.push_back(0);
    temp_vec.push_back(0);
    for(int i=0; i!=3; i++){
	average_dudx_Bagi.appendRow(temp_vec);
	average_dudx_Eulerian.appendRow(temp_vec);
	average_dudx_Lagrangian.appendRow(temp_vec);
    }

    // pre_2. create particles and boundaries from files
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    initVolume = W0*L0*H0;
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    matrix granularStress;
    matrix granularStrain;	// granular strain
    matrix previousStrain(3,3);
    matrix previousEuler(3,3);
    matrix finiteStrain(3,3);	// finite granular strain
    matrix lag_HOT(3,3);	// mixed finite granular strain
    matrix eulerianStrain(3,3);	// eulerian granular strain
    matrix euler_HOT(3,3);
    matrix previousEuler_HOT(3,3);
    matrix previousEuler_dudx(3,3);	// used to test quadratic terms, April 22, 2013
    matrix previousBagi_dudx(3,3);
    matrix spatial_dvdx(3,3);	// average spatial velocity gradient tensor

    matrix initialPosition(3,1);	// initial position of each particle, used for rotation
    matrix currentPosition(3,1);	// current position of each particle, used for rotation
    vec con_position;
    vec curr_velocity;
    REAL prev_x, prev_y, prev_z;

    matrix curr_rate(3,3);	// current spatial deformation rate tensor, July 8, 2013
    matrix curr_dvdx(3,3);	// current spatial velocity gradient tensor
    matrix prev_strain_rate(3,3);
    matrix curr_strain_rate(3,3);	// current strain based on deformation rate tensor


    // rotation tensor
    matrix Rx(3,3);	// rotation tensor in x direction
    matrix Ry(3,3);	// rotation tensor in y direction
    matrix Rz(3,3);	// rotation tensor in z direction

    // initial rotation tensor into identity matrix
    Rx(1,1) = 1;
    Rx(1,2) = 0;
    Rx(1,3) = 0;
    Rx(2,1) = 0;
    Rx(2,2) = 1;
    Rx(2,3) = 0;
    Rx(3,1) = 0;
    Rx(3,2) = 0;
    Rx(3,3) = 1;

    Ry(1,1) = 1;
    Ry(1,2) = 0;
    Ry(1,3) = 0;
    Ry(2,1) = 0;
    Ry(2,2) = 1;
    Ry(2,3) = 0;
    Ry(3,1) = 0;
    Ry(3,2) = 0;
    Ry(3,3) = 1;

    Rz(1,1) = 1;
    Rz(1,2) = 0;
    Rz(1,3) = 0;
    Rz(2,1) = 0;
    Rz(2,2) = 1;
    Rz(2,3) = 0;
    Rz(3,1) = 0;
    Rz(3,2) = 0;
    Rz(3,3) = 1;

    long double theta_rot;	// current rotation angle

    int total_steps;	// total steps will be determined by TIMESTEP and angle that you want to rotate
    total_steps = angle/TIMESTEP;
    std::cout << "total steps: " << total_steps << std::endl;


    // initialize previousStrain(3,3)
    for(int i_ps=0; i_ps!=2; i_ps++){
	for(int j_ps=0; j_ps!=2; j_ps++){
		previousStrain(i_ps+1,j_ps+1) = 0;
		previousEuler(i_ps+1,j_ps+1) = 0;
		previousEuler_HOT(i_ps+1,j_ps+1) = 0;
		previousEuler_dudx(i_ps+1,j_ps+1) = 0;	// used to test quadratic terms
		previousBagi_dudx(i_ps+1,j_ps+1) = 0;
		curr_strain_rate(i_ps+1,j_ps+1) = 0;	// calculate strain by rate, July 8, 2013
	}
    }
    
    // first tessellation
    createInputForQhull();
    callQhull();
    readTesse("tess_info");
    readTesse_finite("tess_info");
    cellVec_init = cellVec;	// keep cellVec_init for lagrangian strain unchanged 

    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL epsilon_w_log, epsilon_l_log, epsilon_h_log;	// lograthmic strain
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
    int  bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    // iterations start here...
    g_iteration=0;
    gettimeofday(&time_w1,NULL);
    do
    {
	theta_rot = (g_iteration+1)*TIMESTEP;
	// calculate rotaion tensor
	if(axis_flag.getx()){	// if x is 1, means rotate along x
		Rx(2,2) = cos(theta_rot);
		Rx(2,3) = -sin(theta_rot);
		Rx(3,2) = sin(theta_rot);
		Rx(3,3) = cos(theta_rot);	
	}

	if(axis_flag.gety()){	// if y is 1, means rotate along y
		Ry(1,1) = cos(theta_rot);
		Ry(1,3) = sin(theta_rot);
		Ry(3,1) = -sin(theta_rot);
		Ry(3,3) = cos(theta_rot);	
	}

	if(axis_flag.getz()){	// if x=z is 1, means rotate along z
		Rz(1,1) = cos(theta_rot);
		Rz(1,2) = -sin(theta_rot);
		Rz(2,1) = sin(theta_rot);
		Rz(2,2)	= cos(theta_rot);
	}

	// apply rotation to each particle
	for(std::vector<particle*>::iterator it=ParticleVec.begin();it!=ParticleVec.end();++it){				
		// get previous position
		prev_x = (*it)->getCurrPosition().getx();
		prev_y = (*it)->getCurrPosition().gety();
		prev_z = (*it)->getCurrPosition().getz();

		initialPosition(1,1) = (*it)->getInitPosition().getx();
		initialPosition(2,1) = (*it)->getInitPosition().gety();
		initialPosition(3,1) = (*it)->getInitPosition().getz();
		
		currentPosition = Rz*Ry*Rx*initialPosition;	// rotation

		con_position.setx(currentPosition(1,1));
		con_position.sety(currentPosition(2,1));
		con_position.setz(currentPosition(3,1));
		(*it)->setCurrPosition(con_position);
		
		// set current velocity
		curr_velocity.setx((currentPosition(1,1)-prev_x)/TIMESTEP);
		curr_velocity.sety((currentPosition(2,1)-prev_y)/TIMESTEP);
		curr_velocity.setz((currentPosition(3,1)-prev_z)/TIMESTEP);
		(*it)->setCurrVelocity(curr_velocity);
    	}
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;

	granularStress.clear();
	granularStress = getGranularStress();	// calculate granular stress

	void_ratio=Volume/getParticleVolume()-1;

	// control the rotate angle at every time step
	// initialize
	minctl[0].rote = TIMESTEP*vec(0,0,0);
	minctl[1].rote = TIMESTEP*vec(0,0,0);

	midctl[0].rote = TIMESTEP*vec(0,0,0);
	midctl[1].rote = TIMESTEP*vec(0,0,0);	

	maxctl[0].rote = TIMESTEP*vec(0,0,0);
	maxctl[1].rote = TIMESTEP*vec(0,0,0);	

	if(axis_flag.getx()){	// if x is 1, means rotate along x
		minctl[0].rote += TIMESTEP*vec(1.0,0,0);
		minctl[1].rote += TIMESTEP*vec(1.0,0,0);

		midctl[0].rote += TIMESTEP*vec(1.0,0,0);
		midctl[1].rote += TIMESTEP*vec(1.0,0,0);	

		maxctl[0].rote += TIMESTEP*vec(1.0,0,0);
		maxctl[1].rote += TIMESTEP*vec(1.0,0,0);	
	}

	if(axis_flag.gety()){	// if y is 1, means rotate along y
		minctl[0].rote += TIMESTEP*vec(0,1.0,0);
		minctl[1].rote += TIMESTEP*vec(0,1.0,0);

		midctl[0].rote += TIMESTEP*vec(0,1.0,0);
		midctl[1].rote += TIMESTEP*vec(0,1.0,0);	

		maxctl[0].rote += TIMESTEP*vec(0,1.0,0);
		maxctl[1].rote += TIMESTEP*vec(0,1.0,0);	
	}

	if(axis_flag.getz()){	// if x=z is 1, means rotate along z
		minctl[0].rote += TIMESTEP*vec(0,0,1.0);
		minctl[1].rote += TIMESTEP*vec(0,0,1.0);

		midctl[0].rote += TIMESTEP*vec(0,0,1.0);
		midctl[1].rote += TIMESTEP*vec(0,0,1.0);	

		maxctl[0].rote += TIMESTEP*vec(0,0,1.0);
		maxctl[1].rote += TIMESTEP*vec(0,0,1.0);	
	}

	updateRB(min,minctl,2);
	updateRB(mid,midctl,2);
	updateRB(max,maxctl,2);
	updateRB6();
	
	// calculate strain by rate
	prev_strain_rate.clear();
	prev_strain_rate = curr_strain_rate;	// set previous strain
	
	// calculate spatial velocity gradient tensor
	spatial_dvdx.clear();
	spatial_dvdx = getAverage_dvdx();
	curr_dvdx.clear();
	curr_dvdx = spatial_dvdx;
	
	curr_rate = 0.5*(curr_dvdx+curr_dvdx.getTrans());	// current spatial deformation rate

	curr_strain_rate = prev_strain_rate+(curr_rate-curr_dvdx.getTrans()*prev_strain_rate-prev_strain_rate*curr_dvdx)*TIMESTEP;

	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
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
	epsilon_w_log = log(l24/W0); epsilon_l_log = log(l13/L0); epsilon_h_log = log(l56/H0);
	if (g_iteration % interval == 0 ){


//	    if(eulerianStrain(1,1)-previousEuler(1,1)>strainThreshold ||
//	       eulerianStrain(2,2)-previousEuler(2,2)>strainThreshold ||
//	       eulerianStrain(3,3)-previousEuler(3,3)>strainThreshold ){	// need to tessellate again
//	    if (g_iteration % interval == 0){	// control tessellation in another way, for test rotation. Since during rotation, Eulerian finite strain will be zero, to test the tessellation, we need to update tessellation then use this method to control update

			// tessellate again
//			createInputForQhull();
//			callQhull();
//			readTesse("tess_info");
//			readTesse_finite("tess_info");

//	    }


	    // calculate granular strain and Euler strain
	    average_dudx_Bagi.clear();	// used to test quadratic terms
	    average_dudx_Eulerian.clear();

	    granularStrain.clear();
	    eulerianStrain.clear();
	    euler_HOT.clear();
//	    granularStrain = getGranularStrain()+previousStrain;
//	    eulerianStrain = getEulerianStrain()+previousEuler;
//	    euler_HOT = getEuler_HOT()+previousEuler_HOT;
//
//	    average_dudx_Bagi = average_dudx_Bagi+previousBagi_dudx;	// used to test quadratic terms
//	    average_dudx_Eulerian = average_dudx_Eulerian+previousEuler_dudx;

	    // calculate granular strain based on intial position directly, instead of previous strain. July 1, 2013
	    granularStrain = getGranularStrain();
	    eulerianStrain = getEulerianStrain();
	    euler_HOT = getEuler_HOT();

	    average_dudx_Bagi = average_dudx_Bagi;	// used to test quadratic terms
	    average_dudx_Eulerian = average_dudx_Eulerian;

//	    // calculate spatial velocity gradient tensor
//	    spatial_dvdx.clear();
//	    spatial_dvdx = getAverage_dvdx();

//	    resetStartPosition();	// reset initial position
	    // update previousStrain
	    previousStrain = granularStrain;
	    previousEuler = eulerianStrain;
	    previousEuler_HOT = euler_HOT;
	    previousBagi_dudx = average_dudx_Bagi;	// used to test quadratic terms
	    previousEuler_dudx = average_dudx_Eulerian;


	    // calculate finite granular strain
	    finiteStrain.clear();
	    average_dudx_Lagrangian.clear();	// used to test quadratic terms
	    finiteStrain = getFiniteStrain();
	    // calculate mixed finite granular strain
	    lag_HOT.clear();
	    lag_HOT = getHigherStrain();

	    gettimeofday(&time_w2,NULL);
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << getDensity()
		        << setw(OWID) << sigma1_1 << setw(OWID) << sigma1_2 << setw(OWID) << sigma2_1
			<< setw(OWID) << sigma2_2 << setw(OWID) << sigma3_1 << setw(OWID) << sigma3_2
		        << setw(OWID) << granularStress(1,1) << setw(OWID) << granularStress(1,2) << setw(OWID) << granularStress(1,3)
		        << setw(OWID) << granularStress(2,1) << setw(OWID) << granularStress(2,2) << setw(OWID) << granularStress(2,3)
		        << setw(OWID) << granularStress(3,1) << setw(OWID) << granularStress(3,2) << setw(OWID) << granularStress(3,3)
		        << setw(OWID) << getAverageRigidPressure()
		        << setw(OWID) << l24 << setw(OWID) << l13 << setw(OWID) << l56
		        << setw(OWID) << Volume
		        << setw(OWID) << epsilon_w
		        << setw(OWID) << epsilon_l
		        << setw(OWID) << epsilon_h
		        << setw(OWID) << (epsilon_w+epsilon_l+epsilon_h)
			<< setw(OWID) << granularStrain(1,1) << setw(OWID) << granularStrain(1,2) << setw(OWID) << granularStrain(1,3)
		        << setw(OWID) << granularStrain(2,1) << setw(OWID) << granularStrain(2,2) << setw(OWID) << granularStrain(2,3)
		        << setw(OWID) << granularStrain(3,1) << setw(OWID) << granularStrain(3,2) << setw(OWID) << granularStrain(3,3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
	                << setw(OWID) << getVibraTimeStep()
	                << setw(OWID) << getImpactTimeStep()
		        << setw(OWID) << timediffsec(time_w1,time_w2)
			<< setw(OWID) << finiteStrain(1,1) << setw(OWID) << finiteStrain(1,2) << setw(OWID) << finiteStrain(1,3)
		        << setw(OWID) << finiteStrain(2,1) << setw(OWID) << finiteStrain(2,2) << setw(OWID) << finiteStrain(2,3)
		        << setw(OWID) << finiteStrain(3,1) << setw(OWID) << finiteStrain(3,2) << setw(OWID) << finiteStrain(3,3)
			<< setw(OWID) << lag_HOT(1,1) << setw(OWID) << lag_HOT(1,2) << setw(OWID) << lag_HOT(1,3)
		        << setw(OWID) << lag_HOT(2,1) << setw(OWID) << lag_HOT(2,2) << setw(OWID) << lag_HOT(2,3)
		        << setw(OWID) << lag_HOT(3,1) << setw(OWID) << lag_HOT(3,2) << setw(OWID) << lag_HOT(3,3)
			<< setw(OWID) << eulerianStrain(1,1) << setw(OWID) << eulerianStrain(1,2) << setw(OWID) << eulerianStrain(1,3)
		        << setw(OWID) << eulerianStrain(2,1) << setw(OWID) << eulerianStrain(2,2) << setw(OWID) << eulerianStrain(2,3)
		        << setw(OWID) << eulerianStrain(3,1) << setw(OWID) << eulerianStrain(3,2) << setw(OWID) << eulerianStrain(3,3)
			<< setw(OWID) << euler_HOT(1,1) << setw(OWID) << euler_HOT(1,2) << setw(OWID) << euler_HOT(1,3)
		        << setw(OWID) << euler_HOT(2,1) << setw(OWID) << euler_HOT(2,2) << setw(OWID) << euler_HOT(2,3)
		        << setw(OWID) << euler_HOT(3,1) << setw(OWID) << euler_HOT(3,2) << setw(OWID) << euler_HOT(3,3)
		        << setw(OWID) << epsilon_w_log
		        << setw(OWID) << epsilon_l_log
		        << setw(OWID) << epsilon_h_log
		        << setw(OWID) << (epsilon_w_log+epsilon_l_log+epsilon_h_log)
			<< setw(OWID) << average_dudx_Bagi(1,1) << setw(OWID) << average_dudx_Bagi(1,2) << setw(OWID) << average_dudx_Bagi(1,3)
		        << setw(OWID) << average_dudx_Bagi(2,1) << setw(OWID) << average_dudx_Bagi(2,2) << setw(OWID) << average_dudx_Bagi(2,3)
		        << setw(OWID) << average_dudx_Bagi(3,1) << setw(OWID) << average_dudx_Bagi(3,2) << setw(OWID) << average_dudx_Bagi(3,3)
		        
			<< setw(OWID) << average_dudx_Lagrangian(1,1) << setw(OWID) << average_dudx_Lagrangian(1,2) << setw(OWID) << average_dudx_Lagrangian(1,3)
		        << setw(OWID) << average_dudx_Lagrangian(2,1) << setw(OWID) << average_dudx_Lagrangian(2,2) << setw(OWID) << average_dudx_Lagrangian(2,3)
		        << setw(OWID) << average_dudx_Lagrangian(3,1) << setw(OWID) << average_dudx_Lagrangian(3,2) << setw(OWID) << average_dudx_Lagrangian(3,3)
			<< setw(OWID) << average_dudx_Eulerian(1,1) << setw(OWID) << average_dudx_Eulerian(1,2) << setw(OWID) << average_dudx_Eulerian(1,3)
		        << setw(OWID) << average_dudx_Eulerian(2,1) << setw(OWID) << average_dudx_Eulerian(2,2) << setw(OWID) << average_dudx_Eulerian(2,3)
		        << setw(OWID) << average_dudx_Eulerian(3,1) << setw(OWID) << average_dudx_Eulerian(3,2) << setw(OWID) << average_dudx_Eulerian(3,3)
			<< setw(OWID) << spatial_dvdx(1,1) << setw(OWID) << spatial_dvdx(1,2) << setw(OWID) << spatial_dvdx(1,3)
		        << setw(OWID) << spatial_dvdx(2,1) << setw(OWID) << spatial_dvdx(2,2) << setw(OWID) << spatial_dvdx(2,3)
		        << setw(OWID) << spatial_dvdx(3,1) << setw(OWID) << spatial_dvdx(3,2) << setw(OWID) << spatial_dvdx(3,3) 
			<< setw(OWID) << curr_strain_rate(1,1) << setw(OWID) << curr_strain_rate(1,2) << setw(OWID) << curr_strain_rate(1,3)
		        << setw(OWID) << curr_strain_rate(2,1) << setw(OWID) << curr_strain_rate(2,2) << setw(OWID) << curr_strain_rate(2,3)
		        << setw(OWID) << curr_strain_rate(3,1) << setw(OWID) << curr_strain_rate(3,2) << setw(OWID) << curr_strain_rate(3,3) 
		        << endl;
/*
	    g_debuginf << setw(OWID) << g_iteration
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
*/
	}
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    g_timeinf << setw(4) << "end" << " " << ctime(&timeStamp) << flush;

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}

  // this is actually a deposit process, but since we cannot calculate granular stress and strain during normal deposition, however in cavity expansion we need to calculate granular stress and strain. Then we write a new function cavityExpand to do cavity expansion and in this function, we will calculate granular stress and strain. June 17, 2013
// the output progress file format is the same as triaxial
void assembly::cavityExpand(int   total_steps,  
		       int   snapshots,
		       int   interval,
		       const char* iniptclfile,   
		       const char* inibdryfile,
		       const char* particlefile, 
		       const char* contactfile,
		       const char* progressfile, 
		       const char* debugfile)
{
    // pre_1: open streams for output.
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
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
		<< setw(OWID) << "sample" // for new 3 sigma's
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
		<< setw(OWID) << "sample"
	        << setw(OWID) << "void"
	        << setw(OWID) << "sample"
	        << setw(OWID) << "coordinate"
	        << setw(OWID) << "vibra"
	        << setw(OWID) << "impact"
	        << setw(OWID) << "wall-clock" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "finite" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag" 
		<< setw(OWID) << "H.O.T._Lag"
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "H.O.T._Euler" 
		<< setw(OWID) << "Logarithmic" 	// for logarithmic strain
		<< setw(OWID) << "Logarithmic"
		<< setw(OWID) << "Logarithmic"
		<< setw(OWID) << "Logarithmic" 
		<< setw(OWID) << "Bagi" 	// for average_dudx_Bagi, used to test quadratic terms, April 22, 2013
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Bagi" 
		<< setw(OWID) << "Lagrangian" 	// for average_dudx_Lagrangian, used to test quadratic terms
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian" 
		<< setw(OWID) << "Lagrangian"
		<< setw(OWID) << "Eulerian" 	// for average_dudx_Eulerian, used to test quadratic terms
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "Eulerian" 
		<< setw(OWID) << "spatial" 	// for spatial velocity gradient tensor, deformation rate tensor
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial"
		<< setw(OWID) << "total volume" << endl
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
	        << setw(OWID) << "sigma_11"
	        << setw(OWID) << "sigma_12"
	        << setw(OWID) << "sigma_13"
	        << setw(OWID) << "sigma_21"
	        << setw(OWID) << "sigma_22"
	        << setw(OWID) << "sigma_23"
	        << setw(OWID) << "sigma_31"
	        << setw(OWID) << "sigma_32"
	        << setw(OWID) << "sigma_33"
	        << setw(OWID) << "mean_stress"
	        << setw(OWID) << "width"
	        << setw(OWID) << "length"
	        << setw(OWID) << "height"
	        << setw(OWID) << "volume"
	        << setw(OWID) << "epsilon_w"
	        << setw(OWID) << "epsilon_l"
	        << setw(OWID) << "epsilon_h"
	        << setw(OWID) << "epsilon_v"
	        << setw(OWID) << "epsilon_11"
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "ratio"
	        << setw(OWID) << "porosity"
	        << setw(OWID) << "number"
	        << setw(OWID) << "t_step"
	        << setw(OWID) << "t_step"
	        << setw(OWID) << "time" 
	        << setw(OWID) << "epsilon_11"	// for finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "epsilon_11"	// for mixed finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
	        << setw(OWID) << "epsilon_11"	// for eulerian finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33" 
	        << setw(OWID) << "epsilon_11"	// for mixed eulerian finite strain
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33" 
		<< setw(OWID) << "epsilon_w"	// for logranthmic strain
	        << setw(OWID) << "epsilon_l"
	        << setw(OWID) << "epsilon_h"
	        << setw(OWID) << "epsilon_v" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Bagi, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Lagrangian, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" 
	        << setw(OWID) << "dudx_11"	// for average_dudx_Eulerian, used to test quadratic terms
	        << setw(OWID) << "dudx_12"
	        << setw(OWID) << "dudx_13"
	        << setw(OWID) << "dudx_21"
	        << setw(OWID) << "dudx_22"
	        << setw(OWID) << "dudx_23"
	        << setw(OWID) << "dudx_31"
	        << setw(OWID) << "dudx_32"
	        << setw(OWID) << "dudx_33" 	        
		<< setw(OWID) << "dvdx_11"	// for spatial velocity gradient tensor, deformation rate tensor
	        << setw(OWID) << "dvdx_12"
	        << setw(OWID) << "dvdx_13"
	        << setw(OWID) << "dvdx_21"
	        << setw(OWID) << "dvdx_22"
	        << setw(OWID) << "dvdx_23"
	        << setw(OWID) << "dvdx_31"
	        << setw(OWID) << "dvdx_32"
	        << setw(OWID) << "dvdx_33"
		<< setw(OWID) << "of cells" << endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1); }
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);


    std::vector<REAL> temp_vec;	// used to test quadratic terms, April 22, 2013
    temp_vec.push_back(0);
    temp_vec.push_back(0);
    temp_vec.push_back(0);
    for(int i=0; i!=3; i++){
	average_dudx_Bagi.appendRow(temp_vec);
	average_dudx_Eulerian.appendRow(temp_vec);
	average_dudx_Lagrangian.appendRow(temp_vec);
    }


    // pre_2. create particles and boundaries from existing files.
    readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    readBoundary(inibdryfile);   // create boundaries.

    // pre_3: define variables used in iterations.
    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    initVolume = W0*L0*H0;
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    matrix granularStress(3,3);
    matrix granularStrain(3,3);	// granular strain
    matrix previousStrain(3,3);
    matrix previousEuler(3,3);
    matrix finiteStrain(3,3);	// finite granular strain
    matrix lag_HOT(3,3);	// mixed finite granular strain
    matrix eulerianStrain(3,3);	// eulerian granular strain
    matrix euler_HOT(3,3);
    matrix previousEuler_HOT(3,3);
    matrix previousEuler_dudx(3,3);	// used to test quadratic terms, April 22, 2013
    matrix previousBagi_dudx(3,3);
    matrix spatial_dvdx(3,3);	// average spatial velocity gradient tensor

    // initialize previousStrain(3,3)
    for(int i_ps=0; i_ps!=2; i_ps++){
	for(int j_ps=0; j_ps!=2; j_ps++){
		previousStrain(i_ps+1,j_ps+1) = 0;
		previousEuler(i_ps+1,j_ps+1) = 0;
		previousEuler_HOT(i_ps+1,j_ps+1) = 0;
		previousEuler_dudx(i_ps+1,j_ps+1) = 0;	// used to test quadratic terms
		previousBagi_dudx(i_ps+1,j_ps+1) = 0;
	}
    }
    
    // first tessellation
    createInputForQhull();
    callQhull();
    readTesse("tess_info");
    readTesse_finite("tess_info");
    cellVec_init = cellVec;	// keep cellVec_init for lagrangian strain unchanged 

    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL epsilon_w_log, epsilon_l_log, epsilon_h_log;	// lograthmic strain


    REAL avgNormal=0;
    REAL avgTangt=0;
    int  stepsnum=0;
    char stepsstr[4];
    char stepsfp[50];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
//    matrix granularStress;	// granular stress
    int  bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    // iterations starting ...
    g_iteration=0; 
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
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();
	//gettimeofday(&time_w2,NULL);
	//cout << setw(OWID) << timediffsec(time_w1,time_w2) << endl;

	// 6. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx(); Volume=l13*l24*l56;
	
	granularStress.clear();
	granularStress = getGranularStress();	// calculate granular stress

	void_ratio=Volume/getParticleVolume()-1;

	// 7. (1) output particles and contacts information as snapshots.
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    time(&timeStamp);
	    g_timeinf << setw(4) << stepsnum << " " << ctime(&timeStamp) << flush;
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info.
	if (g_iteration % interval == 0) {
	    // calculate granular strain and Euler strain
	    average_dudx_Bagi.clear();	// used to test quadratic terms
	    average_dudx_Eulerian.clear();

	    granularStrain.clear();
	    eulerianStrain.clear();
	    euler_HOT.clear();
	    granularStrain = getGranularStrain()+previousStrain;
	    eulerianStrain = getEulerianStrain()+previousEuler;
	    euler_HOT = getEuler_HOT()+previousEuler_HOT;

	    average_dudx_Bagi = average_dudx_Bagi+previousBagi_dudx;	// used to test quadratic terms
	    average_dudx_Eulerian = average_dudx_Eulerian+previousEuler_dudx;

	    if(eulerianStrain(1,1)-previousEuler(1,1)>strainThreshold ||
	       eulerianStrain(2,2)-previousEuler(2,2)>strainThreshold ||
	       eulerianStrain(3,3)-previousEuler(3,3)>strainThreshold ){	// need to tessellate again
			resetStartPosition();	// reset initial position
			// tessellate again
			createInputForQhull();
			callQhull();
			readTesse("tess_info");
			readTesse_finite("tess_info");
			// update previousStrain
			previousStrain = granularStrain;
			previousEuler = eulerianStrain;
			previousEuler_HOT = euler_HOT;
			previousBagi_dudx = average_dudx_Bagi;	// used to test quadratic terms
			previousEuler_dudx = average_dudx_Eulerian;
	    }

	    // calculate spatial velocity gradient tensor
	    spatial_dvdx.clear();
	    spatial_dvdx = getAverage_dvdx();

	    // calculate finite granular strain
	    finiteStrain.clear();
	    average_dudx_Lagrangian.clear();	// used to test quadratic terms
	    finiteStrain = getFiniteStrain();
	    // calculate mixed finite granular strain
	    lag_HOT.clear();
	    lag_HOT = getHigherStrain();

	    gettimeofday(&time_w2,NULL);
	    progressinf << setw(OWID) << g_iteration
		        << setw(OWID) << getPossCntctNum()
		        << setw(OWID) << getActualCntctNum()
		        << setw(OWID) << getAveragePenetration()
		        << setw(OWID) << avgNormal
		        << setw(OWID) << avgTangt
		        << setw(OWID) << getAverageVelocity() 
		        << setw(OWID) << getAverageOmga()
		        << setw(OWID) << getAverageForce()   
		        << setw(OWID) << getAverageMoment()
		        << setw(OWID) << 0
		        << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 0
			<< setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 0
		        << setw(OWID) << granularStress(1,1) << setw(OWID) << granularStress(1,2) << setw(OWID) << granularStress(1,3)
		        << setw(OWID) << granularStress(2,1) << setw(OWID) << granularStress(2,2) << setw(OWID) << granularStress(2,3)
		        << setw(OWID) << granularStress(3,1) << setw(OWID) << granularStress(3,2) << setw(OWID) << granularStress(3,3)
		        << setw(OWID) << getAverageRigidPressure()
		        << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 0
			<< setw(OWID) << granularStrain(1,1) << setw(OWID) << granularStrain(1,2) << setw(OWID) << granularStrain(1,3)
		        << setw(OWID) << granularStrain(2,1) << setw(OWID) << granularStrain(2,2) << setw(OWID) << granularStrain(2,3)
		        << setw(OWID) << granularStrain(3,1) << setw(OWID) << granularStrain(3,2) << setw(OWID) << granularStrain(3,3)
		        << setw(OWID) << void_ratio
		        << setw(OWID) << void_ratio/(1+void_ratio)
		        << setw(OWID) << 2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
	                << setw(OWID) << getVibraTimeStep()
	                << setw(OWID) << getImpactTimeStep()
		        << setw(OWID) << timediffsec(time_w1,time_w2)
			<< setw(OWID) << finiteStrain(1,1) << setw(OWID) << finiteStrain(1,2) << setw(OWID) << finiteStrain(1,3)
		        << setw(OWID) << finiteStrain(2,1) << setw(OWID) << finiteStrain(2,2) << setw(OWID) << finiteStrain(2,3)
		        << setw(OWID) << finiteStrain(3,1) << setw(OWID) << finiteStrain(3,2) << setw(OWID) << finiteStrain(3,3)
			<< setw(OWID) << lag_HOT(1,1) << setw(OWID) << lag_HOT(1,2) << setw(OWID) << lag_HOT(1,3)
		        << setw(OWID) << lag_HOT(2,1) << setw(OWID) << lag_HOT(2,2) << setw(OWID) << lag_HOT(2,3)
		        << setw(OWID) << lag_HOT(3,1) << setw(OWID) << lag_HOT(3,2) << setw(OWID) << lag_HOT(3,3)
			<< setw(OWID) << eulerianStrain(1,1) << setw(OWID) << eulerianStrain(1,2) << setw(OWID) << eulerianStrain(1,3)
		        << setw(OWID) << eulerianStrain(2,1) << setw(OWID) << eulerianStrain(2,2) << setw(OWID) << eulerianStrain(2,3)
		        << setw(OWID) << eulerianStrain(3,1) << setw(OWID) << eulerianStrain(3,2) << setw(OWID) << eulerianStrain(3,3)
			<< setw(OWID) << euler_HOT(1,1) << setw(OWID) << euler_HOT(1,2) << setw(OWID) << euler_HOT(1,3)
		        << setw(OWID) << euler_HOT(2,1) << setw(OWID) << euler_HOT(2,2) << setw(OWID) << euler_HOT(2,3)
		        << setw(OWID) << euler_HOT(3,1) << setw(OWID) << euler_HOT(3,2) << setw(OWID) << euler_HOT(3,3)
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 0
		        << setw(OWID) << 0
			<< setw(OWID) << average_dudx_Bagi(1,1) << setw(OWID) << average_dudx_Bagi(1,2) << setw(OWID) << average_dudx_Bagi(1,3)
		        << setw(OWID) << average_dudx_Bagi(2,1) << setw(OWID) << average_dudx_Bagi(2,2) << setw(OWID) << average_dudx_Bagi(2,3)
		        << setw(OWID) << average_dudx_Bagi(3,1) << setw(OWID) << average_dudx_Bagi(3,2) << setw(OWID) << average_dudx_Bagi(3,3)
		        
			<< setw(OWID) << average_dudx_Lagrangian(1,1) << setw(OWID) << average_dudx_Lagrangian(1,2) << setw(OWID) << average_dudx_Lagrangian(1,3)
		        << setw(OWID) << average_dudx_Lagrangian(2,1) << setw(OWID) << average_dudx_Lagrangian(2,2) << setw(OWID) << average_dudx_Lagrangian(2,3)
		        << setw(OWID) << average_dudx_Lagrangian(3,1) << setw(OWID) << average_dudx_Lagrangian(3,2) << setw(OWID) << average_dudx_Lagrangian(3,3)
			<< setw(OWID) << average_dudx_Eulerian(1,1) << setw(OWID) << average_dudx_Eulerian(1,2) << setw(OWID) << average_dudx_Eulerian(1,3)
		        << setw(OWID) << average_dudx_Eulerian(2,1) << setw(OWID) << average_dudx_Eulerian(2,2) << setw(OWID) << average_dudx_Eulerian(2,3)
		        << setw(OWID) << average_dudx_Eulerian(3,1) << setw(OWID) << average_dudx_Eulerian(3,2) << setw(OWID) << average_dudx_Eulerian(3,3)
			<< setw(OWID) << spatial_dvdx(1,1) << setw(OWID) << spatial_dvdx(1,2) << setw(OWID) << spatial_dvdx(1,3)
		        << setw(OWID) << spatial_dvdx(2,1) << setw(OWID) << spatial_dvdx(2,2) << setw(OWID) << spatial_dvdx(2,3)
		        << setw(OWID) << spatial_dvdx(3,1) << setw(OWID) << spatial_dvdx(3,2) << setw(OWID) << spatial_dvdx(3,3) 
			<< setw(OWID) << getCurrCellVolume()
			<< endl;

	    /*
	    g_debuginf << setw(OWID) << g_iteration
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
	    */

	}

	// 8. loop break conditions.

    } while (++g_iteration < total_steps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    g_timeinf << setw(4) << "end" << " " << ctime(&timeStamp) << flush;

    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}

  // this is used to calculate volume of the assembly based on particle input file, July 14, 2013
  void assembly::calculateVolume(const char* iniptclefile){

    readSample(iniptclefile); // create container and particles, velocity and omga are set zero. 
    // first tessellation
    createInputForQhull();
    callQhull();
    readTesse("tess_info");
    readTesse_finite("tess_info");

	REAL curr_totalVolume = 0;	// the current summed volume of all cells

	REAL temp_volume;
	for(std::vector<cell*>::const_iterator iter=cellVec.begin(); iter!=cellVec.end(); iter++){
		temp_volume = getCellVolume(**iter);
		curr_totalVolume += temp_volume;
	}

    std::cout << "the volume of the assembly is: " << std::endl;
    std::cout << curr_totalVolume << std::endl;

}

  // ------------------------------------------- SPH Part --------------------------------------------------
    void assembly::readSPHData(const char *InputFile){
    	
	std::ifstream ifs(InputFile);
	if(!ifs) {
	    std::cout << "stream error!" << std::endl; exit(-1);
	}
		
   	std::string tmp;
	ifs >> nSPHPoint >> SPHxmin >> SPHxmax;
//	ifs >> dem::TIMESTEP >> nsteps >> snapshots >> rampStep;	// these will be defined in DEM
		
	REAL totalMass, SPHInitialDensity;
	ifs >> totalMass >> SPHInitialDensity;
	REAL SPHParticleMass = totalMass/nSPHPoint;

	ifs >> dem::nu >> dem::smoothLength >> dem::soundSpeed;

//	// parameter of the cylinder boundary
//	ifs >> radius >> center_x >> center_y;	// this will be defined by DEM particle

	// free the spaces of SPH points
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    delete (*pt);
	}
	SPHParticleVec.clear();
	// free the spaces of ghost points in each dem particle
	for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
	    for(std::vector<sph::SPHParticle*>::iterator spht=(*pt)->SPHGhostParticleVec.begin(); spht!=(*pt)->SPHGhostParticleVec.end(); spht++){
		delete (*spht);
	   }
	   (*pt)->SPHGhostParticleVec.clear();
	}

	int isGhost = 0;	// is Ghost particle
	REAL radius_a, radius_b, radius_c;
  	dem::vec pt_position;
	REAL pt_x, pt_y, pt_z;	// the position of the DEM particle
	dem::vec local_x, tmp_xyz;	// local position of ghost point in dem particle
	dem::vec tmp_local;
	// read particle information, create and store PeriParticle objects into periParticleVec
	for(int ip=0; ip<nSPHPoint; ip++){
	    REAL tmp_x, tmp_y, tmp_z;
	    int tmp_int;
	    ifs >> tmp_int >> tmp_x >> tmp_y >> tmp_z;
	    isGhost = 0;
	    tmp_xyz = dem::vec(tmp_x, tmp_y, tmp_z);	
	    for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
		radius_a = (*pt)->getA(); radius_b = (*pt)->getB(); radius_c = (*pt)->getC();
		pt_position = (*pt)->getCurrPosition();
		tmp_local = (*pt)->localVec(tmp_xyz-pt_position);
		if( tmp_local.getx()*tmp_local.getx()/(radius_a*radius_a)+tmp_local.gety()*tmp_local.gety()/(radius_b*radius_b)+tmp_local.getz()*tmp_local.getz()/(radius_c*radius_c) <= 1 ){
		    isGhost = 1;	// is Ghost particle
		    local_x = (*pt)->localVec( dem::vec(tmp_x, tmp_y, tmp_z)-pt_position );
		    sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHParticleMass, SPHInitialDensity, tmp_x, tmp_y, tmp_z, local_x, (*pt), 2); // 2 is ghost SPH particle
		    (*pt)->SPHGhostParticleVec.push_back(tmp_pt);
		    break;	// if this sph point is ghost for dem particle 1, then it cannot be ghost for any others,
				// should avoid the initial overlap of the different dem particles
		}
	    }
	    
	    if(isGhost==0){	// not Ghost particle
	    	sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHParticleMass, SPHInitialDensity, tmp_x, tmp_y, tmp_z, 1);	// 1 is free SPH particle
	    	SPHParticleVec.push_back(tmp_pt);
	    }

	}

    } // end readSPHData()

    void assembly::writeSPHParticleTecplot(std::ofstream &ofs, const int iframe) {
	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(10);
	if(iframe == 0) {
	    ofs << "Title = \"SPH Particle Information\"" << std::endl;
	    ofs << "VARIABLES = \"x\", \"y\",\"z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" \"Vz\" \"Pressure\" \"a_x\" \"a_y\" \"a_z\" \"density_dot\" \"density\" " << std::endl;
	}
	ofs << "ZONE T =\" " << iframe << "-th Load Step\" "<< std::endl;
	// Output the coordinates and the array information
	for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleVec.begin(); pt!= SPHParticleVec.end(); pt++) {

//// not print the most right layer of SPH free particles, 2015/05/19
//if((*pt)->getInitPosition().getx()==25){
//continue;
//}

	    ofs << std::setw(20) << (*pt)->getCurrPosition().getx()
		<< std::setw(20) << (*pt)->getCurrPosition().gety()
		<< std::setw(20) << (*pt)->getCurrPosition().getz()
		<< std::setw(20) << (*pt)->getDisplacement().getx()
		<< std::setw(20) << (*pt)->getDisplacement().gety() 
		<< std::setw(20) << (*pt)->getDisplacement().getz() 
		<< std::setw(20) << (*pt)->getVelocity().getx()
		<< std::setw(20) << (*pt)->getVelocity().gety() 
		<< std::setw(20) << (*pt)->getVelocity().getz();

	    (*pt)->calculateParticlePressure();
 	    ofs << std::setw(20) <<(*pt)->getParticlePressure();

	    ofs << std::setw(20) << (*pt)->getVelocityDot().getx()
		<< std::setw(20) << (*pt)->getVelocityDot().gety() 
		<< std::setw(20) << (*pt)->getVelocityDot().getz()
		<< std::setw(20) << (*pt)->getDensityDot() 
		<< std::setw(20) << (*pt)->getParticleDensity() << std::endl;

	}
	// output sph boundary particle
	for(std::vector<sph::SPHParticle*>::iterator pt = SPHBoundaryParticleVec.begin(); pt!= SPHBoundaryParticleVec.end(); pt++) {
	    ofs << std::setw(20) << (*pt)->getCurrPosition().getx()
		<< std::setw(20) << (*pt)->getCurrPosition().gety()
		<< std::setw(20) << (*pt)->getCurrPosition().getz()
		<< std::setw(20) << (*pt)->getDisplacement().getx()
		<< std::setw(20) << (*pt)->getDisplacement().gety() 
		<< std::setw(20) << (*pt)->getDisplacement().getz() 
		<< std::setw(20) << (*pt)->getVelocity().getx()
		<< std::setw(20) << (*pt)->getVelocity().gety() 
		<< std::setw(20) << (*pt)->getVelocity().getz();

	    (*pt)->calculateParticlePressure();
 	    ofs << std::setw(20) <<(*pt)->getParticlePressure();

	    ofs << std::setw(20) << (*pt)->getVelocityDot().getx()
		<< std::setw(20) << (*pt)->getVelocityDot().gety() 
		<< std::setw(20) << (*pt)->getVelocityDot().getz()
		<< std::setw(20) << (*pt)->getDensityDot()
		<< std::setw(20) << (*pt)->getParticleDensity() << std::endl;

	}

	// output the SPH ghost particles
	std::vector<sph::SPHParticle*>::iterator pt;
	for(std::vector<particle*>::iterator demt=ParticleVec.begin(); demt!=ParticleVec.end(); demt++){
	    for(pt = (*demt)->SPHGhostParticleVec.begin(); pt!= (*demt)->SPHGhostParticleVec.end(); pt++) {
	    	ofs << std::setw(20) << (*pt)->getCurrPosition().getx()
		    << std::setw(20) << (*pt)->getCurrPosition().gety()
		    << std::setw(20) << (*pt)->getCurrPosition().getz()
		    << std::setw(20) << (*pt)->getDisplacement().getx()
		    << std::setw(20) << (*pt)->getDisplacement().gety() 
		    << std::setw(20) << (*pt)->getDisplacement().getz() 
		    << std::setw(20) << (*pt)->getVelocity().getx()
		    << std::setw(20) << (*pt)->getVelocity().gety() 
		    << std::setw(20) << (*pt)->getVelocity().getz();

	    	(*pt)->calculateParticlePressure();
 	    	ofs << std::setw(20) <<(*pt)->getParticlePressure();

	    	ofs << std::setw(20) << (*pt)->getVelocityDot().getx()
		    << std::setw(20) << (*pt)->getVelocityDot().gety() 
		    << std::setw(20) << (*pt)->getVelocityDot().getz()
		    << std::setw(20) << (*pt)->getDensityDot()
		    << std::setw(20) << (*pt)->getParticleDensity() << std::endl;

	    }
	}
    } // end writeSPHParticleTecplot()


    void assembly::printSPHVelocity(const char* str){
  	std::ofstream ofs(str);
  	if(!ofs) {
    	    std::cout << "stream error!" << std::endl; exit(-1);
  	}
  	ofs.setf(std::ios::scientific, std::ios::floatfield);
  	ofs.precision(dem::OPREC);
  	ofs << std::setw(dem::OWID) << "current_x"
      	    << std::setw(dem::OWID) << "current_y"
            << std::setw(dem::OWID) << "current_z"
	    << std::setw(dem::OWID) << "initial_X"
      	    << std::setw(dem::OWID) << "initial_Y"
            << std::setw(dem::OWID) << "initial_Z"
	    << std::setw(dem::OWID) << "Vx"
      	    << std::setw(dem::OWID) << "Vy"
            << std::setw(dem::OWID) << "Vz"
      	    << std::endl;
 
  	dem::vec tmp;
	for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleVec.begin(); pt!= SPHParticleVec.end(); pt++) {
	    tmp = (*pt)->getCurrPosition();
	    ofs << std::setw(dem::OWID) << tmp.getx()
		<< std::setw(dem::OWID) << tmp.gety()
		<< std::setw(dem::OWID) << tmp.getz();

	    tmp = (*pt)->getInitPosition();
	    ofs << std::setw(dem::OWID) << tmp.getx()
		<< std::setw(dem::OWID) << tmp.gety()
		<< std::setw(dem::OWID) << tmp.getz();

	    tmp = (*pt)->getVelocity();
	    ofs << std::setw(dem::OWID) << tmp.getx()
		<< std::setw(dem::OWID) << tmp.gety()
		<< std::setw(dem::OWID) << tmp.getz() << std::endl;
	}
	// output boundary particles
	for(std::vector<sph::SPHParticle*>::iterator pt = SPHBoundaryParticleVec.begin(); pt!= SPHBoundaryParticleVec.end(); pt++) {
	    tmp = (*pt)->getCurrPosition();
	    ofs << std::setw(dem::OWID) << tmp.getx()
		<< std::setw(dem::OWID) << tmp.gety()
		<< std::setw(dem::OWID) << tmp.getz();

	    tmp = (*pt)->getInitPosition();
	    ofs << std::setw(dem::OWID) << tmp.getx()
		<< std::setw(dem::OWID) << tmp.gety()
		<< std::setw(dem::OWID) << tmp.getz();

	    tmp = (*pt)->getVelocity();
	    ofs << std::setw(dem::OWID) << tmp.getx()
		<< std::setw(dem::OWID) << tmp.gety()
		<< std::setw(dem::OWID) << tmp.getz() << std::endl;
	}
	// output the ghost sph particle
	std::vector<sph::SPHParticle*>::iterator pt;
	for(std::vector<particle*>::iterator demt=ParticleVec.begin(); demt!=ParticleVec.end(); demt++){
	    for(pt = (*demt)->SPHGhostParticleVec.begin(); pt!= (*demt)->SPHGhostParticleVec.end(); pt++) {
	    	tmp = (*pt)->getCurrPosition();
	    	ofs << std::setw(dem::OWID) << tmp.getx()
		    << std::setw(dem::OWID) << tmp.gety()
		    << std::setw(dem::OWID) << tmp.getz();

	    	tmp = (*pt)->getInitPosition();
	    	ofs << std::setw(dem::OWID) << tmp.getx()
		    << std::setw(dem::OWID) << tmp.gety()
		    << std::setw(dem::OWID) << tmp.getz();

	    	tmp = (*pt)->getVelocity();
	    	ofs << std::setw(dem::OWID) << tmp.getx()
		    << std::setw(dem::OWID) << tmp.gety()
		    << std::setw(dem::OWID) << tmp.getz() << std::endl;
	    }
	}
  
  	ofs.close();
    } // end printSPHVelocity

    void assembly::printSPHFreeParticle(const char* str){
  	std::ofstream ofs(str);
  	if(!ofs) {
    	    std::cout << "stream error!" << std::endl; exit(-1);
  	}
  	ofs.setf(std::ios::scientific, std::ios::floatfield);
  	ofs.precision(dem::OPREC);

	ofs << SPHParticleVec.size() << std::endl;
  	ofs << std::setw(dem::OWID) << "current_x"
      	    << std::setw(dem::OWID) << "current_y"
            << std::setw(dem::OWID) << "current_z"
	    << std::setw(dem::OWID) << "mass"
      	    << std::setw(dem::OWID) << "density"
      	    << std::endl;
 
  	dem::vec tmp;
	for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleVec.begin(); pt!= SPHParticleVec.end(); pt++) {
	    tmp = (*pt)->getCurrPosition();
	    ofs << std::setw(dem::OWID) << tmp.getx()
		<< std::setw(dem::OWID) << tmp.gety()
		<< std::setw(dem::OWID) << tmp.getz()
		<< std::setw(dem::OWID) << (*pt)->getParticleMass()
		<< std::setw(dem::OWID) << (*pt)->getParticleDensity() << std::endl;
	}
  
  	ofs.close();
    } // end printSPHFreeParticle


    void assembly::printSPHVelocityPressure(const char* str){
  	std::ofstream ofs(str);
  	if(!ofs) {
    	    std::cout << "stream error!" << std::endl; exit(-1);
  	}
  	ofs.setf(std::ios::scientific, std::ios::floatfield);
  	ofs.precision(dem::OPREC);
  	ofs << std::setw(dem::OWID) << "current_x"
      	    << std::setw(dem::OWID) << "current_y"
            << std::setw(dem::OWID) << "current_z"
	    << std::setw(dem::OWID) << "Vx"
      	    << std::setw(dem::OWID) << "Vy"
            << std::setw(dem::OWID) << "Vz"
	    << std::setw(dem::OWID) << "pressure"
      	    << std::endl;
 
  	dem::vec tmp;
	for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleVec.begin(); pt!= SPHParticleVec.end(); pt++) {
	    tmp = (*pt)->getCurrPosition();
	    ofs << std::setw(dem::OWID) << tmp.getx()
		<< std::setw(dem::OWID) << tmp.gety()
		<< std::setw(dem::OWID) << tmp.getz();

	    tmp = (*pt)->getVelocity();
	    ofs << std::setw(dem::OWID) << tmp.getx()
		<< std::setw(dem::OWID) << tmp.gety()
		<< std::setw(dem::OWID) << tmp.getz();
	   
	    (*pt)->calculateParticlePressure();
 	    ofs << std::setw(dem::OWID) << (*pt)->getParticlePressure() << std::endl;
	}
	// output boundary sph partilces
	for(std::vector<sph::SPHParticle*>::iterator pt = SPHBoundaryParticleVec.begin(); pt!= SPHBoundaryParticleVec.end(); pt++) {
	    tmp = (*pt)->getCurrPosition();
	    ofs << std::setw(dem::OWID) << tmp.getx()
		<< std::setw(dem::OWID) << tmp.gety()
		<< std::setw(dem::OWID) << tmp.getz();

	    tmp = (*pt)->getVelocity();
	    ofs << std::setw(dem::OWID) << tmp.getx()
		<< std::setw(dem::OWID) << tmp.gety()
		<< std::setw(dem::OWID) << tmp.getz();
	   
	    (*pt)->calculateParticlePressure();
 	    ofs << std::setw(dem::OWID) << (*pt)->getParticlePressure() << std::endl;
	}
	// output sph ghost particle
	std::vector<sph::SPHParticle*>::iterator pt;
	for(std::vector<particle*>::iterator demt=ParticleVec.begin(); demt!=ParticleVec.end(); demt++){
	    for(pt = (*demt)->SPHGhostParticleVec.begin(); pt!= (*demt)->SPHGhostParticleVec.end(); pt++) {
	    	tmp = (*pt)->getCurrPosition();
	    	ofs << std::setw(dem::OWID) << tmp.getx()
		    << std::setw(dem::OWID) << tmp.gety()
		    << std::setw(dem::OWID) << tmp.getz();

	   	tmp = (*pt)->getVelocity();
	    	ofs << std::setw(dem::OWID) << tmp.getx()
		    << std::setw(dem::OWID) << tmp.gety()
		    << std::setw(dem::OWID) << tmp.getz();
	   
	    	(*pt)->calculateParticlePressure();
 	    	ofs << std::setw(dem::OWID) << (*pt)->getParticlePressure() << std::endl;
	    }
	} 
  
  	ofs.close();
    } // end printSPHVelocityPressure


    //// the momentum equilibrium equation and state equation are implemented as in Monaghan's paper (1994), simulate free surface flow using sph
    //// here I would like to use two ways to sum over the smooth kernel for all the SPH particles
    //// (1) implemented below, it is the easy one which is just checking the SPH particles after the home particle
    void assembly::calculateSPHDensityDotVelocityDot(){
	// initialize the densityDot and velocityDot of all the SPH particles 
	dem::vec tmp_vec = dem::vec(0,0, -(dem::G)*(dem::GRVT_SCL) );
	dem::vec zero_vec = dem::vec(0,0,0);
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    (*pt)->setDensityDotVelocityDotZero();
	    (*pt)->addVelocityDot(tmp_vec);
	    (*pt)->calculateParticleViscosity();	// everytime when use getParticleViscolisity(), make sure that pressure has been calculated!!!!
	    (*pt)->calculateParticlePressure();		// everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
	}
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	    (*pt)->setDensityDotVelocityDotZero();
	    (*pt)->addVelocityDot(zero_vec);
	    (*pt)->calculateParticleViscosity();	// everytime when use getParticleViscolisity(), make sure that pressure has been calculated!!!!
	    (*pt)->calculateParticlePressure();		// everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
	}
	std::vector<sph::SPHParticle*>::iterator gt;
	for(std::vector<particle*>::iterator pdt=ParticleVec.begin(); pdt!=ParticleVec.end(); pdt++){
	    for(gt=(*pdt)->SPHGhostParticleVec.begin(); gt!=(*pdt)->SPHGhostParticleVec.end(); gt++){
	    	(*gt)->setDensityDotVelocityDotZero();
	    	(*gt)->addVelocityDot(tmp_vec);
	    	(*gt)->calculateParticleViscosity();	// everytime when use getParticleViscolisity(), make sure that pressure has been calculated!!!!
	    	(*gt)->calculateParticlePressure();	// everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
	    }
	}

	// temporary variables used in the loop
	dem::vec pta_position;
	dem::vec ptb_position;
	dem::vec delta_aWab;
	dem::vec delta_bWba;
	dem::vec vab;
	dem::vec vba;
	dem::vec vdem;
	dem::vec dva_dt;
	dem::vec dvb_dt;
	dem::vec delta_a;
	dem::vec delta_b;
	REAL pa, pb, rhoa, rhob, mua, mub;
	REAL rab;
        REAL dWab_dra;
	REAL dWba_drb;
	REAL Wab, Wba;
	REAL da, dB;
	REAL beta;
	REAL xa, ya, xB, yB, k, sqrt_xaya;	// variables for Morris' method to calculate da/dB
	std::vector<sph::SPHParticle*>::iterator ptb;
	std::vector<particle*>::iterator demt;
	REAL ra, rb, rc;	// the geometry of the dem particle
	dem::vec pta_local, ptb_local;	// local position of pta and ptb in the dem particle
	dem::vec demt_curr;
	REAL Gamma_ab, mu_ab, vr_dot;
	REAL alpha = 0.3;
	REAL epsilon = 0.0;	// parameter for velocity correction
	REAL Wq, Ra, Rb, phi_4, coefficient_a, coefficient_b;

/*
dem::vec sum_delta_aWab;
dem::vec sum_delta_bWba;
REAL sum_dWab_dra;
REAL sum_dWba_drb;
REAL sum_Wab;
REAL sum_Wba;
REAL delta_v;
*/
	for(std::vector<sph::SPHParticle*>::iterator pta=SPHParticleVec.begin(); pta!=SPHParticleVec.end(); pta++){
	    pa = (*pta)->getParticlePressure();
	    if(pa>=0){
		Ra = 0.006;
	    }
	    else{
		Ra = 0.6;
	    }
//	    mua = (*pta)->getParticleViscosity();
	    rhoa = (*pta)->getParticleDensity();
	    pta_position = (*pta)->getCurrPosition();

//std::cout << std::endl << std::endl;
//std::cout << "pta position: "; pta_position.print(); std::cout << std::endl;
/*
sum_delta_aWab = 0;
sum_delta_bWba = 0;
sum_dWab_dra = 0;
sum_dWba_drb = 0;
sum_Wab = 0;
sum_Wba = 0;
*/
	    for(ptb=pta+1; ptb!=SPHParticleVec.end(); ptb++){
		ptb_position = (*ptb)->getCurrPosition();
		rab = dem::vfabs(pta_position-ptb_position);
		if(rab<=3*dem::smoothLength){ // ptb is in the smooth kernel
		    // add the density dot for pta and ptb
		    delta_aWab = gradientKernelFunction(pta_position, ptb_position);	// this is to add SPH pta 
		    delta_bWba = - delta_aWab;
//		    delta_bWba = gradientKernelFunction(ptb_position, pta_position;	// this is to add SPH ptb
		    dWab_dra = partialKernelFunction(pta_position, ptb_position);	// this is to add SPH pta
		    dWba_drb = dWab_dra;
//		    dWba_drb = partialKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
		    Wab = kernelFunction(pta_position, ptb_position);	// this is to add SPH pta
		    Wba = Wab; 
//		    Wba = kernelFunction(ptb_position, pta_position);	// this is to add SPH ptb

/*
delta_v = (*ptb)->getParticleMass()/(*ptb)->getParticleDensity();
sum_delta_aWab = sum_delta_aWab+delta_aWab*delta_v;
sum_delta_bWba = sum_delta_bWba+delta_bWba*delta_v;
sum_dWab_dra = sum_dWab_dra + dWab_dra*delta_v;
sum_dWba_drb = sum_dWba_drb + dWba_drb*delta_v;
sum_Wab = sum_Wab + Wab*delta_v;
sum_Wba = sum_Wba + Wba*delta_v;
*/


		    vab = (*pta)->getVelocity()-(*ptb)->getVelocity();
		    vba = -vab;

		    (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    pb = (*ptb)->getParticlePressure();
//		    mub = (*ptb)->getParticleViscosity();
		    rhob = (*ptb)->getParticleDensity();

		    vr_dot = vab%(pta_position-ptb_position);
		    if(vr_dot<0){
			mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    	Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    }
		    else{
		 	Gamma_ab = 0;
		    }

		    Wq = kernelFunction(rab/dem::smoothLength);
		    phi_4 = pow(Wq/Wqmin,4);
		    if(pb>=0){
			Rb = 0.006;
		    }
		    else{
			Rb = 0.6;
		    }
		    coefficient_a = 1+Ra*phi_4;
		    coefficient_b = 1+Rb*phi_4;

		    dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
//std::cout << "ptb position: "; ptb_position.print();
//std::cout << "ptb mass: " << (*ptb)->getParticleMass() << "  pa: " << pa << "  rhoa: " << rhoa << "  coefficient_a: " << coefficient_a << "  pb: " << pb << "  rhob: " << rhob << "  coefficient_b: " << coefficient_b << "  Gamma: " << Gamma_ab << std::endl;
//std::cout << "delta_aWab: "; delta_aWab.print();
//std::cout << "pta->dva_dt: "; dva_dt.print(); 
		    (*pta)->addVelocityDot(dva_dt);
		    dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
//std::cout << "ptb->dvb_dt: "; dvb_dt.print(); 
		    (*ptb)->addVelocityDot(dvb_dt);

		    delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
		    (*pta)->addVelocityCorrection(delta_a);
		    delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
		    (*pta)->addVelocityCorrection(delta_b);

		} // end if
	    } // end ptb
	    // interact with the Ghost particles
	    for(demt=ParticleVec.begin(); demt!=ParticleVec.end(); demt++){
		demt_curr = (*demt)->getCurrPosition();
		pta_local = (*demt)->localVec(pta_position-demt_curr);	// the local position of sph point pta
		ra = (*demt)->getA(); rb = (*demt)->getB(); rc = (*demt)->getC();
		k = 1.0/(sqrt(pta_local.getx()*pta_local.getx()/(ra*ra)+pta_local.gety()*pta_local.gety()/(rb*rb)+pta_local.getz()*pta_local.getz()/(rc*rc) ) );
		da = dem::vfabs(pta_local-k*pta_local);	// the distance is the same in rotated coordinates
		if(da>3*dem::smoothLength){	// this dem particle will not have ghost points interacting with this sph point	
		    continue;	// move to next dem particle
		}

	    	for(ptb=(*demt)->SPHGhostParticleVec.begin(); ptb!=(*demt)->SPHGhostParticleVec.end(); ptb++){	
		    ptb_position = (*ptb)->getCurrPosition();
		    rab = dem::vfabs(pta_position-ptb_position);
		    if(rab<=3*dem::smoothLength){ // ptb is in the smooth kernel
		    	// add the density dot for pta and ptb
			delta_aWab = gradientKernelFunction(pta_position, ptb_position);	// this is to add SPH pta 
			delta_bWba = - delta_aWab;
//		    	delta_bWba = gradientKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
		    	dWab_dra = partialKernelFunction(pta_position, ptb_position);	// this is to add SPH pta
			dWba_drb = dWab_dra;
//		    	dWba_drb = partialKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
		    	Wab = kernelFunction(pta_position, ptb_position);	// this is to add SPH pta
			Wba = Wab; 
//		    	Wba = kernelFunction(ptb_position, pta_position);	// this is to add SPH ptb

/*
delta_v = (*ptb)->getParticleMass()/(*ptb)->getParticleDensity();
sum_delta_aWab = sum_delta_aWab+delta_aWab*delta_v;
sum_delta_bWba = sum_delta_bWba+delta_bWba*delta_v;
sum_dWab_dra = sum_dWab_dra + dWab_dra*delta_v;
sum_dWba_drb = sum_dWba_drb + dWba_drb*delta_v;
sum_Wab = sum_Wab + Wab*delta_v;
sum_Wba = sum_Wba + Wba*delta_v;
*/
		    	// calculate Vab as the method shown in Morris's paper, 1996

		    	// (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
			ptb_local = (*demt)->localVec(ptb_position-demt_curr);
			k = 1.0/(sqrt(ptb_local.getx()*ptb_local.getx()/(ra*ra)+ptb_local.gety()*ptb_local.gety()/(rb*rb)+ptb_local.getz()*ptb_local.getz()/(rc*rc) ) );
		    	dB = dem::vfabs(ptb_local-k*ptb_local);

//		    	// (2) here the Morris's method is used 
//		    	xa = pta_position.getx(); ya = pta_position.gety();
//		    	xB = ptb_position.getx(); yB = ptb_position.gety();
//		    	if(ya==0) {da = xa-radius; dB = radius-xB;}
//		    	else if(xa==0) {da = ya-radius; dB = radius-yB;}
//		    	else {
//			    sqrt_xaya = sqrt(xa*xa+ya*ya);
//			    k = radius/sqrt_xaya;
//			    da = radius/k - radius;
//			    dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//		    	}


		    	beta = 1+dB/da;
		    	if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

			vdem = (*demt)->getCurrVelocity()+(*demt)->getCurrOmga()*( ptb_position-demt_curr );
		    	vab = beta*((*pta)->getVelocity()-vdem);
		    	vba = -vab;

		    	(*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    	(*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    	pb = (*ptb)->getParticlePressure();
//		    	mub = (*ptb)->getParticleViscosity();
		    	rhob = (*ptb)->getParticleDensity();

		    	vr_dot = vab%(pta_position-ptb_position);
		    	if(vr_dot<0){
			    mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    	    Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	}
		    	else{
		 	    Gamma_ab = 0;
		    	}

		      	Wq = kernelFunction(rab/dem::smoothLength);
		    	phi_4 = pow(Wq/Wqmin,4);
		    	if(pb>=0){
			    Rb = 0.006;
		    	}
		    	else{
			    Rb = 0.6;
		    	}
		    	coefficient_a = 1+Ra*phi_4;
		    	coefficient_b = 1+Rb*phi_4;

		    	dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
		    	(*pta)->addVelocityDot(dva_dt);
		    	dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
			(*demt)->addForce((*ptb)->getParticleMass()*dvb_dt);
			(*demt)->addMoment( (ptb_position-demt_curr) * (*ptb)->getParticleMass()*dvb_dt );
//		    	(*ptb)->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

		        delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
		        (*pta)->addVelocityCorrection(delta_a);
//		        delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//		        (*pta)->addVelocityCorrection(delta_b);

		    } // end if
	    	} // end ptb
	    } // end demt
	    // interact with boundary particles
	    da = pta_position.getz();	// assume with the bottom boundary
	    if(pta_position.getx() < da){	// with left boundary
		da = pta_position.getx();
	    } 
	    else if(length_x-pta_position.getx() < da){	// with right boundary
		da = length_x-pta_position.getx();
	    }
	    else if(height_z-pta_position.getz() < da){	// with top boundary
		da = height_z-pta_position.getz();
	    }


	    for(ptb=SPHBoundaryParticleVec.begin(); ptb!=SPHBoundaryParticleVec.end(); ptb++){	
		ptb_position = (*ptb)->getCurrPosition();
		rab = dem::vfabs(pta_position-ptb_position);
		if(rab<=3*dem::smoothLength){ // ptb is in the smooth kernel
		    // add the density dot for pta and ptb
		    delta_aWab = gradientKernelFunction(pta_position, ptb_position);	// this is to add SPH pta 
		    delta_bWba = - delta_aWab;
//		    delta_bWba = gradientKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
		    dWab_dra = partialKernelFunction(pta_position, ptb_position);	// this is to add SPH pta
		    dWba_drb = dWab_dra;
//		    dWba_drb = partialKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
		    Wab = kernelFunction(pta_position, ptb_position);	// this is to add SPH pta
		    Wba = Wab; 
//		    Wba = kernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
/*
delta_v = (*ptb)->getParticleMass()/(*ptb)->getParticleDensity();
sum_delta_aWab = sum_delta_aWab+delta_aWab*delta_v;
sum_delta_bWba = sum_delta_bWba+delta_bWba*delta_v;
sum_dWab_dra = sum_dWab_dra + dWab_dra*delta_v;
sum_dWba_drb = sum_dWba_drb + dWba_drb*delta_v;
sum_Wab = sum_Wab + Wab*delta_v;
sum_Wba = sum_Wba + Wba*delta_v;
*/
		    // calculate Vab as the method shown in Morris's paper, 1996

	    	    dB = 0-ptb_position.getz();	// assume with the bottom boundary
	    	    if(ptb_position.getx() < dB){	// with left boundary
		    	dB = 0-ptb_position.getx();
	    	    } 
	    	    else if(ptb_position.getx()-length_x < dB){	// with right boundary
		    	dB = ptb_position.getx()-length_x;
	    	    }
	    	    else if(ptb_position.getz()-height_z < dB){	// with top boundary
		    	dB = ptb_position.getz()-height_z;
	    	    }

		    beta = 1+dB/da;
		    if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

		    vab = beta*(*pta)->getVelocity();
		    vba = -vab;

		    (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    pb = (*ptb)->getParticlePressure();
//		    mub = (*ptb)->getParticleViscosity();
		    rhob = (*ptb)->getParticleDensity();

		    vr_dot = vab%(pta_position-ptb_position);
		    if(vr_dot<0){
			mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    	Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    }
		    else{
		 	Gamma_ab = 0;
		    }

		    Wq = kernelFunction(rab/dem::smoothLength);
		    phi_4 = pow(Wq/Wqmin,4);
		    if(pb>=0){
			Rb = 0.006;
		    }
		    else{
			Rb = 0.6;
		    }
		    coefficient_a = 1+Ra*phi_4;
		    coefficient_b = 1+Rb*phi_4;

		    dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
		    (*pta)->addVelocityDot(dva_dt);
//		    dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
//		    (*ptb)->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

		    delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
		    (*pta)->addVelocityCorrection(delta_a);
//		    delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//		    (*pta)->addVelocityCorrection(delta_b);

		    // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
		    if(rab<=space_interval){ // ptb is in the smooth kernel
		    	dva_dt = D*(pow(space_interval/rab, p1)-pow(space_interval/rab, p2))*(pta_position-ptb_position)/(rab*rab);
		    	(*pta)->addVelocityDot(dva_dt);
		    } // end if
		} // end if
	    } // end ptb


//	    // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
//	    for(ptb=SPHBoundaryParticleVec.begin(); ptb!=SPHBoundaryParticleVec.end(); ptb++){	
//		ptb_position = (*ptb)->getCurrPosition();
//		rab = dem::vfabs(pta_position-ptb_position);
//		if(rab<=space_interval){ // ptb is in the smooth kernel
//		    dva_dt = D*(pow(space_interval/rab, p1)-pow(space_interval/rab, p2))*(pta_position-ptb_position)/(rab*rab);
//		    (*pta)->addVelocityDot(dva_dt);
//		} // end if
//	    } // end ptb

/*
std::cout << std::endl;
std::cout << "particle position: "; (*pta)->getCurrPosition().print(); std::cout << std::endl;
std::cout << "sum_delta_aWab: "; sum_delta_aWab.print(); std::cout << std::endl;
std::cout << "sum_delta_bWba: "; sum_delta_bWba.print(); std::cout << std::endl;
std::cout << "sum_dWab_dra: " << sum_dWab_dra << std::endl;
std::cout << "sum_dWba_drb: " << sum_dWba_drb << std::endl;
std::cout << "sum_Wab: " << sum_Wab << std::endl;
std::cout << "sum_Wba: " << sum_Wba << std::endl;
*/

	} // end pta

    } // end calculateSPHDensityDotVelocityDot()

/*  // misunderstanding of boundary conditions in Monaghan's paper(1994)
    void assembly::applySPHBoundaryForce(){
	REAL p1 = 4; 
	REAL p2 = 2;
	dem::vec pt_position;
	REAL ptx, ptz;
	REAL tmp_a;
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    pt_position = (*pt)->getCurrPosition();
	    ptx = pt_position.getx();
	    ptz = pt_position.getz();
	    if(ptx<=space_interval){	// left boundary
		tmp_a = D*(pow(space_interval/ptx, p1)-pow(space_interval/ptx, p2))/(ptx*ptx);
		(*pt)->addVelocityDot(dem::vec(tmp_a,0,0));
	    }
	    if((length_x-ptx)<=space_interval){	// right boundary
		tmp_a = D*(pow(space_interval/(length_x-ptx), p1)-pow(space_interval/(length_x-ptx), p2))/((length_x-ptx)*(length_x-ptx));
		(*pt)->addVelocityDot(dem::vec(-tmp_a,0,0));
	    }
	    if(ptz<=space_interval){	// bottom boundary
		tmp_a = D*(pow(space_interval/ptz, p1)-pow(space_interval/ptz, p2))/(ptz*ptz);
		(*pt)->addVelocityDot(dem::vec(0,0,tmp_a));
	    }
	    if((height_z-ptz)<=space_interval){	// top boundary
		tmp_a = D*(pow(space_interval/(height_z-ptz), p1)-pow(space_interval/(height_z-ptz), p2))/((height_z-ptz)*(height_z-ptz));
		(*pt)->addVelocityDot(dem::vec(0,0,-tmp_a));
	    }
	}
    } // end applySPHBoundaryForce
*/

    //// here I would like to use two ways to sum over the smooth kernel for all the SPH particles
    //// (1) implemented below, it is the easy one which is just checking the SPH particles after the home particle
    void assembly::calculateSPHDensityDotVelocityDotPeriodic(REAL L){
	// initialize the densityDot and velocityDot of all the SPH particles 
	dem::vec tmp_vec = dem::vec(0,0, -(dem::G)*(dem::GRVT_SCL) );
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    (*pt)->setDensityDotVelocityDotZero();
	    (*pt)->addVelocityDot(tmp_vec);
	    (*pt)->calculateParticleViscosity();	// everytime when use getParticleViscolisity(), make sure that pressure has been calculated!!!!
	    (*pt)->calculateParticlePressure();		// everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
	}
	std::vector<sph::SPHParticle*>::iterator gt;
	for(std::vector<particle*>::iterator pdt=ParticleVec.begin(); pdt!=ParticleVec.end(); pdt++){
	    for(gt=(*pdt)->SPHGhostParticleVec.begin(); gt!=(*pdt)->SPHGhostParticleVec.end(); gt++){
	    	(*gt)->setDensityDotVelocityDotZero();
	    	(*gt)->addVelocityDot(tmp_vec);
	    	(*gt)->calculateParticleViscosity();	// everytime when use getParticleViscolisity(), make sure that pressure has been calculated!!!!
	    	(*gt)->calculateParticlePressure();	// everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
	    }
	}

	// temporary variables used in the loop
	dem::vec pta_position;
	dem::vec ptb_position;
	dem::vec delta_aWab;
	dem::vec delta_bWba;
	dem::vec vab;
	dem::vec vba;
	dem::vec vdem;
	dem::vec dva_dt;
	dem::vec dvb_dt;
	REAL pa, pb, rhoa, rhob, mua, mub;
	REAL rab;
        REAL dWab_dra;
	REAL dWba_drb;
	REAL dB, da;
	REAL beta;
	REAL xa, ya, xB, yB, k, sqrt_xaya;	// variables for Morris' method to calculate da/dB
	std::vector<sph::SPHParticle*>::iterator ptb;
	std::vector<particle*>::iterator demt;
	REAL ra, rb, rc;	// the geometry of the dem particle
	dem::vec pta_local, ptb_local;	// local position of pta and ptb in the dem particle
	dem::vec demt_curr, demt_curr_origin;
	for(std::vector<sph::SPHParticle*>::iterator pta=SPHParticleVec.begin(); pta!=SPHParticleVec.end(); pta++){
	    pa = (*pta)->getParticlePressure();
	    mua = (*pta)->getParticleViscosity();
	    rhoa = (*pta)->getParticleDensity();
	    pta_position = (*pta)->getCurrPosition();
	    for(ptb=pta+1; ptb!=SPHParticleVec.end(); ptb++){
		ptb_position = (*ptb)->getCurrPosition();
		// check if the ptb should be projected to the other side
		// apply periodic boundary condition, http://en.wikipedia.org/wiki/Periodic_boundary_conditions
		// Second, every distance between objects should obey the minimum image criterion.
		if( fabs(ptb_position.getx()-pta_position.getx()) > 0.5*L ){	// should project ptb.getx
		    if(ptb_position.getx() > 0) {
		    	ptb_position.setx( -L+ptb_position.getx() );	// only valid for box centered at the origin
		    }
		    else{
			ptb_position.setx( L+ptb_position.getx() );
		    }
		}
		if( fabs(ptb_position.gety()-pta_position.gety()) > 0.5*L ){	// should project ptb.gety
		    if(ptb_position.gety() > 0) {
		    	ptb_position.sety( -L+ptb_position.gety() );	// only valid for box centered at the origin
		    }
		    else{
			ptb_position.sety( L+ptb_position.gety() );
		    }
		}

		// this is to avoid that a particle in the SPH boundary being projected to the same particle in the other boundary
		if(pta_position==ptb_position){
		    ptb_position = (*ptb)->getCurrPosition();	// then not project it
		}

		rab = dem::vfabs(pta_position-ptb_position);
		if(rab<=3*dem::smoothLength){ // ptb is in the smooth kernel
		    // add the density dot for pta and ptb
		    delta_aWab = gradientKernelFunction(pta_position, ptb_position);	// this is to add SPH pta 
		    delta_bWba = - delta_aWab;
//		    delta_bWba = gradientKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
		    dWab_dra = partialKernelFunction(pta_position, ptb_position);	// this is to add SPH pta
		    dWba_drb = dWab_dra;
//		    dWba_drb = partialKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb

		    vab = (*pta)->getVelocity()-(*ptb)->getVelocity();
		    vba = -vab;

		    (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    pb = (*ptb)->getParticlePressure();
		    mub = (*ptb)->getParticleViscosity();
		    rhob = (*ptb)->getParticleDensity();

		    dva_dt = (*ptb)->getParticleMass()*( -(pa/(rhoa*rhoa)+pb/(rhob*rhob))*delta_aWab + (mua+mub)*vab/(rhoa*rhob)/rab*dWab_dra );
		    (*pta)->addVelocityDot(dva_dt);
		    dvb_dt = (*pta)->getParticleMass()*( -(pb/(rhob*rhob)+pa/(rhoa*rhoa))*delta_bWba + (mub+mua)*vba/(rhob*rhoa)/rab*dWba_drb );
		    (*ptb)->addVelocityDot(dvb_dt);

//if(dem::vfabs(dva_dt)>100){
//std::cout << "dva_dt in free particles: "; dva_dt.print();
//std::cout << "dvb_dt in free particles: "; dvb_dt.print();
//std::cout << "rab: " << rab << "  rhoa: " << rhoa << "  rhob: " << rhob << std::endl;
//std::cout << "pa: " << pa << "   pb: " << pb << std::endl; 
//std::cout << "delta_aWab: "; delta_aWab.print();
//std::cout << "vab: "; vab.print(); std::cout << std::endl;
//}

//std::cout << "dva_dt in free particles: "; dva_dt.print();
//std::cout << "dvb_dt in free particles: "; dvb_dt.print();

		} // end if
	    } // end ptb
	    // interact with the Ghost particles
	    for(demt=ParticleVec.begin(); demt!=ParticleVec.end(); demt++){
		demt_curr = (*demt)->getCurrPosition();
		demt_curr_origin = demt_curr;
		// check if the dem particle should be projected to the other side
		// apply periodic boundary condition, http://en.wikipedia.org/wiki/Periodic_boundary_conditions
		// Second, every distance between objects should obey the minimum image criterion.
		if( fabs(demt_curr.getx()-pta_position.getx()) > 0.5*L ){	// should project ptb.getx
		    if(demt_curr.getx() > 0) {
		    	demt_curr.setx( -L+demt_curr.getx() );	// only valid for box centered at the origin
		     }
		    else{
			demt_curr.setx( L+demt_curr.getx() );
		    }
		}
		if( fabs(demt_curr.gety()-pta_position.gety()) > 0.5*L ){	// should project ptb.gety
		    if(demt_curr.gety() > 0) {
		    	demt_curr.sety( -L+demt_curr.gety() );	// only valid for box centered at the origin
		    }
		    else{
			demt_curr.sety( L+demt_curr.gety() );
		    }
		}
		pta_local = (*demt)->localVec(pta_position-demt_curr);	// the local position of sph point pta
		ra = (*demt)->getA(); rb = (*demt)->getB(); rc = (*demt)->getC();
		k = 1.0/(sqrt(pta_local.getx()*pta_local.getx()/(ra*ra)+pta_local.gety()*pta_local.gety()/(rb*rb)+pta_local.getz()*pta_local.getz()/(rc*rc) ) );
		da = dem::vfabs(pta_local-k*pta_local);	// the distance is the same in rotated coordinates
		if(da>3*dem::smoothLength){	// this dem particle will not have ghost points interacting with this sph point	
		    continue;	// move to next dem particle
		}

	    	for(ptb=(*demt)->SPHGhostParticleVec.begin(); ptb!=(*demt)->SPHGhostParticleVec.end(); ptb++){	
		    ptb_position = (*ptb)->getCurrPosition();
		    // check if the ptb should be projected to the other side
		    // apply periodic boundary condition, http://en.wikipedia.org/wiki/Periodic_boundary_conditions
		    // Second, every distance between objects should obey the minimum image criterion.
		    if( fabs(ptb_position.getx()-pta_position.getx()) > 0.5*L ){	// should project ptb.getx
		    	if(ptb_position.getx() > 0) {
		    	    ptb_position.setx( -L+ptb_position.getx() );	// only valid for box centered at the origin
		     	}
		    	else{
			    ptb_position.setx( L+ptb_position.getx() );
		    	}
		    }
		    if( fabs(ptb_position.gety()-pta_position.gety()) > 0.5*L ){	// should project ptb.gety
		    	if(ptb_position.gety() > 0) {
		    	    ptb_position.sety( -L+ptb_position.gety() );	// only valid for box centered at the origin
		    	}
		    	else{
			    ptb_position.sety( L+ptb_position.gety() );
		    	}
		    }

		    // this is to avoid that a particle in the SPH boundary being projected to the same particle in the other boundary
		    if(pta_position==ptb_position){
		    	ptb_position = (*ptb)->getCurrPosition();	// then not project it
		    }

		    rab = dem::vfabs(pta_position-ptb_position);
		    if(rab<=3*dem::smoothLength){ // ptb is in the smooth kernel
		    	// add the density dot for pta and ptb
			delta_aWab = gradientKernelFunction(pta_position, ptb_position);	// this is to add SPH pta 
			delta_bWba = - delta_aWab;
//		    	delta_bWba = gradientKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
		    	dWab_dra = partialKernelFunction(pta_position, ptb_position);	// this is to add SPH pta
			dWba_drb = dWab_dra;
//		    	dWba_drb = partialKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb

		    	// calculate Vab as the method shown in Morris's paper, 1996
		    	// (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
			ptb_local = (*demt)->localVec((*ptb)->getCurrPosition()-demt_curr_origin);
			k = 1.0/(sqrt(ptb_local.getx()*ptb_local.getx()/(ra*ra)+ptb_local.gety()*ptb_local.gety()/(rb*rb)+ptb_local.getz()*ptb_local.getz()/(rc*rc) ) );
		    	dB = dem::vfabs(ptb_local-k*ptb_local);

//		    	// (2) here the Morris's method is used 
//		    	xa = pta_position.getx(); ya = pta_position.gety();
//		    	xB = ptb_position.getx(); yB = ptb_position.gety();
//		    	if(ya==0) {da = xa-radius; dB = radius-xB;}
//		    	else if(xa==0) {da = ya-radius; dB = radius-yB;}
//		    	else {
//			    sqrt_xaya = sqrt(xa*xa+ya*ya);
//			    k = radius/sqrt_xaya;
//			    da = radius/k - radius;
//			    dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//		    	}

		    	beta = 1+dB/da;
		    	if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

			vdem = (*demt)->getCurrVelocity()+(*demt)->getCurrOmga()*( (*ptb)->getCurrPosition()-demt_curr_origin );
		    	vab = beta*((*pta)->getVelocity()-vdem);
		    	vba = -vab;

		    	(*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    	(*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    	pb = (*ptb)->getParticlePressure();
		    	mub = (*ptb)->getParticleViscosity();
		        rhob = (*ptb)->getParticleDensity();

		        dva_dt = (*ptb)->getParticleMass()*( -(pa/(rhoa*rhoa)+pb/(rhob*rhob))*delta_aWab + (mua+mub)*vab/(rhoa*rhob)/rab*dWab_dra );
		        (*pta)->addVelocityDot(dva_dt);
		        dvb_dt = (*pta)->getParticleMass()*( -(pb/(rhob*rhob)+pa/(rhoa*rhoa))*delta_bWba + (mub+mua)*vba/(rhob*rhoa)/rab*dWba_drb );
//if(dem::vfabs(dva_dt)>100){
//std::cout << "dva_dt in ghost particles: "; dva_dt.print(); 
//std::cout << "dvb_dt in ghost particles: "; dvb_dt.print();
//std::cout << "rab: " << rab << "  rhoa: " << rhoa << "  rhob: " << rhob << std::endl;
//std::cout << "pa: " << pa << "   pb: " << pb << std::endl; 
//std::cout << "delta_aWab: "; delta_aWab.print();
//std::cout << "vab: "; vab.print(); std::cout << std::endl;
//}
			(*demt)->addForce((*ptb)->getParticleMass()*dvb_dt);
			(*demt)->addMoment( ((*ptb)->getCurrPosition()-demt_curr_origin) * (*ptb)->getParticleMass()*dvb_dt );
//		    (*ptb)->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

		    } // end if
	    	} // end ptb
	    } // end demt
	} // end pta

    } // end calculateSPHDensityDotVelocityDotPeriodic()


/*
    //// here I would like to use two ways to sum over the smooth kernel for all the SPH particles
    //// (2) implemented below, partition the computational domain into several parts in 1D. The method
    //// is described in Morris's Ph.D. thesis
    void SPHDomain::calculateDensityDotVelocityDot(){


	// suppose we partition the domain space into 4 parts in x dimension
	std::vector<sph::SPHParticle*> SPHParticleVecPart1;	// the very left one
	std::vector<sph::SPHParticle*> SPHParticleVecPart2;	// the second one
	std::vector<sph::SPHParticle*> SPHParticleVecPart3;	// the third one
	std::vector<sph::SPHParticle*> SPHParticleVecPart4;	// the very right one
	REAL xsecond, xthird, xfourth;
	xthird = 0.5*(SPHXmin+SPHXmax);	// the middle coordinate
	xsecond = 0.5*(SPHXmin+xthird);
	xfourth = 0.5*(xthird+SPHXmax);	

	REAL pt_x;
	dem::vec tmp_vec = dem::vec(0,0, -(dem::G)*(dem::GRVT_SCL) );
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    // initialize the densityDot and velocityDot of all the SPH particles 
	    (*pt)->setDensityDotVelocityDotZero();
	    (*pt)->addVelocityDot(tmp_vec);

	    pt_x = (*pt)->getCurrPosition().getx();
	    if(pt_x<=xsecond){
		SPHParticleVecPart1.push_back(*pt);
	    }
	    else if(pt_x<=xthird){
		SPHParticleVecPart2.push_back(*pt);
	    }
	    else if(pt_x<=xfourth){
		SPHParticleVecPart3.push_back(*pt);
	    }
	    else{
		SPHParticleVecPart4.push_back(*pt);
	    }
	}


	// temporary variables used in the loop
	dem::vec pta_position;
	dem::vec ptb_position;
	dem::vec delta_aWab;
	dem::vec delta_bWba;
	dem::vec vab;
	dem::vec vba;
	dem::vec dva_dt;
	dem::vec dvb_dt;
	REAL pa, pb, rhoa, rhob, mua, mub;
	REAL rab;
        REAL dWab_dra;
	REAL dWba_drb;

	std::vector<sph::SPHParticle*>::iterator pta;
	std::vector<sph::SPHParticle*>::iterator ptb;
	for(pta=SPHParticleVecPart1.begin(); pta!=SPHParticleVecPart1.end(); pta++){
	    pta_position = (*pta)->getCurrPosition();
	    for(ptb=pta+1; ptb!=SPHParticleVecPart1.end(); ptb++){
		ptb_position = (*ptb)->getCurrPosition();
		rab = dem::vfabs(pta_position-ptb_position);
		if(rab<=3*dem::smoothLength){ // ptb is in the smooth kernel
		    // add the density dot for pta and ptb
		    delta_aWab = periDynamics::gradientKernelFunction(pta_position, ptb_position, dem::smoothLength);	// this is to add SPH pta 
		    delta_bWba = periDynamics::gradientKernelFunction(ptb_position, pta_position, dem::smoothLength);	// this is to add SPH ptb
		    dWab_dra = periDynamics::partialKernelFunction(pta_position, ptb_position, dem::smoothLength);	// this is to add SPH pta
		    dWba_drb = periDynamics::partialKernelFunction(ptb_position, pta_position, dem::smoothLength);	// this is to add SPH ptb

		    vab = (*pta)->getVelocity()-(*ptb)->getVelocity();
		    vba = -vab;

		    (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    pa = (*pta)->getParticlePressure();
		    pb = (*ptb)->getParticlePressure();
		    mua = (*pta)->getParticleViscosity();
		    mub = (*ptb)->getParticleViscosity();
		    rhoa = (*pta)->getParticleDensity();
		    rhob = (*ptb)->getParticleDensity();

		    dva_dt = (*ptb)->getParticleMass()*( -(pa/(rhoa*rhoa)+pb/(rhob*rhob))*delta_aWab + (mua+mub)*vab/(rhoa*rhob)/rab*dWab_dra );
		    (*pta)->addVelocityDot(dva_dt);
		    dvb_dt = (*pta)->getParticleMass()*( -(pb/(rhob*rhob)+pa/(rhoa*rhoa))*delta_bWba + (mub+mua)*vba/(rhob*rhoa)/rab*dWba_drb );
		    (*ptb)->addVelocityDot(dvb_dt);
		}
	    }
	    // sum over the rightside partition
	    for(ptb=SPHParticleVecPart2.begin(); ptb!=SPHParticleVecPart2.end(); ptb++){
		ptb_position = (*ptb)->getCurrPosition();
		rab = dem::vfabs(pta_position-ptb_position);
		if(rab<=3*dem::smoothLength){ // ptb is in the smooth kernel
		    // add the density dot for pta and ptb
		    delta_aWab = periDynamics::gradientKernelFunction(pta_position, ptb_position, dem::smoothLength);	// this is to add SPH pta 
		    delta_bWba = periDynamics::gradientKernelFunction(ptb_position, pta_position, dem::smoothLength);	// this is to add SPH ptb
		    dWab_dra = periDynamics::partialKernelFunction(pta_position, ptb_position, dem::smoothLength);	// this is to add SPH pta
		    dWba_drb = periDynamics::partialKernelFunction(ptb_position, pta_position, dem::smoothLength);	// this is to add SPH ptb

		    vab = (*pta)->getVelocity()-(*ptb)->getVelocity();
		    vba = -vab;

		    (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    pa = (*pta)->getParticlePressure();
		    pb = (*ptb)->getParticlePressure();
		    mua = (*pta)->getParticleViscosity();
		    mub = (*ptb)->getParticleViscosity();
		    rhoa = (*pta)->getParticleDensity();
		    rhob = (*ptb)->getParticleDensity();

		    dva_dt = (*ptb)->getParticleMass()*( -(pa/(rhoa*rhoa)+pb/(rhob*rhob))*delta_aWab + (mua+mub)*vab/(rhoa*rhob)/rab*dWab_dra );
		    (*pta)->addVelocityDot(dva_dt);
		    dvb_dt = (*pta)->getParticleMass()*( -(pb/(rhob*rhob)+pa/(rhoa*rhoa))*delta_bWba + (mub+mua)*vba/(rhob*rhoa)/rab*dWba_drb );
		    (*ptb)->addVelocityDot(dvb_dt);
		}
	    }
	}
	for(pta=SPHParticleVecPart2.begin(); pta!=SPHParticleVecPart2.end(); pta++){
	    pta_position = (*pta)->getCurrPosition();
	    for(ptb=pta+1; ptb!=SPHParticleVecPart2.end(); ptb++){
		ptb_position = (*ptb)->getCurrPosition();
		rab = dem::vfabs(pta_position-ptb_position);
		if(rab<=3*dem::smoothLength){ // ptb is in the smooth kernel
		    // add the density dot for pta and ptb
		    delta_aWab = periDynamics::gradientKernelFunction(pta_position, ptb_position, dem::smoothLength);	// this is to add SPH pta 
		    delta_bWba = periDynamics::gradientKernelFunction(ptb_position, pta_position, dem::smoothLength);	// this is to add SPH ptb
		    dWab_dra = periDynamics::partialKernelFunction(pta_position, ptb_position, dem::smoothLength);	// this is to add SPH pta
		    dWba_drb = periDynamics::partialKernelFunction(ptb_position, pta_position, dem::smoothLength);	// this is to add SPH ptb

		    vab = (*pta)->getVelocity()-(*ptb)->getVelocity();
		    vba = -vab;

		    (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    pa = (*pta)->getParticlePressure();
		    pb = (*ptb)->getParticlePressure();
		    mua = (*pta)->getParticleViscosity();
		    mub = (*ptb)->getParticleViscosity();
		    rhoa = (*pta)->getParticleDensity();
		    rhob = (*ptb)->getParticleDensity();

		    dva_dt = (*ptb)->getParticleMass()*( -(pa/(rhoa*rhoa)+pb/(rhob*rhob))*delta_aWab + (mua+mub)*vab/(rhoa*rhob)/rab*dWab_dra );
		    (*pta)->addVelocityDot(dva_dt);
		    dvb_dt = (*pta)->getParticleMass()*( -(pb/(rhob*rhob)+pa/(rhoa*rhoa))*delta_bWba + (mub+mua)*vba/(rhob*rhoa)/rab*dWba_drb );
		    (*ptb)->addVelocityDot(dvb_dt);
		}
	    }
	    // sum over the rightside partition
	    for(ptb=SPHParticleVecPart3.begin(); ptb!=SPHParticleVecPart3.end(); ptb++){
		ptb_position = (*ptb)->getCurrPosition();
		rab = dem::vfabs(pta_position-ptb_position);
		if(rab<=3*dem::smoothLength){ // ptb is in the smooth kernel
		    // add the density dot for pta and ptb
		    delta_aWab = periDynamics::gradientKernelFunction(pta_position, ptb_position, dem::smoothLength);	// this is to add SPH pta 
		    delta_bWba = periDynamics::gradientKernelFunction(ptb_position, pta_position, dem::smoothLength);	// this is to add SPH ptb
		    dWab_dra = periDynamics::partialKernelFunction(pta_position, ptb_position, dem::smoothLength);	// this is to add SPH pta
		    dWba_drb = periDynamics::partialKernelFunction(ptb_position, pta_position, dem::smoothLength);	// this is to add SPH ptb

		    vab = (*pta)->getVelocity()-(*ptb)->getVelocity();
		    vba = -vab;

		    (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    pa = (*pta)->getParticlePressure();
		    pb = (*ptb)->getParticlePressure();
		    mua = (*pta)->getParticleViscosity();
		    mub = (*ptb)->getParticleViscosity();
		    rhoa = (*pta)->getParticleDensity();
		    rhob = (*ptb)->getParticleDensity();

		    dva_dt = (*ptb)->getParticleMass()*( -(pa/(rhoa*rhoa)+pb/(rhob*rhob))*delta_aWab + (mua+mub)*vab/(rhoa*rhob)/rab*dWab_dra );
		    (*pta)->addVelocityDot(dva_dt);
		    dvb_dt = (*pta)->getParticleMass()*( -(pb/(rhob*rhob)+pa/(rhoa*rhoa))*delta_bWba + (mub+mua)*vba/(rhob*rhoa)/rab*dWba_drb );
		    (*ptb)->addVelocityDot(dvb_dt);
		}
	    }
	}
	for(pta=SPHParticleVecPart3.begin(); pta!=SPHParticleVecPart3.end(); pta++){
	    pta_position = (*pta)->getCurrPosition();
	    for(ptb=pta+1; ptb!=SPHParticleVecPart3.end(); ptb++){
		ptb_position = (*ptb)->getCurrPosition();
		rab = dem::vfabs(pta_position-ptb_position);
		if(rab<=3*dem::smoothLength){ // ptb is in the smooth kernel
		    // add the density dot for pta and ptb
		    delta_aWab = periDynamics::gradientKernelFunction(pta_position, ptb_position, dem::smoothLength);	// this is to add SPH pta 
		    delta_bWba = periDynamics::gradientKernelFunction(ptb_position, pta_position, dem::smoothLength);	// this is to add SPH ptb
		    dWab_dra = periDynamics::partialKernelFunction(pta_position, ptb_position, dem::smoothLength);	// this is to add SPH pta
		    dWba_drb = periDynamics::partialKernelFunction(ptb_position, pta_position, dem::smoothLength);	// this is to add SPH ptb

		    vab = (*pta)->getVelocity()-(*ptb)->getVelocity();
		    vba = -vab;

		    (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    pa = (*pta)->getParticlePressure();
		    pb = (*ptb)->getParticlePressure();
		    mua = (*pta)->getParticleViscosity();
		    mub = (*ptb)->getParticleViscosity();
		    rhoa = (*pta)->getParticleDensity();
		    rhob = (*ptb)->getParticleDensity();

		    dva_dt = (*ptb)->getParticleMass()*( -(pa/(rhoa*rhoa)+pb/(rhob*rhob))*delta_aWab + (mua+mub)*vab/(rhoa*rhob)/rab*dWab_dra );
		    (*pta)->addVelocityDot(dva_dt);
		    dvb_dt = (*pta)->getParticleMass()*( -(pb/(rhob*rhob)+pa/(rhoa*rhoa))*delta_bWba + (mub+mua)*vba/(rhob*rhoa)/rab*dWba_drb );
		    (*ptb)->addVelocityDot(dvb_dt);
		}
	    }
	    // sum over the rightside partition
	    for(ptb=SPHParticleVecPart4.begin(); ptb!=SPHParticleVecPart4.end(); ptb++){
		ptb_position = (*ptb)->getCurrPosition();
		rab = dem::vfabs(pta_position-ptb_position);
		if(rab<=3*dem::smoothLength){ // ptb is in the smooth kernel
		    // add the density dot for pta and ptb
		    delta_aWab = periDynamics::gradientKernelFunction(pta_position, ptb_position, dem::smoothLength);	// this is to add SPH pta 
		    delta_bWba = periDynamics::gradientKernelFunction(ptb_position, pta_position, dem::smoothLength);	// this is to add SPH ptb
		    dWab_dra = periDynamics::partialKernelFunction(pta_position, ptb_position, dem::smoothLength);	// this is to add SPH pta
		    dWba_drb = periDynamics::partialKernelFunction(ptb_position, pta_position, dem::smoothLength);	// this is to add SPH ptb

		    vab = (*pta)->getVelocity()-(*ptb)->getVelocity();
		    vba = -vab;

		    (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    pa = (*pta)->getParticlePressure();
		    pb = (*ptb)->getParticlePressure();
		    mua = (*pta)->getParticleViscosity();
		    mub = (*ptb)->getParticleViscosity();
		    rhoa = (*pta)->getParticleDensity();
		    rhob = (*ptb)->getParticleDensity();

		    dva_dt = (*ptb)->getParticleMass()*( -(pa/(rhoa*rhoa)+pb/(rhob*rhob))*delta_aWab + (mua+mub)*vab/(rhoa*rhob)/rab*dWab_dra );
		    (*pta)->addVelocityDot(dva_dt);
		    dvb_dt = (*pta)->getParticleMass()*( -(pb/(rhob*rhob)+pa/(rhoa*rhoa))*delta_bWba + (mub+mua)*vba/(rhob*rhoa)/rab*dWba_drb );
		    (*ptb)->addVelocityDot(dvb_dt);
		}
	    }
	}
	for(pta=SPHParticleVecPart4.begin(); pta!=SPHParticleVecPart4.end(); pta++){
	    pta_position = (*pta)->getCurrPosition();
	    for(ptb=pta+1; ptb!=SPHParticleVecPart4.end(); ptb++){
		ptb_position = (*ptb)->getCurrPosition();
		rab = dem::vfabs(pta_position-ptb_position);
		if(rab<=3*dem::smoothLength){ // ptb is in the smooth kernel
		    // add the density dot for pta and ptb
		    delta_aWab = periDynamics::gradientKernelFunction(pta_position, ptb_position, dem::smoothLength);	// this is to add SPH pta 
		    delta_bWba = periDynamics::gradientKernelFunction(ptb_position, pta_position, dem::smoothLength);	// this is to add SPH ptb
		    dWab_dra = periDynamics::partialKernelFunction(pta_position, ptb_position, dem::smoothLength);	// this is to add SPH pta
		    dWba_drb = periDynamics::partialKernelFunction(ptb_position, pta_position, dem::smoothLength);	// this is to add SPH ptb

		    vab = (*pta)->getVelocity()-(*ptb)->getVelocity();
		    vba = -vab;

		    (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    pa = (*pta)->getParticlePressure();
		    pb = (*ptb)->getParticlePressure();
		    mua = (*pta)->getParticleViscosity();
		    mub = (*ptb)->getParticleViscosity();
		    rhoa = (*pta)->getParticleDensity();
		    rhob = (*ptb)->getParticleDensity();

		    dva_dt = (*ptb)->getParticleMass()*( -(pa/(rhoa*rhoa)+pb/(rhob*rhob))*delta_aWab + (mua+mub)*vab/(rhoa*rhob)/rab*dWab_dra );
		    (*pta)->addVelocityDot(dva_dt);
		    dvb_dt = (*pta)->getParticleMass()*( -(pb/(rhob*rhob)+pa/(rhoa*rhoa))*delta_bWba + (mub+mua)*vba/(rhob*rhoa)/rab*dWba_drb );
		    (*ptb)->addVelocityDot(dvb_dt);
		}
	    }
	    // no rightside partition
	}

    } // end calculateDensityDotVelocityDot()

*/

    void assembly::checkDivision(){	// can be used for 2D and 3D cell division
	//print out the SPH particles in different cells
	int tmp_n = 0;
	char stepsstr[5];
	char filename[50];
//	for(std::vector< std::vector<sph::SPHParticle*> >::iterator pvec=SPHParticleCellVec.begin(); pvec=SPHParticleCellVec.end(); pvec++, tmp_n++){
	for(int pvec=0; pvec<SPHParticleCellVec.size(); pvec++, tmp_n++){

	    if(SPHParticleCellVec[pvec].size()==0){
		continue;
	    }
    	    sprintf(stepsstr,"%4d",tmp_n);
    	    strcpy(filename, "cell_");
    	    strcat(filename, stepsstr);
    	    sprintf(stepsstr,"%4d",SPHParticleCellVec[pvec].size());
	    strcat(filename,"_num");
	    strcat(filename,stepsstr);
	    strcat(filename,".dat");
    	    std::ofstream ofs(filename);
    	    if(!ofs) {
    		std::cout << "stream error!" << std::endl; exit(-1);
    	    }
    	    ofs.setf(std::ios::scientific, std::ios::floatfield);
    	    ofs.precision(dem::OPREC);

//    	    ofs << SPHParticleCellVec[pvec].size() << std::endl;
    	    ofs << std::setw(dem::OWID) << "current_x"
      		<< std::setw(dem::OWID) << "current_y"
        	<< std::setw(dem::OWID) << "current_z"
		<< std::setw(dem::OWID) << "mass"
      		<< std::setw(dem::OWID) << "density"
      		<< std::endl;
 
    	    dem::vec tmp;
    	    for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pvec].begin(); pt!= SPHParticleCellVec[pvec].end(); pt++) {
		tmp = (*pt)->getCurrPosition();
		ofs << std::setw(dem::OWID) << tmp.getx()
	    	    << std::setw(dem::OWID) << tmp.gety()
	    	    << std::setw(dem::OWID) << tmp.getz()
	    	    << std::setw(dem::OWID) << (*pt)->getParticleMass()
	    	    << std::setw(dem::OWID) << (*pt)->getParticleDensity() << std::endl;
    	    }
  
    	    ofs.close();

	}

    } // checkDivision


    void assembly::checkNeighborCells2D(){	// can only be used for 2D 
	//print out the SPH particles in different cells
	int tmp_n = 0;
	char stepsstr[5];
	char filename[50];
//	for(std::vector< std::vector<sph::SPHParticle*> >::iterator pvec=SPHParticleCellVec.begin(); pvec=SPHParticleCellVec.end(); pvec++, tmp_n++){
	for(int pvec=0; pvec<SPHParticleCellVec.size(); pvec++, tmp_n++){

	    if(SPHParticleCellVec[pvec].size()==0){
		continue;
	    }
    	    sprintf(stepsstr,"%4d",tmp_n);
    	    strcpy(filename, "neighbors_");
    	    strcat(filename, stepsstr);
    	    sprintf(stepsstr,"%4d",SPHParticleCellVec[pvec].size());
	    strcat(filename,"_num");
	    strcat(filename,stepsstr);
	    strcat(filename,".dat");
    	    std::ofstream ofs(filename);
    	    if(!ofs) {
    		std::cout << "stream error!" << std::endl; exit(-1);
    	    }
    	    ofs.setf(std::ios::scientific, std::ios::floatfield);
    	    ofs.precision(dem::OPREC);

	    ofs << "Title = \"SPH Particle Information\"" << std::endl;
	    ofs << "VARIABLES = \"x\", \"y\",\"z\" \"mass\" \"density\" " << std::endl;

	    ofs << "ZONE T =\" 0-th Load Step\" "<< std::endl;
    	    dem::vec tmp;
    	    for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pvec].begin(); pt!= SPHParticleCellVec[pvec].end(); pt++) {
		tmp = (*pt)->getCurrPosition();
		ofs << std::setw(dem::OWID) << tmp.getx()
	    	    << std::setw(dem::OWID) << tmp.gety()
	    	    << std::setw(dem::OWID) << tmp.getz()
	    	    << std::setw(dem::OWID) << (*pt)->getParticleMass()
	    	    << std::setw(dem::OWID) << (*pt)->getParticleDensity() << std::endl;
    	    }
	    if(SPHParticleCellVec[pvec+1].size()!=0   || SPHParticleCellVec[pvec+Nx-1].size()!=0 
	    || SPHParticleCellVec[pvec+Nx].size()!=0  || SPHParticleCellVec[pvec+Nx+1].size()!=0){
	    	ofs << "ZONE T =\" 1-th Load Step\" "<< std::endl;
	    }
	    if(pvec+1<numCell && SPHParticleCellVec[pvec+1].size()!=0){
    	    	for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pvec+1].begin(); pt!= SPHParticleCellVec[pvec+1].end(); pt++) {
		    tmp = (*pt)->getCurrPosition();
		    ofs << std::setw(dem::OWID) << tmp.getx()
	    	    	<< std::setw(dem::OWID) << tmp.gety()
	    	    	<< std::setw(dem::OWID) << tmp.getz()
	    	    	<< std::setw(dem::OWID) << (*pt)->getParticleMass()
	    	    	<< std::setw(dem::OWID) << (*pt)->getParticleDensity() << std::endl;
    	    	}
	    }
	    if(pvec+Nx-1<numCell && SPHParticleCellVec[pvec+Nx-1].size()!=0){
    	    	for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pvec+Nx-1].begin(); pt!= SPHParticleCellVec[pvec+Nx-1].end(); pt++) {
		    tmp = (*pt)->getCurrPosition();
		    ofs << std::setw(dem::OWID) << tmp.getx()
	    	    	<< std::setw(dem::OWID) << tmp.gety()
	    	    	<< std::setw(dem::OWID) << tmp.getz()
	    	    	<< std::setw(dem::OWID) << (*pt)->getParticleMass()
	    	    	<< std::setw(dem::OWID) << (*pt)->getParticleDensity() << std::endl;
    	    	}
	    }
	    if(pvec+Nx<numCell && SPHParticleCellVec[pvec+Nx].size()!=0){
    	    	for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pvec+Nx].begin(); pt!= SPHParticleCellVec[pvec+Nx].end(); pt++) {
		    tmp = (*pt)->getCurrPosition();
		    ofs << std::setw(dem::OWID) << tmp.getx()
	    	    	<< std::setw(dem::OWID) << tmp.gety()
	    	    	<< std::setw(dem::OWID) << tmp.getz()
	    	    	<< std::setw(dem::OWID) << (*pt)->getParticleMass()
	    	    	<< std::setw(dem::OWID) << (*pt)->getParticleDensity() << std::endl;
    	    	}
	    }
	    if(pvec+Nx+1<numCell && SPHParticleCellVec[pvec+Nx+1].size()!=0){
    	    	for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pvec+Nx+1].begin(); pt!= SPHParticleCellVec[pvec+Nx+1].end(); pt++) {
		    tmp = (*pt)->getCurrPosition();
		    ofs << std::setw(dem::OWID) << tmp.getx()
	    	    	<< std::setw(dem::OWID) << tmp.gety()
	    	    	<< std::setw(dem::OWID) << tmp.getz()
	    	    	<< std::setw(dem::OWID) << (*pt)->getParticleMass()
	    	    	<< std::setw(dem::OWID) << (*pt)->getParticleDensity() << std::endl;
    	    	}
	    }
  
    	    ofs.close();

	}

    } // checkNeighborCells


    // divide SPH domain into different cells in 2D in xz plane.
    // see the notes 5/20/2015 and 5/21/2015 or Simpson's paper "Numerical techniques for three-dimensional Smoothed Particle Hydrodynamics"
    void assembly::divideSPHDomain2D(){

	// clear the std::vector< std::vector<sph::SPHParticle*> > SPHParticleCellVec
	for(int pvec=0; pvec<SPHParticleCellVec.size(); pvec++){
	    SPHParticleCellVec[pvec].clear();	    
	}

	dem::vec tmp_xyz;
	int num;
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){	// set all free SPH particles into their cells
	    tmp_xyz = (*pt)->getCurrPosition();
	    num = int( (tmp_xyz.getz()-SPHzmin)/(3.0*dem::smoothLength) )*Nx + int( (tmp_xyz.getx()-SPHxmin)/(3.0*dem::smoothLength) );
	    SPHParticleCellVec[num].push_back(*pt);
	}
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){	// set all SPH boundary particles into their cells
	    tmp_xyz = (*pt)->getCurrPosition();
	    num = int( (tmp_xyz.getz()-SPHzmin)/(3.0*dem::smoothLength) )*Nx + int( (tmp_xyz.getx()-SPHxmin)/(3.0*dem::smoothLength) );
	    SPHParticleCellVec[num].push_back(*pt);
	}
	std::vector<sph::SPHParticle*>::iterator gt;
	for(std::vector<particle*>::iterator pdt=ParticleVec.begin(); pdt!=ParticleVec.end(); pdt++){
	    for(gt=(*pdt)->SPHGhostParticleVec.begin(); gt!=(*pdt)->SPHGhostParticleVec.end(); gt++){	// set all SPH ghost particles into their cells
	    	tmp_xyz = (*gt)->getCurrPosition();
	    	num = int( (tmp_xyz.getz()-SPHzmin)/(3.0*dem::smoothLength) )*Nx + int( (tmp_xyz.getx()-SPHxmin)/(3.0*dem::smoothLength) );
	    	SPHParticleCellVec[num].push_back(*gt);
	    }
	}

    } // divideSPHDomain2D


    // divide SPH domain into different cells in 3D in xzy plane. y is the third direction.
    // see the notes 5/20/2015 and 5/21/2015 or Simpson's paper "Numerical techniques for three-dimensional Smoothed Particle Hydrodynamics"
    void assembly::divideSPHDomain3D(){

	// clear the std::vector< std::vector<sph::SPHParticle*> > SPHParticleCellVec
	for(int pvec=0; pvec<SPHParticleCellVec.size(); pvec++){
	    SPHParticleCellVec[pvec].clear();	    
	}

	dem::vec tmp_xyz;
	int num;
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){	// set all free SPH particles into their cells
	    tmp_xyz = (*pt)->getCurrPosition();
	    num = int( (tmp_xyz.gety()-SPHymin)/(3*dem::smoothLength) )*Nx*Nz + int( (tmp_xyz.getz()-SPHzmin)/(3.0*dem::smoothLength) )*Nx + int( (tmp_xyz.getx()-SPHxmin)/(3.0*dem::smoothLength) );
	    SPHParticleCellVec[num].push_back(*pt);
	}
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){	// set all SPH boundary particles into their cells
	    tmp_xyz = (*pt)->getCurrPosition();
	    num = int( (tmp_xyz.gety()-SPHymin)/(3*dem::smoothLength) )*Nx*Nz + int( (tmp_xyz.getz()-SPHzmin)/(3.0*dem::smoothLength) )*Nx + int( (tmp_xyz.getx()-SPHxmin)/(3.0*dem::smoothLength) );
	    SPHParticleCellVec[num].push_back(*pt);
	}
	std::vector<sph::SPHParticle*>::iterator gt;
	for(std::vector<particle*>::iterator pdt=ParticleVec.begin(); pdt!=ParticleVec.end(); pdt++){
	    for(gt=(*pdt)->SPHGhostParticleVec.begin(); gt!=(*pdt)->SPHGhostParticleVec.end(); gt++){	// set all SPH ghost particles into their cells
	    	tmp_xyz = (*gt)->getCurrPosition();
	    num = int( (tmp_xyz.gety()-SPHymin)/(3*dem::smoothLength) )*Nx*Nz + int( (tmp_xyz.getz()-SPHzmin)/(3.0*dem::smoothLength) )*Nx + int( (tmp_xyz.getx()-SPHxmin)/(3.0*dem::smoothLength) );
	    	SPHParticleCellVec[num].push_back(*gt);
	    }
	}

    } // divideSPHDomain3D


#ifdef OPENMP

    //// the momentum equilibrium equation and state equation are implemented as in Monaghan's paper (1994), simulate free surface flow using sph
    //// here the neighboring list of SPH particles is searched by the cells, vector< vector<sph::SPHParticle*> > SPHParticleCellVec;
    void assembly::calculateSPHDensityDotVelocityDotLinkedList2D(){

	// divide the SPH domain into different cells, each cell will contain SPH particles within it
	divideSPHDomain2D();

//	checkDivision();	// pass, May 22, 2015
//	checkNeighborCells2D();	// pass, May 22, 2015

	// initialize the densityDot and velocityDot of all the SPH particles 
	dem::vec tmp_vec = dem::vec(0,0, -(dem::G)*(dem::GRVT_SCL) );
	dem::vec zero_vec = dem::vec(0,0,0);
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    (*pt)->setDensityDotVelocityDotZero();
	    (*pt)->addVelocityDot(tmp_vec);
	    (*pt)->calculateParticleViscosity();	// everytime when use getParticleViscolisity(), make sure that pressure has been calculated!!!!
	    (*pt)->calculateParticlePressure();		// everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
	}
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	    (*pt)->setDensityDotVelocityDotZero();
	    (*pt)->addVelocityDot(zero_vec);
	    (*pt)->calculateParticleViscosity();	// everytime when use getParticleViscolisity(), make sure that pressure has been calculated!!!!
	    (*pt)->calculateParticlePressure();		// everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
	}
	std::vector<sph::SPHParticle*>::iterator gt;
	for(std::vector<particle*>::iterator pdt=ParticleVec.begin(); pdt!=ParticleVec.end(); pdt++){
	    for(gt=(*pdt)->SPHGhostParticleVec.begin(); gt!=(*pdt)->SPHGhostParticleVec.end(); gt++){
	    	(*gt)->setDensityDotVelocityDotZero();
	    	(*gt)->addVelocityDot(tmp_vec);
	    	(*gt)->calculateParticleViscosity();	// everytime when use getParticleViscolisity(), make sure that pressure has been calculated!!!!
	    	(*gt)->calculateParticlePressure();	// everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
	    }
	}

	// temporary variables used in the loop
	dem::vec pta_position;	// ...
	REAL pa, rhoa, mua;	// pa, mua, rhoa ...
	REAL alpha = 0.0;	// ...
	REAL epsilon = 0.0; //...	// parameter for velocity correction
	REAL Ra;	// Ra ...

	int numThisCell;	// number of SPH particles in the current pta cell
	int numNeighborCell;	// number of SPH particles in the neighboring cells

	REAL ptaDensityDot;
	REAL ptaVelocityDotX, ptaVelocityDotY, ptaVelocityDotZ;
	REAL ptaVelocityCorrectionX, ptaVelocityCorrectionY, ptaVelocityCorrectionZ;

	std::vector<sph::SPHParticle*> tmp_particleVec;	// sph particles in neighboring cells
	for(int pvec=0; pvec<SPHParticleCellVec.size(); pvec++){

	    // store all SPH particles in pvec's neighboring cells
	    tmp_particleVec.clear();	// it is the same for the particles in the same cell
	    int pnum;
	    if(pvec+1<numCell){
		pnum = pvec+1;
    	    	for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pnum].begin(); pt!= SPHParticleCellVec[pnum].end(); pt++) {
		    tmp_particleVec.push_back(*pt);
    	    	}
	    }
	    if(pvec+Nx-1<numCell){
		pnum = pvec+Nx-1;
    	    	for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pnum].begin(); pt!= SPHParticleCellVec[pnum].end(); pt++) {
		    tmp_particleVec.push_back(*pt);
    	    	}
	    }
	    if(pvec+Nx<numCell){
		pnum = pvec+Nx;
    	    	for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pnum].begin(); pt!= SPHParticleCellVec[pnum].end(); pt++) {
		    tmp_particleVec.push_back(*pt);
    	    	}
	    }
	    if(pvec+Nx+1<numCell){
		pnum = pvec+Nx+1;
    	    	for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pnum].begin(); pt!= SPHParticleCellVec[pnum].end(); pt++) {
		    tmp_particleVec.push_back(*pt);
    	    	}
	    }

	    numThisCell = SPHParticleCellVec[pvec].size();

	    for(int pta=0; pta<SPHParticleCellVec[pvec].size(); pta++){	// SPH particles in cell pvec
		pa = SPHParticleCellVec[pvec][pta]->getParticlePressure();
	    	if(pa>=0){
		    Ra = 0.006;
	    	}
	    	else{
		    Ra = 0.6;
	    	}
//	    	mua = SPHParticleCellVec[pvec][pta]->getParticleViscosity();
	    	rhoa = SPHParticleCellVec[pvec][pta]->getParticleDensity();
	    	pta_position = SPHParticleCellVec[pvec][pta]->getCurrPosition();

//		if(SPHParticleCellVec[pvec][pta]->getType()!=1){	// pta is not free sph particles, i.e. pta is ghost or boundary particles
//		    continue;	// do not consider pta, as we treat before
//		}

		// initialize the pta variables
	    	ptaDensityDot = 0;
	    	ptaVelocityDotX = 0; ptaVelocityDotY = 0; ptaVelocityDotZ = 0;
	    	ptaVelocityCorrectionX = 0; ptaVelocityCorrectionY = 0; ptaVelocityCorrectionZ = 0;
	    	numNeighborCell = tmp_particleVec.size();		

#pragma omp parallel num_threads(NUM_THREADS) reduction(+:ptaDensityDot, ptaVelocityDotX, ptaVelocityDotY, ptaVelocityDotZ, ptaVelocityCorrectionX, ptaVelocityCorrectionY, ptaVelocityCorrectionZ)

{



		// temporary variables used in the loop

		dem::vec ptb_position;	
		dem::vec delta_aWab;	
		dem::vec delta_bWba;	
		dem::vec vab;	
		dem::vec vba;	
		dem::vec vdem;	
		dem::vec dva_dt;	
		dem::vec dvb_dt;	
		dem::vec delta_a;	
		dem::vec delta_b;	
		REAL pb, rhob, mub;	
		REAL rab;	
        	REAL dWab_dra;	
		REAL dWba_drb;	
		REAL Wab, Wba;	
		REAL da, dB;	
		REAL beta;	
		REAL xa, ya, xB, yB, k, sqrt_xaya; 	// variables for Morris' method to calculate da/dB	
		int ptb;	
		REAL ra, rb, rc; 	// the geometry of the dem particle
		dem::vec pta_local, ptb_local; 	// local position of pta and ptb in the dem particle
		dem::vec demt_curr;	
		REAL Gamma_ab, mu_ab, vr_dot;	
		REAL Wq, Rb, phi_4, coefficient_a, coefficient_b;	
		dem::particle* demt;	

		int thread_num  = omp_get_thread_num();

		int start_num = BLOCK_LOW(thread_num, NUM_THREADS, numThisCell-pta-1);
		int end_num   = BLOCK_HIGH(thread_num, NUM_THREADS, numThisCell-pta-1);
		for(ptb=start_num+pta+1; ptb<end_num+pta+1+1; ptb++){	// sum over the SPH particles in the same cell pvec
		    ptb_position = SPHParticleCellVec[pvec][ptb]->getCurrPosition();
		    rab = dem::vfabs(pta_position-ptb_position);
		    if(rab<=3*dem::smoothLength){	// ptb is in the smooth kernel
		    	pb = SPHParticleCellVec[pvec][ptb]->getParticlePressure();
//		    	mub = SPHParticleCellVec[pvec][ptb]->getParticleViscosity();
		    	rhob = SPHParticleCellVec[pvec][ptb]->getParticleDensity();

		    	Wq = kernelFunction(rab/dem::smoothLength);
		    	phi_4 = pow(Wq/Wqmin,4);
		    	if(pb>=0){
			    Rb = 0.006;
		    	}
		    	else{
			    Rb = 0.6;
		    	}
		    	coefficient_a = 1+Ra*phi_4;
		    	coefficient_b = 1+Rb*phi_4;

			// we have three types of SPH particles: 1, free particle; 2, ghost particle; 3, boundary particle
			// Then we have 3x3 = 9 different types of interactions with the three types of particles
			// so we cannot judge the interaction type only by the type of particle ptb, we need to consider pta also
			// pta          ptb          need to consider or not
			// 1            1              V			free with free
			// 1            2              V			free with ghost
			// 1            3              V			free with boundary

			// 2            1              V			ghost with free
			// 2            2              X			ghost with ghost
			// 2            3              X			ghost with boundary

			// 3            1              V			boundary with free
			// 3            2              X			boundary with ghost
			// 3            3              X 			boundary with boundary


			// add the density dot for pta and ptb
			delta_aWab = gradientKernelFunction(pta_position, ptb_position);	// this is to add SPH pta 
		    	delta_bWba = gradientKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
		    	dWab_dra = partialKernelFunction(pta_position, ptb_position);	// this is to add SPH pta
		    	dWba_drb = partialKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
		    	Wab = kernelFunction(pta_position, ptb_position);	// this is to add SPH pta 
		    	Wba = kernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
			switch(SPHParticleCellVec[pvec][ptb]->getType()){
			  case 1:	// ptb is free SPH particle
			    switch(SPHParticleCellVec[pvec][pta]->getType()){
			      case 1:	// free with free
		    	    	vab = SPHParticleCellVec[pvec][pta]->getVelocity()-SPHParticleCellVec[pvec][ptb]->getVelocity();
		    	   	vba = -vab;

		    	    	ptaDensityDot += SPHParticleCellVec[pvec][ptb]->getParticleMass()*(vab%delta_aWab);
		    	    	SPHParticleCellVec[pvec][ptb]->addDensityDot( SPHParticleCellVec[pvec][pta]->getParticleMass()*(vba%delta_bWba) );

		    	    	vr_dot = vab%(pta_position-ptb_position);
		    	    	if(vr_dot<0){
				    mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    		    Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    	}
		    	    	else{
		 		    Gamma_ab = 0;
		    	    	}


		    	    	dva_dt = -SPHParticleCellVec[pvec][ptb]->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
				ptaVelocityDotX += dva_dt.getx(); ptaVelocityDotY += dva_dt.gety(); ptaVelocityDotZ += dva_dt.getz();
		    	    	dvb_dt = -SPHParticleCellVec[pvec][pta]->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
		    	    	SPHParticleCellVec[pvec][ptb]->addVelocityDot(dvb_dt);

		    	    	delta_a = epsilon*SPHParticleCellVec[pvec][ptb]->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
				ptaVelocityCorrectionX += delta_a.getx(); ptaVelocityCorrectionY += delta_a.gety(); ptaVelocityCorrectionZ += delta_a.getz();
		    	    	delta_b = epsilon*SPHParticleCellVec[pvec][pta]->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
		    	    	SPHParticleCellVec[pvec][ptb]->addVelocityCorrection(delta_b);
			        break;
			      case 2:	// ghost with free
				demt = SPHParticleCellVec[pvec][pta]->getDemParticle();
				demt_curr = demt->getCurrPosition();
			    	ptb_local = demt->localVec(ptb_position-demt_curr);	// the local position of sph point ptb
			    	ra = demt->getA(); rb = demt->getB(); rc = demt->getC();
			    	k = 1.0/(sqrt(ptb_local.getx()*ptb_local.getx()/(ra*ra)+ptb_local.gety()*ptb_local.gety()/(rb*rb)+ptb_local.getz()*ptb_local.getz()/(rc*rc) ) );
			    	da = dem::vfabs(ptb_local-k*ptb_local);	// the distance is the same in rotated coordinates

		    	    	// calculate Vab as the method shown in Morris's paper, 1996

		    	    	// (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
			    	pta_local = demt->localVec(pta_position-demt_curr);
			    	k = 1.0/(sqrt(pta_local.getx()*pta_local.getx()/(ra*ra)+pta_local.gety()*pta_local.gety()/(rb*rb)+pta_local.getz()*pta_local.getz()/(rc*rc) ) );
		    	    	dB = dem::vfabs(pta_local-k*pta_local);

//		    	    // (2) here the Morris's method is used 
//		    	    xa = pta_position.getx(); ya = pta_position.gety();
//		    	    xB = ptb_position.getx(); yB = ptb_position.gety();
//		    	    if(ya==0) {da = xa-radius; dB = radius-xB;}
//		    	    else if(xa==0) {da = ya-radius; dB = radius-yB;}
//		    	    else {
//			    sqrt_xaya = sqrt(xa*xa+ya*ya);
//			    k = radius/sqrt_xaya;
//			    da = radius/k - radius;
//			    dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//		    	    }


		    	    	beta = 1+dB/da;
		    	    	if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

			    	vdem = demt->getCurrVelocity()+demt->getCurrOmga()*( pta_position-demt_curr );
		    	    	vba = beta*(SPHParticleCellVec[pvec][ptb]->getVelocity()-vdem);
		    	    	vab = -vba;

		    	    	ptaDensityDot += SPHParticleCellVec[pvec][ptb]->getParticleMass()*(vab%delta_aWab);
		    	    	SPHParticleCellVec[pvec][ptb]->addDensityDot( SPHParticleCellVec[pvec][pta]->getParticleMass()*(vba%delta_bWba) );

		    	    	vr_dot = vab%(pta_position-ptb_position);
		    	    	if(vr_dot<0){
			    	    mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    	    	    Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    	}
		    	    	else{
		 	    	    Gamma_ab = 0;
		    	    	}

		    	    	dva_dt = -SPHParticleCellVec[pvec][ptb]->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
				#pragma omp critical (demtCritical)
			    	{
				demt->addForce(SPHParticleCellVec[pvec][pta]->getParticleMass()*dva_dt);
			    	demt->addMoment( (pta_position-demt_curr) * SPHParticleCellVec[pvec][pta]->getParticleMass()*dva_dt );
				}
//		    	    	SPHParticleCellVec[pvec][pta]->addVelocityDot(dva_dt);

		    	    	dvb_dt = -SPHParticleCellVec[pvec][pta]->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
		    	    	SPHParticleCellVec[pvec][ptb]->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

//		            	delta_a = epsilon*SPHParticleCellVec[pvec][ptb]->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
//		            	SPHParticleCellVec[pvec][pta]->addVelocityCorrection(delta_a);
		            	delta_b = epsilon*SPHParticleCellVec[pvec][pta]->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
		            	SPHParticleCellVec[pvec][ptb]->addVelocityCorrection(delta_b);
				break;
			      case 3:	// boundary with free

	    		    	// interact with boundary particles
	    		    	da = ptb_position.getz();	// assume with the bottom boundary
	    		    	if(ptb_position.getx() < da){	// with left boundary
				    da = ptb_position.getx();
	    		    	} 
	    		    	else if(length_x-ptb_position.getx() < da){	// with right boundary
				    da = length_x-ptb_position.getx();
	    		    	}
	    		    	else if(height_z-ptb_position.getz() < da){	// with top boundary
				    da = height_z-ptb_position.getz();
	    		    	}

		    	    	// calculate Vab as the method shown in Morris's paper, 1996

	    	    	    	dB = 0-pta_position.getz();	// assume with the bottom boundary
	    	    	    	if(pta_position.getx() < dB){	// with left boundary
		    		    dB = 0-pta_position.getx();
	    	    	    	} 
	    	    	    	else if(pta_position.getx()-length_x < dB){	// with right boundary
		    		    dB = pta_position.getx()-length_x;
	    	    	    	}
	    	    	    	else if(pta_position.getz()-height_z < dB){	// with top boundary
		    		    dB = pta_position.getz()-height_z;
	    	    	    	}

		    	    	beta = 1+dB/da;
		    	    	if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

		    	    	vba = beta*SPHParticleCellVec[pvec][ptb]->getVelocity();
		    	    	vab = -vba;

		    	    	ptaDensityDot += SPHParticleCellVec[pvec][ptb]->getParticleMass()*(vab%delta_aWab);
		    	    	SPHParticleCellVec[pvec][ptb]->addDensityDot( SPHParticleCellVec[pvec][pta]->getParticleMass()*(vba%delta_bWba) );

		    	    	vr_dot = vab%(pta_position-ptb_position);
		    	    	if(vr_dot<0){
				    mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    		    Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    	}
		    	    	else{
		 		    Gamma_ab = 0;
		    	    	}

//		    	    	dva_dt = -SPHParticleCellVec[pvec][ptb]->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
//		    	    	SPHParticleCellVec[pvec][pta]->addVelocityDot(dva_dt);
		    	    	dvb_dt = -SPHParticleCellVec[pvec][pta]->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
		    	    	SPHParticleCellVec[pvec][ptb]->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

//			    	delta_a = epsilon*SPHParticleCellVec[pvec][ptb]->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
//		    	    	SPHParticleCellVec[pvec][pta]->addVelocityCorrection(delta_a);
		    	    	delta_b = epsilon*SPHParticleCellVec[pvec][pta]->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
		    	    	SPHParticleCellVec[pvec][ptb]->addVelocityCorrection(delta_b);
			
//			    	// apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
//		    	    	if(rab<=space_interval){ // ptb is in the smooth kernel
//		    		    dvb_dt = D*(pow(space_interval/rab, p1)-pow(space_interval/rab, p2))*(ptb_position-pta_position)/(rab*rab);
//		    		    ptaVelocityDotX += dva_dt.getx(); ptaVelocityDotY += dva_dt.gety(); ptaVelocityDotZ += dva_dt.getz();
//		    	    	} // end if

				break;
			      default:
			    	std::cout << "SPH particle type of pta should be 1, 2 or 3!" << std::endl;
			    	exit(-1);
			      } // end switch pta

			    break;
			  case 2:	// ptb is ghost particle

			    if(SPHParticleCellVec[pvec][pta]->getType()!=1){	// pta is not free sph particles, i.e. pta is ghost or boundary particles
		    		break;	// do not consider pta, as we treat before
			    }
			    demt = SPHParticleCellVec[pvec][ptb]->getDemParticle();
			    demt_curr = demt->getCurrPosition();
			    pta_local = demt->localVec(pta_position-demt_curr);	// the local position of sph point pta
			    ra = demt->getA(); rb = demt->getB(); rc = demt->getC();
			    k = 1.0/(sqrt(pta_local.getx()*pta_local.getx()/(ra*ra)+pta_local.gety()*pta_local.gety()/(rb*rb)+pta_local.getz()*pta_local.getz()/(rc*rc) ) );
			    da = dem::vfabs(pta_local-k*pta_local);	// the distance is the same in rotated coordinates

		    	    // calculate Vab as the method shown in Morris's paper, 1996

		    	    // (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
			    ptb_local = demt->localVec(ptb_position-demt_curr);
			    k = 1.0/(sqrt(ptb_local.getx()*ptb_local.getx()/(ra*ra)+ptb_local.gety()*ptb_local.gety()/(rb*rb)+ptb_local.getz()*ptb_local.getz()/(rc*rc) ) );
		    	    dB = dem::vfabs(ptb_local-k*ptb_local);

//		    	    // (2) here the Morris's method is used 
//		    	    xa = pta_position.getx(); ya = pta_position.gety();
//		    	    xB = ptb_position.getx(); yB = ptb_position.gety();
//		    	    if(ya==0) {da = xa-radius; dB = radius-xB;}
//		    	    else if(xa==0) {da = ya-radius; dB = radius-yB;}
//		    	    else {
//			    sqrt_xaya = sqrt(xa*xa+ya*ya);
//			    k = radius/sqrt_xaya;
//			    da = radius/k - radius;
//			    dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//		    	    }


		    	    beta = 1+dB/da;
		    	    if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

			    vdem = demt->getCurrVelocity()+demt->getCurrOmga()*( ptb_position-demt_curr );
		    	    vab = beta*(SPHParticleCellVec[pvec][pta]->getVelocity()-vdem);
		    	    vba = -vab;

		    	    ptaDensityDot += SPHParticleCellVec[pvec][ptb]->getParticleMass()*(vab%delta_aWab);
		    	    SPHParticleCellVec[pvec][ptb]->addDensityDot( SPHParticleCellVec[pvec][pta]->getParticleMass()*(vba%delta_bWba) );

		    	    vr_dot = vab%(pta_position-ptb_position);
		    	    if(vr_dot<0){
			    	mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    	    	Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    }
		    	    else{
		 	    	Gamma_ab = 0;
		    	    }

		    	    dva_dt = -SPHParticleCellVec[pvec][ptb]->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
		    	    ptaVelocityDotX += dva_dt.getx(); ptaVelocityDotY += dva_dt.gety(); ptaVelocityDotZ += dva_dt.getz();
		    	    dvb_dt = -SPHParticleCellVec[pvec][pta]->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
			    #pragma omp critical (demtCritical)
			    {
			    demt->addForce(SPHParticleCellVec[pvec][ptb]->getParticleMass()*dvb_dt);
			    demt->addMoment( (ptb_position-demt_curr) * SPHParticleCellVec[pvec][ptb]->getParticleMass()*dvb_dt );
			    }
//		    	    SPHParticleCellVec[pvec][ptb]->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

		            delta_a = epsilon*SPHParticleCellVec[pvec][ptb]->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
		            ptaVelocityCorrectionX += delta_a.getx(); ptaVelocityCorrectionY += delta_a.gety(); ptaVelocityCorrectionZ += delta_a.getz();
//		            delta_b = epsilon*SPHParticleCellVec[pvec][pta]->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//		            SPHParticleCellVec[pvec][ptb]->addVelocityCorrection(delta_b);

			    break;
			  case 3:	// ptb is boundary particle

			    if(SPHParticleCellVec[pvec][pta]->getType()!=1){	// pta is not free sph particles, i.e. pta is ghost or boundary particles
		    		break;	// do not consider pta, as we treat before
			    }

	    		    // interact with boundary particles
	    		    da = pta_position.getz();	// assume with the bottom boundary
	    		    if(pta_position.getx() < da){	// with left boundary
				da = pta_position.getx();
	    		    } 
	    		    else if(length_x-pta_position.getx() < da){	// with right boundary
				da = length_x-pta_position.getx();
	    		    }
	    		    else if(height_z-pta_position.getz() < da){	// with top boundary
				da = height_z-pta_position.getz();
	    		    }

		    	    // calculate Vab as the method shown in Morris's paper, 1996

	    	    	    dB = 0-ptb_position.getz();	// assume with the bottom boundary
	    	    	    if(ptb_position.getx() < dB){	// with left boundary
		    		dB = 0-ptb_position.getx();
	    	    	    } 
	    	    	    else if(ptb_position.getx()-length_x < dB){	// with right boundary
		    		dB = ptb_position.getx()-length_x;
	    	    	    }
	    	    	    else if(ptb_position.getz()-height_z < dB){	// with top boundary
		    		dB = ptb_position.getz()-height_z;
	    	    	    }

		    	    beta = 1+dB/da;
		    	    if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

		    	    vab = beta*SPHParticleCellVec[pvec][pta]->getVelocity();
		    	    vba = -vab;

		    	    ptaDensityDot += SPHParticleCellVec[pvec][ptb]->getParticleMass()*(vab%delta_aWab);
		    	    SPHParticleCellVec[pvec][ptb]->addDensityDot( SPHParticleCellVec[pvec][pta]->getParticleMass()*(vba%delta_bWba) );

		    	    vr_dot = vab%(pta_position-ptb_position);
		    	    if(vr_dot<0){
				mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    		Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    }
		    	    else{
		 		Gamma_ab = 0;
		    	    }

		    	    dva_dt = -SPHParticleCellVec[pvec][ptb]->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
		    	    ptaVelocityDotX += dva_dt.getx(); ptaVelocityDotY += dva_dt.gety(); ptaVelocityDotZ += dva_dt.getz();
//		    	    dvb_dt = -SPHParticleCellVec[pvec][pta]->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
//		    	    SPHParticleCellVec[pvec][ptb]->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

		    	    delta_a = epsilon*SPHParticleCellVec[pvec][ptb]->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
		    	    ptaVelocityCorrectionX += delta_a.getx(); ptaVelocityCorrectionY += delta_a.gety(); ptaVelocityCorrectionZ += delta_a.getz(); 
// //		    	    delta_b = epsilon*SPHParticleCellVec[pvec][pta]->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//		    	    SPHParticleCellVec[pvec][ptb]->addVelocityCorrection(delta_b);
			
//			    // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
//		    	    if(rab<=space_interval){ // ptb is in the smooth kernel
//		    		dva_dt = D*(pow(space_interval/rab, p1)-pow(space_interval/rab, p2))*(pta_position-ptb_position)/(rab*rab);
//		    		ptaVelocityDotX += dva_dt.getx(); ptaVelocityDotY += dva_dt.gety(); ptaVelocityDotZ += dva_dt.getz();
//		    	    } // end if
			    break;
			  default:
			    std::cout << "SPH particle type should be 1, 2 or 3!" << std::endl;
			    exit(-1);
		
			} //end swtich type
		    } // end if 3h
		} // end for ptb in the same cell


		start_num = BLOCK_LOW(thread_num, NUM_THREADS, numNeighborCell);
		end_num   = BLOCK_HIGH(thread_num, NUM_THREADS, numNeighborCell);
		for(ptb=start_num; ptb<end_num+1; ptb++){	// all particles in pvec's neighboring cells 
		    ptb_position = tmp_particleVec[ptb]->getCurrPosition();
		    rab = dem::vfabs(pta_position-ptb_position);
		    if(rab<=3*dem::smoothLength){	// ptb is in the smooth kernel
		    	pb = tmp_particleVec[ptb]->getParticlePressure();
//		    	mub = tmp_particleVec[ptb]->getParticleViscosity();
		    	rhob = tmp_particleVec[ptb]->getParticleDensity();

		    	Wq = kernelFunction(rab/dem::smoothLength);
		    	phi_4 = pow(Wq/Wqmin,4);
		    	if(pb>=0){
			    Rb = 0.006;
		    	}
		    	else{
			    Rb = 0.6;
		    	}
		    	coefficient_a = 1+Ra*phi_4;
		    	coefficient_b = 1+Rb*phi_4;

			// we have three types of SPH particles: 1, free particle; 2, ghost particle; 3, boundary particle
			// Then we have 3x3 = 9 different types of interactions with the three types of particles
			// so we cannot judge the interaction type only by the type of particle ptb, we need to consider pta also
			// pta          ptb          need to consider or not
			// 1            1              V			free with free
			// 1            2              V			free with ghost
			// 1            3              V			free with boundary

			// 2            1              V			ghost with free
			// 2            2              X			ghost with ghost
			// 2            3              X			ghost with boundary

			// 3            1              V			boundary with free
			// 3            2              X			boundary with ghost
			// 3            3              X 			boundary with boundary


			// add the density dot for pta and ptb
			delta_aWab = gradientKernelFunction(pta_position, ptb_position);	// this is to add SPH pta 
		    	delta_bWba = gradientKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
		    	dWab_dra = partialKernelFunction(pta_position, ptb_position);	// this is to add SPH pta
		    	dWba_drb = partialKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
		    	Wab = kernelFunction(pta_position, ptb_position);	// this is to add SPH pta 
		    	Wba = kernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
			switch(tmp_particleVec[ptb]->getType()){
			  case 1:	// ptb is free SPH particle
			    switch(SPHParticleCellVec[pvec][pta]->getType()){
			      case 1:	// free with free
		    	    	vab = SPHParticleCellVec[pvec][pta]->getVelocity()-tmp_particleVec[ptb]->getVelocity();
		    	   	vba = -vab;

		    	    	ptaDensityDot += tmp_particleVec[ptb]->getParticleMass()*(vab%delta_aWab);
		    	    	tmp_particleVec[ptb]->addDensityDot( SPHParticleCellVec[pvec][pta]->getParticleMass()*(vba%delta_bWba) );

		    	    	vr_dot = vab%(pta_position-ptb_position);
		    	    	if(vr_dot<0){
				    mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    		    Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    	}
		    	    	else{
		 		    Gamma_ab = 0;
		    	    	}


		    	    	dva_dt = -tmp_particleVec[ptb]->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
			        ptaVelocityDotX += dva_dt.getx(); ptaVelocityDotY += dva_dt.gety(); ptaVelocityDotZ += dva_dt.getz();
		    	    	dvb_dt = -SPHParticleCellVec[pvec][pta]->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
		    	    	tmp_particleVec[ptb]->addVelocityDot(dvb_dt);

		    	    	delta_a = epsilon*tmp_particleVec[ptb]->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
				ptaVelocityCorrectionX += delta_a.getx(); ptaVelocityCorrectionY += delta_a.gety(); ptaVelocityCorrectionZ += delta_a.getz(); 
		    	    	delta_b = epsilon*SPHParticleCellVec[pvec][pta]->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
		    	    	tmp_particleVec[ptb]->addVelocityCorrection(delta_b);
			        break;
			      case 2:	// ghost with free
				demt = SPHParticleCellVec[pvec][pta]->getDemParticle();
				demt_curr = demt->getCurrPosition();
			    	ptb_local = demt->localVec(ptb_position-demt_curr);	// the local position of sph point ptb
			    	ra = demt->getA(); rb = demt->getB(); rc = demt->getC();
			    	k = 1.0/(sqrt(ptb_local.getx()*ptb_local.getx()/(ra*ra)+ptb_local.gety()*ptb_local.gety()/(rb*rb)+ptb_local.getz()*ptb_local.getz()/(rc*rc) ) );
			    	da = dem::vfabs(ptb_local-k*ptb_local);	// the distance is the same in rotated coordinates

		    	    	// calculate Vab as the method shown in Morris's paper, 1996

		    	    	// (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
			    	pta_local = demt->localVec(pta_position-demt_curr);
			    	k = 1.0/(sqrt(pta_local.getx()*pta_local.getx()/(ra*ra)+pta_local.gety()*pta_local.gety()/(rb*rb)+pta_local.getz()*pta_local.getz()/(rc*rc) ) );
		    	    	dB = dem::vfabs(pta_local-k*pta_local);

//		    	    // (2) here the Morris's method is used 
//		    	    xa = pta_position.getx(); ya = pta_position.gety();
//		    	    xB = ptb_position.getx(); yB = ptb_position.gety();
//		    	    if(ya==0) {da = xa-radius; dB = radius-xB;}
//		    	    else if(xa==0) {da = ya-radius; dB = radius-yB;}
//		    	    else {
//			    sqrt_xaya = sqrt(xa*xa+ya*ya);
//			    k = radius/sqrt_xaya;
//			    da = radius/k - radius;
//			    dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//		    	    }


		    	    	beta = 1+dB/da;
		    	    	if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

			    	vdem = demt->getCurrVelocity()+demt->getCurrOmga()*( pta_position-demt_curr );
		    	    	vba = beta*(tmp_particleVec[ptb]->getVelocity()-vdem);
		    	    	vab = -vba;

		    	    	ptaDensityDot += tmp_particleVec[ptb]->getParticleMass()*(vab%delta_aWab);
		    	    	tmp_particleVec[ptb]->addDensityDot( SPHParticleCellVec[pvec][pta]->getParticleMass()*(vba%delta_bWba) );

		    	    	vr_dot = vab%(pta_position-ptb_position);
		    	    	if(vr_dot<0){
			    	    mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    	    	    Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    	}
		    	    	else{
		 	    	    Gamma_ab = 0;
		    	    	}

		    	    	dva_dt = -tmp_particleVec[ptb]->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
				#pragma omp critical (demtCritical)
			    	{
				demt->addForce(SPHParticleCellVec[pvec][pta]->getParticleMass()*dva_dt);
			    	demt->addMoment( (pta_position-demt_curr) * SPHParticleCellVec[pvec][pta]->getParticleMass()*dva_dt );
				}
//		    	    	SPHParticleCellVec[pvec][pta]->addVelocityDot(dva_dt);

		    	    	dvb_dt = -SPHParticleCellVec[pvec][pta]->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
		    	    	tmp_particleVec[ptb]->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

//		            	delta_a = epsilon*tmp_particleVec[ptb]->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
//		            	SPHParticleCellVec[pvec][pta]->addVelocityCorrection(delta_a);
		            	delta_b = epsilon*SPHParticleCellVec[pvec][pta]->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
		            	tmp_particleVec[ptb]->addVelocityCorrection(delta_b);
				break;
			      case 3:	// boundary with free

	    		    	// interact with boundary particles
	    		    	da = ptb_position.getz();	// assume with the bottom boundary
	    		    	if(ptb_position.getx() < da){	// with left boundary
				    da = ptb_position.getx();
	    		    	} 
	    		    	else if(length_x-ptb_position.getx() < da){	// with right boundary
				    da = length_x-ptb_position.getx();
	    		    	}
	    		    	else if(height_z-ptb_position.getz() < da){	// with top boundary
				    da = height_z-ptb_position.getz();
	    		    	}

		    	    	// calculate Vab as the method shown in Morris's paper, 1996

	    	    	    	dB = 0-pta_position.getz();	// assume with the bottom boundary
	    	    	    	if(pta_position.getx() < dB){	// with left boundary
		    		    dB = 0-pta_position.getx();
	    	    	    	} 
	    	    	    	else if(pta_position.getx()-length_x < dB){	// with right boundary
		    		    dB = pta_position.getx()-length_x;
	    	    	    	}
	    	    	    	else if(pta_position.getz()-height_z < dB){	// with top boundary
		    		    dB = pta_position.getz()-height_z;
	    	    	    	}

		    	    	beta = 1+dB/da;
		    	    	if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

		    	    	vba = beta*tmp_particleVec[ptb]->getVelocity();
		    	    	vab = -vba;

		    	    	ptaDensityDot += tmp_particleVec[ptb]->getParticleMass()*(vab%delta_aWab);
		    	    	tmp_particleVec[ptb]->addDensityDot( SPHParticleCellVec[pvec][pta]->getParticleMass()*(vba%delta_bWba) );

		    	    	vr_dot = vab%(pta_position-ptb_position);
		    	    	if(vr_dot<0){
				    mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    		    Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    	}
		    	    	else{
		 		    Gamma_ab = 0;
		    	    	}

//		    	    	dva_dt = -tmp_particleVec[ptb]->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
//		    	    	SPHParticleCellVec[pvec][pta]->addVelocityDot(dva_dt);
		    	    	dvb_dt = -SPHParticleCellVec[pvec][pta]->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
		    	    	tmp_particleVec[ptb]->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

//			    	delta_a = epsilon*tmp_particleVec[ptb]->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
//		    	    	SPHParticleCellVec[pvec][pta]->addVelocityCorrection(delta_a);
		    	    	delta_b = epsilon*SPHParticleCellVec[pvec][pta]->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
		    	    	tmp_particleVec[ptb]->addVelocityCorrection(delta_b);
			
//			    	// apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
//		    	    	if(rab<=space_interval){ // ptb is in the smooth kernel
//		    		    dvb_dt = D*(pow(space_interval/rab, p1)-pow(space_interval/rab, p2))*(ptb_position-pta_position)/(rab*rab);
//		    		    ptaVelocityDotX += dva_dt.getx(); ptaVelocityDotY += dva_dt.gety(); ptaVelocityDotZ += dva_dt.getz();
//		    	    	} // end if

				break;
			      default:
			    	std::cout << "SPH particle type of pta should be 1, 2 or 3!" << std::endl;
			    	exit(-1);
			      } // end switch pta

			    break;
			  case 2:	// ptb is ghost particle

			    if(SPHParticleCellVec[pvec][pta]->getType()!=1){	// pta is not free sph particles, i.e. pta is ghost or boundary particles
		    		break;	// do not consider pta, as we treat before
			    }
			    demt = tmp_particleVec[ptb]->getDemParticle();
			    demt_curr = demt->getCurrPosition();
			    pta_local = demt->localVec(pta_position-demt_curr);	// the local position of sph point pta
			    ra = demt->getA(); rb = demt->getB(); rc = demt->getC();
			    k = 1.0/(sqrt(pta_local.getx()*pta_local.getx()/(ra*ra)+pta_local.gety()*pta_local.gety()/(rb*rb)+pta_local.getz()*pta_local.getz()/(rc*rc) ) );
			    da = dem::vfabs(pta_local-k*pta_local);	// the distance is the same in rotated coordinates

		    	    // calculate Vab as the method shown in Morris's paper, 1996

		    	    // (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
			    ptb_local = demt->localVec(ptb_position-demt_curr);
			    k = 1.0/(sqrt(ptb_local.getx()*ptb_local.getx()/(ra*ra)+ptb_local.gety()*ptb_local.gety()/(rb*rb)+ptb_local.getz()*ptb_local.getz()/(rc*rc) ) );
		    	    dB = dem::vfabs(ptb_local-k*ptb_local);

//		    	    // (2) here the Morris's method is used 
//		    	    xa = pta_position.getx(); ya = pta_position.gety();
//		    	    xB = ptb_position.getx(); yB = ptb_position.gety();
//		    	    if(ya==0) {da = xa-radius; dB = radius-xB;}
//		    	    else if(xa==0) {da = ya-radius; dB = radius-yB;}
//		    	    else {
//			    sqrt_xaya = sqrt(xa*xa+ya*ya);
//			    k = radius/sqrt_xaya;
//			    da = radius/k - radius;
//			    dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//		    	    }


		    	    beta = 1+dB/da;
		    	    if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

			    vdem = demt->getCurrVelocity()+demt->getCurrOmga()*( ptb_position-demt_curr );
		    	    vab = beta*(SPHParticleCellVec[pvec][pta]->getVelocity()-vdem);
		    	    vba = -vab;

		    	    ptaDensityDot += tmp_particleVec[ptb]->getParticleMass()*(vab%delta_aWab);
		    	    tmp_particleVec[ptb]->addDensityDot( SPHParticleCellVec[pvec][pta]->getParticleMass()*(vba%delta_bWba) );

		    	    vr_dot = vab%(pta_position-ptb_position);
		    	    if(vr_dot<0){
			    	mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    	    	Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    }
		    	    else{
		 	    	Gamma_ab = 0;
		    	    }

		    	    dva_dt = -tmp_particleVec[ptb]->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
			    ptaVelocityDotX += dva_dt.getx(); ptaVelocityDotY += dva_dt.gety(); ptaVelocityDotZ += dva_dt.getz();
		    	    dvb_dt = -SPHParticleCellVec[pvec][pta]->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
			    #pragma omp critical (demtCritical)
			    {
			    demt->addForce(tmp_particleVec[ptb]->getParticleMass()*dvb_dt);
			    demt->addMoment( (ptb_position-demt_curr) * tmp_particleVec[ptb]->getParticleMass()*dvb_dt );
			    }
//		    	    tmp_particleVec[ptb]->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

		            delta_a = epsilon*tmp_particleVec[ptb]->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
			    ptaVelocityCorrectionX += delta_a.getx(); ptaVelocityCorrectionY += delta_a.gety(); ptaVelocityCorrectionZ += delta_a.getz(); 
//		            delta_b = epsilon*SPHParticleCellVec[pvec][pta]->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//		            tmp_particleVec[ptb]->addVelocityCorrection(delta_b);

			    break;
			  case 3:	// ptb is boundary particle

			    if(SPHParticleCellVec[pvec][pta]->getType()!=1){	// pta is not free sph particles, i.e. pta is ghost or boundary particles
		    		break;	// do not consider pta, as we treat before
			    }

	    		    // interact with boundary particles
	    		    da = pta_position.getz();	// assume with the bottom boundary
	    		    if(pta_position.getx() < da){	// with left boundary
				da = pta_position.getx();
	    		    } 
	    		    else if(length_x-pta_position.getx() < da){	// with right boundary
				da = length_x-pta_position.getx();
	    		    }
	    		    else if(height_z-pta_position.getz() < da){	// with top boundary
				da = height_z-pta_position.getz();
	    		    }

		    	    // calculate Vab as the method shown in Morris's paper, 1996

	    	    	    dB = 0-ptb_position.getz();	// assume with the bottom boundary
	    	    	    if(ptb_position.getx() < dB){	// with left boundary
		    		dB = 0-ptb_position.getx();
	    	    	    } 
	    	    	    else if(ptb_position.getx()-length_x < dB){	// with right boundary
		    		dB = ptb_position.getx()-length_x;
	    	    	    }
	    	    	    else if(ptb_position.getz()-height_z < dB){	// with top boundary
		    		dB = ptb_position.getz()-height_z;
	    	    	    }

		    	    beta = 1+dB/da;
		    	    if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

		    	    vab = beta*SPHParticleCellVec[pvec][pta]->getVelocity();
		    	    vba = -vab;

		    	    ptaDensityDot += tmp_particleVec[ptb]->getParticleMass()*(vab%delta_aWab);
		    	    tmp_particleVec[ptb]->addDensityDot( SPHParticleCellVec[pvec][pta]->getParticleMass()*(vba%delta_bWba) );

		    	    vr_dot = vab%(pta_position-ptb_position);
		    	    if(vr_dot<0){
				mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    		Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    }
		    	    else{
		 		Gamma_ab = 0;
		    	    }

		    	    dva_dt = -tmp_particleVec[ptb]->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
			    ptaVelocityDotX += dva_dt.getx(); ptaVelocityDotY += dva_dt.gety(); ptaVelocityDotZ += dva_dt.getz();
//		    	    dvb_dt = -SPHParticleCellVec[pvec][pta]->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
//		    	    tmp_particleVec[ptb]->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

		    	    delta_a = epsilon*tmp_particleVec[ptb]->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
			    ptaVelocityCorrectionX += delta_a.getx(); ptaVelocityCorrectionY += delta_a.gety(); ptaVelocityCorrectionZ += delta_a.getz(); 
// //		    	    delta_b = epsilon*SPHParticleCellVec[pvec][pta]->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//		    	    tmp_particleVec[ptb]->addVelocityCorrection(delta_b);
			
//			    // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
//		    	    if(rab<=space_interval){ // ptb is in the smooth kernel
//		    		dva_dt = D*(pow(space_interval/rab, p1)-pow(space_interval/rab, p2))*(pta_position-ptb_position)/(rab*rab);
//		    		ptaVelocityDotX += dva_dt.getx(); ptaVelocityDotY += dva_dt.gety(); ptaVelocityDotZ += dva_dt.getz();
//		    	    } // end if
			    break;
			  default:
			    std::cout << "SPH particle type should be 1, 2 or 3!" << std::endl;
			    exit(-1);
		
			} //end swtich type
		    } // end if 3h
		} // end for ptb in neighbor cells
} // #pragma omp parallel

		SPHParticleCellVec[pvec][pta]->addDensityDot(ptaDensityDot);
		SPHParticleCellVec[pvec][pta]->addVelocityDot(dem::vec(ptaVelocityDotX,ptaVelocityDotY,ptaVelocityDotZ));
		SPHParticleCellVec[pvec][pta]->addVelocityCorrection(dem::vec(ptaVelocityCorrectionX,ptaVelocityCorrectionY,ptaVelocityCorrectionZ));

	    } // end for pta

	    tmp_particleVec.clear();	// clear elements in tmp-vector for particles neighboring cells, it is important

	} // end for pvec, different cells


//	    // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
//	    for(ptb=SPHBoundaryParticleVec.begin(); ptb!=SPHBoundaryParticleVec.end(); ptb++){	
//		ptb_position = (*ptb)->getCurrPosition();
//		rab = dem::vfabs(pta_position-ptb_position);
//		if(rab<=space_interval){ // ptb is in the smooth kernel
//		    dva_dt = D*(pow(space_interval/rab, p1)-pow(space_interval/rab, p2))*(pta_position-ptb_position)/(rab*rab);
//		    (*pta)->addVelocityDot(dva_dt);
//		} // end if
//	    } // end ptb


    } // end calculateSPHDensityDotVelocityDotLinkedList2D()


// else def OPENMP
#else

    //// the momentum equilibrium equation and state equation are implemented as in Monaghan's paper (1994), simulate free surface flow using sph
    //// here the neighboring list of SPH particles is searched by the cells, vector< vector<sph::SPHParticle*> > SPHParticleCellVec;
    void assembly::calculateSPHDensityDotVelocityDotLinkedList2D(){

	// divide the SPH domain into different cells, each cell will contain SPH particles within it
	divideSPHDomain2D();

//	checkDivision();	// pass, May 22, 2015
//	checkNeighborCells2D();	// pass, May 22, 2015

	// initialize the densityDot and velocityDot of all the SPH particles 
	dem::vec tmp_vec = dem::vec(0,0, -(dem::G)*(dem::GRVT_SCL) );
	dem::vec zero_vec = dem::vec(0,0,0);
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    (*pt)->setDensityDotVelocityDotZero();
	    (*pt)->addVelocityDot(tmp_vec);
	    (*pt)->calculateParticleViscosity();	// everytime when use getParticleViscolisity(), make sure that pressure has been calculated!!!!
	    (*pt)->calculateParticlePressure();		// everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
	}
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	    (*pt)->setDensityDotVelocityDotZero();
	    (*pt)->addVelocityDot(zero_vec);
	    (*pt)->calculateParticleViscosity();	// everytime when use getParticleViscolisity(), make sure that pressure has been calculated!!!!
	    (*pt)->calculateParticlePressure();		// everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
	}
	std::vector<sph::SPHParticle*>::iterator gt;
	for(std::vector<particle*>::iterator pdt=ParticleVec.begin(); pdt!=ParticleVec.end(); pdt++){
	    for(gt=(*pdt)->SPHGhostParticleVec.begin(); gt!=(*pdt)->SPHGhostParticleVec.end(); gt++){
	    	(*gt)->setDensityDotVelocityDotZero();
	    	(*gt)->addVelocityDot(tmp_vec);
	    	(*gt)->calculateParticleViscosity();	// everytime when use getParticleViscolisity(), make sure that pressure has been calculated!!!!
	    	(*gt)->calculateParticlePressure();	// everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
	    }
	}

	// temporary variables used in the loop
	dem::vec pta_position;
	dem::vec ptb_position;
	dem::vec delta_aWab;
	dem::vec delta_bWba;
	dem::vec vab;
	dem::vec vba;
	dem::vec vdem;
	dem::vec dva_dt;
	dem::vec dvb_dt;
	dem::vec delta_a;
	dem::vec delta_b;
	REAL pa, pb, rhoa, rhob, mua, mub;
	REAL rab;
        REAL dWab_dra;
	REAL dWba_drb;
	REAL Wab, Wba;
	REAL da, dB;
	REAL beta;
	REAL xa, ya, xB, yB, k, sqrt_xaya;	// variables for Morris' method to calculate da/dB
	std::vector<sph::SPHParticle*>::iterator ptb;
	REAL ra, rb, rc;	// the geometry of the dem particle
	dem::vec pta_local, ptb_local;	// local position of pta and ptb in the dem particle
	dem::vec demt_curr;
	REAL Gamma_ab, mu_ab, vr_dot;
	REAL alpha = 0.0;
	REAL epsilon = 0.0;	// parameter for velocity correction
	REAL Wq, Ra, Rb, phi_4, coefficient_a, coefficient_b;
	dem::particle* demt;

	std::vector<sph::SPHParticle*> tmp_particleVec;	// sph particles in neighboring cells
	for(int pvec=0; pvec<SPHParticleCellVec.size(); pvec++){

	    // store all SPH particles in pvec's neighboring cells
	    tmp_particleVec.clear();	// it is the same for the particles in the same cell
	    int pnum;
	    if(pvec+1<numCell){
		pnum = pvec+1;
    	    	for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pnum].begin(); pt!= SPHParticleCellVec[pnum].end(); pt++) {
		    tmp_particleVec.push_back(*pt);
    	    	}
	    }
	    if(pvec+Nx-1<numCell){
		pnum = pvec+Nx-1;
    	    	for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pnum].begin(); pt!= SPHParticleCellVec[pnum].end(); pt++) {
		    tmp_particleVec.push_back(*pt);
    	    	}
	    }
	    if(pvec+Nx<numCell){
		pnum = pvec+Nx;
    	    	for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pnum].begin(); pt!= SPHParticleCellVec[pnum].end(); pt++) {
		    tmp_particleVec.push_back(*pt);
    	    	}
	    }
	    if(pvec+Nx+1<numCell){
		pnum = pvec+Nx+1;
    	    	for(std::vector<sph::SPHParticle*>::iterator pt = SPHParticleCellVec[pnum].begin(); pt!= SPHParticleCellVec[pnum].end(); pt++) {
		    tmp_particleVec.push_back(*pt);
    	    	}
	    }

	    for(std::vector<sph::SPHParticle*>::iterator pta=SPHParticleCellVec[pvec].begin(); pta!=SPHParticleCellVec[pvec].end(); pta++){	// SPH particles in cell pvec
		pa = (*pta)->getParticlePressure();
	    	if(pa>=0){
		    Ra = 0.006;
	    	}
	    	else{
		    Ra = 0.6;
	    	}
//	    	mua = (*pta)->getParticleViscosity();
	    	rhoa = (*pta)->getParticleDensity();
	    	pta_position = (*pta)->getCurrPosition();

//		if((*pta)->getType()!=1){	// pta is not free sph particles, i.e. pta is ghost or boundary particles
//		    continue;	// do not consider pta, as we treat before
//		}		

		for(ptb=pta+1; ptb!=SPHParticleCellVec[pvec].end(); ptb++){	// sum over the SPH particles in the same cell pvec
		    ptb_position = (*ptb)->getCurrPosition();
		    rab = dem::vfabs(pta_position-ptb_position);
		    if(rab<=3*dem::smoothLength){	// ptb is in the smooth kernel
		    	pb = (*ptb)->getParticlePressure();
//		    	mub = (*ptb)->getParticleViscosity();
		    	rhob = (*ptb)->getParticleDensity();

		    	Wq = kernelFunction(rab/dem::smoothLength);
		    	phi_4 = pow(Wq/Wqmin,4);
		    	if(pb>=0){
			    Rb = 0.006;
		    	}
		    	else{
			    Rb = 0.6;
		    	}
		    	coefficient_a = 1+Ra*phi_4;
		    	coefficient_b = 1+Rb*phi_4;

			// we have three types of SPH particles: 1, free particle; 2, ghost particle; 3, boundary particle
			// Then we have 3x3 = 9 different types of interactions with the three types of particles
			// so we cannot judge the interaction type only by the type of particle ptb, we need to consider pta also
			// pta          ptb          need to consider or not
			// 1            1              V			free with free
			// 1            2              V			free with ghost
			// 1            3              V			free with boundary

			// 2            1              V			ghost with free
			// 2            2              X			ghost with ghost
			// 2            3              X			ghost with boundary

			// 3            1              V			boundary with free
			// 3            2              X			boundary with ghost
			// 3            3              X 			boundary with boundary


			// add the density dot for pta and ptb
			delta_aWab = gradientKernelFunction(pta_position, ptb_position);	// this is to add SPH pta 
		    	delta_bWba = gradientKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
		    	dWab_dra = partialKernelFunction(pta_position, ptb_position);	// this is to add SPH pta
		    	dWba_drb = partialKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
		    	Wab = kernelFunction(pta_position, ptb_position);	// this is to add SPH pta 
		    	Wba = kernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
			switch((*ptb)->getType()){
			  case 1:	// ptb is free SPH particle
			    switch((*pta)->getType()){
			      case 1:	// free with free
		    	    	vab = (*pta)->getVelocity()-(*ptb)->getVelocity();
		    	   	vba = -vab;

		    	    	(*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    	    	(*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    	    	vr_dot = vab%(pta_position-ptb_position);
		    	    	if(vr_dot<0){
				    mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    		    Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    	}
		    	    	else{
		 		    Gamma_ab = 0;
		    	    	}


		    	    	dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
		    	    	(*pta)->addVelocityDot(dva_dt);
		    	    	dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
		    	    	(*ptb)->addVelocityDot(dvb_dt);

		    	    	delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
		    	    	(*pta)->addVelocityCorrection(delta_a);
		    	    	delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
		    	    	(*ptb)->addVelocityCorrection(delta_b);
			        break;
			      case 2:	// ghost with free
				demt = (*pta)->getDemParticle();
				demt_curr = demt->getCurrPosition();
			    	ptb_local = demt->localVec(ptb_position-demt_curr);	// the local position of sph point ptb
			    	ra = demt->getA(); rb = demt->getB(); rc = demt->getC();
			    	k = 1.0/(sqrt(ptb_local.getx()*ptb_local.getx()/(ra*ra)+ptb_local.gety()*ptb_local.gety()/(rb*rb)+ptb_local.getz()*ptb_local.getz()/(rc*rc) ) );
			    	da = dem::vfabs(ptb_local-k*ptb_local);	// the distance is the same in rotated coordinates

		    	    	// calculate Vab as the method shown in Morris's paper, 1996

		    	    	// (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
			    	pta_local = demt->localVec(pta_position-demt_curr);
			    	k = 1.0/(sqrt(pta_local.getx()*pta_local.getx()/(ra*ra)+pta_local.gety()*pta_local.gety()/(rb*rb)+pta_local.getz()*pta_local.getz()/(rc*rc) ) );
		    	    	dB = dem::vfabs(pta_local-k*pta_local);

//		    	    // (2) here the Morris's method is used 
//		    	    xa = pta_position.getx(); ya = pta_position.gety();
//		    	    xB = ptb_position.getx(); yB = ptb_position.gety();
//		    	    if(ya==0) {da = xa-radius; dB = radius-xB;}
//		    	    else if(xa==0) {da = ya-radius; dB = radius-yB;}
//		    	    else {
//			    sqrt_xaya = sqrt(xa*xa+ya*ya);
//			    k = radius/sqrt_xaya;
//			    da = radius/k - radius;
//			    dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//		    	    }


		    	    	beta = 1+dB/da;
		    	    	if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

			    	vdem = demt->getCurrVelocity()+demt->getCurrOmga()*( pta_position-demt_curr );
		    	    	vba = beta*((*ptb)->getVelocity()-vdem);
		    	    	vab = -vba;

		    	    	(*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    	    	(*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    	    	vr_dot = vab%(pta_position-ptb_position);
		    	    	if(vr_dot<0){
			    	    mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    	    	    Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    	}
		    	    	else{
		 	    	    Gamma_ab = 0;
		    	    	}

		    	    	dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
			    	demt->addForce((*pta)->getParticleMass()*dva_dt);
			    	demt->addMoment( (pta_position-demt_curr) * (*pta)->getParticleMass()*dva_dt );
//		    	    	(*pta)->addVelocityDot(dva_dt);

		    	    	dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
		    	    	(*ptb)->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

//		            	delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
//		            	(*pta)->addVelocityCorrection(delta_a);
		            	delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
		            	(*ptb)->addVelocityCorrection(delta_b);
				break;
			      case 3:	// boundary with free

	    		    	// interact with boundary particles
	    		    	da = ptb_position.getz();	// assume with the bottom boundary
	    		    	if(ptb_position.getx() < da){	// with left boundary
				    da = ptb_position.getx();
	    		    	} 
	    		    	else if(length_x-ptb_position.getx() < da){	// with right boundary
				    da = length_x-ptb_position.getx();
	    		    	}
	    		    	else if(height_z-ptb_position.getz() < da){	// with top boundary
				    da = height_z-ptb_position.getz();
	    		    	}

		    	    	// calculate Vab as the method shown in Morris's paper, 1996

	    	    	    	dB = 0-pta_position.getz();	// assume with the bottom boundary
	    	    	    	if(pta_position.getx() < dB){	// with left boundary
		    		    dB = 0-pta_position.getx();
	    	    	    	} 
	    	    	    	else if(pta_position.getx()-length_x < dB){	// with right boundary
		    		    dB = pta_position.getx()-length_x;
	    	    	    	}
	    	    	    	else if(pta_position.getz()-height_z < dB){	// with top boundary
		    		    dB = pta_position.getz()-height_z;
	    	    	    	}

		    	    	beta = 1+dB/da;
		    	    	if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

		    	    	vba = beta*(*ptb)->getVelocity();
		    	    	vab = -vba;

		    	    	(*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    	    	(*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    	    	vr_dot = vab%(pta_position-ptb_position);
		    	    	if(vr_dot<0){
				    mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    		    Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    	}
		    	    	else{
		 		    Gamma_ab = 0;
		    	    	}

//		    	    	dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
//		    	    	(*pta)->addVelocityDot(dva_dt);
		    	    	dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
		    	    	(*ptb)->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

//			    	delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
//		    	    	(*pta)->addVelocityCorrection(delta_a);
		    	    	delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
		    	    	(*ptb)->addVelocityCorrection(delta_b);
			
//			    	// apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
//		    	    	if(rab<=space_interval){ // ptb is in the smooth kernel
//		    		    dvb_dt = D*(pow(space_interval/rab, p1)-pow(space_interval/rab, p2))*(ptb_position-pta_position)/(rab*rab);
//		    		    (*ptb)->addVelocityDot(dvb_dt);
//		    	    	} // end if

				break;
			      default:
			    	std::cout << "SPH particle type of pta should be 1, 2 or 3!" << std::endl;
			    	exit(-1);
			      } // end switch pta

			    break;
			  case 2:	// ptb is ghost particle

			    if((*pta)->getType()!=1){	// pta is not free sph particles, i.e. pta is ghost or boundary particles
		    		break;	// do not consider pta, as we treat before
			    }
			    demt = (*ptb)->getDemParticle();
			    demt_curr = demt->getCurrPosition();
			    pta_local = demt->localVec(pta_position-demt_curr);	// the local position of sph point pta
			    ra = demt->getA(); rb = demt->getB(); rc = demt->getC();
			    k = 1.0/(sqrt(pta_local.getx()*pta_local.getx()/(ra*ra)+pta_local.gety()*pta_local.gety()/(rb*rb)+pta_local.getz()*pta_local.getz()/(rc*rc) ) );
			    da = dem::vfabs(pta_local-k*pta_local);	// the distance is the same in rotated coordinates

		    	    // calculate Vab as the method shown in Morris's paper, 1996

		    	    // (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
			    ptb_local = demt->localVec(ptb_position-demt_curr);
			    k = 1.0/(sqrt(ptb_local.getx()*ptb_local.getx()/(ra*ra)+ptb_local.gety()*ptb_local.gety()/(rb*rb)+ptb_local.getz()*ptb_local.getz()/(rc*rc) ) );
		    	    dB = dem::vfabs(ptb_local-k*ptb_local);

//		    	    // (2) here the Morris's method is used 
//		    	    xa = pta_position.getx(); ya = pta_position.gety();
//		    	    xB = ptb_position.getx(); yB = ptb_position.gety();
//		    	    if(ya==0) {da = xa-radius; dB = radius-xB;}
//		    	    else if(xa==0) {da = ya-radius; dB = radius-yB;}
//		    	    else {
//			    sqrt_xaya = sqrt(xa*xa+ya*ya);
//			    k = radius/sqrt_xaya;
//			    da = radius/k - radius;
//			    dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//		    	    }


		    	    beta = 1+dB/da;
		    	    if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

			    vdem = demt->getCurrVelocity()+demt->getCurrOmga()*( ptb_position-demt_curr );
		    	    vab = beta*((*pta)->getVelocity()-vdem);
		    	    vba = -vab;

		    	    (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    	    (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    	    vr_dot = vab%(pta_position-ptb_position);
		    	    if(vr_dot<0){
			    	mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    	    	Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    }
		    	    else{
		 	    	Gamma_ab = 0;
		    	    }

		    	    dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
		    	    (*pta)->addVelocityDot(dva_dt);
		    	    dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
			    demt->addForce((*ptb)->getParticleMass()*dvb_dt);
			    demt->addMoment( (ptb_position-demt_curr) * (*ptb)->getParticleMass()*dvb_dt );
//		    	    (*ptb)->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

		            delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
		            (*pta)->addVelocityCorrection(delta_a);
//		            delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//		            (*ptb)->addVelocityCorrection(delta_b);

			    break;
			  case 3:	// ptb is boundary particle

			    if((*pta)->getType()!=1){	// pta is not free sph particles, i.e. pta is ghost or boundary particles
		    		break;	// do not consider pta, as we treat before
			    }

	    		    // interact with boundary particles
	    		    da = pta_position.getz();	// assume with the bottom boundary
	    		    if(pta_position.getx() < da){	// with left boundary
				da = pta_position.getx();
	    		    } 
	    		    else if(length_x-pta_position.getx() < da){	// with right boundary
				da = length_x-pta_position.getx();
	    		    }
	    		    else if(height_z-pta_position.getz() < da){	// with top boundary
				da = height_z-pta_position.getz();
	    		    }

		    	    // calculate Vab as the method shown in Morris's paper, 1996

	    	    	    dB = 0-ptb_position.getz();	// assume with the bottom boundary
	    	    	    if(ptb_position.getx() < dB){	// with left boundary
		    		dB = 0-ptb_position.getx();
	    	    	    } 
	    	    	    else if(ptb_position.getx()-length_x < dB){	// with right boundary
		    		dB = ptb_position.getx()-length_x;
	    	    	    }
	    	    	    else if(ptb_position.getz()-height_z < dB){	// with top boundary
		    		dB = ptb_position.getz()-height_z;
	    	    	    }

		    	    beta = 1+dB/da;
		    	    if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

		    	    vab = beta*(*pta)->getVelocity();
		    	    vba = -vab;

		    	    (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    	    (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    	    vr_dot = vab%(pta_position-ptb_position);
		    	    if(vr_dot<0){
				mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    		Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    }
		    	    else{
		 		Gamma_ab = 0;
		    	    }

		    	    dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
		    	    (*pta)->addVelocityDot(dva_dt);
//		    	    dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
//		    	    (*ptb)->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

		    	    delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
		    	    (*pta)->addVelocityCorrection(delta_a);
//		    	    delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//		    	    (*ptb)->addVelocityCorrection(delta_b);
			
//			    // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
//		    	    if(rab<=space_interval){ // ptb is in the smooth kernel
//		    		dva_dt = D*(pow(space_interval/rab, p1)-pow(space_interval/rab, p2))*(pta_position-ptb_position)/(rab*rab);
//		    		(*pta)->addVelocityDot(dva_dt);
//		    	    } // end if
			    break;
			  default:
			    std::cout << "SPH particle type should be 1, 2 or 3!" << std::endl;
			    exit(-1);
		
			} //end swtich type
		    } // end if 3h
		} // end for ptb in the same cell

		for(ptb=tmp_particleVec.begin(); ptb!=tmp_particleVec.end(); ptb++){	// all particles in pvec's neighboring cells 
		    ptb_position = (*ptb)->getCurrPosition();
		    rab = dem::vfabs(pta_position-ptb_position);
		    if(rab<=3*dem::smoothLength){	// ptb is in the smooth kernel
		    	pb = (*ptb)->getParticlePressure();
//		    	mub = (*ptb)->getParticleViscosity();
		    	rhob = (*ptb)->getParticleDensity();

		    	Wq = kernelFunction(rab/dem::smoothLength);
		    	phi_4 = pow(Wq/Wqmin,4);
		    	if(pb>=0){
			    Rb = 0.006;
		    	}
		    	else{
			    Rb = 0.6;
		    	}
		    	coefficient_a = 1+Ra*phi_4;
		    	coefficient_b = 1+Rb*phi_4;

			// add the density dot for pta and ptb
			delta_aWab = gradientKernelFunction(pta_position, ptb_position);	// this is to add SPH pta 
		    	delta_bWba = gradientKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
		    	dWab_dra = partialKernelFunction(pta_position, ptb_position);	// this is to add SPH pta
		    	dWba_drb = partialKernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
		    	Wab = kernelFunction(pta_position, ptb_position);	// this is to add SPH pta 
		    	Wba = kernelFunction(ptb_position, pta_position);	// this is to add SPH ptb
			switch((*ptb)->getType()){
			  case 1:	// ptb is free SPH particle
			    switch((*pta)->getType()){
			      case 1:	// free with free
		    	    	vab = (*pta)->getVelocity()-(*ptb)->getVelocity();
		    	   	vba = -vab;

		    	    	(*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    	    	(*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    	    	vr_dot = vab%(pta_position-ptb_position);
		    	    	if(vr_dot<0){
				    mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    		    Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    	}
		    	    	else{
		 		    Gamma_ab = 0;
		    	    	}


		    	    	dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
		    	    	(*pta)->addVelocityDot(dva_dt);
		    	    	dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
		    	    	(*ptb)->addVelocityDot(dvb_dt);

		    	    	delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
		    	    	(*pta)->addVelocityCorrection(delta_a);
		    	    	delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
		    	    	(*ptb)->addVelocityCorrection(delta_b);
			        break;
			      case 2:	// ghost with free
				demt = (*pta)->getDemParticle();
				demt_curr = demt->getCurrPosition();
			    	ptb_local = demt->localVec(ptb_position-demt_curr);	// the local position of sph point ptb
			    	ra = demt->getA(); rb = demt->getB(); rc = demt->getC();
			    	k = 1.0/(sqrt(ptb_local.getx()*ptb_local.getx()/(ra*ra)+ptb_local.gety()*ptb_local.gety()/(rb*rb)+ptb_local.getz()*ptb_local.getz()/(rc*rc) ) );
			    	da = dem::vfabs(ptb_local-k*ptb_local);	// the distance is the same in rotated coordinates

		    	    	// calculate Vab as the method shown in Morris's paper, 1996

		    	    	// (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
			    	pta_local = demt->localVec(pta_position-demt_curr);
			    	k = 1.0/(sqrt(pta_local.getx()*pta_local.getx()/(ra*ra)+pta_local.gety()*pta_local.gety()/(rb*rb)+pta_local.getz()*pta_local.getz()/(rc*rc) ) );
		    	    	dB = dem::vfabs(pta_local-k*pta_local);

//		    	    // (2) here the Morris's method is used 
//		    	    xa = pta_position.getx(); ya = pta_position.gety();
//		    	    xB = ptb_position.getx(); yB = ptb_position.gety();
//		    	    if(ya==0) {da = xa-radius; dB = radius-xB;}
//		    	    else if(xa==0) {da = ya-radius; dB = radius-yB;}
//		    	    else {
//			    sqrt_xaya = sqrt(xa*xa+ya*ya);
//			    k = radius/sqrt_xaya;
//			    da = radius/k - radius;
//			    dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//		    	    }


		    	    	beta = 1+dB/da;
		    	    	if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

			    	vdem = demt->getCurrVelocity()+demt->getCurrOmga()*( pta_position-demt_curr );
		    	    	vba = beta*((*ptb)->getVelocity()-vdem);
		    	    	vab = -vba;

		    	    	(*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    	    	(*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    	    	vr_dot = vab%(pta_position-ptb_position);
		    	    	if(vr_dot<0){
			    	    mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    	    	    Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    	}
		    	    	else{
		 	    	    Gamma_ab = 0;
		    	    	}

		    	    	dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
			    	demt->addForce((*pta)->getParticleMass()*dva_dt);
			    	demt->addMoment( (pta_position-demt_curr) * (*pta)->getParticleMass()*dva_dt );
//		    	    	(*pta)->addVelocityDot(dva_dt);

		    	    	dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
		    	    	(*ptb)->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

//		            	delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
//		            	(*pta)->addVelocityCorrection(delta_a);
		            	delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
		            	(*ptb)->addVelocityCorrection(delta_b);
				break;
			      case 3:	// boundary with free

	    		    	// interact with boundary particles
	    		    	da = ptb_position.getz();	// assume with the bottom boundary
	    		    	if(ptb_position.getx() < da){	// with left boundary
				    da = ptb_position.getx();
	    		    	} 
	    		    	else if(length_x-ptb_position.getx() < da){	// with right boundary
				    da = length_x-ptb_position.getx();
	    		    	}
	    		    	else if(height_z-ptb_position.getz() < da){	// with top boundary
				    da = height_z-ptb_position.getz();
	    		    	}

		    	    	// calculate Vab as the method shown in Morris's paper, 1996

	    	    	    	dB = 0-pta_position.getz();	// assume with the bottom boundary
	    	    	    	if(pta_position.getx() < dB){	// with left boundary
		    		    dB = 0-pta_position.getx();
	    	    	    	} 
	    	    	    	else if(pta_position.getx()-length_x < dB){	// with right boundary
		    		    dB = pta_position.getx()-length_x;
	    	    	    	}
	    	    	    	else if(pta_position.getz()-height_z < dB){	// with top boundary
		    		    dB = pta_position.getz()-height_z;
	    	    	    	}

		    	    	beta = 1+dB/da;
		    	    	if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

		    	    	vba = beta*(*ptb)->getVelocity();
		    	    	vab = -vba;

		    	    	(*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    	    	(*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    	    	vr_dot = vab%(pta_position-ptb_position);
		    	    	if(vr_dot<0){
				    mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    		    Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    	}
		    	    	else{
		 		    Gamma_ab = 0;
		    	    	}

//		    	    	dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
//		    	    	(*pta)->addVelocityDot(dva_dt);
		    	    	dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
		    	    	(*ptb)->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

//			    	delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
//		    	    	(*pta)->addVelocityCorrection(delta_a);
		    	    	delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
		    	    	(*ptb)->addVelocityCorrection(delta_b);
			
			    	// apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
		    	    	if(rab<=space_interval){ // ptb is in the smooth kernel
		    		    dvb_dt = D*(pow(space_interval/rab, p1)-pow(space_interval/rab, p2))*(ptb_position-pta_position)/(rab*rab);
		    		    (*ptb)->addVelocityDot(dvb_dt);
		    	    	} // end if

				break;
			      default:
			    	std::cout << "SPH particle type of pta should be 1, 2 or 3!" << std::endl;
			    	exit(-1);
			      } // end switch pta

			    break;
			  case 2:	// ptb is ghost particle

			    if((*pta)->getType()!=1){	// pta is not free sph particles, i.e. pta is ghost or boundary particles
		    		break;	// do not consider pta, as we treat before
			    }

			    demt = (*ptb)->getDemParticle();
			    demt_curr = demt->getCurrPosition();
			    pta_local = demt->localVec(pta_position-demt_curr);	// the local position of sph point pta
			    ra = demt->getA(); rb = demt->getB(); rc = demt->getC();
			    k = 1.0/(sqrt(pta_local.getx()*pta_local.getx()/(ra*ra)+pta_local.gety()*pta_local.gety()/(rb*rb)+pta_local.getz()*pta_local.getz()/(rc*rc) ) );
			    da = dem::vfabs(pta_local-k*pta_local);	// the distance is the same in rotated coordinates

		    	    // calculate Vab as the method shown in Morris's paper, 1996

		    	    // (1) here I wanna use the distance from the point a/b to the surface of the ellipsoid to simplify the problem
			    ptb_local = demt->localVec(ptb_position-demt_curr);
			    k = 1.0/(sqrt(ptb_local.getx()*ptb_local.getx()/(ra*ra)+ptb_local.gety()*ptb_local.gety()/(rb*rb)+ptb_local.getz()*ptb_local.getz()/(rc*rc) ) );
		    	    dB = dem::vfabs(ptb_local-k*ptb_local);

//		    	    // (2) here the Morris's method is used 
//		    	    xa = pta_position.getx(); ya = pta_position.gety();
//		    	    xB = ptb_position.getx(); yB = ptb_position.gety();
//		    	    if(ya==0) {da = xa-radius; dB = radius-xB;}
//		    	    else if(xa==0) {da = ya-radius; dB = radius-yB;}
//		    	    else {
//			    sqrt_xaya = sqrt(xa*xa+ya*ya);
//			    k = radius/sqrt_xaya;
//			    da = radius/k - radius;
//			    dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
//		    	    }


		    	    beta = 1+dB/da;
		    	    if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

			    vdem = demt->getCurrVelocity()+demt->getCurrOmga()*( ptb_position-demt_curr );
		    	    vab = beta*((*pta)->getVelocity()-vdem);
		    	    vba = -vab;

		    	    (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    	    (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    	    vr_dot = vab%(pta_position-ptb_position);
		    	    if(vr_dot<0){
			    	mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    	    	Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    }
		    	    else{
		 	    	Gamma_ab = 0;
		    	    }

		    	    dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
		    	    (*pta)->addVelocityDot(dva_dt);
		    	    dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)*coefficient_b+pa/(rhoa*rhoa)*coefficient_a+Gamma_ab)*delta_bWba;
			    demt->addForce((*ptb)->getParticleMass()*dvb_dt);
			    demt->addMoment( (ptb_position-demt_curr) * (*ptb)->getParticleMass()*dvb_dt );
//		    	    (*ptb)->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

		            delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
		            (*pta)->addVelocityCorrection(delta_a);
//		            delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//		            (*ptb)->addVelocityCorrection(delta_b);

			    break;
			  case 3:	// ptb is boundary particle

			    if((*pta)->getType()!=1){	// pta is not free sph particles, i.e. pta is ghost or boundary particles
		    		break;	// do not consider pta, as we treat before
			    }

	    		    // interact with boundary particles
	    		    da = pta_position.getz();	// assume with the bottom boundary
	    		    if(pta_position.getx() < da){	// with left boundary
				da = pta_position.getx();
	    		    } 
	    		    else if(length_x-pta_position.getx() < da){	// with right boundary
				da = length_x-pta_position.getx();
	    		    }
	    		    else if(height_z-pta_position.getz() < da){	// with top boundary
				da = height_z-pta_position.getz();
	    		    }

		    	    // calculate Vab as the method shown in Morris's paper, 1996

	    	    	    dB = 0-ptb_position.getz();	// assume with the bottom boundary
	    	    	    if(ptb_position.getx() < dB){	// with left boundary
		    		dB = 0-ptb_position.getx();
	    	    	    } 
	    	    	    else if(ptb_position.getx()-length_x < dB){	// with right boundary
		    		dB = ptb_position.getx()-length_x;
	    	    	    }
	    	    	    else if(ptb_position.getz()-height_z < dB){	// with top boundary
		    		dB = ptb_position.getz()-height_z;
	    	    	    }

		    	    beta = 1+dB/da;
		    	    if(beta>1.5 || beta<0 || isnan(beta)){ beta = 1.5; }

		    	    vab = beta*(*pta)->getVelocity();
		    	    vba = -vab;

		    	    (*pta)->addDensityDot( (*ptb)->getParticleMass()*(vab%delta_aWab) );
		    	    (*ptb)->addDensityDot( (*pta)->getParticleMass()*(vba%delta_bWba) );

		    	    vr_dot = vab%(pta_position-ptb_position);
		    	    if(vr_dot<0){
				mu_ab = dem::smoothLength*vr_dot/(rab*rab+0.01*dem::smoothLength*dem::smoothLength);
		    		Gamma_ab = (-alpha*dem::soundSpeed*mu_ab)/(rhoa+rhob)*2;
		    	    }
		    	    else{
		 		Gamma_ab = 0;
		    	    }

		    	    dva_dt = -(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
		    	    (*pta)->addVelocityDot(dva_dt);
//		    	    dvb_dt = -(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
//		    	    (*ptb)->addVelocityDot(dvb_dt);	// the velocities of the ghost particles will not envolve as the same way as others

		    	    delta_a = epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
		    	    (*pta)->addVelocityCorrection(delta_a);
//		    	    delta_b = epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
//		    	    (*ptb)->addVelocityCorrection(delta_b);

		    	    // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
		    	    if(rab<=space_interval){ // ptb is in the smooth kernel
		    		dva_dt = D*(pow(space_interval/rab, p1)-pow(space_interval/rab, p2))*(pta_position-ptb_position)/(rab*rab);
		    		(*pta)->addVelocityDot(dva_dt);
		    	    } // end if
			
			    break;
			  default:
			    std::cout << "SPH particle type should be 1, 2 or 3!" << std::endl;
			    exit(-1);
		
			} //end swtich type
		    } // end if 3h
		} // end for ptb in neighbor cells
	    } // end for pta

	    tmp_particleVec.clear();	// clear elements in tmp-vector for particles neighboring cells, it is important

	} // end for pvec, different cells


//	    // apply the boundary forces by Lennard-Jones potential as in Monaghan's paper(1994)
//	    for(ptb=SPHBoundaryParticleVec.begin(); ptb!=SPHBoundaryParticleVec.end(); ptb++){	
//		ptb_position = (*ptb)->getCurrPosition();
//		rab = dem::vfabs(pta_position-ptb_position);
//		if(rab<=space_interval){ // ptb is in the smooth kernel
//		    dva_dt = D*(pow(space_interval/rab, p1)-pow(space_interval/rab, p2))*(pta_position-ptb_position)/(rab*rab);
//		    (*pta)->addVelocityDot(dva_dt);
//		} // end if
//	    } // end ptb


    } // end calculateSPHDensityDotVelocityDotLinkedList2D()

// end def OPENMP
#endif



    void assembly::updateSPH(){

#ifdef OPENMP

#pragma omp parallel num_threads(NUM_THREADS)
{

	int thread_num = omp_get_thread_num();



	int start_num = BLOCK_LOW(thread_num, NUM_THREADS, numFreeSPH);

	int end_num   = BLOCK_HIGH(thread_num, NUM_THREADS, numFreeSPH);
	for(int i=start_num; i<end_num+1; i++){
	    SPHParticleVec[i]->updateParticle();
	}

	start_num = BLOCK_LOW(thread_num, NUM_THREADS, numBoundarySPH);

	end_num   = BLOCK_HIGH(thread_num, NUM_THREADS, numBoundarySPH);
	for(int i=start_num; i<end_num+1; i++){

	    SPHBoundaryParticleVec[i]->updateParticle();

	}

	start_num = BLOCK_LOW(thread_num, NUM_THREADS, TotalNum);

	end_num   = BLOCK_HIGH(thread_num, NUM_THREADS, TotalNum);

	for(int i=start_num; i<end_num+1; i++){

	    ParticleVec[i]->updateSPHGhostParticle();

	}
}
//else of def OPENMP
#else
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    (*pt)->updateParticle();
	}
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	    (*pt)->updateParticle();
	}
	std::vector<sph::SPHParticle*>::iterator gt;
	for(std::vector<particle*>::iterator demt=ParticleVec.begin(); demt!=ParticleVec.end(); demt++){
	    (*demt)->updateSPHGhostParticle();	// update position by dem particle and sph particle density
	}
// end def OPENMP
#endif

    } // end updateSPH

    REAL assembly::calculateSPHKineticEnergy(){
	REAL tmp_K = 0;
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    tmp_K = tmp_K+(*pt)->getKineticEnergy();
	}
	return tmp_K;
    } // end calculateSPHKineticEnergy

    REAL assembly::getWaterHeight(){
	REAL tmp_H = 0;
	REAL tmp_z;
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    tmp_z = (*pt)->getCurrPosition().getz();
	    if(tmp_H<tmp_z){
		tmp_H = tmp_z;
	    }
	}
	return tmp_H;
    } // end getWaterHeight

    REAL assembly::getWaterWidth(){
	REAL tmp_W = 0;
	REAL tmp_x;
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    tmp_x = (*pt)->getCurrPosition().getx();
	    if(tmp_W<tmp_x){
		tmp_W = tmp_x;
	    }
	}
	return tmp_W;
    } // end getWaterHeight

    void assembly::initialSPHVelocity(){
	// initialization SPH velocity based on equation(4.3) in http://www.artcompsci.org/vol_1/v1_web/node34.html

	calculateSPHDensityDotVelocityDotLinkedList2D();	// calculate velocityDot and densityDot
	// fix the y for all SPH points
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    (*pt)->fixY();
	}
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	    (*pt)->fixXYZ();
	}
	initialSPHLeapFrogVelocity();	// initial velocity only for free SPH particles based on equation (4.3)

    } // initialSPHVelocity()

    void assembly::updateSPHLeapFrogPositionDensity(){	// update particle position and density based on equation (4.1)

#ifdef OPENMP



#pragma omp parallel num_threads(NUM_THREADS)

{

	int thread_num = omp_get_thread_num();



	int start_num = BLOCK_LOW(thread_num, NUM_THREADS, numFreeSPH);

	int end_num   = BLOCK_HIGH(thread_num, NUM_THREADS, numFreeSPH);

	for(int i=start_num; i<end_num+1; i++){

	    SPHParticleVec[i]->updateParticlePositionDensityLeapFrog();

	}



	start_num = BLOCK_LOW(thread_num, NUM_THREADS, numBoundarySPH);

	end_num   = BLOCK_HIGH(thread_num, NUM_THREADS, numBoundarySPH);

	for(int i=start_num; i<end_num+1; i++){

	    SPHBoundaryParticleVec[i]->updateParticleDensity();

	}



	start_num = BLOCK_LOW(thread_num, NUM_THREADS, TotalNum);

	end_num   = BLOCK_HIGH(thread_num, NUM_THREADS, TotalNum);

	for(int i=start_num; i<end_num+1; i++){

	    ParticleVec[i]->updateSPHGhostParticle();

	}

}

// else of def OPENMP

#else

	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){

	    (*pt)->updateParticlePositionDensityLeapFrog();

	}

	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){

	    (*pt)->updateParticleDensity();	// only update density for boundary particles

	}

	std::vector<sph::SPHParticle*>::iterator gt;

	for(std::vector<particle*>::iterator demt=ParticleVec.begin(); demt!=ParticleVec.end(); demt++){

	    (*demt)->updateSPHGhostParticle();	// update position and density by dem particle and sph particle density

	}



// end of def OPENMP

#endif


    } // end updateSPHLeapFrogPositionDensity

    void assembly::updateSPHLeapFrogVelocity(){	// update particle velocity only for free SPH particles based on equation (4.2)

#ifdef OPENMP



#pragma omp parallel num_threads(NUM_THREADS)

{

	int thread_num = omp_get_thread_num();



	int start_num = BLOCK_LOW(thread_num, NUM_THREADS, numFreeSPH);

	int end_num   = BLOCK_HIGH(thread_num, NUM_THREADS, numFreeSPH);

	for(int i=start_num; i<end_num+1; i++){

	    SPHParticleVec[i]->updateParticleVelocityLeapFrog();

	}



}

// else of def OPENMP

#else

	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){

	    (*pt)->updateParticleVelocityLeapFrog();

	}



//end of def OPENMP

#endif

    } // end updateSPHLeapFrogVelocity

    void assembly::initialSPHLeapFrogVelocity(){	// initial particle velocity only for free SPH particles based on equation (4.3)
#ifdef OPENMP



#pragma omp parallel num_threads(NUM_THREADS)

{

	int thread_num = omp_get_thread_num();



	int start_num = BLOCK_LOW(thread_num, NUM_THREADS, numFreeSPH);

	int end_num   = BLOCK_HIGH(thread_num, NUM_THREADS, numFreeSPH);

	for(int i=start_num; i<end_num+1; i++){

	    SPHParticleVec[i]->initialParticleVelocityLeapFrog();

	}



}

// else of def OPENMP

#else

	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){

	    (*pt)->initialParticleVelocityLeapFrog();

	}

// end of def OPENMP

#endif
    } // end initialSPHLeapFrogVelocity

//	REAL factor = 1.0/(120.0*dem::PI*h*h*h);	// 3D quintic kernel factor
//	REAL factor = 7.0/(478.0*dem::PI*h*h);		// 2D quintic kernel factor
    // kernel function, vec is the position of a b, h is smoothing length	
    inline REAL assembly::kernelFunction(const dem::vec& a, const dem::vec& b){
	REAL rab = dem::vfabs(a-b);
	REAL s = rab*one_devide_h;
	REAL item1_5 = (3-s)*(3-s)*(3-s)*(3-s)*(3-s);	// (3-s)^5
	REAL item2_5 = (2-s)*(2-s)*(2-s)*(2-s)*(2-s); 	// (2-s)^5
	REAL item3_5 = (1-s)*(1-s)*(1-s)*(1-s)*(1-s); 	// (1-s)^5
	if(s<1){
	    return factor_kernel*(item1_5-6*item2_5+15*item3_5);
	}
	else if(s<2){
	    return factor_kernel*(item1_5-6*item2_5);
	}
	else if(s<3){
	    return factor_kernel*item1_5;
	}
	else{
	    return 0.0;
	}

    } // end kernelFunction

    // kernel function, s is rab/h, h is smoothing length	
    inline REAL assembly::kernelFunction(REAL s){
	REAL item1_5 = (3-s)*(3-s)*(3-s)*(3-s)*(3-s);	// (3-s)^5
	REAL item2_5 = (2-s)*(2-s)*(2-s)*(2-s)*(2-s); 	// (2-s)^5
	REAL item3_5 = (1-s)*(1-s)*(1-s)*(1-s)*(1-s); 	// (1-s)^5
	if(s<1){
	    return factor_kernel*(item1_5-6*item2_5+15*item3_5);
	}
	else if(s<2){
	    return factor_kernel*(item1_5-6*item2_5);
	}
	else if(s<3){
	    return factor_kernel*item1_5;
	}
	else{
	    return 0.0;
	}

    } // end kernelFunction

//	REAL factor = 1.0/(120.0*dem::PI*h*h*h*h);	// 3D quintic kernel factor
//	REAL factor = 1.0/(120.0*h*h);	// 1D quintic kernel factor
//	REAL factor = 7.0/(478.0*dem::PI*h*h*h);	// 2D quintic kernel factor
    // to calculate delta_aWab, where a is the position of the first particle
    inline dem::vec assembly::gradientKernelFunction(const dem::vec& a, const dem::vec& b){
	REAL rab = dem::vfabs(a-b);
	REAL s = rab*one_devide_h;
	REAL item1_4 = (3-s)*(3-s)*(3-s)*(3-s);	// (3-s)^4
	REAL item2_4 = (2-s)*(2-s)*(2-s)*(2-s); // (2-s)^4
	REAL item3_4 = (1-s)*(1-s)*(1-s)*(1-s); // (1-s)^4
	if(s<1){
	    return factor_kernel_gradient *(-5*item1_4+30*item2_4-75*item3_4)*(a-b)/rab;	// it is strange that s is devided here. compare to Liu's SPH code (National University of Singapore) 
	}
	else if(s<2){
	    return factor_kernel_gradient *(-5*item1_4+30*item2_4)*(a-b)/rab;
	}
	else if(s<3){
	    return factor_kernel_gradient *(-5*item1_4)*(a-b)/rab;
	}
	else{
	    return dem::vec(0.0);
	}

    } // end gradientKernelFunction

    // to calculate partial differential dWab_dra
    inline REAL assembly::partialKernelFunction(const dem::vec& a, const dem::vec& b){
	REAL rab = dem::vfabs(a-b);
	REAL s = rab*one_devide_h;
	REAL item1_4 = (3-s)*(3-s)*(3-s)*(3-s);	// (3-s)^4
	REAL item2_4 = (2-s)*(2-s)*(2-s)*(2-s); // (2-s)^4
	REAL item3_4 = (1-s)*(1-s)*(1-s)*(1-s); // (1-s)^4
	if(s<1){
	    return factor_kernel_gradient *(-5*item1_4+30*item2_4-75*item3_4);
	}
	else if(s<2){
	    return factor_kernel_gradient *(-5*item1_4+30*item2_4);
	}
	else if(s<3){
	    return factor_kernel_gradient *(-5*item1_4);
	}
	else{
	    return 0.0;
	}

    } // end partialKernelFunction




    void assembly::Cylinder2DProblem(int nsteps,              // total_steps
	       			     int snapshots,           // number of snapshots   
	          		     int interval,
				     const char* iniptclfile,   
		       		     const char* particlefile, 
		       		     const char* contactfile,
		       		     const char* progressfile){

   	// pre_1: open streams for output.
    	// particlefile and contactfile are used for snapshots at the end.
    	progressinf.open(progressfile); 
    	if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
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
	            << setw(OWID) << "average" << endl
	            << setw(OWID) << "number"
	            << setw(OWID) << "contacts"
	            << setw(OWID) << "contacts"
	            << setw(OWID) << "penetration"
	            << setw(OWID) << "contact_normal"
	            << setw(OWID) << "contact_tangt"
	            << setw(OWID) << "velocity"
	            << setw(OWID) << "omga"
	            << setw(OWID) << "force"
	            << setw(OWID) << "moment" << endl;

    	// pre_2. create particles and boundaries from existing files.
    	readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 

	REAL L = 0.2;	// unit m
	REAL SPHInitialDensity = 1000;	// unit kg/m^3
//	REAL a0 = 1.5e-7;		// unit m/s^2, the constant acceleration applied to the particles
	REAL a0 = 0;		// unit m/s^2, the constant acceleration applied to the particles
//	dem::nu = 1e-6;	// unit m^2/s, dynamic viscosity
	dem::nu = 1e-4;	// unit m^2/s, dynamic viscosity
	nSPHPoint = 50;	// number of points in each direction, x and y

	space_interval = L/(nSPHPoint);
	dem::smoothLength = 1.5*(space_interval);	// refer to Morris's paper Sect.2.9
	// calculate paramters used in kernel funciton
	one_devide_h = 1.0/dem::smoothLength;
	factor_kernel = 7.0/(478.0*dem::PI*dem::smoothLength*dem::smoothLength);
	factor_kernel_gradient = 7.0/(478.0*dem::PI*dem::smoothLength*dem::smoothLength*dem::smoothLength);

//	dem::soundSpeed = 5.77e-4;
	dem::soundSpeed = 1e-2;
        dem::vec initVelocity = dem::vec(-1.0e-4, 0, 0);	// the initial velocity for dem particle
	dem::vec initOmega = dem::vec(0,0,-2e-2);

	// generate the sph particle and store them in SPHParticleVec
	REAL SPHmass = SPHInitialDensity*L*L*1.0/(nSPHPoint*nSPHPoint);

	// free the spaces of SPH points
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    delete (*pt);
	}
	SPHParticleVec.clear();
		// free the spaces of ghost points in each dem particle
	for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
	    for(std::vector<sph::SPHParticle*>::iterator spht=(*pt)->SPHGhostParticleVec.begin(); spht!=(*pt)->SPHGhostParticleVec.end(); spht++){
		delete (*spht);
	   }
	   (*pt)->SPHGhostParticleVec.clear();
	}

	int isGhost = 0;	// is Ghost particle
	REAL radius_a, radius_b;
  	dem::vec pt_position;
	REAL pt_x, pt_y;	// the position of the DEM particle
	dem::vec local_x;	// local position of ghost point in dem particle
	REAL tmp_x, tmp_y;
	dem::vec tmp_xyz, tmp_local;
	// create and store PeriParticle objects into periParticleVec
	for(int num_x=0; num_x!=nSPHPoint; num_x++){
	    tmp_x = num_x*space_interval-0.5*L+0.5*space_interval;
	    for(int num_y=0; num_y!=nSPHPoint; num_y++){
		tmp_y = num_y*space_interval-0.5*L+0.5*space_interval;
		isGhost = 0;
		tmp_xyz = dem::vec(tmp_x, tmp_y, 0);
	    	for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
		    radius_a = (*pt)->getA(); radius_b = (*pt)->getB();
		    pt_position = (*pt)->getCurrPosition();
		    tmp_local = (*pt)->localVec(tmp_xyz-pt_position);
		    if( tmp_local.getx()*tmp_local.getx()/(radius_a*radius_a)+tmp_local.gety()*tmp_local.gety()/(radius_b*radius_b) <= 1 ){
		    	isGhost = 1;	// is Ghost particle
		    	local_x = (*pt)->localVec( dem::vec(tmp_x, tmp_y, 0)-pt_position );
		    	sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, tmp_y, 0, local_x, (*pt), 2);	// 2 is ghost SPH particle
		    	(*pt)->SPHGhostParticleVec.push_back(tmp_pt);
		    	break;	// if this sph point is ghost for dem particle 1, then it cannot be ghost for any others,
				// should avoid the initial overlap of the different dem particles
		    }
	    	} // end dem particle
	    
	    	if(isGhost==0){	// not Ghost particle
	    	    sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, tmp_y, 0, 1);	// 1 is free SHP particle
	    	    SPHParticleVec.push_back(tmp_pt);
	    	}
	    }
	}

	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    (*pt)->initial();
	}

	int iframe = 0;

    	REAL avgNormal=0;
    	REAL avgTangt=0;
    	int  stepsnum=0;
    	char stepsstr[4];
    	char stepsfp[50];
    	REAL void_ratio=0;
    	REAL bdry_penetr[7];
//    	matrix granularStress;	// granular stress
    	int  bdry_cntnum[7];

//for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
//    (*pt)->setType(1);	// fix the particle
//}

     	dem::vec acceleration = dem::vec(a0,0,0);	// drive the particles in x direction
	dem::vec zeroVec = dem::vec(0,0,0);
	std::ofstream ofs("total_results.dat");
	for(int istep = 0; istep < nsteps; istep++){

//if(istep==nsteps/2){
//for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
//    (*pt)->setType(0);	// free the particle
//}
//}
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
	    rigidBoundaryForce(bdry_penetr, bdry_cntnum);

	    calculateSPHDensityDotVelocityDotPeriodic(L);
	    // add constant acceleration for all SPH points
	    for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
		(*pt)->addVelocityDot(acceleration);
	    }

	    // fix the z for all SPH points
	    for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	      	(*pt)->fixZ();
	    }
	    // should fix the dem particle in z, i.e. force_k = 0; moment_i=moment_j = 0
	    for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
		(*pt)->fixZ();
	    }



//if(istep<0.2*nsteps){
for(std::vector<particle*>::iterator demt=ParticleVec.begin(); demt!=ParticleVec.end(); demt++){
//	   (*demt)->setPrevVelocity(initVelocity);
if(demt==ParticleVec.begin()){
	   (*demt)->setPrevOmga(initOmega);
}
else{
	   (*demt)->setPrevOmga(initOmega);
}
	   (*demt)->setMoment(zeroVec);
}
//}


	    // 5. update particles' velocity/omga/position/orientation based on force/moment.
	    updateParticle();
	    updateSPH();

	    // apply periodic boundary conditions
	    // http://en.wikipedia.org/wiki/Periodic_boundary_conditions
	    // First, an object which has passed through one face of the simulation box should re-enter the box through the opposite face
	    REAL tmp_x, tmp_y;
	    for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
		tmp_x = (*pt)->getCurrPosition().getx();
		tmp_y = (*pt)->getCurrPosition().gety();
		if(tmp_x < -L*0.5) (*pt)->setCurrPositionX(tmp_x+L);	// only valid for box centered at the origin
		if(tmp_x >= L*0.5) (*pt)->setCurrPositionX(tmp_x-L);

		if(tmp_y < -L*0.5) (*pt)->setCurrPositionY(tmp_y+L);
		if(tmp_y >= L*0.5) (*pt)->setCurrPositionY(tmp_y-L);
	    }
	    // should also periodic the dem particle	
	    for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
		tmp_x = (*pt)->getCurrPosition().getx();
		tmp_y = (*pt)->getCurrPosition().gety();
		if(tmp_x < -L*0.5) (*pt)->setCurrPositionX(tmp_x+L);	// only valid for box centered at the origin
		if(tmp_x >= L*0.5) (*pt)->setCurrPositionX(tmp_x-L);

		if(tmp_y < -L*0.5) (*pt)->setCurrPositionY(tmp_y+L);
		if(tmp_y >= L*0.5) (*pt)->setCurrPositionY(tmp_y-L);
	    }

	    if(istep % interval == 0){
		sprintf(stepsstr, "%03d", istep); 
	    	strcpy(stepsfp,"velocity"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
		printSPHVelocityPressure(stepsfp);	

	    	progressinf << setw(OWID) << istep
		            << setw(OWID) << getPossCntctNum()
		            << setw(OWID) << getActualCntctNum()
		            << setw(OWID) << getAveragePenetration()
		            << setw(OWID) << avgNormal
		            << setw(OWID) << avgTangt
		            << setw(OWID) << getAverageVelocity() 
		            << setw(OWID) << getAverageOmga()
		            << setw(OWID) << getAverageForce()   
		            << setw(OWID) << getAverageMoment()
			    << std::endl;
	    } 

	    // output SPH particles information
	    if (istep % (nsteps/snapshots) == 0){
	    	writeSPHParticleTecplot(ofs, iframe);

	    	sprintf(stepsstr, "%03d", iframe); 
	    	strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    	printParticle(stepsfp);
	    
	    	sprintf(stepsstr, "%03d", iframe); 
	    	strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    	printContact(stepsfp);
	   	++iframe;
	    }
	}

	ofs.close();

	sprintf(stepsstr, "%03d", nsteps); 
	strcpy(stepsfp,"velocity"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	printSPHVelocityPressure(stepsfp);

    } // end Cylinder2DProblem


    void assembly::BurstingDam2D(int nsteps,              // total_steps
	       			     int snapshots,           // number of snapshots   
	          		     int interval,
				     const char* iniptclfile,   
				     const char* inibdryfile,
		       		     const char* particlefile, 
		       		     const char* contactfile,
		       		     const char* progressfile){

   	// pre_1: open streams for output.
    	// particlefile and contactfile are used for snapshots at the end.
    	progressinf.open(progressfile); 
    	if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
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
		    << setw(OWID) << "current"
		    << setw(OWID) << "current"
		    << setw(OWID) << "current" << endl
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
		    << setw(OWID) << "time"
		    << setw(OWID) << "height"
		    << setw(OWID) << "width" << endl;

    	// pre_2. create particles and boundaries from existing files.
    	readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    	readBoundary(inibdryfile);   // create boundaries.

	REAL L = 25;	// unit m
	length_x = 100;
	height_z = 50;
  	D = 5*(dem::G)*(dem::GRVT_SCL)*L;
	dem::SPHInitialDensity = 1000;	// unit kg/m^3
	dem::nu = 1e-6;	// unit m^2/s, dynamic viscosity
	nSPHPoint = 101;	// number of points in each direction, x and z
	int num_layers = 5;	// number of layers of boundary sph particles

	space_interval = L/(nSPHPoint-1);
	dem::smoothLength = 1.5*(space_interval);	// refer to Morris's paper Sect.2.9
	// calculate paramters used in kernel funciton
	one_devide_h = 1.0/dem::smoothLength;
	factor_kernel = 7.0/(478.0*dem::PI*dem::smoothLength*dem::smoothLength);
	factor_kernel_gradient = 7.0/(478.0*dem::PI*dem::smoothLength*dem::smoothLength*dem::smoothLength);

//	dem::soundSpeed = 1e-2;
  	dem::gamma = 7;
	dem::soundSpeed = sqrt(200*dem::G*dem::GRVT_SCL*L);
	dem::P0 = dem::soundSpeed*dem::soundSpeed*dem::SPHInitialDensity/dem::gamma;
	dem::P0 = 1.013e5;	// air pressure
	REAL qmin = 0.7593;	// this is fixed for quintic kernel function
	Wqmin = kernelFunction(qmin);
	p1 = 4;	// for the Lennard-Jones boundary forces 
	p2 = 2;
//        dem::vec initVelocity = dem::vec(-1.0e-4, 0, 0);	// the initial velocity for dem particle
//	dem::vec initOmega = dem::vec(0,0,-2e-2);

	// generate the sph particle and store them in SPHParticleVec
	REAL SPHmass = SPHInitialDensity*L*L*1.0/(nSPHPoint*nSPHPoint);

	// free the spaces of SPH points
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    delete (*pt);
	}
	SPHParticleVec.clear();
		// free the spaces of ghost points in each dem particle
	for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
	    for(std::vector<sph::SPHParticle*>::iterator spht=(*pt)->SPHGhostParticleVec.begin(); spht!=(*pt)->SPHGhostParticleVec.end(); spht++){
		delete (*spht);
	   }
	   (*pt)->SPHGhostParticleVec.clear();
	}
	// free the spaces of boundary SPH points
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	    delete (*pt);
	}
	SPHBoundaryParticleVec.clear();

	// create and store PeriParticle objects into periParticleVec
	int isGhost = 0;	// is Ghost particle
	REAL radius_a, radius_b, radius_c;
  	dem::vec pt_position;
	dem::vec local_x;	// local position of ghost point in dem particle
	dem::vec tmp_xyz, tmp_local;
	REAL small_value = 0.01*space_interval;
	for(REAL tmp_z=0-space_interval*(num_layers); tmp_z<height_z+space_interval*num_layers; tmp_z=tmp_z+space_interval){
	    for(REAL tmp_x=0-space_interval*(num_layers); tmp_x<length_x+space_interval*num_layers; tmp_x=tmp_x+space_interval){
		isGhost = 0;
		tmp_xyz = dem::vec(tmp_x, 0, tmp_z);
	    	for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
		    radius_a = (*pt)->getA(); radius_b = (*pt)->getB(); radius_c = (*pt)->getC();
		    pt_position = (*pt)->getCurrPosition();
		    pt_position.sety(0);
		    tmp_local = (*pt)->localVec(tmp_xyz-pt_position);
		    if( tmp_local.getx()*tmp_local.getx()/(radius_a*radius_a)+tmp_local.gety()*tmp_local.gety()/(radius_b*radius_b)+tmp_local.getz()*tmp_local.getz()/(radius_c*radius_c) <= 1 ){
		    	isGhost = 1;	// is Ghost particle
		    	local_x = (*pt)->localVec( dem::vec(tmp_x, 0, tmp_z)-pt_position );
		    	sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, local_x, (*pt), 2);	// 2 is ghost SPH particle
		    	(*pt)->SPHGhostParticleVec.push_back(tmp_pt);
		    	break;	// if this sph point is ghost for dem particle 1, then it cannot be ghost for any others,
				// should avoid the initial overlap of the different dem particles
		    }
	    	} // end dem particle

		if(isGhost==0){	// free/boundary particles
		    if(tmp_x<=-space_interval+small_value || tmp_x>=length_x-small_value || tmp_z<=-space_interval+small_value || tmp_z>=height_z-small_value){	// boundary sph particles
			if(tmp_z>=L+space_interval){
		    	    sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, 3);	// 3 is boundary SPH particle
			    SPHBoundaryParticleVec.push_back(tmp_pt);
			}
			else{
		    	    sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity*pow(1+SPHInitialDensity*(dem::G)*(dem::GRVT_SCL)*(L+space_interval-tmp_z)/(dem::P0), 1.0/(dem::gamma)), tmp_x, 0, tmp_z, 3);	// 3 is boundary SPH particle
			    SPHBoundaryParticleVec.push_back(tmp_pt);
			}
		    } 
		    else if(tmp_x<=L && tmp_z<=L){	// free sph particles
		    	sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity*pow(1+SPHInitialDensity*(dem::G)*(dem::GRVT_SCL)*(L+space_interval-tmp_z)/(dem::P0), 1.0/(dem::gamma)), tmp_x, 0, tmp_z, 1);	// 1 is free SPH particle
//			sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, 1);	// 1 is free SPH particle
	    	    	SPHParticleVec.push_back(tmp_pt);
		    }
		}
	    }
	}

	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    (*pt)->initial();
	}
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	    (*pt)->initial();
	}

#ifdef OPENMP
 	numFreeSPH = SPHParticleVec.size();
	numBoundarySPH = SPHBoundaryParticleVec.size();
#endif

	// generate the vector of cells for SPH particles in order to use linked-list search
	SPHxmax = length_x+space_interval*num_layers;
	SPHxmin = 0-space_interval*(num_layers);
	SPHzmax = height_z+space_interval*num_layers;
	SPHzmin = 0-space_interval*(num_layers);
	Nx = (SPHxmax-SPHxmin)/(3.0*dem::smoothLength)+1;	// (xmax-xmin)/(3h)+1
	Nz = (SPHzmax-SPHzmin)/(3.0*dem::smoothLength)+1;	// (zmax-zmin)/(3h)+1
	numCell = Nx*Nz;
	SPHParticleCellVec.resize(numCell);	// at this point, the cellVec contains numCell vectors,
						// each vector is empty but ready to store the pointer of SPH particles (this process will be in the calculation of SPH forces)
	

/*	// move the right boundary just besides the water particles at x=L, i.e. x=x-(length_x-L-2*space_inter)
	dem::sphDamping = 0.0005;
    	dem::DMP_F = 2.0/dem::TIMESTEP;
    	dem::DMP_M = 2.0/dem::TIMESTEP;
	int breakStep = 100000;	// the step at which the dam will be breaking
if(length_x > L){
	REAL delta_x = length_x-L-2*space_interval;
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	    if( (*pt)->getCurrPosition().getx()>=length_x-small_value ){	// is the right boundary particle
		(*pt)->addCurrPositionX(-delta_x);
	    }
	}
}
*/

  	std::ofstream ofs_top("boundarySPHParticles.dat");
  	if(!ofs_top) {
    	    cout << "stream error!" << endl; exit(-1);
  	}
	ofs_top.setf(std::ios::scientific, std::ios::floatfield);
	ofs_top.precision(10);
	for(std::vector<sph::SPHParticle*>::const_iterator pt = SPHBoundaryParticleVec.begin(); pt!= SPHBoundaryParticleVec.end(); pt++) {
	    ofs_top << std::setw(20) << (*pt)->getCurrPosition().getx()
		    << std::setw(20) << (*pt)->getCurrPosition().gety()
		    << std::setw(20) << (*pt)->getCurrPosition().getz()
		    << std::endl;
	    ofs_top.flush();
	}
	ofs_top.close();

  	std::ofstream ofs_bot("freeSPHParticles.dat");
  	if(!ofs_bot) {
    	    cout << "stream error!" << endl; exit(-1);
  	}
	ofs_bot.setf(std::ios::scientific, std::ios::floatfield);
	ofs_bot.precision(10);
	for(std::vector<sph::SPHParticle*>::const_iterator pt = SPHParticleVec.begin(); pt!= SPHParticleVec.end(); pt++) {
	    ofs_bot << std::setw(20) << (*pt)->getCurrPosition().getx()
		    << std::setw(20) << (*pt)->getCurrPosition().gety()
		    << std::setw(20) << (*pt)->getCurrPosition().getz()
		    << std::endl;
	    ofs_bot.flush();
	}
	ofs_bot.close();

  	std::ofstream ofs_gho("ghostSPHParticles.dat");
  	if(!ofs_gho) {
    	    cout << "stream error!" << endl; exit(-1);
  	}
	ofs_gho.setf(std::ios::scientific, std::ios::floatfield);
	ofs_gho.precision(10);
	for(std::vector<particle*>::const_iterator demt=ParticleVec.begin(); demt!=ParticleVec.end(); demt++){
	    for(std::vector<sph::SPHParticle*>::const_iterator pt = (*demt)->SPHGhostParticleVec.begin(); pt!= (*demt)->SPHGhostParticleVec.end(); pt++) {
	    	ofs_gho << std::setw(20) << (*pt)->getCurrPosition().getx()
		    << std::setw(20) << (*pt)->getCurrPosition().gety()
		    << std::setw(20) << (*pt)->getCurrPosition().getz()
		    << std::endl;
	    	ofs_gho.flush();
	    }
	}
	ofs_gho.close();

	int iframe = 0;

	REAL waterHeight;
	REAL waterWidth;

    	REAL avgNormal=0;
    	REAL avgTangt=0;
    	int  stepsnum=0;
    	char stepsstr[4];
    	char stepsfp[50];
    	REAL void_ratio=0;
    	REAL bdry_penetr[7];
//    	matrix granularStress;	// granular stress
    	int  bdry_cntnum[7];
	REAL sph_K;	// sph kinetic energy

	dem::vec zeroVec = dem::vec(0,0,0);
	std::ofstream ofs("total_results.dat");
	// initialization SPH velocity based on equation(4.3) in http://www.artcompsci.org/vol_1/v1_web/node34.html
	initialSPHVelocity();
	for(int istep = 0; istep < nsteps; istep++){

/*	    // break the dam, let free SPH particles go
	    if(istep == breakStep){
		dem::sphDamping = 0;
		dem::DMP_F = 0.0/dem::TIMESTEP;
    		dem::DMP_M = 0.0/dem::TIMESTEP;
		dem::TIMESTEP = dem::TIMESTEP*0.02;	// change timestep to 2e-4
		for(std::vector<sph::SPHParticle*>::iterator ptt=SPHBoundaryParticleVec.begin(); ptt!=SPHBoundaryParticleVec.end(); ptt++){
	    	    (*ptt)->setCurrPositionInitial();
		}
	    }
*/
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
	    rigidBoundaryForce(bdry_penetr, bdry_cntnum);

// update position and density of SPH particles based on equation (4.1)
	    updateSPHLeapFrogPositionDensity();

	    calculateSPHDensityDotVelocityDotLinkedList2D();

//	    // add constant acceleration for all SPH points
//	    for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
//		(*pt)->addVelocityDot(acceleration);
//	    }
//	    applySPHBoundaryForce();

	    // fix the y for all SPH points
	    for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	      	(*pt)->fixY();
	    }
	    for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	      	(*pt)->fixXYZ();
	    }

	    // should fix the dem particle in y, i.e. force_j = 0; moment_i=moment_k = 0
	    for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
		(*pt)->fixY();
	    }

	    // 5. update particles' velocity/omga/position/orientation based on force/moment.
	    updateParticle();

// update velocity of SPH particles based on equation (4.2)
	    updateSPHLeapFrogVelocity();

	    if(istep % interval == 0){
		sprintf(stepsstr, "%03d", istep); 
	    	strcpy(stepsfp,"velocity"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
		printSPHVelocityPressure(stepsfp);	
		waterHeight = getWaterHeight();
		waterWidth = getWaterWidth();
	    	progressinf << setw(OWID) << istep
		            << setw(OWID) << getPossCntctNum()
		            << setw(OWID) << getActualCntctNum()
		            << setw(OWID) << getAveragePenetration()
		            << setw(OWID) << avgNormal
		            << setw(OWID) << avgTangt
		            << setw(OWID) << getAverageVelocity() 
		            << setw(OWID) << getAverageOmga()
		            << setw(OWID) << getAverageForce()   
		            << setw(OWID) << getAverageMoment()
			    << setw(OWID) << istep*dem::TIMESTEP
			    << setw(OWID) << waterHeight
			    << setw(OWID) << waterWidth
			    << std::endl;

		sph_K = calculateSPHKineticEnergy();
		std::cout << sph_K << std::endl;
	    } 

	    // output SPH particles information
	    if (istep % (nsteps/snapshots) == 0){
	    	writeSPHParticleTecplot(ofs, iframe);

	    	sprintf(stepsstr, "%03d", iframe); 
	    	strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    	printParticle(stepsfp);
	    
	    	sprintf(stepsstr, "%03d", iframe); 
	    	strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    	printContact(stepsfp);

	    	sprintf(stepsstr, "%03d", iframe); 
	    	strcpy(stepsfp, "SPH_free_particles"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    	printSPHFreeParticle(stepsfp);
	   	++iframe;
	    }
	}

	ofs.close();

	sprintf(stepsstr, "%03d", nsteps); 
	strcpy(stepsfp,"velocity"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	printSPHVelocityPressure(stepsfp);
	printSPHFreeParticle("SPH_free_particles_end");

    } // end BurstingDam2D



    void assembly::BurstingDam2DAfterEquilibrium(int nsteps,              // total_steps
	       			     int snapshots,           // number of snapshots   
	          		     int interval,
				     const char* iniptclfile,   
				     const char* inibdryfile,
				     const char* iniSPHptcl,
		       		     const char* particlefile, 
		       		     const char* contactfile,
		       		     const char* progressfile){

   	// pre_1: open streams for output.
    	// particlefile and contactfile are used for snapshots at the end.
    	progressinf.open(progressfile); 
    	if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
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
	            << setw(OWID) << "average" << endl
	            << setw(OWID) << "number"
	            << setw(OWID) << "contacts"
	            << setw(OWID) << "contacts"
	            << setw(OWID) << "penetration"
	            << setw(OWID) << "contact_normal"
	            << setw(OWID) << "contact_tangt"
	            << setw(OWID) << "velocity"
	            << setw(OWID) << "omga"
	            << setw(OWID) << "force"
	            << setw(OWID) << "moment" << endl;

    	// pre_2. create particles and boundaries from existing files.
    	readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    	readBoundary(inibdryfile);   // create boundaries.

	REAL L = 25;	// unit m
	length_x = 25;
	height_z = 50;
  	D = 5*(dem::G)*(dem::GRVT_SCL)*L;
	dem::SPHInitialDensity = 1000;	// unit kg/m^3
	dem::nu = 1e-6;	// unit m^2/s, dynamic viscosity
	nSPHPoint = 51;	// number of points in each direction, x and z
	int num_layers = 5;	// number of layers of boundary sph particles

	space_interval = L/(nSPHPoint-1);
	dem::smoothLength = 1.5*(space_interval);	// refer to Morris's paper Sect.2.9
	// calculate paramters used in kernel funciton
	one_devide_h = 1.0/dem::smoothLength;
	factor_kernel = 7.0/(478.0*dem::PI*dem::smoothLength*dem::smoothLength);
	factor_kernel_gradient = 7.0/(478.0*dem::PI*dem::smoothLength*dem::smoothLength*dem::smoothLength);

//	dem::soundSpeed = 1e-2;
  	dem::gamma = 7;
	dem::soundSpeed = sqrt(200*dem::G*dem::GRVT_SCL*L);
	dem::P0 = dem::soundSpeed*dem::soundSpeed*dem::SPHInitialDensity/dem::gamma*0.005;
	REAL qmin = 0.7593;	// this is fixed for quintic kernel function
	Wqmin = kernelFunction(qmin);
	p1 = 4;	// for the Lennard-Jones boundary forces 
	p2 = 2;
//        dem::vec initVelocity = dem::vec(-1.0e-4, 0, 0);	// the initial velocity for dem particle
//	dem::vec initOmega = dem::vec(0,0,-2e-2);

	// generate the sph particle and store them in SPHParticleVec
	REAL SPHmass = SPHInitialDensity*L*L*1.0/(nSPHPoint*nSPHPoint);

	// free the spaces of SPH points
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    delete (*pt);
	}
	SPHParticleVec.clear();
		// free the spaces of ghost points in each dem particle
	for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
	    for(std::vector<sph::SPHParticle*>::iterator spht=(*pt)->SPHGhostParticleVec.begin(); spht!=(*pt)->SPHGhostParticleVec.end(); spht++){
		delete (*spht);
	   }
	   (*pt)->SPHGhostParticleVec.clear();
	}
	// free the spaces of boundary SPH points
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	    delete (*pt);
	}
	SPHBoundaryParticleVec.clear();

	// create and store PeriParticle objects into periParticleVec
	int isGhost = 0;	// is Ghost particle
	REAL radius_a, radius_b, radius_c;
  	dem::vec pt_position;
	dem::vec local_x;	// local position of ghost point in dem particle
	dem::vec tmp_xyz, tmp_local;
	REAL small_value = 0.01*space_interval;
	for(REAL tmp_z=0-space_interval*(num_layers-1); tmp_z<height_z+space_interval*num_layers; tmp_z=tmp_z+space_interval){
	    for(REAL tmp_x=0-space_interval*(num_layers-1); tmp_x<length_x+space_interval*num_layers; tmp_x=tmp_x+space_interval){
		isGhost = 0;
		tmp_xyz = dem::vec(tmp_x, 0, tmp_z);
	    	for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
		    radius_a = (*pt)->getA(); radius_b = (*pt)->getB(); radius_c = (*pt)->getC();
		    pt_position = (*pt)->getCurrPosition();
		    pt_position.sety(0);
		    tmp_local = (*pt)->localVec(tmp_xyz-pt_position);
		    if( tmp_local.getx()*tmp_local.getx()/(radius_a*radius_a)+tmp_local.gety()*tmp_local.gety()/(radius_b*radius_b)+tmp_local.getz()*tmp_local.getz()/(radius_c*radius_c) <= 1 ){
		    	isGhost = 1;	// is Ghost particle
		    	local_x = (*pt)->localVec( dem::vec(tmp_x, 0, tmp_z)-pt_position );
		    	sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, local_x, (*pt), 2); 	// 2 is ghost SPH particle
		    	(*pt)->SPHGhostParticleVec.push_back(tmp_pt);
		    	break;	// if this sph point is ghost for dem particle 1, then it cannot be ghost for any others,
				// should avoid the initial overlap of the different dem particles
		    }
	    	} // end dem particle

		if(isGhost==0){	// free/boundary particles
		    if(tmp_x<=0+small_value || tmp_x>=length_x-small_value || tmp_z<=0+small_value || tmp_z>=height_z-small_value){	// boundary sph particles
			if(tmp_z>=L+space_interval){
		    	    sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, 3);	// 3 is boundary SPH particle
			    SPHBoundaryParticleVec.push_back(tmp_pt);
			}
			else{
		    	    sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity*pow(1+SPHInitialDensity*(dem::G)*(dem::GRVT_SCL)*(L+space_interval-tmp_z)/(dem::P0), 1.0/(dem::gamma)), tmp_x, 0, tmp_z, 3);	// 3 is boundary SPH particle
			    SPHBoundaryParticleVec.push_back(tmp_pt);
			}
		    } 
		}
	    }
	}

	// read SPH free particles
	std::ifstream ifs_sph(iniSPHptcl);
	if(!ifs_sph) {
	    std::cout << "stream error!" << std::endl; exit(-1);
	}

   	std::string tmp;
	ifs_sph >> nSPHPoint;	// here nSPHPoint is the total number of free SPH particles
	ifs_sph >> tmp >> tmp >> tmp >> tmp >> tmp;	// read header 

	REAL tmp_x, tmp_y, tmp_z;
	REAL tmp_mass, tmp_density;	
	isGhost = 0;	// is Ghost particle
	// read particle information, create and store PeriParticle objects into periParticleVec
	for(int ip=0; ip<nSPHPoint; ip++){
	    ifs_sph >> tmp_x >> tmp_y >> tmp_z >> tmp_mass >> tmp_density;
	    isGhost = 0;
	    tmp_xyz = dem::vec(tmp_x, 0, tmp_z);	
	    for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
	    	radius_a = (*pt)->getA(); radius_b = (*pt)->getB(); radius_c = (*pt)->getC();
		pt_position = (*pt)->getCurrPosition();
		pt_position.sety(0);
		tmp_local = (*pt)->localVec(tmp_xyz-pt_position);
		if( tmp_local.getx()*tmp_local.getx()/(radius_a*radius_a)+tmp_local.gety()*tmp_local.gety()/(radius_b*radius_b)+tmp_local.getz()*tmp_local.getz()/(radius_c*radius_c) <= 1 ){
		    isGhost = 1;	// is Ghost particle
		    local_x = (*pt)->localVec( dem::vec(tmp_x, 0, tmp_z)-pt_position );
		    sph::SPHParticle* tmp_pt = new sph::SPHParticle(tmp_mass, tmp_density, tmp_x, 0, tmp_z, local_x, (*pt), 2);	// 2 is ghost SPH particle
		    (*pt)->SPHGhostParticleVec.push_back(tmp_pt);
		    break;	// if this sph point is ghost for dem particle 1, then it cannot be ghost for any others,
				// should avoid the initial overlap of the different dem particles
		}
	    } // end dem particle
	    
	    if(isGhost==0){	// not Ghost particle
	    	sph::SPHParticle* tmp_pt = new sph::SPHParticle(tmp_mass, tmp_density, tmp_x, tmp_y, tmp_z, 1);		// 1 is free SPH particle
	    	SPHParticleVec.push_back(tmp_pt);
	    }

	}


	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    (*pt)->initial();
	}
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	    (*pt)->initial();
	}

  	std::ofstream ofs_top("boundarySPHParticles.dat");
  	if(!ofs_top) {
    	    cout << "stream error!" << endl; exit(-1);
  	}
	ofs_top.setf(std::ios::scientific, std::ios::floatfield);
	ofs_top.precision(10);
	for(std::vector<sph::SPHParticle*>::const_iterator pt = SPHBoundaryParticleVec.begin(); pt!= SPHBoundaryParticleVec.end(); pt++) {
	    ofs_top << std::setw(20) << (*pt)->getCurrPosition().getx()
		    << std::setw(20) << (*pt)->getCurrPosition().gety()
		    << std::setw(20) << (*pt)->getCurrPosition().getz()
		    << std::endl;
	    ofs_top.flush();
	}
	ofs_top.close();

  	std::ofstream ofs_bot("freeSPHParticles.dat");
  	if(!ofs_bot) {
    	    cout << "stream error!" << endl; exit(-1);
  	}
	ofs_bot.setf(std::ios::scientific, std::ios::floatfield);
	ofs_bot.precision(10);
	for(std::vector<sph::SPHParticle*>::const_iterator pt = SPHParticleVec.begin(); pt!= SPHParticleVec.end(); pt++) {
	    ofs_bot << std::setw(20) << (*pt)->getCurrPosition().getx()
		    << std::setw(20) << (*pt)->getCurrPosition().gety()
		    << std::setw(20) << (*pt)->getCurrPosition().getz()
		    << std::endl;
	    ofs_bot.flush();
	}
	ofs_bot.close();

  	std::ofstream ofs_gho("ghostSPHParticles.dat");
  	if(!ofs_gho) {
    	    cout << "stream error!" << endl; exit(-1);
  	}
	ofs_gho.setf(std::ios::scientific, std::ios::floatfield);
	ofs_gho.precision(10);
	for(std::vector<particle*>::const_iterator demt=ParticleVec.begin(); demt!=ParticleVec.end(); demt++){
	    for(std::vector<sph::SPHParticle*>::const_iterator pt = (*demt)->SPHGhostParticleVec.begin(); pt!= (*demt)->SPHGhostParticleVec.end(); pt++) {
	    	ofs_gho << std::setw(20) << (*pt)->getCurrPosition().getx()
		    << std::setw(20) << (*pt)->getCurrPosition().gety()
		    << std::setw(20) << (*pt)->getCurrPosition().getz()
		    << std::endl;
	    	ofs_gho.flush();
	    }
	}
	ofs_gho.close();

	int iframe = 0;

    	REAL avgNormal=0;
    	REAL avgTangt=0;
    	int  stepsnum=0;
    	char stepsstr[4];
    	char stepsfp[50];
    	REAL void_ratio=0;
    	REAL bdry_penetr[7];
//    	matrix granularStress;	// granular stress
    	int  bdry_cntnum[7];
	REAL sph_K;	// sph kinetic energy

	dem::vec zeroVec = dem::vec(0,0,0);
	std::ofstream ofs("total_results.dat");
	for(int istep = 0; istep < nsteps; istep++){

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
	    rigidBoundaryForce(bdry_penetr, bdry_cntnum);

	    calculateSPHDensityDotVelocityDot();
//	    // add constant acceleration for all SPH points
//	    for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
//		(*pt)->addVelocityDot(acceleration);
//	    }
//	    applySPHBoundaryForce();

	    // fix the y for all SPH points
	    for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	      	(*pt)->fixY();
	    }
	    for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	      	(*pt)->fixXYZ();
	    }

	    // should fix the dem particle in y, i.e. force_j = 0; moment_i=moment_k = 0
	    for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
		(*pt)->fixY();
	    }

	    // 5. update particles' velocity/omga/position/orientation based on force/moment.
	    updateParticle();
	    updateSPH();

	    if(istep % interval == 0){
		sprintf(stepsstr, "%03d", istep); 
	    	strcpy(stepsfp,"velocity"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
		printSPHVelocityPressure(stepsfp);	

	    	progressinf << setw(OWID) << istep
		            << setw(OWID) << getPossCntctNum()
		            << setw(OWID) << getActualCntctNum()
		            << setw(OWID) << getAveragePenetration()
		            << setw(OWID) << avgNormal
		            << setw(OWID) << avgTangt
		            << setw(OWID) << getAverageVelocity() 
		            << setw(OWID) << getAverageOmga()
		            << setw(OWID) << getAverageForce()   
		            << setw(OWID) << getAverageMoment()
			    << std::endl;

		sph_K = calculateSPHKineticEnergy();
		std::cout << sph_K << std::endl;
	    } 

	    // output SPH particles information
	    if (istep % (nsteps/snapshots) == 0){
	    	writeSPHParticleTecplot(ofs, iframe);

	    	sprintf(stepsstr, "%03d", iframe); 
	    	strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    	printParticle(stepsfp);
	    
	    	sprintf(stepsstr, "%03d", iframe); 
	    	strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    	printContact(stepsfp);

	    	sprintf(stepsstr, "%03d", iframe); 
	    	strcpy(stepsfp, "SPH_free_particles"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    	printSPHFreeParticle(stepsfp);
	   	++iframe;
	    }
	}

	ofs.close();

	sprintf(stepsstr, "%03d", nsteps); 
	strcpy(stepsfp,"velocity"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	printSPHVelocityPressure(stepsfp);
	printSPHFreeParticle("SPH_free_particles_end");

    } // end BurstingDam2DAfterEquilibrium


    // connecting containers
    void assembly::Connectivity2DAfterEquilibrium(int nsteps,              // total_steps
	       			     int snapshots,           // number of snapshots   
	          		     int interval,
				     const char* iniptclfile,   
				     const char* inibdryfile,
				     const char* iniSPHptcl,
		       		     const char* particlefile, 
		       		     const char* contactfile,
		       		     const char* progressfile){

   	// pre_1: open streams for output.
    	// particlefile and contactfile are used for snapshots at the end.
    	progressinf.open(progressfile); 
    	if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
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
	            << setw(OWID) << "average" << endl
	            << setw(OWID) << "number"
	            << setw(OWID) << "contacts"
	            << setw(OWID) << "contacts"
	            << setw(OWID) << "penetration"
	            << setw(OWID) << "contact_normal"
	            << setw(OWID) << "contact_tangt"
	            << setw(OWID) << "velocity"
	            << setw(OWID) << "omga"
	            << setw(OWID) << "force"
	            << setw(OWID) << "moment" << endl;

    	// pre_2. create particles and boundaries from existing files.
    	readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    	readBoundary(inibdryfile);   // create boundaries.

	REAL L = 25;	// unit m
	length_x = 50;
	height_z = 50;
  	D = 5*(dem::G)*(dem::GRVT_SCL)*L;
	dem::SPHInitialDensity = 1000;	// unit kg/m^3
	dem::nu = 1e-6;	// unit m^2/s, dynamic viscosity
	nSPHPoint = 51;	// number of points in each direction, x and z
	int num_layers = 5;	// number of layers of boundary sph particles
	REAL middle_plate_zmin = 5;	// the zmin of the middle plate
	REAL middle_plate_zmax = 20;

	space_interval = L/(nSPHPoint-1);
	dem::smoothLength = 1.5*(space_interval);	// refer to Morris's paper Sect.2.9
	// calculate paramters used in kernel funciton
	one_devide_h = 1.0/dem::smoothLength;
	factor_kernel = 7.0/(478.0*dem::PI*dem::smoothLength*dem::smoothLength);
	factor_kernel_gradient = 7.0/(478.0*dem::PI*dem::smoothLength*dem::smoothLength*dem::smoothLength);

//	dem::soundSpeed = 1e-2;
  	dem::gamma = 7;
	dem::soundSpeed = sqrt(200*dem::G*dem::GRVT_SCL*L);
	dem::P0 = dem::soundSpeed*dem::soundSpeed*dem::SPHInitialDensity/dem::gamma*0.005;
	REAL qmin = 0.7593;	// this is fixed for quintic kernel function
	Wqmin = kernelFunction(qmin);
	p1 = 4;	// for the Lennard-Jones boundary forces 
	p2 = 2;
//        dem::vec initVelocity = dem::vec(-1.0e-4, 0, 0);	// the initial velocity for dem particle
//	dem::vec initOmega = dem::vec(0,0,-2e-2);

	// generate the sph particle and store them in SPHParticleVec
	REAL SPHmass = SPHInitialDensity*L*L*1.0/(nSPHPoint*nSPHPoint);

	// free the spaces of SPH points
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    delete (*pt);
	}
	SPHParticleVec.clear();
		// free the spaces of ghost points in each dem particle
	for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
	    for(std::vector<sph::SPHParticle*>::iterator spht=(*pt)->SPHGhostParticleVec.begin(); spht!=(*pt)->SPHGhostParticleVec.end(); spht++){
		delete (*spht);
	   }
	   (*pt)->SPHGhostParticleVec.clear();
	}
	// free the spaces of boundary SPH points
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	    delete (*pt);
	}
	SPHBoundaryParticleVec.clear();

	// create and store PeriParticle objects into periParticleVec
	int isGhost = 0;	// is Ghost particle
	REAL radius_a, radius_b, radius_c;
  	dem::vec pt_position;
	dem::vec local_x;	// local position of ghost point in dem particle
	dem::vec tmp_xyz, tmp_local;
	REAL small_value = 0.01*space_interval;
	for(REAL tmp_z=0-space_interval*(num_layers-1); tmp_z<height_z+space_interval*num_layers; tmp_z=tmp_z+space_interval){
	    for(REAL tmp_x=0-space_interval*(num_layers-1); tmp_x<length_x+space_interval*num_layers; tmp_x=tmp_x+space_interval){
		isGhost = 0;
		tmp_xyz = dem::vec(tmp_x, 0, tmp_z);
	    	for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
		    radius_a = (*pt)->getA(); radius_b = (*pt)->getB(); radius_c = (*pt)->getC();
		    pt_position = (*pt)->getCurrPosition();
		    pt_position.sety(0);
		    tmp_local = (*pt)->localVec(tmp_xyz-pt_position);
		    if( tmp_local.getx()*tmp_local.getx()/(radius_a*radius_a)+tmp_local.gety()*tmp_local.gety()/(radius_b*radius_b)+tmp_local.getz()*tmp_local.getz()/(radius_c*radius_c) <= 1 ){
		    	isGhost = 1;	// is Ghost particle
		    	local_x = (*pt)->localVec( dem::vec(tmp_x, 0, tmp_z)-pt_position );
		    	sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, local_x, (*pt), 2);	// 2 is ghost SPH particle
		    	(*pt)->SPHGhostParticleVec.push_back(tmp_pt);
		    	break;	// if this sph point is ghost for dem particle 1, then it cannot be ghost for any others,
				// should avoid the initial overlap of the different dem particles
		    }
	    	} // end dem particle

		if(isGhost==0){	// free/boundary particles
		    if(tmp_x<=0+small_value || tmp_x>=length_x-small_value || tmp_z<=0+small_value || tmp_z>=height_z-small_value){	// boundary sph particles
			if(tmp_z>=L+space_interval){
		    	    sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, 3);	// 3 is boundary SPH particle
			    SPHBoundaryParticleVec.push_back(tmp_pt);
			}
			else{
		    	    sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity*pow(1+SPHInitialDensity*(dem::G)*(dem::GRVT_SCL)*(L+space_interval-tmp_z)/(dem::P0), 1.0/(dem::gamma)), tmp_x, 0, tmp_z, 3);	// 3 is boundary SPH particle
			    SPHBoundaryParticleVec.push_back(tmp_pt);
			}
		    }
		    else if(tmp_x<=L+space_interval && tmp_z<=L+space_interval){	// free sph particles
		    	sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity*pow(1+SPHInitialDensity*(dem::G)*(dem::GRVT_SCL)*(L+space_interval-tmp_z)/(dem::P0), 1.0/(dem::gamma)), tmp_x, 0, tmp_z, 1);	// 1 is free SPH particle
//			sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, 1);	// 1 is free SPH particle
	    	    	SPHParticleVec.push_back(tmp_pt);
		    }

		    if(tmp_x>=L+2*space_interval-small_value && tmp_x<=L+(2+num_layers-1)*space_interval+small_value && tmp_z>=middle_plate_zmin-small_value && tmp_z<=middle_plate_zmax+small_value ){	// the middle plate boundary
			if(tmp_z>=L+space_interval){
		    	    sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, 3);	// 3 is boundary SPH particle
			    SPHBoundaryParticleVec.push_back(tmp_pt);
			}
			else{
		    	    sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity*pow(1+SPHInitialDensity*(dem::G)*(dem::GRVT_SCL)*(L+space_interval-tmp_z)/(dem::P0), 1.0/(dem::gamma)), tmp_x, 0, tmp_z, 3);	// 3 is boundary SPH particle
			    SPHBoundaryParticleVec.push_back(tmp_pt);
			}
		    }
		}
	    }
	}

/*
	// read SPH free particles
	std::ifstream ifs_sph(iniSPHptcl);
	if(!ifs_sph) {
	    std::cout << "stream error!" << std::endl; exit(-1);
	}

   	std::string tmp;
	ifs_sph >> nSPHPoint;	// here nSPHPoint is the total number of free SPH particles
	ifs_sph >> tmp >> tmp >> tmp >> tmp >> tmp;	// read header 

	REAL tmp_x, tmp_y, tmp_z;
	REAL tmp_mass, tmp_density;	
	isGhost = 0;	// is Ghost particle
	// read particle information, create and store PeriParticle objects into periParticleVec
	for(int ip=0; ip<nSPHPoint; ip++){
	    ifs_sph >> tmp_x >> tmp_y >> tmp_z >> tmp_mass >> tmp_density;
	    isGhost = 0;
	    tmp_xyz = dem::vec(tmp_x, 0, tmp_z);	
	    for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
	    	radius_a = (*pt)->getA(); radius_b = (*pt)->getB(); radius_c = (*pt)->getC();
		pt_position = (*pt)->getCurrPosition();
		pt_position.sety(0);
		tmp_local = (*pt)->localVec(tmp_xyz-pt_position);
		if( tmp_local.getx()*tmp_local.getx()/(radius_a*radius_a)+tmp_local.gety()*tmp_local.gety()/(radius_b*radius_b)+tmp_local.getz()*tmp_local.getz()/(radius_c*radius_c) <= 1 ){
		    isGhost = 1;	// is Ghost particle
		    local_x = (*pt)->localVec( dem::vec(tmp_x, 0, tmp_z)-pt_position );
		    sph::SPHParticle* tmp_pt = new sph::SPHParticle(tmp_mass, tmp_density, tmp_x, 0, tmp_z, local_x, (*pt), 2);	// 2 is ghost SPH particle
		    (*pt)->SPHGhostParticleVec.push_back(tmp_pt);
		    break;	// if this sph point is ghost for dem particle 1, then it cannot be ghost for any others,
				// should avoid the initial overlap of the different dem particles
		}
	    } // end dem particle
	    
	    if(isGhost==0){	// not Ghost particle
	    	sph::SPHParticle* tmp_pt = new sph::SPHParticle(tmp_mass, tmp_density, tmp_x, tmp_y, tmp_z, 1);	// 1 is free SPH particle
	    	SPHParticleVec.push_back(tmp_pt);
	    }

	}
*/


	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    (*pt)->initial();
	}
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	    (*pt)->initial();
	}

  	std::ofstream ofs_top("boundarySPHParticles.dat");
  	if(!ofs_top) {
    	    cout << "stream error!" << endl; exit(-1);
  	}
	ofs_top.setf(std::ios::scientific, std::ios::floatfield);
	ofs_top.precision(10);
	for(std::vector<sph::SPHParticle*>::const_iterator pt = SPHBoundaryParticleVec.begin(); pt!= SPHBoundaryParticleVec.end(); pt++) {
	    ofs_top << std::setw(20) << (*pt)->getCurrPosition().getx()
		    << std::setw(20) << (*pt)->getCurrPosition().gety()
		    << std::setw(20) << (*pt)->getCurrPosition().getz()
		    << std::endl;
	    ofs_top.flush();
	}
	ofs_top.close();

  	std::ofstream ofs_bot("freeSPHParticles.dat");
  	if(!ofs_bot) {
    	    cout << "stream error!" << endl; exit(-1);
  	}
	ofs_bot.setf(std::ios::scientific, std::ios::floatfield);
	ofs_bot.precision(10);
	for(std::vector<sph::SPHParticle*>::const_iterator pt = SPHParticleVec.begin(); pt!= SPHParticleVec.end(); pt++) {
	    ofs_bot << std::setw(20) << (*pt)->getCurrPosition().getx()
		    << std::setw(20) << (*pt)->getCurrPosition().gety()
		    << std::setw(20) << (*pt)->getCurrPosition().getz()
		    << std::endl;
	    ofs_bot.flush();
	}
	ofs_bot.close();

  	std::ofstream ofs_gho("ghostSPHParticles.dat");
  	if(!ofs_gho) {
    	    cout << "stream error!" << endl; exit(-1);
  	}
	ofs_gho.setf(std::ios::scientific, std::ios::floatfield);
	ofs_gho.precision(10);
	for(std::vector<particle*>::const_iterator demt=ParticleVec.begin(); demt!=ParticleVec.end(); demt++){
	    for(std::vector<sph::SPHParticle*>::const_iterator pt = (*demt)->SPHGhostParticleVec.begin(); pt!= (*demt)->SPHGhostParticleVec.end(); pt++) {
	    	ofs_gho << std::setw(20) << (*pt)->getCurrPosition().getx()
		    << std::setw(20) << (*pt)->getCurrPosition().gety()
		    << std::setw(20) << (*pt)->getCurrPosition().getz()
		    << std::endl;
	    	ofs_gho.flush();
	    }
	}
	ofs_gho.close();

	int iframe = 0;

    	REAL avgNormal=0;
    	REAL avgTangt=0;
    	int  stepsnum=0;
    	char stepsstr[4];
    	char stepsfp[50];
    	REAL void_ratio=0;
    	REAL bdry_penetr[7];
//    	matrix granularStress;	// granular stress
    	int  bdry_cntnum[7];
	REAL sph_K;	// sph kinetic energy

	dem::vec zeroVec = dem::vec(0,0,0);
	std::ofstream ofs("total_results.dat");
	for(int istep = 0; istep < nsteps; istep++){

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
	    rigidBoundaryForce(bdry_penetr, bdry_cntnum);

	    calculateSPHDensityDotVelocityDot();
//	    // add constant acceleration for all SPH points
//	    for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
//		(*pt)->addVelocityDot(acceleration);
//	    }
//	    applySPHBoundaryForce();

	    // fix the y for all SPH points
	    for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	      	(*pt)->fixY();
	    }
	    for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	      	(*pt)->fixXYZ();
	    }

	    // should fix the dem particle in y, i.e. force_j = 0; moment_i=moment_k = 0
	    for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
		(*pt)->fixY();
	    }

	    // 5. update particles' velocity/omga/position/orientation based on force/moment.
	    updateParticle();
	    updateSPH();

	    if(istep % interval == 0){
		sprintf(stepsstr, "%03d", istep); 
	    	strcpy(stepsfp,"velocity"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
		printSPHVelocityPressure(stepsfp);	

	    	progressinf << setw(OWID) << istep
		            << setw(OWID) << getPossCntctNum()
		            << setw(OWID) << getActualCntctNum()
		            << setw(OWID) << getAveragePenetration()
		            << setw(OWID) << avgNormal
		            << setw(OWID) << avgTangt
		            << setw(OWID) << getAverageVelocity() 
		            << setw(OWID) << getAverageOmga()
		            << setw(OWID) << getAverageForce()   
		            << setw(OWID) << getAverageMoment()
			    << std::endl;

		sph_K = calculateSPHKineticEnergy();
		std::cout << sph_K << std::endl;
	    } 

	    // output SPH particles information
	    if (istep % (nsteps/snapshots) == 0){
	    	writeSPHParticleTecplot(ofs, iframe);

	    	sprintf(stepsstr, "%03d", iframe); 
	    	strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    	printParticle(stepsfp);
	    
	    	sprintf(stepsstr, "%03d", iframe); 
	    	strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    	printContact(stepsfp);

	    	sprintf(stepsstr, "%03d", iframe); 
	    	strcpy(stepsfp, "SPH_free_particles"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    	printSPHFreeParticle(stepsfp);
	   	++iframe;
	    }
	}

	ofs.close();

	sprintf(stepsstr, "%03d", nsteps); 
	strcpy(stepsfp,"velocity"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	printSPHVelocityPressure(stepsfp);
	printSPHFreeParticle("SPH_free_particles_end");

    } // end Connectivity2DAfterEquilibrium


    void assembly::BurstingDam3D(int nsteps,              // total_steps
	       			     int snapshots,           // number of snapshots   
	          		     int interval,
				     const char* iniptclfile,   
				     const char* inibdryfile,
		       		     const char* particlefile, 
		       		     const char* contactfile,
		       		     const char* progressfile){

   	// pre_1: open streams for output.
    	// particlefile and contactfile are used for snapshots at the end.
    	progressinf.open(progressfile); 
    	if(!progressinf) { cout << "stream error!" << endl; exit(-1); }
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
		    << setw(OWID) << "current"
		    << setw(OWID) << "current"
		    << setw(OWID) << "current" << endl
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
		    << setw(OWID) << "time"
		    << setw(OWID) << "height"
		    << setw(OWID) << "width" << endl;

    	// pre_2. create particles and boundaries from existing files.
    	readSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    	readBoundary(inibdryfile);   // create boundaries.

	REAL L = 25;	// unit m
	length_x = 100;
	height_z = 50;
  	D = 5*(dem::G)*(dem::GRVT_SCL)*L;
	dem::SPHInitialDensity = 1000;	// unit kg/m^3
	dem::nu = 1e-6;	// unit m^2/s, dynamic viscosity
	nSPHPoint = 51;	// number of points in each direction, x and z
	int num_layers = 5;	// number of layers of boundary sph particles

	space_interval = L/(nSPHPoint-1);	// calculate paramters used in kernel funciton
	dem::smoothLength = 1.5*(space_interval);	// refer to Morris's paper Sect.2.9
	one_devide_h = 1.0/dem::smoothLength;
	factor_kernel =1.0/(120.0*dem::PI*dem::smoothLength*dem::smoothLength*dem::smoothLength);
	factor_kernel_gradient = 1.0/(120.0*dem::PI*dem::smoothLength*dem::smoothLength*dem::smoothLength*dem::smoothLength);


//	dem::soundSpeed = 1e-2;
  	dem::gamma = 7;
	dem::soundSpeed = sqrt(200*dem::G*dem::GRVT_SCL*L);
	dem::P0 = dem::soundSpeed*dem::soundSpeed*dem::SPHInitialDensity/dem::gamma;
	dem::P0 = 1.013e5;	// air pressure
	REAL qmin = 0.7593;	// this is fixed for quintic kernel function
	Wqmin = kernelFunction(qmin);
	p1 = 4;	// for the Lennard-Jones boundary forces 
	p2 = 2;
//        dem::vec initVelocity = dem::vec(-1.0e-4, 0, 0);	// the initial velocity for dem particle
//	dem::vec initOmega = dem::vec(0,0,-2e-2);

	// generate the sph particle and store them in SPHParticleVec
	REAL SPHmass = SPHInitialDensity*L*L*L*1.0/(nSPHPoint*nSPHPoint*nSPHPoint);

	// free the spaces of SPH points
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    delete (*pt);
	}
	SPHParticleVec.clear();
		// free the spaces of ghost points in each dem particle
	for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
	    for(std::vector<sph::SPHParticle*>::iterator spht=(*pt)->SPHGhostParticleVec.begin(); spht!=(*pt)->SPHGhostParticleVec.end(); spht++){
		delete (*spht);
	   }
	   (*pt)->SPHGhostParticleVec.clear();
	}
	// free the spaces of boundary SPH points
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	    delete (*pt);
	}
	SPHBoundaryParticleVec.clear();

	// create and store PeriParticle objects into periParticleVec
	int isGhost = 0;	// is Ghost particle
	REAL radius_a, radius_b, radius_c;
  	dem::vec pt_position;
	dem::vec local_x;	// local position of ghost point in dem particle
	dem::vec tmp_xyz, tmp_local;
	REAL small_value = 0.01*space_interval;
    for(REAL tmp_y=0-space_interval*(num_layers); tmp_y<L+space_interval*num_layers; tmp_y=tmp_y+space_interval){
	for(REAL tmp_x=0-space_interval*(num_layers); tmp_x<length_x+space_interval*num_layers; tmp_x=tmp_x+space_interval){
	    for(REAL tmp_z=0-space_interval*(num_layers); tmp_z<height_z+space_interval*num_layers; tmp_z=tmp_z+space_interval){
		isGhost = 0;
		tmp_xyz = dem::vec(tmp_x, tmp_y, tmp_z);
	    	for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
		    radius_a = (*pt)->getA(); radius_b = (*pt)->getB(); radius_c = (*pt)->getC();
		    pt_position = (*pt)->getCurrPosition();
		    tmp_local = (*pt)->localVec(tmp_xyz-pt_position);
		    if( tmp_local.getx()*tmp_local.getx()/(radius_a*radius_a)+tmp_local.gety()*tmp_local.gety()/(radius_b*radius_b)+tmp_local.getz()*tmp_local.getz()/(radius_c*radius_c) <= 1 ){
		    	isGhost = 1;	// is Ghost particle
		    	local_x = (*pt)->localVec( dem::vec(tmp_x, tmp_y, tmp_z)-pt_position );
		    	sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, tmp_y, tmp_z, local_x, (*pt), 2);	// 2 is ghost SPH particle
		    	(*pt)->SPHGhostParticleVec.push_back(tmp_pt);
		    	break;	// if this sph point is ghost for dem particle 1, then it cannot be ghost for any others,
				// should avoid the initial overlap of the different dem particles
		    }
	    	} // end dem particle

		if(isGhost==0){	// free/boundary particles
		    if(tmp_x<=-space_interval+small_value || tmp_x>=length_x-small_value || tmp_z<=-space_interval+small_value || tmp_z>=height_z-small_value){	// boundary sph particles
			if(tmp_z>=L+space_interval){
		    	    sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, 3);	// 3 is boundary SPH particle
			    SPHBoundaryParticleVec.push_back(tmp_pt);
			}
			else{
		    	    sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity*pow(1+SPHInitialDensity*(dem::G)*(dem::GRVT_SCL)*(L+space_interval-tmp_z)/(dem::P0), 1.0/(dem::gamma)), tmp_x, 0, tmp_z, 3);	// 3 is boundary SPH particle
			    SPHBoundaryParticleVec.push_back(tmp_pt);
			}
		    } 
		    else if(tmp_x<=L && tmp_z<=L){	// free sph particles
		    	sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity*pow(1+SPHInitialDensity*(dem::G)*(dem::GRVT_SCL)*(L+space_interval-tmp_z)/(dem::P0), 1.0/(dem::gamma)), tmp_x, 0, tmp_z, 1);	// 1 is free SPH particle
//			sph::SPHParticle* tmp_pt = new sph::SPHParticle(SPHmass, SPHInitialDensity, tmp_x, 0, tmp_z, 1);	// 1 is free SPH particle
	    	    	SPHParticleVec.push_back(tmp_pt);
		    }
		}
	    }
	}
    }

	for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	    (*pt)->initial();
	}
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	    (*pt)->initial();
	}

	// generate the vector of cells for SPH particles in order to use linked-list search
	SPHxmax = length_x+space_interval*num_layers;
	SPHxmin = 0-space_interval*(num_layers);
	SPHzmax = height_z+space_interval*num_layers;
	SPHzmin = 0-space_interval*(num_layers);
	Nx = (SPHxmax-SPHxmin)/(3.0*dem::smoothLength)+1;	// (xmax-xmin)/(3h)+1
	Nz = (SPHzmax-SPHzmin)/(3.0*dem::smoothLength)+1;	// (zmax-zmin)/(3h)+1
	numCell = Nx*Nz;
	SPHParticleCellVec.resize(numCell);	// at this point, the cellVec contains numCell vectors,
						// each vector is empty but ready to store the pointer of SPH particles (this process will be in the calculation of SPH forces)
	

/*	// move the right boundary just besides the water particles at x=L, i.e. x=x-(length_x-L-2*space_inter)
	dem::sphDamping = 0.0005;
    	dem::DMP_F = 2.0/dem::TIMESTEP;
    	dem::DMP_M = 2.0/dem::TIMESTEP;
	int breakStep = 100000;	// the step at which the dam will be breaking
if(length_x > L){
	REAL delta_x = length_x-L-2*space_interval;
	for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	    if( (*pt)->getCurrPosition().getx()>=length_x-small_value ){	// is the right boundary particle
		(*pt)->addCurrPositionX(-delta_x);
	    }
	}
}
*/

  	std::ofstream ofs_top("boundarySPHParticles.dat");
  	if(!ofs_top) {
    	    cout << "stream error!" << endl; exit(-1);
  	}
	ofs_top.setf(std::ios::scientific, std::ios::floatfield);
	ofs_top.precision(10);
	for(std::vector<sph::SPHParticle*>::const_iterator pt = SPHBoundaryParticleVec.begin(); pt!= SPHBoundaryParticleVec.end(); pt++) {
	    ofs_top << std::setw(20) << (*pt)->getCurrPosition().getx()
		    << std::setw(20) << (*pt)->getCurrPosition().gety()
		    << std::setw(20) << (*pt)->getCurrPosition().getz()
		    << std::endl;
	    ofs_top.flush();
	}
	ofs_top.close();

  	std::ofstream ofs_bot("freeSPHParticles.dat");
  	if(!ofs_bot) {
    	    cout << "stream error!" << endl; exit(-1);
  	}
	ofs_bot.setf(std::ios::scientific, std::ios::floatfield);
	ofs_bot.precision(10);
	for(std::vector<sph::SPHParticle*>::const_iterator pt = SPHParticleVec.begin(); pt!= SPHParticleVec.end(); pt++) {
	    ofs_bot << std::setw(20) << (*pt)->getCurrPosition().getx()
		    << std::setw(20) << (*pt)->getCurrPosition().gety()
		    << std::setw(20) << (*pt)->getCurrPosition().getz()
		    << std::endl;
	    ofs_bot.flush();
	}
	ofs_bot.close();

  	std::ofstream ofs_gho("ghostSPHParticles.dat");
  	if(!ofs_gho) {
    	    cout << "stream error!" << endl; exit(-1);
  	}
	ofs_gho.setf(std::ios::scientific, std::ios::floatfield);
	ofs_gho.precision(10);
	for(std::vector<particle*>::const_iterator demt=ParticleVec.begin(); demt!=ParticleVec.end(); demt++){
	    for(std::vector<sph::SPHParticle*>::const_iterator pt = (*demt)->SPHGhostParticleVec.begin(); pt!= (*demt)->SPHGhostParticleVec.end(); pt++) {
	    	ofs_gho << std::setw(20) << (*pt)->getCurrPosition().getx()
		    << std::setw(20) << (*pt)->getCurrPosition().gety()
		    << std::setw(20) << (*pt)->getCurrPosition().getz()
		    << std::endl;
	    	ofs_gho.flush();
	    }
	}
	ofs_gho.close();

	int iframe = 0;

	REAL waterHeight;
	REAL waterWidth;

    	REAL avgNormal=0;
    	REAL avgTangt=0;
    	int  stepsnum=0;
    	char stepsstr[4];
    	char stepsfp[50];
    	REAL void_ratio=0;
    	REAL bdry_penetr[7];
//    	matrix granularStress;	// granular stress
    	int  bdry_cntnum[7];
	REAL sph_K;	// sph kinetic energy

	dem::vec zeroVec = dem::vec(0,0,0);
	std::ofstream ofs("total_results.dat");
	// initialization SPH velocity based on equation(4.3) in http://www.artcompsci.org/vol_1/v1_web/node34.html
	initialSPHVelocity();
	for(int istep = 0; istep < nsteps; istep++){

/*	    // break the dam, let free SPH particles go
	    if(istep == breakStep){
		dem::sphDamping = 0;
		dem::DMP_F = 0.0/dem::TIMESTEP;
    		dem::DMP_M = 0.0/dem::TIMESTEP;
		dem::TIMESTEP = dem::TIMESTEP*0.02;	// change timestep to 2e-4
		for(std::vector<sph::SPHParticle*>::iterator ptt=SPHBoundaryParticleVec.begin(); ptt!=SPHBoundaryParticleVec.end(); ptt++){
	    	    (*ptt)->setCurrPositionInitial();
		}
	    }
*/
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
	    rigidBoundaryForce(bdry_penetr, bdry_cntnum);

// update position and density of SPH particles based on equation (4.1)
	    updateSPHLeapFrogPositionDensity();

	    calculateSPHDensityDotVelocityDotLinkedList2D();
//	    // add constant acceleration for all SPH points
//	    for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
//		(*pt)->addVelocityDot(acceleration);
//	    }
//	    applySPHBoundaryForce();

	    // fix the y for all SPH points
	    for(std::vector<sph::SPHParticle*>::iterator pt=SPHParticleVec.begin(); pt!=SPHParticleVec.end(); pt++){
	      	(*pt)->fixY();
	    }
	    for(std::vector<sph::SPHParticle*>::iterator pt=SPHBoundaryParticleVec.begin(); pt!=SPHBoundaryParticleVec.end(); pt++){
	      	(*pt)->fixXYZ();
	    }

	    // should fix the dem particle in y, i.e. force_j = 0; moment_i=moment_k = 0
	    for(std::vector<particle*>::iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
		(*pt)->fixY();
	    }

	    // 5. update particles' velocity/omga/position/orientation based on force/moment.
	    updateParticle();

// update velocity of SPH particles based on equation (4.2)
	    updateSPHLeapFrogVelocity();

	    if(istep % interval == 0){
		sprintf(stepsstr, "%03d", istep); 
	    	strcpy(stepsfp,"velocity"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
		printSPHVelocityPressure(stepsfp);	
		waterHeight = getWaterHeight();
		waterWidth = getWaterWidth();
	    	progressinf << setw(OWID) << istep
		            << setw(OWID) << getPossCntctNum()
		            << setw(OWID) << getActualCntctNum()
		            << setw(OWID) << getAveragePenetration()
		            << setw(OWID) << avgNormal
		            << setw(OWID) << avgTangt
		            << setw(OWID) << getAverageVelocity() 
		            << setw(OWID) << getAverageOmga()
		            << setw(OWID) << getAverageForce()   
		            << setw(OWID) << getAverageMoment()
			    << setw(OWID) << istep*dem::TIMESTEP
			    << setw(OWID) << waterHeight
			    << setw(OWID) << waterWidth
			    << std::endl;

		sph_K = calculateSPHKineticEnergy();
		std::cout << sph_K << std::endl;
	    } 

	    // output SPH particles information
	    if (istep % (nsteps/snapshots) == 0){
	    	writeSPHParticleTecplot(ofs, iframe);

	    	sprintf(stepsstr, "%03d", iframe); 
	    	strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    	printParticle(stepsfp);
	    
	    	sprintf(stepsstr, "%03d", iframe); 
	    	strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    	printContact(stepsfp);

	    	sprintf(stepsstr, "%03d", iframe); 
	    	strcpy(stepsfp, "SPH_free_particles"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    	printSPHFreeParticle(stepsfp);
	   	++iframe;
	    }
	}

	ofs.close();

	sprintf(stepsstr, "%03d", nsteps); 
	strcpy(stepsfp,"velocity"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	printSPHVelocityPressure(stepsfp);
	printSPHFreeParticle("SPH_free_particles_end");

    } // end BurstingDam3D



} // namespace dem
