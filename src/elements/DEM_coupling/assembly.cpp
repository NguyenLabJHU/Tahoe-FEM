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
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <cassert>
#include <utility>

#include "GhostParticleT.h"

#ifdef OPENMP
#include <omp.h>
#define OPENMP_IMPL 4
#endif

// OPENMP_IMPL:
// 0: OpenMP implementation 0, ts partitions, based on linked list
// 1: OpenMP implementation 1, ts partitions, based on vector
// 2: OpenMP implementation 2, no partition, each thread leaps by ts until completed
// 3: OpenMP implementation 3, no partition, each thread leaps by ts until num/2 and handles two particles.
// 4: OpenMP implementation 4, no partition, parallel for, various loop scheduling: (static), (static,1), (dynamic), (dynamic,1)

//#define BINNING
//#define TIME_PROFILE

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
	MaxRadius=pt->getA()>MaxRadius?pt->getA():MaxRadius;

	ParticleVec.push_back(pt);
    }
    ifs.close();
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
	else       // cylindrical boundary
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
      FaceVec.clear();

      ifs>>RgdBdryNum;

      for(int i=0;i<RgdBdryNum;i++){
	ifs>>type;
	if(type==1){ // plane boundary
	  rbptr=new plnrgd_bdry<particle>(ifs);
	  FaceVec.push_back(rbptr);}
	else        // cylindrical boundary -by cupescathy
	  rbptr=new cylrgd_bdry<particle>(ifs);
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

/*void assembly::findParticleOnTunnelFace(){  // find particle on tunnel face which is part of cavity boundaries --by cupescathy
    std::vector<RGDBDRY*>::iterator rt;
    for(rt = FaceVec.begin(); rt != FaceVec.end(); ++rt)
	(*rt)->findParticleOnBoundary(ParticleVec);
}*/

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
  for(rt = CavityRBVec.begin(); rt != CavityRBVec.end(); ++rt){
	  if ((*rt)->getOrder()==1)
		  (*rt)->rigidBF(BdryTgtMap);
	  else
		  (*rt)->rigidBF();
  }
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

// continue deposition with deposited particles and additional generated particles -by cupescathy
void assembly::deposit_Continue(int   freetype,
			       int   total_steps,
			       int   snapshots,
			       int   interval,
			       REAL  rFloHeight,
			       const char* iniptclfile,
			       const char* inibdryfile,
			       const char* layer_input,// additional generation
			       const char* particlefile,
			       const char* contactfile,
			       const char* progressfile,
			       const char* trmparticle,
			       const char* trmboundary,
			       const char* debugfile)
{
  generateAll(iniptclfile, layer_input, freetype, rFloHeight);

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
  REAL diameter = gradInfo.getMaxPtclRadius()*2;

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

// generate particles by combining deposited and additional particles -by cupescathy
void assembly::generateAll(const char* particlefile,
		    const char* layer_input,
			int freetype,
			REAL rFloHeight)
{
  readSample(layer_input);
  REAL x,y,z;
  particle* newptcl;
  REAL dimx     = container.getDimx();
  REAL dimy     = container.getDimy();
  REAL dimz     = container.getDimz();
  REAL diameter = gradInfo.getMaxPtclRadius()*2;

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

// create a specimen from discrete particles through floating and then gravitation,
// boundaries are composed of fixed particles.
void assembly::deposit_PtclBdry(gradation& grad,
				int   freetype,
				int   total_steps,
				int   snapshots,
				int   interval,
				REAL  rFloHeight,
				REAL  cof,
				const char* iniptclfile,
				const char* particlefile,
				const char* contactfile,
				const char* progressfile,
				const char* debugfile)
{
	generate_p(grad, iniptclfile, freetype, rFloHeight,cof);

	deposit_p(total_steps,        // total_steps
		  snapshots,          // number of snapshots
		  interval,           // print interval
		  iniptclfile,        // input file, initial particles
		  particlefile,       // output file, resulted particles, including snapshots
		  contactfile,        // output file, resulted contacts, including snapshots
		  progressfile,       // output file, statistical info
		  debugfile);         // output file, debug info
    }

// freetype:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
// ht- how many times of size would be the floating height
void assembly::generate_p(gradation&  grad,
			 const char* particlefile,
			 int freetype,
			 REAL rFloHeight,
			 REAL cof)  // ratio of diameter of gost particle to maximum diameter of gradation
{

	// get container's information. This will be used for particle generation including fixed particles
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
	REAL z2 = container.getMaxCorner().getz() + diameter;

	REAL x0 = container.getCenter().getx();
	REAL y0 = container.getCenter().gety();
	REAL z0 = container.getCenter().getz();

    REAL wall=2.2; // wall - wall height; ht - free particle height
    REAL est =1.05;
    REAL gdiameter=diameter*cof;

    // particle boundary 1
    x=container.getMaxCorner().getx();
    for (y=container.getMinCorner().gety(); y<2*y0*est; y+=gdiameter)
	for(z=container.getMinCorner().getz()+gdiameter; z< 2*z0*wall; z+=gdiameter){
	    newptcl = new particle(TotalNum+1, 1, vec(x,y,z), gdiameter/2, YOUNG, POISSON);
	    ParticleVec.push_back(newptcl);
	    TotalNum++;
	}

    // particle boundary 2
    y=container.getMaxCorner().gety();
    for (x=container.getMinCorner().getx(); x<2*x0*est; x+=gdiameter)
	for(z=container.getMinCorner().getz()+gdiameter; z<2*z0*wall; z+=gdiameter){
	    newptcl = new particle(TotalNum+1, 1, vec(x,y,z), gdiameter/2, YOUNG, POISSON);
	    ParticleVec.push_back(newptcl);
	    TotalNum++;
	}

    // particle boundary 3
    x=container.getMinCorner().getx();
    for (y=container.getMinCorner().gety(); y<2*y0*est; y+=gdiameter)
	for(z=container.getMinCorner().getz()+gdiameter; z<2*z0*wall; z+=gdiameter){
	    newptcl = new particle(TotalNum+1, 1, vec(x,y,z), gdiameter/2, YOUNG, POISSON);
	    ParticleVec.push_back(newptcl);
	    TotalNum++;
	}

    // particle boundary 4
    y=container.getMinCorner().gety();
    for (x=container.getMinCorner().getx(); x<2*x0*est; x+=gdiameter)
	for(z=container.getMinCorner().getz()+gdiameter; z<2*z0*wall; z+=gdiameter){
	    newptcl = new particle(TotalNum+1, 1, vec(x,y,z), gdiameter/2, YOUNG, POISSON);
	    ParticleVec.push_back(newptcl);
	    TotalNum++;
	}

    // particle boundary 6
    z=container.getMinCorner().getz();
    for (y=container.getMinCorner().gety(); y<2*y0*est; y+=gdiameter)
	for( x=container.getMinCorner().getx(); x<2*x0*est; x+=gdiameter){
	    newptcl = new particle(TotalNum+1, 1, vec(x,y,z), gdiameter/2, YOUNG, POISSON);
	    ParticleVec.push_back(newptcl);
	    TotalNum++;
	}

    // generate free particles
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
    //readSample(iniptclfile);   // create container and particles, velocity and omga are set zero.
    //readBoundary(inibdryfile); // create boundaries.
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
	        << "sigma3_1        sigma3_2           p             width          length           "
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
	l56=getTopFreeParticlePosition().getz() - getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx(); Volume=l13*l24*l56;
	void_ratio=Volume/getParticleVolume()-1;

	/*
	l56=getTopFreeParticlePosition().getz() - (-dimn/2);
	l24=dimn*rsize;
	l13=dimn*rsize;
	Volume=l13*l24*l56;
	void_ratio=Volume/getParticleVolume()-1;
	*/

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

	// 8. loop break condition
	if (   fabs(sigma1_1-sigma)/sigma < STRESS_ERROR && fabs(sigma1_2-sigma)/sigma < STRESS_ERROR
	    && fabs(sigma2_1-sigma)/sigma < STRESS_ERROR && fabs(sigma2_2-sigma)/sigma < STRESS_ERROR
	    && fabs(sigma3_1-sigma)/sigma < STRESS_ERROR && fabs(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
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
	if (g_iteration % interval == 0 ){
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
	                << setw(OWID) << getVibraTimeStep()
	                << setw(OWID) << getImpactTimeStep()
		        << setw(OWID) << timediffsec(time_w1,time_w2)
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

	// Most time it is balanced, so use progressinf instead.
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

//trim particles by given value of height-by Boning
void assembly::trimByHeight(REAL Height,
			const char* inptclefile,	// initial particle file
			const char* inbdryfile,		// initial boundary file
			const char* resultptcle){	// result particle file

	readSample(inptclefile);	// create container and particles, velocity and omga are set zero.
	readBoundary(inbdryfile);	// create boundaries

	// pre_3. define variables used in iterations
	REAL bdryW = getApt(2).gety()-getApt(4).gety();
    	REAL bdryL = getApt(1).getx()-getApt(3).getx();
	REAL bdryH = getApt(5).getz()-getApt(6).getz();	// height of boundary box

	trimHistoryNum = TotalNum;

	//update container information
	REAL dx,dy,dz,cx,cy,cz;
	dx=container.getDimx();
	dy=container.getDimy();
	dz=Height;
	cx=dx/2;
	cy=dy/2;
	cz=dz/2;
	container.set(container.getDimx(),container.getDimy(),Height,vec(cx,cy,cz));

 	REAL z1, z2;
 	z1 = getApt(6).getz();
	z2 = z1+Height;

  	std::vector<particle*>::iterator itr;
  	vec center;
  	REAL mass = 0;
	int num_iter = 0;
  	for (itr = ParticleVec.begin(); itr != ParticleVec.end(); ){
    		center=(*itr)->getCurrPosition();
    		if(center.getz() > z2)
      		{
			delete (*itr); // release memory
			itr = ParticleVec.erase(itr);
	      	}
    		else
      			++itr;
		num_iter++;
//		std::cout << num_iter << std::endl;
 	 }

  	for(itr=ParticleVec.begin();itr!=ParticleVec.end();++itr)
    		mass += (*itr)->getMass();

  	Volume = container.getDimx() * container.getDimy() * container.getDimz();
  	BulkDensity = mass/Volume;

  	TotalNum = ParticleVec.size();
  	printParticle(resultptcle);

}

//generate new sample by an existing deposited sample-by cupescathy
void assembly::gnrtPtclBySample(const char* samplefile,
		REAL row,
		REAL column,
		REAL layer,
		const char* resultsample,
		REAL offset,
		int bdrynum,
		const char* boundaryfile)
{
	std::ifstream ifs(samplefile);
	REAL x,y,z,x0,y0,z0,cx,cy,cz,dx,dy,dz;
	int k=0;//total number of basic sample
	int RORC;
	int ID, type, TN;
	int kx=0,ky=0,kz=0;// counter of sample number in each direction
	REAL a, b, c, px,py,pz,dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz;
	REAL vx,vy,vz,omx,omy,omz,fx,fy,fz,mx,my,mz;

	ParticleVec.clear();//clear existing particle information

	for (z=1;z<=layer;z++){
		kz++;
		for (y=1;y<=column;y++){
			ky++;
			for (x=1;x<=row;x++){
				kx++;
				k++;
				std::ifstream ifs(samplefile);

				ifs >> TotalNum >> RORC;

				//generate new container
				ifs >> cx >> cy >> cz >> dx >> dy >> dz;
				container.set(dx,dy,dz,vec(cx,cy,cz));

				x0=container.getDimx();
				y0=container.getDimy();
				z0=container.getDimz();

				dx=dx+x0*(kx-1);
				dy=dy+y0*(ky-1);
				dz=dz+z0*(kz-1);
				cx=cx+x0*(kx-1)/2;
				cy=cy+y0*(ky-1)/2;
				cz=cz+z0*(kz-1)/2;
				container.set(dx,dy,dz,vec(cx,cy,cz));//update information of container
				Volume = container.getDimx() * container.getDimy() * container.getDimz();
				TN=TotalNum;
				char s[20];
				ifs>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
				   >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s;
				for (int i=0;i<TotalNum;i++){
					//generate new particles
					ifs>>ID>>type>>a>>b>>c>>px>>py>>pz>>dax>>day>>daz>>dbx>>dby>>dbz>>dcx>>dcy>>dcz
					   >>vx>>vy>>vz>>omx>>omy>>omz>>fx>>fy>>fz>>mx>>my>>mz;
					ID=ID+TN*(k-1);
					px=px+x0*(kx-1);
					py=py+y0*(ky-1);
					pz=pz+(z0+offset)*(kz-1);
					particle* pt= new particle(ID,type,vec(a,b,c),vec(px,py,pz),vec(dax,day,daz),vec(dbx,dby,dbz),vec(dcx,dcy,dcz),YOUNG,POISSON);
					ParticleVec.push_back(pt);
		       }
				ifs.close();
			}
		kx=0;
		}
	ky=0;
	}
	//BulkDensity = mass/Volume;
 	TotalNum = ParticleVec.size();
	printParticle(resultsample);
	buildBoundary(bdrynum,boundaryfile);
}

// make a cylinder cavity (which can be treated as a tunnel) inside the sample and remove particles in the cavity-by cupescathy (tested)
void assembly::trimCavCyl(bool toRebuild,
			  const char* particlefile,
			  const char* cavparticlefile)
{
  if (toRebuild) readSample(particlefile);
  trimHistoryNum = TotalNum;

  REAL r0,h0,r,h,x_0,y_0,z_0,x_1,y_1,z_1;
  r0 = cavcyl.getRadius();
  h0 = cavcyl.getHeight();
  x_0 = cavcyl.getCenter().getx();
  y_0 = cavcyl.getCenter().gety();
  z_0 = cavcyl.getCenter().getz();

  /*r=0;
  h=0;
  x_1=0;
  y_1=0;
  z_1=0;*/

  std::vector<particle*>::iterator itr;
  vec center;
  REAL mass = 0;
  REAL delta;

  for (itr = ParticleVec.begin(); itr != ParticleVec.end(); ){
    center=(*itr)->getCurrPosition();
    x_1=center.getx();
    y_1=center.gety();
    z_1=center.getz();
    delta=(*itr)->getA();
    r=sqrt((y_0-y_1)*(y_0-y_1)+(z_0-z_1)*(z_0-z_1)); // a cylinder domain with an axis along x direction
    h=fabs(x_1-x_0);
    if((r-delta) < r0 && (h-delta)< h0/2)
      {
	delete (*itr); // release memory
	itr=ParticleVec.erase(itr);
      }
    else
      ++itr;
  }

  for(itr=ParticleVec.begin();itr!=ParticleVec.end();++itr)
    mass += (*itr)->getMass();

  Volume = container.getDimx() * container.getDimy() * container.getDimz()
          -cavcyl.getVolume();
  BulkDensity = mass/Volume;

  TotalNum = ParticleVec.size();
  printParticle(cavparticlefile);
}


// build cylindrical cavity boundary -by cupescathy
void assembly::buildCavCylBoundary(int existMaxId, const char* boundaryfile)
{
  std::ofstream ofs(boundaryfile);
  if(!ofs) { cout << "stream error!" << endl; exit(-1);}

  REAL r,h,x0,y0,z0;
  r = cavcyl.getRadius();
  h = cavcyl.getHeight();
  x0 = cavcyl.getCenter().getx();
  y0 = cavcyl.getCenter().gety();
  z0 = cavcyl.getCenter().getz();

  int bdrynum = 2;

  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs << setw(OWID) << 0
      << setw(OWID) << bdrynum << endl << endl;

  // boundary 1 (cylindrical)
  ofs << setw(OWID) << 2 << endl // cylindrical boundary -by cupescathy
      << setw(OWID) << existMaxId + 1
      << setw(OWID) << 2
      << setw(OWID) << 0 << endl

      << setw(OWID) << 2
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << x0
      << setw(OWID) << y0
      << setw(OWID) << z0
      << setw(OWID) << r
      << setw(OWID) << -1 << endl

      << setw(OWID) << 1
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << x0-h/2
      << setw(OWID) << y0
      << setw(OWID) << z0
      << setw(OWID) << 0
      << setw(OWID) << 0 << endl<< endl

   // boundary 2 (plane)
      << setw(OWID) << 1 << endl           // type
      << setw(OWID) << existMaxId + 2      // ID
      << setw(OWID) << 2                   // limnum
      << setw(OWID) << 0 << endl           // area

      << setw(OWID) << 1                   // order
      << setw(OWID) << 1                   // dirc.x
      << setw(OWID) << 0                   // dirc.y
      << setw(OWID) << 0                   // dirc.z
      << setw(OWID) << x0-h/2              // apt.x
      << setw(OWID) << y0                  // apt.y
      << setw(OWID) << z0                  // apt.z
      << setw(OWID) << 0                   // radius
      << setw(OWID) << 0 << endl           // side

      << setw(OWID) << 2
      << setw(OWID) << 1
      << setw(OWID) << 0
      << setw(OWID) << 0
      << setw(OWID) << x0
      << setw(OWID) << y0
      << setw(OWID) << z0
      << setw(OWID) << r
      << setw(OWID) << -1 << endl<< endl;

      ofs.close();
}

// tunneling sample generation -by cupescathy
void assembly::TunnelingSample(vec cntrposi, // tunnel center position which locates at certain boundary
		const char* iniprtls, // initial particle file
		const char* inibdrys, // initial boundary file
		const char* modelprtls, //model particle file
		const char* tunnelbdrys, // tunnel boundary file
		REAL D, // tunnel diameter
		REAL lgth, // initial tunneling distance
		int   total_steps,  // used for initial deposition after particles are trimmed
		int   snapshots,
		int   interval,
		const char* particlefile,
		const char* contactfile,
		const char* progressfile,
		const char* debugfile
){

	readSample(iniprtls);
	readBoundary(inibdrys); //create initial boundary which does not include the cylindrical boundary

	// step_1: generate geometry information of the tunnel

	vec tmp;
	tmp=vec(cntrposi.getx()-lgth/2, cntrposi.gety(), cntrposi.getz()); // get the center of cylinder
	cavcyl.setRadius(D/2);
	cavcyl.setHeight(lgth);
	cavcyl.setCenter(tmp);

	// step_2: generate the tunnel inside particles by cavity trimming function
	trimCavCyl(true, iniprtls, modelprtls);

	// step_3: create boundary of tunnel
	buildCavCylBoundary(5, tunnelbdrys);

	// step_4: build initial model through deposition

	depositAfterCavity(total_steps,
			snapshots,
			interval,
			modelprtls,
			inibdrys,
			tunnelbdrys,
			particlefile,
			contactfile,
			progressfile,
			debugfile);

}

// tunneling function. Typically tunneling is a deposition process -by cupescathy
void assembly::Tunneling (const char* samplefile, // generated tunneling sample
		const char* bdryfile, // boundary file without cavity file
		const char* cavitybdry, // tunnel boundary which treated as a cylinder with outside particles
		REAL exspeed,    // excavation speed
		REAL Height, // initial tunnel length
		REAL Radius, // tunnel radius
		REAL Cof,  // cofficient of over-excavation and under-excavation
		int   total_steps,
		int   snapshots,
		int   interval,
	    const char* particlefile,
		const char* tunnelfile, // tunneling information file including face pressure, particles deleted, etc.
		const char* cavityupdates,
		const char* contactfile,
		const char* progressfile,
		const char* debugfile
		){

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

    std::ofstream tunnelinginf(tunnelfile);
    if(!tunnelinginf) { cout << "stream error!" << endl; exit(-1);}
    tunnelinginf.setf(std::ios::scientific, std::ios::floatfield);
    tunnelinginf.precision(OPREC);
    tunnelinginf << setw(OWID) << "itr"
	        << setw(OWID) << "cur_itr"
	        << setw(OWID) << "tol_tulgth"
	        << setw(OWID) << "cur_tulgth"
	        << setw(OWID) << "time"
	        << setw(OWID) << "cur_tMass"
	        << setw(OWID) << "cur_Volume"
	        << setw(OWID) << "tol_Mass"
	        << setw(OWID) << "tol_Volume"
	        << setw(OWID) << "FacePressure"
	        << setw(OWID) << "tol_deleted"
	        << setw(OWID) << "del_times" << endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout << "stream error!" << endl; exit(-1); }
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);



    // pre_2. create particles and boundaries from existing files.
    readSample(samplefile); // create container and particles, velocity and omga are set zero.
    readBoundary(bdryfile);   // create boundaries.
	readCavityBoundary(cavitybdry); // create cylindrical boundary file

    // pre_3: define variables used in iterations.
    REAL l13, l24, l56, H, r;
    REAL l=0; // excavation length varies every excavation step, not time step, ideal condition
    REAL L=0; // actual excavation length of every excavation step
    REAL L0=0; // current total excavation length
    REAL L1=exspeed*TIMESTEP; // excavation length every time step, distinguished from L
    REAL AvePres=0; // average supporting pressure
    REAL TolPres=0; // total supporting pressure
	REAL TotalExMass=0; // total mass of particles excavated
	REAL TotalExVolume=0; // total volume of particles excavated
    REAL avgNormal=0;
    REAL avgTangt=0;
    int  stepsnum=0;
    char stepsstr[4];
    char stepsfp[50];
    REAL void_ratio=0;
    REAL void_ratio_ini=0;
    REAL bdry_penetr[7];
    int delnum=0;  // particle deletion times
    int prtlnum=0; // deleted particles number
    int  bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    int         tunnel[1]={6};    // boundary 6, cylindrical boundary
    int         face[1]={7};    // boundary 7, tunnel face
    UPDATECTL   tunnelctl[1];
    UPDATECTL   facectl[1];


    // iterations starting ...
    g_iteration=0;
    int c_iteration=0; //record steps from last particle deletion
    gettimeofday(&time_w1,NULL);

    do
    {
	// 1. create possible boundary particles and contacts between particles.
    findContact();
    findParticleOnBoundary();
	findParticleOnCavity();
	//findParticleOnTunnelFace();

	// 2. set particles' forces/moments as zero before each re-calculation,
	clearForce();

	// 3. calculate contact forces/moments and apply them to particles.
	internalForce(avgNormal, avgTangt);

	// 4. calculate boundary forces/moments and apply them to particles.
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	cavityBoundaryForce();

	// 5. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();

	// 6. calculate specimen initial void ratio.
	l56=getTopFreeParticlePosition().getz() - getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();
	H=Height;
	r=Radius;
	Volume=l13*l24*l56-PI*r*r*H;
	if (g_iteration==1){
		void_ratio_ini=Volume/getParticleVolume()-1;
	}


	//displacement control, advancement length
	tunnelctl[0].tran=vec(-TIMESTEP*exspeed/2,0,0);
	facectl[0].tran=vec(-TIMESTEP*exspeed,0,0);


	// 7. delete particles in contact with tunnel face
	REAL exMass=0; // mass of particles excavated every time
	REAL exVolume=0; // volume of particles excavated every time
	TolPres=0;
	//c_iteration=g_iteration-c_iteration;
	if (fabs(L-L0)<L1){
		L0=0;
		c_iteration=0;
		delnum +=1;
	    std::vector<RGDBDRY*>::const_iterator it;
		for (it = FaceVec.begin(); it != FaceVec.end(); ++it){
			    std::cout << (*it)->getBdryID();
				TolPres = (*it)->getNormalForce().getx();
				deleteParticleOnTunnelFace();
				exMass = (*it)->getDelMass();
				exVolume = (*it)->getDelVolume();
				prtlnum = (*it)->getCntnum();
				l=(void_ratio_ini+1)*exVolume/(PI*r*r); // l is determined basing on current average void ratio
				TotalExMass += exMass;
				TotalExVolume += exVolume;
			}
	    tunnelinginf << setw(OWID) << g_iteration
		        << setw(OWID) << c_iteration
		        << setw(OWID) << exspeed*TIMESTEP*g_iteration
		        << setw(OWID) << l
		        << setw(OWID) << timediffsec(time_w1,time_w2)
		        << setw(OWID) << exMass
		        << setw(OWID) << exVolume
		        << setw(OWID) << TotalExMass
		        << setw(OWID) << TotalExVolume
		        << setw(OWID) << TolPres
		        << setw(OWID) << prtlnum
		        << setw(OWID) << delnum << endl;
	    TotalNum -= prtlnum; // count particles remaining in the container
	}
	L=l*Cof;
	c_iteration +=1;
	L0=exspeed*TIMESTEP*c_iteration; // current total excavation length since last particle deletion

	// 8. update boundaries' position and orientation
	H = H-TIMESTEP*exspeed;
	updateCB(tunnel,tunnelctl,1);
	updateCB(face,facectl,1);
	updateCB2();


	// 9. (1) output particles and contacts information as snapshots.
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum);
	    strcpy(stepsfp, particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum);
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum);
	    strcpy(stepsfp, cavityupdates); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printCavityBoundary(stepsfp);

	    time(&timeStamp);
	    g_timeinf << setw(4) << stepsnum << " " << ctime(&timeStamp) << flush;
	    ++stepsnum;
	}

	// 9. (2) output stress and strain info.
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
    tunnelinginf.close();


}

// update CoefOfLimits[0] for cavity boundary --by cupescathy
void assembly::updateCB(int bn[], UPDATECTL cbctl[], int num){
    for(int i=0;i<num;i++){
	for(std::vector<RGDBDRY*>::iterator rt=CavityRBVec.begin();rt!=CavityRBVec.end();++rt){
	    if((*rt)->getBdryID()==bn[i]){
		(*rt)->update(cbctl[i]);
		break;
	    }
	}
    }
}

// update CoefOfLimits[1,..] for all cavity boundaries, one cylindrical and one plane boundary in this case --by cupescathy
void assembly::updateCB2(){
    for(std::vector<RGDBDRY*>::iterator rt=CavityRBVec.begin();rt!=CavityRBVec.end();++rt){
	if((*rt)->getBdryID()==6){
	    for(std::vector<RGDBDRY*>::iterator lt=CavityRBVec.begin();lt!=CavityRBVec.end();++lt){
		if((*lt)->getBdryID()==7)
		    (*rt)->CoefOfLimits[1].apt=(*lt)->CoefOfLimits[0].apt;
	    }
	}
	else if((*rt)->getBdryID()==7){
	    for(std::vector<RGDBDRY*>::iterator lt=CavityRBVec.begin();lt!=CavityRBVec.end();++lt){
		if((*lt)->getBdryID()==6)
		    (*rt)->CoefOfLimits[1].apt=(*lt)->CoefOfLimits[0].apt;
	    }
	  }
    }
}

void assembly::deleteParticleOnTunnelFace(){  // delete particle on tunnel face which is part of cavity boundaries --by cupescathy
    std::vector<RGDBDRY*>::iterator rt;
    for(rt = FaceVec.begin(); rt != FaceVec.end(); ++rt)
	(*rt)->deleteParticleOnBoundary(ParticleVec);
}

// DE part for DE-FE coupling

void assembly::ReadGhostSample(char* str,
			  double& TimeStep,
			  int& NumStep,
			  Tahoe::iArray2DT& fGhostElemSet,
			  std::vector<particle*>& fGhostParticleList){
    ifstream ifs(str);
    if(!ifs) {
	cout<<"stream error in createSample!"<<endl; exit(-1);
    }

    char s[50];
    int Major, Minor;

    ifs >> s

	>> s >> s >> NUM_STEP
	>> s >> s >> TIMESTEP
	>> s >> s >> MASS_SCL
	>> s >> s >> MNT_SCL
	>> s >> s >> GRVT_SCL
	>> s >> s >> DMP_F
	>> s >> s >> DMP_M

	>> s >> s >> DMP_CNT
	>> s >> s >> FRICTION
	>> s >> s >> BDRYFRIC
	>> s >> s >> COHESION

	>> s >> s >> COMPRESS_RATE
	>> s >> s >> RELEASE_RATE
	>> s >> s >> PILE_RATE
	>> s >> s >> STRESS_ERROR

	>> s >> s >> idum

	>> s
	>> TotalNum
	>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
	>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
	>>s>>s;

    RORC = 1;

    ParticleVec.clear();
    int ID, type, elemgrp, elemnum;
    REAL a, b, c, px,py,pz,dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz;
    REAL vx,vy,vz,omx,omy,omz,fx,fy,fz,mx,my,mz;

    for (int i = 0; i < TotalNum; i++)
    {
	ifs >> ID >> type
	    >> a >> b>> c >> px >> py >> pz >> dax >> day >> daz >> dbx >> dby >> dbz >> dcx >> dcy >> dcz
	    >> vx >> vy >> vz >> omx >> omy >> omz >> fx >> fy >> fz >> mx >> my >> mz;

	particle* pt = NULL;
	if (type == 10) {
	    ifs >> elemgrp >> elemnum;
	    pt = new Tahoe::GhostParticleT(ID,type,vec(a,b,c),vec(px,py,pz),vec(dax,day,daz),vec(dbx,dby,dbz),vec(dcx,dcy,dcz), elemgrp,elemnum-1);
	    fGhostParticleList.push_back(pt);
	}
	else
		pt = new particle(ID,type,vec(a,b,c),vec(px,py,pz),vec(dax,day,daz),vec(dbx,dby,dbz),vec(dcx,dcy,dcz));

//      optional settings for a particle's initial status
//	pt->setPrevVelocity(vec(vx,vy,vz));
//	pt->setCurrVelocity(vec(vx,vy,vz));
//	pt->setPrevOmga(vec(omx,omy,omz));
//	pt->setCurrOmga(vec(omx,omy,omz));
//	pt->setConstForce(vec(fx,fy,fz));  // constant force, not initial force
//	pt->setConstMoment(vec(mx,my,mz)); // constant moment, not initial moment

	ParticleVec.push_back(pt);
    }

    TimeStep = TIMESTEP;
    NumStep  = NUM_STEP;

    ifs >> s
	>> Major
	>> Minor;

    fGhostElemSet.Dimension(Major, Minor);
    if (Major > 0) fGhostElemSet.ReadNumbered(ifs);

    for (int i = 0; i < fGhostElemSet.MajorDim(); i++)
	fGhostElemSet(i,1) -= 1;

    ifs.close();
}

void assembly::Run(int total_steps, int PrintNum)
{
    // pre_1. define variables used in iterations.
    REAL avgNormal=0;
    REAL avgTangt=0;

    // iterations starting ...
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles.
	if (g_iteration%UPDATE_CNT==0)
	    findContact();

	// 2. set particles' forces/moments as zero before each re-calculation,
	clearForce(); // apply constant force in 1 step
	//setForceZero(PrintNum); // apply constant force in 100 steps

	// 3. calculate contact forces/moments and apply them to particles.
	internalForce(avgNormal, avgTangt);

	// 4. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();

	// 5. calculate specimen void ratio.

	// 6. loop break conditions.

    } while (++g_iteration < total_steps);

}


} // namespace dem ends

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
