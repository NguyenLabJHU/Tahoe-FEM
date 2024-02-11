#include "Vec.h"
#include "Tetra.h"
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/Qhull.h"

using namespace std;

const double PI = 3.1415927;
const double EPS = 1.0E-12;
static const int OWID = 15;
static const int OPREC= 6;

int main(int argc, char *argv[])
{
  if(argc < 2) {
    cout << endl 
	 << "-- Plot particle centroid tetrahedrons --"<<endl
	 << "Usage:" << endl
	 << "1) process a single file:  plot_tetra_particle particle_file" << endl
	 << "   --example: plot_tetra_particle triaxial_particle_008" << endl
	 << "2) process multiple files: plot_tetra_particle particle_file_prefix first_suffix last_suffix suffix_increment" << endl
	 << "   --example: plot_tetra_particle triaxial_particle 1 100 5" << endl << endl;
    return -1;
  }	

  int first, last, incre;
  if(argc == 2) {
    first = 0;
    last  = 1;
    incre = 2;
  }
  else {
    first = atoi(argv[2]);
    last  = atoi(argv[3]);
    incre = atoi(argv[4]);
  }

  ifstream ifs;
  ofstream ofs, ofs2;
  char filein[50];
  char fileout[50], fileout2[50];
  char num[4], s[20], snum[20];

  int id, type;       
  double radius_a,       radius_b,       radius_c;     
  double position_x,     position_y,     position_z;       
  double axle_a_x,       axle_a_y,       axle_a_z;      
  double axle_b_x,       axle_b_y,       axle_b_z;      
  double axle_c_x,       axle_c_y,       axle_c_z;    
  double velocity_x,     velocity_y,     velocity_z;         
  double omga_x,         omga_y,         omga_z;       
  double force_x,        force_y,        force_z;       
  double moment_x,       moment_y,       moment_z;

  ofstream ofs3;
  ofs3.open("tetra_particle_stats");
  if(!ofs3)  { cout<<"stream error 4!"<<endl; exit(-1);}
  ofs3.setf(ios::scientific, ios::floatfield);
  ofs3 << setw(OWID) << "snap"
       << setw(OWID) << "avgLength"
       << setw(OWID) << "avgSideLen"
       << setw(OWID) << "avgBottLen"
       << setw(OWID) << "sideBottRatio"
       << setw(OWID) << "topSolidAng_sr"
       << setw(OWID) << "topConeAng_deg"
       << setw(OWID) << "avgVolume"
       << std::endl;

  for(int step=first; step<=last; step+=incre) {
    if(argc == 2) {
      strcpy(filein, argv[1]);
    } else {
      sprintf(num, "%03d", step);
      strcpy(filein, argv[1]);
      strcat(filein, "_");
      strcat(filein, num);
    }
    strcpy(fileout, "tetra_block_");
    strcat(fileout, filein);
    strcat(fileout, ".dat");

    strcpy(fileout2, "tetra_point_");
    strcat(fileout2, filein);
    strcat(fileout2, ".dat");

    cout << "generating file " << fileout << " " << fileout2 << " ......" <<endl;

    ifs.open(filein);
    if(!ifs)  { cout<<"stream error 1!"<<endl; exit(-1);}
    ofs.open(fileout);
    if(!ofs)  { cout<<"stream error 2!"<<endl; exit(-1);}
    ofs.setf(ios::scientific, ios::floatfield);
    ofs2.open(fileout2);
    if(!ofs2)  { cout<<"stream error 3!"<<endl; exit(-1);}
    ofs2.setf(ios::scientific, ios::floatfield);

    int totalParticle;
    ifs >> totalParticle;
    ifs >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
        >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
        >>s>>s>>s>>s>>s>>s>>s>>s>>s;

    // prepare data for Qhull
    std::stringstream coordStream;
    coordStream << 3 << " " << totalParticle << " ";
    std::vector<Vec> spaceCoords(totalParticle);
    std::size_t index = 0;
    for(int k = 0; k < totalParticle; ++k) {

      ifs >> id >> type
	  >> radius_a   >>  radius_b   >>  radius_c     
	  >> position_x >>  position_y >>  position_z       
	  >> axle_a_x   >>  axle_a_y   >>  axle_a_z       
	  >> axle_b_x   >>  axle_b_y   >>  axle_b_z       
	  >> axle_c_x   >>  axle_c_y   >>  axle_c_z     
	  >> velocity_x >>  velocity_y >>  velocity_z         
	  >> omga_x     >>  omga_y     >>  omga_z        
	  >> force_x    >>  force_y    >>  force_z       
	  >> moment_x   >>  moment_y   >>  moment_z;

      coordStream << position_x << " "
		  << position_y << " "
		  << position_z << " ";

      // store coordindates
      spaceCoords[index++] = Vec(position_x, position_y, position_z);     
    }
    ifs.close();

    // run Qhull and output
    orgQhull::RboxPoints rbox;
    rbox.appendPoints(coordStream); 
    orgQhull::Qhull qhull;
    qhull.runQhull(rbox, "d Qbb Qt i"); // "d Qbb Qt i"; // "d Qt Qz Qs i"
    std::stringstream connectStream;
    qhull.setOutputStream(&connectStream);
    qhull.outputQhull();

    int totalTetra;
    connectStream >> totalTetra; // totalTetra generated by Qhull, may include coplanar tetrahedrons
    std::vector<Tetra> tetraVec;
    int i, j, k, l;
    for(int it=0; it!=totalTetra; it++) {
      connectStream >> i >> j >> k >> l;
      ++i;
      ++j;
      ++k;
      ++l; // Qhull IDs start from 0; FEM IDs start from 1
      Tetra tmpTetra(i, j, k, l, spaceCoords[i-1], spaceCoords[j-1], spaceCoords[k-1], spaceCoords[l-1]);
      if (fabs(tmpTetra.getVolume()) > EPS) // only store non-coplanar tetrahedrons
	tetraVec.push_back(tmpTetra);
    }
    for (int i = 0; i < tetraVec.size(); ++i)
      tetraVec[i].setNodeOrder();
    
    // print total average of each variable
    double avgLength = 0;
    double avgSideLen = 0;
    double avgBottLen = 0;
    double sideBottRatio = 0;
    double topSolidAng_sr = 0;
    double avgVolume = 0;
    for (int i = 0; i < tetraVec.size(); ++i) {
      avgLength += tetraVec[i].getInfo()[0];
      avgSideLen += tetraVec[i].getInfo()[1];
      avgBottLen += tetraVec[i].getInfo()[2];
      sideBottRatio += tetraVec[i].getInfo()[3];
      topSolidAng_sr += tetraVec[i].getInfo()[4];
      avgVolume += tetraVec[i].getInfo()[6];
    }
    avgLength /= tetraVec.size();
    avgSideLen /= tetraVec.size();
    avgBottLen /= tetraVec.size();
    sideBottRatio /= tetraVec.size();
    topSolidAng_sr /= tetraVec.size();
    avgVolume /= tetraVec.size();
    ofs3 << setw(OWID) << step
         << setw(OWID) << avgLength
	 << setw(OWID) << avgSideLen
	 << setw(OWID) << avgBottLen
	 << setw(OWID) << sideBottRatio
	 << setw(OWID) << topSolidAng_sr
	 << setw(OWID) << acos(1 - topSolidAng_sr/(2*PI)) *2 *180/PI
	 << setw(OWID) << avgVolume
	 << std::endl;

    // print Tecplot Block format
    ofs << "VARIABLES=" << endl
	<< setw(OWID) << "x"
	<< setw(OWID) << "y"
	<< setw(OWID) << "z"
	<< setw(OWID) << "length"
	<< setw(OWID) << "sideLen"
	<< setw(OWID) << "bottLen"
	<< setw(OWID) << "sideBottRatio"
	<< setw(OWID) << "topSolidAng_sr"
	<< setw(OWID) << "topConeAng_deg"
	<< setw(OWID) << "volume"
	<< endl;
    ofs << "ZONE N=" << totalParticle << ", E=" << tetraVec.size()  <<", DATAPACKING=BLOCK, \
VARLOCATION=([4,5,6,7,8,9,10]=CELLCENTERED), ZONETYPE=FETETRAHEDRON, STRANDID=2, SOLUTIONTIME=" << step << endl;

    // Tecplot: 
    // BLOCK format must be used for cell-centered data.
    // For nodal variables, provide the values for each variable in nodal order. 
    // Similarly, for cell-centered values, provide the variable values in cell order.

    // Implement A: for current Tecplot line character limit 32,000
    int lineLen = 32000;
    int valNum = lineLen / OWID - 10; // leave room for 10 values
    int kut;

    // print variables
    // coordinates
    kut = 0;
    for (std::size_t i = 0; i < spaceCoords.size(); ++i) {
      ofs << std::setw(OWID) << spaceCoords[i].getX();
      ++kut; if (kut >= valNum) {ofs << std::endl; kut = 0;}
    }
    ofs << std::endl;

    kut = 0;
    for (std::size_t i = 0; i < spaceCoords.size(); ++i) {
      ofs << std::setw(OWID) << spaceCoords[i].getY();
      ++kut; if (kut >= valNum) {ofs << std::endl; kut = 0;}
    }
    ofs << std::endl;

    kut = 0;
    for (std::size_t i = 0; i < spaceCoords.size(); ++i) {
      ofs << std::setw(OWID) << spaceCoords[i].getZ();
      ++kut; if (kut >= valNum) {ofs << std::endl; kut = 0;}
    }
    ofs << std::endl;

    // variables
    kut = 0;
    for (std::size_t i = 0; i < tetraVec.size(); ++i) {
      ofs << std::setw(OWID) << tetraVec[i].getInfo()[0];
      ++kut; if (kut >= valNum) {ofs << std::endl; kut = 0;}
    }
    ofs << std::endl;

    kut = 0;
    for (std::size_t i = 0; i < tetraVec.size(); ++i) {
      ofs << std::setw(OWID) << tetraVec[i].getInfo()[1];
      ++kut; if (kut >= valNum) {ofs << std::endl; kut = 0;}
    }
    ofs << std::endl;

    kut = 0;
    for (std::size_t i = 0; i < tetraVec.size(); ++i) {
      ofs << std::setw(OWID) << tetraVec[i].getInfo()[2];
      ++kut; if (kut >= valNum) {ofs << std::endl; kut = 0;}
    }
    ofs << std::endl;

    kut = 0;
    for (std::size_t i = 0; i < tetraVec.size(); ++i) {
      ofs << std::setw(OWID) << tetraVec[i].getInfo()[3];
      ++kut; if (kut >= valNum) {ofs << std::endl; kut = 0;}
    }
    ofs << std::endl;

    kut = 0;
    for (std::size_t i = 0; i < tetraVec.size(); ++i) {
      ofs << std::setw(OWID) << tetraVec[i].getInfo()[4];
      ++kut; if (kut >= valNum) {ofs << std::endl; kut = 0;}
    }
    ofs << std::endl;

    kut = 0;
    for (std::size_t i = 0; i < tetraVec.size(); ++i) {
      ofs << std::setw(OWID) << tetraVec[i].getInfo()[5];
      ++kut; if (kut >= valNum) {ofs << std::endl; kut = 0;}
    }
    ofs << std::endl;

    kut = 0;
    for (std::size_t i = 0; i < tetraVec.size(); ++i) {
      ofs << std::setw(OWID) << tetraVec[i].getInfo()[6];
      ++kut; if (kut >= valNum) {ofs << std::endl; kut = 0;}
    }
    ofs << std::endl;

    // print connections
    for (std::size_t i = 0; i < tetraVec.size(); ++i)
      ofs << std::setw(OWID) << tetraVec[i].getI()
	  << std::setw(OWID) << tetraVec[i].getJ()
	  << std::setw(OWID) << tetraVec[i].getK()
	  << std::setw(OWID) << tetraVec[i].getL()
	  << std::endl;
    ofs.close();

    // print Tecplot POINT format
    ofs2<< "VARIABLES=" << endl
	<< setw(OWID) << "x"
	<< setw(OWID) << "y"
	<< setw(OWID) << "z"
	<< endl;
    ofs2<< "ZONE N=" << totalParticle << ", E=" << tetraVec.size()  <<", DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON, STRANDID=2, SOLUTIONTIME=" << step << endl;
    for (int i = 0; i < spaceCoords.size(); ++i)
      ofs2 << std::setw(OWID) << spaceCoords[i].getX() << std::setw(OWID) << spaceCoords[i].getY() << std::setw(OWID) << spaceCoords[i].getZ() << std::endl;
    for (std::size_t i = 0; i < tetraVec.size(); ++i)
      ofs2<< std::setw(OWID) << tetraVec[i].getI()
	  << std::setw(OWID) << tetraVec[i].getJ()
	  << std::setw(OWID) << tetraVec[i].getK()
	  << std::setw(OWID) << tetraVec[i].getL()
	  << std::endl;
    ofs2.close();

  } // end of for(int step=first; step<=last; step+=incre)

  ofs3.close();
  return 0;
}
