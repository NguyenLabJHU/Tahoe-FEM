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
	 << "-- Plot contact tetrahedrons using data processed by plotcontact --"<<endl
	 << "Usage:" << endl
	 << "1) process a single file:  plot_tetra_contact contact_file" << endl
	 << "   --example: plot_tetra_contact triaxial_contact_008 (not triaxial_contact_008.dat)" << endl
	 << "2) process multiple files: plot_tetra_contact contact_file_prefix  first_suffix  last_suffix  suffix_increment" << endl
	 << "   --example: plot_tetra_contact triaxial_contact  1  100  5" << endl << endl;
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

  int ptcl_1;
  int ptcl_2;
  double point1_x;
  double point1_y;
  double point1_z;
  double point2_x;
  double point2_y;
  double point2_z;
  double radius_1;
  double radius_2;
  double penetration;
  double tangt_disp;
  double contact_radius;
  double R0;
  double E0;
  double normal_force;
  double tangt_force;
  double normal_all;
  double x;
  double y;
  double z;
  double normal_x;
  double normal_y;
  double normal_z;
  double tangt_x;
  double tangt_y;
  double tangt_z;
  double total_x;
  double total_y;
  double total_z;
  double normal;
  double shear;
  double total;
  double vibra_t_step;
  double impact_t_step;

  ofstream ofs3;
  ofs3.open("tetra_contact_stats");
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

    strcat(filein, ".dat");

    ifs.open(filein);
    if(!ifs)  { cout<<"stream error 1!"<<endl; exit(-1);}

    // read data and obtain lines, snum reads "I=1216252," from "ZONE I=1216252, DATAPACKING=POINT"
    ifs >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
        >>s>>s>>s>>s>>s>>s>>s>>s>>snum>>s;
    std::string strall(snum);
    strall.erase(strall.end()-1);
    strall.erase(0,2);
    //std::cout << strall << std::endl;
    int totalContact = atol(strall.c_str());
    //std::cout << totalContact << std::endl;

    if (totalContact < 5) { // Qhull needs at least 5 points in 3D
      ifs.close();
    } else {
      cout << "generating file " << fileout << " " << fileout2 << " ......" <<endl;
      ofs.open(fileout);
      if(!ofs)  { cout<<"stream error 2!"<<endl; exit(-1);}
      ofs.setf(ios::scientific, ios::floatfield);
      ofs2.open(fileout2);
      if(!ofs2)  { cout<<"stream error 3!"<<endl; exit(-1);}
      ofs2.setf(ios::scientific, ios::floatfield);

      // prepare data for Qhull
      std::stringstream coordStream;
      coordStream << 3 << " " << totalContact << " ";
      std::vector<Vec> spaceCoords(totalContact);
      std::size_t index = 0;
      for(int k = 0; k < totalContact; ++k) {
	ifs >> x >> y >> z
	    >> normal_x >> normal_y >> normal_z
	    >> tangt_x >> tangt_y >> tangt_z
	    >> total_x >> total_y >> total_z
	    >> normal >> shear >> total >> penetration;

	coordStream << x << " "
		    << y << " "
		    << z << " ";     

	// store coordindates
	spaceCoords[index++] = Vec(x, y, z); 
      }
      ifs.close();

      // run Qhull and output
      orgQhull::RboxPoints rbox;
      rbox.appendPoints(coordStream); 
      orgQhull::Qhull qhull;
      qhull.runQhull(rbox, "d Qbb Qt i"); // "d Qt Qz Qs i"
      std::stringstream connectStream;
      qhull.setOutputStream(&connectStream);
      qhull.outputQhull();

      int totalTetra;
      connectStream >> totalTetra;
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
      double topConeAng_deg = 0;
      double avgVolume = 0;
      for (int i = 0; i < tetraVec.size(); ++i) {
	avgLength += tetraVec[i].getInfo()[0];
	avgSideLen += tetraVec[i].getInfo()[1];
	avgBottLen += tetraVec[i].getInfo()[2];
	sideBottRatio += tetraVec[i].getInfo()[3];
	topSolidAng_sr += tetraVec[i].getInfo()[4];
	topConeAng_deg += tetraVec[i].getInfo()[5];
	avgVolume += tetraVec[i].getInfo()[6];
      }
      avgLength /= tetraVec.size();
      avgSideLen /= tetraVec.size();
      avgBottLen /= tetraVec.size();
      sideBottRatio /= tetraVec.size();
      topSolidAng_sr /= tetraVec.size();
      topConeAng_deg /= tetraVec.size();
      avgVolume /= tetraVec.size();
      ofs3 << setw(OWID) << step
	   << setw(OWID) << avgLength
	   << setw(OWID) << avgSideLen
	   << setw(OWID) << avgBottLen
	   << setw(OWID) << sideBottRatio
	   << setw(OWID) << topSolidAng_sr
	   << setw(OWID) << topConeAng_deg
	   << setw(OWID) << avgVolume
	   << std::endl;

      // print Tecplot Block format
      ofs << "VARIABLES=" << endl
	  << setw(OWID) << "x"
	  << setw(OWID) << "y"
	  << setw(OWID) << "z"
	  << setw(OWID) << "avgLength"
	  << setw(OWID) << "avgSideLen"
	  << setw(OWID) << "avgBottLen"
	  << setw(OWID) << "sideBottRatio"
	  << setw(OWID) << "topSolidAng_sr"
	  << setw(OWID) << "topConeAng_deg"
	  << setw(OWID) << "avgVolume"
	  << endl;
      ofs << "ZONE N=" << totalContact << ", E=" << tetraVec.size()  <<", DATAPACKING=BLOCK, \
VARLOCATION=([4,5,6,7,8,9,10]=CELLCENTERED), ZONETYPE=FETETRAHEDRON" << endl;

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
      ofs2<< "ZONE N=" << totalContact << ", E=" << tetraVec.size()  <<", DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON" << endl;
      for (int i = 0; i < spaceCoords.size(); ++i)
	ofs2 << std::setw(OWID) << spaceCoords[i].getX() << std::setw(OWID) << spaceCoords[i].getY() << std::setw(OWID) << spaceCoords[i].getZ() << std::endl;
      for (std::size_t i = 0; i < tetraVec.size(); ++i)
	ofs2<< std::setw(OWID) << tetraVec[i].getI()
	    << std::setw(OWID) << tetraVec[i].getJ()
	    << std::setw(OWID) << tetraVec[i].getK()
	    << std::setw(OWID) << tetraVec[i].getL()
	    << std::endl;
      ofs2.close();

    } // end of if (totalContact < 5)

  } // end of for(int step=first; step<=last; step+=incre)
  
  ofs3.close();
  return 0;
}
