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
  ofstream ofs;
  char filein[50];
  char fileout[50];
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

  for(int step=first; step<=last; step+=incre) {
    if(argc == 2) {
      strcpy(filein, argv[1]);
    } else {
      sprintf(num, "%03d", step);
      strcpy(filein, argv[1]);
      strcat(filein, "_");
      strcat(filein, num);
    }
    strcpy(fileout, filein);
    strcat(filein, ".dat");
    strcat(fileout, "_tetra.dat");      

    ifs.open(filein);
    if(!ifs)  { cout<<"stream error!"<<endl; exit(-1);}

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
    }
    else {
      cout << "generating file " << fileout << " ......" <<endl;
      ofs.open(fileout);
      if(!ofs)  { cout<<"stream error!"<<endl; exit(-1);}
      ofs.setf(ios::scientific, ios::floatfield);

      // prepare data for Qhull
      std::stringstream coordStream;
      coordStream << 3 << " " << totalContact << " ";
      for(int k = 0; k < totalContact; ++k) {
	ifs >> x >> y >> z
	    >> normal_x >> normal_y >> normal_z
	    >> tangt_x >> tangt_y >> tangt_z
	    >> total_x >> total_y >> total_z
	    >> normal >> shear >> total >> penetration;
	coordStream << x << " "
		    << y << " "
		    << z << " ";     
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
      ofs << "VARIABLES=" << endl
	  << setw(OWID) << "x"
	  << setw(OWID) << "y"
	  << setw(OWID) << "z"
	  << endl;
      ofs << "ZONE N=" << totalContact << ", E=" << totalTetra  <<", DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON" << endl;

      // open again
      ifs.open(filein);
      if(!ifs)  { cout<<"stream error!"<<endl; exit(-1);}
      ifs >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
	  >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s;
      for(int k = 0; k < totalContact; ++k) {
	ifs >> x >> y >> z
	    >> normal_x >> normal_y >> normal_z
	    >> tangt_x >> tangt_y >> tangt_z
	    >> total_x >> total_y >> total_z
	    >> normal >> shear >> total >> penetration;

	ofs << setw(OWID) << x
	    << setw(OWID) << y
	    << setw(OWID) << z << std::endl;

      }
      ifs.close();

      int m, n, i, j;
      for(int it=0; it!=totalTetra; it++) {
	connectStream >> m >> n >> i >> j;
	++m;
	++n;
	++i;
	++j; // Qhull IDs start from 0; FEM IDs start from 1
	ofs << setw(OWID) << m << setw(OWID) << n << setw(OWID) << i << setw(OWID) << j << std::endl;
      }

      ofs.close();
    } // end of if (totalContact >= 5)
  }
  
  return 0;
}
