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
	 << "-- Plot particle centroid meshes --"<<endl
	 << "Usage:" << endl
	 << "1) process a single file:  plot_tetra_particle particle_file" << endl
	 << "   --example: plot_tetra_particle triaxial_particle_008" << endl
	 << "2) process multiple files: plot_tetra_particle contact_file_prefix  first_suffix  last_suffix  suffix_increment" << endl
	 << "   --example: plot_tetra_particle triaxial_particle  1  100  5" << endl << endl;
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

  for(int step=first; step<=last; step+=incre) {
    if(argc == 2) {
      strcpy(filein, argv[1]);
      strcpy(fileout, filein);
    } else {
      sprintf(num, "%03d", step);
      strcpy(filein, argv[1]);
      strcat(filein, "_");
      strcat(filein, num);
    }
    strcpy(fileout, filein);
    strcat(fileout, "_tetra.dat");      
    cout << "generating file " << fileout << " ......" <<endl;

    ifs.open(filein);
    if(!ifs)  { cout<<"stream error 1!"<<endl; exit(-1);}
    ofs.open(fileout);
    if(!ofs)  { cout<<"stream error 2!"<<endl; exit(-1);}
    ofs.setf(ios::scientific, ios::floatfield);

    int totalParticle;
    ifs >> totalParticle;
    ifs >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
        >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
        >>s>>s>>s>>s>>s>>s>>s>>s>>s;

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

    // prepare data for Qhull
    std::stringstream coordStream;
    coordStream << 3 << " " << totalParticle << " ";
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
    ofs << "ZONE N=" << totalParticle << ", E=" << totalTetra  <<", DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON" << endl;

    // open again
    ifs.open(filein);
    if(!ifs)  { cout<<"stream error 3!"<<endl; exit(-1);}
    ifs >> totalParticle;
    ifs >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
        >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
        >>s>>s>>s>>s>>s>>s>>s>>s>>s;

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

      ofs << setw(OWID) << position_x
	  << setw(OWID) << position_y
	  << setw(OWID) << position_z << std::endl;

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
  }

  return 0;
}
