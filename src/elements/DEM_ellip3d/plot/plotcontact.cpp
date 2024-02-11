#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <cmath>

using namespace std;

const double PI = 3.1415927;
static const int OWID = 15;

int main(int argc, char *argv[])
{
  if(argc < 2) {
    cout << endl 
	 << "-- Plot contact force --"<<endl
	 << "Usage:" << endl
	 << "1) process a single file:  plotcontact contact_file" << endl
	 << "   --example: plotcontact triaxial_contact_008" << endl
	 << "2) process multiple files: plotcontact contact_file_prefix  first_suffix  last_suffix  suffix_increment" << endl
	 << "   --example: plotcontact triaxial_contact  1  100  5" << endl << endl;
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
  char num[4], s[20];

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
  double contact_x;
  double contact_y;
  double contact_z;
  double normal_x;
  double normal_y;
  double normal_z;
  double tangt_x;
  double tangt_y;
  double tangt_z;
  double vibra_t_step;
  double impact_t_step;
  double dir_x;
  double dir_y;
  double dir_z;

  for(int snapshot = first; snapshot <= last; snapshot += incre) {
    if(argc == 2)
      strcpy(filein, argv[1]);
    else {
      sprintf(num, "%03d", snapshot);
      strcpy(filein, argv[1]);
      strcat(filein, "_");
      strcat(filein, num);
    }
      
    strcpy(fileout, filein);
    strcat(fileout, ".dat");
    cout << "generating file " << fileout << " ......" <<endl;

    ifs.open(filein);
    if(!ifs)  { cout<<"stream error!"<<endl; exit(-1);}
    ofs.open(fileout);
    if(!ofs)  { cout<<"stream error!"<<endl; exit(-1);}
    ofs.setf(ios::scientific, ios::floatfield);

    int totalNum;
    ifs >> totalNum;
    ifs>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
       >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
       >>s>>s>>s; // add dir_x, dir_y, dir_z

    ofs << "VARIABLES=" << endl
	<< setw(OWID) << "x"
	<< setw(OWID) << "y"
	<< setw(OWID) << "z"
	<< setw(OWID) << "normal_x"
	<< setw(OWID) << "normal_y"
	<< setw(OWID) << "normal_z"
	<< setw(OWID) << "tangt_x"
	<< setw(OWID) << "tangt_y"
	<< setw(OWID) << "tangt_z"
	<< setw(OWID) << "total_x"
	<< setw(OWID) << "total_y"
	<< setw(OWID) << "total_z"
	<< setw(OWID) << "normal"
	<< setw(OWID) << "shear"
	<< setw(OWID) << "total"
	<< setw(OWID) << "penetration"
	<< std::setw(OWID) << "dir_x"
	<< std::setw(OWID) << "dir_y"
	<< std::setw(OWID) << "dir_z"
	<< endl;

    ofs << "ZONE I=" << totalNum <<", DATAPACKING=POINT, STRANDID=1, SOLUTIONTIME=" << snapshot << endl;
    for(int it = 0; it < totalNum; ++it) {
      ifs >> ptcl_1 >> ptcl_2 >> point1_x >> point1_y >> point1_z >> point2_x >> point2_y >> point2_z 
	  >> radius_1 >> radius_2 >> penetration >> tangt_disp >> contact_radius >> R0 >> E0 >> normal_force >> tangt_force >> normal_all
	  >> contact_x >> contact_y >> contact_z >> normal_x >> normal_y >> normal_z >> tangt_x >> tangt_y >> tangt_z >> vibra_t_step >> impact_t_step >> dir_x >> dir_y >> dir_z;
  
      ofs << setw(OWID) << contact_x
	  << setw(OWID) << contact_y
	  << setw(OWID) << contact_z
	  << setw(OWID) << normal_x
	  << setw(OWID) << normal_y
	  << setw(OWID) << normal_z
	  << setw(OWID) << tangt_x
	  << setw(OWID) << tangt_y
	  << setw(OWID) << tangt_z
	  << setw(OWID) << normal_x + tangt_x
	  << setw(OWID) << normal_y + tangt_y
	  << setw(OWID) << normal_z + tangt_z

	  << setw(OWID) << sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z)
	  << setw(OWID) << sqrt(tangt_x*tangt_x   + tangt_y*tangt_y   + tangt_z*tangt_z)
	  << setw(OWID) << sqrt( (normal_x+tangt_x)*(normal_x+tangt_x) + (normal_y+tangt_y)*(normal_y+tangt_y) + (normal_z+tangt_z)*(normal_z+tangt_z) )

	  << setw(OWID) << penetration
	  << std::setw(OWID) << dir_x
	  << std::setw(OWID) << dir_y
	  << std::setw(OWID) << dir_z
	  << endl;   
    }
	
    ifs.close();
    ofs.close();

  }

  return 0;
}
