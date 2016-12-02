#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <cmath>

using namespace std;
static const int OWID = 15;

int main(int argc, char *argv[])
{
  if(argc < 5) {
    cout << endl
	 << "-- Colletct Dual Particles' Contact Information --" << endl
	 << "Usage: probedual contact_file_prefix  first_suffix  last_suffix  suffix_increment" << endl
	 << "   --example: probedual dep_contact 0 100 5" << endl << endl;
    return -1;
  }	

  int first, last, incre, n, ic;
  first   = atoi(argv[2]);
  last    = atoi(argv[3]);
  incre   = atoi(argv[4]);

  ifstream ifs;
  ofstream ofs;
  char filein1[50];
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

  strcpy(fileout, argv[1]);
  strcat(fileout, "_dual");
  ofs.open(fileout);
  if(!ofs)   { cout<<"ofs stream error!"<<endl; exit(-1);}
  ofs.setf(ios::scientific, ios::floatfield);
  ofs << std::setw(OWID) << "snapshot"
      << std::setw(OWID) << "ptcl_1"
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
      << std::setw(OWID) << "normal_all_z"
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
 
  for(n=first; n<=last; n+=incre) {
    sprintf(num, "%03d", n);
    strcpy(filein1, argv[1]);
    strcat(filein1, "_");
    strcat(filein1, num);

    ifs.open(filein1);
    if(!ifs)  { cout<<"ifs stream error!"<<endl; exit(-1);}
    std::string line;

    while (getline(ifs, line)) {
      std::istringstream ss(line);
      ss >> ptcl_1
	 >> ptcl_2
	 >> point1_x
	 >> point1_y
	 >> point1_z
	 >> point2_x
	 >> point2_y
	 >> point2_z
	 >> radius_1
	 >> radius_2
	 >> penetration
	 >> tangt_disp
	 >> contact_radius
	 >> R0
	 >> E0
	 >> normal_force
	 >> tangt_force
	 >> normal_all
	 >> contact_x
	 >> contact_y
	 >> contact_z
	 >> normal_x
	 >> normal_y
	 >> normal_z
	 >> tangt_x
	 >> tangt_y
	 >> tangt_z
	 >> vibra_t_step
	 >> impact_t_step;

      ofs << std::setw(OWID) << n
	  << std::setw(OWID) << ptcl_1       
	  << std::setw(OWID) << ptcl_2       
	  << std::setw(OWID) << point1_x     
	  << std::setw(OWID) << point1_y     
	  << std::setw(OWID) << point1_z     
	  << std::setw(OWID) << point2_x     
	  << std::setw(OWID) << point2_y     
	  << std::setw(OWID) << point2_z     
	  << std::setw(OWID) << radius_1     
	  << std::setw(OWID) << radius_2     
	  << std::setw(OWID) << penetration  
	  << std::setw(OWID) << tangt_disp   
	  << std::setw(OWID) << contact_radius
	  << std::setw(OWID) << R0           
	  << std::setw(OWID) << E0           
	  << std::setw(OWID) << normal_force 
	  << std::setw(OWID) << tangt_force  
	  << std::setw(OWID) << normal_all 
	  << std::setw(OWID) << contact_x    
	  << std::setw(OWID) << contact_y    
	  << std::setw(OWID) << contact_z    
	  << std::setw(OWID) << normal_x     
	  << std::setw(OWID) << normal_y     
	  << std::setw(OWID) << normal_z     
	  << std::setw(OWID) << tangt_x      
	  << std::setw(OWID) << tangt_y      
	  << std::setw(OWID) << tangt_z      
	  << std::setw(OWID) << vibra_t_step 
	  << std::setw(OWID) << impact_t_step
	  << std::endl;
    }
      
    ifs.close();
  }
  ofs.close();
  
  return 0;
}
