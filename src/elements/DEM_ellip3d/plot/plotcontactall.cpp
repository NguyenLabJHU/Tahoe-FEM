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

using namespace std;

const double PI = 3.1415927;
static const int OWID = 15;
static const int OPREC= 6;

int main(int argc, char *argv[])
{
  if(argc == 1) {
    cout << endl 
	 << "-- Merge internal and boundary contact forces from Tecplot .dat files --"<<endl
	 << "Usage:" << endl
	 << "1) process a single snapshot:  plotcontactall contact_file boundary_file" << endl
	 << "   --example: plotcontactall triaxial_contact_008 triaxial_bdrycntc_008 (NOT triaxial_contact_008.dat triaxial_bdrycntc_008.dat)" << endl
	 << "2) process multiple snapshots: plotcontactall contact_file_prefix boundary_file_prefix  first_suffix  last_suffix  suffix_increment" << endl
	 << "   --example: plotcontactall triaxial_contact triaxial_bdrycntc  1  100  5" << endl << endl;
    return -1;
  }	

  int first, last, incre;
  if(argc == 3) {
    first = 0;
    last  = 1;
    incre = 2;
  } else if (argc == 6) {
    first = atoi(argv[3]);
    last  = atoi(argv[4]);
    incre = atoi(argv[5]);
  }

  ifstream ifs, ifs2;
  ofstream ofs;
  char filein[50], filein2[50];
  char fileout[50];
  char num[4], s[20], snum[20];

  for(int snapshot = first; snapshot <= last; snapshot += incre) {
    if(argc == 3) {
      strcpy(filein, argv[1]);
      strcpy(filein2, argv[2]);
    } else {
      sprintf(num, "%03d", snapshot);
      strcpy(filein, argv[1]);
      strcat(filein, "_");
      strcat(filein, num);
      strcpy(filein2, argv[2]);
      strcat(filein2, "_");
      strcat(filein2, num);
    }
      
    strcpy(fileout, filein);
    strcat(fileout, "_all.dat");
    strcat(filein, ".dat");
    strcat(filein2, ".dat");
    cout << "generating file " << fileout << " ......" <<endl;

    ifs.open(filein);
    if(!ifs)  { cout<<"stream error!"<<endl; exit(-1);}
    ifs2.open(filein2);
    if(!ifs2) { cout<<"stream error!"<<endl; exit(-1);}
    ofs.open(fileout);
    if(!ofs)  { cout<<"stream error!"<<endl; exit(-1);}
    ofs.setf(ios::scientific, ios::floatfield);

    ifs >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s 
        >> s >> s >> s >> s >> s >> s >> s >> s >> snum >> s >> s >> s;
    std::string strall(snum);
    strall.erase(strall.end()-1);
    strall.erase(0,2);
    //std::cout << strall << std::endl;
    int internalContact = atol(strall.c_str());
    //std::cout << internalContact << std::endl;

    ifs2 >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s >> snum >> s >> s >> s;
    std::string strall2(snum);
    strall2.erase(strall2.end()-1);
    strall2.erase(0,2);
    std::cout << strall2 << std::endl;
    int bdryContact = atol(strall2.c_str());
    std::cout << bdryContact << std::endl;

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
	<< endl;

    ofs << "ZONE I=" << internalContact + bdryContact <<", DATAPACKING=POINT, STRANDID=3, SOLUTIONTIME=" << snapshot << endl;

    double x, y, z, normal_x, normal_y, normal_z, tangt_x, tangt_y, tangt_z, total_x, total_y, total_z, normal, shear, total, penetration;
    for(int it = 0; it < internalContact; ++it) {
      ifs >> x >> y >> z >> normal_x >> normal_y >> normal_z >> tangt_x >> tangt_y >> tangt_z >> total_x >> total_y >> total_z >> normal >> shear >> total >> penetration;
      ofs << setw(OWID) << x
	  << setw(OWID) << y
	  << setw(OWID) << z
	  << setw(OWID) << normal_x
	  << setw(OWID) << normal_y
	  << setw(OWID) << normal_z
	  << setw(OWID) << tangt_x
	  << setw(OWID) << tangt_y
	  << setw(OWID) << tangt_z
	  << setw(OWID) << total_x
	  << setw(OWID) << total_y
	  << setw(OWID) << total_z
	  << setw(OWID) << normal
	  << setw(OWID) << shear
	  << setw(OWID) << total
	  << setw(OWID) << penetration
	  << endl;
    }

    for(int it = 0; it < bdryContact; ++it) {
      ifs2 >> x >> y >> z >> normal_x >> normal_y >> normal_z >> tangt_x >> tangt_y >> tangt_z >> penetration;
      normal= sqrt(normal_x * normal_x + normal_y * normal_y + normal_z * normal_z);
      shear = sqrt(tangt_x * tangt_x + tangt_y * tangt_y + tangt_z * tangt_z);
      ofs << setw(OWID) << x
	  << setw(OWID) << y
	  << setw(OWID) << z
	  << setw(OWID) << -normal_x
	  << setw(OWID) << -normal_y
	  << setw(OWID) << -normal_z
	  << setw(OWID) << -tangt_x
	  << setw(OWID) << -tangt_y
	  << setw(OWID) << -tangt_z
	  << setw(OWID) << -normal_x - tangt_x
	  << setw(OWID) << -normal_y - tangt_y
	  << setw(OWID) << -normal_z - tangt_z
	  << setw(OWID) << normal
	  << setw(OWID) << shear
	  << setw(OWID) << sqrt(normal*normal + shear*shear)
	  << setw(OWID) << penetration
	  << endl;
    }
	
    ifs.close();
    ifs2.close();
    ofs.close();

  }

  return 0;
}
