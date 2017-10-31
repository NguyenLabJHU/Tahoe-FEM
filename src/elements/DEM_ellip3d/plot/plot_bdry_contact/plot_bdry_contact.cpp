#include "Vec.h"
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
static const int OPREC= 6;

class BdryContact {
public:
  Vec point;
  Vec normal;
  Vec tangt;
  REAL penetr;

public:
  BdryContact()
    :point(0), normal(0), tangt(0), penetr(0) 
  {}
    
  BdryContact(Vec pt, Vec nm, Vec tg, REAL pntr)
    :point(pt), normal(nm), tangt(tg), penetr(pntr) 
  {}
    
  void print(std::ostream &os) {
    os << std::setw(OWID) << point.getX()
       << std::setw(OWID) << point.getY()
       << std::setw(OWID) << point.getZ()
       << std::setw(OWID) << normal.getX()
       << std::setw(OWID) << normal.getY()
       << std::setw(OWID) << normal.getZ() 
       << std::setw(OWID) << tangt.getX()
       << std::setw(OWID) << tangt.getY()
       << std::setw(OWID) << tangt.getZ() 
       << std::setw(OWID) << penetr << std::endl;
  }
};

int main(int argc, char *argv[])
{
  if(argc < 2) {
    cout << endl 
	 << "-- Plot contact force on boundaries --"<<endl
	 << "Usage:" << endl
	 << "1) process a single file:  plot_bdry_contact boundary_contact_file" << endl
	 << "   --example: plot_bdry_contact triaxial_bdrycntc_008" << endl
	 << "2) process multiple files: plot_bdry_contact boudnary_contact_file_prefix  first_suffix  last_suffix  suffix_increment" << endl
	 << "   --example: plot_bdry_contact triaxial_bdrycntc  1  100  5" << endl << endl;
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

    std::vector<BdryContact> contactInfo;
  for(int n=first; n<=last; n+=incre) {
    if(argc == 2)
      strcpy(filein, argv[1]);
    else {
      sprintf(num, "%03d", n);
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

    contactInfo.clear();
    int bdryIndex, bdryPointNum;
    double pos_x, pos_y, pos_z, normal_x, normal_y, normal_z, tangt_x, tangt_y, tangt_z, pentr;

    int totalBdryNum; // = 5;
    ifs >> totalBdryNum;
    for (int bdryNum = 1; bdryNum <= totalBdryNum; ++bdryNum) {

      ifs >> bdryIndex >> bdryPointNum >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s;
      for (int pointNum = 0; pointNum < bdryPointNum; ++pointNum) {
	ifs >> pos_x >> pos_y >> pos_z >> normal_x >> normal_y >> normal_z >> tangt_x >> tangt_y >> tangt_z >> pentr;
	contactInfo.push_back(BdryContact(Vec(pos_x, pos_y, pos_z), Vec(normal_x, normal_y, normal_z), Vec(tangt_x, tangt_y, tangt_z), pentr ));
      }
      ifs >> normal_x >> normal_y >> normal_z >> tangt_x >> tangt_y >> tangt_z >> pentr;

    }

    ofs	<< std::setw(OWID) << "VARIABLES =" << std::endl
	<< std::setw(OWID) << "pos_x"
	<< std::setw(OWID) << "pos_y"
	<< std::setw(OWID) << "pos_z"
	<< std::setw(OWID) << "normal_x"
	<< std::setw(OWID) << "normal_y"
	<< std::setw(OWID) << "normal_z"
	<< std::setw(OWID) << "tangt_x"
	<< std::setw(OWID) << "tangt_y"
	<< std::setw(OWID) << "tangt_z"
	<< std::setw(OWID) << "pentr" << std::endl;

    ofs << "ZONE I=" << contactInfo.size() <<", DATAPACKING=POINT" << std::endl;
    for (std::vector<BdryContact>::iterator it = contactInfo.begin(); it != contactInfo.end(); ++it)
      it->print(ofs);

    ifs.close();
    ofs.close();

  }

  return 0;
}
