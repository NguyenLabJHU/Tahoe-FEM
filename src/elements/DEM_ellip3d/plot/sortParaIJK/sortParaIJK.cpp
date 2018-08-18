#include "Vec.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cstdlib>
#include <map>

class IJK {
public:
  std::size_t i;
  std::size_t j;
  std::size_t k;

public:
  IJK() : i(0), j(0), k(0) {}
  IJK(int _i, int _j, int _k) : i(_i), j(_j), k(_k) {}
};

struct cmpIJK {
  bool operator()(const IJK &v1, const IJK &v2) const {

    /*
    // compare in i,j,k order
    if (v1.i < v2.i) {
      return true;
    } else if (v1.i == v2.i) {

      if (v1.j < v2.j) {
	return true;
      } else if (v1.j == v2.j) {

	if (v1.k < v2.k) {
	  return true;
	} else {
	  return false;
	}

      } else {
	return false;
      }

    } else {
      return false;
    }

    */
    // compare in k,j,i order, only in this order does Tecplot plot correctly (Tecplot requires I vary fastest for IJK-ordered data).
    if (v1.k < v2.k) {
      return true;
    } else if (v1.k == v2.k) {

      if (v1.j < v2.j) {
	return true;
      } else if (v1.j == v2.j) {

	if (v1.i < v2.i) {
	  return true;
	} else {
	  return false;
	}

      } else {
	return false;
      }

    } else {
      return false;
    }

  }
}; 

class GasVar {
public:
  Vec    coords;
  IJK    ijk;
  double Mach;
  double density;
  Vec    momentum;
  double energy;
  Vec    velocity;
  double pressure;
  double temperature;
  double mask;
  Vec    penalForce;
  Vec    pressureForce;

public:
  GasVar() : coords(0), ijk(0,0,0), Mach(0), density(0), momentum(0), energy(0), velocity(0), pressure(0), temperature(0), mask(0), penalForce(0), pressureForce(0) {}
  GasVar(Vec co, IJK i1, double m, double d, Vec mom, double en, Vec ve, double pr, double te, int ma, Vec pena, Vec pres) 
    : coords(co), ijk(i1), Mach(m), density(d), momentum(mom), energy(en), velocity(ve), pressure(pr), temperature(te), mask(ma), penalForce(pena), pressureForce(pres) {}

};

char *combineString(char *cstr, const char *str, int num, int width);

int main(int argc, char *argv[])
{
  if(argc < 5) {
    std::cout << std::endl 
	      << "-- sort unordered IJK data obtained from parallel IO --" << std::endl
	      << "-- note the original files are to be replaced --" << std::endl
	      << "  Usage:  sortParaIJK  file_prefix  start_snap  end_snap  increment" << std::endl
	      << "Example:  sortParaIJK  couple_fluidplot  0  100  1" << std::endl << std::endl;
    return -1;
  }

  const int OWID = 15; 
  const int OPREC = 6; 
  
  int startSnap = atoi(argv[2]);
  int endSnap   = atoi(argv[3]);
  int incSnap   = atoi(argv[4]);

  std::map<IJK, GasVar, cmpIJK> gasVarMap;

  ///////////////////////////////////////////////////////////
  for (int snapLoop = startSnap; snapLoop <= endSnap; snapLoop += incSnap) {
    gasVarMap.clear();

    ///////////////////////////////////////////////////////////
    char str[50], argsuf[50];
    strcpy(argsuf, argv[1]); strcat(argsuf, "_");
    combineString(str, argsuf, snapLoop, 3);
    strcat(str, ".dat");
    std::cout << "processing and replacing " << str << " ......" << std::endl;

    std::ifstream ifs(str);//, std::ifstream::binary);
    if(!ifs) { std::cout << "ifstream error." << std::endl; exit(-1); }
    std::string line, line1, line2;
    getline(ifs, line1);
    getline(ifs, line2);

    int lineNum = 0;
    ///////////////////////////////////////////////////////////
    while (getline(ifs, line)) {
      std::istringstream ss(line);
      double x, y, z;
      int i, j, k;
      double Mach;
      double density;
      double m1, m2, m3;
      double energy;
      double v1, v2, v3;
      double pressure;
      double temperature;
      double mask;
      double pen1, pen2, pen3;
      double pre1, pre2, pre3;
      ss >> x
	 >> y
	 >> z
	 >> i
	 >> j
	 >> k
	 >> Mach
	 >> density
	 >> m1
	 >> m2
	 >> m3
	 >> energy
	 >> v1
	 >> v2
	 >> v3
	 >> pressure
	 >> temperature
	 >> mask
	 >> pen1
	 >> pen2
	 >> pen3
	 >> pre1
	 >> pre2
	 >> pre3;
 
      gasVarMap[IJK(i,j,k)] = GasVar( Vec(x,y,z), IJK(i,j,k), Mach, density, Vec(m1,m2,m3), energy, Vec(v1,v2,v3), pressure, temperature, mask, Vec(pen1,pen2,pen3), Vec(pre1,pre2,pre3) );
      ++lineNum;
    }
    ifs.close();
    remove(str); // remove file

    // Tecplot example:
    //VARIABLES = x              y              z        iGlobal        jGlobal        kGlobal           Mach        density      momentumX      momentumY      momentumZ         energy      velocityX      velocityY      velocityZ       pressure    temperature           mask        penalFx        penalFy        penalFz     pressureFx     pressureFy     pressureFz
    //ZONE I=             17, J=             17, K=             17, DATAPACKING=POINT, STRANDID=1, SOLUTIONTIME=             10

    std::ofstream ofs(str);
    //std::string nstr = "new";
    //std::ofstream ofs((nstr+std::string(str)).c_str());
    if(!ofs) { std::cout << "ofstream error." << std::endl; exit(-1); }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);

    //ofs << line1 << std::endl;
    ofs << "  VARIABLES = x              y              z           Mach        density      momentumX      momentumY      momentumZ         energy      velocityX      velocityY      velocityZ       pressure    temperature           mask        penalFx        penalFy        penalFz     pressureFx     pressureFy     pressureFz" << std::endl;
    ofs << line2 << std::endl;
    
    //std::cout << "lineNum=" << lineNum << " mapSize=" << gasVarMap.size() << std::endl;
    std::map<IJK, GasVar, cmpIJK>::const_iterator it;    
    for (it = gasVarMap.begin(); it != gasVarMap.end(); ++it)
      ofs << std::setw(OWID) << it->second.coords.getX()       
	  << std::setw(OWID) << it->second.coords.getY()
	  << std::setw(OWID) << it->second.coords.getZ()
	//<< std::setw(OWID) << it->first.i     
	//<< std::setw(OWID) << it->first.j    
	//<< std::setw(OWID) << it->first.k     
	  << std::setw(OWID) << it->second.Mach     
	  << std::setw(OWID) << it->second.density     
	  << std::setw(OWID) << it->second.momentum.getX()
	  << std::setw(OWID) << it->second.momentum.getY()
	  << std::setw(OWID) << it->second.momentum.getZ()
	  << std::setw(OWID) << it->second.energy   
	  << std::setw(OWID) << it->second.velocity.getX()
	  << std::setw(OWID) << it->second.velocity.getY()
	  << std::setw(OWID) << it->second.velocity.getZ()
	  << std::setw(OWID) << it->second.pressure 
	  << std::setw(OWID) << it->second.temperature  
	  << std::setw(OWID) << it->second.mask
	  << std::setw(OWID) << it->second.penalForce.getX()
	  << std::setw(OWID) << it->second.penalForce.getY()
	  << std::setw(OWID) << it->second.penalForce.getZ()
	  << std::setw(OWID) << it->second.pressureForce.getX()
	  << std::setw(OWID) << it->second.pressureForce.getY()
	  << std::setw(OWID) << it->second.pressureForce.getZ()
	  << std::endl;

    ofs.close();
  } // snapLoop

  return 0;
}

char *combineString(char *cstr, const char *str, int num, int width) {
  std::string obj(str);
  std::stringstream ss;
  ss << std::setw(width) << std::setfill('0') << std::right << num;
  obj += ss.str();
  return strcpy( cstr, obj.c_str() );
}


