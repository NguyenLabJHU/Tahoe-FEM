#include "Parameter.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstddef>
#include <cstdlib>

namespace dem {

  void Parameter::readIn(const char *input) {
    std::ifstream ifs;
    ifs.open(input);
    if (!ifs) { std::cout << "stream error: Parameter.cpp" << std::endl; exit(-1); }
    std::string line;
    std::istringstream ssline;
    std::string str, str2;
    REAL val;

    // 28 generic parameters
    for (std::size_t i = 0; i < 28; ++i) {
      while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
      ssline.clear(); ssline.str(line);
      ssline >> str >> val;
      parameter[str] = val;
    }

    // for different types of simulation
    std::size_t simuType = static_cast<std::size_t> (parameter["simuType"]);
    switch (simuType) {
    case 001: // proceed from preset state
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 1; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 002: // tuneMassPercentage
      for (std::size_t i = 0; i < 12; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      for (std::size_t i = 0; i < static_cast<std::size_t> (parameter["sieveNum"]); ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	REAL percent, size;
	ssline >> percent >> size;
	gradation.push_back(std::make_pair(percent, size));
      } 
      break;

    case 003: // trimOnly
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 1; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 004: // removeBySphere
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 5; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 005: // calculate mass percentage
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 1; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 101: // depositIntoContainer  
      for (std::size_t i = 0; i < 12; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      for (std::size_t i = 0; i < static_cast<std::size_t> (parameter["sieveNum"]); ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	REAL percent, size;
	ssline >> percent >> size;
	gradation.push_back(std::make_pair(percent, size));
      } 
      break;
  
    case 102: // resumeDepositIntoContainer 
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 201: // isotropic 1
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 6; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 202: // isotropic 2
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 7; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 203: // isotropic 3
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 3; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      for (std::size_t i = 0; i < static_cast<std::size_t> (parameter["sigmaPoints"]); ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	REAL sigma;
	ssline >> sigma;
	sigmaPath.push_back(sigma);
      } 
      for (std::size_t i = 0; i < 3; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 301: // odometer 1
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 7; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 302: // odometer 2
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 3; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      for (std::size_t i = 0; i < static_cast<std::size_t> (parameter["sigmaPoints"]); ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	REAL sigma;
	ssline >> sigma;
	sigmaPath.push_back(sigma);
      } 
      for (std::size_t i = 0; i < 3; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 401: // triaxial 1
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 5; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 402: // triaxial 2
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 6; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 411: // plane strain 1
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 6; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 412: // plain strain 2
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 7; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 501: // true triaxial 1
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 9; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 502: // true triaxial 2
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 10; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 601: // expandCavityParticle
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 7; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 602: // resumeExpandCavityParticle
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 1; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 701: // couple with sonic fluid flow, bottom "left" part, R-H conditions
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 25; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 702: // couple with sonic fluid flow, bottom "left" part
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 27; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 703: // couple with sonic fluid flow, rectangular "left" part
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 32; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 704: // couple with sonic fluid flow, spherical "left" part
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 30; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    case 705: // couple with sonic fluid flow, rectangular "left" part with a zone below
      for (std::size_t i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (std::size_t i = 0; i < 35; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      break;

    }
    ifs.close();
  
  }
  
  void Parameter::writeOut() {
    std::map<std::string, REAL> &param = dem::Parameter::getSingleton().parameter;
    std::vector<std::pair<REAL, REAL> > &grada = dem::Parameter::getSingleton().gradation;
    std::map<std::string, std::string> &file = dem::Parameter::getSingleton().datafile;
    std::vector<REAL> &sigma = dem::Parameter::getSingleton().sigmaPath;
  
    for (std::map<std::string, REAL>::const_iterator it = param.begin(); it != param.end(); ++it)
      std::cout << it->first << "  " << it->second << std::endl;
  
    for (std::size_t i = 0; i < grada.size(); ++i) 
      std::cout << grada[i].first << "  " << grada[i].second << std::endl;

    for (std::map<std::string, std::string>::const_iterator it = file.begin(); it != file.end(); ++it)
      std::cout << it->first << "  " << it->second << std::endl;  

    for (std::vector<REAL>::const_iterator it = sigma.begin(); it != sigma.end(); ++it)
      std::cout << (*it) << std::endl;
  }
  
}
