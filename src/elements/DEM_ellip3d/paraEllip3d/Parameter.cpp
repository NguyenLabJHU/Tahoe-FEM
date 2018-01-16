#include "Parameter.h"
#include "const.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstddef>
#include <cstdlib>
#include <iomanip>

namespace dem {

  void Parameter::readIn(const char *input) {
    std::ifstream ifs;
    ifs.open(input);
    if (!ifs) { debugInf << "stream error: Parameter.cpp" << std::endl; exit(-1); }
    std::string line;
    std::istringstream ssline;
    std::string str, str2;
    REAL val;

    // 25 generic parameters
    int geneNum = 25;
    std::size_t i;
    for (i = 0; i < geneNum; ++i) {
      while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
      ssline.clear(); ssline.str(line);
      ssline >> str >> val;
      parameter[str] = val;
    }
    if (i != geneNum) {std::cout << "generic parameters missing ..." << std::endl; exit(-1);}

    // for different types of simulation
    std::size_t simuType = static_cast<std::size_t> (parameter["simuType"]);
    switch (simuType) {
    case 001: { // proceed from preset state
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 4; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 002: { // tuneMassPercentage
      std::size_t i;
      for (i = 0; i < 12; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      for (i = 0; i < static_cast<std::size_t> (parameter["sieveNum"]); ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	REAL percent, size;
	ssline >> percent >> size;
	gradation.push_back(std::make_pair(percent, size));
      } 
    }
      break;

    case 003: { // trimOnly
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 1; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 004: { // removeBySphere
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 5; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 005: { // calculate mass percentage
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 1; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 101: { // depositIntoContainer  
      std::size_t i;
      int typeNum = 12;
      for (i = 0; i < typeNum; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	if (ifs.eof()) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      if (i != typeNum) {std::cout << "type parameters part (1) missing ..." << std::endl; exit(-1);}

      typeNum = static_cast<std::size_t> (parameter["sieveNum"]);
      for (i = 0; i < typeNum; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	if (ifs.eof()) break;
	ssline.clear(); ssline.str(line);
	REAL percent, size;
	ssline >> percent >> size;
	gradation.push_back(std::make_pair(percent, size));
      }
      if (i != typeNum) {std::cout << "type parameters part (2) missing ..." << std::endl; exit(-1);}

      typeNum = 7;
      for (i = 0; i < typeNum; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	if (ifs.eof()) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      } 
      if (i != typeNum) {std::cout << "type parameters part (3) missing ..." << std::endl; exit(-1);}

    }
      break;
  
    case 102: { // resumeDepositIntoContainer 
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 5; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 201: { // isotropic 1
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 7; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 202: { // isotropic 2
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 8; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 203: { // isotropic 3
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 4; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      for (i = 0; i < static_cast<std::size_t> (parameter["sigmaPoints"]); ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	REAL sigma;
	ssline >> sigma;
	sigmaPath.push_back(sigma);
      } 
      for (i = 0; i < 3; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 301: { // odometer 1
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 8; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 302: { // odometer 2
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 4; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      for (i = 0; i < static_cast<std::size_t> (parameter["sigmaPoints"]); ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	REAL sigma;
	ssline >> sigma;
	sigmaPath.push_back(sigma);
      } 
      for (i = 0; i < 3; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 401: { // triaxial 1
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 6; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 402: { // triaxial 2
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 7; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 411: { // plane strain 1
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 7; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 412: {// plain strain 2
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 8; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 501: {// true triaxial 1
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 10; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 502: { // true triaxial 2
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 11; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 601: { // expandCavityParticle
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 11; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 602: {// resumeExpandCavityParticle
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 4; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 603: { // oedometer impact from top
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 3; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
    }
      break;

    case 701: { // couple with gas flow, bottom "left" part, R-H conditions
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 28; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      std::size_t printPtclNum = static_cast<std::size_t> (parameter["printPtclNum"]);
      cfdPrintPtcls.resize(printPtclNum);;
      for (i = 0; i < printPtclNum; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> cfdPrintPtcls[i];
      } 
    }
      break;

    case 702: { // couple with gas flow, bottom "left" part
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 30; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      std::size_t printPtclNum = static_cast<std::size_t> (parameter["printPtclNum"]);
      cfdPrintPtcls.resize(printPtclNum);;
      for (i = 0; i < printPtclNum; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> cfdPrintPtcls[i];
      }
    }
      break;

    case 703: { // couple with gas flow, rectangular "left" part
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 35; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      std::size_t printPtclNum = static_cast<std::size_t> (parameter["printPtclNum"]);
      cfdPrintPtcls.resize(printPtclNum);;
      for (i = 0; i < printPtclNum; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> cfdPrintPtcls[i];
      }
    }
      break;

    case 704: { // couple with gas flow, spherical "left" part
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 33; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      std::size_t printPtclNum = static_cast<std::size_t> (parameter["printPtclNum"]);
      cfdPrintPtcls.resize(printPtclNum);;
      for (i = 0; i < printPtclNum; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> cfdPrintPtcls[i];
      }
    }
      break;

    case 705: { // couple with gas flow, rectangular "left" part with a zone below
      std::size_t i;
      for (i = 0; i < 2; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> str2;
	datafile[str] = str2;
      }
      for (i = 0; i < 38; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> str >> val;
	parameter[str] = val;
      }
      std::size_t printPtclNum = static_cast<std::size_t> (parameter["printPtclNum"]);
      cfdPrintPtcls.resize(printPtclNum);;
      for (i = 0; i < printPtclNum; ++i) {
	while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
	ssline.clear(); ssline.str(line);
	ssline >> cfdPrintPtcls[i];
      }
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
    std::vector<std::size_t> &ptcl = dem::Parameter::getSingleton().cfdPrintPtcls;
  
    for (std::map<std::string, REAL>::const_iterator it = param.begin(); it != param.end(); ++it)
      debugInf << std::setw(OWID) << it->first << std::setw(OWID) << it->second << std::endl;
  
    if (grada.size() != 0)
      debugInf << std::setw(OWID) << "gradation:" << std::endl;
    for (std::size_t i = 0; i < grada.size(); ++i) 
      debugInf << std::setw(OWID) << grada[i].first << std::setw(OWID) << grada[i].second << std::endl;

    for (std::map<std::string, std::string>::const_iterator it = file.begin(); it != file.end(); ++it)
      debugInf << std::setw(OWID) << it->first << std::setw(OWID) << it->second << std::endl;  

    if (sigma.size() != 0)
      debugInf << std::setw(OWID) << "sigma:" << std::endl;
    for (std::vector<REAL>::const_iterator it = sigma.begin(); it != sigma.end(); ++it)
      debugInf << std::setw(OWID) << (*it) << std::endl;

    if (ptcl.size() != 0)
      debugInf << std::setw(OWID) << "particles:" << std::endl;
    for (std::vector<std::size_t>::const_iterator it = ptcl.begin(); it != ptcl.end(); ++it)
      debugInf << std::setw(OWID) << (*it) << std::endl;

    debugInf << std::endl;
  }
  
}
