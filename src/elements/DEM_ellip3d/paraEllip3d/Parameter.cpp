#include "Parameter.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

namespace dem {

void Parameter::readIn(const char *input) {
  std::ifstream ifs;
  ifs.open(input);
  if (!ifs) { std::cout << "stream error!" << std::endl; exit(-1); }
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
  std::size_t simuType = (std::size_t) parameter["simuType"];
  switch (simuType) {
  case 1: // depositIntoContainer  
    for (std::size_t i = 0; i < 11; ++i) {
      while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
      ssline.clear(); ssline.str(line);
      ssline >> str >> val;
      parameter[str] = val;
    }
    for (std::size_t i = 0; i < (std::size_t) parameter["sieveNum"]; ++i) {
      while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
      ssline.clear(); ssline.str(line);
      REAL percent, size;
      ssline >> percent >> size;
      gradation.push_back(std::make_pair(percent, size));
    } 
    break;
  
  case 2: // resumeDepositIntoContainer 
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

  case 3: // expandCavityParticle
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

  case 4: // resumeExpandCavityParticle
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

  case 5: // trimOnly
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

  case 6: // isotropic 1
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

  case 7: // isotropic 2
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

  case 8: // isotropic 3
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
    for (std::size_t i = 0; i < (std::size_t) parameter["sigmaPoints"]; ++i) {
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

  case 100: // couple with sonic fluid flow
    for (std::size_t i = 0; i < 2; ++i) {
      while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
      ssline.clear(); ssline.str(line);
      ssline >> str >> str2;
      datafile[str] = str2;
    }
    for (std::size_t i = 0; i < 8; ++i) {
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
