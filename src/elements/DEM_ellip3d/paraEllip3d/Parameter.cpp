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

  // 26 generic parameters
  for (int i = 0; i < 26; ++i) {
    while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
    ssline.clear(); ssline.str(line);
    ssline >> str >> val;
    parameter[str] = val;
  }

  if ((int) parameter["simuType"] == 0) {
    for (int i = 0; i < 11; ++i) {
      while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
      ssline.clear(); ssline.str(line);
      ssline >> str >> val;
      parameter[str] = val;
    }
    
    for (int i = 0; i < (int) parameter["sieveNum"]; ++i) {
      while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
      ssline.clear(); ssline.str(line);
      REAL percent, size;
      ssline >> percent >> size;
      gradation.push_back(std::make_pair(percent, size));
    }
    
  }
  else if ((int) parameter["simuType"] == 1) {
    for (int i = 0; i < 2; ++i) {
      while (getline(ifs, line) ) if (line[0] != '#' && line.compare("") != 0 ) break;
      ssline.clear(); ssline.str(line);
      ssline >> str >> str2;
      datafile[str] = str2;
    }
  } 
  else if ((int) parameter["simuType"] == 2) {
  }
  
  ifs.close();
  
}
  
void Parameter::writeOut() {
  std::map<std::string, REAL> &param = dem::Parameter::getSingleton().parameter;
  std::vector<std::pair<REAL, REAL> > &grada = dem::Parameter::getSingleton().gradation;
  std::map<std::string, std::string> &file = dem::Parameter::getSingleton().datafile;
  
  for (std::map<std::string, REAL>::const_iterator it = param.begin(); it != param.end(); ++it)
    std::cout << it->first << "  " << it->second << std::endl;
  
  for (int i = 0; i < grada.size(); ++i) 
    std::cout << grada[i].first << "  " << grada[i].second << std::endl;

  for (std::map<std::string, std::string>::const_iterator it = file.begin(); it != file.end(); ++it)
    std::cout << it->first << "  " << it->second << std::endl;  
}
  
}
