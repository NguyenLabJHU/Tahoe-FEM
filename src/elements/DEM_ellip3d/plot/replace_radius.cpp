#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>

const std::size_t OWID  = 15;   // output width
const std::size_t OPREC = 6;    // output precision, number of digits after decimal dot

int main(int argc, char *argv[])
{
  if(argc < 4) {
    std::cout << std::endl 
	      << "-- Replace particle radius --" << std::endl
	      << "Usage:" << std::endl
	      << "replace_radius  particle_file  radius_b  radius_c" << std::endl
	      << "Example:" << std::endl
	      << "replace_radius  deposit_particle_end  2.1E-3  1.6E-3" << std::endl << std::endl;
    return -1;
  }	

  double new_b = atof(argv[2]);
  double new_c = atof(argv[3]);

  std::ifstream ifs;
  std::ofstream ofs;
  char filein[50];
  char fileout[50];

  int totalNum, id, type;
  double a, b, c, x, y, z, l1, l2, l3, m1, m2, m3, n1, n2, n3;
  double vx, vy, vz, omx, omy, omz, fx, fy, fz, mx, my, mz;

  strcpy(filein, argv[1]);
  strcpy(fileout, filein);
  strcat(fileout, ".new");

  ifs.open(filein);
  if(!ifs) { std::cout << "stream error!" << std::endl; exit(-1); }
  ofs.open(fileout);
  if(!ofs) { std::cout << "stream error!" << std::endl; exit(-1); }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  ifs >> totalNum;
  std::string s[29];
  for (int i = 0; i < 29; ++i)
    ifs >> s[i];

  ofs << std::setw(OWID) << totalNum << std::endl;
  for (int i = 0; i < 29; ++i)
    ofs << std::setw(OWID) << s[i];
  ofs << std::endl;

  int nid = 0;
  for(int ptcl = 0; ptcl < totalNum; ++ptcl) {
    ifs >> id >> type >> a >> b >> c >> x >> y >> z >> l1 >> m1 >> n1 >> l2 >> m2 >> n2 >> l3 >> m3 >> n3
	>> vx >> vy >> vz >> omx >> omy >> omz >> fx >> fy >> fz >> mx >> my >> mz;
    ofs << std::setw(OWID) << id << std::setw(OWID) << type << std::setw(OWID) 
	<< std::setw(OWID) << a  << std::setw(OWID) << new_b << std::setw(OWID) << new_c 
	<< std::setw(OWID) << x  << std::setw(OWID) << y  << std::setw(OWID) << z 
	<< std::setw(OWID) << l1 << std::setw(OWID) << m1 << std::setw(OWID) << n1 
	<< std::setw(OWID) << l2 << std::setw(OWID) << m2 << std::setw(OWID) << n2 
	<< std::setw(OWID) << l3 << std::setw(OWID) << m3 << std::setw(OWID) << n3
	<< std::setw(OWID) << vx << std::setw(OWID) << vy << std::setw(OWID) << vz 
	<< std::setw(OWID) << omx<< std::setw(OWID) << omy<< std::setw(OWID) << omz 
	<< std::setw(OWID) << fx << std::setw(OWID) << fy << std::setw(OWID) << fz 
	<< std::setw(OWID) << mx << std::setw(OWID) << my << std::setw(OWID) << mz << std::endl;
  }

  std::size_t sieveNum;
  ifs >> sieveNum;
  std::vector<double> percent(sieveNum), size(sieveNum);
  for (std::size_t i = 0; i < sieveNum; ++i)
    ifs >> percent[i] >> size[i];
  double ratio_ba, ratio_ca;
  ifs >> ratio_ba >> ratio_ca;

  ofs << std::endl << std::setw(OWID) << sieveNum << std::endl;
  for (std::size_t i = 0; i < sieveNum; ++i)
    ofs << std::setw(OWID) << percent[i] << std::setw(OWID) << size[i] << std::endl;
  ofs << std::endl << std::setw(OWID) << ratio_ba << std::setw(OWID) << ratio_ca << std::endl;


  ifs.close();
  ofs.close();

  return 0;
}
