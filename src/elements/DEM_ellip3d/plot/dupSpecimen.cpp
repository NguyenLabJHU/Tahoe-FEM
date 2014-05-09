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
  if(argc < 6) {
    std::cout << std::endl 
	      << "-- duplicate specimen in x, y, z directions --" << std::endl
	      << "Usage:" << std::endl
	      << "dupSpecimen  boundary_file  particle_file  2  2  2" << std::endl
	      << "where boundary_file is a closed container"<< std::endl << std::endl;
    return -1;
  }	

  int dupX = atoi(argv[3]);
  int dupY = atoi(argv[4]);
  int dupZ = atoi(argv[5]);

  std::ifstream ifs;
  std::ofstream ofs;
  char filein[50];
  char fileout[50];

  ////////////////////////////////// boundary file
  strcpy(filein, argv[1]);
  strcpy(fileout, filein);
  strcat(fileout, ".dup");
  ifs.open(filein);
  if(!ifs) { std::cout << "stream error!" << std::endl; exit(-1); }
  ofs.open(fileout);
  if(!ofs) { std::cout << "stream error!" << std::endl; exit(-1); }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  double x1, y1, z1, x2, y2, z2;
  ifs >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z1
      << std::setw(OWID) << x1 + (x2-x1)*dupX 
      << std::setw(OWID) << y1 + (y2-y1)*dupY 
      << std::setw(OWID) << z1 + (z2-z1)*dupZ << std::endl << std::endl
      << std::setw(OWID) << 6 << std::endl << std::endl;

  double x0 = x1 + (x2-x1)*dupX/2;
  double y0 = y1 + (y2-y1)*dupY/2;
  double z0 = z1 + (z2-z1)*dupZ/2;

  // boundary 1
  ofs << std::setw(OWID) << 1
      << std::setw(OWID) << 0 << std::endl

      << std::setw(OWID) << 1      
      << std::setw(OWID) << -1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 
      << std::setw(OWID) << x1
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0 << std::endl << std::endl

    // boundary 2
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 << std::endl

      << std::setw(OWID) << 2
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0 
      << std::setw(OWID) << x2     
      << std::setw(OWID) << y0
      << std::setw(OWID) << z0 << std::endl << std::endl
      
    // boundary 3
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 << std::endl

      << std::setw(OWID) << 3
      << std::setw(OWID) << 0
      << std::setw(OWID) << -1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << x0
      << std::setw(OWID) << y1
      << std::setw(OWID) << z0 << std::endl << std::endl
      
    // boundary 4
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 << std::endl

      << std::setw(OWID) << 4
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y2
      << std::setw(OWID) << z0 << std::endl << std::endl
      
    // boundary 5
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 << std::endl

      << std::setw(OWID) << 5
      << std::setw(OWID) << 0
      << std::setw(OWID) << 0
      << std::setw(OWID) << -1
      << std::setw(OWID) << x0
      << std::setw(OWID) << y0
      << std::setw(OWID) << z1 << std::endl << std::endl
      
    // boundary 6
      << std::setw(OWID) << 1
      << std::setw(OWID) << 0 << std::endl

      << std::setw(OWID) << 6
      << std::setw(OWID) << 0 
      << std::setw(OWID) << 0
      << std::setw(OWID) << 1    
      << std::setw(OWID) << x0      
      << std::setw(OWID) << y0
      << std::setw(OWID) << z2 << std::endl << std::endl;
    
  ifs.close();
  ofs.close();
  
  ////////////////////////////////// boundary file
  int totalNum, id, type;
  double a, b, c, x, y, z, l1, l2, l3, m1, m2, m3, n1, n2, n3;
  double vx, vy, vz, omx, omy, omz, fx, fy, fz, mx, my, mz;

  strcpy(filein, argv[2]);
  strcpy(fileout, filein);
  strcat(fileout, ".dup");

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
  ofs << std::setw(OWID) << totalNum *dupX *dupY *dupZ << std::endl;
  for (int i = 0; i < 29; ++i)
    ofs << std::setw(OWID) << s[i];
  ofs << std::endl;

  int nid = 0;
  for(int ptcl = 0; ptcl < totalNum; ++ptcl) {
    ifs >> id >> type >> a >> b >> c >> x >> y >> z >> l1 >> m1 >> n1 >> l2 >> m2 >> n2 >> l3 >> m3 >> n3
	>> vx >> vy >> vz >> omx >> omy >> omz >> fx >> fy >> fz >> mx >> my >> mz;

    double xnew, ynew, znew;
    for (int i = 0; i < dupX; ++i) {
      xnew = x + i*(x2-x1);
      for (int j= 0; j < dupY; ++j) {
	ynew = y + j*(y2-y1);
	for (int k = 0; k < dupZ; ++k) {
	  znew = z + k*(z2-z1);
	  ofs << std::setw(OWID) << ++nid << std::setw(OWID) << type 
	      << std::setw(OWID) << a << std::setw(OWID) << b << std::setw(OWID) << c 
	      << std::setw(OWID) << xnew << std::setw(OWID) << ynew << std::setw(OWID) << znew 
	      << std::setw(OWID) << l1 << std::setw(OWID) << m1 << std::setw(OWID) << n1 
	      << std::setw(OWID) << l2 << std::setw(OWID) << m2 << std::setw(OWID) << n2 
	      << std::setw(OWID) << l3 << std::setw(OWID) << m3 << std::setw(OWID) << n3
	      << std::setw(OWID) << vx << std::setw(OWID) << vy << std::setw(OWID) << vz 
	      << std::setw(OWID) << omx << std::setw(OWID) << omy << std::setw(OWID) << omz 
	      << std::setw(OWID) << fx << std::setw(OWID) << fy << std::setw(OWID) << fz 
	      << std::setw(OWID) << mx << std::setw(OWID) << my << std::setw(OWID) << mz
	      << std::endl;
	}
      }
    }
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
