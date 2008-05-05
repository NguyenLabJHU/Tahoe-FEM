#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

int main(int argc, char *argv[])
{
    if(argc < 5)
    {
	cout << endl
	     << "-- Colletct a Particle's Information --" << endl
	     << "Usage: probe1p particle_number particle_file_prefix  first_suffix  last_suffix  suffix_increment" << endl
	     << "   --example: probe1p 25 dep_particle 0 100 5" << endl << endl;
	return -1;
    }	

    int first, last, incre, pnumber, n;
    first   = atoi(argv[3]);
    last    = atoi(argv[4]);
    incre   = atoi(argv[5]);
    pnumber = atoi(argv[1]);

    ifstream ifs1;
    ofstream ofs;
    char filein1[50];
    char fileout[50];
    char num[4], s[20];

    int ID, type, TotalNum, RORS;
    long double cx, cy, cz, rd, wd, lt, ht;
    long double a, b, c, x0, y0, z0, l1, l2, l3, m1, m2, m3, n1, n2, n3;
    long double vx, vy, vz, ox, oy, oz, fx, fy, fz, mx, my, mz;

    strcpy(fileout, "particle_");
    char pnum[6];
    sprintf(pnum, "%05d", pnumber);
    strcat(fileout, pnum);

    ofs.open(fileout);
    if(!ofs)   { cout<<"ofs stream error!"<<endl; exit(-1);}
    ofs.setf(ios::scientific, ios::floatfield);
    ofs<<"  snapshot      ID         radius_a        radius_b        radius_c   "
       <<"   position_x      position_y      position_z   "
       <<"    axle_a_x         axle_a_y        axle_a_z    "
       <<"    axle_b_x        axle_b_y        axle_b_z    "
       <<"    axle_c_x        axle_c_y        axle_c_z    "
       <<"   velocity_x      velocity_y      velocity_z   "
       <<"     omga_x          omga_y          omga_z     "
       <<"    force_x         force_y         force_z     "
       <<"    moment_x        moment_y        moment_z    "<<endl;
 
    for(n=first; n<=last; n+=incre)
    {
	sprintf(num, "%03d", n);
	strcpy(filein1, argv[2]);
	strcat(filein1, "_");
	strcat(filein1, num);

	ifs1.open(filein1);
	if(!ifs1)  { cout<<"ifs1 stream error!"<<endl; exit(-1);}
	ifs1 >> TotalNum >> RORS;
	if(RORS==0)
	    ifs1 >> cx >> cy >> cz >> rd >> ht;
	else
	    ifs1 >> cx >> cy >> cz >> wd >> lt >> ht;
	ifs1>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
	    >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s;
	
	do {
	    ifs1 >> ID >> type >> a >> b >> c >> x0 >> y0 >> z0 >> l1 >> m1 >> n1 >> l2 >> m2 >> n2 >> l3 >> m3 >> n3
		 >> vx >> vy >> vz >> ox >> oy >> oz >> fx >> fy >> fz >> mx >> my >> mz;
	} while(ID!=pnumber);


	ofs << setw(10) << n  << setw(10) << ID 
	    << setw(16) << a  << setw(16) << b  << setw(16) << c
	    << setw(16) << x0 << setw(16) << y0 << setw(16) << z0
	    << setw(16) << l1 << setw(16) << m1 << setw(16) << n1 
	    << setw(16) << l2 << setw(16) << m2 << setw(16) << n2 
	    << setw(16) << l3 << setw(16) << m3 << setw(16) << n3
	    << setw(16) << vx << setw(16) << vy << setw(16) << vz
	    << setw(16) << ox << setw(16) << oy << setw(16) << oz
	    << setw(16) << fx << setw(16) << fy << setw(16) << fz
	    << setw(16) << mx << setw(16) << my << setw(16) << mz << endl;

	ifs1.close();
    }
    ofs.close();

    return 0;
}
