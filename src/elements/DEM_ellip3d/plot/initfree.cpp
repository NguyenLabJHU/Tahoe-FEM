#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

int main(int argc, char *argv[])
{
    if(argc < 3)
    {
	cout << endl
	     << "-- Colletct Free Particles' Information --" << endl
	     << "Usage:" << endl
	     << "1) process a single file:  initfree fixed_particle_number particle_file" << endl
	     << "   --example: initfree 900 pile_particle_008" << endl
	     << "2) process multiple files: initfree fixed_particle_number particle_file_prefix  first_suffix  last_suffix  suffix_increment" << endl
	     << "   --example: initfree 900 pile_particle 0 100 5" << endl << endl;
	return -1;
    }	

    int first, last, incre, pnumber, n;
    if(argc == 3) {
	first = 0;
	last  = 1;
	incre = 2;
    }
    else
    {
	first   = atoi(argv[3]);
	last    = atoi(argv[4]);
	incre   = atoi(argv[5]);
    }

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

    for(n=first; n<=last; n+=incre)
    {
	if(argc == 3){
	    strcpy(filein1, argv[2]);
	    strcpy(fileout, filein1);
	    strcat(fileout, "_free");
	}
	else
	{
	    sprintf(num, "%03d", n);
	    strcpy(filein1, argv[2]);
	    strcat(filein1, "_");
	    strcpy(fileout, filein1);
	    strcat(fileout,"free_");
	    strcat(filein1, num);
	    strcat(fileout, num);
	}

	ifs1.open(filein1);
	if(!ifs1)  { cout<<"ifs1 stream error!"<<endl; exit(-1);}

	ofs.open(fileout);
	if(!ofs)   { cout<<"ofs stream error!"<<endl; exit(-1);}
	ofs.setf(ios::scientific, ios::floatfield);
	cout<<fileout<<"..."<<endl;

	ifs1 >> TotalNum >> RORS;
	ofs<<setw(10)<<TotalNum-pnumber<<setw(10)<<RORS<<endl;

	if(RORS==0) {
	    ifs1 >> cx >> cy >> cz >> rd >> ht;
	    ofs<<setw(16)<<cx
	       <<setw(16)<<cy
	       <<setw(16)<<cz
	       <<setw(16)<<rd
	       <<setw(16)<<ht<<endl;
	}
	else {
	    ifs1 >> cx >> cy >> cz >> wd >> lt >> ht;
	    ofs<<setw(16)<<cx
	       <<setw(16)<<cy
	       <<setw(16)<<cz
	       <<setw(16)<<wd
	       <<setw(16)<<lt
	       <<setw(16)<<ht<<endl;
	}
	ifs1>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
	    >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s;

	ofs<<"     ID         type      radius_a        radius_b        radius_c   "
	   <<"    position_x            position_y            position_z        "
	   <<"    axle_a_x              axle_a_y              axle_a_z          "
	   <<"    axle_b_x              axle_b_y              axle_b_z          "
	   <<"    axle_c_x              axle_c_y              axle_c_z          "
	   <<"   velocity_x      velocity_y      velocity_z   "
	   <<"     omga_x          omga_y          omga_z     "
	   <<"    force_x         force_y         force_z     "
	   <<"    moment_x        moment_y        moment_z    "
	   <<endl;
	
	for (int k = 0; k < TotalNum; ++k) {
	    ifs1 >> ID >> type >> a >> b >> c >> x0 >> y0 >> z0 >> l1 >> m1 >> n1 >> l2 >> m2 >> n2 >> l3 >> m3 >> n3
		 >> vx >> vy >> vz >> ox >> oy >> oz >> fx >> fy >> fz >> mx >> my >> mz;
	    if (ID > pnumber) {
		ofs  << setw(10) << ID << setw(10) << type << setw(16) << a << setw(16) << b << setw(16) << c;
		ofs.precision(12);
		ofs  << setw(22) << x0 << setw(22) << y0 << setw(22) << z0
		     << setw(22) << l1 << setw(22) << m1 << setw(22) << n1
		     << setw(22) << l2 << setw(22) << m2 << setw(22) << n2
		     << setw(22) << l3 << setw(22) << m3 << setw(22) << n3;
		ofs.precision(6);
		ofs  << setw(16) << vx << setw(16) << vy << setw(16) << vz
		     << setw(16) << ox << setw(16) << oy << setw(16) << oz
		     << setw(16) << fx << setw(16) << fy << setw(16) << fz
		     << setw(16) << mx << setw(16) << my << setw(16) << mz << endl;
	    }

	} 

	char add;
	while (ifs1) {
	    ifs1.get(add);
	    ofs.put(add);
	}

	ifs1.close();
	ofs.close();
    }


    return 0;
}
