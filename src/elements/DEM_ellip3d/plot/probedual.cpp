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
	     << "-- Colletct Dual Particles' Contact Information --" << endl
	     << "Usage: probedual contact_file_prefix  first_suffix  last_suffix  suffix_increment" << endl
	     << "   --example: probedual dep_contact 0 100 5" << endl << endl;
	return -1;
    }	

    int first, last, incre, n, ic;
    first   = atoi(argv[2]);
    last    = atoi(argv[3]);
    incre   = atoi(argv[4]);

    ifstream ifs1;
    ofstream ofs;
    char filein1[50];
    char fileout[50];
    char num[4], s[20];

    int TotalCntct, ptcl_1, ptcl_2;
    long double point1_x, point1_y, point1_z, point2_x, point2_y, point2_z;
    long double radius_1, radius_2, penetration, tang_dispmt, contact_radius, R0, E0, normal_force, shear_force;

    strcpy(fileout, argv[1]);
    strcat(fileout, "_dual");
    ofs.open(fileout);
    if(!ofs)   { cout<<"ofs stream error!"<<endl; exit(-1);}
    ofs.setf(ios::scientific, ios::floatfield);
    ofs<<"  snapshot     ptcl_1     ptcl_2    point1_x        point1_y        point1_z        "
       <<"point2_x       point2_y        point2_z         radius_1        radius_2      "
       <<"penetration     tang_dispmt     contact_radius       R0              E0         normal_force     tang_force"<<endl;
 
    for(n=first; n<=last; n+=incre)
    {
	sprintf(num, "%03d", n);
	strcpy(filein1, argv[1]);
	strcat(filein1, "_");
	strcat(filein1, num);

	ifs1.open(filein1);
	if(!ifs1)  { cout<<"ifs1 stream error!"<<endl; exit(-1);}
	ifs1 >> TotalCntct;
	ifs1 >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s;
	if (TotalCntct==0)
	    ofs <<setw(10)<<n<<setw(10)<<0<<setw(10)<<0<<setw(16)<<0<<setw(16)<<0<<setw(16)<<0
		<<setw(16)<<0<<setw(16)<<0<<setw(16)<<0<<setw(16)<<0<<setw(16)<<0<<setw(16)<<0
		<<setw(16)<<0<<setw(16)<<0<<setw(16)<<0<<setw(16)<<0<<setw(16)<<0<<setw(16)<<0<<endl;
	else{
	    for (ic=0;ic<TotalCntct;++ic) {
		ifs1>>ptcl_1>>ptcl_2>>point1_x>>point1_y>>point1_z>>point2_x>>point2_y>>point2_z
		    >>radius_1>>radius_2>>penetration>>tang_dispmt>>contact_radius>>R0>>E0>>normal_force>>shear_force;
		ofs <<setw(10)<<n<<setw(10)<<ptcl_1<<setw(10)<<ptcl_2<<setw(16)<<point1_x<<setw(16)<<point1_y<<setw(16)<<point1_z
		    <<setw(16)<<point2_x<<setw(16)<<point2_y<<setw(16)<<point2_z<<setw(16)<<radius_1<<setw(16)<<radius_2
		    <<setw(16)<<penetration<<setw(16)<<tang_dispmt<<setw(16)<<contact_radius<<setw(16)<<R0<<setw(16)<<E0
		    <<setw(16)<<normal_force<<setw(16)<<shear_force<<endl;
	    }
	    
	}
	ifs1.close();
    }
    ofs.close();

    return 0;
}
