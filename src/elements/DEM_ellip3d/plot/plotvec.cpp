#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

const long double PI = 3.1415927;

int main(int argc, char *argv[])
{
    if(argc < 2)
    {
	cout << endl 
	     << "-- Plot Vectors such as Force or Velocity --"<<endl
	     << "Usage:" << endl
	     << "1) process a single file:  plotvec particle_file" << endl
	     << "   --example: plotvec iso_particle_008" << endl
	     << "2) process multiple files: plotvec particle_file_prefix  first_suffix  last_suffix  suffix_increment" << endl
	     << "   --example: plotvec iso_particle  0  100  5" << endl << endl;
	return -1;
    }	

    int first, last, incre;
    if(argc == 2) {
	first = 0;
	last  = 1;
	incre = 2;
    }
    else
    {
	first = atoi(argv[2]);
	last  = atoi(argv[3]);
	incre = atoi(argv[4]);
    }

    ifstream ifs;
    ofstream ofs;
    char filein[50];
    char fileout[50];
    char num[4], s[20];

    int ID, type, TotalNum, RORS;
    long double cx, cy, cz, rd, wd, lt, ht;
    long double a, b, c, x0, y0, z0, l1, l2, l3, m1, m2, m3, n1, n2, n3, v1, v2, v3, w1, w2, w3, f1, f2, f3, mt1, mt2, mt3;
    int n, k;

    for(n=first; n<=last; n+=incre)
    {
	if(argc == 2)
	    strcpy(filein, argv[1]);
	else
	{
	    sprintf(num, "%03d", n);
	    strcpy(filein, argv[1]);
	    strcat(filein, "_");
	    strcat(filein, num);
	}

	strcpy(fileout, filein);
	strcat(fileout, "_vec.dat");
        cout << "generating file " << fileout << " ......" <<endl;

	ifs.open(filein);
	if(!ifs)  { cout<<"stream error!"<<endl; exit(-1);}
	ofs.open(fileout);
	if(!ofs)  { cout<<"stream error!"<<endl; exit(-1);}
	ofs.setf(ios::scientific, ios::floatfield);

	ifs >> TotalNum >> RORS;
	if(RORS==0)
	    ifs >> cx >> cy >> cz >> rd >> ht;
	else
	    ifs >> cx >> cy >> cz >> wd >> lt >> ht;

	ifs>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
	   >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s;

	ofs << "VARIABLES=" << endl
	    << setfill(' ') << setw(16) << "x"
	    << setfill(' ') << setw(16) << "y"
	    << setfill(' ') << setw(16) << "z"
	    << setfill(' ') << setw(16) << "velocity_x"
	    << setfill(' ') << setw(16) << "velocity_y"
	    << setfill(' ') << setw(16) << "velocity_z"
	    << setfill(' ') << setw(16) << "force_x"
	    << setfill(' ') << setw(16) << "force_y"
	    << setfill(' ') << setw(16) << "force_z"
	    << setfill(' ') << setw(16) << "omga_x"
	    << setfill(' ') << setw(16) << "omga_y"
	    << setfill(' ') << setw(16) << "omga_z"
	    << setfill(' ') << setw(16) << "moment_x"
	    << setfill(' ') << setw(16) << "moment_y"
	    << setfill(' ') << setw(16) << "moment_z"
	    << endl;

	ofs << "ZONE I=" << TotalNum <<", DATAPACKING=POINT" << endl;
	for(k = 0; k < TotalNum; ++k)
	{
	    ifs >> ID >> type >> a >> b >> c >> x0 >> y0 >> z0 >> l1 >> m1 >> n1 >> l2 >> m2 >> n2 >> l3 >> m3 >> n3
		>>v1>>v2>>v3>>w1>>w2>>w3>>f1>>f2>>f3>>mt1>>mt2>>mt3;
	    
	    ofs << setfill(' ') << setw(16) << x0
		<< setfill(' ') << setw(16) << y0
		<< setfill(' ') << setw(16) << z0
		<< setfill(' ') << setw(16) << v1
		<< setfill(' ') << setw(16) << v2
		<< setfill(' ') << setw(16) << v3
		<< setfill(' ') << setw(16) << f1
		<< setfill(' ') << setw(16) << f2
		<< setfill(' ') << setw(16) << f3
		<< setfill(' ') << setw(16) << w1
		<< setfill(' ') << setw(16) << w2
		<< setfill(' ') << setw(16) << w3
		<< setfill(' ') << setw(16) << mt1
		<< setfill(' ') << setw(16) << mt2
		<< setfill(' ') << setw(16) << mt3
		<< endl;   
	}
	
	ifs.close();
	ofs.close();
    }

    return 0;
}
