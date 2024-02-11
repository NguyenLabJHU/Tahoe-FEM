#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>

// this file is to generate boundary for deposition, 
// no inputs, specific dimensions will be defined in the code directly
// July 25, 2014, Boning Zhang

using namespace std;
int main(){

    // dimensions of boundary, need to be changed for specific boundary
    long double x_max = 0.0015;
    long double x_min =-0.0015;
    long double y_max = 0.0015;
    long double y_min =-0.0015;
    long double z_min = 0;
    long double z_max = 0.003;


    long double x_mid = (x_max+x_min)*0.5;
    long double y_mid = (y_max+y_min)*0.5;
    long double z_mid = (z_max+z_min)*0.5;

    long double area_x = (z_max-z_min)*(y_max-y_min);
    long double area_y = (z_max-z_min)*(x_max-x_min);
    long double area_z = (x_max-x_min)*(y_max-y_min);

    int OWID = 16;

    std::ofstream ofs("flo_boundary_end");
    if(!ofs) {
	std::cout << "stream error!" << std::endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);

    ofs << setw(OWID) << 0
        << setw(OWID) << 6 << std::endl;	// 5 boundaries since deposition
    
    // the boundary number 1
    ofs << setw(OWID) << 1 << std::endl
	<< setw(OWID) << 1 << setw(OWID) << 5 << setw(OWID) << area_x << std::endl
	<< setw(OWID) << 1 << setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << x_max << setw(OWID) << y_mid << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive x boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) <<-1 << setw(OWID) << 0 << setw(OWID) << x_mid << setw(OWID) << y_min << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative y boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << x_mid << setw(OWID) << y_max << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive y boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << x_mid << setw(OWID) << y_mid << setw(OWID) << z_max << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive z boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) <<-1 << setw(OWID) << x_mid << setw(OWID) << y_mid << setw(OWID) << z_min << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative z boundary
	<< std::endl

    // the boundary nubmer 2
        << setw(OWID) << 1 << std::endl
	<< setw(OWID) << 2 << setw(OWID) << 5 << setw(OWID) << area_y << std::endl
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << x_mid << setw(OWID) << y_max << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive y boundary
	<< setw(OWID) << 1 << setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << x_max << setw(OWID) << y_mid << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive x boundary
	<< setw(OWID) << 1 << setw(OWID) <<-1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << x_min << setw(OWID) << y_mid << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative x boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << x_mid << setw(OWID) << y_mid << setw(OWID) << z_max << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive z boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) <<-1 << setw(OWID) << x_mid << setw(OWID) << y_mid << setw(OWID) << z_min << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative z boundary
	<< std::endl

    // the boundary number 3
        << setw(OWID) << 1 << std::endl
	<< setw(OWID) << 3 << setw(OWID) << 5 << setw(OWID) << area_x << std::endl
	<< setw(OWID) << 1 << setw(OWID) <<-1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << x_min << setw(OWID) << y_mid << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative x boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) <<-1 << setw(OWID) << 0 << setw(OWID) << x_mid << setw(OWID) << y_min << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative y boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << x_mid << setw(OWID) << y_max << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive y boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << x_mid << setw(OWID) << y_mid << setw(OWID) << z_max << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive z boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) <<-1 << setw(OWID) << x_mid << setw(OWID) << y_mid << setw(OWID) << z_min << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative z boundary
	<< std::endl

    // the boundary number 4
        << setw(OWID) << 1 << std::endl
	<< setw(OWID) << 4 << setw(OWID) << 5 << setw(OWID) << area_y << std::endl
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) <<-1 << setw(OWID) << 0 << setw(OWID) << x_mid << setw(OWID) << y_min << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative y boundary
	<< setw(OWID) << 1 << setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << x_max << setw(OWID) << y_mid << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive x boundary
	<< setw(OWID) << 1 << setw(OWID) <<-1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << x_min << setw(OWID) << y_mid << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative x boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << x_mid << setw(OWID) << y_mid << setw(OWID) << z_max << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive z boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) <<-1 << setw(OWID) << x_mid << setw(OWID) << y_mid << setw(OWID) << z_min << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative z boundary
	<< std::endl

    // the boundary number 5
        << setw(OWID) << 1 << std::endl
	<< setw(OWID) << 5 << setw(OWID) << 5 << setw(OWID) << area_z << std::endl
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << x_mid << setw(OWID) << y_mid << setw(OWID) << z_max << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive z boundary
	<< setw(OWID) << 1 << setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << x_max << setw(OWID) << y_mid << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive x boundary
	<< setw(OWID) << 1 << setw(OWID) <<-1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << x_min << setw(OWID) << y_mid << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative x boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << x_mid << setw(OWID) << y_max << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive y boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) <<-1 << setw(OWID) << 0 << setw(OWID) << x_mid << setw(OWID) << y_min << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative y boundary
	<< std::endl

    // the boundary number 6
        << setw(OWID) << 1 << std::endl
	<< setw(OWID) << 6 << setw(OWID) << 5 << setw(OWID) << area_z << std::endl
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << -1 << setw(OWID) << x_mid << setw(OWID) << y_mid << setw(OWID) << z_min << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative z boundary
	<< setw(OWID) << 1 << setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << x_max << setw(OWID) << y_mid << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive x boundary
	<< setw(OWID) << 1 << setw(OWID) <<-1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << x_min << setw(OWID) << y_mid << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative x boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << x_mid << setw(OWID) << y_max << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive y boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) <<-1 << setw(OWID) << 0 << setw(OWID) << x_mid << setw(OWID) << y_min << setw(OWID) << z_mid << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative y boundary
	<< std::endl;
    

    ofs.close();


    return 0;

}
