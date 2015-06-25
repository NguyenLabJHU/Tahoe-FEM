#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>

// this file is to generate boundary for deposition, 
// no inputs, specific dimensions will be defined in the code directly
// June 25, 2015, Boning Zhang

using namespace std;
int main(){

    long double center_x = 0;	// centroid on the bottom boundary
    long double center_y = 0;
    long double center_z = 0;

    long double radius = 6.35e-3;	// radius of the cylinder, unit m
    long double height = 10.3e-3;	// height of the cylinder, unit m
    long double PI = 3.141592653589;

    int isInside = -1;	// -1 means particles are inside the cylinder, 1 means particles are outside of cylinder

    int OWID = 16;
    std::ofstream ofs("ini_boundary_tri");
    if(!ofs) {
	std::cout << "stream error!" << std::endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);

    ofs << setw(OWID) << 0
        << setw(OWID) << 3 << std::endl;	// 3 boundaries since trixial
    
    // the boundary number 1, bottom boundary
    ofs << setw(OWID) << 1 << std::endl
	<< setw(OWID) << 1 << setw(OWID) << 2 << setw(OWID) << PI*radius*radius << std::endl
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) <<-1 << setw(OWID) << center_x << setw(OWID) << center_y << setw(OWID) << center_z << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative z boundary
	<< setw(OWID) << 2 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << center_x << setw(OWID) << center_y << setw(OWID) << center_z << setw(OWID) << radius << setw(OWID) << isInside << std::endl	// cylindrical side boundary
	<< std::endl

    // the boundary number 2, top boundary
        << setw(OWID) << 1 << std::endl
	<< setw(OWID) << 2 << setw(OWID) << 2 << setw(OWID) << PI*radius*radius << std::endl
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << center_x << setw(OWID) << center_y << setw(OWID) << center_z+height << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive z boundary
	<< setw(OWID) << 2 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << center_x << setw(OWID) << center_y << setw(OWID) << center_z << setw(OWID) << radius << setw(OWID) << isInside << std::endl	// cylindrical side boundary
	<< std::endl

    // the boundary nubmer 3, side boundary
        << setw(OWID) << 2 << std::endl
	<< setw(OWID) << 3 << setw(OWID) << 3 << setw(OWID) << 2.0*PI*radius*height << std::endl
	<< setw(OWID) << 2 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << center_x << setw(OWID) << center_y << setw(OWID) << center_z << setw(OWID) << radius << setw(OWID) << isInside << std::endl	// cylindrical side boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) <<-1 << setw(OWID) << center_x << setw(OWID) << center_y << setw(OWID) << center_z << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative z boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << center_x << setw(OWID) << center_y << setw(OWID) << center_z+height << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// positive z boundary
	<< std::endl;
    

    ofs.close();


    return 0;

}
