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
    int isInside = -1;	// -1 means particles are inside the cylinder, 1 means particles are outside of cylinder

    int OWID = 16;
    std::ofstream ofs("ini_boundary_dep");
    if(!ofs) {
	std::cout << "stream error!" << std::endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);

    ofs << setw(OWID) << 0
        << setw(OWID) << 2 << std::endl;	// 2 boundaries since deposition
    
    // the boundary number 1, bottom boundary
    ofs << setw(OWID) << 1 << std::endl
	<< setw(OWID) << 1 << setw(OWID) << 2 << setw(OWID) << 0 << std::endl
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) <<-1 << setw(OWID) << center_x << setw(OWID) << center_y << setw(OWID) << center_z << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative z boundary
	<< setw(OWID) << 2 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << center_x << setw(OWID) << center_y << setw(OWID) << center_z << setw(OWID) << radius << setw(OWID) << isInside << std::endl	// cylindrical side boundary
	<< std::endl

    // the boundary nubmer 2, side boundary
        << setw(OWID) << 2 << std::endl
	<< setw(OWID) << 2 << setw(OWID) << 2 << setw(OWID) << 0 << std::endl
	<< setw(OWID) << 2 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) << 1 << setw(OWID) << center_x << setw(OWID) << center_y << setw(OWID) << center_z << setw(OWID) << radius << setw(OWID) << isInside << std::endl	// cylindrical side boundary
	<< setw(OWID) << 1 << setw(OWID) << 0 << setw(OWID) << 0 << setw(OWID) <<-1 << setw(OWID) << center_x << setw(OWID) << center_y << setw(OWID) << center_z << setw(OWID) << 0 << setw(OWID) << 0 << std::endl	// negative z boundary
	<< std::endl;
    

    ofs.close();


    return 0;

}
