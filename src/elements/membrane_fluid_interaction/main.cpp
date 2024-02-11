#include "realtypes.h"
#include "surfacetension.h" 
#include <cstring>
#include <stdlib.h>
#include <ctime>

// main program
int main(int argc, char* argv[])
{
    // Part 0: command line arguments and timestamps
    time_t time1, time2;
    time(&time1);
    if (argc == 2) {
//	dem::NUM_THREADS = atoi(argv[1]);
	std::cout << "command line: " << argv[0] << " " << argv[1] << std::endl; 
    }


    // Part 1: start simulation
    memFluid::surfacetension axi1;
    axi1.solve("axi1", "Axi1", 100);

    //axi1.solve(const std::string& plain_axi, const char* filename, int interval);


    // Part 2: record run time
    time(&time2);
    std::cout << std::endl 
	      << "simulation start time: " << ctime(&time1);
    std::cout << "simulation  end  time: " << ctime(&time2);
    return 0;
}

