#include <iostream>
#include <string>
#include <vector>

// user defined head files
#include "SPHDomain.h"

using namespace std;
using namespace sph;

// Main Function
int main(int argc, char * argv[]) {

    SPHDomain sphCase;
    sphCase.CouetteFlow();
    return 0;
}
