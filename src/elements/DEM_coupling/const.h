// Typically, the paramenters here do not need to change

#ifndef CONST_H
#define CONST_H

namespace dem { 

// PI value
const long double PI      = 3.141592653589793;

// absolute numerical precision
const long double PREC    = 1.0e-12;

// random shape for each particle
//#define RANDOM_SHAPE

// particle material property (quartz sand E=29GPa, v= 0.25)
const long double YOUNG   = 2.90e+10;  
const long double POISSON = 0.25;      
const long double Gs      = 2.65;     


// step interval to update contacts between particles
const int UPDATE_CNT      = 1; //50; 

}

#endif
