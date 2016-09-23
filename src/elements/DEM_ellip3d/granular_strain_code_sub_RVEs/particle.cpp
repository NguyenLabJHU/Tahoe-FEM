#include "particle.h"
#include <iostream>



particle::particle(int n, REAL* posi_x, REAL* posi_y, REAL* posi_z, REAL tmp_v){
    ID = n;
    position_x = posi_x;
    position_y = posi_y;
    position_z = posi_z;

    curr_position.setx(position_x[0]);
    curr_position.sety(position_y[0]);
    curr_position.setz(position_z[0]);
    prev_position = curr_position;

    volume = tmp_v;


} // end particle

void particle::update(int step, REAL time_interval){

    prev_position = curr_position;
    curr_position.setx(position_x[step]);
    curr_position.sety(position_y[step]);
    curr_position.setz(position_z[step]);

    curr_velocity = (curr_position-prev_position)/time_interval;

} // end update
