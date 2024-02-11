#include "particle.h"
#include "parameter.h"
#include "ran.h"
#include "root6.h"
#include "eig3.h"
#include <iostream>
#include <math.h>
#include <algorithm>
#include <iomanip>

//#define MOMENT
#ifdef MOMENT
const int START = 10000;  // at which time step to apply moment? for moment rotation test only.
#define SLIP  // if defined, stick and slip; otherwise slide.
#endif
//main.cpp: dem::TIMESTEP = 5.0e-07; A.deposit(12000,... 

//#define MINDLIN_ASSUMED

namespace dem {

void particle::init() {	// August 19, 2013
    // generate orientation of axle a/b/c
    REAL l1, m1, n1, l2, m2, n2, mod, A, B, C;
    
    while (1) {
      l1 = ran(&idum)*2 - 1;
      m1 = ran(&idum)*2 - 1;
      n1 = ran(&idum)*2 - 1;
      mod = sqrt(l1*l1 + m1*m1 + n1*n1);
      l1 /= mod;
      m1 /= mod;
      n1 /= mod;
      
      l2 = ran(&idum)*2 - 1;
      // solve n2
      A = m1*m1 + n1*n1;
      B = 2 * l1 * l2 * n1;
      C = l1*l1*l2*l2 + m1*m1*l2*l2 - m1*m1;
      if (B*B - 4*A*C > EPS)
	break;
    }

    int sign = 2*ran(&idum)-1 > 0 ? 1:-1;
    n2 = (-B + sign*sqrt(B*B - 4*A*C) ) / (2*A);
    m2 = - (l1*l2 + n1*n2) / m1;
    curr_direction_a=vec(acos(l1), acos(m1), acos(n1));
    curr_direction_b=vec(acos(l2), acos(m2), acos(n2));
    curr_direction_c=vacos(normalize(vcos(curr_direction_a)*vcos(curr_direction_b)));

    // local coordinate of center_geo of every octant
    REAL x1, x2, x3, x4, x5, x6, x7, x8;
    REAL y1, y2, y3, y4, y5, y6, y7, y8;
    REAL z1, z2, z3, z4, z5, z6, z7, z8;
    x1 = 3.0/8.0*aplus; y1 = 3.0/8.0*bplus; z1 = 3.0/8.0*cplus;
    x2 = -3.0/8.0*aminus; y2 = 3.0/8.0*bplus; z2 = 3.0/8.0*cplus;
    x3 = -3.0/8.0*aminus; y3 = -3.0/8.0*bminus; z3 = 3.0/8.0*cplus;
    x4 = 3.0/8.0*aplus; y4 = -3.0/8.0*bminus; z4 = 3.0/8.0*cplus;
    x5 = 3.0/8.0*aplus; y5 = 3.0/8.0*bplus; z5 = -3.0/8.0*cminus;
    x6 = -3.0/8.0*aminus; y6 = 3.0/8.0*bplus; z6 = -3.0/8.0*cminus;
    x7 = -3.0/8.0*aminus; y7 = -3.0/8.0*bminus; z7 = -3.0/8.0*cminus;
    x8 = 3.0/8.0*aplus; y8 = -3.0/8.0*bminus; z8 = -3.0/8.0*cminus;

    prev_position=curr_position;
    prev_direction_a=curr_direction_a;
    prev_direction_b=curr_direction_b;
    prev_direction_c=curr_direction_c;
    prev_velocity=curr_velocity=0;
    prev_omga=curr_omga=0;
    curr_acceleration=prev_acceleration=0;
    curr_acce_rotate = prev_acce_rotate = 0;
    force=pre_force=0;moment=pre_moment=0;mres=0;
    const_force=const_moment=0;
    density=Gs*1.0e3;
    REAL v1, v2, v3, v4, v5, v6, v7, v8;	// volumes of eight octants
    v1 = 1/6.0*PI*aplus*bplus*cplus;
    v2 = 1/6.0*PI*aminus*bplus*cplus;
    v3 = 1/6.0*PI*aminus*bminus*cplus;
    v4 = 1/6.0*PI*aplus*bminus*cplus;
    v5 = 1/6.0*PI*aplus*bplus*cminus;
    v6 = 1/6.0*PI*aminus*bplus*cminus;
    v7 = 1/6.0*PI*aminus*bminus*cminus;
    v8 = 1/6.0*PI*aplus*bminus*cminus;

    volume= v1+v2+v3+v4+v5+v6+v7+v8;
    // local coordinate of center of mass
    REAL xl_center_mass, yl_center_mass, zl_center_mass;
    xl_center_mass = (v1*x1+v2*x2+v3*x3+v4*x4+v5*x5+v6*x6+v7*x7+v8*x8)/volume;
    yl_center_mass = (v1*y1+v2*y2+v3*y3+v4*y4+v5*y5+v6*y6+v7*y7+v8*y8)/volume;
    zl_center_mass = (v1*z1+v2*z2+v3*z3+v4*z4+v5*z5+v6*z6+v7*z7+v8*z8)/volume;
    local_center_mass = vec(xl_center_mass, yl_center_mass, zl_center_mass);

    curr_center_mass = curr_position + globalVec(local_center_mass);
    prev_center_mass = curr_center_mass;	// August 19, 2013
    init_center_mass = curr_center_mass;	// initial center for granular strain
    start_center_mass = curr_center_mass;

    mass=density*volume;
    REAL Ixx, Iyy, Izz;
    // moment of inertia with respect to center of geometry
    Ixx = v1*density/5.0*(bplus*bplus+cplus*cplus)+v2*density/5.0*(bplus*bplus+cplus*cplus)
	 +v3*density/5.0*(bminus*bminus+cplus*cplus)+v4*density/5.0*(bminus*bminus+cplus*cplus)
         +v5*density/5.0*(bplus*bplus+cminus*cminus)+v6*density/5.0*(bplus*bplus+cminus*cminus)
	 +v7*density/5.0*(bminus*bminus+cminus*cminus)+v8*density/5.0*(bminus*bminus+cminus*cminus);

    Iyy = v1*density/5.0*(aplus*aplus+cplus*cplus)+v2*density/5.0*(aminus*aminus+cplus*cplus)
	 +v3*density/5.0*(aminus*aminus+cplus*cplus)+v4*density/5.0*(aplus*aplus+cplus*cplus)
         +v5*density/5.0*(aplus*aplus+cminus*cminus)+v6*density/5.0*(aminus*aminus+cminus*cminus)
	 +v7*density/5.0*(aminus*aminus+cminus*cminus)+v8*density/5.0*(aplus*aplus+cminus*cminus);

    Izz = v1*density/5.0*(aplus*aplus+bplus*bplus)+v2*density/5.0*(aminus*aminus+bplus*bplus)
	 +v3*density/5.0*(aminus*aminus+bminus*bminus)+v4*density/5.0*(aplus*aplus+bminus*bminus)
         +v5*density/5.0*(aplus*aplus+bplus*bplus)+v6*density/5.0*(aminus*aminus+bplus*bplus)
	 +v7*density/5.0*(aminus*aminus+bminus*bminus)+v8*density/5.0*(aplus*aplus+bminus*bminus);
    // moment of inertial with respect to center of mass
    Ixx = Ixx+xl_center_mass*xl_center_mass*mass;
    Iyy = Iyy+yl_center_mass*yl_center_mass*mass;
    Izz = Izz+zl_center_mass*zl_center_mass*mass;

    J=vec(Ixx,Iyy,Izz);
    cntnum=0;
    inContact=false;
    isAbleDivide = true;
    candidacy = 0;	// for particle fracture, October 17, 2013
    numBroken = 0;
    typeBroken= 0;
    delta_a = 0;
    delta_b = 0;	// the increment in b axle
    delta_c = 0;	// -1 means not transit

    weibullPhi = ran(&idum);
    REAL tmp_parameter = log(pow(1.0/weibullPhi, 1.0/dem::weibullModulus))
			*pow(getAverageRadius()/dem::basicRadius, -2*dem::weibullModulus);
    strengthHoek = dem::sigmaCompress*tmp_parameter;
    strengthContact = dem::ContactTensileCritical*tmp_parameter;

    matrix zero3x3(3,3);
    average_stress = zero3x3;
//    GlobCoef(1, 1, 1);	// calculate coefficients when it is needed, August 19, 2013
    GlobCoef();	// September 6, 2013

    // the integrals over poly-ellipsoids for average stress calculation
    int_x  = PI*0.0625*( aplus*aplus*bplus*cplus     - aminus*aminus*bplus*cplus 
			    - aminus*aminus*bminus*cplus  + aplus*aplus*bminus*cplus
			    + aplus*aplus*bplus*cminus    - aminus*aminus*bplus*cminus
			    - aminus*aminus*bminus*cminus + aplus*aplus*bminus*cminus );

    int_y  = PI*0.0625*( aplus*bplus*bplus*cplus     + aminus*bplus*bplus*cplus
			    - aminus*bminus*bminus*cplus  - aplus*bminus*bminus*cplus
			    + aplus*bplus*bplus*cminus    + aminus*bplus*bplus*cminus
			    - aminus*bminus*bminus*cminus - aplus*bminus*bminus*cminus );

    int_z  = PI*0.0625*( aplus*bplus*cplus*cplus     + aminus*bplus*cplus*cplus
			    + aminus*bminus*cplus*cplus   + aplus*bminus*cplus*cplus
			    - aplus*bplus*cminus*cminus   - aminus*bplus*cminus*cminus
			    - aminus*bminus*cminus*cminus - aplus*bminus*cminus*cminus ); 

    int_xy = 1.0/15.0*(  aplus*aplus*bplus*bplus*cplus      - aminus*aminus*bplus*bplus*cplus
			    + aminus*aminus*bminus*bminus*cplus  - aplus*aplus*bminus*bminus*cplus
			    + aplus*aplus*bplus*bplus*cminus     - aminus*aminus*bplus*bplus*cminus
			    + aminus*aminus*bminus*bminus*cminus - aplus*aplus*bminus*bminus*cminus );

    int_xz = 1.0/15.0*(  aplus*aplus*bplus*cplus*cplus      - aminus*aminus*bplus*cplus*cplus
			    - aminus*aminus*bminus*cplus*cplus   + aplus*aplus*bminus*cplus*cplus
			    - aplus*aplus*bplus*cminus*cminus    + aminus*aminus*bplus*cminus*cminus 
			    + aminus*aminus*bminus*cminus*cminus - aplus*aplus*bminus*cminus*cminus );

    int_yz = 1.0/15.0*(  aplus*bplus*bplus*cplus*cplus      + aminus*bplus*bplus*cplus*cplus 
			    - aminus*bminus*bminus*cplus*cplus   - aplus*bminus*bminus*cplus*cplus
			    - aplus*bplus*bplus*cminus*cminus    - aminus*bplus*bplus*cminus*cminus 
			    + aminus*bminus*bminus*cminus*cminus + aplus*bminus*bminus*cminus*cminus );

    int_x2 = PI/30.0*(   aplus*aplus*aplus*bplus*cplus      + aminus*aminus*aminus*bplus*cplus
			    + aminus*aminus*aminus*bminus*cplus  + aplus*aplus*aplus*bminus*cplus
			    + aplus*aplus*aplus*bplus*cminus     + aminus*aminus*aminus*bplus*cminus
			    + aminus*aminus*aminus*bminus*cminus + aplus*aplus*aplus*bminus*cminus );

    int_y2 = PI/30.0*(   aplus*bplus*bplus*bplus*cplus      + aminus*bplus*bplus*bplus*cplus
			    + aminus*bminus*bminus*bminus*cplus  + aplus*bminus*bminus*bminus*cplus
			    + aplus*bplus*bplus*bplus*cminus     + aminus*bplus*bplus*bplus*cminus
			    + aminus*bminus*bminus*bminus*cminus + aplus*bminus*bminus*bminus*cminus );

    int_z2 = PI/30.0*(   aplus*bplus*cplus*cplus*cplus      + aminus*bplus*cplus*cplus*cplus
			    + aminus*bminus*cplus*cplus*cplus    + aplus*bminus*cplus*cplus*cplus
			    + aplus*bplus*cminus*cminus*cminus   + aminus*bplus*cminus*cminus*cminus
			    + aminus*bminus*cminus*cminus*cminus + aplus*bminus*cminus*cminus*cminus );
}

particle::particle(int n, int tp, vec center_geo, REAL r, REAL yng, REAL poi)	// August 16, 2013
  :ID(n), type(tp), curr_position(center_geo), aplus(r), aminus(r), bplus(r), bminus(r), cplus(r), cminus(r), young(yng), poisson(poi) {
  init ();
}


particle::particle(int n, int tp, vec center_geo, REAL raplus, REAL raminus, REAL rbplus, REAL rbminus, REAL rcplus, REAL rcminus, REAL yng, REAL poi)	// August 16, 2013
  :ID(n), type(tp), curr_position(center_geo), aplus(raplus), aminus(raminus), bplus(rbplus), bminus(rbminus), cplus(rcplus), cminus(rcminus), young(yng), poisson(poi) {
  init ();
}


particle::particle(int n, int tp, vec center, gradation& grad, REAL yng, REAL poi)
  :ID(n), type(tp), curr_position(center), young(yng), poisson(poi) {	// August 19, 2013

  REAL atotal, btotal, ctotal;	// aplus+aminus, bplus+bminus, cplus+cminus
    // generate a particle based on gradation distribution
  REAL sievenum = grad.getSieveNum();
  for (int k=0;k<sievenum;k++){
    if (ran(&idum) <= grad.getPercent()[sievenum-1-k]){
      atotal=grad.getSize()[sievenum-1-k];
      break;
    }
  }
  
#ifdef RANDOM_SHAPE
  grad.setPtclRatioBA(ran(&idum));
  grad.setPtclRatioCA(ran(&idum));
#endif
  
  btotal=atotal*grad.getPtclRatioBA();
  ctotal=atotal*grad.getPtclRatioCA();
  
  REAL ratio_plus = ran(&idum);	// ratio of plus to total
  if (ratio_plus > 0.8)
	ratio_plus = 0.8;	// to avoid the poly-ellipsoid be too sharp
  else if (ratio_plus < 0.2)
	ratio_plus = 0.2;	// to avoid it be too sharp
  aplus = atotal*ratio_plus;
  aminus = atotal-aplus;

  ratio_plus = ran(&idum);	// ratio of plus to total
  if (ratio_plus > 0.8)
	ratio_plus = 0.8;	// to avoid the poly-ellipsoid be too sharp
  else if (ratio_plus < 0.2)
	ratio_plus = 0.2;	// to avoid it be too sharp
  bplus = btotal*ratio_plus;
  bminus = btotal-bplus;

  ratio_plus = ran(&idum);	// ratio of plus to total
  if (ratio_plus > 0.8)
	ratio_plus = 0.8;	// to avoid the poly-ellipsoid be too sharp
  else if (ratio_plus < 0.2)
	ratio_plus = 0.2;	// to avoid it be too sharp
  cplus = ctotal*ratio_plus;
  cminus = ctotal-cplus;  

  init();
}


particle::particle(int id, int tp, REAL raplus, REAL raminus, REAL rbplus, REAL rbminus, REAL rcplus, REAL rcminus, vec position_geo, vec dirca, vec dircb, vec dircc, REAL yng, REAL poi)	// August 16, 2013
  :ID(id), type(tp), young(yng), poisson(poi) {
    aplus=raplus;
    aminus=raminus;
    bplus=rbplus;
    bminus=rbminus;
    cplus=rcplus;
    cminus=rcminus;
    // local coordinate of center_geo of every octant
    REAL x1, x2, x3, x4, x5, x6, x7, x8;
    REAL y1, y2, y3, y4, y5, y6, y7, y8;
    REAL z1, z2, z3, z4, z5, z6, z7, z8;
    x1 = 3.0/8.0*aplus; y1 = 3.0/8.0*bplus; z1 = 3.0/8.0*cplus;
    x2 = -3.0/8.0*aminus; y2 = 3.0/8.0*bplus; z2 = 3.0/8.0*cplus;
    x3 = -3.0/8.0*aminus; y3 = -3.0/8.0*bminus; z3 = 3.0/8.0*cplus;
    x4 = 3.0/8.0*aplus; y4 = -3.0/8.0*bminus; z4 = 3.0/8.0*cplus;
    x5 = 3.0/8.0*aplus; y5 = 3.0/8.0*bplus; z5 = -3.0/8.0*cminus;
    x6 = -3.0/8.0*aminus; y6 = 3.0/8.0*bplus; z6 = -3.0/8.0*cminus;
    x7 = -3.0/8.0*aminus; y7 = -3.0/8.0*bminus; z7 = -3.0/8.0*cminus;
    x8 = 3.0/8.0*aplus; y8 = -3.0/8.0*bminus; z8 = -3.0/8.0*cminus;
    
    curr_position=prev_position=position_geo;	// initial position for granular strain
    curr_direction_a=prev_direction_a=dirca;
    curr_direction_b=prev_direction_b=dircb;
    curr_direction_c=prev_direction_c=dircc;
    curr_velocity=prev_velocity=0;
    curr_omga=prev_omga=0;
    curr_acceleration=prev_acceleration=0;
    curr_acce_rotate = prev_acce_rotate = 0;
    force=pre_force=0;
    moment=pre_moment=0;mres=0;
    const_force=const_moment=0;
    cntnum=0;
    density=Gs*1.0e3;
    REAL v1, v2, v3, v4, v5, v6, v7, v8;	// volumes of eight octants
    v1 = 1/6.0*PI*aplus*bplus*cplus;
    v2 = 1/6.0*PI*aminus*bplus*cplus;
    v3 = 1/6.0*PI*aminus*bminus*cplus;
    v4 = 1/6.0*PI*aplus*bminus*cplus;
    v5 = 1/6.0*PI*aplus*bplus*cminus;
    v6 = 1/6.0*PI*aminus*bplus*cminus;
    v7 = 1/6.0*PI*aminus*bminus*cminus;
    v8 = 1/6.0*PI*aplus*bminus*cminus;

    volume= v1+v2+v3+v4+v5+v6+v7+v8;
    // local coordinate of center of mass
    REAL xl_center_mass, yl_center_mass, zl_center_mass;
    xl_center_mass = (v1*x1+v2*x2+v3*x3+v4*x4+v5*x5+v6*x6+v7*x7+v8*x8)/volume;
    yl_center_mass = (v1*y1+v2*y2+v3*y3+v4*y4+v5*y5+v6*y6+v7*y7+v8*y8)/volume;
    zl_center_mass = (v1*z1+v2*z2+v3*z3+v4*z4+v5*z5+v6*z6+v7*z7+v8*z8)/volume;
    local_center_mass = vec(xl_center_mass, yl_center_mass, zl_center_mass);

    curr_center_mass = curr_position + globalVec(local_center_mass);
    prev_center_mass = curr_center_mass;	// August 19, 2013
    init_center_mass = curr_center_mass;
    start_center_mass = curr_center_mass;   

    mass=density*volume;
    REAL Ixx, Iyy, Izz;
    // moment of inertia with respect to center of geometry
    Ixx = v1*density/5.0*(bplus*bplus+cplus*cplus)+v2*density/5.0*(bplus*bplus+cplus*cplus)
	 +v3*density/5.0*(bminus*bminus+cplus*cplus)+v4*density/5.0*(bminus*bminus+cplus*cplus)
         +v5*density/5.0*(bplus*bplus+cminus*cminus)+v6*density/5.0*(bplus*bplus+cminus*cminus)
	 +v7*density/5.0*(bminus*bminus+cminus*cminus)+v8*density/5.0*(bminus*bminus+cminus*cminus);

    Iyy = v1*density/5.0*(aplus*aplus+cplus*cplus)+v2*density/5.0*(aminus*aminus+cplus*cplus)
	 +v3*density/5.0*(aminus*aminus+cplus*cplus)+v4*density/5.0*(aplus*aplus+cplus*cplus)
         +v5*density/5.0*(aplus*aplus+cminus*cminus)+v6*density/5.0*(aminus*aminus+cminus*cminus)
	 +v7*density/5.0*(aminus*aminus+cminus*cminus)+v8*density/5.0*(aplus*aplus+cminus*cminus);

    Izz = v1*density/5.0*(aplus*aplus+bplus*bplus)+v2*density/5.0*(aminus*aminus+bplus*bplus)
	 +v3*density/5.0*(aminus*aminus+bminus*bminus)+v4*density/5.0*(aplus*aplus+bminus*bminus)
         +v5*density/5.0*(aplus*aplus+bplus*bplus)+v6*density/5.0*(aminus*aminus+bplus*bplus)
	 +v7*density/5.0*(aminus*aminus+bminus*bminus)+v8*density/5.0*(aplus*aplus+bminus*bminus);
    // moment of inertial with respect to center of mass
    Ixx = Ixx+xl_center_mass*xl_center_mass*mass;
    Iyy = Iyy+yl_center_mass*yl_center_mass*mass;
    Izz = Izz+zl_center_mass*zl_center_mass*mass;

    J=vec(Ixx,Iyy,Izz);
    inContact=false;
    isAbleDivide = true;
    candidacy =0;	// for particle fracture. Initial candidacy=0; October 17, 2013
    numBroken =0;
    typeBroken=0;

    delta_a = 0;
    delta_b = 0;	// the increment in b axle
    delta_c = 0;	// -1 means not transit

    weibullPhi = ran(&idum);
    REAL tmp_parameter = log(pow(1.0/weibullPhi, 1.0/dem::weibullModulus))
			*pow(getAverageRadius()/dem::basicRadius, -2*dem::weibullModulus);
    strengthHoek = dem::sigmaCompress*tmp_parameter;
    strengthContact = dem::ContactTensileCritical*tmp_parameter;

    matrix zero3x3(3,3);
    average_stress = zero3x3;

//    GlobCoef(1, 1, 1);	// calculate coefficients when it is needed, August 16, 2013
    GlobCoef();	// September 6, 2013

    // the integrals over poly-ellipsoids for average stress calculation
    int_x  = PI*0.0625*( aplus*aplus*bplus*cplus     - aminus*aminus*bplus*cplus 
			    - aminus*aminus*bminus*cplus  + aplus*aplus*bminus*cplus
			    + aplus*aplus*bplus*cminus    - aminus*aminus*bplus*cminus
			    - aminus*aminus*bminus*cminus + aplus*aplus*bminus*cminus );

    int_y  = PI*0.0625*( aplus*bplus*bplus*cplus     + aminus*bplus*bplus*cplus
			    - aminus*bminus*bminus*cplus  - aplus*bminus*bminus*cplus
			    + aplus*bplus*bplus*cminus    + aminus*bplus*bplus*cminus
			    - aminus*bminus*bminus*cminus - aplus*bminus*bminus*cminus );

    int_z  = PI*0.0625*( aplus*bplus*cplus*cplus     + aminus*bplus*cplus*cplus
			    + aminus*bminus*cplus*cplus   + aplus*bminus*cplus*cplus
			    - aplus*bplus*cminus*cminus   - aminus*bplus*cminus*cminus
			    - aminus*bminus*cminus*cminus - aplus*bminus*cminus*cminus ); 

    int_xy = 1.0/15.0*(  aplus*aplus*bplus*bplus*cplus      - aminus*aminus*bplus*bplus*cplus
			    + aminus*aminus*bminus*bminus*cplus  - aplus*aplus*bminus*bminus*cplus
			    + aplus*aplus*bplus*bplus*cminus     - aminus*aminus*bplus*bplus*cminus
			    + aminus*aminus*bminus*bminus*cminus - aplus*aplus*bminus*bminus*cminus );

    int_xz = 1.0/15.0*(  aplus*aplus*bplus*cplus*cplus      - aminus*aminus*bplus*cplus*cplus
			    - aminus*aminus*bminus*cplus*cplus   + aplus*aplus*bminus*cplus*cplus
			    - aplus*aplus*bplus*cminus*cminus    + aminus*aminus*bplus*cminus*cminus 
			    + aminus*aminus*bminus*cminus*cminus - aplus*aplus*bminus*cminus*cminus );

    int_yz = 1.0/15.0*(  aplus*bplus*bplus*cplus*cplus      + aminus*bplus*bplus*cplus*cplus 
			    - aminus*bminus*bminus*cplus*cplus   - aplus*bminus*bminus*cplus*cplus
			    - aplus*bplus*bplus*cminus*cminus    - aminus*bplus*bplus*cminus*cminus 
			    + aminus*bminus*bminus*cminus*cminus + aplus*bminus*bminus*cminus*cminus );

    int_x2 = PI/30.0*(   aplus*aplus*aplus*bplus*cplus      + aminus*aminus*aminus*bplus*cplus
			    + aminus*aminus*aminus*bminus*cplus  + aplus*aplus*aplus*bminus*cplus
			    + aplus*aplus*aplus*bplus*cminus     + aminus*aminus*aminus*bplus*cminus
			    + aminus*aminus*aminus*bminus*cminus + aplus*aplus*aplus*bminus*cminus );

    int_y2 = PI/30.0*(   aplus*bplus*bplus*bplus*cplus      + aminus*bplus*bplus*bplus*cplus
			    + aminus*bminus*bminus*bminus*cplus  + aplus*bminus*bminus*bminus*cplus
			    + aplus*bplus*bplus*bplus*cminus     + aminus*bplus*bplus*bplus*cminus
			    + aminus*bminus*bminus*bminus*cminus + aplus*bminus*bminus*bminus*cminus );

    int_z2 = PI/30.0*(   aplus*bplus*cplus*cplus*cplus      + aminus*bplus*cplus*cplus*cplus
			    + aminus*bminus*cplus*cplus*cplus    + aplus*bminus*cplus*cplus*cplus
			    + aplus*bplus*cminus*cminus*cminus   + aminus*bplus*cminus*cminus*cminus
			    + aminus*bminus*cminus*cminus*cminus + aplus*bminus*cminus*cminus*cminus );

}

// this will create the sub-particle which is the sub-poly-ellipsoid in the negative axle of the two sub-poly-ellipsoids that are 
// broken from the original particle along the break plane
// the sub-poly-ellipsoid which is in the positive axle is used to change the original particle, this process is defined
// in the function particle::breakItSelf(int break_plane)
particle::particle(const particle & pt, int break_plane){ // break_plane 1--ab plane, 2--ac plane, 3--bc plane
    ID = pt.ID;
    type = pt.type;

    young = pt.young; poisson = pt.poisson;
    prev_position = pt.prev_position;
    prev_center_mass = pt.prev_center_mass;
    init_center_mass = pt.init_center_mass;	
    start_center_mass = pt.start_center_mass;
    curr_direction_a = pt.curr_direction_a; curr_direction_b = pt.curr_direction_b; curr_direction_c = pt.curr_direction_c;
    prev_direction_a = pt.prev_direction_a; prev_direction_b = pt.prev_direction_b; prev_direction_c = pt.prev_direction_c;

    curr_velocity = pt.curr_velocity; prev_velocity = pt.prev_velocity;
    curr_omga = pt.curr_omga; prev_omga = pt.prev_omga;
    curr_acceleration = pt.curr_acceleration; prev_acceleration = pt.prev_acceleration;
    curr_acce_rotate  = pt.curr_acce_rotate;  prev_acce_rotate  = pt.prev_acce_rotate;
    
    pre_force = pt.pre_force; force = pt.force;
    pre_moment = pt.pre_moment; moment = pt.moment;
    
    const_force = pt.const_force; const_moment = pt.const_moment;
    
    flb_force = pt.flb_force; flb_moment = pt.flb_moment; mres = pt.mres;           
    density = pt.density; 
    //mass = pt.mass;
    //REAL volume = pt.volume;
    //vec  J = pt.J;              // moment of inertia in local body-fixed frame

    kinetEnergy = pt.kinetEnergy; // kinetic energy
    cntnum = 0;
    inContact = false; // in contact with other particle or boundary
    isAbleDivide = true;
    candidacy = 1;	/////////// use 2 instead of 1 is to avoid fractured particles to be subdivided again.
    numBroken = pt.numBroken+1;	// it is important that we create a sub-poly-ellipsoid first and then break itself
    typeBroken = pt.typeBroken;	// typeBroken has been determined in calculateBreakPlane()

//    weibullPhi = ran(&idum);
    weibullPhi = pt.weibullPhi;

    matrix zero3x3(3,3);
    average_stress = zero3x3;
    IsFBP = pt.IsFBP;

    // get new center, here the lx,mx,nx and so on are the angles, not the cooridinates of the directional vectors
    REAL lx = curr_direction_a.getx(); REAL mx = curr_direction_a.gety(); REAL nx = curr_direction_a.getz();
    REAL ly = curr_direction_b.getx(); REAL my = curr_direction_b.gety(); REAL ny = curr_direction_b.getz();
    REAL lz = curr_direction_c.getx(); REAL mz = curr_direction_c.gety(); REAL nz = curr_direction_c.getz();

    // calculate the geometries of the sub-divided particles based on oct_num
    REAL shift_ratio = 0.1;
    REAL left_ratio = 0.948683305522212;	// this is calculated as in week summary 2014_09_19_mason_deposition_results_and_fracture_model
    REAL x_oct, y_oct, z_oct;
    switch (break_plane)
    {
	case 1:	// break along ab-plane
	z_oct = -pt.cminus;
	aplus = left_ratio*pt.aplus; aminus = left_ratio*pt.aminus;
	bplus = left_ratio*pt.bplus; bminus = left_ratio*pt.bminus;
	cplus = shift_ratio*pt.cminus; cminus = 0.9*pt.cminus;

	delta_c = -(0.2*cminus-0.8*cplus)/dem::numStepTransition;

    	curr_position.setx(pt.curr_position.getx() + cosl(lz)*shift_ratio*z_oct);
    	curr_position.sety(pt.curr_position.gety() + cosl(mz)*shift_ratio*z_oct);
    	curr_position.setz(pt.curr_position.getz() + cosl(nz)*shift_ratio*z_oct);

	break;

	case 2:	// break along ac-plane
	y_oct = -pt.bminus;
	aplus = left_ratio*pt.aplus; aminus = left_ratio*pt.aminus;
	bplus = shift_ratio*pt.bminus; bminus = 0.9*pt.bminus;
	cplus = left_ratio*pt.cplus; cminus = left_ratio*pt.cminus;

	delta_b = -(0.2*bminus-0.8*bplus)/dem::numStepTransition;

    	curr_position.setx(pt.curr_position.getx() + cosl(ly)*shift_ratio*y_oct);
    	curr_position.sety(pt.curr_position.gety() + cosl(my)*shift_ratio*y_oct);
    	curr_position.setz(pt.curr_position.getz() + cosl(ny)*shift_ratio*y_oct);

	break;

	case 3: // break along bc-plane
	x_oct = -pt.aminus;
	aplus = shift_ratio*pt.aminus; aminus = 0.9*pt.aminus;
	bplus = left_ratio*pt.bplus; bminus = left_ratio*pt.bminus;
	cplus = left_ratio*pt.cplus; cminus = left_ratio*pt.cminus;

	delta_a = -(0.2*aminus-0.8*aplus)/dem::numStepTransition;

    	curr_position.setx(pt.curr_position.getx() + cosl(lx)*shift_ratio*x_oct);
    	curr_position.sety(pt.curr_position.gety() + cosl(mx)*shift_ratio*x_oct);
    	curr_position.setz(pt.curr_position.getz() + cosl(nx)*shift_ratio*x_oct);
		
	break;

	default:
	std::cout << "break plane should be 1/2/3 in particle.cpp..." << std::endl;
	exit(-1);
    }

    // calculate mass, volume, J
   
    // local coordinate of center_geo of every octant
    REAL x1, x2, x3, x4, x5, x6, x7, x8;
    REAL y1, y2, y3, y4, y5, y6, y7, y8;
    REAL z1, z2, z3, z4, z5, z6, z7, z8;
    x1 = 3.0/8.0*aplus; y1 = 3.0/8.0*bplus; z1 = 3.0/8.0*cplus;
    x2 = -3.0/8.0*aminus; y2 = 3.0/8.0*bplus; z2 = 3.0/8.0*cplus;
    x3 = -3.0/8.0*aminus; y3 = -3.0/8.0*bminus; z3 = 3.0/8.0*cplus;
    x4 = 3.0/8.0*aplus; y4 = -3.0/8.0*bminus; z4 = 3.0/8.0*cplus;
    x5 = 3.0/8.0*aplus; y5 = 3.0/8.0*bplus; z5 = -3.0/8.0*cminus;
    x6 = -3.0/8.0*aminus; y6 = 3.0/8.0*bplus; z6 = -3.0/8.0*cminus;
    x7 = -3.0/8.0*aminus; y7 = -3.0/8.0*bminus; z7 = -3.0/8.0*cminus;
    x8 = 3.0/8.0*aplus; y8 = -3.0/8.0*bminus; z8 = -3.0/8.0*cminus;

    REAL v1, v2, v3, v4, v5, v6, v7, v8;	// volumes of eight octants
    v1 = 1.0/6.0*PI*aplus*bplus*cplus;
    v2 = 1.0/6.0*PI*aminus*bplus*cplus;
    v3 = 1.0/6.0*PI*aminus*bminus*cplus;
    v4 = 1.0/6.0*PI*aplus*bminus*cplus;
    v5 = 1.0/6.0*PI*aplus*bplus*cminus;
    v6 = 1.0/6.0*PI*aminus*bplus*cminus;
    v7 = 1.0/6.0*PI*aminus*bminus*cminus;
    v8 = 1.0/6.0*PI*aplus*bminus*cminus;

    volume= v1+v2+v3+v4+v5+v6+v7+v8;
    // local coordinate of center of mass
    REAL xl_center_mass, yl_center_mass, zl_center_mass;
    xl_center_mass = (v1*x1+v2*x2+v3*x3+v4*x4+v5*x5+v6*x6+v7*x7+v8*x8)/volume;
    yl_center_mass = (v1*y1+v2*y2+v3*y3+v4*y4+v5*y5+v6*y6+v7*y7+v8*y8)/volume;
    zl_center_mass = (v1*z1+v2*z2+v3*z3+v4*z4+v5*z5+v6*z6+v7*z7+v8*z8)/volume;
    local_center_mass = vec(xl_center_mass, yl_center_mass, zl_center_mass);

    curr_center_mass = curr_position + globalVec(local_center_mass);

    // prev mass center will be used to determine the next position, which is needed to be changed. 
    // and for each subparticle, we need to figure out the corresponding previous mass center such that the vector pointing from
    // previous mass center to curr mass center of the subparticle should be the same as the vector from previous mass center to 
    // curr mass center of the original particle.
    prev_center_mass = pt.prev_center_mass - pt.curr_center_mass + curr_center_mass;

    mass=density*volume;
    REAL Ixx, Iyy, Izz;
    // moment of inertia with respect to center of geometry
    Ixx = v1*density/5.0*(bplus*bplus+cplus*cplus)+v2*density/5.0*(bplus*bplus+cplus*cplus)
	 +v3*density/5.0*(bminus*bminus+cplus*cplus)+v4*density/5.0*(bminus*bminus+cplus*cplus)
         +v5*density/5.0*(bplus*bplus+cminus*cminus)+v6*density/5.0*(bplus*bplus+cminus*cminus)
	 +v7*density/5.0*(bminus*bminus+cminus*cminus)+v8*density/5.0*(bminus*bminus+cminus*cminus);

    Iyy = v1*density/5.0*(aplus*aplus+cplus*cplus)+v2*density/5.0*(aminus*aminus+cplus*cplus)
	 +v3*density/5.0*(aminus*aminus+cplus*cplus)+v4*density/5.0*(aplus*aplus+cplus*cplus)
         +v5*density/5.0*(aplus*aplus+cminus*cminus)+v6*density/5.0*(aminus*aminus+cminus*cminus)
	 +v7*density/5.0*(aminus*aminus+cminus*cminus)+v8*density/5.0*(aplus*aplus+cminus*cminus);

    Izz = v1*density/5.0*(aplus*aplus+bplus*bplus)+v2*density/5.0*(aminus*aminus+bplus*bplus)
	 +v3*density/5.0*(aminus*aminus+bminus*bminus)+v4*density/5.0*(aplus*aplus+bminus*bminus)
         +v5*density/5.0*(aplus*aplus+bplus*bplus)+v6*density/5.0*(aminus*aminus+bplus*bplus)
	 +v7*density/5.0*(aminus*aminus+bminus*bminus)+v8*density/5.0*(aplus*aplus+bminus*bminus);
    // moment of inertial with respect to center of mass
    Ixx = Ixx+xl_center_mass*xl_center_mass*mass;
    Iyy = Iyy+yl_center_mass*yl_center_mass*mass;
    Izz = Izz+zl_center_mass*zl_center_mass*mass;

    J=vec(Ixx,Iyy,Izz);    

    REAL tmp_parameter = log(pow(1.0/weibullPhi, 1.0/dem::weibullModulus))
			*pow(getAverageRadius()/dem::basicRadius, -2*dem::weibullModulus);
    strengthHoek = dem::sigmaCompress*tmp_parameter;
    strengthContact = dem::ContactTensileCritical*tmp_parameter;

    GlobCoef();	

    // the integrals over poly-ellipsoids for average stress calculation
    int_x  = PI*0.0625*( aplus*aplus*bplus*cplus     - aminus*aminus*bplus*cplus 
			    - aminus*aminus*bminus*cplus  + aplus*aplus*bminus*cplus
			    + aplus*aplus*bplus*cminus    - aminus*aminus*bplus*cminus
			    - aminus*aminus*bminus*cminus + aplus*aplus*bminus*cminus );

    int_y  = PI*0.0625*( aplus*bplus*bplus*cplus     + aminus*bplus*bplus*cplus
			    - aminus*bminus*bminus*cplus  - aplus*bminus*bminus*cplus
			    + aplus*bplus*bplus*cminus    + aminus*bplus*bplus*cminus
			    - aminus*bminus*bminus*cminus - aplus*bminus*bminus*cminus );

    int_z  = PI*0.0625*( aplus*bplus*cplus*cplus     + aminus*bplus*cplus*cplus
			    + aminus*bminus*cplus*cplus   + aplus*bminus*cplus*cplus
			    - aplus*bplus*cminus*cminus   - aminus*bplus*cminus*cminus
			    - aminus*bminus*cminus*cminus - aplus*bminus*cminus*cminus ); 

    int_xy = 1.0/15.0*(  aplus*aplus*bplus*bplus*cplus      - aminus*aminus*bplus*bplus*cplus
			    + aminus*aminus*bminus*bminus*cplus  - aplus*aplus*bminus*bminus*cplus
			    + aplus*aplus*bplus*bplus*cminus     - aminus*aminus*bplus*bplus*cminus
			    + aminus*aminus*bminus*bminus*cminus - aplus*aplus*bminus*bminus*cminus );

    int_xz = 1.0/15.0*(  aplus*aplus*bplus*cplus*cplus      - aminus*aminus*bplus*cplus*cplus
			    - aminus*aminus*bminus*cplus*cplus   + aplus*aplus*bminus*cplus*cplus
			    - aplus*aplus*bplus*cminus*cminus    + aminus*aminus*bplus*cminus*cminus 
			    + aminus*aminus*bminus*cminus*cminus - aplus*aplus*bminus*cminus*cminus );

    int_yz = 1.0/15.0*(  aplus*bplus*bplus*cplus*cplus      + aminus*bplus*bplus*cplus*cplus 
			    - aminus*bminus*bminus*cplus*cplus   - aplus*bminus*bminus*cplus*cplus
			    - aplus*bplus*bplus*cminus*cminus    - aminus*bplus*bplus*cminus*cminus 
			    + aminus*bminus*bminus*cminus*cminus + aplus*bminus*bminus*cminus*cminus );

    int_x2 = PI/30.0*(   aplus*aplus*aplus*bplus*cplus      + aminus*aminus*aminus*bplus*cplus
			    + aminus*aminus*aminus*bminus*cplus  + aplus*aplus*aplus*bminus*cplus
			    + aplus*aplus*aplus*bplus*cminus     + aminus*aminus*aminus*bplus*cminus
			    + aminus*aminus*aminus*bminus*cminus + aplus*aplus*aplus*bminus*cminus );

    int_y2 = PI/30.0*(   aplus*bplus*bplus*bplus*cplus      + aminus*bplus*bplus*bplus*cplus
			    + aminus*bminus*bminus*bminus*cplus  + aplus*bminus*bminus*bminus*cplus
			    + aplus*bplus*bplus*bplus*cminus     + aminus*bplus*bplus*bplus*cminus
			    + aminus*bminus*bminus*bminus*cminus + aplus*bminus*bminus*bminus*cminus );

    int_z2 = PI/30.0*(   aplus*bplus*cplus*cplus*cplus      + aminus*bplus*cplus*cplus*cplus
			    + aminus*bminus*cplus*cplus*cplus    + aplus*bminus*cplus*cplus*cplus
			    + aplus*bplus*cminus*cminus*cminus   + aminus*bplus*cminus*cminus*cminus
			    + aminus*bminus*cminus*cminus*cminus + aplus*bminus*cminus*cminus*cminus );
}


// change the particle itself to the fractured particle which is breaked by the break_plane of the original particle
// and also in the positive part, September 22, 2014. 
void particle::breakItSelf(int break_plane){	

    cntnum = 0;
    inContact = false; // in contact with other particle or boundary
    isAbleDivide = true;
    candidacy = 1;	/////////// use 2 instead of 1 is to not want fractured particles to be subdivided again.
    numBroken++;

    matrix zero3x3(3,3);
    average_stress = zero3x3;

    // get new center, here the lx,mx,nx and so on are the angles, not the cooridinates of the directional vectors
    REAL lx = curr_direction_a.getx(); REAL mx = curr_direction_a.gety(); REAL nx = curr_direction_a.getz();
    REAL ly = curr_direction_b.getx(); REAL my = curr_direction_b.gety(); REAL ny = curr_direction_b.getz();
    REAL lz = curr_direction_c.getx(); REAL mz = curr_direction_c.gety(); REAL nz = curr_direction_c.getz();

    // calculate the geometries of the sub-divided particles based on oct_num
    REAL shift_ratio = 0.1;
    REAL left_ratio = 0.948683305522212;	// this is calculated as in week summary 2014_09_19_mason_deposition_results_and_fracture_model
    REAL x_oct, y_oct, z_oct;
    switch (break_plane)
    {
	case 1:	// break along ab-plane
	z_oct = cplus;
	aplus = left_ratio*aplus; aminus = left_ratio*aminus;
	bplus = left_ratio*bplus; bminus = left_ratio*bminus;
	cminus = shift_ratio*cplus; cplus = 0.9*cplus;	// cminus has to be before cplus

	delta_c = (0.2*cplus-0.8*cminus)/dem::numStepTransition;

    	curr_position.setx(curr_position.getx() + cosl(lz)*shift_ratio*z_oct);
    	curr_position.sety(curr_position.gety() + cosl(mz)*shift_ratio*z_oct);
    	curr_position.setz(curr_position.getz() + cosl(nz)*shift_ratio*z_oct);

	break;

	case 2:	// break along ac-plane
	y_oct = bplus;
	aplus = left_ratio*aplus; aminus = left_ratio*aminus;
	bminus = shift_ratio*bplus; bplus = 0.9*bplus;	// bminus has to be before bplus
	cplus = left_ratio*cplus; cminus = left_ratio*cminus;

	delta_b = (0.2*bplus-0.8*bminus)/dem::numStepTransition;

    	curr_position.setx(curr_position.getx() + cosl(ly)*shift_ratio*y_oct);
    	curr_position.sety(curr_position.gety() + cosl(my)*shift_ratio*y_oct);
    	curr_position.setz(curr_position.getz() + cosl(ny)*shift_ratio*y_oct);

	break;

	case 3: // break along bc-plane
	x_oct = aplus;
	aminus = shift_ratio*aplus; aplus = 0.9*aplus; 	// aminus has to be before aplus
	bplus = left_ratio*bplus; bminus = left_ratio*bminus;
	cplus = left_ratio*cplus; cminus = left_ratio*cminus;

	delta_a = (0.2*aplus-0.8*aminus)/dem::numStepTransition;

    	curr_position.setx(curr_position.getx() + cosl(lx)*shift_ratio*x_oct);
    	curr_position.sety(curr_position.gety() + cosl(mx)*shift_ratio*x_oct);
    	curr_position.setz(curr_position.getz() + cosl(nx)*shift_ratio*x_oct);
		
	break;

	default:
	std::cout << "break plane should be 1/2/3 in particle.cpp..." << std::endl;
	exit(-1);
    }

    // calculate mass, volume, J
   
    // local coordinate of center_geo of every octant
    REAL x1, x2, x3, x4, x5, x6, x7, x8;
    REAL y1, y2, y3, y4, y5, y6, y7, y8;
    REAL z1, z2, z3, z4, z5, z6, z7, z8;
    x1 = 3.0/8.0*aplus; y1 = 3.0/8.0*bplus; z1 = 3.0/8.0*cplus;
    x2 = -3.0/8.0*aminus; y2 = 3.0/8.0*bplus; z2 = 3.0/8.0*cplus;
    x3 = -3.0/8.0*aminus; y3 = -3.0/8.0*bminus; z3 = 3.0/8.0*cplus;
    x4 = 3.0/8.0*aplus; y4 = -3.0/8.0*bminus; z4 = 3.0/8.0*cplus;
    x5 = 3.0/8.0*aplus; y5 = 3.0/8.0*bplus; z5 = -3.0/8.0*cminus;
    x6 = -3.0/8.0*aminus; y6 = 3.0/8.0*bplus; z6 = -3.0/8.0*cminus;
    x7 = -3.0/8.0*aminus; y7 = -3.0/8.0*bminus; z7 = -3.0/8.0*cminus;
    x8 = 3.0/8.0*aplus; y8 = -3.0/8.0*bminus; z8 = -3.0/8.0*cminus;

    REAL v1, v2, v3, v4, v5, v6, v7, v8;	// volumes of eight octants
    v1 = 1.0/6.0*PI*aplus*bplus*cplus;
    v2 = 1.0/6.0*PI*aminus*bplus*cplus;
    v3 = 1.0/6.0*PI*aminus*bminus*cplus;
    v4 = 1.0/6.0*PI*aplus*bminus*cplus;
    v5 = 1.0/6.0*PI*aplus*bplus*cminus;
    v6 = 1.0/6.0*PI*aminus*bplus*cminus;
    v7 = 1.0/6.0*PI*aminus*bminus*cminus;
    v8 = 1.0/6.0*PI*aplus*bminus*cminus;

    volume= v1+v2+v3+v4+v5+v6+v7+v8;
    // local coordinate of center of mass
    REAL xl_center_mass, yl_center_mass, zl_center_mass;
    xl_center_mass = (v1*x1+v2*x2+v3*x3+v4*x4+v5*x5+v6*x6+v7*x7+v8*x8)/volume;
    yl_center_mass = (v1*y1+v2*y2+v3*y3+v4*y4+v5*y5+v6*y6+v7*y7+v8*y8)/volume;
    zl_center_mass = (v1*z1+v2*z2+v3*z3+v4*z4+v5*z5+v6*z6+v7*z7+v8*z8)/volume;
    local_center_mass = vec(xl_center_mass, yl_center_mass, zl_center_mass);

    prev_center_mass = prev_center_mass - curr_center_mass;
    curr_center_mass = curr_position + globalVec(local_center_mass);

    // prev mass center will be used to determine the next position, which is needed to be changed. 
    // and for each subparticle, we need to figure out the corresponding previous mass center such that the vector pointing from
    // previous mass center to curr mass center of the subparticle should be the same as the vector from previous mass center to 
    // curr mass center of the original particle.
    prev_center_mass += curr_center_mass;

    mass=density*volume;
    REAL Ixx, Iyy, Izz;
    // moment of inertia with respect to center of geometry
    Ixx = v1*density/5.0*(bplus*bplus+cplus*cplus)+v2*density/5.0*(bplus*bplus+cplus*cplus)
	 +v3*density/5.0*(bminus*bminus+cplus*cplus)+v4*density/5.0*(bminus*bminus+cplus*cplus)
         +v5*density/5.0*(bplus*bplus+cminus*cminus)+v6*density/5.0*(bplus*bplus+cminus*cminus)
	 +v7*density/5.0*(bminus*bminus+cminus*cminus)+v8*density/5.0*(bminus*bminus+cminus*cminus);

    Iyy = v1*density/5.0*(aplus*aplus+cplus*cplus)+v2*density/5.0*(aminus*aminus+cplus*cplus)
	 +v3*density/5.0*(aminus*aminus+cplus*cplus)+v4*density/5.0*(aplus*aplus+cplus*cplus)
         +v5*density/5.0*(aplus*aplus+cminus*cminus)+v6*density/5.0*(aminus*aminus+cminus*cminus)
	 +v7*density/5.0*(aminus*aminus+cminus*cminus)+v8*density/5.0*(aplus*aplus+cminus*cminus);

    Izz = v1*density/5.0*(aplus*aplus+bplus*bplus)+v2*density/5.0*(aminus*aminus+bplus*bplus)
	 +v3*density/5.0*(aminus*aminus+bminus*bminus)+v4*density/5.0*(aplus*aplus+bminus*bminus)
         +v5*density/5.0*(aplus*aplus+bplus*bplus)+v6*density/5.0*(aminus*aminus+bplus*bplus)
	 +v7*density/5.0*(aminus*aminus+bminus*bminus)+v8*density/5.0*(aplus*aplus+bminus*bminus);
    // moment of inertial with respect to center of mass
    Ixx = Ixx+xl_center_mass*xl_center_mass*mass;
    Iyy = Iyy+yl_center_mass*yl_center_mass*mass;
    Izz = Izz+zl_center_mass*zl_center_mass*mass;

    J=vec(Ixx,Iyy,Izz);    

//    weibullPhi = ran(&idum);
    REAL tmp_parameter = log(pow(1.0/weibullPhi, 1.0/dem::weibullModulus))
			*pow(getAverageRadius()/dem::basicRadius, -2*dem::weibullModulus);
    strengthHoek = dem::sigmaCompress*tmp_parameter;
    strengthContact = dem::ContactTensileCritical*tmp_parameter;

    GlobCoef();	

    // the integrals over poly-ellipsoids for average stress calculation
    int_x  = PI*0.0625*( aplus*aplus*bplus*cplus     - aminus*aminus*bplus*cplus 
			    - aminus*aminus*bminus*cplus  + aplus*aplus*bminus*cplus
			    + aplus*aplus*bplus*cminus    - aminus*aminus*bplus*cminus
			    - aminus*aminus*bminus*cminus + aplus*aplus*bminus*cminus );

    int_y  = PI*0.0625*( aplus*bplus*bplus*cplus     + aminus*bplus*bplus*cplus
			    - aminus*bminus*bminus*cplus  - aplus*bminus*bminus*cplus
			    + aplus*bplus*bplus*cminus    + aminus*bplus*bplus*cminus
			    - aminus*bminus*bminus*cminus - aplus*bminus*bminus*cminus );

    int_z  = PI*0.0625*( aplus*bplus*cplus*cplus     + aminus*bplus*cplus*cplus
			    + aminus*bminus*cplus*cplus   + aplus*bminus*cplus*cplus
			    - aplus*bplus*cminus*cminus   - aminus*bplus*cminus*cminus
			    - aminus*bminus*cminus*cminus - aplus*bminus*cminus*cminus ); 

    int_xy = 1.0/15.0*(  aplus*aplus*bplus*bplus*cplus      - aminus*aminus*bplus*bplus*cplus
			    + aminus*aminus*bminus*bminus*cplus  - aplus*aplus*bminus*bminus*cplus
			    + aplus*aplus*bplus*bplus*cminus     - aminus*aminus*bplus*bplus*cminus
			    + aminus*aminus*bminus*bminus*cminus - aplus*aplus*bminus*bminus*cminus );

    int_xz = 1.0/15.0*(  aplus*aplus*bplus*cplus*cplus      - aminus*aminus*bplus*cplus*cplus
			    - aminus*aminus*bminus*cplus*cplus   + aplus*aplus*bminus*cplus*cplus
			    - aplus*aplus*bplus*cminus*cminus    + aminus*aminus*bplus*cminus*cminus 
			    + aminus*aminus*bminus*cminus*cminus - aplus*aplus*bminus*cminus*cminus );

    int_yz = 1.0/15.0*(  aplus*bplus*bplus*cplus*cplus      + aminus*bplus*bplus*cplus*cplus 
			    - aminus*bminus*bminus*cplus*cplus   - aplus*bminus*bminus*cplus*cplus
			    - aplus*bplus*bplus*cminus*cminus    - aminus*bplus*bplus*cminus*cminus 
			    + aminus*bminus*bminus*cminus*cminus + aplus*bminus*bminus*cminus*cminus );

    int_x2 = PI/30.0*(   aplus*aplus*aplus*bplus*cplus      + aminus*aminus*aminus*bplus*cplus
			    + aminus*aminus*aminus*bminus*cplus  + aplus*aplus*aplus*bminus*cplus
			    + aplus*aplus*aplus*bplus*cminus     + aminus*aminus*aminus*bplus*cminus
			    + aminus*aminus*aminus*bminus*cminus + aplus*aplus*aplus*bminus*cminus );

    int_y2 = PI/30.0*(   aplus*bplus*bplus*bplus*cplus      + aminus*bplus*bplus*bplus*cplus
			    + aminus*bminus*bminus*bminus*cplus  + aplus*bminus*bminus*bminus*cplus
			    + aplus*bplus*bplus*bplus*cminus     + aminus*bplus*bplus*bplus*cminus
			    + aminus*bminus*bminus*bminus*cminus + aplus*bminus*bminus*bminus*cminus );

    int_z2 = PI/30.0*(   aplus*bplus*cplus*cplus*cplus      + aminus*bplus*cplus*cplus*cplus
			    + aminus*bminus*cplus*cplus*cplus    + aplus*bminus*cplus*cplus*cplus
			    + aplus*bplus*cminus*cminus*cminus   + aminus*bplus*cminus*cminus*cminus
			    + aminus*bminus*cminus*cminus*cminus + aplus*bminus*cminus*cminus*cminus );

} // end breakItSelf()


// transit bad-shaped poly-ellipsoids generated in fracture
// to normal shaped poly-ellipsoids, suggested by Dr. Fu at LLNL
// make sure this transition (movement of centroid) will not introduce 
// dditional velocity to the system
void particle::shapeTransition(){

  bool isTransit = false;

  if(delta_a!=0){	// need transition
    // shape transition
    // new centroid
    curr_position.setx(curr_position.getx() + cosl(curr_direction_a.getx())*delta_a);
    curr_position.sety(curr_position.gety() + cosl(curr_direction_a.gety())*delta_a);
    curr_position.setz(curr_position.getz() + cosl(curr_direction_a.getz())*delta_a);
    // new aplus and aminus, b and c will be the same
    aplus = aplus-delta_a; aminus = aminus+delta_a;	// if move toward aplus, then delta_a itself is positive
    isTransit = true;
    // check
    if(aplus>0.25*aminus && aminus>0.25*aplus)	// in normal shape
	delta_a = 0;
  } // end if delta_a

  if(delta_b!=0){	// need transition
    // shape transition
    // new centroid
    curr_position.setx(curr_position.getx() + cosl(curr_direction_b.getx())*delta_b);
    curr_position.sety(curr_position.gety() + cosl(curr_direction_b.gety())*delta_b);
    curr_position.setz(curr_position.getz() + cosl(curr_direction_b.getz())*delta_b);
    // new bplus and bminus, a and c will be the same
    bplus = bplus-delta_b; bminus = bminus+delta_b;	// if move toward aplus, then delta_a itself is positive
    isTransit = true;
    // check
    if(bplus>0.25*bminus && bminus>0.25*bplus)	// in normal shape
	delta_b = 0;
  } // end if delta_b

  if(delta_c!=0){	// need transition
    // shape transition
    // new centroid
    curr_position.setx(curr_position.getx() + cosl(curr_direction_c.getx())*delta_c);
    curr_position.sety(curr_position.gety() + cosl(curr_direction_c.gety())*delta_c);
    curr_position.setz(curr_position.getz() + cosl(curr_direction_c.getz())*delta_c);
    // new cplus and cminus, a and b will be the same
    cplus = cplus-delta_c; cminus = cminus+delta_c;	// if move toward aplus, then delta_a itself is positive
    isTransit = true;
    // check
    if(cplus>0.25*cminus && cminus>0.25*cplus)	// in normal shape
	delta_c = 0;
  } // end if delta_c

  if(isTransit){
    calculateGeometry();	// calculate geometry of the new poly-ellipsoid
  }

} // shapeTransition()


void particle::calculateGeometry(){
    // calculate mass, volume, J
   
    // local coordinate of center_geo of every octant
    REAL x1, x2, x3, x4, x5, x6, x7, x8;
    REAL y1, y2, y3, y4, y5, y6, y7, y8;
    REAL z1, z2, z3, z4, z5, z6, z7, z8;
    x1 = 3.0/8.0*aplus; y1 = 3.0/8.0*bplus; z1 = 3.0/8.0*cplus;
    x2 = -3.0/8.0*aminus; y2 = 3.0/8.0*bplus; z2 = 3.0/8.0*cplus;
    x3 = -3.0/8.0*aminus; y3 = -3.0/8.0*bminus; z3 = 3.0/8.0*cplus;
    x4 = 3.0/8.0*aplus; y4 = -3.0/8.0*bminus; z4 = 3.0/8.0*cplus;
    x5 = 3.0/8.0*aplus; y5 = 3.0/8.0*bplus; z5 = -3.0/8.0*cminus;
    x6 = -3.0/8.0*aminus; y6 = 3.0/8.0*bplus; z6 = -3.0/8.0*cminus;
    x7 = -3.0/8.0*aminus; y7 = -3.0/8.0*bminus; z7 = -3.0/8.0*cminus;
    x8 = 3.0/8.0*aplus; y8 = -3.0/8.0*bminus; z8 = -3.0/8.0*cminus;

    REAL v1, v2, v3, v4, v5, v6, v7, v8;	// volumes of eight octants
    v1 = 1.0/6.0*PI*aplus*bplus*cplus;
    v2 = 1.0/6.0*PI*aminus*bplus*cplus;
    v3 = 1.0/6.0*PI*aminus*bminus*cplus;
    v4 = 1.0/6.0*PI*aplus*bminus*cplus;
    v5 = 1.0/6.0*PI*aplus*bplus*cminus;
    v6 = 1.0/6.0*PI*aminus*bplus*cminus;
    v7 = 1.0/6.0*PI*aminus*bminus*cminus;
    v8 = 1.0/6.0*PI*aplus*bminus*cminus;

    volume= v1+v2+v3+v4+v5+v6+v7+v8;
    // local coordinate of center of mass
    REAL xl_center_mass, yl_center_mass, zl_center_mass;
    xl_center_mass = (v1*x1+v2*x2+v3*x3+v4*x4+v5*x5+v6*x6+v7*x7+v8*x8)/volume;
    yl_center_mass = (v1*y1+v2*y2+v3*y3+v4*y4+v5*y5+v6*y6+v7*y7+v8*y8)/volume;
    zl_center_mass = (v1*z1+v2*z2+v3*z3+v4*z4+v5*z5+v6*z6+v7*z7+v8*z8)/volume;
    local_center_mass = vec(xl_center_mass, yl_center_mass, zl_center_mass);

    prev_center_mass = prev_center_mass - curr_center_mass;
    curr_center_mass = curr_position + globalVec(local_center_mass);

    // prev mass center will be used to determine the next position, which is needed to be changed. 
    // and for each subparticle, we need to figure out the corresponding previous mass center such that the vector pointing from
    // previous mass center to curr mass center of the subparticle should be the same as the vector from previous mass center to 
    // curr mass center of the original particle.
    prev_center_mass += curr_center_mass;

    mass=density*volume;
    REAL Ixx, Iyy, Izz;
    // moment of inertia with respect to center of geometry
    Ixx = v1*density/5.0*(bplus*bplus+cplus*cplus)+v2*density/5.0*(bplus*bplus+cplus*cplus)
	 +v3*density/5.0*(bminus*bminus+cplus*cplus)+v4*density/5.0*(bminus*bminus+cplus*cplus)
         +v5*density/5.0*(bplus*bplus+cminus*cminus)+v6*density/5.0*(bplus*bplus+cminus*cminus)
	 +v7*density/5.0*(bminus*bminus+cminus*cminus)+v8*density/5.0*(bminus*bminus+cminus*cminus);

    Iyy = v1*density/5.0*(aplus*aplus+cplus*cplus)+v2*density/5.0*(aminus*aminus+cplus*cplus)
	 +v3*density/5.0*(aminus*aminus+cplus*cplus)+v4*density/5.0*(aplus*aplus+cplus*cplus)
         +v5*density/5.0*(aplus*aplus+cminus*cminus)+v6*density/5.0*(aminus*aminus+cminus*cminus)
	 +v7*density/5.0*(aminus*aminus+cminus*cminus)+v8*density/5.0*(aplus*aplus+cminus*cminus);

    Izz = v1*density/5.0*(aplus*aplus+bplus*bplus)+v2*density/5.0*(aminus*aminus+bplus*bplus)
	 +v3*density/5.0*(aminus*aminus+bminus*bminus)+v4*density/5.0*(aplus*aplus+bminus*bminus)
         +v5*density/5.0*(aplus*aplus+bplus*bplus)+v6*density/5.0*(aminus*aminus+bplus*bplus)
	 +v7*density/5.0*(aminus*aminus+bminus*bminus)+v8*density/5.0*(aplus*aplus+bminus*bminus);
    // moment of inertial with respect to center of mass
    Ixx = Ixx+xl_center_mass*xl_center_mass*mass;
    Iyy = Iyy+yl_center_mass*yl_center_mass*mass;
    Izz = Izz+zl_center_mass*zl_center_mass*mass;

    J=vec(Ixx,Iyy,Izz);    

    GlobCoef();	

} // calculateGeometry()


REAL particle::getMaxRadius() const{	// August 21, 2013
	REAL temp;
	temp = aplus;
	if(temp < aminus)
		temp = aminus;
	if(temp < bplus)
		temp = bplus;
	if(temp < bminus)
		temp = bminus;
	if(temp < cplus)
		temp = cplus;
	if(temp < cminus)
		temp = cminus;

	return temp;
}

REAL particle::getMinRadius() const{	// August 27, 2013
	REAL temp;
	temp = aplus;
	if(temp > aminus)
		temp = aminus;
	if(temp > bplus)
		temp = bplus;
	if(temp > bminus)
		temp = bminus;
	if(temp > cplus)
		temp = cplus;
	if(temp > cminus)
		temp = cminus;

	return temp;
}

// 1: rotational energy is 1/2(I1*w1^2+I2*w2^2+I3*w3^2), where each term
//    is expressed in local frame.
// 2. angular velocities in global frame needs to be converted to those
//    in local frame.
REAL particle::getTransEnergy() const{
    return mass*pow(vfabs(curr_velocity),2)/2;
}


REAL particle::getRotatEnergy() const{
    vec curr_local_omga, tmp;

    tmp=vcos(curr_direction_a); curr_local_omga.setx(tmp%curr_omga);
    tmp=vcos(curr_direction_b); curr_local_omga.sety(tmp%curr_omga);
    tmp=vcos(curr_direction_c); curr_local_omga.setz(tmp%curr_omga);

    return J.getx()*pow(curr_local_omga.getx(),2)/2 +
	   J.gety()*pow(curr_local_omga.gety(),2)/2 +
	   J.getz()*pow(curr_local_omga.getz(),2)/2;
}


REAL particle::getKinetEnergy() const{
    return getTransEnergy() + getRotatEnergy();
}


REAL particle::getPotenEnergy(REAL ref) const{	// August 16, 2013
    return G*mass*(curr_center_mass.getz() - ref);
}


void   particle::getGlobCoef(REAL coef[], int num_oct) const{	// September 6, 2013
    switch (num_oct){
    	case 1:
	    for (int i=0;i<10;i++)
                coef[i]=this->coef1[i];
	    break;
    	case 2:
	    for (int i=0;i<10;i++)
                coef[i]=this->coef2[i];
	    break;
    	case 3:
	    for (int i=0;i<10;i++)
                coef[i]=this->coef3[i];
	    break;
    	case 4:
	    for (int i=0;i<10;i++)
                coef[i]=this->coef4[i];
	    break;
    	case 5:
	    for (int i=0;i<10;i++)
                coef[i]=this->coef5[i];
	    break;
    	case 6:
	    for (int i=0;i<10;i++)
                coef[i]=this->coef6[i];
	    break;
    	case 7:
	    for (int i=0;i<10;i++)
                coef[i]=this->coef7[i];
	    break;
    	case 8:
	    for (int i=0;i<10;i++)
                coef[i]=this->coef8[i];
	    break;
	default:
	    std::cout << "number of octant is larger then 8 in getGlobCoef()!" << std::cout;
	    exit(-1);
	    break;
    }
}


void particle::print() const{	//August 16, 2013
         std::cout<<"aplus="<<aplus<<'\t'<<"aminus="<<aminus<<'\t'<<"bplus="<<bplus<<'\t'<<"bminus="<<bminus<<'\t'<<"cplus="<<cplus<<'\t'<<"cminus="<<cminus<<std::endl;
    std::cout<<"curr_direction_a=";
    vcos(curr_direction_a).print();
    std::cout<<"curr_direction_b=";
    vcos(curr_direction_b).print();
    std::cout<<"curr_direction_c=";
    vcos(curr_direction_c).print();
    std::cout<<"curr_position=";
    curr_position.print();
}


REAL particle::surfaceError(vec pt) const{
    REAL x=pt.getx();
    REAL y=pt.gety();
    REAL z=pt.getz();
    return coef1[0]*x*x+coef1[1]*y*y+coef1[2]*z*z+coef1[3]*x*y+coef1[4]*y*z+coef1[5]*z*x+
	coef1[6]*x+coef1[7]*y+coef1[8]*z+coef1[9];
}


//void particle::GlobCoef(int asign, int bsign, int csign){	// August 21, 2013. September 6, 2013
void particle::GlobCoef(){	// September 6, 2013
    //coef[0]--x^2, coef[1]--y^2, coef[2]--z^2, coef[3]--xy, coef[4]--yz, coef[5]--xz
    //coef[6]--x, coef[7]--y, coef[8]--z, coef[9]--const
    REAL a, b, c;
/*    switch (asign){
	case 1:
	   a = aplus;
	   break;
	case -1:
	   a = aminus;
	   break;
	default:	// approximation
	   a = (aplus+aminus)*0.5;
	   break;
    }
    switch (bsign){
	case 1:
	   b = bplus;
	   break;
	case -1:
	   b = bminus;
	   break;
	default:	// approximation
	   b = (bplus+bminus)*0.5;
	   break;
    }
    switch (csign){
	case 1:
	   c = cplus;
	   break;
	case -1:
	   c = cminus;
	   break;
	default:	// approximation
	   c = (cplus+cminus)*0.5;
	   break;
    }
*/
    vec v1=vcos(curr_direction_a);
    vec v2=vcos(curr_direction_b);
    vec v3=vcos(curr_direction_c);
    REAL X0=curr_position.getx();
    REAL Y0=curr_position.gety();
    REAL Z0=curr_position.getz();
    REAL l1=v1.getx();
    REAL m1=v1.gety();
    REAL n1=v1.getz();
    REAL l2=v2.getx();
    REAL m2=v2.gety();
    REAL n2=v2.getz();
    REAL l3=v3.getx();
    REAL m3=v3.gety();
    REAL n3=v3.getz();

    REAL divd;


    // the first octant
    a=aplus; b=bplus; c=cplus;
    if(a==b&&b==c){	// this octant is a sphere
	coef1[0]=1;
	coef1[1]=1;
	coef1[2]=1;
	coef1[3]=0;
	coef1[4]=0;
	coef1[5]=0;
	coef1[6]=-2*curr_position.getx();
	coef1[7]=-2*curr_position.gety();
	coef1[8]=-2*curr_position.getz();
	coef1[9]=pow(vfabs(curr_position),2)-a*a;
    }
    else{
    	coef1[0]=l1*l1/a/a+l2*l2/b/b+l3*l3/c/c;
    	coef1[1]=m1*m1/a/a+m2*m2/b/b+m3*m3/c/c;
    	coef1[2]=n1*n1/a/a+n2*n2/b/b+n3*n3/c/c;
    	coef1[3]=(2*l1*m1)/a/a + (2*l2*m2)/b/b + (2*l3*m3)/c/c;
    	coef1[4]=(2*m1*n1)/a/a + (2*m2*n2)/b/b + (2*m3*n3)/c/c;
    	coef1[5]=(2*l1*n1)/a/a + (2*l2*n2)/b/b + (2*l3*n3)/c/c;
    	coef1[6]=
		-2*l1*m1*Y0*pow(a,-2) - 2*l1*n1*Z0*pow(a,-2) - 
		2*l2*m2*Y0*pow(b,-2) - 2*l2*n2*Z0*pow(b,-2) - 
		2*l3*m3*Y0*pow(c,-2) - 2*l3*n3*Z0*pow(c,-2) - 
		2*X0*pow(a,-2)*pow(l1,2) - 2*X0*pow(b,-2)*pow(l2,2) - 
		2*X0*pow(c,-2)*pow(l3,2);
    	coef1[7]=
		(-2*l1*m1*X0)/a/a - (2*l2*m2*X0)/b/b - 
		(2*l3*m3*X0)/c/c - (2*m1*m1*Y0)/a/a - 
		(2*m2*m2*Y0)/b/b - (2*m3*m3*Y0)/c/c - 
		(2*m1*n1*Z0)/a/a - (2*m2*n2*Z0)/b/b - 
		(2*m3*n3*Z0)/c/c;
    	coef1[8]=
		(-2*l1*n1*X0)/a/a - (2*l2*n2*X0)/b/b - 
		(2*l3*n3*X0)/c/c - (2*m1*n1*Y0)/a/a - 
		(2*m2*n2*Y0)/b/b - (2*m3*n3*Y0)/c/c - 
		(2*n1*n1*Z0)/a/a - (2*n2*n2*Z0)/b/b - 
		(2*n3*n3*Z0)/c/c;
    	coef1[9]=
		-1 + 2*l1*m1*X0*Y0*pow(a,-2) + 2*l1*n1*X0*Z0*pow(a,-2) + 
		2*m1*n1*Y0*Z0*pow(a,-2) + 2*l2*m2*X0*Y0*pow(b,-2) + 
		2*l2*n2*X0*Z0*pow(b,-2) + 2*m2*n2*Y0*Z0*pow(b,-2) + 
		2*l3*m3*X0*Y0*pow(c,-2) + 2*l3*n3*X0*Z0*pow(c,-2) + 
		2*m3*n3*Y0*Z0*pow(c,-2) + 
		pow(a,-2)*pow(l1,2)*pow(X0,2) + 
		pow(b,-2)*pow(l2,2)*pow(X0,2) + 
		pow(c,-2)*pow(l3,2)*pow(X0,2) + 
		pow(a,-2)*pow(m1,2)*pow(Y0,2) + 
		pow(b,-2)*pow(m2,2)*pow(Y0,2) + 
		pow(c,-2)*pow(m3,2)*pow(Y0,2) + 
		pow(a,-2)*pow(n1,2)*pow(Z0,2) + 
		pow(b,-2)*pow(n2,2)*pow(Z0,2) + 
		pow(c,-2)*pow(n3,2)*pow(Z0,2);
    	divd=coef1[0];
    	for (int kk=0;kk<10;kk++){  // when a particle is initialized or updated, coef[0] is set as 1.0.
	    coef1[kk]/=divd;
    	}
    } // end if

    // the second octant
    a=aminus; b=bplus; c=cplus;
    if(a==b&&b==c){	// this octant is a sphere
	coef2[0]=1;
	coef2[1]=1;
	coef2[2]=1;
	coef2[3]=0;
	coef2[4]=0;
	coef2[5]=0;
	coef2[6]=-2*curr_position.getx();
	coef2[7]=-2*curr_position.gety();
	coef2[8]=-2*curr_position.getz();
	coef2[9]=pow(vfabs(curr_position),2)-a*a;
    }
    else{
    	coef2[0]=l1*l1/a/a+l2*l2/b/b+l3*l3/c/c;
    	coef2[1]=m1*m1/a/a+m2*m2/b/b+m3*m3/c/c;
    	coef2[2]=n1*n1/a/a+n2*n2/b/b+n3*n3/c/c;
    	coef2[3]=(2*l1*m1)/a/a + (2*l2*m2)/b/b + (2*l3*m3)/c/c;
    	coef2[4]=(2*m1*n1)/a/a + (2*m2*n2)/b/b + (2*m3*n3)/c/c;
    	coef2[5]=(2*l1*n1)/a/a + (2*l2*n2)/b/b + (2*l3*n3)/c/c;
    	coef2[6]=
		-2*l1*m1*Y0*pow(a,-2) - 2*l1*n1*Z0*pow(a,-2) - 
		2*l2*m2*Y0*pow(b,-2) - 2*l2*n2*Z0*pow(b,-2) - 
		2*l3*m3*Y0*pow(c,-2) - 2*l3*n3*Z0*pow(c,-2) - 
		2*X0*pow(a,-2)*pow(l1,2) - 2*X0*pow(b,-2)*pow(l2,2) - 
		2*X0*pow(c,-2)*pow(l3,2);
    	coef2[7]=
		(-2*l1*m1*X0)/a/a - (2*l2*m2*X0)/b/b - 
		(2*l3*m3*X0)/c/c - (2*m1*m1*Y0)/a/a - 
		(2*m2*m2*Y0)/b/b - (2*m3*m3*Y0)/c/c - 
		(2*m1*n1*Z0)/a/a - (2*m2*n2*Z0)/b/b - 
		(2*m3*n3*Z0)/c/c;
	coef2[8]=
		(-2*l1*n1*X0)/a/a - (2*l2*n2*X0)/b/b - 
		(2*l3*n3*X0)/c/c - (2*m1*n1*Y0)/a/a - 
		(2*m2*n2*Y0)/b/b - (2*m3*n3*Y0)/c/c - 
		(2*n1*n1*Z0)/a/a - (2*n2*n2*Z0)/b/b - 
		(2*n3*n3*Z0)/c/c;
    	coef2[9]=
		-1 + 2*l1*m1*X0*Y0*pow(a,-2) + 2*l1*n1*X0*Z0*pow(a,-2) + 
		2*m1*n1*Y0*Z0*pow(a,-2) + 2*l2*m2*X0*Y0*pow(b,-2) + 
		2*l2*n2*X0*Z0*pow(b,-2) + 2*m2*n2*Y0*Z0*pow(b,-2) + 
		2*l3*m3*X0*Y0*pow(c,-2) + 2*l3*n3*X0*Z0*pow(c,-2) + 
		2*m3*n3*Y0*Z0*pow(c,-2) + 
		pow(a,-2)*pow(l1,2)*pow(X0,2) + 
		pow(b,-2)*pow(l2,2)*pow(X0,2) + 
		pow(c,-2)*pow(l3,2)*pow(X0,2) + 
		pow(a,-2)*pow(m1,2)*pow(Y0,2) + 
		pow(b,-2)*pow(m2,2)*pow(Y0,2) + 
		pow(c,-2)*pow(m3,2)*pow(Y0,2) + 
		pow(a,-2)*pow(n1,2)*pow(Z0,2) + 
		pow(b,-2)*pow(n2,2)*pow(Z0,2) + 
		pow(c,-2)*pow(n3,2)*pow(Z0,2);
    	divd=coef2[0];
    	for (int kk=0;kk<10;kk++){  // when a particle is initialized or updated, coef[0] is set as 1.0.
	    coef2[kk]/=divd;
    	}
    } // end if

    // the third octant
    a=aminus; b=bminus; c=cplus;
    if(a==b&&b==c){	// this octant is a sphere
	coef3[0]=1;
	coef3[1]=1;
	coef3[2]=1;
	coef3[3]=0;
	coef3[4]=0;
	coef3[5]=0;
	coef3[6]=-2*curr_position.getx();
	coef3[7]=-2*curr_position.gety();
	coef3[8]=-2*curr_position.getz();
	coef3[9]=pow(vfabs(curr_position),2)-a*a;
    }
    else{
    	coef3[0]=l1*l1/a/a+l2*l2/b/b+l3*l3/c/c;
    	coef3[1]=m1*m1/a/a+m2*m2/b/b+m3*m3/c/c;
    	coef3[2]=n1*n1/a/a+n2*n2/b/b+n3*n3/c/c;
    	coef3[3]=(2*l1*m1)/a/a + (2*l2*m2)/b/b + (2*l3*m3)/c/c;
    	coef3[4]=(2*m1*n1)/a/a + (2*m2*n2)/b/b + (2*m3*n3)/c/c;
    	coef3[5]=(2*l1*n1)/a/a + (2*l2*n2)/b/b + (2*l3*n3)/c/c;
    	coef3[6]=
		-2*l1*m1*Y0*pow(a,-2) - 2*l1*n1*Z0*pow(a,-2) - 
		2*l2*m2*Y0*pow(b,-2) - 2*l2*n2*Z0*pow(b,-2) - 
		2*l3*m3*Y0*pow(c,-2) - 2*l3*n3*Z0*pow(c,-2) - 
		2*X0*pow(a,-2)*pow(l1,2) - 2*X0*pow(b,-2)*pow(l2,2) - 
		2*X0*pow(c,-2)*pow(l3,2);
    	coef3[7]=
		(-2*l1*m1*X0)/a/a - (2*l2*m2*X0)/b/b - 
		(2*l3*m3*X0)/c/c - (2*m1*m1*Y0)/a/a - 
		(2*m2*m2*Y0)/b/b - (2*m3*m3*Y0)/c/c - 
		(2*m1*n1*Z0)/a/a - (2*m2*n2*Z0)/b/b - 
		(2*m3*n3*Z0)/c/c;
    	coef3[8]=
		(-2*l1*n1*X0)/a/a - (2*l2*n2*X0)/b/b - 
		(2*l3*n3*X0)/c/c - (2*m1*n1*Y0)/a/a - 
		(2*m2*n2*Y0)/b/b - (2*m3*n3*Y0)/c/c - 
		(2*n1*n1*Z0)/a/a - (2*n2*n2*Z0)/b/b - 
		(2*n3*n3*Z0)/c/c;
    	coef3[9]=
		-1 + 2*l1*m1*X0*Y0*pow(a,-2) + 2*l1*n1*X0*Z0*pow(a,-2) + 
		2*m1*n1*Y0*Z0*pow(a,-2) + 2*l2*m2*X0*Y0*pow(b,-2) + 
		2*l2*n2*X0*Z0*pow(b,-2) + 2*m2*n2*Y0*Z0*pow(b,-2) + 
		2*l3*m3*X0*Y0*pow(c,-2) + 2*l3*n3*X0*Z0*pow(c,-2) + 
		2*m3*n3*Y0*Z0*pow(c,-2) + 
		pow(a,-2)*pow(l1,2)*pow(X0,2) + 
		pow(b,-2)*pow(l2,2)*pow(X0,2) + 
		pow(c,-2)*pow(l3,2)*pow(X0,2) + 
		pow(a,-2)*pow(m1,2)*pow(Y0,2) + 
		pow(b,-2)*pow(m2,2)*pow(Y0,2) + 
		pow(c,-2)*pow(m3,2)*pow(Y0,2) + 
		pow(a,-2)*pow(n1,2)*pow(Z0,2) + 
		pow(b,-2)*pow(n2,2)*pow(Z0,2) + 
		pow(c,-2)*pow(n3,2)*pow(Z0,2);
    	divd=coef3[0];
    	for (int kk=0;kk<10;kk++){  // when a particle is initialized or updated, coef[0] is set as 1.0.
	    coef3[kk]/=divd;
    	}
    } // end if

    // the fourth octant
    a=aplus; b=bminus; c=cplus;
    if(a==b&&b==c){	// this octant is a sphere
	coef4[0]=1;
	coef4[1]=1;
	coef4[2]=1;
	coef4[3]=0;
	coef4[4]=0;
	coef4[5]=0;
	coef4[6]=-2*curr_position.getx();
	coef4[7]=-2*curr_position.gety();
	coef4[8]=-2*curr_position.getz();
	coef4[9]=pow(vfabs(curr_position),2)-a*a;
    }
    else{
    	coef4[0]=l1*l1/a/a+l2*l2/b/b+l3*l3/c/c;
    	coef4[1]=m1*m1/a/a+m2*m2/b/b+m3*m3/c/c;
    	coef4[2]=n1*n1/a/a+n2*n2/b/b+n3*n3/c/c;
    	coef4[3]=(2*l1*m1)/a/a + (2*l2*m2)/b/b + (2*l3*m3)/c/c;
    	coef4[4]=(2*m1*n1)/a/a + (2*m2*n2)/b/b + (2*m3*n3)/c/c;
    	coef4[5]=(2*l1*n1)/a/a + (2*l2*n2)/b/b + (2*l3*n3)/c/c;
    	coef4[6]=
		-2*l1*m1*Y0*pow(a,-2) - 2*l1*n1*Z0*pow(a,-2) - 
		2*l2*m2*Y0*pow(b,-2) - 2*l2*n2*Z0*pow(b,-2) - 
		2*l3*m3*Y0*pow(c,-2) - 2*l3*n3*Z0*pow(c,-2) - 
		2*X0*pow(a,-2)*pow(l1,2) - 2*X0*pow(b,-2)*pow(l2,2) - 
		2*X0*pow(c,-2)*pow(l3,2);
    	coef4[7]=
		(-2*l1*m1*X0)/a/a - (2*l2*m2*X0)/b/b - 
		(2*l3*m3*X0)/c/c - (2*m1*m1*Y0)/a/a - 
		(2*m2*m2*Y0)/b/b - (2*m3*m3*Y0)/c/c - 
		(2*m1*n1*Z0)/a/a - (2*m2*n2*Z0)/b/b - 
		(2*m3*n3*Z0)/c/c;
    	coef4[8]=
		(-2*l1*n1*X0)/a/a - (2*l2*n2*X0)/b/b - 
		(2*l3*n3*X0)/c/c - (2*m1*n1*Y0)/a/a - 
		(2*m2*n2*Y0)/b/b - (2*m3*n3*Y0)/c/c - 
		(2*n1*n1*Z0)/a/a - (2*n2*n2*Z0)/b/b - 
		(2*n3*n3*Z0)/c/c;
    	coef4[9]=
		-1 + 2*l1*m1*X0*Y0*pow(a,-2) + 2*l1*n1*X0*Z0*pow(a,-2) + 
		2*m1*n1*Y0*Z0*pow(a,-2) + 2*l2*m2*X0*Y0*pow(b,-2) + 
		2*l2*n2*X0*Z0*pow(b,-2) + 2*m2*n2*Y0*Z0*pow(b,-2) + 
		2*l3*m3*X0*Y0*pow(c,-2) + 2*l3*n3*X0*Z0*pow(c,-2) + 
		2*m3*n3*Y0*Z0*pow(c,-2) + 
		pow(a,-2)*pow(l1,2)*pow(X0,2) + 
		pow(b,-2)*pow(l2,2)*pow(X0,2) + 
		pow(c,-2)*pow(l3,2)*pow(X0,2) + 
		pow(a,-2)*pow(m1,2)*pow(Y0,2) + 
		pow(b,-2)*pow(m2,2)*pow(Y0,2) + 
		pow(c,-2)*pow(m3,2)*pow(Y0,2) + 
		pow(a,-2)*pow(n1,2)*pow(Z0,2) + 
		pow(b,-2)*pow(n2,2)*pow(Z0,2) + 
		pow(c,-2)*pow(n3,2)*pow(Z0,2);
    	divd=coef4[0];
    	for (int kk=0;kk<10;kk++){  // when a particle is initialized or updated, coef[0] is set as 1.0.
	    coef4[kk]/=divd;
    	}
    } // end if

    // the fifth octant
    a=aplus; b=bplus; c=cminus;
    if(a==b&&b==c){	// this octant is a sphere
	coef5[0]=1;
	coef5[1]=1;
	coef5[2]=1;
	coef5[3]=0;
	coef5[4]=0;
	coef5[5]=0;
	coef5[6]=-2*curr_position.getx();
	coef5[7]=-2*curr_position.gety();
	coef5[8]=-2*curr_position.getz();
	coef5[9]=pow(vfabs(curr_position),2)-a*a;
    }
    else{
    	coef5[0]=l1*l1/a/a+l2*l2/b/b+l3*l3/c/c;
    	coef5[1]=m1*m1/a/a+m2*m2/b/b+m3*m3/c/c;
    	coef5[2]=n1*n1/a/a+n2*n2/b/b+n3*n3/c/c;
    	coef5[3]=(2*l1*m1)/a/a + (2*l2*m2)/b/b + (2*l3*m3)/c/c;
    	coef5[4]=(2*m1*n1)/a/a + (2*m2*n2)/b/b + (2*m3*n3)/c/c;
    	coef5[5]=(2*l1*n1)/a/a + (2*l2*n2)/b/b + (2*l3*n3)/c/c;
    	coef5[6]=
		-2*l1*m1*Y0*pow(a,-2) - 2*l1*n1*Z0*pow(a,-2) - 
		2*l2*m2*Y0*pow(b,-2) - 2*l2*n2*Z0*pow(b,-2) - 
		2*l3*m3*Y0*pow(c,-2) - 2*l3*n3*Z0*pow(c,-2) - 
		2*X0*pow(a,-2)*pow(l1,2) - 2*X0*pow(b,-2)*pow(l2,2) - 
		2*X0*pow(c,-2)*pow(l3,2);
    	coef5[7]=
		(-2*l1*m1*X0)/a/a - (2*l2*m2*X0)/b/b - 
		(2*l3*m3*X0)/c/c - (2*m1*m1*Y0)/a/a - 
		(2*m2*m2*Y0)/b/b - (2*m3*m3*Y0)/c/c - 
		(2*m1*n1*Z0)/a/a - (2*m2*n2*Z0)/b/b - 
		(2*m3*n3*Z0)/c/c;
    	coef5[8]=
		(-2*l1*n1*X0)/a/a - (2*l2*n2*X0)/b/b - 
		(2*l3*n3*X0)/c/c - (2*m1*n1*Y0)/a/a - 
		(2*m2*n2*Y0)/b/b - (2*m3*n3*Y0)/c/c - 
		(2*n1*n1*Z0)/a/a - (2*n2*n2*Z0)/b/b - 
		(2*n3*n3*Z0)/c/c;
    	coef5[9]=
		-1 + 2*l1*m1*X0*Y0*pow(a,-2) + 2*l1*n1*X0*Z0*pow(a,-2) + 
		2*m1*n1*Y0*Z0*pow(a,-2) + 2*l2*m2*X0*Y0*pow(b,-2) + 
		2*l2*n2*X0*Z0*pow(b,-2) + 2*m2*n2*Y0*Z0*pow(b,-2) + 
		2*l3*m3*X0*Y0*pow(c,-2) + 2*l3*n3*X0*Z0*pow(c,-2) + 
		2*m3*n3*Y0*Z0*pow(c,-2) + 
		pow(a,-2)*pow(l1,2)*pow(X0,2) + 
		pow(b,-2)*pow(l2,2)*pow(X0,2) + 
		pow(c,-2)*pow(l3,2)*pow(X0,2) + 
		pow(a,-2)*pow(m1,2)*pow(Y0,2) + 
		pow(b,-2)*pow(m2,2)*pow(Y0,2) + 
		pow(c,-2)*pow(m3,2)*pow(Y0,2) + 
		pow(a,-2)*pow(n1,2)*pow(Z0,2) + 
		pow(b,-2)*pow(n2,2)*pow(Z0,2) + 
		pow(c,-2)*pow(n3,2)*pow(Z0,2);
    	divd=coef5[0];
    	for (int kk=0;kk<10;kk++){  // when a particle is initialized or updated, coef[0] is set as 1.0.
	    coef5[kk]/=divd;
    	}
    } // end if

    // the sixth octant
    a=aminus; b=bplus; c=cminus;
    if(a==b&&b==c){	// this octant is a sphere
	coef6[0]=1;
	coef6[1]=1;
	coef6[2]=1;
	coef6[3]=0;
	coef6[4]=0;
	coef6[5]=0;
	coef6[6]=-2*curr_position.getx();
	coef6[7]=-2*curr_position.gety();
	coef6[8]=-2*curr_position.getz();
	coef6[9]=pow(vfabs(curr_position),2)-a*a;
    }
    else{
    	coef6[0]=l1*l1/a/a+l2*l2/b/b+l3*l3/c/c;
    	coef6[1]=m1*m1/a/a+m2*m2/b/b+m3*m3/c/c;
    	coef6[2]=n1*n1/a/a+n2*n2/b/b+n3*n3/c/c;
    	coef6[3]=(2*l1*m1)/a/a + (2*l2*m2)/b/b + (2*l3*m3)/c/c;
    	coef6[4]=(2*m1*n1)/a/a + (2*m2*n2)/b/b + (2*m3*n3)/c/c;
    	coef6[5]=(2*l1*n1)/a/a + (2*l2*n2)/b/b + (2*l3*n3)/c/c;
    	coef6[6]=
		-2*l1*m1*Y0*pow(a,-2) - 2*l1*n1*Z0*pow(a,-2) - 
		2*l2*m2*Y0*pow(b,-2) - 2*l2*n2*Z0*pow(b,-2) - 
		2*l3*m3*Y0*pow(c,-2) - 2*l3*n3*Z0*pow(c,-2) - 
		2*X0*pow(a,-2)*pow(l1,2) - 2*X0*pow(b,-2)*pow(l2,2) - 
		2*X0*pow(c,-2)*pow(l3,2);
    	coef6[7]=
		(-2*l1*m1*X0)/a/a - (2*l2*m2*X0)/b/b - 
		(2*l3*m3*X0)/c/c - (2*m1*m1*Y0)/a/a - 
		(2*m2*m2*Y0)/b/b - (2*m3*m3*Y0)/c/c - 
		(2*m1*n1*Z0)/a/a - (2*m2*n2*Z0)/b/b - 
		(2*m3*n3*Z0)/c/c;
    	coef6[8]=
		(-2*l1*n1*X0)/a/a - (2*l2*n2*X0)/b/b - 
		(2*l3*n3*X0)/c/c - (2*m1*n1*Y0)/a/a - 
		(2*m2*n2*Y0)/b/b - (2*m3*n3*Y0)/c/c - 
		(2*n1*n1*Z0)/a/a - (2*n2*n2*Z0)/b/b - 
		(2*n3*n3*Z0)/c/c;
    	coef6[9]=
		-1 + 2*l1*m1*X0*Y0*pow(a,-2) + 2*l1*n1*X0*Z0*pow(a,-2) + 
		2*m1*n1*Y0*Z0*pow(a,-2) + 2*l2*m2*X0*Y0*pow(b,-2) + 
		2*l2*n2*X0*Z0*pow(b,-2) + 2*m2*n2*Y0*Z0*pow(b,-2) + 
		2*l3*m3*X0*Y0*pow(c,-2) + 2*l3*n3*X0*Z0*pow(c,-2) + 
		2*m3*n3*Y0*Z0*pow(c,-2) + 
		pow(a,-2)*pow(l1,2)*pow(X0,2) + 
		pow(b,-2)*pow(l2,2)*pow(X0,2) + 
		pow(c,-2)*pow(l3,2)*pow(X0,2) + 
		pow(a,-2)*pow(m1,2)*pow(Y0,2) + 
		pow(b,-2)*pow(m2,2)*pow(Y0,2) + 
		pow(c,-2)*pow(m3,2)*pow(Y0,2) + 
		pow(a,-2)*pow(n1,2)*pow(Z0,2) + 
		pow(b,-2)*pow(n2,2)*pow(Z0,2) + 
		pow(c,-2)*pow(n3,2)*pow(Z0,2);
    	divd=coef6[0];
    	for (int kk=0;kk<10;kk++){  // when a particle is initialized or updated, coef[0] is set as 1.0.
	    coef6[kk]/=divd;
    	}
    } // end if

    // the seventh octant
    a=aminus; b=bminus; c=cminus;
    if(a==b&&b==c){	// this octant is a sphere
	coef7[0]=1;
	coef7[1]=1;
	coef7[2]=1;
	coef7[3]=0;
	coef7[4]=0;
	coef7[5]=0;
	coef7[6]=-2*curr_position.getx();
	coef7[7]=-2*curr_position.gety();
	coef7[8]=-2*curr_position.getz();
	coef7[9]=pow(vfabs(curr_position),2)-a*a;
    }
    else{
    	coef7[0]=l1*l1/a/a+l2*l2/b/b+l3*l3/c/c;
    	coef7[1]=m1*m1/a/a+m2*m2/b/b+m3*m3/c/c;
    	coef7[2]=n1*n1/a/a+n2*n2/b/b+n3*n3/c/c;
    	coef7[3]=(2*l1*m1)/a/a + (2*l2*m2)/b/b + (2*l3*m3)/c/c;
    	coef7[4]=(2*m1*n1)/a/a + (2*m2*n2)/b/b + (2*m3*n3)/c/c;
    	coef7[5]=(2*l1*n1)/a/a + (2*l2*n2)/b/b + (2*l3*n3)/c/c;
    	coef7[6]=
		-2*l1*m1*Y0*pow(a,-2) - 2*l1*n1*Z0*pow(a,-2) - 
		2*l2*m2*Y0*pow(b,-2) - 2*l2*n2*Z0*pow(b,-2) - 
		2*l3*m3*Y0*pow(c,-2) - 2*l3*n3*Z0*pow(c,-2) - 
		2*X0*pow(a,-2)*pow(l1,2) - 2*X0*pow(b,-2)*pow(l2,2) - 
		2*X0*pow(c,-2)*pow(l3,2);
    	coef7[7]=
		(-2*l1*m1*X0)/a/a - (2*l2*m2*X0)/b/b - 
		(2*l3*m3*X0)/c/c - (2*m1*m1*Y0)/a/a - 
		(2*m2*m2*Y0)/b/b - (2*m3*m3*Y0)/c/c - 
		(2*m1*n1*Z0)/a/a - (2*m2*n2*Z0)/b/b - 
		(2*m3*n3*Z0)/c/c;
    	coef7[8]=
		(-2*l1*n1*X0)/a/a - (2*l2*n2*X0)/b/b - 
		(2*l3*n3*X0)/c/c - (2*m1*n1*Y0)/a/a - 
		(2*m2*n2*Y0)/b/b - (2*m3*n3*Y0)/c/c - 
		(2*n1*n1*Z0)/a/a - (2*n2*n2*Z0)/b/b - 
		(2*n3*n3*Z0)/c/c;
    	coef7[9]=
		-1 + 2*l1*m1*X0*Y0*pow(a,-2) + 2*l1*n1*X0*Z0*pow(a,-2) + 
		2*m1*n1*Y0*Z0*pow(a,-2) + 2*l2*m2*X0*Y0*pow(b,-2) + 
		2*l2*n2*X0*Z0*pow(b,-2) + 2*m2*n2*Y0*Z0*pow(b,-2) + 
		2*l3*m3*X0*Y0*pow(c,-2) + 2*l3*n3*X0*Z0*pow(c,-2) + 
		2*m3*n3*Y0*Z0*pow(c,-2) + 
		pow(a,-2)*pow(l1,2)*pow(X0,2) + 
		pow(b,-2)*pow(l2,2)*pow(X0,2) + 
		pow(c,-2)*pow(l3,2)*pow(X0,2) + 
		pow(a,-2)*pow(m1,2)*pow(Y0,2) + 
		pow(b,-2)*pow(m2,2)*pow(Y0,2) + 
		pow(c,-2)*pow(m3,2)*pow(Y0,2) + 
		pow(a,-2)*pow(n1,2)*pow(Z0,2) + 
		pow(b,-2)*pow(n2,2)*pow(Z0,2) + 
		pow(c,-2)*pow(n3,2)*pow(Z0,2);
    	divd=coef7[0];
    	for (int kk=0;kk<10;kk++){  // when a particle is initialized or updated, coef[0] is set as 1.0.
	    coef7[kk]/=divd;
    	}
    } // end if

    // the eighth octant
    a=aplus; b=bminus; c=cminus;
    if(a==b&&b==c){	// this octant is a sphere
	coef8[0]=1;
	coef8[1]=1;
	coef8[2]=1;
	coef8[3]=0;
	coef8[4]=0;
	coef8[5]=0;
	coef8[6]=-2*curr_position.getx();
	coef8[7]=-2*curr_position.gety();
	coef8[8]=-2*curr_position.getz();
	coef8[9]=pow(vfabs(curr_position),2)-a*a;
    }
    else{
    	coef8[0]=l1*l1/a/a+l2*l2/b/b+l3*l3/c/c;
    	coef8[1]=m1*m1/a/a+m2*m2/b/b+m3*m3/c/c;
    	coef8[2]=n1*n1/a/a+n2*n2/b/b+n3*n3/c/c;
    	coef8[3]=(2*l1*m1)/a/a + (2*l2*m2)/b/b + (2*l3*m3)/c/c;
    	coef8[4]=(2*m1*n1)/a/a + (2*m2*n2)/b/b + (2*m3*n3)/c/c;
    	coef8[5]=(2*l1*n1)/a/a + (2*l2*n2)/b/b + (2*l3*n3)/c/c;
    	coef8[6]=
		-2*l1*m1*Y0*pow(a,-2) - 2*l1*n1*Z0*pow(a,-2) - 
		2*l2*m2*Y0*pow(b,-2) - 2*l2*n2*Z0*pow(b,-2) - 
		2*l3*m3*Y0*pow(c,-2) - 2*l3*n3*Z0*pow(c,-2) - 
		2*X0*pow(a,-2)*pow(l1,2) - 2*X0*pow(b,-2)*pow(l2,2) - 
		2*X0*pow(c,-2)*pow(l3,2);
    	coef8[7]=
		(-2*l1*m1*X0)/a/a - (2*l2*m2*X0)/b/b - 
		(2*l3*m3*X0)/c/c - (2*m1*m1*Y0)/a/a - 
		(2*m2*m2*Y0)/b/b - (2*m3*m3*Y0)/c/c - 
		(2*m1*n1*Z0)/a/a - (2*m2*n2*Z0)/b/b - 
		(2*m3*n3*Z0)/c/c;
    	coef8[8]=
		(-2*l1*n1*X0)/a/a - (2*l2*n2*X0)/b/b - 
		(2*l3*n3*X0)/c/c - (2*m1*n1*Y0)/a/a - 
		(2*m2*n2*Y0)/b/b - (2*m3*n3*Y0)/c/c - 
		(2*n1*n1*Z0)/a/a - (2*n2*n2*Z0)/b/b - 
		(2*n3*n3*Z0)/c/c;
    	coef8[9]=
		-1 + 2*l1*m1*X0*Y0*pow(a,-2) + 2*l1*n1*X0*Z0*pow(a,-2) + 
		2*m1*n1*Y0*Z0*pow(a,-2) + 2*l2*m2*X0*Y0*pow(b,-2) + 
		2*l2*n2*X0*Z0*pow(b,-2) + 2*m2*n2*Y0*Z0*pow(b,-2) + 
		2*l3*m3*X0*Y0*pow(c,-2) + 2*l3*n3*X0*Z0*pow(c,-2) + 
		2*m3*n3*Y0*Z0*pow(c,-2) + 
		pow(a,-2)*pow(l1,2)*pow(X0,2) + 
		pow(b,-2)*pow(l2,2)*pow(X0,2) + 
		pow(c,-2)*pow(l3,2)*pow(X0,2) + 
		pow(a,-2)*pow(m1,2)*pow(Y0,2) + 
		pow(b,-2)*pow(m2,2)*pow(Y0,2) + 
		pow(c,-2)*pow(m3,2)*pow(Y0,2) + 
		pow(a,-2)*pow(n1,2)*pow(Z0,2) + 
		pow(b,-2)*pow(n2,2)*pow(Z0,2) + 
		pow(c,-2)*pow(n3,2)*pow(Z0,2);
    	divd=coef8[0];
    	for (int kk=0;kk<10;kk++){  // when a particle is initialized or updated, coef[0] is set as 1.0.
	    coef8[kk]/=divd;
    	}
    } // end if

}


bool particle::intersectWithLine(vec v, vec dirc, vec rt[], int num_oct) const{	// September 6, 2013
    REAL x0=v.getx();
    REAL y0=v.gety();
    REAL z0=v.getz();
    REAL p=dirc.getx();
    REAL q=dirc.gety();
    REAL r=dirc.getz();
    REAL a, b, c, d, e, f, g, h, i, j;
    switch (num_oct){
	case 1:
   	    a=coef1[0];
    	    b=coef1[1];
    	    c=coef1[2];
    	    d=coef1[3];
    	    e=coef1[4];
    	    f=coef1[5];
    	    g=coef1[6];
    	    h=coef1[7];
    	    i=coef1[8];
     	    j=coef1[9];
	    break;
	case 2:
   	    a=coef2[0];
    	    b=coef2[1];
    	    c=coef2[2];
    	    d=coef2[3];
    	    e=coef2[4];
    	    f=coef2[5];
    	    g=coef2[6];
    	    h=coef2[7];
    	    i=coef2[8];
     	    j=coef2[9];
	    break;
	case 3:
   	    a=coef3[0];
    	    b=coef3[1];
    	    c=coef3[2];
    	    d=coef3[3];
    	    e=coef3[4];
    	    f=coef3[5];
    	    g=coef3[6];
    	    h=coef3[7];
    	    i=coef3[8];
     	    j=coef3[9];
	    break;
	case 4:
   	    a=coef4[0];
    	    b=coef4[1];
    	    c=coef4[2];
    	    d=coef4[3];
    	    e=coef4[4];
    	    f=coef4[5];
    	    g=coef4[6];
    	    h=coef4[7];
    	    i=coef4[8];
     	    j=coef4[9];
	    break;
	case 5:
   	    a=coef5[0];
    	    b=coef5[1];
    	    c=coef5[2];
    	    d=coef5[3];
    	    e=coef5[4];
    	    f=coef5[5];
    	    g=coef5[6];
    	    h=coef5[7];
    	    i=coef5[8];
     	    j=coef5[9];
	    break;
	case 6:
   	    a=coef6[0];
    	    b=coef6[1];
    	    c=coef6[2];
    	    d=coef6[3];
    	    e=coef6[4];
    	    f=coef6[5];
    	    g=coef6[6];
    	    h=coef6[7];
    	    i=coef6[8];
     	    j=coef6[9];
	    break;
	case 7:
   	    a=coef7[0];
    	    b=coef7[1];
    	    c=coef7[2];
    	    d=coef7[3];
    	    e=coef7[4];
    	    f=coef7[5];
    	    g=coef7[6];
    	    h=coef7[7];
    	    i=coef7[8];
     	    j=coef7[9];
	    break;
	case 8:
   	    a=coef8[0];
    	    b=coef8[1];
    	    c=coef8[2];
    	    d=coef8[3];
    	    e=coef8[4];
    	    f=coef8[5];
    	    g=coef8[6];
    	    h=coef8[7];
    	    i=coef8[8];
     	    j=coef8[9];
	    break;
	default:
	    std::cout << "number of octant exceeds 8 in intersectWithLine()!" << std::cout;
	    exit(-1);
	    break;
    }

    REAL A = a*p*p + b*q*q + c*r*r + d*p*q + e*q*r + f*r*p;
    REAL B = 2*a*p*x0 + 2*b*q*y0 + 2*c*r*z0
	            + d*p*y0 + d*q*x0 + e*q*z0 + e*r*y0 + f*p*z0 + f*r*x0
	            + g*p + h*q + i*r;
    REAL C = a*x0*x0 + b*y0*y0 + c*z0*z0 + d*x0*y0 + e*y0*z0 +f*z0*x0
	            + g*x0 + h*y0 + i*z0 + j;
    
    REAL delta=B*B-4*A*C;

    if (delta < 0){
	g_debuginf<<"particle.cpp: g_iteration="<<g_iteration
		  <<" delta < 0 in intersectWithLine()"<<std::endl;
	return false;
    }
    else{
	REAL t1=(-B+sqrt(delta))/(2*A);
	REAL t2=(-B-sqrt(delta))/(2*A);
	
	rt[0].setx(t1*p + x0);
	rt[0].sety(t1*q + y0);
	rt[0].setz(t1*r + z0);
	rt[1].setx(t2*p + x0);
	rt[1].sety(t2*q + y0);
	rt[1].setz(t2*r + z0);   
	return true;
    }
}


//    1. This member function is coded based on Mathematical equations
//       in local frame, x^2/a^2 + y^2/b^2 + z^2/c^2 =1, in
//       seeking appropriate osculating circle among an infinite number of
//       osculating circles passing through the contact point.
//    2. r = 2*r1*r2/(r1+r2)
//    3. It is important to eliminate float exceptions in computations, that is, 
//       when dz/dx == infinite, coordinate x & z are switched to use dx/dz == 0.
//    4. When a point is close to the equator, for example, fabs(z)==0,
//       float exception is prone to occurring, then a switch is needed
//       as above.
REAL particle::getRadius(vec v, int asign, int bsign, int csign) const{
//REAL particle::getRadius(vec v){
    REAL a, b, c;
    switch (asign){
	case 1:
	   a = aplus;
	   break;
	case -1:
	   a = aminus;
	   break;
	default:	// approximation
	   a = (aplus+aminus)*0.5;
	   break;
    }
    switch (bsign){
	case 1:
	   b = bplus;
	   break;
	case -1:
	   b = bminus;
	   break;
	default:	// approximation
	   b = (bplus+bminus)*0.5;
	   break;
    }
    switch (csign){
	case 1:
	   c = cplus;
	   break;
	case -1:
	   c = cminus;
	   break;
	default:	// approximation
	   c = (cplus+cminus)*0.5;
	   break;
    }

//    a = aplus; b = bplus; c = cplus;
    if(a==b&&b==c)	// this octant is a sphere
	return a;

    REAL per=1.0e-4; // define when a point is close to equator
    REAL ra=a;       // semi-axles of ellipsoid
    REAL rb=b;
    REAL rc=c;

    // get the local coodinates of vector v, the point on the particle's surface
    vec v1=vcos(curr_direction_a);
    vec v2=vcos(curr_direction_b);
    vec v3=vcos(curr_direction_c);
    REAL X0=curr_position.getx();
    REAL Y0=curr_position.gety();
    REAL Z0=curr_position.getz();
    REAL x1=v.getx()-X0;
    REAL y1=v.gety()-Y0;
    REAL z1=v.getz()-Z0;
    REAL l1=v1.getx();
    REAL m1=v1.gety();
    REAL n1=v1.getz();
    REAL l2=v2.getx();
    REAL m2=v2.gety();
    REAL n2=v2.getz();
    REAL l3=v3.getx();
    REAL m3=v3.gety();
    REAL n3=v3.getz();
    REAL x=l1*x1 + m1*y1 + n1*z1;
    REAL y=l2*x1 + m2*y1 + n2*z1;
    REAL z=l3*x1 + m3*y1 + n3*z1;

    REAL tmp;
    if (fabs(z)<=c*per) {     // switch x & z, use 0 instead of infinity
	tmp=ra; ra=rc; rc=tmp;
	tmp=x; x=z; z=tmp; 
	if (fabs(z)<=a*per) { // switch y & z, use 0 instead of infinity
	    tmp=ra; ra=rb; rb=tmp;
	    tmp=y; y=z; z=tmp; 
	}     
    }

    REAL p=-rc*rc/ra/ra*x/z;
    REAL q=-rc*rc/rb/rb*y/z;
    REAL r=-rc*rc/ra/ra*(1/z+rc*rc/ra/ra*x*x/pow(z,3));
    REAL t=-rc*rc/rb/rb*(1/z+rc*rc/rb/rb*y*y/pow(z,3));
    REAL s=-pow(rc,4)/ra/ra/rb/rb*x*y/pow(z,3);
    REAL n  = sqrt(1+p*p+q*q);

    REAL A,B,C;
    A=r*t-s*s;
    B=n*(2*p*q*s-(1+p*p)*t-(1+q*q)*r);
    C=n*n*n*n;


    // if delta < 0, then it is usually -1.0e-20, caused by computational precision.
    /*
    if (B*B-4*A*C<0){
	g_debuginf<<"particle.cpp: g_iteration="<<g_iteration
		  <<" delta < 0 in getRadius()"
		  <<" delta="<<B*B-4*A*C
		  <<" -C/B="<<-C/B
		  <<std::endl;
    }
    */
    return fabs(-C/B*2.0); // 2*r1*r2/(r1+r2)
}

/////////////////////////////////////
//REAL particle::getRadius(vec v, int asign, int bsign, int csign) const{
//REAL particle::getRadius(vec v){
//    REAL a, b, c;
/*    switch (asign){
	case 1:
	   a = aplus;
	   break;
	case -1:
	   a = aminus;
	   break;
	default:	// approximation
	   a = (aplus+aminus)*0.5;
	   break;
    }
    switch (bsign){
	case 1:
	   b = bplus;
	   break;
	case -1:
	   b = bminus;
	   break;
	default:	// approximation
	   b = (bplus+bminus)*0.5;
	   break;
    }
    switch (csign){
	case 1:
	   c = cplus;
	   break;
	case -1:
	   c = cminus;
	   break;
	default:	// approximation
	   c = (cplus+cminus)*0.5;
	   break;
    }
*/
/*    a = aplus; b = bplus; c = cplus;
    if(a==b&&b==c)	// this octant is a sphere
	return a;

    REAL per=1.0e-4; // define when a point is close to equator
    REAL ra=a;       // semi-axles of ellipsoid
    REAL rb=b;
    REAL rc=c;

    // get the local coodinates of vector v, the point on the particle's surface
    vec v1=vcos(curr_direction_a);
    vec v2=vcos(curr_direction_b);
    vec v3=vcos(curr_direction_c);
    REAL X0=curr_position.getx();
    REAL Y0=curr_position.gety();
    REAL Z0=curr_position.getz();
    REAL x1=v.getx()-X0;
    REAL y1=v.gety()-Y0;
    REAL z1=v.getz()-Z0;
    REAL l1=v1.getx();
    REAL m1=v1.gety();
    REAL n1=v1.getz();
    REAL l2=v2.getx();
    REAL m2=v2.gety();
    REAL n2=v2.getz();
    REAL l3=v3.getx();
    REAL m3=v3.gety();
    REAL n3=v3.getz();
    REAL x=l1*x1 + m1*y1 + n1*z1;
    REAL y=l2*x1 + m2*y1 + n2*z1;
    REAL z=l3*x1 + m3*y1 + n3*z1;

    REAL tmp;
    if (fabs(z)<=c*per) {     // switch x & z, use 0 instead of infinity
	tmp=ra; ra=rc; rc=tmp;
	tmp=x; x=z; z=tmp; 
	if (fabs(z)<=a*per) { // switch y & z, use 0 instead of infinity
	    tmp=ra; ra=rb; rb=tmp;
	    tmp=y; y=z; z=tmp; 
	}     
    }

    REAL p=-rc*rc/ra/ra*x/z;
    REAL q=-rc*rc/rb/rb*y/z;
    REAL r=-rc*rc/ra/ra*(1/z+rc*rc/ra/ra*x*x/pow(z,3));
    REAL t=-rc*rc/rb/rb*(1/z+rc*rc/rb/rb*y*y/pow(z,3));
    REAL s=-pow(rc,4)/ra/ra/rb/rb*x*y/pow(z,3);
    REAL n  = sqrt(1+p*p+q*q);

    REAL A,B,C;
    A=r*t-s*s;
    B=n*(2*p*q*s-(1+p*p)*t-(1+q*q)*r);
    C=n*n*n*n;

*/
    // if delta < 0, then it is usually -1.0e-20, caused by computational precision.
    /*
    if (B*B-4*A*C<0){
	g_debuginf<<"particle.cpp: g_iteration="<<g_iteration
		  <<" delta < 0 in getRadius()"
		  <<" delta="<<B*B-4*A*C
		  <<" -C/B="<<-C/B
		  <<std::endl;
    }
    */
//    return fabs(-C/B*2.0); // 2*r1*r2/(r1+r2)
//}
/////////////////////////////////////////////


void particle::clearForce(){
    force=const_force;
    moment=const_moment;

    force += vec(0,0,-G*mass*GRVT_SCL); // Unit is Newton, GRVT_SCL is for amplification.
    inContact = false;

    if (getType()==3){ // pile
	force -= vec(0,0,-G*mass*GRVT_SCL); 
    }

#ifdef MOMENT
	REAL m[20]={ 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
		       80, 70, 60, 50, 40, 30, 20, 10, 0};
#ifdef SLIP
	for (int i=0;i<20;++i) m[i] *= 1.0e-8;
#else
	for (int i=0;i<20;++i) m[i] *= 2.0e-8; 
#endif
	int s[20];
	for (int i=0;i<20;++i)
	    s[i] = START + i*100;
	
	for (int i=0;i<19;++i)
	    if (g_iteration>=s[i] && g_iteration<s[i+1] )
		moment += vec(0,m[i],0);
	if (g_iteration>=s[19] )
	    moment += vec(0,m[19],0);
#endif
}

dem::vec particle::calculateInitialCohesiveForce(){
    dem::vec fc_tmp;

    REAL atf=DMP_F*TIMESTEP; 
//    REAL atm=DMP_M*TIMESTEP;  
    fc_tmp = mass*MASS_SCL*( TIMESTEP*curr_acceleration+prev_velocity*(1-(2-atf)/(2+atf)) )/(TIMESTEP*2.0/(2+atf))-force;

    return fc_tmp;
} // calculateInitialCohesiveForce()


// central difference integration method
void particle::update() {	// August 16, 2013

    if (getType()==0 || getType()==5 || getType()==7) { // 0-free, 1-fixed, 5-free bounary particle, 7-impacting ellipsoidal bullet
	// It is important to distinguish global frame from local frame!
	vec prev_local_omga;
	vec curr_local_omga;
	vec local_moment;
	vec tmp;
	REAL atf=DMP_F*TIMESTEP; 
	REAL atm=DMP_M*TIMESTEP; 

//	// get global previous center of mass before the change of curr_directions, August 16, 2013
//	vec pre_center_mass, curr_center_mass;	//global coordinates of previous and current center of mass
//	pre_center_mass = pre_position + globalVec(local_center_mass);

	// moment: angular kinetics (rotational) equations are in local frame,
	// so global values need to be converted to those in local frame when applying equations
	tmp=vcos(getCurrDirecA()); local_moment.setx(tmp%moment); prev_local_omga.setx(tmp%prev_omga); // l1,m1,n1
	tmp=vcos(getCurrDirecB()); local_moment.sety(tmp%moment); prev_local_omga.sety(tmp%prev_omga); // l2,m2,n2
	tmp=vcos(getCurrDirecC()); local_moment.setz(tmp%moment); prev_local_omga.setz(tmp%prev_omga); // l3,m3,n3
	
	curr_local_omga.setx( prev_local_omga.getx()*(2-atm)/(2+atm) + local_moment.getx()/(J.getx()*MNT_SCL)*TIMESTEP*2/(2+atm) ); 
	curr_local_omga.sety( prev_local_omga.gety()*(2-atm)/(2+atm) + local_moment.gety()/(J.gety()*MNT_SCL)*TIMESTEP*2/(2+atm) );
	curr_local_omga.setz( prev_local_omga.getz()*(2-atm)/(2+atm) + local_moment.getz()/(J.getz()*MNT_SCL)*TIMESTEP*2/(2+atm) );
	
	// convert local angular velocities to those in global frame in order to rotate a particle in global space
	tmp=vcos( vec(curr_direction_a.getx(),curr_direction_b.getx(),curr_direction_c.getx()) ); // l1,l2,l3
	curr_omga.setx(tmp%curr_local_omga);
	
	tmp=vcos( vec(curr_direction_a.gety(),curr_direction_b.gety(),curr_direction_c.gety()) ); // m1,m2,m3
	curr_omga.sety(tmp%curr_local_omga);   
	
	tmp=vcos( vec(curr_direction_a.getz(),curr_direction_b.getz(),curr_direction_c.getz()) ); // n1,n2,n3
	curr_omga.setz(tmp%curr_local_omga);	

	curr_direction_a=vacos(normalize(rotateVec(vcos(prev_direction_a),curr_omga*TIMESTEP)));
	curr_direction_b=vacos(normalize(rotateVec(vcos(prev_direction_b),curr_omga*TIMESTEP)));
	curr_direction_c=vacos(normalize(rotateVec(vcos(prev_direction_c),curr_omga*TIMESTEP)));

	// force: translational kinetics equations are in global frame
	curr_velocity=prev_velocity*(2-atf)/(2+atf)+force/(mass*MASS_SCL)*TIMESTEP*2/(2+atf);
//curr_position = prev_position + curr_velocity*TIMESTEP;
	curr_center_mass = prev_center_mass + curr_velocity*TIMESTEP;	// August 16, 2013
//REAL move_dist = vfabs(curr_center_mass-prev_center_mass);
//if(move_dist>0.02*getMinRadius()){
//  curr_center_mass = prev_center_mass + 0.02*getMinRadius()*normalize(curr_center_mass-prev_center_mass);
//}

	curr_position = curr_center_mass + globalVec(-local_center_mass);

    }
#ifdef MOMENT
    else if (getType()==2) { //special case 2 (moment): translate first, then rotate
	vec prev_local_omga;
	vec curr_local_omga;
	vec local_moment;
	vec tmp;
	REAL atf=DMP_F*TIMESTEP; 
	REAL atm=DMP_M*TIMESTEP; 
	curr_velocity=prev_velocity*(2-atf)/(2+atf)+force/(mass*MASS_SCL)*TIMESTEP*2/(2+atf);
	if (g_iteration < START){	// only translate, then displacements of center_mass and center_geo are the same
	    curr_center_mass = prev_center_mass + curr_velocity*TIMESTEP;	
	    curr_position = prev_position + curr_velocity*TIMESTEP;
	}

	tmp=vcos(getCurrDirecA()); local_moment.setx(tmp%moment); prev_local_omga.setx(tmp%prev_omga); // l1,m1,n1
	tmp=vcos(getCurrDirecB()); local_moment.sety(tmp%moment); prev_local_omga.sety(tmp%prev_omga); // l2,m2,n2
	tmp=vcos(getCurrDirecC()); local_moment.setz(tmp%moment); prev_local_omga.setz(tmp%prev_omga); // l3,m3,n3
	
	curr_local_omga.setx( prev_local_omga.getx()*(2-atm)/(2+atm) + local_moment.getx()/(J.getx()*MNT_SCL)*TIMESTEP*2/(2+atm) ); 
	curr_local_omga.sety( prev_local_omga.gety()*(2-atm)/(2+atm) + local_moment.gety()/(J.gety()*MNT_SCL)*TIMESTEP*2/(2+atm) );
	curr_local_omga.setz( prev_local_omga.getz()*(2-atm)/(2+atm) + local_moment.getz()/(J.getz()*MNT_SCL)*TIMESTEP*2/(2+atm) );

	if (g_iteration >= START) {	
	    tmp=vcos( vec(curr_direction_a.getx(),curr_direction_b.getx(),curr_direction_c.getx()) ); // l1,l2,l3
	    curr_omga.setx(tmp%curr_local_omga);
	    
	    tmp=vcos( vec(curr_direction_a.gety(),curr_direction_b.gety(),curr_direction_c.gety()) ); // m1,m2,m3
	    curr_omga.sety(tmp%curr_local_omga);   
	    
	    tmp=vcos( vec(curr_direction_a.getz(),curr_direction_b.getz(),curr_direction_c.getz()) ); // n1,n2,n3
	    curr_omga.setz(tmp%curr_local_omga);
	    
	    curr_direction_a=vacos(normalize(rotateVec(vcos(prev_direction_a),curr_omga*TIMESTEP)));
	    curr_direction_b=vacos(normalize(rotateVec(vcos(prev_direction_b),curr_omga*TIMESTEP)));
	    curr_direction_c=vacos(normalize(rotateVec(vcos(prev_direction_c),curr_omga*TIMESTEP)));
	}
    }
#endif
    else if (getType()==3) { //special case 3 (displacemental ellipsoidal pile): translate in vertical direction only
	curr_velocity.setx(0);	
	curr_velocity.sety(0);
	curr_velocity.setz(-PILE_RATE);
	curr_center_mass = prev_center_mass + curr_velocity*TIMESTEP;	
        curr_position = prev_position + curr_velocity*TIMESTEP;
    }
    else if (getType()==4) { //special case 4 (impacting ellipsoidal penetrator): impact with inital velocity in vertical direction only 
	REAL atf=DMP_F*TIMESTEP; 
	curr_velocity=prev_velocity*(2-atf)/(2+atf)+force/(mass*MASS_SCL)*TIMESTEP*2/(2+atf);
	curr_velocity.setx(0);	
	curr_velocity.sety(0);
	curr_center_mass = prev_center_mass + curr_velocity*TIMESTEP;	
        curr_position = prev_position + curr_velocity*TIMESTEP;
	
    }

    // Below is needed for all cases
    // ensure three axles perpendicular to each other, and being unit vector
    if(curr_direction_a==0)
	curr_direction_a=vacos(normalize(vcos(curr_direction_b)*vcos(curr_direction_c)));
    if(curr_direction_b==0)
	curr_direction_b=vacos(normalize(vcos(curr_direction_c)*vcos(curr_direction_a)));
    if(curr_direction_c==0)
	curr_direction_c=vacos(normalize(vcos(curr_direction_a)*vcos(curr_direction_b)));


    // caclulate curr_acceleration and curr_acce_rotate
    // April 23, 2014
    curr_acceleration = (curr_velocity-prev_velocity)/TIMESTEP;
    curr_acce_rotate  = (curr_omga-prev_omga)/TIMESTEP;

    // update acceleration, April 23, 2014
    prev_acceleration = curr_acceleration;
    prev_acce_rotate = curr_acce_rotate;
    prev_position=curr_position;
    prev_center_mass=curr_center_mass;
    prev_direction_a=curr_direction_a;
    prev_direction_b=curr_direction_b;
    prev_direction_c=curr_direction_c;
    prev_velocity=curr_velocity;
    prev_omga=curr_omga;
    pre_force=force; 
    pre_moment=moment;

    cntnum=0;
//    GlobCoef(1, 1, 1);   // every time the particle is updated, the algebra expression is also updated. 
    GlobCoef();	// September 6, 2013
// update coefficients when it is needed, August 16, 2013
}


// calculate the average stress for the particle
// April 25, 2014
void particle::calcStress(){

    // calculate the average stress caused by body force, i.e. gravity here
    vec grav_global = vec(0,0,-G*mass*GRVT_SCL);
    vec grav_local = localVec(grav_global);
    matrix grav_mat(3,1);
    // since integrals are calculated in the local coordinates, then we need to 
    // use the vector in local coordinates
    grav_mat(1,1) = grav_local.getx();
    grav_mat(2,1) = grav_local.gety();
    grav_mat(3,1) = grav_local.getz();

    matrix int_mat(1,3);
    int_mat(1,1) = int_x; int_mat(1,2) = int_y; int_mat(1,3) = int_z;

    matrix sigma_g(3,3);
    sigma_g = density*grav_mat*int_mat;

    // calculate the average stress caused by acceleration
    matrix sigma_a(3,3);
    
    vec local_acceleration = localVec(curr_acceleration); 
    REAL a1 = local_acceleration.getx();	// acceleration at the 
    REAL a2 = local_acceleration.gety();
    REAL a3 = local_acceleration.getz();

    vec local_acce_rotate = localVec(curr_acce_rotate);
    REAL alpha1 = local_acce_rotate.getx();
    REAL alpha2 = local_acce_rotate.gety();
    REAL alpha3 = local_acce_rotate.getz();

    vec local_omga = localVec(curr_omga);
    REAL omga1 = local_omga.getx();
    REAL omga2 = local_omga.gety();
    REAL omga3 = local_omga.getz();

    sigma_a(1,1) = density*( a1*int_x + alpha2*int_xz - alpha3*int_xy 
			   + omga2*omga1*int_xy - omga2*omga2*int_x2
			   + omga3*omga1*int_xz - omga3*omga3*int_x2 );

    sigma_a(1,2) = density*( a1*int_y + alpha2*int_yz - alpha3*int_y2
			   + omga2*omga1*int_y2 - omga2*omga2*int_x2 
			   + omga3*omga1*int_yz - omga3*omga3*int_xy );

    sigma_a(1,3) = density*( a1*int_z + alpha2*int_xz - alpha3*int_yz
			   + omga2*omga1*int_yz - omga2*omga2*int_xz 
			   + omga3*omga1*int_z2 - omga3*omga3*int_xz );

    sigma_a(2,1) = density*( a2*int_x - alpha1*int_xz + alpha3*int_x2
			   - omga1*omga1*int_xy + omga1*omga2*int_x2 
			   + omga3*omga2*int_xz - omga3*omga3*int_xy );

    sigma_a(2,2) = density*( a2*int_y - alpha1*int_yz + alpha3*int_xy
			   - omga1*omga1*int_x2 + omga1*omga2*int_xy 
			   + omga3*omga2*int_yz - omga3*omga3*int_y2 );

    sigma_a(2,3) = density*( a2*int_z - alpha1*int_z2 + alpha3*int_xz
			   - omga1*omga1*int_yz + omga1*omga2*int_xz 
			   + omga3*omga2*int_z2 - omga3*omga3*int_yz );

    sigma_a(3,1) = density*( a3*int_x + alpha1*int_xy - alpha2*int_x2
			   - omga1*omga1*int_xz + omga1*omga3*int_x2 
			   - omga2*omga2*int_xz + omga2*omga3*int_xy );

    sigma_a(3,2) = density*( a3*int_y + alpha1*int_y2 - alpha2*int_xy
			   - omga1*omga1*int_yz + omga1*omga3*int_xy 
			   - omga2*omga2*int_yz + omga2*omga3*int_y2 );

    sigma_a(3,3) = density*( a3*int_z + alpha1*int_yz - alpha2*int_xz
			   - omga1*omga1*int_z2 + omga1*omga3*int_xz 
			   - omga2*omga2*int_z2 + omga2*omga3*int_yz );

    average_stress = (average_stress - sigma_a + sigma_g)/volume;

} // calcStress()


int particle::calculateBreakPlane(){

    if(isAbleDivide==false) return -1; // cannot be sub-divided

    int break_plane = -1;
    average_stress = 0.5*(average_stress+average_stress.getTrans());	// in case it is slightly not symmetric, it is local coordinate
    // calculate principal tensile stress and principal directions, 
    // refer to http://barnesc.blogspot.com/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
    REAL sigma[3][3];
    sigma[0][0]=average_stress(1,1); sigma[0][1]=average_stress(1,2); sigma[0][2]=average_stress(1,3);
    sigma[1][0]=average_stress(2,1); sigma[1][1]=average_stress(2,2); sigma[1][2]=average_stress(2,3);
    sigma[2][0]=average_stress(3,1); sigma[2][1]=average_stress(3,2); sigma[2][2]=average_stress(3,3);
    // allocate memory for the eigenvalues/vectors
    REAL V[3][3];	// eigenvectors
    REAL d[3];		// eigenvalues

    dem::eigen_decomposition(sigma,V,d);	// d[0]<d[1]<d[2], d[2] is the maximum tensile stress

    // compare the maximum tensile stress d[2] with the critical stress
    // the critical stress should be a function related to grain size, it is not established yet,
//    if(d[2]<dem::sigma_critical) return -1;	// not break
// our stress formulation will give sigma_max = 0 for particles under unixial compression experiment, we noticed that this not correct, since from Uintah simulations, we know that there is actually tensile volume inside the particle. But there is no other way to calculate stress inside the particle. Then we change the fracture criterion from sigma_max > sigma_critical to 2tau-p
// where tau = 0.5*(sigma_max-sigma_min) & p = 0.5*(sigma_max+sigma_min)

//    //-------------- below is the maximum shear stress criterion --------------
//    REAL tau = 0.5*(d[2]-d[0]); REAL p = 0.5*(d[2]+d[0]);
//    if(2*tau-p<dem::sigma_critical) return -1;	// not break
//    // ------------- above is the maximum shear stress criterion --------------

    //-------------- below is the Hoek-Brown criterion ------------------------
    REAL sigma1 = -d[0]; REAL sigma3 = -d[2];	// in their notation, compression is positive, and sigma1 is the major stress
  if( sigma1 > sigma3+strengthHoek*sqrt(mi*sigma3/strengthHoek+1) ){	// break
    //-------------- agove is the Hoek-Brown criterion ------------------------

    // if this is a spherical particle, then set its principal directions to be the same as the principal directions of stress
    if(aplus==aminus && bplus==bminus && cplus==cminus
    && aplus==bplus  && bplus==cplus ){	// means this particle is sphere, do not write as a==b==c
    	rotatePrincipalDirections(V);
	rotateBreakPlaneRandomly();	// rotate the break plane randomly a little
	// now the principal directions are the same as the that of stress, so the local stress are principal stresses
	average_stress(1,1) = d[0]; average_stress(1,2) = 0;    average_stress(1,3) = 0; 
	average_stress(2,1) = 0;    average_stress(2,2) = d[1]; average_stress(2,3) = 0; 
	average_stress(3,1) = 0;    average_stress(3,2) = 0;    average_stress(3,3) = d[2]; 
    }

    // at present the break plane should be along the principal directions of the poly-ellipsoid
    REAL a_min = getAMin(); REAL a_max = getAMax();
    REAL b_min = getBMin(); REAL b_max = getBMax();
    REAL c_min = getCMin(); REAL c_max = getCMax();
    break_plane = -1;
    if(a_min>0.2*a_max && b_min>0.2*b_max && c_min>0.2*c_max){  // case (1)
	// case (1): calculate the principle stress along three principle axles and then
	// determine which plane is the break plane
	break_plane = getBreakPlane();	// 1 means ab-plane, 2 means ac-plane, 3 means bc-plane
        numBrokenType1++;
	typeBroken = -1;
	if(break_plane!=1 && break_plane!=2 && break_plane!=3){
	    std::cout << "breaking plane should be one of xy, xz or yz!" << std::endl;
	    exit(-1);
	}
    }
    else if(a_min<0.2*a_max && b_min<0.2*b_max && c_min<0.2*c_max){  // case (2)
	isAbleDivide = false;	// do not break
    }
    else if(a_min<0.2*a_max && b_min<0.2*b_max && c_min>0.2*c_max){  // case (3)(4)_1
	if(isCLongEnough()){ // c^+ + c^- > 1.2max(a_max,b_max)
	    break_plane = 1;	// break along plane-ab
	    numBrokenType1++;
	    typeBroken = -1;
	}
	else{
	    isAbleDivide = false;
	}
    }
    else if(a_min<0.2*a_max && b_min>0.2*b_max && c_min<0.2*c_max){  // case (3)(4)_2
	if(isBLongEnough()){ // b^+ + b^- > 1.2max(a_max,c_max)
	    break_plane = 2;	// break along plane-ac
	    numBrokenType1++;
	    typeBroken = -1;
	}
	else{
	    isAbleDivide = false;;
	}
    }
    else if(a_min>0.2*a_max && b_min<0.2*b_max && c_min<0.2*c_max){  // case (3)(4)_3
	if(isALongEnough()){ // a^+ + a^- > 1.2max(b_max,c_max)
	    break_plane = 3;	// break along plane-bc
	    numBrokenType1++;
	    typeBroken = -1;
	}
     	else{
	    isAbleDivide = false;;
    	}
    }
    else if(a_min<0.2*a_max && b_min>0.2*b_max && c_min>0.2*c_max){  // case (5)_1
	break_plane = getBreakPlane(1);
	numBrokenType1++;
	typeBroken = -1;
    }
    else if(a_min>0.2*a_max && b_min<0.2*b_max && c_min>0.2*c_max){  // case (5)_2
	break_plane = getBreakPlane(2);
	numBrokenType1++;
	typeBroken = -1;
    }
    else if(a_min>0.2*a_max && b_min>0.2*b_max && c_min<0.2*c_max){  // case (5)_3
	break_plane = getBreakPlane(3);
	numBrokenType1++;
	typeBroken = -1;
    }
    else {
	std::cout << "Error: it is more than 5 cases in fracture model..." << std::endl;
	exit(-1);
    }
  } // the Hoek-Brown sub-division criterion
  else
 if(numCriticalContacts>=3){ // break
    dem::vec unitNormal = dem::normalize(dem::vec(contact1-contact2)*dem::vec(contact1-contact3)); // contact is global coordinate
    // if this is a spherical particle, then set its principal directions to be the same as the principal directions of stress
    if(aplus==aminus && bplus==bminus && cplus==cminus
    && aplus==bplus  && bplus==cplus ){	// means this particle is sphere, do not write as a==b==c
    	rotatePrincipalDirections(unitNormal); 
	rotateBreakPlaneRandomly();	// rotate the break plane randomly a little
	numBrokenType2++;
	typeBroken = 1;
	return 3;	// if rotate, means (1) this is sphere, it is normal shape, so its break plane is determined by contacts
			// (2) after rotation, curr_direction_a is orthogonal to the break plane.
			// 1 means ab-plane, 2 means ac-plane, 3 means bc-plane
    }

    // at present the break plane should be along the principal directions of the poly-ellipsoid
    REAL a_min = getAMin(); REAL a_max = getAMax();
    REAL b_min = getBMin(); REAL b_max = getBMax();
    REAL c_min = getCMin(); REAL c_max = getCMax();
    break_plane = -1;
	
    dem::vec direcA=vcos(curr_direction_a);
    dem::vec direcB=vcos(curr_direction_b);
    dem::vec direcC=vcos(curr_direction_c);
    REAL dotAN=fabs(direcA%unitNormal);	// absolute fabs() is very important 
    REAL dotBN=fabs(direcB%unitNormal);
    REAL dotCN=fabs(direcC%unitNormal);
    if(a_min>0.2*a_max && b_min>0.2*b_max && c_min>0.2*c_max){  // case (1)
	// case (1): determine the break plane based on the contacts only
	if(dotAN>=dotBN && dotAN>=dotCN)
	    break_plane = 3;	// 1 means ab-plane, 2 means ac-plane, 3 means bc-plane
	else if(dotBN>=dotAN && dotBN>=dotCN)
	    break_plane = 2;	// 1 means ab-plane, 2 means ac-plane, 3 means bc-plane
	else if(dotCN>=dotAN && dotCN>=dotBN)
	    break_plane = 1;	// 1 means ab-plane, 2 means ac-plane, 3 means bc-plane

	if(break_plane!=1 && break_plane!=2 && break_plane!=3){
	    std::cout << "breaking plane should be one of xy, xz or yz!" << std::endl;
	    exit(-1);
	}

	numBrokenType2++;
	typeBroken = 1;
    }
    else if(a_min<0.2*a_max && b_min<0.2*b_max && c_min<0.2*c_max){  // case (2)
	isAbleDivide = false;	// do not break
    }
    else if(a_min<0.2*a_max && b_min<0.2*b_max && c_min>0.2*c_max){  // case (3)(4)_1
	if(isCLongEnough()){ // c^+ + c^- > 1.2max(a_max,b_max)
	    break_plane = 1;	// break along plane-ab
	    numBrokenType2++;
	    typeBroken = 1;
	}
	else{
	    isAbleDivide = false;
	}
    }
    else if(a_min<0.2*a_max && b_min>0.2*b_max && c_min<0.2*c_max){  // case (3)(4)_2
	if(isBLongEnough()){ // b^+ + b^- > 1.2max(a_max,c_max)
	    break_plane = 2;	// break along plane-ac
	    numBrokenType2++;
	    typeBroken = 1;
	}
	else{
	    isAbleDivide = false;;
	}
    }
    else if(a_min>0.2*a_max && b_min<0.2*b_max && c_min<0.2*c_max){  // case (3)(4)_3
	if(isALongEnough()){ // a^+ + a^- > 1.2max(b_max,c_max)
	    break_plane = 3;	// break along plane-bc
	    numBrokenType2++;
	    typeBroken = 1;
	}
     	else{
	    isAbleDivide = false;;
    	}
    }
    else if(a_min<0.2*a_max && b_min>0.2*b_max && c_min>0.2*c_max){  // case (5)_1
	if(isBLongerThanC()){  // bplus+bminus>1.2*(cplus+cminus)
	    break_plane = 2;	// break along plane-ac
	    numBrokenType2++;
	    typeBroken = 1;
	}
	else if(isCLongerThanB()){  // cplus+cminus>1.2*(bplus+bminus)
	    break_plane = 1;	// break along plane-ab
	    numBrokenType2++;
	    typeBroken = 1;
	}
	else{
	    numBrokenType2++;
	    typeBroken = 1;
	    if(dotBN>dotCN)	// axle c should be in
		return 2;	// break along plane-ac
	    else
		return 1;	// break along plane-ab
	}
    }
    else if(a_min>0.2*a_max && b_min<0.2*b_max && c_min>0.2*c_max){  // case (5)_2
	if(isALongerThanC()){  // Aplus+Aminus>1.2*(cplus+cminus)
	    typeBroken = 1;
	    numBrokenType2++;
	    return 3;
	}
	else if(isCLongerThanA()){  // cplus+cminus>1.2*(aplus+aminus)
	    typeBroken = 1;
	    numBrokenType2++;   
	    return 1;
	}
	else{
	    typeBroken = 1;
	    numBrokenType2++;
	    if(dotAN>dotCN)	// axle c should be in
		return 3;
	    else
		return 1;
 	}
    }
    else if(a_min>0.2*a_max && b_min>0.2*b_max && c_min<0.2*c_max){  // case (5)_3
	if(isALongerThanB()){  // Aplus+Aminus>1.2*(bplus+bminus)
	    typeBroken = 1;
	    numBrokenType2++;
	    return 3;	// break along plane-bc
	}
	else if(isBLongerThanA()){  // bplus+bminus>1.2*(aplus+aminus)
	    typeBroken = 1;
	    numBrokenType2++;
	    return 2;	// break along plane-ac
	}
	else{
	    typeBroken = 1;
	    numBrokenType2++;
	    if(dotAN>dotBN)	// axle b should be in
		return 3;
	    else
		return 2;
	}
    }
    else {
	std::cout << "Error: it is more than 5 cases in fracture model..." << std::endl;
	std::cout << "a_min: " << a_min << ", a_max: " << a_max << ", b_min: " << b_min <<
		   ", b_max: " << b_max << ", c_min: " << c_min << ", c_max: " << c_max << std::endl;
	exit(-1);
    }

  } // maximum tensile stress at 3 contacts
  else if(numCriticalContacts==2){ // break, with only two critical contacts
    if(stress2<1.5*strengthContact)
	return -1;

    dem::vec unitNormal = dem::normalize(dem::vec(contact1-contact2)); // contact is global coordinate
    // if this is a spherical particle, then set its principal directions to be the same as the principal directions of stress
    if(aplus==aminus && bplus==bminus && cplus==cminus
    && aplus==bplus  && bplus==cplus ){	// means this particle is sphere, do not write as a==b==c
    	rotatePrincipalDirections(unitNormal); 
	rotateBreakPlaneRandomly();	// rotate the break plane randomly a little
	numBrokenType2++;
	typeBroken = 2;
	return 1;	// if rotate, means (1) this is sphere, it is normal shape, so its break plane is determined by contacts
			// (2) after rotation, curr_direction_a is the vector connecting contact1 and contact2
			// 1 means ab-plane, 2 means ac-plane, 3 means bc-plane
    }

    // at present the break plane should be along the principal directions of the poly-ellipsoid
    REAL a_min = getAMin(); REAL a_max = getAMax();
    REAL b_min = getBMin(); REAL b_max = getBMax();
    REAL c_min = getCMin(); REAL c_max = getCMax();
    break_plane = -1;
	
    dem::vec direcA=vcos(curr_direction_a);	// direcA is unit vector
    dem::vec direcB=vcos(curr_direction_b);
    dem::vec direcC=vcos(curr_direction_c);
    REAL dotAN=fabs(direcA%unitNormal);	// absolute fabs() is very important 
    REAL dotBN=fabs(direcB%unitNormal);
    REAL dotCN=fabs(direcC%unitNormal);
    if(a_min>0.2*a_max && b_min>0.2*b_max && c_min>0.2*c_max){  // case (1)
	// case (1): determine the break plane based on the contacts only
	if(dotAN>=dotBN && dotAN>=dotCN){
	    if(b_min+b_max<c_min+c_max)
		break_plane = 1;	// 1 means ab-plane, 2 means ac-plane, 3 means bc-plane
	    else
		break_plane = 2;	// 1 means ab-plane, 2 means ac-plane, 3 means bc-plane
	}	
	else if(dotBN>=dotAN && dotBN>=dotCN){
	    if(a_min+a_max<c_min+c_max)
		break_plane = 1;	// 1 means ab-plane, 2 means ac-plane, 3 means bc-plane
	    else
		break_plane = 3;	// 1 means ab-plane, 2 means ac-plane, 3 means bc-plane
	}
	else if(dotCN>=dotAN && dotCN>=dotBN){
	    if(a_min+a_max<b_min+b_max)
		break_plane = 2;	// 1 means ab-plane, 2 means ac-plane, 3 means bc-plane
	    else
		break_plane = 3;	// 1 means ab-plane, 2 means ac-plane, 3 means bc-plane
	}

	if(break_plane!=1 && break_plane!=2 && break_plane!=3){
	    std::cout << "breaking plane should be one of xy, xz or yz!" << std::endl;
	    exit(-1);
	}

	numBrokenType2++;
	typeBroken = 2;
    }
    else if(a_min<0.2*a_max && b_min<0.2*b_max && c_min<0.2*c_max){  // case (2)
	isAbleDivide = false;	// do not break
    }
    else if(a_min<0.2*a_max && b_min<0.2*b_max && c_min>0.2*c_max){  // case (3)(4)_1
	if(isCLongEnough()){ // c^+ + c^- > 1.2max(a_max,b_max)
	    break_plane = 1;	// break along plane-ab
	    numBrokenType2++;
	    typeBroken = 2;
	}
	else{
	    isAbleDivide = false;
	}
    }
    else if(a_min<0.2*a_max && b_min>0.2*b_max && c_min<0.2*c_max){  // case (3)(4)_2
	if(isBLongEnough()){ // b^+ + b^- > 1.2max(a_max,c_max)
	    break_plane = 2;	// break along plane-ac
	    numBrokenType2++;
	    typeBroken = 2;
	}
	else{
	    isAbleDivide = false;;
	}
    }
    else if(a_min>0.2*a_max && b_min<0.2*b_max && c_min<0.2*c_max){  // case (3)(4)_3
	if(isALongEnough()){ // a^+ + a^- > 1.2max(b_max,c_max)
	    break_plane = 3;	// break along plane-bc
	    numBrokenType2++;
	    typeBroken = 2;
	}
     	else{
	    isAbleDivide = false;;
    	}
    }
    else if(a_min<0.2*a_max && b_min>0.2*b_max && c_min>0.2*c_max){  // case (5)_1
	if(isBLongerThanC()){  // bplus+bminus>1.2*(cplus+cminus)
	    break_plane = 2;	// break along plane-ac
	    numBrokenType2++;
	    typeBroken = 2;
	}
	else if(isCLongerThanB()){  // cplus+cminus>1.2*(bplus+bminus)
	    break_plane = 1;	// break along plane-ab
	    numBrokenType2++;
	    typeBroken = 2;
	}
	else{
	    numBrokenType2++;
	    typeBroken = 2;
	    if(dotCN>dotBN)	// axle c should be in
		return 2;	// break along plane-ac
	    else
		return 1;	// break along plane-ab
	}
    }
    else if(a_min>0.2*a_max && b_min<0.2*b_max && c_min>0.2*c_max){  // case (5)_2
	if(isALongerThanC()){  // Aplus+Aminus>1.2*(cplus+cminus)
	    typeBroken = 2;
	    numBrokenType2++;
	    return 3;
	}
	else if(isCLongerThanA()){  // cplus+cminus>1.2*(aplus+aminus)
	    typeBroken = 2;
	    numBrokenType2++;   
	    return 1;
	}
	else{
	    typeBroken = 2;
	    numBrokenType2++;
	    if(dotCN>dotAN)	// axle c should be in
		return 3;
	    else
		return 1;
 	}
    }
    else if(a_min>0.2*a_max && b_min>0.2*b_max && c_min<0.2*c_max){  // case (5)_3
	if(isALongerThanB()){  // Aplus+Aminus>1.2*(bplus+bminus)
	    typeBroken = 2;
	    numBrokenType2++;
	    return 3;	// break along plane-bc
	}
	else if(isBLongerThanA()){  // bplus+bminus>1.2*(aplus+aminus)
	    typeBroken = 2;
	    numBrokenType2++;
	    return 2;	// break along plane-ac
	}
	else{
	    typeBroken = 2;
	    numBrokenType2++;
	    if(dotBN>dotAN)	// axle b should be in
		return 3;
	    else
		return 2;
	}
    }
    else {
	std::cout << "Error: it is more than 5 cases in fracture model..." << std::endl;
	std::cout << "a_min: " << a_min << ", a_max: " << a_max << ", b_min: " << b_min <<
		   ", b_max: " << b_max << ", c_min: " << c_min << ", c_max: " << c_max << std::endl;
	exit(-1);
    }

  } // maximum tensile stress at 2 contacts

    if(isAbleDivide==false)  // cannot be divided
	return -1;

    return break_plane;

} // calculateBreakPlane


// rotate the principal directions of sphere to be the same as the principal directions of stress
void particle::rotatePrincipalDirections(REAL V[3][3]){
    // the principal directions of stress calculated as V[3][3] are currently in local coordinate
    dem::vec v1 = globalVec( vec(V[0][0],V[1][0],V[2][0]) );	// principal direction of d[0]
    dem::vec v2 = globalVec( vec(V[0][1],V[1][1],V[2][1]) );	// principal direction of d[1]
    dem::vec v3 = globalVec( vec(V[0][2],V[1][2],V[2][2]) );	// principal direction of d[2]

    // set current direction as the principal directions of stress
    // onething needs to notice is that the curr_direction_a is the angle of a between x,y,z axles
    dem::vec e1=vec(1,0,0); dem::vec e2=vec(0,1,0); dem::vec e3=vec(0,0,1);
    dem::vec v1_tmp = normalize(v1);
    curr_direction_a.setx( acos(v1_tmp%e1) );
    curr_direction_a.sety( acos(v1_tmp%e2) );
    curr_direction_a.setz( acos(v1_tmp%e3) );

    dem::vec v2_tmp = normalize(v2);
    curr_direction_b.setx( acos(v2_tmp%e1) );
    curr_direction_b.sety( acos(v2_tmp%e2) );
    curr_direction_b.setz( acos(v2_tmp%e3) );

    dem::vec v3_tmp = normalize(v3);
    curr_direction_c.setx( acos(v3_tmp%e1 ) );
    curr_direction_c.sety( acos(v3_tmp%e2 ) );
    curr_direction_c.setz( acos(v3_tmp%e3 ) );
    // since prev_directions are still the original geometry, so we set them to be the same as current ones
    // to avoid the virtual rotation acceleration calculated
    prev_direction_a = curr_direction_a; 
    prev_direction_b = curr_direction_b; 
    prev_direction_c = curr_direction_c;
} // rotatePrincipalDirections


// rotate the principal directions of sphere to be one of the normal vector of the break plane
// determined by the maximum tensile stress in contacts
// after this rotation, the curr_direction_a is the same as unitN, which is orthogonal to the 
// break plane.
void particle::rotatePrincipalDirections(dem::vec unitN){
    // generate the other two directions orthogonal to the unitN randomly
    REAL l1=unitN.getx();
    REAL m1=unitN.gety();
    REAL n1=unitN.getz();
    REAL l2,m2,n2,A,B,C;
    REAL tmp_EPS = 1e-3;
    if( (fabs(fabs(l1)-1)<tmp_EPS &&fabs(m1)<tmp_EPS && fabs(n1)<tmp_EPS)
     || (fabs(l1)<tmp_EPS &&fabs(fabs(m1)-1)<tmp_EPS && fabs(n1)<tmp_EPS)
     || (fabs(l1)<tmp_EPS &&fabs(m1)<tmp_EPS && fabs(fabs(n1)-1)<tmp_EPS) ){	// the ideal break plane is very close to the 
	return;								// principle directions, no need to rotate
    }
    else{
    	while (1) {
      	    l2 = ran(&idum)*2 - 1;
      	    // solve n2
      	    A = m1*m1 + n1*n1;
      	    B = 2 * l1 * l2 * n1;
      	    C = l1*l1*l2*l2 + m1*m1*l2*l2 - m1*m1;
      	    if (B*B - 4*A*C > EPS)
		break;
    	}

    	int sign = 2*ran(&idum)-1 > 0 ? 1:-1;
    	n2 = (-B + sign*sqrt(B*B - 4*A*C) ) / (2*A);
    	m2 = - (l1*l2 + n1*n2) / m1;
    }

    curr_direction_a=vec(acos(l1), acos(m1), acos(n1));
    curr_direction_b=vec(acos(l2), acos(m2), acos(n2));
    curr_direction_c=vacos(normalize(vcos(curr_direction_a)*vcos(curr_direction_b)));

    // since prev_directions are still the original geometry, so we set them to be the same as current ones
    // to avoid the virtual rotation acceleration calculated
    prev_direction_a = curr_direction_a; 
    prev_direction_b = curr_direction_b; 
    prev_direction_c = curr_direction_c;
} // rotatePrincipalDirections


void particle::rotateToXYZDirections(){
    curr_direction_a=vec(acos(1), acos(0), acos(0));
    curr_direction_b=vec(acos(0), acos(1), acos(0));
    curr_direction_c=vec(acos(0), acos(0), acos(1));

    // since prev_directions are still the original geometry, so we set them to be the same as current ones
    // to avoid the virtual rotation acceleration calculated
    prev_direction_a = curr_direction_a; 
    prev_direction_b = curr_direction_b; 
    prev_direction_c = curr_direction_c;
} // rotatePrincipalDirections

// rotate the break plane of spheres randomly a little,
// rotation of break plane actually is converted to the rotation of principel directions,
// since currently we can only subdivide the poly-ellipsoid along the principel directions
void particle::rotateBreakPlaneRandomly(){

    REAL theta_max = 0.1745;	// unit rad, the maximum angle that we can rotate the break plane in x,y,z directions
    REAL num = ran(&idum);	// num belongs (0, 1)
    num = num*0.5;		// num belongs (0, 0.5)
    num = num+0.5;		// num belongs (0.5, 1)
    int sign = 2*ran(&idum)-1 > 0 ? 1:-1;
    num = num*sign;
    REAL alpha = num*theta_max;	// rotation angle along x direction

    num = ran(&idum)*0.5+0.5;	// num belongs (0.5,1)
    sign = 2*ran(&idum)-1 > 0 ? 1:-1;
    REAL beta  = num*sign*theta_max;	// rotation angle along y direction

    num = ran(&idum)*0.5+0.5;	// num belongs (0.5,1)
    sign = 2*ran(&idum)-1 > 0 ? 1:-1;
    REAL gamma = num*sign*theta_max;	// rotation angle along z direction

//    REAL alpha = 0.5*dem::PI;
//    REAL beta  = 0.5*dem::PI;
//    REAL gamma = 0.5*dem::PI;

    dem::vec a_prime = dem::vec(cos(beta)*cos(gamma), 
			        cos(gamma)*sin(alpha)*sin(beta)-cos(alpha)*sin(gamma),
			        sin(alpha)*sin(gamma)+cos(alpha)*cos(gamma)*sin(beta)  );
    dem::vec b_prime = dem::vec(cos(beta)*sin(gamma),
			        cos(alpha)*cos(gamma)+sin(alpha)*sin(beta)*sin(gamma),
			        cos(alpha)*sin(beta)*sin(gamma)-cos(gamma)*sin(alpha)  );
    dem::vec c_prime = dem::vec(-sin(beta), cos(beta)*sin(alpha), cos(alpha)*cos(beta) );

    // convert the new rotated prinple directions in local coordinate to global coordinates
    dem::vec a_global = globalVec(a_prime);
    dem::vec b_global = globalVec(b_prime);
    dem::vec c_global = globalVec(c_prime);

    curr_direction_a=vec(acos(a_global.getx()), acos(a_global.gety()), acos(a_global.getz()) );
    curr_direction_b=vec(acos(b_global.getx()), acos(b_global.gety()), acos(b_global.getz()) );
    curr_direction_c=vec(acos(c_global.getx()), acos(c_global.gety()), acos(c_global.getz()) );

    // since prev_directions are still the original geometry, so we set them to be the same as current ones
    // to avoid the virtual rotation acceleration calculated
    prev_direction_a = curr_direction_a; 
    prev_direction_b = curr_direction_b; 
    prev_direction_c = curr_direction_c;    

} // rotateBreakPlaneRandomly


vec particle::localVec(vec v) const{
    // v is a vector in global coordinates, it is to be transformed into local coordinates
    vec l=vcos(vec(curr_direction_a.getx(),curr_direction_b.getx(),curr_direction_c.getx()));
    vec m=vcos(vec(curr_direction_a.gety(),curr_direction_b.gety(),curr_direction_c.gety()));
    vec n=vcos(vec(curr_direction_a.getz(),curr_direction_b.getz(),curr_direction_c.getz()));
    return l*v.getx()+m*v.gety()+n*v.getz();
}


vec particle::globalVec(vec v) const{
    vec l=vcos(curr_direction_a);
    vec m=vcos(curr_direction_b);
    vec n=vcos(curr_direction_c);
    return l*v.getx()+m*v.gety()+n*v.getz();
}


bool particle::nearestPTOnPlane(REAL p, REAL q, REAL r, REAL s, vec& ptnp, int num_oct) const {	// August 19, 2013. September 6, 2013
    if(aplus==aminus && aminus==bplus && bplus==bminus && bminus==cplus && cplus==cminus){
      vec tnm=vec(p,q,r)/sqrt(p*p+q*q+r*r);
      // signed distance from particle center to plane
      REAL l_nm=(curr_position.getx()*p+curr_position.gety()*q+curr_position.getz()*r+s)/sqrt(p*p+q*q+r*r); 
      ptnp=curr_position-l_nm*tnm;
      if( (aplus-fabs(l_nm)) / (2.0*aplus) > MINOVERLAP) // intersect
	return true;
      else              // no intersect,
	return false;
    }
    REAL a, b, c, d, e, f, g, h, i, j;
    switch (num_oct){
	case 1: 	
   	    a=coef1[0];
    	    b=coef1[1];
    	    c=coef1[2];
    	    d=coef1[3];
    	    e=coef1[4];
    	    f=coef1[5];
    	    g=coef1[6];
    	    h=coef1[7];
    	    i=coef1[8];
     	    j=coef1[9];
	    break;
	case 2: 	
   	    a=coef2[0];
    	    b=coef2[1];
    	    c=coef2[2];
    	    d=coef2[3];
    	    e=coef2[4];
    	    f=coef2[5];
    	    g=coef2[6];
    	    h=coef2[7];
    	    i=coef2[8];
     	    j=coef2[9];
	    break;
	case 3: 
   	    a=coef3[0];
    	    b=coef3[1];
    	    c=coef3[2];
    	    d=coef3[3];
    	    e=coef3[4];
    	    f=coef3[5];
    	    g=coef3[6];
    	    h=coef3[7];
    	    i=coef3[8];
     	    j=coef3[9];
	    break;
	case 4: 	
   	    a=coef4[0];
    	    b=coef4[1];
    	    c=coef4[2];
    	    d=coef4[3];
    	    e=coef4[4];
    	    f=coef4[5];
    	    g=coef4[6];
    	    h=coef4[7];
    	    i=coef4[8];
     	    j=coef4[9];
	    break;
	case 5: 
   	    a=coef5[0];
    	    b=coef5[1];
    	    c=coef5[2];
    	    d=coef5[3];
    	    e=coef5[4];
    	    f=coef5[5];
    	    g=coef5[6];
    	    h=coef5[7];
    	    i=coef5[8];
     	    j=coef5[9];
	    break;
	case 6: 	
   	    a=coef6[0];
    	    b=coef6[1];
    	    c=coef6[2];
    	    d=coef6[3];
    	    e=coef6[4];
    	    f=coef6[5];
    	    g=coef6[6];
    	    h=coef6[7];
    	    i=coef6[8];
     	    j=coef6[9];
	    break;
	case 7: 	
   	    a=coef7[0];
    	    b=coef7[1];
    	    c=coef7[2];
    	    d=coef7[3];
    	    e=coef7[4];
    	    f=coef7[5];
    	    g=coef7[6];
    	    h=coef7[7];
    	    i=coef7[8];
     	    j=coef7[9];
	    break;
	case 8: 	
   	    a=coef8[0];
    	    b=coef8[1];
    	    c=coef8[2];
    	    d=coef8[3];
    	    e=coef8[4];
    	    f=coef8[5];
    	    g=coef8[6];
    	    h=coef8[7];
    	    i=coef8[8];
     	    j=coef8[9];
	    break;
	default:
	    std::cout << "number of octant exceeds 8 in intersectWithLine()!" << std::cout;
	    exit(-1);
	    break;
    }

    REAL domi=
	 e*e*p*p + 4*c*d*p*q - 4*a*c*q*q + 
	 f*f*q*q - 2*d*f*q*r + 
	 d*d*r*r - 
	 2*e*(f*p*q + d*p*r - 2*a*q*r) - 
	 4*b*(c*p*p + r*(-f*p + a*r));
    REAL x=
	(-(f*i*q*q) - 2*b*i*p*r + f*h*q*r + d*i*q*r + 
	 2*b*g*r*r - d*h*r*r - e*e*p*s - 
	 2*b*f*r*s - 2*c*(h*p*q - g*q*q - 2*b*p*s + 
	 d*q*s) + e*(i*p*q + h*p*r - 2*g*q*r + f*q*s + 
	 d*r*s))/domi;
    REAL y=
	(f*i*p*q - 2*f*h*p*r + d*i*p*r + f*g*q*r - 
	 2*a*i*q*r - d*g*r*r + 2*a*h*r*r - 
	 f*f*q*s + d*f*r*s + 
	 2*c*(h*p*p - g*p*q - d*p*s + 2*a*q*s) + 
	 e*(-i*p*p + g*p*r + f*p*s - 2*a*r*s))/domi;
    REAL z=
	(f*h*p*q - 2*d*i*p*q - f*g*q*q + 
	 2*a*i*q*q + d*h*p*r + d*g*q*r - 2*a*h*q*r + 
	 d*f*q*s - d*d*r*s + 
	 e*(-h*p*p + g*p*q + d*p*s - 2*a*q*s) + 
	 2*b*(i*p*p - g*p*r - f*p*s + 2*a*r*s))/domi;
    ptnp=vec(x,y,z);

    REAL val=a*x*x+b*y*y+c*z*z+d*x*y+e*y*z+f*x*z+g*x+h*y+i*z+j;

    if (val >= 0) // not intersect
	return false;
    else  // intersect
	return true;
}


void particle::planeRBForce(plnrgd_bdry<particle>* plb,	// August 19, 2013. apply moment on the mass center. August 22, 2013. September 6, 2013
			    std::map<int,std::vector<boundarytgt> >& BdryTgtMap,
			    std::vector<boundarytgt>& vtmp,
			    REAL &penetr){
	// (p,q,r) are in the same direction as the outward normal vector,
	// hence it is not necessary to provide information about which side the particle is about the plane.
	REAL p,q,r,s;
	BdryCoef tmp1=*((plb->CoefOfLimits).begin());
	p=tmp1.dirc.getx();
	q=tmp1.dirc.gety();
	r=tmp1.dirc.getz();
	s=-tmp1.dirc%tmp1.apt;  // plane equation: p(x-x0)+q(y-y0)+r(z-z0)=0, that is, px+qy+rz+s=0

	// judge which octant of this particle should be used. The result will be sued in GlobCoef() and getRadius()
	int xsign, ysign, zsign;	// xsign: -1, aminus; 0, invalid; 1, aplus
	xsign=0; ysign=0; zsign=0;
	// get the global coordinates of the six points of this particle
	vec local_aplus = vec(aplus, 0, 0);
	vec local_aminus = vec(-aminus, 0, 0);
	vec local_bplus = vec(0, bplus, 0);
	vec local_bminus = vec(0, -bminus, 0);
	vec local_cplus = vec(0, 0, cplus);
	vec local_cminus = vec(0, 0, -cminus);
	vec global_aplus = curr_position+globalVec(local_aplus);
	vec global_aminus = curr_position+globalVec(local_aminus);
	vec global_bplus = curr_position+globalVec(local_bplus);
	vec global_bminus = curr_position+globalVec(local_bminus);
	vec global_cplus = curr_position+globalVec(local_cplus);
	vec global_cminus = curr_position+globalVec(local_cminus);

	// cases for six different boundaries
	// use plb->distToBdry to deal with the boundaries that are not normal to the axles
	if( plb->distToBdry(global_aplus) > plb->distToBdry(global_aminus) )
	    xsign = 1;	// notice that the dist will be negative if inside, and positive outside,
	else		// so dist(aplus)>dist(aminus) means aplus is closer to the boundary than aminus
	    xsign =-1;	// that is why ">" is used
	if( plb->distToBdry(global_bplus) > plb->distToBdry(global_bminus) )
	    ysign = 1;
	else
	    ysign =-1;
	if( plb->distToBdry(global_cplus) > plb->distToBdry(global_cminus) )
	    zsign = 1;
	else
	    zsign =-1;

/*
	if(p>0 && q==0 && r==0){	// right boundary, positive x
		if(global_aplus.getx()>global_aminus.getx())
			xsign = 1;
		else
			xsign = -1;
		if(global_bplus.getx()>global_bminus.getx())
			ysign = 1;
		else 
			ysign = -1;
		if(global_cplus.getx()>global_cminus.getx())
			zsign = 1;
		else 
			zsign = -1;
	}
	else if(p<0 && q==0 && r==0){	// left boundary, negative x
		if(global_aplus.getx()<global_aminus.getx())
			xsign = 1;
		else
			xsign = -1;
		if(global_bplus.getx()<global_bminus.getx())
			ysign = 1;
		else 
			ysign = -1;
		if(global_cplus.getx()<global_cminus.getx())
			zsign = 1;
		else 
			zsign = -1;
	}
	else if(p==0 && q>0 && r==0){	// back boundary, positive y
		if(global_aplus.gety()>global_aminus.gety())
			xsign = 1;
		else 
			xsign = -1;
		if(global_bplus.gety()>global_bminus.gety())
			ysign = 1;
		else 
			ysign = -1;
		if(global_cplus.gety()>global_cminus.gety())
			zsign = 1;
		else
			zsign = -1;

	}
	else if(p==0 && q<0 && r==0){	// front boundary, negative y
		if(global_aplus.gety()<global_aminus.gety())
			xsign = 1;
		else 
			xsign = -1;
		if(global_bplus.gety()<global_bminus.gety())
			ysign = 1;
		else 
			ysign = -1;
		if(global_cplus.gety()<global_cminus.gety())
			zsign = 1;
		else
			zsign = -1;

	}
	else if(p==0 && q==0 && r>0){	// top boundary, positive z
		if(global_aplus.getz()>global_aminus.getz())
			xsign = 1;
		else
			xsign = -1;
		if(global_bplus.getz()>global_bminus.getz())
			ysign = 1;
		else 
			ysign = -1;
		if(global_cplus.getz()>global_cminus.getz())
			zsign = 1;
		else
			zsign = -1;
	}
	else if(p==0 && q==0 && r<0){	// bottom boundary, negative z
		if(global_aplus.getz()<global_aminus.getz())
			xsign = 1;
		else
			xsign = -1;
		if(global_bplus.getz()<global_bminus.getz())
			ysign = 1;
		else 
			ysign = -1;
		if(global_cplus.getz()<global_cminus.getz())
			zsign = 1;
		else
			zsign = -1;
	}
	else{	// boundary is not parallel to the three axel faces, not deal with this case at present. August 22, 2013
		std::cout << "Boundary is not parallel to the three axel faces" << std::endl;
		exit(-1);
	}
*/
	// get octant number
	int num_oct = 1; 
	if(xsign==1 && ysign==1 && zsign==1)
	    num_oct = 1;
	if(xsign==-1 && ysign==1 && zsign==1)
	    num_oct = 2;
	if(xsign==-1 && ysign==-1 && zsign==1)
	    num_oct = 3;
	if(xsign==1 && ysign==-1 && zsign==1)
	    num_oct = 4;
	if(xsign==1 && ysign==1 && zsign==-1)
	    num_oct = 5;
	if(xsign==-1 && ysign==1 && zsign==-1)
	    num_oct = 6;
	if(xsign==-1 && ysign==-1 && zsign==-1)
	    num_oct = 7;
	if(xsign==1 && ysign==-1 && zsign==-1)
	    num_oct = 8;

	vec pt1;
	if (!nearestPTOnPlane(p, q, r, s, pt1, num_oct)) // the particle and the plane does not intersect
	  return;

	// if particle and plane intersect:
	cntnum++;
	inContact = true;
	vec dirc=normalize(vec(p,q,r));
	vec rt[2];
	if (!intersectWithLine(pt1,dirc,rt, num_oct))    // the line and ellipsoid surface does not intersect
	    return;

	vec pt2;
	/*
	if (p*rt[0].getx()+q*rt[0].gety()+r*rt[0].getz()+s > 0)
		pt2=rt[0];
	else
		pt2=rt[1];
	*/
	if (vfabs(rt[0]-pt1) < vfabs(rt[1]-pt1) )
	    pt2 = rt[0];
	else
	    pt2 = rt[1];


	// obtain normal force
	REAL penetration=vfabs(pt1-pt2);
	if (penetration / (2.0*getRadius(pt2, xsign, ysign, zsign) ) <= MINOVERLAP){
	    return;
	}

/*
std::cout << std::endl;
std::cout << pt2.getx() << "   " << pt2.gety() << "   " << pt2.getz() << " ;" << std::endl << std::setw(16) << aplus << std::setw(16) << aminus << std::setw(16) << bplus << std::setw(16) << bminus << std::setw(16) << cplus << std::setw(16) << cminus << std::setw(16) << getCurrPosition().getx() << std::setw(16) << getCurrPosition().gety() << std::setw(16) << getCurrPosition().getz() << std::setw(16) << getCurrDirecA().getx() << std::setw(16) << getCurrDirecA().gety() << std::setw(16) << getCurrDirecA().getz() << std::setw(16) << getCurrDirecB().getx() << std::setw(16) << getCurrDirecB().gety() << std::setw(16) << getCurrDirecB().getz() << std::setw(16) << getCurrDirecC().getx() << std::setw(16) << getCurrDirecC().gety() << std::setw(16) << getCurrDirecC().getz() << std::endl;
*/

	REAL R0=getRadius(pt2, xsign, ysign, zsign);
//std::cout << "Radius with octant in planeRBForce():    " << R0 << std::endl;
//REAL	R0=getRadius(pt2);
//std::cout << "Radius without octant in planeRBForce(): " << R0 << std::endl;
	REAL E0=YOUNG/(1-POISSON*POISSON); // rigid wall has infinite YOUNG's modulus
	REAL allowedOverlap = 2.0 * R0 * MAXOVERLAP;
	if (penetration > allowedOverlap) {
	  g_debuginf << "particle.cpp: g_iter=" << g_iteration 
		     << " ptclId=" << getID()
		     << " bdryId=" << plb->bdry_id
		     << " penetr=" << penetration 
		     << " allow=" << allowedOverlap << std::endl;
	  penetration = allowedOverlap;
	}

#ifdef MEASURE_EPS
	penetration = nearbyint (penetration/MEPS) * MEPS;
#endif
	REAL contact_radius=sqrt(penetration*R0);
	penetr = penetration;
/*
std::cout << "penetration in planeRBbdryForce(): " << penetration << std::endl;
std::cout << "R0 in planeRBbdryForce(): " << R0 << std::endl;
std::cout << "E0 in planeRBbdryForce(): " << E0 << std::endl;
std::cout << "NormDirc in planeRBbdryForce(): " << -dirc.getx() << ", " << -dirc.gety() << ", " << -dirc.getz() << std::endl;
*/
	vec NormDirc=-dirc; //normalize(pt1-pt2);
	vec NormalForce=sqrt(penetration*penetration*penetration)*sqrt(R0)*4*E0/3*NormDirc; // pow(penetration,1.5), a serious bug


	/*
	g_debuginf<<" "<<g_iteration
		  <<" "<<getID()
		  <<" "<<plb->bdry_id
		  <<" "<<pt1.getx()
		  <<" "<<pt1.gety()
		  <<" "<<pt1.getz()
		  <<" "<<rt[0].getx()
		  <<" "<<rt[0].gety()
		  <<" "<<rt[0].getz()
		  <<" "<<rt[1].getx()
		  <<" "<<rt[1].gety()
		  <<" "<<rt[1].getz()
		  <<" "<<vfabs(rt[0]-pt1)
		  <<" "<<vfabs(rt[1]-pt1)
		  <<" "<<penetration
		  <<std::endl;
	*/

	// apply normal force
	addForce(NormalForce);
	addMoment(((pt1+pt2)*0.5-curr_center_mass)*NormalForce);

        // add force terms to the average stress
	// April 23, 2014
    	// now the average particle stress is calculated in the local coordinate
    	// system initially, refer to my own reserch notes in page 141, October 31, 2014
  	vec localNormalForce = localVec(NormalForce);	// convert NormalForce to local coordinate
	matrix tmp_force(3,1); 
	tmp_force(1,1) = localNormalForce.getx(); tmp_force(2,1) = localNormalForce.gety(); tmp_force(3,1) = localNormalForce.getz();

	vec localPosi = localVec((pt1+pt2)*0.5);
	matrix tmp_posi(1,3);
	tmp_posi(1,1) = localPosi.getx(); tmp_posi(1,2) = localPosi.gety(); tmp_posi(1,3) = localPosi.getz();
	matrix tmp_stress = tmp_force*tmp_posi;
	addStress(tmp_stress);
	
	// Here, calculate the maximum tensile stress for this contact and store to particle p1
	// refer to Theory of Elasticiy, Timoshenko and Goodier
	REAL tmp_tensileStress;	// this is the maximum tensile stress at this contact
	REAL tmp_q0 = 3.0*vfabs(NormalForce)/(2.0*PI*contact_radius*contact_radius);	// q0 is the maximum normal pressure
	tmp_tensileStress = (1.0-2.0*POISSON)/3.0*tmp_q0;

	if(tmp_tensileStress>ContactTensileCritical){	// this contact point is critical point
	    addMaximuContactTensile(tmp_tensileStress, (pt1+pt2)*0.5);
    	}

	// obtain normal damping force
	vec veloc2 = getCurrVelocity() + getCurrOmga()*((pt1+pt2)*0.5-getCurrCenterMass());
	REAL kn = pow(6*vfabs(NormalForce)*R0*pow(E0,2),1.0/3.0);
	REAL DMP_CRTC = 2*sqrt(getMass()*kn); // critical damping
	vec CntDampingForce  = DMP_CNT * DMP_CRTC * ((-veloc2)%NormDirc)*NormDirc;

	// apply normal damping force
	addForce(CntDampingForce);
	addMoment(((pt1+pt2)*0.5-curr_center_mass)*CntDampingForce);

//	vec localCntForce = localVec(CntDampingForce);
//	tmp_force(1,1) = localCntForce.getx(); tmp_force(2,1) = localCntForce.gety(); tmp_force(3,1) = localCntForce.getz();
//	tmp_stress = tmp_force*tmp_posi;
//	addStress(tmp_stress);

	vec TgtForce = 0;
	if (BDRYFRIC != 0){
	    // checkin previous tangential force and displacement
	    vec PreTgtForce;
	    vec PreTgtDisp;
	    bool PreTgtLoading=true;
	    vec  TgtDispStart;
	    REAL TgtPeak=0;

	    bool TgtLoading=true;
	    std::vector<boundarytgt>::iterator it;
	    for (it=BdryTgtMap[plb->bdry_id].begin();it!=BdryTgtMap[plb->bdry_id].end();++it){
		if (ID == it->ptcl) {
		    PreTgtForce  =it->TgtForce;
		    PreTgtDisp   =it->TgtDisp;
		    PreTgtLoading=it->TgtLoading;
		    TgtDispStart =it->TgtDispStart;
		    TgtPeak      =it->TgtPeak;
		    break;
		}
	    }
		
	    // obtain tangtential force
	    REAL G0 = 0.5*YOUNG/(1+POISSON);
	    // Vr = Vb + w (crossdot) r, each item needs to be in either global or local frame; 
	    //      here global frame is used for better convenience.
	    vec RelaDispInc = (curr_velocity+curr_omga*((pt1+pt2)*0.5-curr_center_mass))*TIMESTEP;  
	    vec TgtDispInc = RelaDispInc-(RelaDispInc%NormDirc)*NormDirc;
	    vec TgtDisp    = PreTgtDisp + TgtDispInc; // PreTgtDisp read by checkin
	    vec TgtDirc;

	    if (vfabs(TgtDisp) == 0)
		TgtDirc = 0;
	    else
		TgtDirc = normalize(-TgtDisp); // TgtDirc points along tangential forces exerted on particle 1

	    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	    // linear friction model
	    REAL fP  = BDRYFRIC*vfabs(NormalForce);
	    REAL ks  = 4*G0*contact_radius/(2-POISSON);
	    TgtForce = PreTgtForce + ks*(-TgtDispInc); // PreTgtForce read by checkin

	    vec FricDampingForce = 0;
	    if (vfabs(TgtForce) > fP) // slide case
		TgtForce = fP*TgtDirc;
	    else { // adhered/slip case
		
		// obtain tangential damping force
		vec RelaVel = curr_velocity + curr_omga*((pt1+pt2)*0.5-curr_center_mass);  
		vec TgtVel  = RelaVel - (RelaVel%NormDirc)*NormDirc;
		REAL DMP_CRTC = 2*sqrt(getMass()*ks); // critical damping
		FricDampingForce = 1.0 * DMP_CRTC * (-TgtVel);
	    }

	    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	    // Mindlin's model (loading/unloading condition assumed)
	    // This model is not recommended as it is impossible to strictly determine loading/unloading condition
            // unless load is known (the case of pure moment rotation).
#ifdef MINDLIN_ASSUMED
	    REAL val = 0;
	    fP = FRICTION*vfabs(NormalForce);
	    TgtLoading = (PreTgtDisp%TgtDispInc >= 0); 
	    
	    if (TgtLoading) {              // loading
		if (!PreTgtLoading) {      // pre-step is unloading
		    val = 8*G0*contact_radius*vfabs(TgtDispInc)/(3*(2-POISSON)*fP);
		    TgtDispStart = PreTgtDisp;
		}
		else                       // pre-step is loading
		    val = 8*G0*contact_radius*vfabs(TgtDisp-TgtDispStart)/(3*(2-POISSON)*fP);
		
		if (val > 1.0)              
		    TgtForce = fP*TgtDirc;
		else {
		    ks = 4*G0*contact_radius/(2-POISSON)*sqrt(1-val);
		    //incremental method
		    TgtForce = PreTgtForce + ks*(-TgtDispInc); // TgtDispInc determines signs
		    //total value method: TgtForce = fP*(1-pow(1-val, 1.5))*TgtDirc;
		}
	    }
	    else {                         // unloading
		if (PreTgtLoading) {       // pre-step is loading
		    val = 8*G0*contact_radius*vfabs(TgtDisp-TgtDispStart)/(3*(2-POISSON)*fP);
		    TgtPeak = vfabs(PreTgtForce);
		}
		else                       // pre-step is unloading
		    val = 8*G0*contact_radius*vfabs(TgtDisp-TgtDispStart)/(3*(2-POISSON)*fP);
		
		if (val > 1.0 || TgtPeak > fP)  
		    TgtForce = fP*TgtDirc;
		else {
		    ks = 2*sqrt(2)*G0*contact_radius/(2-POISSON) * sqrt(1+pow(1-TgtPeak/fP,2.0/3.0)+val);
		    //incremental method
		    TgtForce = PreTgtForce + ks*(-TgtDispInc); // TgtDispInc determines signs
		    //total value method: TgtForce = (TgtPeak-2*fP*(1-sqrt(2)/4*pow(1+ pow(1-TgtPeak/fP,2.0/3.0) + val,1.5)))*TgtDirc;
		}
	    }
	    
	    if (vfabs(TgtForce) > fP) // slice case
		TgtForce = fP*TgtDirc;
	    else { // adhered/slip case
		
		// obtain tangential damping force
		vec RelaVel = curr_velocity + curr_omga*((pt1+pt2)*0.5-curr_center_mass);  
		vec TgtVel  = RelaVel - (RelaVel%NormDirc)*NormDirc;
		REAL DMP_CRTC = 2*sqrt(getMass()*ks); // critical damping
		FricDampingForce = 1.0 * DMP_CRTC * (-TgtVel);
	    }

#endif
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////	

	    /*
	    if (g_iteration%100==0)  
	    g_debuginf<<"particle.cpp, g_iteraion="<<g_iteration
		      <<" NormalForce="<<vfabs(NormalForce)
		      <<" CntDampingForce= "<<vfabs(CntDampingForce)
		      <<" kn="<<kn
		      <<" TgtForce="<<vfabs(TgtForce)
		      <<" FricDampingForce="<<vfabs(FricDampingForce)
		      <<" ks="<<ks
		      <<std::endl;
	    */
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////


	    // apply tangential force
	    addForce(TgtForce);
	    addMoment(((pt1+pt2)*0.5-curr_center_mass)*TgtForce); 

	    vec localTgtForce = localVec(TgtForce);
	    tmp_force(1,1) = localTgtForce.getx(); tmp_force(2,1) = localTgtForce.gety(); tmp_force(3,1) = localTgtForce.getz();
	    tmp_stress = tmp_force*tmp_posi;
	    addStress(tmp_stress);

	    // apply tangential damping force for adhered/slip case
	    addForce(FricDampingForce);

	    // update current tangential force and displacement, don't checkout.
	    // checkout in rigidBF() ensures BdryTgtMap update after each particles
            // contacting this boundary is processed.
	    vtmp.push_back(boundarytgt(ID,TgtForce,TgtDisp,TgtLoading,TgtDispStart,TgtPeak));

	}

	plb->normal -= NormalForce;
	plb->tangt  -= TgtForce;
//	plb->moment-=(((pt1+pt2)/2-curr_position)*NormForce+
//		      ((pt1+pt2)/2-curr_position)*TgtForce));
}


vec particle::cylinderRBForce(int bdry_id, const cylinder& S, int side){	// August 19, 2013. August 22, 2013. September 6, 2013
	//side=-1, the particles are inside the cylinder
	//side=+1, the particles are outside the cylinder
	REAL x0=S.getCenter().getx();
	REAL y0=S.getCenter().gety();
	REAL r=S.getRadius();
	REAL coef_bdry[10]={1,1,0,0,0,0,-2*x0,-2*y0,0,x0*x0+y0*y0-r*r};

	// determine which octant of this particle should be used
	int xsign, ysign, zsign;
	xsign=0; ysign=0; zsign=0;
	// get the global coordinates of the six points of this particle
	vec local_aplus = vec(aplus, 0, 0);
	vec local_aminus = vec(-aminus, 0, 0);
	vec local_bplus = vec(0, bplus, 0);
	vec local_bminus = vec(0, -bminus, 0);
	vec local_cplus = vec(0, 0, cplus);
	vec local_cminus = vec(0, 0, -cminus);
	vec global_aplus = curr_position+globalVec(local_aplus);
	vec global_aminus = curr_position+globalVec(local_aminus);
	vec global_bplus = curr_position+globalVec(local_bplus);
	vec global_bminus = curr_position+globalVec(local_bminus);
	vec global_cplus = curr_position+globalVec(local_cplus);
	vec global_cminus = curr_position+globalVec(local_cminus);
	// get the distances of these six points to the project of the cylinder center, (x0, y0, 0)
	REAL dist_aplus = pow(global_aplus.getx()-x0, 2) + pow(global_aplus.gety()-y0,2);
	REAL dist_aminus = pow(global_aminus.getx()-x0, 2) + pow(global_aminus.gety()-y0, 2);
	REAL dist_bplus = pow(global_bplus.getx()-x0, 2) + pow(global_bplus.gety()-y0,2);
	REAL dist_bminus = pow(global_bminus.getx()-x0, 2) + pow(global_bminus.gety()-y0, 2);
	REAL dist_cplus = pow(global_cplus.getx()-x0, 2) + pow(global_cplus.gety()-y0,2);
	REAL dist_cminus = pow(global_cminus.getx()-x0, 2) + pow(global_cminus.gety()-y0, 2);
	// judge which octant based on these distances
	if(dist_aplus > dist_aminus)
		xsign = 1;
	else
		xsign = -1;
	if(dist_bplus > dist_bminus)
		ysign = 1;
	else
		ysign = -1;
	if(dist_cplus > dist_cminus)
		zsign = 1;
	else
		zsign = -1;
//	// calculate the coefficients for intersectWithLine()
//	GlobCoef(xsign, ysign, zsign);
//GlobCoef();

	// get octant number
	int num_oct = 1; 
	vec pt1;
	if(xsign==1 && ysign==1 && zsign==1){
	    num_oct = 1;
	    if (!root6(coef1,coef_bdry,pt1))  //on the cylinder and within the particle
	        return 0; //no contact
	}
	if(xsign==-1 && ysign==1 && zsign==1){
	    num_oct = 2;
	    if (!root6(coef2,coef_bdry,pt1))  //on the cylinder and within the particle
	        return 0; //no contact
	}
	if(xsign==-1 && ysign==-1 && zsign==1){
	    num_oct = 3;
	    if (!root6(coef3,coef_bdry,pt1))  //on the cylinder and within the particle
	        return 0; //no contact
	}
	if(xsign==1 && ysign==-1 && zsign==1){
	    num_oct = 4;
	    if (!root6(coef4,coef_bdry,pt1))  //on the cylinder and within the particle
	        return 0; //no contact
	}
	if(xsign==1 && ysign==1 && zsign==-1){
	    num_oct = 5;
	    if (!root6(coef5,coef_bdry,pt1))  //on the cylinder and within the particle
	        return 0; //no contact
	}
	if(xsign==-1 && ysign==1 && zsign==-1){
	    num_oct = 6;
	    if (!root6(coef6,coef_bdry,pt1))  //on the cylinder and within the particle
	        return 0; //no contact
	}
	if(xsign==-1 && ysign==-1 && zsign==-1){
	    num_oct = 7;
	    if (!root6(coef7,coef_bdry,pt1))  //on the cylinder and within the particle
	        return 0; //no contact
	}
	if(xsign==1 && ysign==-1 && zsign==-1){
	    num_oct = 8;
	    if (!root6(coef8,coef_bdry,pt1))  //on the cylinder and within the particle
	        return 0; //no contact
	}

	cntnum++;
	vec rt[2];
	vec cz=vec(S.getCenter().getx(),S.getCenter().gety(),pt1.getz());
	vec tmp=pt1-cz;
	intersectWithLine(pt1, normalize(tmp),rt, num_oct);
	vec pt2;

	if ((rt[0]-pt1)%tmp*side<0)
		pt2=rt[0];
	else
		pt2=rt[1];
	//vec pt2=vfabs(rt[0]-cz)>vfabs(rt[1]-cz)?rt[0]:rt[1];
	REAL radius=getRadius(pt2, xsign, ysign, zsign);//pt2.print();std::cout<<radius;getchar();
//	REAL radius=getRadius(pt2);
	REAL E0=0.5*YOUNG/(1-POISSON*POISSON);
	REAL R0=(r*radius)/(r+radius);
	REAL rou=vfabs(pt1-pt2);
	vec NormDirc=normalize(pt1-pt2);
	REAL nfc=sqrt(rou*rou*rou)*sqrt(R0)*4*E0/3.0; // pow(rou,1.5), a serious bug
	vec NormalForce=nfc*NormDirc;

	addForce(NormalForce);
	addMoment(((pt1+pt2)*0.5-getCurrCenterMass())*NormalForce);	    

        // add force terms to the average stress
	// April 23, 2014
	vec localNormalForce = localVec(NormalForce);
	matrix tmp_force(3,1); 
	tmp_force(1,1) = localNormalForce.getx(); tmp_force(2,1) = localNormalForce.gety(); tmp_force(3,1) = localNormalForce.getz();

	vec localPosi = localVec((pt1+pt2)*0.5);
	matrix tmp_posi(1,3);
	tmp_posi(1,1) = localPosi.getx(); tmp_posi(1,2) = localPosi.gety(); tmp_posi(1,3) = localPosi.getz();
	matrix tmp_stress = tmp_force*tmp_posi;
	addStress(tmp_stress);
	
	return NormalForce;
}


bool particle::isCLongEnough() const{	// if c^+ + c^- > 1.2max(A_max, B_max)

    REAL amax = (aplus>aminus?aplus:aminus);
    REAL bmax = (bplus>bminus?bplus:bminus);

    REAL right_value = amax>bmax?amax:bmax;
    right_value = 1.2*right_value;

    return cplus+cminus > right_value;

} // end isCLongEnough()

bool particle::isBLongEnough() const{	// if b^+ + b^- > 1.2max(A_max, C_max)

    REAL amax = (aplus>aminus?aplus:aminus);
    REAL cmax = (cplus>cminus?cplus:cminus);

    REAL right_value = amax>cmax?amax:cmax;
    right_value = 1.2*right_value;

    return bplus+bminus > right_value;

} // end isBLongEnough()

bool particle::isALongEnough() const{	// if a^+ + a^- > 1.2max(B_max, C_max)

    REAL bmax = (bplus>bminus?bplus:bminus);
    REAL cmax = (cplus>cminus?cplus:cminus);

    REAL right_value = bmax>cmax?bmax:cmax;
    right_value = 1.2*right_value;

    return aplus+aminus > right_value;

} // end isCLongEnough()


// determine which plane is the break plane based on the particle average stress 
// as the method shown in the week summary of 2014_09_19.
// return value: 1--ab plane, 2--ac plane, 3 --bc plane
int particle::getBreakPlane(){

/*
    matrix trans_Q(3,3); // the transformation tensor
    trans_Q(1,1) = cosl(curr_direction_a.getx()); 
    trans_Q(1,2) = cosl(curr_direction_a.gety()); 
    trans_Q(1,3) = cosl(curr_direction_a.getz()); 

    trans_Q(2,1) = cosl(curr_direction_b.getx()); 
    trans_Q(2,2) = cosl(curr_direction_b.gety()); 
    trans_Q(2,3) = cosl(curr_direction_b.getz()); 

    trans_Q(3,1) = cosl(curr_direction_c.getx()); 
    trans_Q(3,2) = cosl(curr_direction_c.gety()); 
    trans_Q(3,3) = cosl(curr_direction_c.getz()); 

    matrix local_stress = trans_Q*average_stress*trans_Q.getTrans();
*/

    // now the average particle stress is calculated in the local coordinate
    // system initially, refer to my own reserch notes in page 141, October 31, 2014
    REAL sigma_aa = average_stress(1,1);
    REAL sigma_bb = average_stress(2,2);
    REAL sigma_cc = average_stress(3,3);

    if(sigma_aa >= sigma_bb && sigma_aa >= sigma_cc){  // sigma_aa is the largest
	return 3;	// plane-bc
    }
    else if(sigma_bb >= sigma_aa && sigma_bb >= sigma_cc){  // sigma_bb is the largest
	return 2;	// plane-ac
    }
    else if(sigma_cc >= sigma_aa && sigma_cc >= sigma_bb){  // sigma_cc is the largest
	return 1;	// plane-ab
    }

} // end getBreakPlane()


// determine which plane is the break plane
// the input value means that axle is in bad length ratio, which should be
// one axle of the break plane
// input: 1--axle a, 2--axle b, 3--axle c
// return: 1--ab plane, 2--ac plane, 3--bc plane
int particle::getBreakPlane(int abc){

/*
    matrix trans_Q(3,3); // the transformation tensor
    trans_Q(1,1) = cosl(curr_direction_a.getx()); 
    trans_Q(1,2) = cosl(curr_direction_a.gety()); 
    trans_Q(1,3) = cosl(curr_direction_a.getz()); 

    trans_Q(2,1) = cosl(curr_direction_b.getx()); 
    trans_Q(2,2) = cosl(curr_direction_b.gety()); 
    trans_Q(2,3) = cosl(curr_direction_b.getz()); 

    trans_Q(3,1) = cosl(curr_direction_c.getx()); 
    trans_Q(3,2) = cosl(curr_direction_c.gety()); 
    trans_Q(3,3) = cosl(curr_direction_c.getz()); 

    matrix local_stress = trans_Q*average_stress*trans_Q.getTrans();
*/

    // now the average particle stress is calculated in the local coordinate
    // system initially, refer to my own reserch notes in page 141, October 31, 2014
    REAL sigma_aa = average_stress(1,1);
    REAL sigma_bb = average_stress(2,2);
    REAL sigma_cc = average_stress(3,3);

    switch (abc){
	case 1:	// axle a must be in the break plane
	    if(isBLongerThanC()){  // bplus+bminus>1.2*(cplus+cminus)
		return 2;	// break along plane-ac
	    }
	    else if(isCLongerThanB()){  // cplus+cminus>1.2*(bplus+bminus)
		return 1;	// break along plane-ab
	    }
	    else{
		if(sigma_bb>sigma_cc)	// axle c should be in
		    return 2;	// break along plane-ac
		else
		    return 1;	// break along plane-ab
	    }

	    break;
	case 2:	// axle b must be in the break plane
	    if(isALongerThanC()){  // Aplus+Aminus>1.2*(cplus+cminus)
		return 3;
	    }
	    else if(isCLongerThanA()){  // cplus+cminus>1.2*(aplus+aminus)   
		return 1;
	    }
	    else{
	  	if(sigma_aa>sigma_cc)	// axle c should be in
		    return 3;
		else
		    return 1;
 	    }

	    break;
	case 3:	// axle c must be in the break plane
	    if(isALongerThanB()){  // Aplus+Aminus>1.2*(bplus+bminus)
		return 3;	// break along plane-bc
	    }
	    else if(isBLongerThanA()){  // bplus+bminus>1.2*(aplus+aminus)
		return 2;	// break along plane-ac
	    }
	    else{
		if(sigma_aa>sigma_bb)	// axle b should be in
		    return 3;
		else
		    return 2;
	    }
	
	    break;
	default:
	    std::cout << "break plane type should be 1/2/3..." << std::endl;
	    exit(-1);

    }

} // end getBreakPlane()

// clear average stress after the end of each step, April 23, 2014, also clear contact stresses
void particle::clearStress(){
    matrix zero3x3(3,3); 
    average_stress = zero3x3;

    numCriticalContacts = 0;	
    stress1 = 0;	// it is very important to initialize the stress at every timestep
    stress2 = 0;
    stress3 = 0;	

    dem::vec tmp_zero = dem::vec(0,0,0);
    contact1 = tmp_zero;	
    contact2 = tmp_zero;	
    contact3 = tmp_zero;

//    initTypeBroken();	// this makes initialization of typeBroken at each step, which makes
			// typeBroken represent the broken state at this current step

} // clearStress()

void particle::addMaximuContactTensile(REAL tmp_stress, dem::vec tmp_contact){
    // check if this critical contact point is larger than the existing three ones
    if(tmp_stress>stress1){	// this tensile stress is larger than the current largest one
	numCriticalContacts++;
	stress3=stress2; contact3=contact2;
	stress2=stress1; contact2=contact1;
	stress1=tmp_stress; contact1=tmp_contact;
    }
    else if(tmp_stress>stress2){	// this tensile stress is less than the current largest one, but larger than
					// current middle one
	numCriticalContacts++;
	stress3=stress2; contact3=contact2;
	stress2=tmp_stress; contact2=tmp_contact;
    }
    else if(tmp_stress>stress3){	// this tensile stress is less than the current middle one, but larger than 
					// the current smallest one
	numCriticalContacts++;
	stress3=tmp_stress; contact3=tmp_contact;
    }// this algorithm can also handle stress1==stress2==stress3

} // addMaximumContactTensile()


} // namespace dem ends
