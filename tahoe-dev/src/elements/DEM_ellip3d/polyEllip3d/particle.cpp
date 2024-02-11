#include "particle.h"
#include "parameter.h"
#include "ran.h"
#include "root6.h"
#include <iostream>

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
//    GlobCoef(1, 1, 1);	// calculate coefficients when it is needed, August 19, 2013
    GlobCoef();	// September 6, 2013
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



std::cout << "particle ID: " << ID << ", center of mass: " << std::setprecision(16) << curr_center_mass.getx() << ", " << std::setprecision(16) << curr_center_mass.gety() << ", " << std::setprecision(16) << curr_center_mass.getz() << std::endl;


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
//vec J_poly = J; 
//REAL mass_poly = mass;
//REAL volume_poly = volume;
//    inContact=false;
//volume=4/3.0*PI*aplus*bplus*cplus;
//mass=density*volume;
//J=vec(mass/5*(bplus*bplus+cplus*cplus),mass/5*(aplus*aplus+cplus*cplus),mass/5*(aplus*aplus+bplus*bplus));
//if(mass_poly != mass){
//	cout << "mass_ploy != mass" << std::endl;
//}
//if(volume_poly != volume){
//	cout << "volume_ploy != volume" << std::endl;
//}
//if(J_poly != J){
//	cout << "J_ploy != J" << std::endl;
//}
//    GlobCoef(1, 1, 1);	// calculate coefficients when it is needed, August 16, 2013
    GlobCoef();	// September 6, 2013
//std::cout << "test in particle.cpp particle(): " << std::endl;
//std::cout << "mass: " << mass << std::endl;
//std::cout << "volume: " << volume << std::endl;
//std::cout << "J: " <<std::endl; J.print();
//std::cout << "local_center_mass: " << std::endl; local_center_mass.print();


}


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
//std::cout << "max radius is (getMaxRadius() in particle.cpp): " << temp << std::endl;
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
//std::cout << "min radius is (getMinRadius() in particle.cpp): " << temp << std::endl;
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
	return;
    }

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
	return;
    }

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
	return;
    }

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
	return;
    }

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
	return;
    }

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
	return;
    }

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
	return;
    }

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

    // the eightth octant
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
	return;
    }

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


// central difference integration method
void particle::update() {	// August 16, 2013

    if (getType()==0 || getType()==5) { // 0-free, 1-fixed, 5-free bounary particle
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
	
//vec curr_omga_test = globalVec(curr_local_omga);
//if(curr_omga_test != curr_omga) {
//	std::cout << "curr_omga_test != curr_omga in particle.cpp update()!\n The globalVec() is not correct!" << std::endl;
//	std::cout << "curr_omga"; curr_omga.print();
//	std:;cout << "curr_omga_test"; curr_omga_test.print();
//}	

	curr_direction_a=vacos(normalize(rotateVec(vcos(prev_direction_a),curr_omga*TIMESTEP)));
	curr_direction_b=vacos(normalize(rotateVec(vcos(prev_direction_b),curr_omga*TIMESTEP)));
	curr_direction_c=vacos(normalize(rotateVec(vcos(prev_direction_c),curr_omga*TIMESTEP)));

	// force: translational kinetics equations are in global frame
	curr_velocity=prev_velocity*(2-atf)/(2+atf)+force/(mass*MASS_SCL)*TIMESTEP*2/(2+atf);
//curr_position = prev_position + curr_velocity*TIMESTEP;
	curr_center_mass = prev_center_mass + curr_velocity*TIMESTEP;	// August 16, 2013
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
//std::cout << c << "   " << e << "   " << h << "   " << j << std::endl;
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

//std::cout << "x y z of the boundary: " << p << "  " << q << "  " << r << std::endl;

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

// used to debug
//std::cout << "boundary contanct octant in planeRBForce(): " << xsign << "   " << ysign << "    " << zsign << std::endl;
//xsign = 1; ysign = 1; zsign = 1;
//REAL test_coef[10];
//std::cout << "in planeBdryForce() particle.cpp: " << std::endl;
//getGlobCoef(test_coef);
//std::cout << "Coefficients before GlobCoef(): " << test_coef[0] << "  " << test_coef[2] << "  " << test_coef[5] << "  " << test_coef[7] << "  " << test_coef[9] << std::endl; 

//	// calcualte new coefficients for nearestPTOnPlane() and intersectWithLine()
//	GlobCoef(xsign, ysign, zsign);

//getGlobCoef(test_coef);
//std::cout << "Coefficients from GlobCoef(): " << test_coef[0] << "  " << test_coef[2] << "  " << test_coef[5] << "  " << test_coef[7] << "  " << test_coef[9] << std::endl; 
//std::cout << std::endl;

//GlobCoef();

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
//std::cout << pt1.getx() << "   " << pt1.gety() << "   " << pt1.getz() << std::endl;
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

//std::cout << pt2.getx() << "   " << pt2.gety() << "   " << pt2.getz() << std::endl;

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
	addMoment(((pt1+pt2)/2-curr_center_mass)*NormalForce);
	
	// obtain normal damping force
	vec veloc2 = getCurrVelocity() + getCurrOmga()*((pt1+pt2)/2-getCurrCenterMass());
	REAL kn = pow(6*vfabs(NormalForce)*R0*pow(E0,2),1.0/3.0);
	REAL DMP_CRTC = 2*sqrt(getMass()*kn); // critical damping
	vec CntDampingForce  = DMP_CNT * DMP_CRTC * ((-veloc2)%NormDirc)*NormDirc;

	// apply normal damping force
	addForce(CntDampingForce);
	addMoment(((pt1+pt2)/2-curr_center_mass)*CntDampingForce);

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
	    REAL G0 = YOUNG/2/(1+POISSON);
	    // Vr = Vb + w (crossdot) r, each item needs to be in either global or local frame; 
	    //      here global frame is used for better convenience.
	    vec RelaDispInc = (curr_velocity+curr_omga*((pt1+pt2)/2-curr_center_mass))*TIMESTEP;  
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
		vec RelaVel = curr_velocity + curr_omga*((pt1+pt2)/2-curr_center_mass);  
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
		vec RelaVel = curr_velocity + curr_omga*((pt1+pt2)/2-curr_center_mass);  
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
	    addMoment(((pt1+pt2)/2-curr_center_mass)*TgtForce); 

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
	REAL nfc=sqrt(rou*rou*rou)*sqrt(R0)*4*E0/3; // pow(rou,1.5), a serious bug
	vec NormalForce=nfc*NormDirc;

	addForce(NormalForce);
	addMoment(((pt1+pt2)/2-getCurrCenterMass())*NormalForce);	    
	
	return NormalForce;
}

} // namespace dem ends
