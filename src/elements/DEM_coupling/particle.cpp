#include <cmath>
#include <iostream>

#include "const.h"
#include "ran.h"
#include "root6.h"
#include "particle.h"

//#define MOMENT
#ifdef MOMENT
const int START = 10000;  // at which time step to apply moment? for moment rotation test only.
#define SLIP  // if defined, stick and slip; otherwise slide.
#endif
//ellip3d.cxx: dem::TIMESTEP = 5.0e-07; A.deposit(12000,... 

using namespace std;

namespace dem {

extern ofstream g_exceptioninf;
extern int g_iteration;
extern long idum;
extern long double GRVT_SCL;
extern long double DMP_F;
extern long double DMP_M;
extern long double MASS_SCL;
extern long double MNT_SCL;
extern long double PILE_RATE;
extern long double BDRYFRIC;

particle::particle(int n, int tp, vec center, long double r){
    type=tp;
    ID=n;
    long double tmp0=ran(&idum);
    long double tmp1=ran(&idum);
    long double tmp2=ran(&idum);
    long double tmp3=ran(&idum);
    long double tmp4=ran(&idum);
    long double tmp5=ran(&idum);
    long double tmp6=ran(&idum);
    long double tmp7=ran(&idum);  // tmpx is between (0,1), (2*tmpx-1) is between (-1,1)

    a=r;
    b=r;
    c=r;

    // generate orientation of axle a/b/c
    long double l1,m1,n1,l2,m2,n2,x,y,z;
    n1 = 2*tmp1-1;                  // (-1,1)
    z  = acosl(n1);                 
    m1 = (2*tmp2-1)*sqrtl(1-n1*n1); // (-1,1)
    y  = acosl(m1);
    int sign = 2*tmp3-1>0?1:-1;
    l1 = sign*sqrtl(1-m1*m1-n1*n1); // (-1,1)
    x  = acosl(l1);
    curr_direction_a=vec(x,y,z);    // axle a,b,c are in the direction of x,y,z respectively
    
    n2 = 2*tmp4 - 1;                // (-1,1)
    m2 = 2*tmp5 - 1;                // (-1,1)
    l2 = -(m1*m2+n1*n2)/l1;         // ensure a/b/c perpendicular
    long double mod2 = sqrtl(l2*l2+m2*m2+n2*n2);
    n2 /= mod2;
    m2 /= mod2;
    l2 /= mod2;
    z  = acosl(n2);
    y  = acosl(m2);
    x  = acosl(l2);
    curr_direction_b=vec(x,y,z);

    curr_direction_c=vacosl(normalize(vcosl(curr_direction_a)*vcosl(curr_direction_b)));

    curr_position=prev_position=center;
    prev_direction_a=curr_direction_a;
    prev_direction_b=curr_direction_b;
    prev_direction_c=curr_direction_c;
    curr_velocity=prev_velocity=0;
    curr_omga=prev_omga=0;
    curr_acceleration=prev_acceleration=0;
    force=pre_force=0;moment=pre_moment=0;mres=0;
    const_force=const_moment=0;
    density=Gs*1.0e3;
    volume=4/3.0*PI*a*b*c;
    mass=density*volume;
    J=vec(mass/5*(b*b+c*c),mass/5*(a*a+c*c),mass/5*(a*a+b*b));
    cntnum=0;
    GlobCoef();
}


particle::particle(int n, int tp, vec center, long double _a, long double _b, long double _c){
    type=tp;
    ID=n;
    long double tmp0=ran(&idum);
    long double tmp1=ran(&idum);
    long double tmp2=ran(&idum);
    long double tmp3=ran(&idum);
    long double tmp4=ran(&idum);
    long double tmp5=ran(&idum);
    long double tmp6=ran(&idum);
    long double tmp7=ran(&idum);  // tmpx is between (0,1), (2*tmpx-1) is between (-1,1)

    a=_a;
    b=_b;
    c=_c;

    // generate orientation of axle a/b/c
    long double l1,m1,n1,l2,m2,n2,x,y,z;
    n1 = 2*tmp1-1;                  // (-1,1)
    z  = acosl(n1);                 
    m1 = (2*tmp2-1)*sqrtl(1-n1*n1); // (-1,1)
    y  = acosl(m1);
    int sign = 2*tmp3-1>0?1:-1;
    l1 = sign*sqrtl(1-m1*m1-n1*n1); // (-1,1)
    x  = acosl(l1);
    curr_direction_a=vec(x,y,z);    // axle a,b,c are in the direction of x,y,z respectively
    
    n2 = 2*tmp4 - 1;                // (-1,1)
    m2 = 2*tmp5 - 1;                // (-1,1)
    l2 = -(m1*m2+n1*n2)/l1;         // ensure a/b/c perpendicular
    long double mod2 = sqrtl(l2*l2+m2*m2+n2*n2);
    n2 /= mod2;
    m2 /= mod2;
    l2 /= mod2;
    z  = acosl(n2);
    y  = acosl(m2);
    x  = acosl(l2);
    curr_direction_b=vec(x,y,z);

    curr_direction_c=vacosl(normalize(vcosl(curr_direction_a)*vcosl(curr_direction_b)));

    curr_position=prev_position=center;
    prev_direction_a=curr_direction_a;
    prev_direction_b=curr_direction_b;
    prev_direction_c=curr_direction_c;
    curr_velocity=prev_velocity=0;
    curr_omga=prev_omga=0;
    curr_acceleration=prev_acceleration=0;
    force=pre_force=0;moment=pre_moment=0;mres=0;
    const_force=const_moment=0;
    density=Gs*1.0e3;
    volume=4/3.0*PI*a*b*c;
    mass=density*volume;
    J=vec(mass/5*(b*b+c*c),mass/5*(a*a+c*c),mass/5*(a*a+b*b));
    cntnum=0;
    GlobCoef();
}


particle::particle(int n, int tp, vec center, gradation& grad){
    type=tp;
    ID=n;
    long double tmp0=ran(&idum);
    long double tmp1=ran(&idum);
    long double tmp2=ran(&idum);
    long double tmp3=ran(&idum);
    long double tmp4=ran(&idum);
    long double tmp5=ran(&idum);
    long double tmp6=ran(&idum);
    long double tmp7=ran(&idum);
    // tmpx is between (0,1), (2*tmpx-1) is between (-1,1)

    // generate a particle based on gradation distribution
    for (int k=0;k<grad.seivenum;k++){
	if (tmp0 <= grad.percent[grad.seivenum-1-k]){
	    a=grad.ptclsize[grad.seivenum-1-k];
	    break;
	}
    }

#ifdef RANDOM_SHAPE
    grad.ratio_ba=tmp6;
    grad.ratio_ca=tmp7;
#endif

    b=a*grad.ratio_ba;
    c=a*grad.ratio_ca;

    // generate orientation of axle a/b/c
    long double l1,m1,n1,l2,m2,n2,x,y,z;
    n1 = 2*tmp1-1;                  // (-1,1)
    z  = acosl(n1);                 
    m1 = (2*tmp2-1)*sqrtl(1-n1*n1); // (-1,1)
    y  = acosl(m1);
    int sign = 2*tmp3-1>0?1:-1;
    l1 = sign*sqrtl(1-m1*m1-n1*n1); // (-1,1)
    x  = acosl(l1);
    curr_direction_a=vec(x,y,z);    // axle a,b,c are in the direction of x,y,z respectively
    
    n2 = 2*tmp4 - 1;                // (-1,1)
    m2 = 2*tmp5 - 1;                // (-1,1)
    l2 = -(m1*m2+n1*n2)/l1;         // ensure a/b/c perpendicular
    long double mod2 = sqrtl(l2*l2+m2*m2+n2*n2);
    n2 /= mod2;
    m2 /= mod2;
    l2 /= mod2;
    z  = acosl(n2);
    y  = acosl(m2);
    x  = acosl(l2);
    curr_direction_b=vec(x,y,z);

    curr_direction_c=vacosl(normalize(vcosl(curr_direction_a)*vcosl(curr_direction_b)));
    
    curr_position=prev_position=center;
    prev_direction_a=curr_direction_a;
    prev_direction_b=curr_direction_b;
    prev_direction_c=curr_direction_c;
    curr_velocity=prev_velocity=0;
    curr_omga=prev_omga=0;
    curr_acceleration=prev_acceleration=0;
    force=pre_force=0;moment=pre_moment=0;mres=0;
    const_force=const_moment=0;
    density=Gs*1.0e3;
    volume=4/3.0*PI*a*b*c;
    mass=density*volume;
    J=vec(mass/5*(b*b+c*c),mass/5*(a*a+c*c),mass/5*(a*a+b*b));
    cntnum=0;
    GlobCoef();
}


particle::particle(int id, int tp, vec dim, vec position, vec dirca, vec dircb, vec dircc){
    type=tp;
    ID=id;
    a=dim.getx();
    b=dim.gety();
    c=dim.getz();
    curr_position=prev_position=position;
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
    volume=4/3.0*PI*a*b*c;
    mass=density*volume;
    J=vec(mass/5*(b*b+c*c),mass/5*(a*a+c*c),mass/5*(a*a+b*b));
    GlobCoef();
}


int    particle::getID() const {return ID;}
int    particle::getType() const {return type;}
long double particle::getA() const {return a;}
long double particle::getB() const {return b;}
long double particle::getC() const {return c;}
long double particle::getVolume() const {return volume;}
long double particle::getMass() const {return mass;}
long double particle::getDensity() const {return density;}
vec    particle::getCurrPosition() const {return curr_position;}
vec    particle::getPrevPosition() const {return prev_position;}
vec    particle::getCurrDirecA() const {return curr_direction_a;}
vec    particle::getCurrDirecB() const {return curr_direction_b;}
vec    particle::getCurrDirecC() const {return curr_direction_c;}
vec    particle::getPrevDirecA() const {return prev_direction_a;}
vec    particle::getPrevDirecB() const {return prev_direction_b;}
vec    particle::getPrevDirecC() const {return prev_direction_c;}
vec    particle::getCurrVelocity() const {return curr_velocity;}
vec    particle::getPrevVelocity() const {return prev_velocity;}
vec    particle::getCurrOmga() const {return curr_omga;}
vec    particle::getPrevOmga() const {return prev_omga;}
vec    particle::getCurrAcceleration() const {return curr_acceleration;}
vec    particle::getPrevAcceleration() const {return prev_acceleration;}
vec    particle::getForce() const {return force;}
vec    particle::getMoment() const {return moment;}
vec    particle::getConstForce() const {return const_force;}
vec    particle::getConstMoment() const {return const_moment;}
vec    particle::getJ() const {return J;}


// 1: rotational energy is 1/2(I1*w1^2+I2*w2^2+I3*w3^2), where each term
//    is expressed in local frame.
// 2. angular velocities in global frame needs to be converted to those
//    in local frame.
long double particle::getTransEnergy() const{
    return mass*powl(vfabsl(curr_velocity),2)/2;
}


long double particle::getRotaEnergy() const{
    vec curr_local_omga, tmp;

    tmp=vcosl(curr_direction_a); curr_local_omga.setx(tmp%curr_omga);
    tmp=vcosl(curr_direction_b); curr_local_omga.sety(tmp%curr_omga);
    tmp=vcosl(curr_direction_c); curr_local_omga.setz(tmp%curr_omga);

    return J.getx()*powl(curr_local_omga.getx(),2)/2 +
	   J.gety()*powl(curr_local_omga.gety(),2)/2 +
	   J.getz()*powl(curr_local_omga.getz(),2)/2;
}


long double particle::getKinEnergy() const{
    return getTransEnergy() + getRotaEnergy();
}


long double particle::getPotEnergy(long double ref) const{
    return 9.8*mass*(curr_position.getz() - ref);
}


void   particle::setID(int n){ID=n;}
void   particle::setA(long double dd){a=dd;}
void   particle::setB(long double dd){b=dd;}
void   particle::setC(long double dd){c=dd;}
void   particle::setCurrPosition(vec vv){curr_position=vv;}
void   particle::setPrevPosition(vec vv){prev_position=vv;}
void   particle::setCurrDirecA(vec vv){curr_direction_a=vv;}
void   particle::setCurrDirecB(vec vv){curr_direction_b=vv;}
void   particle::setCurrDirecC(vec vv){curr_direction_c=vv;}
void   particle::setPrevDirecA(vec vv){prev_direction_a=vv;}
void   particle::setPrevDirecB(vec vv){prev_direction_b=vv;}
void   particle::setPrevDirecC(vec vv){prev_direction_c=vv;}
void   particle::setCurrVelocity(vec vv){curr_velocity=vv;}
void   particle::setPrevVelocity(vec vv){prev_velocity=vv;}
void   particle::setCurrOmga(vec vv){curr_omga=vv;}
void   particle::setPrevOmga(vec vv){prev_omga=vv;}
void   particle::setCurrAcceleration(vec vv){curr_acceleration=vv;}
void   particle::setPrevAcceleration(vec vv){prev_acceleration=vv;}
void   particle::setForce(vec vv){force=vv;}
void   particle::setMoment(vec vv){moment=vv;}
void   particle::setConstForce(vec vv){const_force=vv;}
void   particle::setConstMoment(vec vv){const_moment=vv;}
void   particle::setJ(vec v){J=v;}
void   particle::setMass(long double d){mass=d;}
void   particle::getGlobCoef(long double coef[]) const{
    for (int i=0;i<10;i++)
	coef[i]=this->coef[i];
}


void particle::print() const{
    cout<<"a="<<a<<'\t'<<"b="<<b<<'\t'<<"c="<<c<<endl;
    cout<<"curr_direction_a=";
    vcosl(curr_direction_a).print();
    cout<<"curr_direction_b=";
    vcosl(curr_direction_b).print();
    cout<<"curr_direction_c=";
    vcosl(curr_direction_c).print();
    cout<<"curr_position=";
    curr_position.print();
}


long double particle::surfaceError(vec pt) const{
    long double x=pt.getx();
    long double y=pt.gety();
    long double z=pt.getz();
    return coef[0]*x*x+coef[1]*y*y+coef[2]*z*z+coef[3]*x*y+coef[4]*y*z+coef[5]*z*x+
	coef[6]*x+coef[7]*y+coef[8]*z+coef[9];
}


void particle::GlobCoef(){
    //coef[0]--x^2, coef[1]--y^2, coef[2]--z^2, coef[3]--xy, coef[4]--yz, coef[5]--xz
    //coef[6]--x, coef[7]--y, coef[8]--z, coef[9]--const
    if(a==b&&b==c){
	coef[0]=1;
	coef[1]=1;
	coef[2]=1;
	coef[3]=0;
	coef[4]=0;
	coef[5]=0;
	coef[6]=-2*curr_position.getx();
	coef[7]=-2*curr_position.gety();
	coef[8]=-2*curr_position.getz();
	coef[9]=powl(vfabsl(curr_position),2)-a*a;
	return;
    }
    vec v1=vcosl(curr_direction_a);
    vec v2=vcosl(curr_direction_b);
    vec v3=vcosl(curr_direction_c);
    long double X0=curr_position.getx();
    long double Y0=curr_position.gety();
    long double Z0=curr_position.getz();
    long double l1=v1.getx();
    long double m1=v1.gety();
    long double n1=v1.getz();
    long double l2=v2.getx();
    long double m2=v2.gety();
    long double n2=v2.getz();
    long double l3=v3.getx();
    long double m3=v3.gety();
    long double n3=v3.getz();
    coef[0]=l1*l1/a/a+l2*l2/b/b+l3*l3/c/c;
    coef[1]=m1*m1/a/a+m2*m2/b/b+m3*m3/c/c;
    coef[2]=n1*n1/a/a+n2*n2/b/b+n3*n3/c/c;
    coef[3]=(2*l1*m1)/a/a + (2*l2*m2)/b/b + (2*l3*m3)/c/c;
    coef[4]=(2*m1*n1)/a/a + (2*m2*n2)/b/b + (2*m3*n3)/c/c;
    coef[5]=(2*l1*n1)/a/a + (2*l2*n2)/b/b + (2*l3*n3)/c/c;
    coef[6]=
	-2*l1*m1*Y0*powl(a,-2) - 2*l1*n1*Z0*powl(a,-2) - 
	2*l2*m2*Y0*powl(b,-2) - 2*l2*n2*Z0*powl(b,-2) - 
	2*l3*m3*Y0*powl(c,-2) - 2*l3*n3*Z0*powl(c,-2) - 
	2*X0*powl(a,-2)*powl(l1,2) - 2*X0*powl(b,-2)*powl(l2,2) - 
	2*X0*powl(c,-2)*powl(l3,2);
    coef[7]=
	(-2*l1*m1*X0)/a/a - (2*l2*m2*X0)/b/b - 
	(2*l3*m3*X0)/c/c - (2*m1*m1*Y0)/a/a - 
	(2*m2*m2*Y0)/b/b - (2*m3*m3*Y0)/c/c - 
	(2*m1*n1*Z0)/a/a - (2*m2*n2*Z0)/b/b - 
	(2*m3*n3*Z0)/c/c;
    coef[8]=
	(-2*l1*n1*X0)/a/a - (2*l2*n2*X0)/b/b - 
	(2*l3*n3*X0)/c/c - (2*m1*n1*Y0)/a/a - 
	(2*m2*n2*Y0)/b/b - (2*m3*n3*Y0)/c/c - 
	(2*n1*n1*Z0)/a/a - (2*n2*n2*Z0)/b/b - 
	(2*n3*n3*Z0)/c/c;
    coef[9]=
	-1 + 2*l1*m1*X0*Y0*powl(a,-2) + 2*l1*n1*X0*Z0*powl(a,-2) + 
	2*m1*n1*Y0*Z0*powl(a,-2) + 2*l2*m2*X0*Y0*powl(b,-2) + 
	2*l2*n2*X0*Z0*powl(b,-2) + 2*m2*n2*Y0*Z0*powl(b,-2) + 
	2*l3*m3*X0*Y0*powl(c,-2) + 2*l3*n3*X0*Z0*powl(c,-2) + 
	2*m3*n3*Y0*Z0*powl(c,-2) + 
	powl(a,-2)*powl(l1,2)*powl(X0,2) + 
	powl(b,-2)*powl(l2,2)*powl(X0,2) + 
	powl(c,-2)*powl(l3,2)*powl(X0,2) + 
	powl(a,-2)*powl(m1,2)*powl(Y0,2) + 
	powl(b,-2)*powl(m2,2)*powl(Y0,2) + 
	powl(c,-2)*powl(m3,2)*powl(Y0,2) + 
	powl(a,-2)*powl(n1,2)*powl(Z0,2) + 
	powl(b,-2)*powl(n2,2)*powl(Z0,2) + 
	powl(c,-2)*powl(n3,2)*powl(Z0,2);
    long double divd=coef[0];
    for (int kk=0;kk<10;kk++){  // when a particle is initialized or updated, coef[0] is set as 1.0.
	coef[kk]/=divd;
    }
}


bool particle::intersectWithLine(vec v, vec dirc, vec rt[]) const{
    long double x0=v.getx();
    long double y0=v.gety();
    long double z0=v.getz();
    long double p=dirc.getx();
    long double q=dirc.gety();
    long double r=dirc.getz();
    long double a=coef[0];
    long double b=coef[1];
    long double c=coef[2];
    long double d=coef[3];
    long double e=coef[4];
    long double f=coef[5];
    long double g=coef[6];
    long double h=coef[7];
    long double i=coef[8];
    long double j=coef[9];

    long double A = a*p*p + b*q*q + c*r*r + d*p*q + e*q*r + f*r*p;
    long double B = 2*a*p*x0 + 2*b*q*y0 + 2*c*r*z0
	            + d*p*y0 + d*q*x0 + e*q*z0 + e*r*y0 + f*p*z0 + f*r*x0
	            + g*p + h*q + i*r;
    long double C = a*x0*x0 + b*y0*y0 + c*z0*z0 + d*x0*y0 + e*y0*z0 +f*z0*x0
	            + g*x0 + h*y0 + i*z0 + j;
    
    long double delta=B*B-4*A*C;
    if (delta < 0){
	g_exceptioninf<<"g_iteration="<<setw(10)<<g_iteration
		    <<", delta < 0 in intersectWithLine() of particle.cxx."<<endl;
	return false;
    }
    else{
	long double t1=(-B+sqrtl(delta))/(2*A);
	long double t2=(-B-sqrtl(delta))/(2*A);
	
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
//    2. r = r1*r2/(r1+r2)
//    3. It is essential to eliminate float exceptions in computations, that is, 
//       when dz/dx == infinite, coordinate x & z are switched to use dx/dz == 0.
//    4. When a point is close to the equator, for example, fabsl(z)==0,
//       float exception is prone to occurring, then a switch is needed
//       as above.
long double particle::getRadius(vec v) const{
    if(a==b&&b==c)
	return a;

    long double per=1.0e-4; // define when a point is close to equator
    long double ra=a;       // semi-axles of ellipsoid
    long double rb=b;
    long double rc=c;

    // get the local coodinates of vector v, the point on the particle's surface
    vec v1=vcosl(curr_direction_a);
    vec v2=vcosl(curr_direction_b);
    vec v3=vcosl(curr_direction_c);
    long double X0=curr_position.getx();
    long double Y0=curr_position.gety();
    long double Z0=curr_position.getz();
    long double x1=v.getx()-X0;
    long double y1=v.gety()-Y0;
    long double z1=v.getz()-Z0;
    long double l1=v1.getx();
    long double m1=v1.gety();
    long double n1=v1.getz();
    long double l2=v2.getx();
    long double m2=v2.gety();
    long double n2=v2.getz();
    long double l3=v3.getx();
    long double m3=v3.gety();
    long double n3=v3.getz();
    long double x=l1*x1 + m1*y1 + n1*z1;
    long double y=l2*x1 + m2*y1 + n2*z1;
    long double z=l3*x1 + m3*y1 + n3*z1;

    long double tmp;
    if (fabsl(z)<=c*per) {     // switch x & z, use 0 instead of infinity
	tmp=ra; ra=rc; rc=tmp;
	tmp=x; x=z; z=tmp; 
	if (fabsl(z)<=a*per) { // switch y & z, use 0 instead of infinity
	    tmp=ra; ra=rb; rb=tmp;
	    tmp=y; y=z; z=tmp; 
	}     
    }

    long double p=-rc*rc/ra/ra*x/z;
    long double q=-rc*rc/rb/rb*y/z;
    long double r=-rc*rc/ra/ra*(1/z+rc*rc/ra/ra*x*x/powl(z,3));
    long double t=-rc*rc/rb/rb*(1/z+rc*rc/rb/rb*y*y/powl(z,3));
    long double s=-powl(rc,4)/ra/ra/rb/rb*x*y/powl(z,3);
    long double n  = sqrtl(1+p*p+q*q);

    long double A,B,C;
    A=r*t-s*s;
    B=n*(2*p*q*s-(1+p*p)*t-(1+q*q)*r);
    C=n*n*n*n;


    // if delta < 0, then it is usually -1.0e-20, caused by computational precision.
/*
    if (B*B-4*A*C<0){
	g_exceptioninf<<"g_iteration="<<setw(10)<<g_iteration
		      <<", delta < 0 in getRadius() of particle.cxx."
		      <<setw(16)<<B*B-4*A*C
		      <<setw(16)<<-C/B
		      <<endl;
    }
*/
    return fabsl(-C/B); // r1*r2/(r1+r2)
}


void particle::setZero(){
    force=const_force;
    moment=const_moment;
    force += vec(0,0,-9.8*mass*GRVT_SCL); // Unit is Newton, GRVT_SCL is for amplification.
    if (getType()==3) // pile
	force -= vec(0,0,-9.8*mass*GRVT_SCL); 

#ifdef MOMENT
	long double m[20]={ 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
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


void particle::setZero(int PrintNum){
    force  = const_force * (PrintNum+1)/100;
    moment = const_moment * (PrintNum+1)/100;
    force += vec(0,0,-9.8*mass*GRVT_SCL); // Unit is Newton, GRVT_SCL is for amplification.
    if (getType()==3) // pile
	force -= vec(0,0,-9.8*mass*GRVT_SCL); 
}


// central difference integration method
void particle::update() {

    if (getType()==0) { // 0-free, 1-fixed
	// It is essential to distinguish global frame from local frame!
	vec prev_local_omga;
	vec curr_local_omga;
	vec local_moment;
	vec tmp;
	long double atf=DMP_F*TIMESTEP; 
	long double atm=DMP_M*TIMESTEP; 
	
	// force: translational kinetics equations are in global frame
	curr_velocity=prev_velocity*(2-atf)/(2+atf)+force/(mass*MASS_SCL)*TIMESTEP*2/(2+atf);
	curr_position = prev_position + curr_velocity*TIMESTEP;

	// moment: angular kinetics (rotational) equations are in local frame,
	// so global values need to be converted to those in local frame when applying equations
	tmp=vcosl(getCurrDirecA()); local_moment.setx(tmp%moment); prev_local_omga.setx(tmp%prev_omga); // l1,m1,n1
	tmp=vcosl(getCurrDirecB()); local_moment.sety(tmp%moment); prev_local_omga.sety(tmp%prev_omga); // l2,m2,n2
	tmp=vcosl(getCurrDirecC()); local_moment.setz(tmp%moment); prev_local_omga.setz(tmp%prev_omga); // l3,m3,n3
	
	curr_local_omga.setx( prev_local_omga.getx()*(2-atm)/(2+atm) + local_moment.getx()/(J.getx()*MNT_SCL)*TIMESTEP*2/(2+atm) ); 
	curr_local_omga.sety( prev_local_omga.gety()*(2-atm)/(2+atm) + local_moment.gety()/(J.gety()*MNT_SCL)*TIMESTEP*2/(2+atm) );
	curr_local_omga.setz( prev_local_omga.getz()*(2-atm)/(2+atm) + local_moment.getz()/(J.getz()*MNT_SCL)*TIMESTEP*2/(2+atm) );
	
	// convert local angular velocities to those in global frame in order to rotate a particle in global space
	tmp=vcosl( vec(curr_direction_a.getx(),curr_direction_b.getx(),curr_direction_c.getx()) ); // l1,l2,l3
	curr_omga.setx(tmp%curr_local_omga);
	
	tmp=vcosl( vec(curr_direction_a.gety(),curr_direction_b.gety(),curr_direction_c.gety()) ); // m1,m2,m3
	curr_omga.sety(tmp%curr_local_omga);   
	
	tmp=vcosl( vec(curr_direction_a.getz(),curr_direction_b.getz(),curr_direction_c.getz()) ); // n1,n2,n3
	curr_omga.setz(tmp%curr_local_omga);
	
	curr_direction_a=vacosl(normalize(rotateVec(vcosl(prev_direction_a),curr_omga*TIMESTEP)));
	curr_direction_b=vacosl(normalize(rotateVec(vcosl(prev_direction_b),curr_omga*TIMESTEP)));
	curr_direction_c=vacosl(normalize(rotateVec(vcosl(prev_direction_c),curr_omga*TIMESTEP)));
    }
#ifdef MOMENT
    else if (getType()==2) { //special case 2 (moment): translate first, then rotate
	vec prev_local_omga;
	vec curr_local_omga;
	vec local_moment;
	vec tmp;
	long double atf=DMP_F*TIMESTEP; 
	long double atm=DMP_M*TIMESTEP; 
	curr_velocity=prev_velocity*(2-atf)/(2+atf)+force/(mass*MASS_SCL)*TIMESTEP*2/(2+atf);
	if (g_iteration < START)
	    curr_position = prev_position + curr_velocity*TIMESTEP;	

	tmp=vcosl(getCurrDirecA()); local_moment.setx(tmp%moment); prev_local_omga.setx(tmp%prev_omga); // l1,m1,n1
	tmp=vcosl(getCurrDirecB()); local_moment.sety(tmp%moment); prev_local_omga.sety(tmp%prev_omga); // l2,m2,n2
	tmp=vcosl(getCurrDirecC()); local_moment.setz(tmp%moment); prev_local_omga.setz(tmp%prev_omga); // l3,m3,n3
	
	curr_local_omga.setx( prev_local_omga.getx()*(2-atm)/(2+atm) + local_moment.getx()/(J.getx()*MNT_SCL)*TIMESTEP*2/(2+atm) ); 
	curr_local_omga.sety( prev_local_omga.gety()*(2-atm)/(2+atm) + local_moment.gety()/(J.gety()*MNT_SCL)*TIMESTEP*2/(2+atm) );
	curr_local_omga.setz( prev_local_omga.getz()*(2-atm)/(2+atm) + local_moment.getz()/(J.getz()*MNT_SCL)*TIMESTEP*2/(2+atm) );

	if (g_iteration >= START) {	
	    tmp=vcosl( vec(curr_direction_a.getx(),curr_direction_b.getx(),curr_direction_c.getx()) ); // l1,l2,l3
	    curr_omga.setx(tmp%curr_local_omga);
	    
	    tmp=vcosl( vec(curr_direction_a.gety(),curr_direction_b.gety(),curr_direction_c.gety()) ); // m1,m2,m3
	    curr_omga.sety(tmp%curr_local_omga);   
	    
	    tmp=vcosl( vec(curr_direction_a.getz(),curr_direction_b.getz(),curr_direction_c.getz()) ); // n1,n2,n3
	    curr_omga.setz(tmp%curr_local_omga);
	    
	    curr_direction_a=vacosl(normalize(rotateVec(vcosl(prev_direction_a),curr_omga*TIMESTEP)));
	    curr_direction_b=vacosl(normalize(rotateVec(vcosl(prev_direction_b),curr_omga*TIMESTEP)));
	    curr_direction_c=vacosl(normalize(rotateVec(vcosl(prev_direction_c),curr_omga*TIMESTEP)));
	}
    }
#endif
    else if (getType()==3) { //special case 3 (displacemental ellipsoidal pile): translate in vertical direction only
	curr_velocity.setx(0);	
	curr_velocity.sety(0);
	curr_velocity.setz(-PILE_RATE);
	curr_position = prev_position + curr_velocity*TIMESTEP;
    }
    else if (getType()==4) { //special case 4 (impacting ellipsoidal penetrator): impact with inital velocity in vertical direction only 
	long double atf=DMP_F*TIMESTEP; 
	curr_velocity=prev_velocity*(2-atf)/(2+atf)+force/(mass*MASS_SCL)*TIMESTEP*2/(2+atf);
	curr_velocity.setx(0);	
	curr_velocity.sety(0);
	curr_position = prev_position + curr_velocity*TIMESTEP;
    }

    // Below is needed for all cases
    // ensure three axles perpendicular to each other, and being unit vector
    if(curr_direction_a==0)
	curr_direction_a=vacosl(normalize(vcosl(curr_direction_b)*vcosl(curr_direction_c)));
    if(curr_direction_b==0)
	curr_direction_b=vacosl(normalize(vcosl(curr_direction_c)*vcosl(curr_direction_a)));
    if(curr_direction_c==0)
	curr_direction_c=vacosl(normalize(vcosl(curr_direction_a)*vcosl(curr_direction_b)));

    prev_position=curr_position;
    prev_direction_a=curr_direction_a;
    prev_direction_b=curr_direction_b;
    prev_direction_c=curr_direction_c;
    prev_velocity=curr_velocity;
    prev_omga=curr_omga;
    pre_force=force; 
    pre_moment=moment;

    cntnum=0;
    GlobCoef();   // every time the particle is updated, the algebra expression is also updated
}


vec particle::localVec(vec v) const{
    // v is a vector in global coordinates, it is to be transformed into local coordinates
    vec l=vcosl(vec(curr_direction_a.getx(),curr_direction_b.getx(),curr_direction_c.getx()));
    vec m=vcosl(vec(curr_direction_a.gety(),curr_direction_b.gety(),curr_direction_c.gety()));
    vec n=vcosl(vec(curr_direction_a.getz(),curr_direction_b.getz(),curr_direction_c.getz()));
    return l*v.getx()+m*v.gety()+n*v.getz();
}


vec particle::globalVec(vec v) const{
    vec l=vcosl(curr_direction_a);
    vec m=vcosl(curr_direction_b);
    vec n=vcosl(curr_direction_c);
    return l*v.getx()+m*v.gety()+n*v.getz();
}


bool particle::nearestPTOnPlane(long double p, long double q, long double r, long double s, vec& ptnp) const {
    if(a==b&&b==c){
      vec tnm=vec(p,q,r)/sqrtl(p*p+q*q+r*r);
      // signed distance from particle center to plane
      long double l_nm=(curr_position.getx()*p+curr_position.gety()*q+curr_position.getz()*r+s)/sqrtl(p*p+q*q+r*r); 
      ptnp=curr_position-l_nm*tnm;
      if(fabsl(l_nm)<a) // intersect
	return true;
      else              // no intersect, including tangent
	return false;
    }

    long double a=coef[0];
    long double b=coef[1];
    long double c=coef[2];
    long double d=coef[3];
    long double e=coef[4];
    long double f=coef[5];
    long double g=coef[6];
    long double h=coef[7];
    long double i=coef[8];
    long double j=coef[9];
    long double domi=
	 e*e*p*p + 4*c*d*p*q - 4*a*c*q*q + 
	 f*f*q*q - 2*d*f*q*r + 
	 d*d*r*r - 
	 2*e*(f*p*q + d*p*r - 2*a*q*r) - 
	 4*b*(c*p*p + r*(-f*p + a*r));
    long double x=
	(-(f*i*q*q) - 2*b*i*p*r + f*h*q*r + d*i*q*r + 
	 2*b*g*r*r - d*h*r*r - e*e*p*s - 
	 2*b*f*r*s - 2*c*(h*p*q - g*q*q - 2*b*p*s + 
	 d*q*s) + e*(i*p*q + h*p*r - 2*g*q*r + f*q*s + 
	 d*r*s))/domi;
    long double y=
	(f*i*p*q - 2*f*h*p*r + d*i*p*r + f*g*q*r - 
	 2*a*i*q*r - d*g*r*r + 2*a*h*r*r - 
	 f*f*q*s + d*f*r*s + 
	 2*c*(h*p*p - g*p*q - d*p*s + 2*a*q*s) + 
	 e*(-i*p*p + g*p*r + f*p*s - 2*a*r*s))/domi;
    long double z=
	(f*h*p*q - 2*d*i*p*q - f*g*q*q + 
	 2*a*i*q*q + d*h*p*r + d*g*q*r - 2*a*h*q*r + 
	 d*f*q*s - d*d*r*s + 
	 e*(-h*p*p + g*p*q + d*p*s - 2*a*q*s) + 
	 2*b*(i*p*p - g*p*r - f*p*s + 2*a*r*s))/domi;
    ptnp=vec(x,y,z);

    long double val=a*x*x+b*y*y+c*z*z+d*x*y+e*y*z+f*x*z+g*x+h*y+i*z+j;

    if (val >= 0) // not intersect
	return false;
    else  // intersect
	return true;
}


void particle::planeRBForce(plnrgd_bdry<particle>* plb,
			    map<int,vector<boundarytgt> >& BdryTgtMap,
			    vector<boundarytgt>& vtmp,
			    long double &penetr){
	// (p,q,r) are in the same direction as the outward normal vector,
	// hence it is not necessary to provide information about which side the particle is about the plane.
	long double p,q,r,s;
	int  bdry_id;
	BdryCoef tmp1=*((plb->CoefOfLimits).begin());
	p=tmp1.dirc.getx();
	q=tmp1.dirc.gety();
	r=tmp1.dirc.getz();
	s=-tmp1.dirc%tmp1.apt;  // plane equation: p(x-x0)+q(y-y0)+r(z-z0)=0, that is, px+qy+rz+s=0
	bdry_id=plb->bdry_id;

	vec pt1;
	if (!nearestPTOnPlane(p, q, r, s, pt1)) // the particle and the plane does not intersect
	  return;

	// if particle and plane intersect:
	cntnum++;
	vec dirc=normalize(vec(p,q,r));
	vec rt[2];
	if (!intersectWithLine(pt1,dirc,rt))    // the line and ellipsoid surface does not intersect
	    return;

	vec pt2;
/*
	if (p*rt[0].getx()+q*rt[0].gety()+r*rt[0].getz()+s > 0)
		pt2=rt[0];
	else
		pt2=rt[1];
*/
	if (vfabsl(rt[0]-pt1) < vfabsl(rt[1]-pt1) )
	    pt2 = rt[0];
	else
	    pt2 = rt[1];

	// obtain normal force
	long double penetration=vfabsl(pt1-pt2);
	if (penetration==0)  // plane is tangent to ellipsoid (this case can not be detected by nearestPTOnPlane because of computational precision)
	    return;
	penetr = penetration;
	long double R0=getRadius(pt2);
	long double contact_radius=sqrtl(penetration*R0);
	long double E0=YOUNG/(1-POISSON*POISSON); // rigid wall has infinite YOUNG's modulus
	vec NormDirc=-dirc; //normalize(pt1-pt2);
	vec NormalForce=sqrtl(penetration*penetration*penetration)*sqrtl(R0)*4*E0/3*NormDirc; // powl(penetration,1.5), a serious bug

/*
	g_exceptioninf<<setw(10)<<g_iteration
		      <<setw(10)<<getID()
		      <<setw(10)<<plb->bdry_id
		      <<setw(16)<<pt1.getx()
		      <<setw(16)<<pt1.gety()
		      <<setw(16)<<pt1.getz()
		      <<setw(16)<<rt[0].getx()
		      <<setw(16)<<rt[0].gety()
		      <<setw(16)<<rt[0].getz()
		      <<setw(16)<<rt[1].getx()
		      <<setw(16)<<rt[1].gety()
		      <<setw(16)<<rt[1].getz()
		      <<setw(16)<<vfabsl(rt[0]-pt1)
		      <<setw(16)<<vfabsl(rt[1]-pt1)
		      <<setw(16)<<penetration
		      <<endl;
*/

	// apply normal force
	addForce(NormalForce);
	addMoment(((pt1+pt2)/2-curr_position)*NormalForce);
	
	// obtain normal damping force
	vec cp = (pt1+pt2)/2;        
	vec veloc2 = getCurrVelocity() + getCurrOmga()*(cp-getCurrPosition());
	long double kn = powl(6*vfabsl(NormalForce)*R0*powl(E0,2),1.0/3.0);
	long double DMP_CRTC = 2*sqrtl(getMass()*kn); // critical damping
	vec CntDampingForce  = DMP_CNT * DMP_CRTC * ((-veloc2)%NormDirc)*NormDirc;

	// apply normal damping force
	addForce(CntDampingForce) ;

	vec TgtForce = 0;
	if (BDRYFRIC != 0){
	    // checkin previous tangential force and displacement
	    vec PreTgtForce;
	    vec PreTgtDisp;
	    bool PreTgtLoading=true;
	    vec  TgtDispStart;
	    long double TgtPeak=0;

	    bool TgtLoading=true;
	    vector<boundarytgt>::iterator it;
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
	    long double G0 = YOUNG/2/(1+POISSON);
	    vec cp = (pt1+pt2)/2;
	    // Vr = Vb + w (crossdot) r, each item needs to be in either global or local frame; 
	    //      here global frame is used for better convenience.
	    vec RelaDispInc = (curr_velocity+curr_omga*(cp-curr_position))*TIMESTEP;  
	    vec TgtDispInc = RelaDispInc-(RelaDispInc%NormDirc)*NormDirc;
	    vec TgtDisp    = PreTgtDisp + TgtDispInc; // PreTgtDisp read by checkin
	    vec TgtDirc;

	    if (vfabsl(TgtDisp) == 0)
		TgtDirc = 0;
	    else
		TgtDirc = normalize(-TgtDisp); // TgtDirc points along tangential forces exerted on particle 1

	    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	    // linear friction model
	    long double fP  = BDRYFRIC*vfabsl(NormalForce);
	    long double ks  = 4*G0*contact_radius/(2-POISSON);
	    TgtForce = PreTgtForce + ks*(-TgtDispInc); // PreTgtForce read by checkin
	    if (vfabsl(TgtForce) > fP)
		TgtForce = fP*TgtDirc;
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////

	    // apply tangential force
	    addForce(TgtForce);
	    addMoment(((pt1+pt2)/2-curr_position)*TgtForce); 

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


vec particle::cylinderRBForce(int bdry_id, const cylinder& S, int side){
	//side=-1, the particles are inside the cylinder
	//side=+1, the particles are outside the cylinder
	long double x0=S.get_center().getx();
	long double y0=S.get_center().gety();
	long double r=S.get_radius();
	long double coef2[10]={1,1,0,0,0,0,-2*x0,-2*y0,0,x0*x0+y0*y0-r*r};
	vec pt1;
	if (!root6(coef,coef2,pt1))  //on the cylinder and within the particle
	    return 0; //no contact
	cntnum++;
	vec rt[2];
	vec cz=vec(S.get_center().getx(),S.get_center().gety(),pt1.getz());
	vec tmp=pt1-cz;
	intersectWithLine(pt1, normalize(tmp),rt);
	vec pt2;

	if ((rt[0]-pt1)%tmp*side<0)
		pt2=rt[0];
	else
		pt2=rt[1];
	//vec pt2=vfabsl(rt[0]-cz)>vfabsl(rt[1]-cz)?rt[0]:rt[1];
	long double radius=getRadius(pt2);//pt2.print();cout<<radius;getchar();
	long double E0=0.5*YOUNG/(1-POISSON*POISSON);
	long double R0=(r*radius)/(r+radius);
	long double rou=vfabsl(pt1-pt2);
	vec NormDirc=normalize(pt1-pt2);
	long double nfc=sqrtl(rou*rou*rou)*sqrtl(R0)*4*E0/3; // powl(rou,1.5), a serious bug
	vec NormalForce=nfc*NormDirc;

	addForce(NormalForce);
	addMoment(((pt1+pt2)/2-getCurrPosition())*NormalForce);	    
	
	return NormalForce;
}

} // namespace dem ends
