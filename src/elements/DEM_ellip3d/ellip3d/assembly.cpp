//                          ---------------------
//                         /                    /|
//                        /                    / |
//                       /                    /  |
//                      /         5          /   |
//                     /                    /    |height
//                    /                    /     |                    z (sigma3)
//                   /                    /      |                    |
//                  |---------------------       |                    |
//                  |                    |   2   |                    |____ y (sigma1)
//                  |                    |       /                   /
//                  |                    |      /                   /
//                  |         1          |     /                   x (sigma2) 
//                  |                    |    /length
//                  |                    |   /
//                  |                    |  /
//                  |                    | /
//                  |                    |/
//                  ----------------------
//                         width
//
//    sigma1_1 & sigma1_2 refers to side 2 & side 4 respectively,
//    sigma2_1 & sigma2_2 refers to side 1 & side 3 respectively,
//    sigma3_1 & sigma3_2 refers to side 5 & side 6 respectively,
//
//    int mid[2]={1,3};    // boundary 1 and 3
//    int max[2]={2,4};    // boundary 2 and 4
//    int min[2]={5,6};    // boundary 5 and 6
//    min/mid/max does not mean actual magnitude of values, just signs

#include "assembly.h"
#include "parameter.h"
#include "timefunc.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <ctime>

#ifdef OPENMP
#include <omp.h>
#endif

//#define TIME_PROFILE

using namespace std;

static time_t timeStamp; // for file timestamping
static struct timeval timew1, timew2; // for wall-clock time record
static struct timeval timep1, timep2; // for internal wall-clock time profiling

namespace dem {

ofstream        progressinf;

void assembly::printParticle(const char* str) const
{
    ofstream ofs(str);
    if(!ofs) {
	cout<<"stream error in printParticle!"<<endl; exit(-1);
    }
    ofs.setf(ios::scientific, ios::floatfield);
    ofs<<setw(10)<<TotalNum<<setw(10)<<RORC<<endl;
    if(RORC==0)
	ofs<<setw(16)<<S.get_center().getx()
	   <<setw(16)<<S.get_center().gety()
	   <<setw(16)<<S.get_center().getz()
	   <<setw(16)<<S.get_radius()
	   <<setw(16)<<S.get_height()<<endl;
    else
	ofs<<setw(16)<<R.get_center().getx()
	   <<setw(16)<<R.get_center().gety()
	   <<setw(16)<<R.get_center().getz()
	   <<setw(16)<<R.get_width()
	   <<setw(16)<<R.get_length()
	   <<setw(16)<<R.get_height()<<endl;
    ofs<<"     ID          type      radius_a        radius_b        radius_c   "
       <<"    position_x            position_y            position_z        "
       <<"    axle_a_x              axle_a_y              axle_a_z          "
       <<"    axle_b_x              axle_b_y              axle_b_z          "
       <<"    axle_c_x              axle_c_y              axle_c_z          "
       <<"   velocity_x      velocity_y      velocity_z   "
       <<"     omga_x          omga_y          omga_z     "
       <<"    force_x         force_y         force_z     "
       <<"    moment_x        moment_y        moment_z    "
       <<endl;

    vec tmp;
    list<particle*>::const_iterator  it;
    for (it=ParticleList.begin();it!=ParticleList.end();++it)
    {
	ofs<<setw(10)<<(*it)->getID()
	   <<setw(10)<<(*it)->getType()
	   <<setw(16)<<(*it)->getA()
	   <<setw(16)<<(*it)->getB()
	   <<setw(16)<<(*it)->getC();
	
	ofs.precision(12);
	tmp=(*it)->getCurrPosition();
	ofs<<setw(22)<<tmp.getx()
	   <<setw(22)<<tmp.gety()
	   <<setw(22)<<tmp.getz();
	
	tmp=(*it)->getCurrDirecA();
	ofs<<setw(22)<<tmp.getx()
	   <<setw(22)<<tmp.gety()
	   <<setw(22)<<tmp.getz();
	
	tmp=(*it)->getCurrDirecB();
	ofs<<setw(22)<<tmp.getx()
	   <<setw(22)<<tmp.gety()
	   <<setw(22)<<tmp.getz();
	
	tmp=(*it)->getCurrDirecC();
	ofs<<setw(22)<<tmp.getx()
	   <<setw(22)<<tmp.gety()
	   <<setw(22)<<tmp.getz();
	ofs.precision(6);
	
	tmp=(*it)->getCurrVelocity();
	ofs<<setw(16)<<tmp.getx()
	   <<setw(16)<<tmp.gety()
	   <<setw(16)<<tmp.getz();
	
	tmp=(*it)->getCurrOmga();
	ofs<<setw(16)<<tmp.getx()
	   <<setw(16)<<tmp.gety()
	   <<setw(16)<<tmp.getz();
	
	tmp=(*it)->getForce();
	ofs<<setw(16)<<tmp.getx()
	   <<setw(16)<<tmp.gety()
	   <<setw(16)<<tmp.getz();
	
	tmp=(*it)->getMoment();
	ofs<<setw(16)<<tmp.getx()
	   <<setw(16)<<tmp.gety()
	   <<setw(16)<<tmp.getz()<<endl;
    }

    ofs.close();
}


void assembly::printRectPile(const char* str)
{
    ofstream ofs(str, ios_base::app);
    if(!ofs) {
	cout<<"stream error in printRectPile!"<<endl; exit(-1);
    }
    ofs.setf(ios::scientific, ios::floatfield);

    ofs<<setw(10)<<8<<setw(10)<<6<<endl;
    vec pos[8];
    for(list<RGDBDRY*>::iterator rt=RBList.begin();rt!=RBList.end();++rt){
	if((*rt)->getBdryID()==7){
	    pos[0]=vec((*rt)->CoefOfLimits[0].apt.getx(),
		       (*rt)->CoefOfLimits[1].apt.gety(),
		       (*rt)->CoefOfLimits[4].apt.getz());
	    pos[1]=vec((*rt)->CoefOfLimits[0].apt.getx(),
		       (*rt)->CoefOfLimits[2].apt.gety(),
		       (*rt)->CoefOfLimits[4].apt.getz());
	    pos[5]=vec((*rt)->CoefOfLimits[0].apt.getx(),
		       (*rt)->CoefOfLimits[2].apt.gety(),
		       (*rt)->CoefOfLimits[3].apt.getz());
	    pos[4]=vec((*rt)->CoefOfLimits[0].apt.getx(),
		       (*rt)->CoefOfLimits[1].apt.gety(),
		       (*rt)->CoefOfLimits[3].apt.getz());
	}
	else if((*rt)->getBdryID()==9) {
	    pos[2]=vec((*rt)->CoefOfLimits[0].apt.getx(),
		       (*rt)->CoefOfLimits[2].apt.gety(),
		       (*rt)->CoefOfLimits[4].apt.getz());
	    pos[3]=vec((*rt)->CoefOfLimits[0].apt.getx(),
		       (*rt)->CoefOfLimits[1].apt.gety(),
		       (*rt)->CoefOfLimits[4].apt.getz());
	    pos[7]=vec((*rt)->CoefOfLimits[0].apt.getx(),
		       (*rt)->CoefOfLimits[1].apt.gety(),
		       (*rt)->CoefOfLimits[3].apt.getz());
	    pos[6]=vec((*rt)->CoefOfLimits[0].apt.getx(),
		       (*rt)->CoefOfLimits[2].apt.gety(),
		       (*rt)->CoefOfLimits[3].apt.getz());
	}
    }

    for (int i=0;i<8;i++)
	ofs<<setw(16)<<pos[i].getx()<<setw(16)<<pos[i].gety()<<setw(16)<<pos[i].getz()<<endl;

    ofs<<setw(10)<<1<<setw(10)<<2<<setw(10)<<6<<setw(10)<<5<<endl
       <<setw(10)<<2<<setw(10)<<3<<setw(10)<<7<<setw(10)<<6<<endl
       <<setw(10)<<3<<setw(10)<<4<<setw(10)<<8<<setw(10)<<7<<endl
       <<setw(10)<<4<<setw(10)<<1<<setw(10)<<5<<setw(10)<<8<<endl
       <<setw(10)<<1<<setw(10)<<4<<setw(10)<<3<<setw(10)<<2<<endl
       <<setw(10)<<5<<setw(10)<<6<<setw(10)<<7<<setw(10)<<8<<endl;

    ofs.close();
}


//  1. it is important and helpful to mark a member function as const
//     if it does NOT change member data.
//  2. when a constant member function traverses member data, it can
//     NOT change the data.
//  3. then if it traverses a member data of a list, it should use a
//     const_iterator, otherwise compiler will give errors.
//  4. a const_iterator such as it also guarantees that (*it) will NOT
//     change any data. if (*it) call a modification function, the 
//     compiler will give errors.
void assembly::printContact(const char* str) const
{
    ofstream ofs(str);
    if(!ofs) {
	cout<<"stream error in printContact!"<<endl; exit(-1);
    }
    ofs.setf(ios::scientific, ios::floatfield);
    ofs<<setw(10)<<ActualCntctNum<<endl;
    ofs<<"    ptcl_1     ptcl_2     point1_x        point1_y        point1_z        "
       <<"point2_x       point2_y        point2_z         radius_1        radius_2      "
       <<"penetration    tangt_dispmt     contact_radius       R0              E0         "
       <<"normal_force     tangt_force"
       <<"         x0              y0              z0     "
       <<"         n1              n2              n3     "
       <<"         t1              t2              t3"
       <<endl;
    list<CONTACT>::const_iterator it;
    for (it=ContactList.begin();it!=ContactList.end();++it)
	ofs<<setw(10)<<(*it).getP1()->getID()
	   <<setw(10)<<(*it).getP2()->getID()
	   <<setw(16)<<(*it).getPoint1().getx()
	   <<setw(16)<<(*it).getPoint1().gety()
	   <<setw(16)<<(*it).getPoint1().getz()
	   <<setw(16)<<(*it).getPoint2().getx()
	   <<setw(16)<<(*it).getPoint2().gety()
	   <<setw(16)<<(*it).getPoint2().getz()
	   <<setw(16)<<(*it).getRadius1()
	   <<setw(16)<<(*it).getRadius2()
	   <<setw(16)<<(*it).getPenetration()
	   <<setw(16)<<(*it).getTgtDisp()
	   <<setw(16)<<(*it).getContactRadius()
	   <<setw(16)<<(*it).getR0()
	   <<setw(16)<<(*it).getE0()
	   <<setw(16)<<(*it).getNormalForce()
	   <<setw(16)<<(*it).getTgtForce()
	   <<setw(16)<<( (*it).getPoint1().getx()+(*it).getPoint2().getx() )/2
	   <<setw(16)<<( (*it).getPoint1().gety()+(*it).getPoint2().gety() )/2
	   <<setw(16)<<( (*it).getPoint1().getz()+(*it).getPoint2().getz() )/2
	   <<setw(16)<<(*it).NormalForceVec().getx()
	   <<setw(16)<<(*it).NormalForceVec().gety()
	   <<setw(16)<<(*it).NormalForceVec().getz()
	   <<setw(16)<<(*it).TgtForceVec().getx()
	   <<setw(16)<<(*it).TgtForceVec().gety()
	   <<setw(16)<<(*it).TgtForceVec().getz()
	   <<endl;
    ofs.close();
}


void assembly::snapshot(const char* str) const{
    ofstream ofs(str);
    if(!ofs) {
	cout<<"stream error in snapshot!"<<endl; exit(-1);
    }
    ofs.setf(ios::scientific, ios::floatfield);
    ofs<<"TotalNum="<<setw(10)<<TotalNum<<endl;
    vec tmp;
    ofs<<"     ID         radius_a        radius_b        radius_c   "
       <<"   position_x     position_y     position_z    "
       <<"      a_cosx          a_cosy          a_cosz     "
       <<"      b_cosx         b_cosy          b_cosz     "
       <<"      c_cosx         c_cosy          c_cosz     "
       <<"   velocity_x      velocity_y      velocity_z   "
       <<"     omga_x          omga_y          omga_z     "
       <<"    force_x         force_y         force_z     "
       <<"    moment_x        moment_y        moment_z    "
       <<"     ID"<<endl;    
    list<particle*>::const_iterator it;
    for (it=ParticleList.begin();it!=ParticleList.end();++it)
    {
	ofs<<setw(10)<<(*it)->getID()
	   <<setw(16)<<(*it)->getA()
	   <<setw(16)<<(*it)->getB()
	   <<setw(16)<<(*it)->getC();

	tmp=(*it)->getCurrPosition();
	ofs<<setw(16)<<tmp.getx()
	   <<setw(16)<<tmp.gety()
	   <<setw(16)<<tmp.getz();
	
	tmp=(*it)->getCurrDirecA();
	ofs<<setw(16)<<tmp.getx()
	   <<setw(16)<<tmp.gety()
	   <<setw(16)<<tmp.getz();

	tmp=(*it)->getCurrDirecB();
	ofs<<setw(16)<<tmp.getx()
	   <<setw(16)<<tmp.gety()
	   <<setw(16)<<tmp.getz();

	tmp=(*it)->getCurrDirecC();
	ofs<<setw(16)<<tmp.getx()
	   <<setw(16)<<tmp.gety()
	   <<setw(16)<<tmp.getz();

	tmp=(*it)->getCurrVelocity();
	ofs<<setw(16)<<tmp.getx()
	   <<setw(16)<<tmp.gety()
	   <<setw(16)<<tmp.getz();

	tmp=(*it)->getCurrOmga();
	ofs<<setw(16)<<tmp.getx()
	   <<setw(16)<<tmp.gety()
	   <<setw(16)<<tmp.getz();

	tmp=(*it)->getForce();
	ofs<<setw(16)<<tmp.getx()
	   <<setw(16)<<tmp.gety()
	   <<setw(16)<<tmp.getz();

	tmp=(*it)->getMoment();
	ofs<<setw(16)<<tmp.getx()
	   <<setw(16)<<tmp.gety()
	   <<setw(16)<<tmp.getz();

	ofs<<setw(10)<<(*it)->getID()<<endl;
    }
    ofs.close();
}

	
void assembly::createSample(const char* str){
    ifstream ifs(str);
    if(!ifs) {
	cout<<"stream error in createSample!"<<endl; exit(-1);
    }
    ifs >> TotalNum >> RORC;

    long double cx,cy,cz,rd,wd,lt,ht;
    if(RORC==0){
	ifs >> cx >> cy >> cz >> rd >> ht;
	S.set_center(vec(cx,cy,cz));
	S.set_radius(rd);
	S.set_height(ht);
	Volume = PI * powl(S.get_radius(),2) * S.get_height();
    }
    else{
	ifs >> cx >> cy >> cz >> wd >> lt >> ht;
	R.set_center(vec(cx,cy,cz));
	R.set_width(wd);
	R.set_length(lt);
	R.set_height(ht);
	Volume = R.get_width() * R.get_length() * R.get_height();
    }
    
    char s[20];
    ifs>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
       >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s;

    ParticleList.clear();
    int ID, type;
    long double a, b, c, px,py,pz,dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz;
    long double vx,vy,vz,omx,omy,omz,fx,fy,fz,mx,my,mz;
    for (int i=0;i<TotalNum;i++){
	ifs>>ID>>type>>a>>b>>c>>px>>py>>pz>>dax>>day>>daz>>dbx>>dby>>dbz>>dcx>>dcy>>dcz
	   >>vx>>vy>>vz>>omx>>omy>>omz>>fx>>fy>>fz>>mx>>my>>mz;
	particle* pt= new particle(ID,type,vec(a,b,c),vec(px,py,pz),vec(dax,day,daz),vec(dbx,dby,dbz),vec(dcx,dcy,dcz));

//      optional settings for a particle's initial status
//	pt->setPrevVelocity(vec(vx,vy,vz));
//	pt->setCurrVelocity(vec(vx,vy,vz));
//	pt->setPrevOmga(vec(omx,omy,omz));
//	pt->setCurrOmga(vec(omx,omy,omz));
//	pt->setConstForce(vec(fx,fy,fz));  // constant force, not initial force
//	pt->setConstMoment(vec(mx,my,mz)); // constant moment, not initial moment

	ParticleList.push_back(pt);
    }
    ifs.close();
}


#ifdef OPENMP	
void assembly::findContact(){ // OpenMP version
    ContactList.clear();
    PossCntctNum = 0;

#ifdef TIME_PROFILE
    gettimeofday(&timep1,NULL); 
#endif
    int iam, nt, i, j, ipoints, npoints, rpoints;
    list<particle*>::iterator ot, it, pt;
    vec u,v;
    npoints = ParticleList.size();
    ot = ParticleList.begin();

#pragma omp parallel num_threads(NUM_THREADS), private(iam, nt, i, j, ipoints, it, pt, u, v)
    {
	iam = omp_get_thread_num();
	nt  = omp_get_num_threads();
	ipoints = npoints/nt;  // divide into nt partitions
	rpoints = npoints%nt;  // remainder of the division
	it = ot;

	// determine starting point and extend of each partition
	if (rpoints == 0) {
	  for (i = 0; i < iam * ipoints; ++i)
	    ++it;              // starting point of each partition
	}
	else {
	  if (iam < rpoints) {
	    ipoints += 1;      // ipoints changed
	    for (i = 0; i < iam * ipoints ; ++i)
	      ++it;
	  }
	  else {
	    for (i = 0; i < rpoints * (ipoints + 1) + (iam - rpoints) * ipoints; ++ i)
	      ++it;
	  }

	}

	// explore each partition
	for (j=0;j<ipoints;++j,++it) { 
	    u=(*it)->getCurrPosition();
	    for (pt=it,++pt;pt!=ParticleList.end();++pt){
		v=(*pt)->getCurrPosition();
		if (   ( vfabsl(v-u) < (*it)->getA() + (*pt)->getA())
		    && ( (*it)->getType() !=  1 || (*pt)->getType() != 1  )      // not both are fixed particles
		    && ( (*it)->getType() !=  5 || (*pt)->getType() != 5  )      // not both are free boundary particles
		    && ( (*it)->getType() != 10 || (*pt)->getType() != 10 )  ) { // not both are ghost particles
		    contact<particle> tmpct(*it, *pt); // a local and temparory object
		    ++PossCntctNum;
		    if(tmpct.isOverlapped())
#pragma omp critical
			ContactList.push_back(tmpct);    // containers use value semantics, so a "copy" is pushed back.
		}
	    }
	}
    }
    
#ifdef TIME_PROFILE
    gettimeofday(&timep2,NULL);
    g_debuginf<<setw(10)<<timediffsec(); 
#endif
	
    ActualCntctNum = ContactList.size();
}

#else
void assembly::findContact(){ // serial version
    ContactList.clear();
    PossCntctNum = 0;

#ifdef TIME_PROFILE
    gettimeofday(&timep1,NULL); 
#endif
    list<particle*>::iterator it, pt;
    vec u,v;
    for (it=ParticleList.begin();it!=ParticleList.end();++it){
	u=(*it)->getCurrPosition();
	for (pt=it,++pt;pt!=ParticleList.end();++pt){
	    v=(*pt)->getCurrPosition();
	    if (   ( vfabsl(v-u) < (*it)->getA() + (*pt)->getA())
		&& ( (*it)->getType() !=  1 || (*pt)->getType() != 1  )      // not both are fixed particles
		&& ( (*it)->getType() !=  5 || (*pt)->getType() != 5  )      // not both are free boundary particles
		&& ( (*it)->getType() != 10 || (*pt)->getType() != 10 )  ) { // not both are ghost particles
		contact<particle> tmpct(*it, *pt); // a local and temparory object
		++PossCntctNum;
		if(tmpct.isOverlapped())
		  ContactList.push_back(tmpct);    // containers use value semantics, so a "copy" is pushed back.
	    }
	}
    }	

#ifdef TIME_PROFILE
    gettimeofday(&timep2,NULL);
    g_debuginf<<setw(10)<<timediffsec(); 
#endif
 
    ActualCntctNum = ContactList.size();
}
#endif

/* another OpenMP version, simpler but slower

#ifdef OPENMP	
void assembly::findContact(){
    ContactList.clear();
    PossCntctNum = 0;

#ifdef TIME_PROFILE
    gettimeofday(&timep1,NULL); 
#endif
    list<particle*>::iterator ot, it, pt;
    vec u,v;
    int i,j,n;
    n =ParticleList.size();
    ot=ParticleList.begin();

#pragma omp parallel for private(j, it, pt, u, v)    
    for (i=0; i < n; i++){
      it=ot;
      for (j=0; j < i; j++)
	++it;
      u=(*it)->getCurrPosition();
      for (pt=it,++pt;pt!=ParticleList.end();++pt){
	v=(*pt)->getCurrPosition();
	if (vfabsl(v-u) < (*it)->getA() + (*pt)->getA() ) {
	  contact<particle> tmpct(*it, *pt); // a local and temparory object
	  ++PossCntctNum;
	  if(tmpct.isOverlapped())
#pragma omp critical
	      ContactList.push_back(tmpct);    // containers use value semantics, so a "copy" is pushed back.
	}
      }
    }
#ifdef TIME_PROFILE
    gettimeofday(&timep2,NULL);
    g_debuginf<<setw(10)<<timediffsec(); 
#endif
	
    ActualCntctNum = ContactList.size();
}

#else

*/

long double assembly::getDensity() const{
    long double dens=0;
    list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it)
	dens+=(*it)->getMass();
    return dens/=Volume;
}


long double assembly::getAveragePenetration() const{
    int totalcntct = ContactList.size();
    if (totalcntct==0)
	return 0;
    else {
	long double pene=0;
	for (list<CONTACT>::const_iterator it=ContactList.begin();it!=ContactList.end();++it)
	    pene += it->getPenetration(); 
	return pene/totalcntct;
    }
}


long double assembly::getMinTimeStep() const {
    int totalcntct = ContactList.size();
    if (totalcntct == 0)
	return 0;
    else {
	list<CONTACT>::const_iterator it=ContactList.begin();
        long double minTimeStep = it->getTimeStep();
	for (++it; it != ContactList.end(); ++it) {
	  long double val = it->getTimeStep(); 
	  minTimeStep =  val < minTimeStep ? val : minTimeStep;
	}
	return minTimeStep;
    }
 }


long double assembly::getAverageVelocity() const{
    long double avgv=0;
    int count=0;
    list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it)
	if ((*it)->getType()==0) {
	    avgv+=vfabsl((*it)->getCurrVelocity());
	    count++;
	}
    return avgv/=count;
}


long double assembly::getAverageOmga() const{
    long double avgv=0;
    int count=0;
    list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it)
	if ((*it)->getType()==0){
	    avgv+=vfabsl((*it)->getCurrOmga());
	    count++;
	}
    return avgv/=count;
}


long double assembly::getAverageForce() const{
    long double avgv=0;
    int count=0;
    list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it)
	if ((*it)->getType()==0){
	    avgv+=vfabsl((*it)->getForce());
	    count++;
	}
    return avgv/count;
}


long double assembly::getAverageMoment() const{
    long double avgv=0;
    int count=0;
    list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it)
	if ((*it)->getType()==0){
	    avgv+=vfabsl((*it)->getMoment());
	    count++;
	}
    return avgv/=count;
}


long double assembly::getParticleVolume() const{
    long double avgv=0;
    list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it)
	if ((*it)->getType()==0)
	    avgv+=(*it)->getVolume();
    return avgv;
}


vec assembly::getTopFreeParticlePosition() const{
    list<particle*>::const_iterator it,jt,kt;
    it=ParticleList.begin();
    while (it!=ParticleList.end() && (*it)->getType()!=0)   // find the 1st free particle
	++it;

    if (it==ParticleList.end())    // no free particles
	return 0;

    jt=it; 
    kt=it;
    
    // two cases:
    // 1: 1st particle is not free
    // 2: 1st particle is free
    if (++kt!=ParticleList.end()){ // case1: more than 2 particles; case 2: more than 1 particle
	for(++it;it!=ParticleList.end();++it){
	    if ((*it)->getType()==0)
		if ((*it)->getCurrPosition().getz() > (*jt)->getCurrPosition().getz())
		    jt=it;
	}
	return (*jt)->getCurrPosition();
    }
    else {
	if ((*it)->getType()==0)  // case1: only 2 particles, the 2nd one is free; case2: only 1 particle
	    return (*it)->getCurrPosition();
	else
	    return 0;
    }

}


long double assembly::ellipPileForce() {
    long double val=0;
    for(list<particle*>::iterator it=ParticleList.begin();it!=ParticleList.end();++it)
	if ((*it)->getType()==3) {
	    val = (*it)->getForce().getz();
	    break;
	}
    return val;
}


vec assembly::ellipPileDimn() {
    vec val;
    for(list<particle*>::iterator it=ParticleList.begin();it!=ParticleList.end();++it)
	if ((*it)->getType()==3) {
	    val = vec((*it)->getA(), (*it)->getB(), (*it)->getC());
	    break;
	}
    return val;
}


long double assembly::ellipPileTipZ() {
    long double val=0;
    for(list<particle*>::iterator it=ParticleList.begin();it!=ParticleList.end();++it)
	if ((*it)->getType()==3) {
	    val = (*it)->getCurrPosition().getz()-(*it)->getA();
	    break;
	}
    return val;
}


long double assembly::ellipPilePeneVol() {
    long double val=0;
    if (getTopFreeParticlePosition().getz()-ellipPileTipZ()<=0)
	val=0;
    else{
	// low: a signed number as lower limit for volumetric integration
	long double low=ellipPileTipZ() + ellipPileDimn().getx() - getTopFreeParticlePosition().getz(); 
	long double lowint=low-powl(low,3)/3.0/powl(ellipPileDimn().getx(),2);
	val = PI * ellipPileDimn().gety() * ellipPileDimn().getz()
	      *(2.0/3*ellipPileDimn().getx()-lowint);
    }
    return val;
}


void assembly::ellipPileUpdate(){
    for(list<particle*>::iterator it=ParticleList.begin();it!=ParticleList.end();++it){
	if ((*it)->getType()==3) {
	    (*it)->curr_velocity.setx(0);	
	    (*it)->curr_velocity.sety(0);
	    (*it)->curr_velocity.setz(-PILE_RATE);
	    (*it)->curr_position = (*it)->prev_position + (*it)->curr_velocity*TIMESTEP;
	}
    }
}


long double assembly::getTransEnergy() const{
    long double engy=0;
    list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it){
	if ((*it)->getType()==0)
	    engy+=(*it)->getTransEnergy();
    }
    return engy;
}


long double assembly::getRotatEnergy() const{
    long double engy=0;
    list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it){
	if ((*it)->getType()==0)
	    engy+=(*it)->getRotatEnergy();
    }
    return engy;
}


long double assembly::getKinetEnergy() const{
    long double engy=0;
    list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it){
	if ((*it)->getType()==0)
	    engy+=(*it)->getKinetEnergy();
    }
    return engy;
}


long double assembly::getPotenEnergy(long double ref) const{
    long double engy=0;
    list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it){
	if ((*it)->getType()==0)
	    engy+=(*it)->getPotenEnergy(ref);
    }
    return engy;
}


void assembly::clearForce(){
    for(list<particle*>::iterator it=ParticleList.begin();it!=ParticleList.end();++it){
	(*it)->clearForce();
    }
}


void assembly::flexiBoundaryForceZero(){
    for(list<particle*>::iterator it=ParticleList.begin();it!=ParticleList.end();++it){
	(*it)->flb_force=0;
	(*it)->flb_moment=0;
    }
}


void assembly::initFBForce(){
    for(list<particle*>::iterator it=ParticleList.begin();it!=ParticleList.end();++it){
	(*it)->force+=(*it)->flb_force;
	(*it)->moment+=(*it)->flb_moment;
    }
}


void assembly::internalForce(long double& avgnm, long double& avgsh){
    avgnm=0;
    avgsh=0;

    int totalcntct = ContactList.size();
    if(totalcntct==0){
	avgnm = 0;
	avgsh = 0;
    }
    else{
	list<CONTACT>::iterator it;
	for (it=ContactList.begin();it!=ContactList.end();++it)
	    it->checkinPreTgt(CntTgtVec); // checkin previous tangential force and displacment    
	
	CntTgtVec.clear(); // CntTgtVec must be cleared before filling in new values.

#ifdef TIME_PROFILE
	gettimeofday(&timep1,NULL); 
#endif 
	for (it=ContactList.begin();it!=ContactList.end();++it){
	    it->contactForce();           // cannot be parallelized as it may change a particle's force simultaneously.
	    it->checkoutTgt(CntTgtVec);   // checkout current tangential force and displacment
	    avgnm += it->getNormalForce();
	    avgsh += it->getTgtForce();
	}
	avgnm /= totalcntct;
	avgsh /= totalcntct;

#ifdef TIME_PROFILE
	gettimeofday(&timep2,NULL);
	g_debuginf<<setw(10)<<timediffsec()<<endl; 
#endif

    }
}


void assembly::updateParticle(){
    for(list<particle*>::iterator it=ParticleList.begin();it!=ParticleList.end();++it){
	(*it)->update();
    }
}


void assembly::createRigidBoundary(ifstream &ifs){
    rgd_bdry<particle>* rbptr;
    int type;
    RBList.clear();
    ifs>>RgdBdryNum;
    
    for(int i=0;i<RgdBdryNum;i++){
	ifs>>type;
	if(type==1) // plane boundary
	    rbptr=new plnrgd_bdry<particle>(ifs);
	else        // cylindrical boundary
	    rbptr=new cylrgd_bdry<particle>(ifs);
	RBList.push_back(rbptr);
    }
}


void assembly::createFlexiBoundary(ifstream &ifs){
    flb_bdry<particle>* fbptr;
    int type;
    FBList.clear();
    ifs>>FlbBdryNum;
    
    for(int i=0;i<FlbBdryNum;i++){
	ifs>>type;
	if(type==1) // plane boundary
	    fbptr=new plnflb_bdry<particle>(ifs);
	else        // cylindrical boundary
	    fbptr=new cylflb_bdry<particle>(ifs);
	FBList.push_back(fbptr);
    }
}

	
void assembly::createBoundary(const char* str){
    ifstream ifs(str);
    if(!ifs) {
	cout<<"stream error in createBoundary!"<<endl; exit(-1);
    }
    ifs >> BdryType;
    if(BdryType==0){      // rigid boundaries
	createRigidBoundary(ifs);
    }
    else if(BdryType==1){ // flexible boundaries
	createFlexiBoundary(ifs);
    }
    ifs.close();
}


void assembly::displayBoundary() const{
    printf("rgd boundary number%5d\n",RgdBdryNum);
    list<RGDBDRY*>::const_iterator rt;
    list<FLBBDRY*>::const_iterator ft;
    for(rt=RBList.begin();rt!=RBList.end();++rt)
	(*rt)->disp();
    for(ft=FBList.begin();ft!=FBList.end();++ft)
	(*ft)->disp();
}


void assembly::printBoundary(const char* str) const
{
    ofstream ofs(str);
    if(!ofs) {
	cout<<"stream error in printBoundary!"<<endl; exit(-1);
    }
    ofs.setf(ios::scientific, ios::floatfield);

    ofs<<setw(10)<<BdryType
       <<setw(10)<<RgdBdryNum<<endl;
    
    list<RGDBDRY*>::const_iterator rt;
    for(rt=RBList.begin();rt!=RBList.end();++rt)
	(*rt)->disp(ofs);
    ofs<<endl;

    ofs.close();
}


void assembly::findParticleOnBoundary(){
    list<RGDBDRY*>::iterator rt;
    list<FLBBDRY*>::iterator ft;
    for(rt=RBList.begin();rt!=RBList.end();++rt)
	(*rt)->findParticleOnBoundary(ParticleList);
    for(ft=FBList.begin();ft!=FBList.end();++ft)
	(*ft)->findParticleOnBoundary(ParticleList);
}


void assembly::findParticleOnLine(){
    list<FLBBDRY*>::iterator ft;
    for(ft=FBList.begin();ft!=FBList.end();++ft)
	(*ft)->findParticleOnLine();
}


void assembly::createFlbNet(){
    list<FLBBDRY*>::iterator ft;
    for(ft=FBList.begin();ft!=FBList.end();++ft)
	(*ft)->createFlbNet();
}


void assembly::rigidBoundaryForce(){
  list<RGDBDRY*>::iterator rt;
  for(rt=RBList.begin();rt!=RBList.end();++rt)
    (*rt)->rigidBF(BdryTgtMap);

  /*
  vector<boundarytgt>::iterator it;
  list<RGDBDRY*>::iterator rt;

  for(rt=RBList.begin();rt!=RBList.end();++rt){	
    (*rt)->rigidBF(BdryTgtMap);
    for (it=BdryTgtMap[(*rt)->bdry_id].begin();it!=BdryTgtMap[(*rt)->bdry_id].end();++it){
      g_debuginf<<setw(10)<<g_iteration
		<<setw(10)<<(*rt)->bdry_id
		<<setw(10)<<BdryTgtMap[(*rt)->bdry_id].size()
		<<setw(16)<<it->TgtForce.getx()
		<<setw(16)<<it->TgtForce.gety()
		<<setw(16)<<it->TgtForce.getz()
		<<endl;
      //<<setw(16)<<it->TgtPeak<<endl;
    }
  }
  */
}


void assembly::rigidBoundaryForce(long double penetr[],int cntnum[]){
  list<RGDBDRY*>::iterator rt;
  for(rt=RBList.begin();rt!=RBList.end();++rt){	
    (*rt)->rigidBF(BdryTgtMap);
    if ((*rt)->getBdryID()==1){
	penetr[1] = (*rt)->getAvgPenetr();
	cntnum[1] = (*rt)->getCntnum();
    }
    else if ((*rt)->getBdryID()==2){
	penetr[2] = (*rt)->getAvgPenetr();
	cntnum[2] = (*rt)->getCntnum();
    }
    else if ((*rt)->getBdryID()==3){
	penetr[3] = (*rt)->getAvgPenetr();
	cntnum[3] = (*rt)->getCntnum();
    }
    else if ((*rt)->getBdryID()==4){
	penetr[4] = (*rt)->getAvgPenetr();
	cntnum[4] = (*rt)->getCntnum();
    }
    else if ((*rt)->getBdryID()==5){
	penetr[5] = (*rt)->getAvgPenetr();
	cntnum[5] = (*rt)->getCntnum();
    }
    else if ((*rt)->getBdryID()==6){
	penetr[6] = (*rt)->getAvgPenetr();
	cntnum[6] = (*rt)->getCntnum();
    }
  }
}


void assembly::flexiBoundaryForce(){
    list<FLBBDRY*>::iterator ft;
    for(ft=FBList.begin();ft!=FBList.end();++ft)
	(*ft)->flxbBF();
}


vec assembly::getNormalForce(int bdry) const{
    list<RGDBDRY*>::const_iterator it;
    for(it=RBList.begin();it!=RBList.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getNormalForce();
    }
    return 0;
}


vec assembly::getShearForce(int bdry) const{
    list<RGDBDRY*>::const_iterator it;
    for(it=RBList.begin();it!=RBList.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getShearForce();
    }
    return 0;
}


long double assembly::getAvgNormal(int bdry) const{
    list<RGDBDRY*>::const_iterator it;
    for(it=RBList.begin();it!=RBList.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getAvgNormal();
    }
    return 0;
}


vec assembly::getApt(int bdry) const{
    list<RGDBDRY*>::const_iterator it;
    for(it=RBList.begin();it!=RBList.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getApt();
    }
    return 0;
}


vec assembly::getDirc(int bdry) const{
    list<RGDBDRY*>::const_iterator it;
    for(it=RBList.begin();it!=RBList.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getDirc();
    }
    return 0;
}


long double assembly::getArea(int n) const{
    list<RGDBDRY*>::const_iterator it;
    for(it=RBList.begin();it!=RBList.end();++it){
	if((*it)->getBdryID()==n)
	    return (*it)->area;
    }
    return 0;
}


void assembly::setArea(int n, long double a){
    list<RGDBDRY*>::iterator it;
    for(it=RBList.begin();it!=RBList.end();++it){
	if((*it)->getBdryID()==n)
	    (*it)->area=a;
    }
}


long double assembly::getAverageRigidPressure() const{
    list<RGDBDRY*>::const_iterator rt;
    long double avgpres=0;
    for(rt=RBList.begin();rt!=RBList.end();++rt)
	avgpres+=vfabsl((*rt)->getNormalForce())/(*rt)->getArea();
    return avgpres/=RgdBdryNum;
}


// only update CoefOfLimits[0] for specified boundaries
void assembly::updateRB(int bn[], UPDATECTL rbctl[], int num){
    for(int i=0;i<num;i++){
	for(list<RGDBDRY*>::iterator rt=RBList.begin();rt!=RBList.end();++rt){
	    if((*rt)->getBdryID()==bn[i]){
		(*rt)->update(rbctl[i]);
		break;
	    }
	}
    }
}


// update CoefOfLimits[1,2,3,4] for all 6 boundaries
void assembly::updateRB6(){
    for(list<RGDBDRY*>::iterator rt=RBList.begin();rt!=RBList.end();++rt){
	if((*rt)->getBdryID()==1 || (*rt)->getBdryID()==3){
	    for(list<RGDBDRY*>::iterator lt=RBList.begin();lt!=RBList.end();++lt){
		if((*lt)->getBdryID()==4)
		    (*rt)->CoefOfLimits[1].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==2)
		    (*rt)->CoefOfLimits[2].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==5)
		    (*rt)->CoefOfLimits[3].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==6)
		    (*rt)->CoefOfLimits[4].apt=(*lt)->CoefOfLimits[0].apt;
	    }
	}
	else if((*rt)->getBdryID()==2 || (*rt)->getBdryID()==4){
	    for(list<RGDBDRY*>::iterator lt=RBList.begin();lt!=RBList.end();++lt){
		if((*lt)->getBdryID()==1)
		    (*rt)->CoefOfLimits[1].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==3)
		    (*rt)->CoefOfLimits[2].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==5)
		    (*rt)->CoefOfLimits[3].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==6)
		    (*rt)->CoefOfLimits[4].apt=(*lt)->CoefOfLimits[0].apt;
	    }

	}
	else if((*rt)->getBdryID()==5 || (*rt)->getBdryID()==6){
	    for(list<RGDBDRY*>::iterator lt=RBList.begin();lt!=RBList.end();++lt){
		if((*lt)->getBdryID()==1)
		    (*rt)->CoefOfLimits[1].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==3)
		    (*rt)->CoefOfLimits[2].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==2)
		    (*rt)->CoefOfLimits[3].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==4)
		    (*rt)->CoefOfLimits[4].apt=(*lt)->CoefOfLimits[0].apt;
	    }

	}
	
    }
}


// upgrade CoefOfLimits[1,2,3,4] for rectangular pile
void assembly::updateRectPile(){
    for(list<RGDBDRY*>::iterator rt=RBList.begin();rt!=RBList.end();++rt){
	if((*rt)->getBdryID()==7 || (*rt)->getBdryID()==9 ){
	    for(list<RGDBDRY*>::iterator lt=RBList.begin();lt!=RBList.end();++lt){
		if((*lt)->getBdryID()==10)
		    (*rt)->CoefOfLimits[1].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==8)
		    (*rt)->CoefOfLimits[2].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==11)
		    (*rt)->CoefOfLimits[3].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==12)
		    (*rt)->CoefOfLimits[4].apt=(*lt)->CoefOfLimits[0].apt;
	    }
	}
	else if((*rt)->getBdryID()==8 || (*rt)->getBdryID()==10){
	    for(list<RGDBDRY*>::iterator lt=RBList.begin();lt!=RBList.end();++lt){
		if((*lt)->getBdryID()==7)
		    (*rt)->CoefOfLimits[1].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==9)
		    (*rt)->CoefOfLimits[2].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==11)
		    (*rt)->CoefOfLimits[3].apt=(*lt)->CoefOfLimits[0].apt;
		else if((*lt)->getBdryID()==12)
		    (*rt)->CoefOfLimits[4].apt=(*lt)->CoefOfLimits[0].apt;
	    }
	}
    }
}


void assembly::updateFB(int bn[], UPDATECTL fbctl[], int num){
    list<FLBBDRY*>::iterator ft;
    int i,k=1;
    for(i=0;i<num;i++){
	for(ft=FBList.begin();ft!=FBList.end();++ft){
	    if(k++==bn[i]){
		UPDATECTL ctl[2];
		ctl[0]=fbctl[bn[i]*2-2];
		ctl[1]=fbctl[bn[i]*2-1];
		(*ft)->update(ctl,2);
		break;
	    }
	}
    }
}


// create a specimen from discreate particles through floating and then gravitation,
// file cre_particle contains the final particle information,
// file cre_boundary contains the final boundary information.
void assembly::deposit_RgdBdry(gradation& grad,
			       int   freetype,
			       int   total_steps,  
			       int   snapshots,
			       int   interval,
			       long double height,
			       const char* iniptclfile,   
			       const char* inibdryfile,
			       const char* particlefile, 
			       const char* contactfile,
			       const char* progressfile, 
			       const char* creparticle,
			       const char* creboundary,
			       const char* debugfile)
{
    if (grad.rorc == 1) {
	RORC = grad.rorc;
	R.set_center(vec(0,0,0));
	R.set_width(grad.dimn);
	R.set_length(grad.dimn);
	R.set_height(grad.dimn);
	
	generate(grad, iniptclfile, freetype, height); 
        // 3.0 for uniform size of (2.5e-3,*0.8,*0.6); if not uniform, it may be larger, say, 4.5
        // 3.0 for uniform spheres of 2.5e-3

	setBoundary(grad.rorc, 5, grad.dimn, inibdryfile);

	deposit(total_steps,        // total_steps
		snapshots,          // number of snapshots
		interval,           // record interval
		iniptclfile,        // input file, initial particles
		inibdryfile,        // input file, initial boundaries
		particlefile,       // output file, resulted particles, including snapshots 
		contactfile,        // output file, resulted contacts, including snapshots 
		progressfile,       // output file, statistical info
		debugfile);         // output file, debug info

	setBoundary(grad.rorc,      // rectangular--1 or cylindrical--0?
		6,                  
		grad.dimn,          // specimen dimension
		"trm_boundary");    // output file, containing boundaries info
	
	trim(grad.rorc,             // rectangular--1 or cylindrical--0?
	     "dep_particle_end",    // input file, particles to be trimmed
	     "trm_boundary",  
	     creparticle,
	     creboundary);

    }
}


// freetype:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
// ht- how many times of size would be the floating height
void assembly::generate(gradation&  grad,
			const char* particlefile,
			int freetype,
			long double ht)
{
    long double x,y,z;
    particle* newptcl;
    TotalNum = 0;
    long double est =1.02;
    int grid=9;  
    // grid: dimension of free particle array.
    // 7 - small dimn container
    // 9 - medium dimn container 
    // 11- large dimn container 

    long double dimn=grad.dimn;
    if (freetype == 0) {      // just one free particle
	newptcl = new particle(TotalNum+1, 0, vec(dimn/2/40,dimn/2/20,dimn/2), grad);
	ParticleList.push_back(newptcl);
	TotalNum++;
    }
    else if (freetype == 1) { // a horizontal layer of free particles
	z=dimn/2;
	for (x=-dimn/2*(grid-1)/10; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
	    for (y=-dimn/2*(grid-1)/10; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5){
		newptcl = new particle(TotalNum+1, 0, vec(x,y,z), grad);
		ParticleList.push_back(newptcl);
		TotalNum++;
	    }
    }
    else if (freetype == 2) { // multiple layers of free particles
	long double offset=0; // 0 for ellipsoids; dimn/2/5/5 for spheres
	if (grad.ratio_ba==1.0 && grad.ratio_ca==1.0)
	    offset = dimn/2/5/5;
	for (z=dimn/2; z<dimn/2 + dimn*ht; z+=dimn/2/5) {
	//for (z=-dimn/2*4/5; z<dimn/2 + dimn*ht; z+=dimn/2/10) { // spheres
	    for (x=-dimn/2*(grid-1)/10+offset; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
		for (y=-dimn/2*(grid-1)/10+offset; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5){
		    newptcl = new particle(TotalNum+1, 0, vec(x,y,z), grad);
		    ParticleList.push_back(newptcl);
		    TotalNum++;
		}	
	    offset *= -1;
	}
    }

    printParticle(particlefile);

}


// create a specimen from discreate particles through floating and then gravitation,
// boundaries are composed of fixed particles.
void assembly::deposit_PtclBdry(gradation& grad,
				int   freetype,
				long double rsize,
				int   total_steps,  
				int   snapshots,
				int   interval,
				const char* iniptclfile,   
				const char* particlefile, 
				const char* contactfile,
				const char* progressfile, 
				const char* debugfile)
{
    if (grad.rorc == 1) {
	RORC = grad.rorc;
	R.set_center(vec(0,0,0));
	R.set_width(grad.dimn);
	R.set_length(grad.dimn);
	R.set_height(grad.dimn);
	
	generate_p(grad, iniptclfile, freetype, rsize, 4.0);
	deposit_p(total_steps,        // total_steps
		  snapshots,          // number of snapshots
		  interval,           // record interval
		  grad.dimn,          // dimension of particle-composed-boundary
		  rsize,              // relative container size
		  iniptclfile,        // input file, initial particles
		  particlefile,       // output file, resulted particles, including snapshots 
		  contactfile,        // output file, resulted contacts, including snapshots 
		  progressfile,       // output file, statistical info
		  debugfile);         // output file, debug info
    }
}


// freetype:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
// ht- how many times of size would be the floating height
void assembly::generate_p(gradation&  grad,
			 const char* particlefile,
			 int freetype,
			 long double rsize,
			 long double ht)
{
    long double x,y,z;
    particle* newptcl;
    TotalNum = 0;
    long double wall=2.2; // wall - wall height; ht - free particle height
    long double est =1.02;
    int grid=static_cast<int> (nearbyint(rsize*10)-1);  

    // grid: dimension of free particle array.
    // 7 - small dimn container
    // 9 - medium dimn container 
    // 11- large dimn container 

    long double dimn=grad.dimn;
    // particle boundary 1
    x=dimn/2*(grid+1)/10;
    for (y=-dimn/2*grid/10; y<dimn/2*grid/10*est; y+=dimn/2/5)
	for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5){
	    newptcl = new particle(TotalNum+1, 1, vec(x,y,z), grad.ptclsize[0]*0.99);
	    ParticleList.push_back(newptcl);
	    TotalNum++;
	}

    // particle boundary 2
    y=dimn/2*(grid+1)/10;
    for (x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5)
	for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5){
	    newptcl = new particle(TotalNum+1, 1, vec(x,y,z), grad.ptclsize[0]*0.99 );
	    ParticleList.push_back(newptcl);
	    TotalNum++;
	}

    // particle boundary 3
    x=-dimn/2*(grid+1)/10;
    for (y=-dimn/2*grid/10; y<dimn/2*grid/10*est; y+=dimn/2/5)
	for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5){
	    newptcl = new particle(TotalNum+1, 1, vec(x,y,z), grad.ptclsize[0]*0.99);
	    ParticleList.push_back(newptcl);
	    TotalNum++;
	}

    // particle boundary 4
    y=-dimn/2*(grid+1)/10;
    for (x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5)
	for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5){
	    newptcl = new particle(TotalNum+1, 1, vec(x,y,z), grad.ptclsize[0]*0.99);
	    ParticleList.push_back(newptcl);
	    TotalNum++;
	}

    // particle boundary 6
    z=-dimn/2;
    for (y=-dimn/2*grid/10; y<dimn/2*grid/10*est; y+=dimn/2/5)
	for( x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5){
	    newptcl = new particle(TotalNum+1, 1, vec(x,y,z), grad.ptclsize[0]*0.99);
	    ParticleList.push_back(newptcl);
	    TotalNum++;
	}

    if (freetype == 0) {      // just one free particle
	newptcl = new particle(TotalNum+1, 0, vec(dimn/2/40,dimn/2/20,dimn/2), grad);
	ParticleList.push_back(newptcl);
	TotalNum++;
    }
    else if (freetype == 1) { // a horizontal layer of free particles
	z=dimn/2;
	for (x=-dimn/2*(grid-1)/10; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
	    for (y=-dimn/2*(grid-1)/10; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5){
		newptcl = new particle(TotalNum+1, 0, vec(x,y,z), grad);
		ParticleList.push_back(newptcl);
		TotalNum++;
	    }
    }
    else if (freetype == 2) { // multiple layers of free particles
	for (z=dimn/2; z<dimn/2 + dimn*ht; z+=dimn/2/5)
	    for (x=-dimn/2*(grid-1)/10; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
		for (y=-dimn/2*(grid-1)/10; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5){
		    newptcl = new particle(TotalNum+1, 0, vec(x,y,z), grad);
		    ParticleList.push_back(newptcl);
		    TotalNum++;
		}	
    }
    
    printParticle(particlefile);
    
}


void assembly::scale_PtclBdry(int   total_steps,  
			      int   snapshots,
			      int   interval,
			      long double dimn,
			      long double rsize,
			      const char* iniptclfile,   
			      const char* particlefile, 
			      const char* contactfile,
			      const char* progressfile, 
			      const char* debugfile)
{
    deposit_p(total_steps,        // total_steps
	      snapshots,          // number of snapshots
	      interval,           // record interval
	      dimn,               // dimension of particle-composed-boundary
	      rsize,              // relative container size
	      iniptclfile,        // input file, initial particles
	      particlefile,       // output file, resulted particles, including snapshots 
	      contactfile,        // output file, resulted contacts, including snapshots 
	      progressfile,       // output file, statistical info
	      debugfile);         // output file, debug info
}


// collapse a deposited specimen through gravitation
void assembly::collapse(int   rors, 
			int   total_steps,  
			int   snapshots,
			int   interval,
			const char* iniptclfile,
			const char* initboundary,
			const char* particlefile,
			const char* contactfile,
			const char* progressfile,
			const char* debugfile)
{
    setBoundary(rors,           // rectangular--1 or cylindrical--0?
	    1,                  // 1-only bottom boundary;5-no top boundary;6-boxed 6 boundaries
	    0.05,               // specimen dimension
	    initboundary);      // output file, containing boundaries info
    
    deposit(total_steps,        // number of iterations
	    snapshots,          // number of snapshots
	    interval,           // record interval
	    iniptclfile,        // input file, initial particles
	    initboundary,       // input file, boundaries
	    particlefile,       // output file, resulted particles, including snapshots 
	    contactfile,        // output file, resulted contacts, including snapshots 
	    progressfile,       // output file, statistical info
	    debugfile);         // output file, debug info
}

  
void assembly::setBoundary(int   rors,
		       int   bdrynum,
		       long double dimn,
		       const char* boundaryfile)
{
    ofstream ofs(boundaryfile);
    if(!ofs) { cout<<"stream error!"<<endl; exit(-1);}
    ofs.setf(ios::scientific, ios::floatfield);
    ofs<<setw(10)<<0
       <<setw(10)<<bdrynum<<endl<<endl;

    if (rors == 1){
	if (bdrynum == 1){   // only a bottom boundary
	    ofs<<setw(10)<<1<<endl
	       <<setw(10)<<6
	       <<setw(10)<<5
	       <<setw(16)<<dimn*dimn<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<dimn/2*50
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<-dimn/2*50
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<dimn/2*50
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<-dimn/2*50
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl;
	}
	else if (bdrynum == 4){ // no top/bottom boundary
	    // boundary 1
	    ofs<<setw(10)<<1<<endl
	       <<setw(10)<<1
	       <<setw(10)<<4
	       <<setw(16)<<dimn*dimn<<endl

	       <<setw(10)<<1
	       <<setw(16)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0 
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0 
	       <<setw(16)<<0
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0     
	       <<setw(16)<<0
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl<<endl

	       // boundary 2
	       <<setw(10)<<1<<endl
	       <<setw(10)<<2
	       <<setw(10)<<4
	       <<setw(16)<<dimn*dimn<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0     
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0 
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0 
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0     
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0      
	       <<setw(16)<<0     
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl<<endl

	       // boundary 3
	       <<setw(10)<<1<<endl
	       <<setw(10)<<3
	       <<setw(10)<<4
	       <<setw(16)<<dimn*dimn<<endl

	       <<setw(10)<<1
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0 
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0 
	       <<setw(16)<<0     
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0      
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0      
	       <<setw(16)<<0     
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl<<endl

	       // boundary 4
	       <<setw(10)<<1<<endl
	       <<setw(10)<<4
	       <<setw(10)<<4
	       <<setw(16)<<dimn*dimn<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<-1
	       <<setw(16)<<0 
	       <<setw(16)<<0      
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0 
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0 
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0      
	       <<setw(16)<<0     
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl<<endl;
	}
	else if (bdrynum == 5){ // no top boundary
	    // boundary 1
	    ofs<<setw(10)<<1<<endl
	       <<setw(10)<<1
	       <<setw(10)<<4
	       <<setw(16)<<dimn*dimn<<endl

	       <<setw(10)<<1
	       <<setw(16)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0 
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0 
	       <<setw(16)<<0
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0     
	       <<setw(16)<<0
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl<<endl

	       // boundary 2
	       <<setw(10)<<1<<endl
	       <<setw(10)<<2
	       <<setw(10)<<4
	       <<setw(16)<<dimn*dimn<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0     
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0 
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0 
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0     
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0      
	       <<setw(16)<<0     
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl<<endl

	       // boundary 3
	       <<setw(10)<<1<<endl
	       <<setw(10)<<3
	       <<setw(10)<<4
	       <<setw(16)<<dimn*dimn<<endl

	       <<setw(10)<<1
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0 
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0 
	       <<setw(16)<<0     
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0      
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0      
	       <<setw(16)<<0     
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl<<endl

	       // boundary 4
	       <<setw(10)<<1<<endl
	       <<setw(10)<<4
	       <<setw(10)<<4
	       <<setw(16)<<dimn*dimn<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<-1
	       <<setw(16)<<0 
	       <<setw(16)<<0      
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0 
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0 
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0      
	       <<setw(16)<<0     
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl<<endl
		
	       // boundary 6
	       <<setw(10)<<1<<endl
	       <<setw(10)<<6
	       <<setw(10)<<5
	       <<setw(16)<<dimn*dimn<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<-dimn/2 
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl;
	}
	else if (bdrynum == 6){ // all 6 boundaries
	       // boundary 1
	    ofs<<setw(10)<<1<<endl
	       <<setw(10)<<1
	       <<setw(10)<<5
	       <<setw(16)<<dimn*dimn<<endl

	       <<setw(10)<<1
	       <<setw(16)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0     
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0     
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0       
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<1
	       <<setw(16)<<0      
	       <<setw(16)<<0     
	       <<setw(16)<<dimn/2 
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0 
	       <<setw(16)<<-1
	       <<setw(16)<<0      
	       <<setw(16)<<0      
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl<<endl

	       // boundary 2
	       <<setw(10)<<1<<endl
	       <<setw(10)<<2
	       <<setw(10)<<5
	       <<setw(16)<<dimn*dimn<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<1
	       <<setw(16)<<0     
	       <<setw(16)<<0     
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<-dimn/2 
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<1
	       <<setw(16)<<0      
	       <<setw(16)<<0     
	       <<setw(16)<<dimn/2 
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0 
	       <<setw(16)<<-1
	       <<setw(16)<<0      
	       <<setw(16)<<0      
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl<<endl

	       // boundary 3
	       <<setw(10)<<1<<endl
	       <<setw(10)<<3
	       <<setw(10)<<5
	       <<setw(16)<<dimn*dimn<<endl

	       <<setw(10)<<1
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0     
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0     
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0  
	       <<setw(16)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0       
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<1
	       <<setw(16)<<0      
	       <<setw(16)<<0     
	       <<setw(16)<<dimn/2 
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0 
	       <<setw(16)<<-1
	       <<setw(16)<<0      
	       <<setw(16)<<0      
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl<<endl

	       // boundary 4
	       <<setw(10)<<1<<endl
	       <<setw(10)<<4
	       <<setw(10)<<5
	       <<setw(16)<<dimn*dimn<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<-1
	       <<setw(16)<<0     
	       <<setw(16)<<0      
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<-1 
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<-dimn/2 
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<1
	       <<setw(16)<<0      
	       <<setw(16)<<0     
	       <<setw(16)<<dimn/2 
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0 
	       <<setw(16)<<-1
	       <<setw(16)<<0      
	       <<setw(16)<<0      
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl<<endl

	       // boundary 5
	       <<setw(10)<<1<<endl
	       <<setw(10)<<5
	       <<setw(10)<<5
	       <<setw(16)<<dimn*dimn<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0
	       <<setw(16)<<1     
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(16)<<dimn/2 
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<-1 
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<-dimn/2 
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl<<endl

	       // boundary 6
	       <<setw(10)<<1<<endl
	       <<setw(10)<<6
	       <<setw(10)<<5
	       <<setw(16)<<dimn*dimn<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0
	       <<setw(16)<<-1    
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<1
	       <<setw(16)<<0 
	       <<setw(16)<<0
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<-1 
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(16)<<-dimn/2 
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<1
	       <<setw(16)<<0
	       <<setw(16)<<0      
	       <<setw(16)<<dimn/2
	       <<setw(16)<<0      
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl

	       <<setw(10)<<1
	       <<setw(16)<<0
	       <<setw(16)<<-1
	       <<setw(16)<<0 
	       <<setw(16)<<0      
	       <<setw(16)<<-dimn/2
	       <<setw(16)<<0
	       <<setw(16)<<0
	       <<setw(10)<<0<<endl<<endl;
	}

    }
    else{
    }
    
    ofs.close();
}


void assembly::trim(int   rors,
		    const char* iniptclfile,
		    const char* inibdryfile,
		    const char* particlefile,
		    const char* boundaryfile)
{
    createSample(iniptclfile);
    createBoundary(inibdryfile);

    list<particle*>::iterator itr,itp;
    vec center;
    long double mass = 0;

    if(rors == 1) {
	if (RgdBdryNum == 1) {
	}
	else if (RgdBdryNum == 4) {
	    long double W0 = getApt(2).gety()-getApt(4).gety();
	    long double L0 = getApt(1).getx()-getApt(3).getx();
	    R.set_width(W0); 
	    R.set_length(L0); 
	    
	    for(itr=ParticleList.begin();itr!=ParticleList.end();++itr){
		center=(*itr)->getCurrPosition();
		if(fabsl(center.getx()) >= L0/2 ||
		   fabsl(center.gety()) >= W0/2 )
		{
		    itp = itr;
		    --itr;
		    delete (*itp); // release memory
		    ParticleList.erase(itp); 
		}
	    }

	}
	else if (RgdBdryNum == 5) {
	    long double W0 = getApt(2).gety()-getApt(4).gety();
	    long double L0 = getApt(1).getx()-getApt(3).getx();
	    long double H0 = -getApt(6).getz()*4;
	    R.set_width(W0); 
	    R.set_length(L0); 
	    R.set_height(H0);
	    R.set_center(0);
	    Volume = W0*L0*H0;
	    
	    for(itr=ParticleList.begin();itr!=ParticleList.end();++itr){
		center=(*itr)->getCurrPosition();
		if(fabsl(center.getx()) >= L0/2 ||
		   fabsl(center.gety()) >= W0/2 ||
		   center.getz() <= getApt(6).getz() )
		{
		    itp = itr;
		    --itr;
		    delete (*itp); // release memory
		    ParticleList.erase(itp); 
		}
	    }

	}
	else if (RgdBdryNum == 6) {
	    long double W0 = getApt(2).gety()-getApt(4).gety();
	    long double L0 = getApt(1).getx()-getApt(3).getx();
	    long double H0 = getApt(5).getz()-getApt(6).getz();
	    R.set_width(W0); 
	    R.set_length(L0); 
	    R.set_height(H0);
	    R.set_center(0);
	    Volume = W0*L0*H0;
	    
	    for(itr=ParticleList.begin();itr!=ParticleList.end();++itr){
		center=(*itr)->getCurrPosition();
		if(fabsl(center.getx()) >= L0/2 ||
		   fabsl(center.gety()) >= W0/2 ||
		   fabsl(center.getz()) >= H0/2 )
		{
		    itp = itr;
		    --itr;
		    delete (*itp); // release memory
		    ParticleList.erase(itp); 
		}
	    }
	}

	TotalNum = ParticleList.size();
	for(itr=ParticleList.begin();itr!=ParticleList.end();++itr)
	    mass += (*itr)->getMass();
	
	BulkDensity = mass/Volume;
	// Gradation =;
    }
    else {
    }
    
    printParticle(particlefile);
    printBoundary(boundaryfile);
}


void assembly::TrimPtclBdryByHeight(double height,
			    const char* iniptclfile,
			    const char* particlefile)
{
    createSample(iniptclfile);

    list<particle*>::iterator itr, itp;
    for(itr=ParticleList.begin();itr!=ParticleList.end();++itr){
	if ( (*itr)->getType() == 1 ) { // 1-fixed
	    vec center=(*itr)->getCurrPosition();
	    if(center.getz() > height)
	    {
		itp = itr;
		--itr;
		delete (*itp); // release memory
		ParticleList.erase(itp); 
	    }
	    else
		(*itr)->setType(10); // 10-ghost
	}
    }

    TotalNum = ParticleList.size();
    
    printParticle(particlefile);
}


// deposit floating particles into a container through applying gravity,
// the container can be as simple as a bottom plate
void assembly::deposit(int   total_steps,  
		       int   snapshots,
		       int   interval,
		       const char* iniptclfile,   
		       const char* inibdryfile,
		       const char* particlefile, 
		       const char* contactfile,
		       const char* progressfile, 
		       const char* debugfile)
{
    // pre_1: open streams for output.
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1); }
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average       translational    rotational       "
	       <<"kinetic        potential         total           void            sample       coordination"
	       <<"       sample           sample          sample          sample          sample          sample"
	       <<"          sample          sample          sample         sample           sample         "
	       <<" sample         sample          sample          sample          sample        minimum        wall-clock" << endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment         energy           energy          "
	       <<"energy         energy           energy          ratio          porosity         number       "
	       <<"   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume        epsilon_w       epsilon_l       epsilon_h       "
	       <<"epsilon-v      time_step          time" << endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1); }
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from existing files.
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries.

    // pre_3: define variables used in iterations.
    long double l13, l24, l56;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    long double void_ratio=0;
    long double bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    // iterations starting ...
    g_iteration=0; 
    gettimeofday(&timew1,NULL);
    do
    {
	// 1. create possible boundary particles and contacts between particles.
        findContact();
        findParticleOnBoundary();

	// 2. set particles' forces/moments as zero before each re-calculation,
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles.
	internalForce(avgNormal, avgTangt);

	// 4. calculate boundary forces/moments and apply them to particles.
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();

	// 6. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx(); Volume=l13*l24*l56;
	void_ratio=Volume/getParticleVolume()-1;

	// 7. (1) output particles and contacts information as snapshots.
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
     	    cout << stepsfp;

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	    time(&timeStamp);
	    cout << "  " << stepsfp << "  " << ctime(&timeStamp);
	}

	// 7. (2) output stress and strain info.
	if (g_iteration % interval == 0) {
	    gettimeofday(&timew2,NULL);
	    long double t1=getTransEnergy();
	    long double t2=getRotatEnergy();
	    long double t3=getPotenEnergy(-0.025);
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<t1
		       <<setw(16)<<t2
		       <<setw(16)<<(t1+t2)
		       <<setw(16)<<t3
		       <<setw(16)<<(t1+t2+t3)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/TotalNum
	               <<setw(16)<<getMinTimeStep()
		       <<setw(16)<<timediffsec(timew1,timew2)
		       <<endl;

	    /*
	    g_debuginf<<setw(10)<<g_iteration
		      <<setw(16)<<bdry_penetr[1]
		      <<setw(16)<<bdry_penetr[2]
		      <<setw(16)<<bdry_penetr[3]
		      <<setw(16)<<bdry_penetr[4]
		      <<setw(16)<<bdry_penetr[6]
		      <<setw(16)<<bdry_cntnum[1]
		      <<setw(16)<<bdry_cntnum[2]
		      <<setw(16)<<bdry_cntnum[3]
		      <<setw(16)<<bdry_cntnum[4]
		      <<setw(16)<<bdry_cntnum[6]
		      <<endl;
	    */

	}

	// 8. loop break conditions.

    } while (++g_iteration < total_steps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);
    cout << stepsfp;

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    cout << "  " << stepsfp << "  " << ctime(&timeStamp);

    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


// actual deposit function for the case of fixed particle boundaries
void assembly::deposit_p(int   total_steps,  
			 int   snapshots,
			 int   interval,
			 long double dimn,
			 long double rsize,
			 const char* iniptclfile,   
			 const char* particlefile, 
			 const char* contactfile,
			 const char* progressfile, 
			 const char* debugfile)
{
    // pre_1: open streams for output.
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1); }
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"deposit..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average       translational    rotational       "
	       <<"kinetic        potential         total           void            sample       coordination"
	       <<"       sample           sample          sample          sample          sample          sample"
	       <<"          sample          sample          sample         sample           sample         "
	       <<" sample          sample          sample          sample          sample"<<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"         omga            force           moment         energy           energy          "
	       <<"energy         energy            energy          ratio          porosity         number       "
	       <<"   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       "
	       <<"epsilon-v"<<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1); }
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from existing files.
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 

    // pre_3: define variables used in iterations.
    long double l13, l24, l56;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    long double void_ratio=0;

    // iterations starting ...
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles.
	findContact();

	// 2. set particles' forces/moments as zero before each re-calculation,
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles.
	internalForce(avgNormal, avgTangt);

	// 4. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();

	// 5. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - (-dimn/2);
	l24=dimn*rsize;
	l13=dimn*rsize;
	Volume=l13*l24*l56;
	void_ratio=Volume/getParticleVolume()-1;

	// 6. (1) output particles and contacts information as snapshots.
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 6. (2) output statistics info.
	if (g_iteration % interval == 0) {
	    long double t1=getTransEnergy();
	    long double t2=getRotatEnergy();
	    long double t3=getPotenEnergy(-0.025);
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<t1
		       <<setw(16)<<t2
		       <<setw(16)<<(t1+t2)
		       <<setw(16)<<t3
		       <<setw(16)<<(t1+t2+t3)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*getActualCntctNum()/TotalNum
		       <<endl;
	}

	// 7. loop break conditions.


    } while (++g_iteration < total_steps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


// squeeze paticles inside a container by moving the boundaries
void assembly::squeeze(int   total_steps,  
		       int   init_steps,
		       int   snapshots,
		       int   interval,
		       int   flag,
		       const char* iniptclfile,   
		       const char* inibdryfile,
		       const char* particlefile, 
		       const char* boundaryfile,
		       const char* contactfile,
		       const char* progressfile, 
		       const char* debugfile)
{
    // pre_1: open streams for output.
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1); }
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"deposit..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average       translational    rotational       "
	       <<"kinetic        potential         total           void            sample       coordination"
	       <<"       sample           sample          sample          sample          sample          sample"
	       <<"          sample          sample          sample         sample           sample         "
	       <<" sample          sample          sample          sample          sample"<<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"         omga            force           moment         energy           energy          "
	       <<"energy         energy            energy          ratio          porosity         number       "
	       <<"   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       "
	       <<"epsilon-v"<<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1); }
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from existing files.
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries.

    // pre_3: define variables used in iterations.
    long double l13, l24, l56;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];

    int         mid[2]={1,3};    // boundary 1 and 3
    UPDATECTL   midctl[2];
    long double void_ratio=0;
    long double bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    // iterations starting ...
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles.
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces/moments as zero before each re-calculation,
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles.
	internalForce(avgNormal, avgTangt);

	// 4. calculate boundary forces/moments and apply them to particles.
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment.
	updateParticle();

	// 6. calculate sample void ratio.
	l56=getTopFreeParticlePosition().getz() -getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx(); Volume=l13*l24*l56;
	void_ratio=Volume/getParticleVolume()-1;

	// displacement control
	if (g_iteration > init_steps) {
	    if (flag==1) // loosen, totally remove the wall
		midctl[0].tran=vec(TIMESTEP*1.0e+0*flag,0,0);
	    else         // squeeze
		midctl[0].tran=vec(TIMESTEP*5.0e-3*flag,0,0);
	    updateRB(mid,midctl,2);
	}

	// 7. (1) output particles and contacts information as snapshots.
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info.
	if (g_iteration % interval == 0) {
	    long double t1=getTransEnergy();
	    long double t2=getRotatEnergy();
	    long double t3=getPotenEnergy(-0.025);
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<t1
		       <<setw(16)<<t2
		       <<setw(16)<<(t1+t2)
		       <<setw(16)<<t3
		       <<setw(16)<<(t1+t2+t3)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    g_debuginf<<setw(10)<<g_iteration
		      <<setw(16)<<bdry_penetr[1]
		      <<setw(16)<<bdry_penetr[2]
		      <<setw(16)<<bdry_penetr[3]
		      <<setw(16)<<bdry_penetr[4]
		      <<setw(16)<<bdry_penetr[6]
		      <<setw(16)<<bdry_cntnum[1]
		      <<setw(16)<<bdry_cntnum[2]
		      <<setw(16)<<bdry_cntnum[3]
		      <<setw(16)<<bdry_cntnum[4]
		      <<setw(16)<<bdry_cntnum[6]
		      <<endl;

	}

	// 8. loop break conditions.

    } while (++g_iteration < total_steps);
    
    // post_1. store the final snapshot of particles & contacts.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);

    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


// Isotropically compress floating particles to a specific confining pressure, which is usually a low
// value in order to create an intial status. Force boundaries are used. This process may be not 
// physically true.
void assembly::isotropic(int   total_steps,
			 int   snapshots, 
			 int   interval,
			 long double sigma,			  
			 const char* iniptclfile,   
			 const char* inibdryfile,
			 const char* particlefile, 
			 const char* boundaryfile,
			 const char* contactfile,  
			 const char* progressfile,
			 const char* balancedfile, 
			 const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1);}
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"isotropic..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(ios::scientific, ios::floatfield);
    balancedinf<<"isotropic..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3: define variables used in iterations
    long double W0 = getApt(2).gety()-getApt(4).gety();
    long double L0 = getApt(1).getx()-getApt(3).getx();
    long double H0 = getApt(5).getz()-getApt(6).getz();
    long double l13, l24, l56, min_area, mid_area, max_area;
    long double sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    long double epsilon_w, epsilon_l, epsilon_h;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         mid[2]={1,3};    // boundary 1 and 3
    int         max[2]={2,4};    // boundary 2 and 4
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   midctl[2];
    UPDATECTL   maxctl[2];
    UPDATECTL   minctl[2];
    long double void_ratio=0;
    long double bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    // iterations start here...
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabsl(getNormalForce(2))/max_area; sigma1_2=vfabsl(getNormalForce(4))/max_area;
	sigma2_1=vfabsl(getNormalForce(1))/mid_area; sigma2_2=vfabsl(getNormalForce(3))/mid_area;
	sigma3_1=vfabsl(getNormalForce(5))/min_area; sigma3_2=vfabsl(getNormalForce(6))/min_area;
	void_ratio=Volume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	if (sigma2_1<sigma)
	    midctl[0].tran=vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma)
	    midctl[1].tran=vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma)
	    maxctl[0].tran=vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma)
	    maxctl[1].tran=vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=vec(0,-TIMESTEP*RELEASE_RATE,0);
	
	updateRB(min,minctl,2);
	updateRB(mid,midctl,2);
	updateRB(max,maxctl,2);
	updateRB6();
	
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0 ){
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    g_debuginf<<setw(10)<<g_iteration
		      <<setw(16)<<bdry_penetr[1]
		      <<setw(16)<<bdry_penetr[2]
		      <<setw(16)<<bdry_penetr[3]
		      <<setw(16)<<bdry_penetr[4]
		      <<setw(16)<<bdry_penetr[5]
		      <<setw(16)<<bdry_penetr[6]
		      <<setw(16)<<bdry_cntnum[1]
		      <<setw(16)<<bdry_cntnum[2]
		      <<setw(16)<<bdry_cntnum[3]
		      <<setw(16)<<bdry_cntnum[4]
		      <<setw(16)<<bdry_cntnum[5]
		      <<setw(16)<<bdry_cntnum[6]
		      <<endl;
	}

	// 8. loop break condition
	if (   fabsl(sigma1_1-sigma)/sigma < STRESS_ERROR && fabsl(sigma1_2-sigma)/sigma < STRESS_ERROR
	    && fabsl(sigma2_1-sigma)/sigma < STRESS_ERROR && fabsl(sigma2_2-sigma)/sigma < STRESS_ERROR
	    && fabsl(sigma3_1-sigma)/sigma < STRESS_ERROR && fabsl(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()    
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()    
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    break;
	}

    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile);  strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_a. This function
// increases confining pressure step by step to sigma_b, making it possible to find equilibrium 
// state where particle pressure equals confining pressure. Force boundaries are used
void assembly::isotropic(int   total_steps,
			 int   snapshots, 
			 int   interval,
			 long double sigma_a,
			 long double sigma_b,
			 int   sigma_division,
			 const char* iniptclfile,   
			 const char* inibdryfile,
			 const char* particlefile, 
			 const char* boundaryfile,
			 const char* contactfile,  
			 const char* progressfile,
			 const char* balancedfile, 
			 const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1);}
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"isotropic..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(ios::scientific, ios::floatfield);
    balancedinf<<"isotropic..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3: define variables used in iterations
    long double W0 = getApt(2).gety()-getApt(4).gety();
    long double L0 = getApt(1).getx()-getApt(3).getx();
    long double H0 = getApt(5).getz()-getApt(6).getz();
    long double l13, l24, l56, min_area, mid_area, max_area;
    long double sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    long double epsilon_w, epsilon_l, epsilon_h;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         mid[2]={1,3};    // boundary 1 and 3
    int         max[2]={2,4};    // boundary 2 and 4
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   midctl[2];
    UPDATECTL   maxctl[2];
    UPDATECTL   minctl[2];
    long double void_ratio=0;
    long double bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    long double sigma=sigma_a;
    long double sigma_inc=(sigma_b-sigma_a)/sigma_division;

    // iterations start here...
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabsl(getNormalForce(2))/max_area; sigma1_2=vfabsl(getNormalForce(4))/max_area;
	sigma2_1=vfabsl(getNormalForce(1))/mid_area; sigma2_2=vfabsl(getNormalForce(3))/mid_area;
	sigma3_1=vfabsl(getNormalForce(5))/min_area; sigma3_2=vfabsl(getNormalForce(6))/min_area;
	void_ratio=Volume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	if (sigma2_1<sigma)
	    midctl[0].tran=vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma)
	    midctl[1].tran=vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma)
	    maxctl[0].tran=vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma)
	    maxctl[1].tran=vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=vec(0,-TIMESTEP*RELEASE_RATE,0);
	
	updateRB(min,minctl,2);
	updateRB(mid,midctl,2);
	updateRB(max,maxctl,2);
	updateRB6();
	
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0 ){
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    g_debuginf<<setw(10)<<g_iteration
		      <<setw(16)<<bdry_penetr[1]
		      <<setw(16)<<bdry_penetr[2]
		      <<setw(16)<<bdry_penetr[3]
		      <<setw(16)<<bdry_penetr[4]
		      <<setw(16)<<bdry_penetr[5]
		      <<setw(16)<<bdry_penetr[6]
		      <<setw(16)<<bdry_cntnum[1]
		      <<setw(16)<<bdry_cntnum[2]
		      <<setw(16)<<bdry_cntnum[3]
		      <<setw(16)<<bdry_cntnum[4]
		      <<setw(16)<<bdry_cntnum[5]
		      <<setw(16)<<bdry_cntnum[6]
		      <<endl;
	}

	// 8. find the balanced status and increase confining pressure
	if (   fabsl(sigma1_1-sigma)/sigma < STRESS_ERROR && fabsl(sigma1_2-sigma)/sigma < STRESS_ERROR
	    && fabsl(sigma2_1-sigma)/sigma < STRESS_ERROR && fabsl(sigma2_2-sigma)/sigma < STRESS_ERROR
	    && fabsl(sigma3_1-sigma)/sigma < STRESS_ERROR && fabsl(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()    
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    sigma += sigma_inc;
	}

	// 9. loop break condition
	if (   fabsl(sigma1_1-sigma_b)/sigma_b < STRESS_ERROR && fabsl(sigma1_2-sigma_b)/sigma_b < STRESS_ERROR
	    && fabsl(sigma2_1-sigma_b)/sigma_b < STRESS_ERROR && fabsl(sigma2_2-sigma_b)/sigma_b < STRESS_ERROR
	    && fabsl(sigma3_1-sigma_b)/sigma_b < STRESS_ERROR && fabsl(sigma3_2-sigma_b)/sigma_b < STRESS_ERROR ) {
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()    
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    break;
	}
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


// loading-unloading-reloading of isotropic compression
// the stress path is defined by sigma_points and sigma_values[]
void assembly::isotropic(int   total_steps,  
			 int   snapshots, 
			 int   interval,
			 int   sigma_points,  
			 long double sigma_values[],  
			 int   sigma_division,	  
			 const char* iniptclfile,  
			 const char* inibdryfile,
			 const char* particlefile, 
			 const char* boundaryfile,
			 const char* contactfile,  
			 const char* progressfile,
			 const char* balancedfile, 
			 const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1);}
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"isotropic..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(ios::scientific, ios::floatfield);
    balancedinf<<"isotropic..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3: define variables used in iterations
    long double W0 = getApt(2).gety()-getApt(4).gety();
    long double L0 = getApt(1).getx()-getApt(3).getx();
    long double H0 = getApt(5).getz()-getApt(6).getz();
    long double l13, l24, l56, min_area, mid_area, max_area;
    long double sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    long double epsilon_w, epsilon_l, epsilon_h;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         mid[2]={1,3};    // boundary 1 and 3
    int         max[2]={2,4};    // boundary 2 and 4
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   midctl[2];
    UPDATECTL   maxctl[2];
    UPDATECTL   minctl[2];
    long double void_ratio=0;
    long double bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    int  i=0;
    long double sigma=sigma_values[i];
    long double sigma_inc=(sigma_values[i+1]-sigma_values[i])/sigma_division;
    long double sigma_b=sigma_values[sigma_points-1];

    // iterations start here...
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();
	
	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabsl(getNormalForce(2))/max_area; sigma1_2=vfabsl(getNormalForce(4))/max_area;
	sigma2_1=vfabsl(getNormalForce(1))/mid_area; sigma2_2=vfabsl(getNormalForce(3))/mid_area;
	sigma3_1=vfabsl(getNormalForce(5))/min_area; sigma3_2=vfabsl(getNormalForce(6))/min_area;
	void_ratio=Volume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	if (sigma2_1<sigma)
	    midctl[0].tran=vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma)
	    midctl[1].tran=vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma)
	    maxctl[0].tran=vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma)
	    maxctl[1].tran=vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=vec(0,-TIMESTEP*RELEASE_RATE,0);
	
	updateRB(min,minctl,2);
	updateRB(mid,midctl,2);
	updateRB(max,maxctl,2);
	updateRB6();
	
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0 ){
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    g_debuginf<<setw(10)<<g_iteration
		      <<setw(16)<<bdry_penetr[1]
		      <<setw(16)<<bdry_penetr[2]
		      <<setw(16)<<bdry_penetr[3]
		      <<setw(16)<<bdry_penetr[4]
		      <<setw(16)<<bdry_penetr[5]
		      <<setw(16)<<bdry_penetr[6]
		      <<setw(16)<<bdry_cntnum[1]
		      <<setw(16)<<bdry_cntnum[2]
		      <<setw(16)<<bdry_cntnum[3]
		      <<setw(16)<<bdry_cntnum[4]
		      <<setw(16)<<bdry_cntnum[5]
		      <<setw(16)<<bdry_cntnum[6]
		      <<endl;
	}

	// 8. find the balanced status and increase confining pressure
	if (   fabsl(sigma1_1-sigma)/sigma < STRESS_ERROR && fabsl(sigma1_2-sigma)/sigma < STRESS_ERROR
	    && fabsl(sigma2_1-sigma)/sigma < STRESS_ERROR && fabsl(sigma2_2-sigma)/sigma < STRESS_ERROR
	    && fabsl(sigma3_1-sigma)/sigma < STRESS_ERROR && fabsl(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()    
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    sigma += sigma_inc;
	    if (sigma==sigma_values[i+1]) {
		i++;
		sigma=sigma_values[i];
		sigma_inc=(sigma_values[i+1]-sigma_values[i])/sigma_division;
	    }

	}

	// 9. loop break condition
	if (   fabsl(sigma1_1-sigma_b)/sigma_b < STRESS_ERROR && fabsl(sigma1_2-sigma_b)/sigma_b < STRESS_ERROR
	    && fabsl(sigma2_1-sigma_b)/sigma_b < STRESS_ERROR && fabsl(sigma2_2-sigma_b)/sigma_b < STRESS_ERROR
	    && fabsl(sigma3_1-sigma_b)/sigma_b < STRESS_ERROR && fabsl(sigma3_2-sigma_b)/sigma_b < STRESS_ERROR ) {
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()    
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    break;
	}
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_3. This function
// increases vertical pressure step by step to sigma_1, thus making it possible to find out
// balanced status where top & bottom particle pressure equals major principle stress. 
// Side boundaries are fixed, top and bottom plates are force-controlled.
void assembly::odometer(int   total_steps,  
			int   snapshots, 
			int   interval,
			long double sigma_3,     
			long double sigma_1,    
			int   sigma_division,			  
			const char* iniptclfile,  
			const char* inibdryfile,
			const char* particlefile,
			const char* boundaryfile,
			const char* contactfile,  
			const char* progressfile,
			const char* balancedfile, 
			const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1);}
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"odometer..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(ios::scientific, ios::floatfield);
    balancedinf<<"odometer..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries
 
    // pre_3. define variables used in iterations
    long double W0 = getApt(2).gety()-getApt(4).gety();
    long double L0 = getApt(1).getx()-getApt(3).getx();
    long double H0 = getApt(5).getz()-getApt(6).getz();
    long double l13, l24, l56, min_area, mid_area, max_area;
    long double sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    long double epsilon_w, epsilon_l, epsilon_h;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];

    int min[2]={5,6};    // minimum stress acting on boundary 5 and 6
    UPDATECTL minctl[2];
    long double void_ratio=0;
    long double bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    long double sigma=sigma_3;
    long double sigma_inc=(sigma_1-sigma_3)/sigma_division;

    // iterations start here...
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();
	
	// 2. set particles' forces and moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces and moments
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/displacement based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabsl(getNormalForce(2))/max_area; sigma1_2=vfabsl(getNormalForce(4))/max_area;
	sigma2_1=vfabsl(getNormalForce(1))/mid_area; sigma2_2=vfabsl(getNormalForce(3))/mid_area;
	sigma3_1=vfabsl(getNormalForce(5))/min_area; sigma3_2=vfabsl(getNormalForce(6))/min_area;
	void_ratio=Volume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	updateRB(min,minctl,2);
	updateRB6();

	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0){
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	}

	// 8. find balanced status of odometer compression
	if (fabsl(sigma3_1-sigma)/sigma < STRESS_ERROR && fabsl(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()    
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    sigma += sigma_inc;
	}

	// 9. loop break condition
	if (fabsl(sigma3_1-sigma_1)/sigma_1 < STRESS_ERROR && fabsl(sigma3_2-sigma_1)/sigma_1 < STRESS_ERROR) {
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()    
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    break;
	}
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_3. This function
// increases vertical pressure step by step to sigma_1, thus making it possible to find out
// balanced status where top & bottom particle pressure equals major principle stress. 
// Side boundaries are fixed, top and bottom plates are force-controlled. Unloading path is
// applied.
void assembly::odometer(int   total_steps,  
			int   snapshots, 
			int   interval,
			int   sigma_points,  
			long double sigma_values[],  
			int   sigma_division,			  
			const char* iniptclfile,  
			const char* inibdryfile,
			const char* particlefile, 
			const char* boundaryfile,
			const char* contactfile,  
			const char* progressfile,
			const char* balancedfile, 
			const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1);}
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"odometer..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(ios::scientific, ios::floatfield);
    balancedinf<<"odometer..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries
 
    // pre_3. define variables used in iterations
    long double W0 = getApt(2).gety()-getApt(4).gety();
    long double L0 = getApt(1).getx()-getApt(3).getx();
    long double H0 = getApt(5).getz()-getApt(6).getz();
    long double l13, l24, l56, min_area, mid_area, max_area;
    long double sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    long double epsilon_w, epsilon_l, epsilon_h;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];

    int min[2]={5,6};    // minimum stress acting on boundary 5 and 6
    UPDATECTL minctl[2];
    long double void_ratio=0;
    long double bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }


    int  i=0;
    long double sigma=sigma_values[i];
    long double sigma_inc=(sigma_values[i+1]-sigma_values[i])/sigma_division;
    long double sigma_b=sigma_values[sigma_points-1];

    // iterations start here...
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces and moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces and moments
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/displacement based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabsl(getNormalForce(2))/max_area; sigma1_2=vfabsl(getNormalForce(4))/max_area;
	sigma2_1=vfabsl(getNormalForce(1))/mid_area; sigma2_2=vfabsl(getNormalForce(3))/mid_area;
	sigma3_1=vfabsl(getNormalForce(5))/min_area; sigma3_2=vfabsl(getNormalForce(6))/min_area;
	void_ratio=Volume/getParticleVolume()-1;

	if (sigma3_1<sigma)
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma)
	    minctl[1].tran=vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	updateRB(min,minctl,2);
	updateRB6();

	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0){
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	}

	// 8. find balanced status of odometer compression
	if (fabsl(sigma3_1-sigma)/sigma < STRESS_ERROR && fabsl(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()    
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    sigma += sigma_inc;
	    if (sigma==sigma_values[i+1]) {
		i++;
		sigma=sigma_values[i];
		sigma_inc=(sigma_values[i+1]-sigma_values[i])/sigma_division;
	    }
	}

	// 9. loop break condition
	if (fabsl(sigma3_1-sigma_b)/sigma_b < STRESS_ERROR && fabsl(sigma3_2-sigma_b)/sigma_b < STRESS_ERROR) {
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()    
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    break;
	}
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


void assembly::unconfined(int   total_steps,  
			  int   snapshots,	
			  int   interval,
			  const char* iniptclfile,  
			  const char* inibdryfile,
			  const char* particlefile,
			  const char* contactfile,  
			  const char* progressfile,
			  const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.  
    progressinf.open(progressfile);
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1);}
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"unconfined..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample"
	       <<"          sample          sample          sample"<<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       "
	       <<"epsilon-v"<<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries
 
    // pre_3. define variables used in iterations
    long double sigma3_1, sigma3_2;
    int    stepsnum=0;
    char   stepsstr[4];
    char   stepsfp[50];
    long double avgNormal=0;
    long double avgTangt=0;
    int    min[2]={5,6};    //  boundary 5 and 6
    UPDATECTL minctl[2];

    // iterations start here...
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces and moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces and moments
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces
	rigidBoundaryForce();

	// 5. update particles' velocity/omga/displacement based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	sigma3_1=vfabsl(getNormalForce(5))/getArea(5); sigma3_2=vfabsl(getNormalForce(6))/getArea(6);
	minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	minctl[1].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);
	updateRB(min,minctl,2);

	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output progress info.
	if (g_iteration % interval == 0)
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<0
		       <<setw(16)<<0<<setw(16)<<0
		       <<setw(16)<<0<<setw(16)<<0
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<endl;
/*
	// 8. loop break condition
	if (getAverageForce() < 1.0) {
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()    
		       <<setw(16)<<getAverageMoment()<<endl;
	    break;
	}
*/
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


// This function initializes triaxial sample to a certain confining pressure.
void assembly::triaxialPtclBdryIni(int   total_steps,  
				   int   snapshots, 
				   int   interval,
				   double sigma,
				   const char* iniptclfile, 
				   const char* inibdryfile,
				   const char* particlefile,
				   const char* boundaryfile,
				   const char* contactfile, 
				   const char* progressfile,
				   const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1);}
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"triaxial..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    long double H0 = getApt(5).getz()-getApt(6).getz();
    long double l56= 0;
    long double sigma3_1, sigma3_2;
    long double epsilon_h;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   minctl[2];

    // iterations start here...
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	rigidBoundaryForce();

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	sigma3_1=vfabsl(getNormalForce(5))/2.5e-3; sigma3_2=vfabsl(getNormalForce(6))/2.5e-3;

	// force control
	if (sigma3_1 < sigma)
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);

	if (sigma3_2 < sigma)
	    minctl[1].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);

	updateRB(min,minctl,2);
	
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	l56=getApt(5).getz()-getApt(6).getz();
	epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0 ){
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<0<<setw(16)<<0
		       <<setw(16)<<0<<setw(16)<<0
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()
		       <<setw(16)<<0<<setw(16)<<0<<setw(16)<<l56
		       <<setw(16)<<0
		       <<setw(16)<<0
		       <<setw(16)<<0
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<0
		       <<setw(16)<<0
		       <<setw(16)<<2.0*getActualCntctNum()/TotalNum
		       <<endl;

	}

	// 9. loop break condition: through displacement control mechanism
	if (   fabsl(sigma3_1-sigma)/sigma < STRESS_ERROR && fabsl(sigma3_2-sigma)/sigma < STRESS_ERROR )
	       break;
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


// This function performs triaxial compression test.
// Displacement boundaries are used in axial direction.
void assembly::triaxialPtclBdry(int   total_steps,  
				int   snapshots, 
				int   interval,
				const char* iniptclfile, 
				const char* inibdryfile,
				const char* particlefile,
				const char* boundaryfile,
				const char* contactfile, 
				const char* progressfile,
				const char* balancedfile,
				const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1);}
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"triaxial..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(ios::scientific, ios::floatfield);
    balancedinf<<"triaxial..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    long double H0 = getApt(5).getz()-getApt(6).getz();
    long double l56= 0;
    long double sigma3_1, sigma3_2;
    long double epsilon_h;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   minctl[2];

    // iterations start here...
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	rigidBoundaryForce();

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	sigma3_1=vfabsl(getNormalForce(5))/2.5e-3; sigma3_2=vfabsl(getNormalForce(6))/2.5e-3;

	// displacement control
	if(g_iteration < 100001) {
	minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	minctl[1].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);

	updateRB(min,minctl,2);
	}
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	l56=getApt(5).getz()-getApt(6).getz();
	epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0 ){
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<0<<setw(16)<<0
		       <<setw(16)<<0<<setw(16)<<0
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()
		       <<setw(16)<<0<<setw(16)<<0<<setw(16)<<l56
		       <<setw(16)<<0
		       <<setw(16)<<0
		       <<setw(16)<<0
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<0
		       <<setw(16)<<0
		       <<setw(16)<<2.0*getActualCntctNum()/TotalNum
		       <<endl;

	}

/* Most time it is balanced, so use progressinf instead.
	// 8. find the balanced status and increase confining pressure
	if (   fabsl(sigma1_1-sigma_a)/sigma_a < STRESS_ERROR && fabsl(sigma1_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabsl(sigma2_1-sigma_a)/sigma_a < STRESS_ERROR && fabsl(sigma2_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabsl(sigma3_1-sigma3_2)/(sigma3_1+sigma3_2)*2<=0.05) {
	    balancedinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()    
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)<<endl;
	}
*/
	// 9. loop break condition: through displacement control mechanism
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_a. This function
// performs triaxial compression test. Displacement boundaries are used in axial direction.
void assembly::triaxial(int   total_steps,  
			int   snapshots, 
			int   interval,
			long double sigma_a,	  
			const char* iniptclfile, 
			const char* inibdryfile,
			const char* particlefile,
			const char* boundaryfile,
			const char* contactfile, 
			const char* progressfile,
			const char* balancedfile,
			const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1);}
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<"       minimum        wall-clock" << endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<"        time_step          time" << endl;

    ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(ios::scientific, ios::floatfield);
    balancedinf<<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<"       minimum        wall-clock"<< endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<"        time_step          time" << endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    long double W0 = getApt(2).gety()-getApt(4).gety();
    long double L0 = getApt(1).getx()-getApt(3).getx();
    long double H0 = getApt(5).getz()-getApt(6).getz();
    long double l13, l24, l56, min_area, mid_area, max_area;
    long double sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    long double epsilon_w, epsilon_l, epsilon_h;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         mid[2]={1,3};    // boundary 1 and 3
    int         max[2]={2,4};    // boundary 2 and 4
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   midctl[2];
    UPDATECTL   maxctl[2];
    UPDATECTL   minctl[2];
    long double void_ratio=0;
    long double bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    // iterations start here...
    g_iteration=0;
    gettimeofday(&timew1,NULL);
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabsl(getNormalForce(2))/max_area; sigma1_2=vfabsl(getNormalForce(4))/max_area;
	sigma2_1=vfabsl(getNormalForce(1))/mid_area; sigma2_2=vfabsl(getNormalForce(3))/mid_area;
	sigma3_1=vfabsl(getNormalForce(5))/min_area; sigma3_2=vfabsl(getNormalForce(6))/min_area;
	void_ratio=Volume/getParticleVolume()-1;

	// displacement control
	minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	minctl[1].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);

	// force control
	if (sigma2_1<sigma_a)
	    midctl[0].tran=vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma_a)
	    midctl[1].tran=vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma_a)
	    maxctl[0].tran=vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma_a)
	    maxctl[1].tran=vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=vec(0,-TIMESTEP*RELEASE_RATE,0);
	
	updateRB(min,minctl,2);
	updateRB(mid,midctl,2);
	updateRB(max,maxctl,2);
	updateRB6();
	
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    cout << stepsfp;
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	    time(&timeStamp);
	    cout << "  " << stepsfp << "  " << ctime(&timeStamp);
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0 ){
	    gettimeofday(&timew2,NULL);
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
	               <<setw(16)<<getMinTimeStep()
		       <<setw(16)<<timediffsec(timew1,timew2)
		       <<endl;
	    g_debuginf<<setw(10)<<g_iteration
		      <<setw(16)<<getTransEnergy()
		      <<setw(16)<<getRotatEnergy()
		      <<setw(16)<<bdry_penetr[1]
		      <<setw(16)<<bdry_penetr[2]
		      <<setw(16)<<bdry_penetr[3]
		      <<setw(16)<<bdry_penetr[4]
		      <<setw(16)<<bdry_penetr[5]
		      <<setw(16)<<bdry_penetr[6]
		      <<setw(16)<<bdry_cntnum[1]
		      <<setw(16)<<bdry_cntnum[2]
		      <<setw(16)<<bdry_cntnum[3]
		      <<setw(16)<<bdry_cntnum[4]
		      <<setw(16)<<bdry_cntnum[5]
		      <<setw(16)<<bdry_cntnum[6]
		      <<endl;
	}

	// Most time it is balanced, so use progressinf instead.
	// 8. find the balanced status and increase confining pressure
	if (   fabsl(sigma1_1-sigma_a)/sigma_a < STRESS_ERROR && fabsl(sigma1_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabsl(sigma2_1-sigma_a)/sigma_a < STRESS_ERROR && fabsl(sigma2_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabsl(sigma3_1-sigma3_2)/(sigma3_1+sigma3_2)*2<=0.05) {
	    balancedinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()    
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
	               <<setw(16)<<getMinTimeStep()
		       <<setw(16)<<timediffsec(timew1,timew2)
		       <<endl;
	}

	// 9. loop break condition: through displacement control mechanism
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);
    cout << stepsfp;

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    cout << "  " << stepsfp << "  " << ctime(&timeStamp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_a. This function
// performs triaxial compression test with unloading. Displacement boundaries are used in 
// axial direction.
void assembly::triaxial(int   total_steps,  
			int   unload_step,
			int   snapshots, 
			int   interval,
			long double sigma_a,	  
			const char* iniptclfile,  
			const char* inibdryfile,
			const char* particlefile,
			const char* boundaryfile,
			const char* contactfile,
			const char* progressfile,
			const char* balancedfile,
			const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1);}
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"triaxial..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(ios::scientific, ios::floatfield);
    balancedinf<<"triaxial..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    long double W0 = getApt(2).gety()-getApt(4).gety();
    long double L0 = getApt(1).getx()-getApt(3).getx();
    long double H0 = getApt(5).getz()-getApt(6).getz();
    long double l13, l24, l56, min_area, mid_area, max_area;
    long double sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    long double epsilon_w, epsilon_l, epsilon_h;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         mid[2]={1,3};    // boundary 1 and 3
    int         max[2]={2,4};    // boundary 2 and 4
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   midctl[2];
    UPDATECTL   maxctl[2];
    UPDATECTL   minctl[2];
    long double void_ratio=0;
    long double bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    // iterations start here...
    bool reload=false;
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();
	
	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabsl(getNormalForce(2))/max_area; sigma1_2=vfabsl(getNormalForce(4))/max_area;
	sigma2_1=vfabsl(getNormalForce(1))/mid_area; sigma2_2=vfabsl(getNormalForce(3))/mid_area;
	sigma3_1=vfabsl(getNormalForce(5))/min_area; sigma3_2=vfabsl(getNormalForce(6))/min_area;
	void_ratio=Volume/getParticleVolume()-1;

	// displacement control
	if (g_iteration <= unload_step){ //loading
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	    minctl[1].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);
	}
	else { 
	    if (reload==false) { // unloading
		if (fabsl(sigma3_1-sigma_a)/sigma_a > STRESS_ERROR && 
		    fabsl(sigma3_2-sigma_a)/sigma_a > STRESS_ERROR){
		    minctl[0].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);
		    minctl[1].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
		}
		else  // reloading
		    reload=true;
	    }
	    else {
		minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
		minctl[1].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);
	    }
	}
	
	// force control
	if (sigma2_1<sigma_a)
	    midctl[0].tran=vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma_a)
	    midctl[1].tran=vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma_a)
	    maxctl[0].tran=vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma_a)
	    maxctl[1].tran=vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=vec(0,-TIMESTEP*RELEASE_RATE,0);
	
	updateRB(min,minctl,2);
	updateRB(mid,midctl,2);
	updateRB(max,maxctl,2);
	updateRB6();
	
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0 ){
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    g_debuginf<<setw(10)<<g_iteration
		      <<setw(16)<<bdry_penetr[1]
		      <<setw(16)<<bdry_penetr[2]
		      <<setw(16)<<bdry_penetr[3]
		      <<setw(16)<<bdry_penetr[4]
		      <<setw(16)<<bdry_penetr[5]
		      <<setw(16)<<bdry_penetr[6]
		      <<setw(16)<<bdry_cntnum[1]
		      <<setw(16)<<bdry_cntnum[2]
		      <<setw(16)<<bdry_cntnum[3]
		      <<setw(16)<<bdry_cntnum[4]
		      <<setw(16)<<bdry_cntnum[5]
		      <<setw(16)<<bdry_cntnum[6]
		      <<endl;
	}

/*
	// 8. find the balanced status and increase confining pressure
	if (   fabsl(sigma1_1-sigma_a)/sigma_a < STRESS_ERROR && fabsl(sigma1_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabsl(sigma2_1-sigma_a)/sigma_a < STRESS_ERROR && fabsl(sigma2_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabsl(sigma3_1-sigma3_2)/(sigma3_1+sigma3_2)*2<=0.05) {
	    balancedinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()    
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)<<endl;
	}
*/
	// 9. loop break condition: through displacement control mechanism
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


// The specimen has been deposited with gravitation within boundaries composed of particles.
// A rectangular pile is then drived into the particles using displacement control.
void assembly::rectPile_Disp(int   total_steps,  
			     int   snapshots, 
			     int   interval,
			     const char* iniptclfile,  
			     const char* inibdryfile,
			     const char* particlefile, 
			     const char* boundaryfile,
			     const char* contactfile,  
			     const char* progressfile,
			     const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1); }
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"pile penetrate..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average       translational    rotational       "
	       <<"kinetic        potential        total           sample           sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample"
	       <<"          sample          sample          sample"<<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"         omga            force           moment         energy           energy          "
	       <<"energy         energy          energy          density         "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       "
	       <<"epsilon-v"<<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(ios::scientific, ios::floatfield);
    g_debuginf<<" iteration    end_bearing     side_friction   total_force"<<endl;

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    int    stepsnum=0;
    char   stepsstr[4];
    char   stepsfp[50];
    long double avgNormal=0;
    long double avgTangt=0;
    
    int pile[2]={11,12}; // top and down boundaries
    UPDATECTL pilectl[2];

    // iterations start here...
    g_iteration=0;
    do
    {
      // 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();
	
	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	rigidBoundaryForce();

	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation

	// displacement control of the pile
	pilectl[0].tran=vec(0,0,-TIMESTEP*PILE_RATE);
	pilectl[1].tran=vec(0,0,-TIMESTEP*PILE_RATE);

	updateRB(pile, pilectl, 2); 
	updateRectPile();
	if (g_iteration % interval == 0) {
	    long double  f7=getShearForce( 7).getz();
	    long double  f8=getShearForce( 8).getz();
	    long double  f9=getShearForce( 9).getz();
	    long double f10=getShearForce(10).getz();
	    long double  fn=getNormalForce(12).getz();
	    g_debuginf<<setw(10)<<g_iteration
		      <<setw(16)<<fn
		      <<setw(16)<<(f7+f8+f9+f10)
		      <<setw(16)<<(fn+f7+f8+f9+f10)
		      <<endl;
	}

	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    printRectPile(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output statistics info.
	if (g_iteration % interval == 0) {
	    long double t1=getTransEnergy();
	    long double t2=getRotatEnergy();
	    long double t3=getPotenEnergy(-0.025);
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<t1
		       <<setw(16)<<t2
		       <<setw(16)<<(t1+t2)
		       <<setw(16)<<t3
		       <<setw(16)<<(t1+t2+t3)<<endl;
	}

	// 8. loop break condition
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);
    printRectPile(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


// The specimen has been deposited with gravitation within boundaries composed of particles.
// An ellipsoidal pile is then drived into the particles using displacement control.
void assembly::ellipPile_Disp(int   total_steps,  
			      int   snapshots, 
			      int   interval,
			      long double dimn,
			      long double rsize,
			      const char* iniptclfile,
			      const char* particlefile, 
			      const char* contactfile,  
			      const char* progressfile,
			      const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1); }
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"pile penetrate..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average       translational    rotational       "
	       <<"kinetic        potential         total           void            sample       coordination"
	       <<"       sample           sample          sample          sample          sample          sample"
	       <<"          sample          sample          sample         sample           sample         "
	       <<" sample          sample          sample          sample          sample"<<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"         omga            force           moment         energy           energy          "
	       <<"energy         energy            energy          ratio          porosity         number       "
	       <<"   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       "
	       <<"epsilon-v"<<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 

    // pre_3. define variables used in iterations
    long double l13, l24, l56;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    long double void_ratio=0;
    
    // iterations start here...
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 5. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - (-dimn/2);
	l24=dimn*rsize;
	l13=dimn*rsize;
	Volume=l13*l24*l56-ellipPilePeneVol();
	void_ratio=Volume/getParticleVolume()-1;

	// 6. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 6. (2) output statistics info.
	if (g_iteration % interval == 0) {
	    long double t1=getTransEnergy();
	    long double t2=getRotatEnergy();
	    long double t3=getPotenEnergy(-0.025);
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<t1
		       <<setw(16)<<t2
		       <<setw(16)<<(t1+t2)
		       <<setw(16)<<t3
		       <<setw(16)<<(t1+t2+t3)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*getActualCntctNum()/TotalNum
		       <<endl;
	    g_debuginf<<setw(10)<<g_iteration
		      <<setw(16)<<getTopFreeParticlePosition().getz()
		      <<setw(16)<<ellipPileTipZ()
		      <<setw(16)<<getTopFreeParticlePosition().getz()-ellipPileTipZ()
		      <<setw(16)<<l13*l24*l56
		      <<setw(16)<<ellipPilePeneVol()
		      <<setw(16)<<Volume
		      <<endl;
	}

	// 7. loop break condition
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


// The specimen has been deposited with gravitation within rigid boundaries.
// An ellipsoidal penetrator is then impacted into the particles with initial velocity.
void assembly::ellipPile_Impact(int   total_steps,  
				int   snapshots, 
				int   interval,
				long double dimn,
				const char* iniptclfile,
				const char* inibdryfile,
				const char* particlefile, 
				const char* contactfile,  
				const char* progressfile,
				const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1); }
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"penetrator impact..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average       translational    rotational       "
	       <<"kinetic        potential         total           void            sample       coordination"
	       <<"       sample           sample          sample          sample          sample          sample"
	       <<"          sample          sample          sample         sample           sample         "
	       <<" sample          sample          sample          sample          sample"<<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"         omga            force           moment         energy           energy          "
	       <<"energy         energy            energy          ratio          porosity         number       "
	       <<"   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       "
	       <<"epsilon-v"<<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles
    createBoundary(inibdryfile);   // create boundaries.

    // pre_3. define variables used in iterations
    long double l13, l24, l56;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    long double void_ratio=0;
    long double bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }
    
    // iterations start here...
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);

	// 4. calculate boundary forces/moments and apply them to particles.
	rigidBoundaryForce(bdry_penetr, bdry_cntnum);
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - (-dimn/2);
	l24=dimn;
	l13=dimn;
	Volume=l13*l24*l56-ellipPilePeneVol();
	void_ratio=Volume/getParticleVolume()-1;

	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output statistics info.
	if (g_iteration % interval == 0) {
	    long double t1=getTransEnergy();
	    long double t2=getRotatEnergy();
	    long double t3=getPotenEnergy(-0.025);
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<t1
		       <<setw(16)<<t2
		       <<setw(16)<<(t1+t2)
		       <<setw(16)<<t3
		       <<setw(16)<<(t1+t2+t3)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    g_debuginf<<setw(10)<<g_iteration
		      <<setw(16)<<bdry_penetr[1]
		      <<setw(16)<<bdry_penetr[2]
		      <<setw(16)<<bdry_penetr[3]
		      <<setw(16)<<bdry_penetr[4]
		      <<setw(16)<<bdry_penetr[6]
		      <<setw(16)<<bdry_cntnum[1]
		      <<setw(16)<<bdry_cntnum[2]
		      <<setw(16)<<bdry_cntnum[3]
		      <<setw(16)<<bdry_cntnum[4]
		      <<setw(16)<<bdry_cntnum[6]
		      <<endl;
	    /*
	    g_debuginf<<setw(10)<<g_iteration
		      <<setw(16)<<getTopFreeParticlePosition().getz()
		      <<setw(16)<<ellipPileTipZ()
		      <<setw(16)<<getTopFreeParticlePosition().getz()-ellipPileTipZ()
		      <<setw(16)<<l13*l24*l56
		      <<setw(16)<<ellipPilePeneVol()
		      <<setw(16)<<Volume
		      <<endl;
	    */
	}

	// 8. loop break condition
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}


// The specimen has been deposited with gravitation within particle boundaries.
// An ellipsoidal penetrator is then impacted into the particles with initial velocity.
void assembly::ellipPile_Impact_p(int   total_steps,  
				  int   snapshots, 
				  int   interval,
				  long double dimn,
				  const char* iniptclfile,
				  const char* particlefile, 
				  const char* contactfile,  
				  const char* progressfile,
				  const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1); }
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"penetrator impact..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average       translational    rotational       "
	       <<"kinetic        potential         total           void            sample       coordination"
	       <<"       sample           sample          sample          sample          sample          sample"
	       <<"          sample          sample          sample         sample           sample         "
	       <<" sample          sample          sample          sample          sample"<<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"         omga            force           moment         energy           energy          "
	       <<"energy         energy            energy          ratio          porosity         number       "
	       <<"   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       "
	       <<"epsilon-v"<<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles

    // pre_3. define variables used in iterations
    long double l13, l24, l56;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    long double void_ratio=0;
    
    // iterations start here...
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 5. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - (-dimn/2);
	l24=dimn;
	l13=dimn;
	Volume=l13*l24*l56-ellipPilePeneVol();
	void_ratio=Volume/getParticleVolume()-1;

	// 6. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 6. (2) output statistics info.
	if (g_iteration % interval == 0) {
	    long double t1=getTransEnergy();
	    long double t2=getRotatEnergy();
	    long double t3=getPotenEnergy(-0.025);
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<t1
		       <<setw(16)<<t2
		       <<setw(16)<<(t1+t2)
		       <<setw(16)<<t3
		       <<setw(16)<<(t1+t2+t3)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*getActualCntctNum()/TotalNum
		       <<endl;
	    g_debuginf<<setw(10)<<g_iteration
		      <<setw(16)<<getTopFreeParticlePosition().getz()
		      <<setw(16)<<ellipPileTipZ()
		      <<setw(16)<<getTopFreeParticlePosition().getz()-ellipPileTipZ()
		      <<setw(16)<<l13*l24*l56
		      <<setw(16)<<ellipPilePeneVol()
		      <<setw(16)<<Volume
		      <<endl;
	}

	// 7. loop break condition
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    g_debuginf.close();
}



// The specimen has been deposited with gravitation within boundaries composed of particles.
// An ellipsoidal pile is then drived into the particles using force control.
// Not recommended.
void assembly::ellipPile_Force(int   total_steps,  
			       int   snapshots,
			       int   interval,
			       long double dimn,
			       long double force,
			       int   division,
			       const char* iniptclfile,
			       const char* particlefile, 
			       const char* contactfile,  
			       const char* progressfile,
			       const char* balancedfile,
			       const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile); 
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1); }
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"pile penetrate..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average       translational    rotational       "
	       <<"kinetic        potential         total           void            sample       coordination"
	       <<"       sample           sample          sample          sample          sample          sample"
	       <<"          sample          sample          sample         sample           sample         "
	       <<" sample          sample          sample          sample          sample"<<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"         omga            force           moment         energy           energy          "
	       <<"energy         energy            energy          ratio          porosity         number       "
	       <<"   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       "
	       <<"epsilon-v"<<endl;

    ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(ios::scientific, ios::floatfield);
    balancedinf<<"pile penetrate..."<<endl
	       <<"   iteration   apply_force    pile_tip_pos     pile_force"<<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 

    // pre_3. define variables used in iterations
    long double l13, l24, l56;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    long double void_ratio=0;

    long double zforce_inc=force/division;
    long double zforce=zforce_inc;

    // iterations start here...
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();

	// 5. calculate specimen void ratio.
	l56=getTopFreeParticlePosition().getz() - (-dimn/2);
	l24=dimn;
	l13=dimn;
	Volume=l13*l24*l56-ellipPilePeneVol();
	void_ratio=Volume/getParticleVolume()-1;
	
	// 6. update pile external force and position
	if(zforce>ellipPileForce())
	    ellipPileUpdate();

	if(fabsl(ellipPileForce()-zforce)/zforce < STRESS_ERROR ){
	    balancedinf<<setw(10)<<g_iteration
		       <<setw(16)<<zforce
		       <<setw(16)<<getTopFreeParticlePosition().getz()-ellipPileTipZ()
		       <<setw(16)<<ellipPileForce()
		       <<endl;
	    zforce += zforce_inc;
	}

	if( g_iteration % interval == 0){
	    g_debuginf<<setw(10)<<g_iteration
		      <<setw(16)<<zforce
		      <<setw(16)<<getTopFreeParticlePosition().getz()-ellipPileTipZ()
		      <<setw(16)<<ellipPileForce()
		      <<endl;
	}

	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);

	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output statistics info.
	if (g_iteration % interval == 0) {
	    long double t1=getTransEnergy();
	    long double t2=getRotatEnergy();
	    long double t3=getPotenEnergy(-0.025);
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<t1
		       <<setw(16)<<t2
		       <<setw(16)<<(t1+t2)
		       <<setw(16)<<t3
		       <<setw(16)<<(t1+t2+t3)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*getActualCntctNum()/TotalNum
		       <<endl;
	}

	// 8. loop break condition
	if (fabsl((zforce-force)/force)<0.001)
	    break;
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}


// The specimen has been isotropically compressed to confining pressure sigma_a. This function
// performs true triaxial test. Force boundaries are used.
void assembly::truetriaxial(int   total_steps,  
			    int   snapshots, 
			    int   interval,
			    long double sigma_a,     
			    long double sigma_w,
			    long double sigma_l,     
			    long double sigma_h,   
			    int   sigma_division,
			    const char* iniptclfile,  
			    const char* inibdryfile,
			    const char* particlefile, 
			    const char* boundaryfile,
			    const char* contactfile,  
			    const char* progressfile,
			    const char* balancedfile, 
			    const char* debugfile) 
{
    // pre_1: open streams for output
    // particlefile and contactfile are used for snapshots at the end.
    progressinf.open(progressfile);
    if(!progressinf) { cout<<"stream error!"<<endl; exit(-1);}
    progressinf.setf(ios::scientific, ios::floatfield);
    progressinf<<"true triaxial..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(ios::scientific, ios::floatfield);
    balancedinf<<"true triaxial..."<<endl
	       <<"     iteration possible  actual      average	    average         average         average"
	       <<"         average         average         average        sample            sample     "
	       <<"     sample          sample          sample          sample          sample          "
	       <<"sample          sample         sample           sample          sample          sample     "
	       <<"     sample          sample          sample          void            sample        coordinate"
	       <<endl
	       <<"       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
	       <<"          omga            force           moment        density          "
	       <<"sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
	       <<"sigma3_1        sigma3_2           p             width          length           "
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon-v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(ios::scientific, ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    long double W0 = getApt(2).gety()-getApt(4).gety();
    long double L0 = getApt(1).getx()-getApt(3).getx();
    long double H0 = getApt(5).getz()-getApt(6).getz();
    long double l13, l24, l56, min_area, mid_area, max_area;
    long double sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    long double epsilon_w, epsilon_l, epsilon_h;
    long double avgNormal=0;
    long double avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int mid[2]={1,3};    // boundary 1 and 3
    int max[2]={2,4};    // boundary 2 and 4
    int min[2]={5,6};    // boundary 5 and 6
    UPDATECTL midctl[2];
    UPDATECTL maxctl[2];
    UPDATECTL minctl[2];
    long double void_ratio=0;
    long double bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    long double sigma_w1=sigma_a;
    long double sigma_l1=sigma_a;
    long double sigma_h1=sigma_a;
    long double sigma_w_inc=(sigma_w-sigma_a)/sigma_division;
    long double sigma_l_inc=(sigma_l-sigma_a)/sigma_division;
    long double sigma_h_inc=(sigma_h-sigma_a)/sigma_division;

    // iterations start here...
    g_iteration=0;
    do
    {
	// 1. create possible boundary particles and contacts between particles
	findContact();
	findParticleOnBoundary();

	// 2. set particles' forces/moments as zero before each re-calculation
	clearForce();	

	// 3. calculate contact forces/moments and apply them to particles
	internalForce(avgNormal, avgTangt);
	
	// 4. calculate boundary forces/moments and apply them to particles
	rigidBoundaryForce();
	
	// 5. update particles' velocity/omga/position/orientation based on force/moment
	updateParticle();
	
	// 6. update boundaries' position and orientation
	l56=getApt(5).getz()-getApt(6).getz();
	l24=getApt(2).gety()-getApt(4).gety();
	l13=getApt(1).getx()-getApt(3).getx();    Volume=l13*l24*l56;
	min_area=l13*l24;    mid_area=l56*l24;    max_area=l56*l13;
	setArea(5,min_area); setArea(6,min_area); setArea(1,mid_area);
	setArea(3,mid_area); setArea(2,max_area); setArea(4,max_area);
	sigma1_1=vfabsl(getNormalForce(2))/max_area; sigma1_2=vfabsl(getNormalForce(4))/max_area;
	sigma2_1=vfabsl(getNormalForce(1))/mid_area; sigma2_2=vfabsl(getNormalForce(3))/mid_area;
	sigma3_1=vfabsl(getNormalForce(5))/min_area; sigma3_2=vfabsl(getNormalForce(6))/min_area;
	void_ratio=Volume/getParticleVolume()-1;

	if (sigma3_1<sigma_h1)
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	else
	    minctl[0].tran=vec(0,0,TIMESTEP*RELEASE_RATE);
	
	if (sigma3_2<sigma_h1)
	    minctl[1].tran=vec(0,0,TIMESTEP*COMPRESS_RATE);
	else
	    minctl[1].tran=vec(0,0,-TIMESTEP*RELEASE_RATE);
	
	if (sigma2_1<sigma_l1)
	    midctl[0].tran=vec(-TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[0].tran=vec(TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma2_2<sigma_l1)
	    midctl[1].tran=vec(TIMESTEP*COMPRESS_RATE,0,0);
	else
	    midctl[1].tran=vec(-TIMESTEP*RELEASE_RATE,0,0);
	
	if (sigma1_1<sigma_w1)
	    maxctl[0].tran=vec(0,-TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[0].tran=vec(0,TIMESTEP*RELEASE_RATE,0);
	
	if (sigma1_2<sigma_w1)
	    maxctl[1].tran=vec(0,TIMESTEP*COMPRESS_RATE,0);
	else
	    maxctl[1].tran=vec(0,-TIMESTEP*RELEASE_RATE,0);
	
	updateRB(min,minctl,2);
	updateRB(mid,midctl,2);
	updateRB(max,maxctl,2);
	updateRB6();
	
	// 7. (1) output particles and contacts information
	if (g_iteration % (total_steps/snapshots) == 0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp,particlefile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printParticle(stepsfp);
	    
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    printContact(stepsfp);
	    ++stepsnum;
	}

	// 7. (2) output stress and strain info
	epsilon_w = (W0-l24)/W0; epsilon_l = (L0-l13)/L0; epsilon_h = (H0-l56)/H0;
	if (g_iteration % interval == 0 ){
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()   
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    g_debuginf<<setw(10)<<g_iteration
		      <<setw(16)<<bdry_penetr[1]
		      <<setw(16)<<bdry_penetr[2]
		      <<setw(16)<<bdry_penetr[3]
		      <<setw(16)<<bdry_penetr[4]
		      <<setw(16)<<bdry_penetr[5]
		      <<setw(16)<<bdry_penetr[6]
		      <<setw(16)<<bdry_cntnum[1]
		      <<setw(16)<<bdry_cntnum[2]
		      <<setw(16)<<bdry_cntnum[3]
		      <<setw(16)<<bdry_cntnum[4]
		      <<setw(16)<<bdry_cntnum[5]
		      <<setw(16)<<bdry_cntnum[6]
		      <<endl;
	}

	// 8. find the balanced status and increase confining pressure
	if (   fabsl(sigma1_1-sigma_w1)/sigma_w1 < STRESS_ERROR && fabsl(sigma1_2-sigma_w1)/sigma_w1 < STRESS_ERROR
	    && fabsl(sigma2_1-sigma_l1)/sigma_l1 < STRESS_ERROR && fabsl(sigma2_2-sigma_l1)/sigma_l1 < STRESS_ERROR
	    && fabsl(sigma3_1-sigma_h1)/sigma_h1 < STRESS_ERROR && fabsl(sigma3_2-sigma_h1)/sigma_h1 < STRESS_ERROR ) {
	    balancedinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()    
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    sigma_w1 += sigma_w_inc;
	    sigma_l1 += sigma_l_inc;
	    sigma_h1 += sigma_h_inc;
	}

	// 9. loop break condition
	if (   fabsl(sigma1_1-sigma_w)/sigma_w < STRESS_ERROR && fabsl(sigma1_2-sigma_w)/sigma_w < STRESS_ERROR
	    && fabsl(sigma2_1-sigma_l)/sigma_l < STRESS_ERROR && fabsl(sigma2_2-sigma_l)/sigma_l < STRESS_ERROR
	    && fabsl(sigma3_1-sigma_h)/sigma_h < STRESS_ERROR && fabsl(sigma3_2-sigma_h)/sigma_h < STRESS_ERROR ) {
	    progressinf<<setw(10)<<g_iteration
		       <<setw(10)<<getPossCntctNum()
		       <<setw(10)<<getActualCntctNum()
		       <<setw(16)<<getAveragePenetration()
		       <<setw(16)<<avgNormal
		       <<setw(16)<<avgTangt
		       <<setw(16)<<getAverageVelocity() 
		       <<setw(16)<<getAverageOmga()
		       <<setw(16)<<getAverageForce()    
		       <<setw(16)<<getAverageMoment()
		       <<setw(16)<<getDensity()
		       <<setw(16)<<sigma1_1<<setw(16)<<sigma1_2
		       <<setw(16)<<sigma2_1<<setw(16)<<sigma2_2
		       <<setw(16)<<sigma3_1<<setw(16)<<sigma3_2
		       <<setw(16)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(16)<<l24<<setw(16)<<l13<<setw(16)<<l56
		       <<setw(16)<<Volume
		       <<setw(16)<<epsilon_w
		       <<setw(16)<<epsilon_l
		       <<setw(16)<<epsilon_h
		       <<setw(16)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(16)<<void_ratio
		       <<setw(16)<<void_ratio/(1+void_ratio)
		       <<setw(16)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    break;
	}
	
    } while (++g_iteration < total_steps);

    // post_1. store the final snapshot of particles, contacts and boundaries.
    strcpy(stepsfp, particlefile); strcat(stepsfp, "_end");
    printParticle(stepsfp);

    strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
    printContact(stepsfp);

    strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
    printBoundary(stepsfp);
    
    // post_2. close streams
    progressinf.close();
    balancedinf.close();
    g_debuginf.close();
}

} // namespace dem ends

/* 
void assembly::dircShear(long double rate, long double roterate,long double stress,const char* iniptclfile,
						 const char* boundaryfile, const char* responsefile, const char* resultfile,
						 const char* trackfile){
	createSample(iniptclfile);//create particles 
	createBoundary(boundaryfile);//create rigid boundaries

	FILE* fprslt=fopen(responsefile,"w");

	FILE* fp;
	fp=fopen(trackfile,"w");

	clearForce();

	int upanddown[2]={5,6};
	UPDATECTL updownctl[2];
	updownctl[0].expnd=1;
	updownctl[0].fixpt=0;
	updownctl[0].rote=0;
	updownctl[0].tran=vec(0,0,rate*TIMESTEP);
	updownctl[1]=updownctl[0];

	int load[2]={2,4};
	UPDATECTL loadctl[2];
	loadctl[0].expnd=1;
	loadctl[0].fixpt=getApt(2);
	loadctl[0].rote=vec(roterate*TIMESTEP,0,0);
	loadctl[0].tran=0;
	loadctl[1]=loadctl[0];
	loadctl[1].fixpt=getApt(4);

	list<RGDBDRY*>::iterator rt;
	fprintf(fprslt,"bdry_1_norm_x  bdry_1_norm_y  bdry_1_norm_z  bdry_1_shar_x  bdry_1_shar_y  bdry_1_shar_z  \
bdry_2_norm_x  bdry_2_norm_y  bdry_2_norm_z  bdry_2_shar_x  bdry_2_shar_y  bdry_2_shar_z  \
bdry_3_norm_x  bdry_3_norm_y  bdry_3_norm_z  bdry_3_shar_x  bdry_3_shar_y  bdry_3_shar_z  \
bdry_4_norm_x  bdry_4_norm_y  bdry_4_norm_z  bdry_4_shar_x  bdry_4_shar_y  bdry_4_shar_z  \
bdry_5_norm_x  bdry_5_norm_y  bdry_5_norm_z  bdry_5_shar_x  bdry_5_shar_y  bdry_5_shar_z  \
bdry_6_norm_x  bdry_6_norm_y  bdry_6_norm_z  bdry_6_shar_x  bdry_6_shar_y  bdry_6_shar_z\n");

	long double avgsigma;
	long double l13, l24, l56, min_area, mid_area, max_area, lead;
	long double sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
	vec tmpnorm, tmpshar;
	long double av=0;
	long double ao=0;
	long double af=0;
	long double am=0;
	long double avgNormal=0;
	long double avgTangt=0;

	progressinf<<"DircShearing..."<<endl
	         <<"iter_num   "
                 <<"init_contact  "
                 <<"contact  "
	         <<"normal force      "
                 <<"velocity        "
	         <<"omga            "
	         <<"force           "
	         <<"moment"<<endl;

	g_iteration=0;
	do{
                cout<<"DircShearing..."<<g_iteration<<endl;
		progressinf<<setw(10)<<g_iteration;

		findParticleOnBoundary();
		findContact();

		internalForce(avgNormal, avgTangt);
		rigidBoundaryForce();
		//track(fp,5);
		progressinf<<setw(16)<<getAverageVelocity()
		         <<setw(16)<<getAverageOmga()
		         <<setw(16)<<getAverageForce()
			 <<setw(16)<<getAverageMoment()<<endl;
		contactUpdate();
		updateParticle();

		l56=getApt(5).getz()-getApt(6).getz();
		l24=getApt(2).gety()-getApt(4).gety();
		l13=getApt(1).getx()-getApt(3).getx();
		min_area=l13*l24;
		mid_area=l56*l24;
		lead=fabsl(normalize(getDirc(2))%vec(0,1,0));
		max_area=l56*l13;
		setArea(5,min_area);
		setArea(6,min_area);
		setArea(1,mid_area);
		setArea(3,mid_area);
		setArea(2,max_area);
		setArea(4,max_area);
		avgsigma=getAverageRigidPressure();
		printf("avgsigma=%15.3lf\n",avgsigma);
		sigma1_1=fabsl(getNormalForce(2))/max_area;
		sigma1_2=fabsl(getNormalForce(4))/max_area;
		sigma2_1=fabsl(getNormalForce(1))/mid_area;
		sigma2_2=fabsl(getNormalForce(3))/mid_area;
		sigma3_1=fabsl(getNormalForce(5))/min_area;
		sigma3_2=fabsl(getNormalForce(6))/min_area;
		if(sigma3_1<stress)
			updownctl[0].tran=vec(0,0,-rate*TIMESTEP);
		else
			updownctl[0].tran=vec(0,0,rate*TIMESTEP);
		if(sigma3_2<stress)
			updownctl[1].tran=vec(0,0,rate*TIMESTEP);
		else
			updownctl[1].tran=vec(0,0,-rate*TIMESTEP);
		updateRB(upanddown,updownctl,2);
		if (1){
			updateRB(load,loadctl,2);
			fprintf(fprslt,"%15.6lf",lead);
			for(rt=RBList.begin();rt!=RBList.end();++rt){
				tmpnorm=(*rt)->getNormalForce()/(*rt)->getArea()/1000;
				tmpshar=(*rt)->getShearForce()/(*rt)->getArea()/1000;
				fprintf(fprslt,"%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf",
					tmpnorm.x,tmpnorm.y,tmpnorm.z,
					tmpshar.x,tmpshar.y,tmpshar.z);
			}
			fprintf(fprslt,"\n");
		}
	}while(++g_iteration<10000);
//	displayBoundary();
	fclose(fp);
	fclose(fprslt);
	printParticle(resultfile);
}
*/

/* 
void assembly::soft_tric(long double _sigma3,long double _b,const char* iniptclfile,
						   const char* boundaryfile,const char* responsefile,
						   const char* resultfile,const char* trackfile){
	createSample(iniptclfile); //create particles 
	createBoundary(boundaryfile);

	FILE* fprslt=fopen(responsefile,"w");
	FILE* fp=fopen(trackfile,"w");
	
	clearForce();

	int pre_it=0;
	int pre_snap=0;
	int snapnum=0;
	char snapfile[80];

	list<RGDBDRY*>::iterator rt;

	int max[2]={1,2};//maximum stress acting on boundary 5 and 6
	UPDATECTL maxctl[2];
	long double loading_rate=0.01;

	long double avgsigma;
	long double af, av, am, ao, adr, pre_af;
	vec disp, tmp;
	av=ao=af=am=adr=pre_af=0;

	progressinf<<"Soft_tric..."<<endl
	         <<"iter_num   "
                 <<"init_contact  "
                 <<"contact  "
	         <<"normal force      "
                 <<"velocity        "
	         <<"omga            "
	         <<"force           "
	         <<"moment          "
		 <<"friction"<<endl;

	g_iteration=0;
	do{
                cout<<"Soft_tric..."<<g_iteration<<endl;
		progressinf<<setw(10)<<g_iteration;

		findParticleOnBoundary();
		findParticleOnLine();
		createFlbNet();
		flexiBoundaryForceZero();
		flexiBoundaryForce();
		findContact();

		initFBForce();
		internalForce();
		rigidBoundaryForce();
		//track(fp,5);
		progressinf<<setw(16)<<(av=getAverageVelocity())
		         <<setw(16)<<(ao=getAverageOmga())
		         <<setw(16)<<(af=getAverageForce())
		         <<setw(16)<<(am=getAverageMoment())
			 <<setw(16)<<(adr=avgDgrFric());
		contactUpdate();
		
		avgsigma=getAverageRigidPressure();
		maxctl[0].tran=TIMESTEP*vec(0,0,-loading_rate);
		maxctl[1].tran=TIMESTEP*vec(0,0,loading_rate);
		if(af<0.03&&g_iteration-pre_it>=20||g_iteration-pre_it>=500){
		//if(g_iteration-pre_it>=50){
			pre_it=g_iteration;
		        updateRB(max,maxctl,2);
			if(g_iteration-pre_snap>=5000){
			    snapnum++;
			    sprintf(snapfile,"%s%d","snap",snapnum);
			    pre_snap=g_iteration;
			    snapshot(snapfile);
			}
			fprintf(fprslt,"%10d%6d%10.3lf%15.6lf%15.6lf%15.6lf%15.6lf",g_iteration,ActualCntctNum,adr,af,am,av,ao);
			for(rt=RBList.begin();rt!=RBList.end();++rt){
				disp=(*rt)->getApt();
				tmp=(*rt)->getNormalForce()/1000/(*rt)->getArea();
				fprintf(fprslt,"%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf",
					disp.getx(),disp.gety(),disp.getz(),tmp.getx(),tmp.gety(),tmp.getz());
			}
			fprintf(fprslt,"\n");
		}
	updateParticle();
	pre_af=af;
	}while(++g_iteration<1000000);
//	displayBoundary();
	fclose(fp);
	fclose(fprslt);
	printParticle(resultfile);
}//end of soft_tric
*/

/* 
void assembly::shallowFoundation(const char* iniptclfile, const char* boundaryfile,const char* responsefile, 
	const char* resultfile, const char* trackfile)
{
	createSample(iniptclfile);//create particles 
	createBoundary(boundaryfile);

	FILE* fprslt=fopen(responsefile,"w");

	FILE* fp;
	fp=fopen(trackfile,"w");

	int pre_it=0;
	int pre_snap=0;
	int snapnum=0;
	char snapfile[80];

	clearForce();

	list<RGDBDRY*>::iterator rt;

//	int mid[2]={2,4};//intermediate stress acting on boundary 2 and 4
//	UPDATECTL midctl[2];
	int max[2]={5,6};//maximum stress acting on boundary 5 and 6
	UPDATECTL maxctl[2];
//	int min[2]={1,3};//minimum stress acting on boundary 1 and 3
//	UPDATECTL minctl[2];
	long double loading_rate=0.01;

	long double avgsigma;
	long double af, av, am, ao, adr, pre_af;
	int nbdry;
	vec disp, tmp, zbdry_velocity_0;
	av=ao=af=am=adr=pre_af=0;

	progressinf<<"Shallow Foundation..."<<endl
	         <<"iter_num   "
                 <<"init_contact  "
                 <<"contact  "
	         <<"normal force      "
                 <<"velocity        "
	         <<"omga            "
	         <<"force           "
	         <<"moment"<<endl;

	g_iteration=0;
	do{
                cout<<"Shallow Foundation..."<<g_iteration<<endl;
		progressinf<<setw(10)<<g_iteration;

		findParticleOnBoundary();
		findParticleOnLine();
		createFlbNet();
		flexiBoundaryForceZero();
		flexiBoundaryForce();
		findContact();
		
		initFBForce();

		internalForce();
		rigidBoundaryForce();
		//track(fp,5);
		progressinf<<setw(16)<<(av=getAverageVelocity())
		         <<setw(16)<<(ao=getAverageOmga())
		         <<setw(16)<<(af=getAverageForce())
		         <<setw(16)<<(am=getAverageMoment())
			 <<setw(16)<<(adr=avgDgrFric());

		contactUpdate();
		
		avgsigma=getAverageRigidPressure();
		zbdry_velocity_0=vec(0,0,-loading_rate);
		maxctl[0].tran=TIMESTEP*zbdry_velocity_0;
		if(1){
			if(af<0.05&&g_iteration-pre_it>=20
				  ||g_iteration-pre_it>=500){
			pre_it=g_iteration;
		        updateRB(&max[0],&maxctl[0],1);
			nbdry=1;

			if(g_iteration-pre_snap>=20000){
			    snapnum++;
			    sprintf(snapfile,"%s%d","snap",snapnum);
			    pre_snap=g_iteration;
			    snapshot(snapfile);
			}
			fprintf(fprslt,"%10d%6d%10.3lf%15.6lf%15.6lf%15.6lf%15.6lf",g_iteration,ActualCntctNum,adr,af,am,av,ao);
			for(rt=RBList.begin();rt!=RBList.end();++rt,++nbdry){
				disp=(*rt)->getApt();
				tmp=(*rt)->getNormalForce()/1000;
				tmp/=getArea(nbdry);
				fprintf(fprslt,"%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf",
					disp.getx(),disp.gety(),disp.getz(),tmp.getx(),tmp.gety(),tmp.getz());
			}
			fprintf(fprslt,"\n");
			}
		}
		updateParticle();
		pre_af=af;
	}while(++g_iteration<1000000);
//	displayBoundary();
	fclose(fp);
	fclose(fprslt);
	printParticle(resultfile);
}
*/

/* 
void assembly::simpleShear(long double _sigma3,long double _b,
			const char* iniptclfile,const char* boundaryfile,
			const char* responsefile,const char* resultfile, const char* trackfile)
{
	createSample(iniptclfile);//create particles 
	createBoundary(boundaryfile);
	FILE* fprslt=fopen(responsefile,"w");

	FILE* fp;
	fp=fopen(trackfile,"w");

	clearForce();

	int pre_it=0;
	int pre_snap=0;
	int snapnum=0;
	char snapfile[80];

	list<RGDBDRY*>::iterator rt;

	int mid[2]={2,4};//intermediate stress acting on boundary 2 and 4
	UPDATECTL midctl[2];
	int max[2]={5,6};//maximum stress acting on boundary 5 and 6
	UPDATECTL maxctl[2];
	int min[2]={1,3};//minimum stress acting on boundary 1 and 3
	UPDATECTL minctl[2];
//	long double loading_rate=0.01;
	long double increment=0.0001;
	long double angular_velocity=0.1;
	vec increment_velocity_x(increment,0,0);
	vec increment_velocity_y(0,increment,0);
	vec increment_velocity_z(0,0,increment);
	vec xbdry_velocity_0,xbdry_velocity_1;
	vec ybdry_velocity_0,ybdry_velocity_1;
	vec zbdry_velocity_0,zbdry_velocity_1;

	long double avgsigma;
	long double af, av, am, ao, adr, pre_af;
	long double sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2, sigma2;
	long double ita1_1, ita1_2, ita2_1, ita2_2, ita3_1, ita3_2;
	long double ar;
	int nbdry;
	vec disp, angl, nm, sh;
	av=ao=af=am=adr=pre_af=0;

	progressinf<<"SimpleShearing..."<<endl
	         <<"iter_num   "
                 <<"init_contact  "
                 <<"contact  "
	         <<"normal force      "
                 <<"velocity        "
	         <<"omga            "
	         <<"force           "
	         <<"moment          "
		 <<"friction"<<endl;

	g_iteration=0;
	do{
                cout<<"SimpleShearinging..."<<g_iteration<<endl;
		progressinf<<setw(10)<<g_iteration;

		findParticleOnBoundary();
		findContact();

		internalForce();
		rigidBoundaryForce();
		flexiBoundaryForce();
		//track(fp,5);
		progressinf<<setw(16)<<(av=getAverageVelocity())
		         <<setw(16)<<(ao=getAverageOmga())
		         <<setw(16)<<(af=getAverageForce())
		         <<setw(16)<<(am=getAverageMoment())
			 <<setw(16)<<(adr=avgDgrFric());
		contactUpdate();
		
		avgsigma=getAverageRigidPressure();
		sigma1_1=fabsl(getNormalForce(5))/getArea(5);
		ita1_1=fabsl(getShearForce(5))/getArea(5);
		sigma1_2=fabsl(getNormalForce(6))/getArea(6);
		ita1_2=fabsl(getShearForce(6))/getArea(6);
		sigma2_1=fabsl(getNormalForce(2))/getArea(2);
		ita2_1=fabsl(getShearForce(2))/getArea(2);
		sigma2_2=fabsl(getNormalForce(4))/getArea(4);
		ita2_2=fabsl(getShearForce(4))/getArea(4);
		sigma3_1=fabsl(getNormalForce(1))/getArea(1);
		ita3_1=fabsl(getShearForce(1))/getArea(1);
		sigma3_2=fabsl(getNormalForce(3))/getArea(3);
		ita3_2=fabsl(getShearForce(3))/getArea(1);
		if(sigma3_1<_sigma3)
			xbdry_velocity_0=-increment_velocity_x;
		else
			xbdry_velocity_0=increment_velocity_x;
		if(sigma3_2<_sigma3)
			xbdry_velocity_1=increment_velocity_x;
		else
			xbdry_velocity_1=-increment_velocity_x;
		sigma2=_sigma3;//_b*sigma1+(1-_b)*_sigma3;
		if(sigma2_1<sigma2)
			ybdry_velocity_0=-increment_velocity_y;
		else
			ybdry_velocity_0=increment_velocity_y;
		if(sigma2_2<sigma2)
			ybdry_velocity_1=increment_velocity_y;
	 	else
			ybdry_velocity_1=-increment_velocity_y;
		if(sigma1_1<_sigma3)
			zbdry_velocity_0=-increment_velocity_z;
		else
			zbdry_velocity_0=increment_velocity_z;
		if(sigma1_2<_sigma3)
			zbdry_velocity_1=increment_velocity_z;
		else
			zbdry_velocity_1=-increment_velocity_z;
		minctl[0].tran=TIMESTEP*xbdry_velocity_0;
		minctl[0].fixpt=getApt(1);
		minctl[0].rote=TIMESTEP*vec(0,0,angular_velocity);
		minctl[1].fixpt=getApt(3);
		minctl[1].tran=TIMESTEP*xbdry_velocity_1;
		minctl[1].rote=TIMESTEP*vec(0,0,angular_velocity);
		midctl[0].fixpt=getApt(2);
		midctl[0].tran=TIMESTEP*ybdry_velocity_0;
		midctl[0].rote=TIMESTEP*vec(0,0,-angular_velocity);
		midctl[1].fixpt=getApt(4);
		midctl[1].tran=TIMESTEP*ybdry_velocity_1;
		midctl[1].rote=TIMESTEP*vec(0,0,-angular_velocity);
		maxctl[0].tran=TIMESTEP*zbdry_velocity_0;
		maxctl[1].tran=TIMESTEP*zbdry_velocity_1;
		//if(af<0.01){
		if(1){
			UPDATECTL tmpctl;

			tmpctl.tran=minctl[0].tran;
			updateRB(&min[0],&tmpctl,1);

			tmpctl.tran=minctl[1].tran;
			updateRB(&min[1],&tmpctl,1);

			tmpctl.tran=midctl[0].tran;
			updateRB(&mid[0],&tmpctl,1);

			tmpctl.tran=midctl[1].tran;
			updateRB(&mid[1],&tmpctl,1);

			tmpctl.tran=maxctl[0].tran;
			updateRB(&max[0],&tmpctl,1);

			tmpctl.tran=maxctl[1].tran;
			updateRB(&max[1],&tmpctl,1);
			
			if(af<0.02&&fabsl(sigma3_1-_sigma3)<0.02*_sigma3
				  &&fabsl(sigma3_2-_sigma3)<0.02*_sigma3
				  &&fabsl(sigma2_1-sigma2)<0.02*_sigma3
				  &&fabsl(sigma2_2-sigma2)<0.02*_sigma3
				  &&fabsl(sigma1_1-_sigma3)<0.02*_sigma3
				  &&fabsl(sigma1_2-_sigma3)<0.02*_sigma3
				  &&g_iteration-pre_it>=20
				  ||g_iteration-pre_it>=500){
			pre_it=g_iteration;
		        updateRB(min,minctl,2);
			updateRB(mid,midctl,2);
			nbdry=1;

			if(g_iteration-pre_snap>=20000){
			    snapnum++;
			    sprintf(snapfile,"%s%d","snap",snapnum);
			    pre_snap=g_iteration;
			    snapshot(snapfile);
			}
			fprintf(fprslt,"%10d%6d%10.3lf%15.6lf%15.6lf%15.6lf%15.6lf",g_iteration,ActualCntctNum,adr,af,am,av,ao);
			for(rt=RBList.begin();rt!=RBList.end();++rt,++nbdry){
				ar=(*rt)->getArea();
				disp=(*rt)->getApt();
				angl=(*rt)->getDirc();
				nm=(*rt)->getNormalForce()/1000/ar;
				sh=(*rt)->getShearForce()/1000/ar;
				fprintf(fprslt,"%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf", 
disp.x,disp.y,disp.z,angl.x,angl.y,angl.z,nm.x,nm.y,nm.z,sh.x,sh.y,sh.z);
			}
			fprintf(fprslt,"\n");
			}
		}
	updateParticle();
	pre_af=af;
	}while(++g_iteration<300000);
//	displayBoundary();
	fclose(fp);
	fclose(fprslt);
	printParticle(resultfile);
}
*/

/* 
void assembly::earthPressure(long double pressure,bool IsPassive, 
				const char* iniptclfile, const char* boundaryfile,
				const char* responsefile, const char* resultfile,
				const char* trackfile)
{
	createSample(iniptclfile);//create particles 
	createBoundary(boundaryfile);

	FILE* fprslt=fopen(responsefile,"w");

	FILE* fp;
	fp=fopen(trackfile,"w");

	int pre_it=0;
	int pre_snap=0;
	int snapnum=0;
	char snapfile[80];
	clearForce();

	list<RGDBDRY*>::iterator rt;

	int wall[1]={1};
	UPDATECTL wallctl[1];
//	long double loading_rate=0.001;

	long double avgsigma;
	long double af, av, am, ao, adr, pre_af;
	int nbdry;
	vec disp, tmp;
	av=ao=af=am=adr=pre_af=0;

	progressinf<<"EarthPressure..."<<endl
	         <<"iter_num   "
                 <<"init_contact  "
                 <<"contact  "
	         <<"normal force      "
                 <<"velocity        "
	         <<"omga            "
	         <<"force           "
	         <<"moment          "
		 <<"friction"<<endl;

	g_iteration=0;
	do{
	        cout<<"EarthPressure..."<<g_iteration<<endl;
		progressinf<<setw(10)<<g_iteration;

		findParticleOnBoundary();
		findParticleOnLine();
		createFlbNet();
		flexiBoundaryForceZero();
		flexiBoundaryForce();
		findContact();
		
                initFBForce();
		internalForce();
		rigidBoundaryForce();
		//track(fp,5);

		progressinf<<setw(16)<<(av=getAverageVelocity())
		         <<setw(16)<<(ao=getAverageOmga())
		         <<setw(16)<<(af=getAverageForce())
		         <<setw(16)<<(am=getAverageMoment())
			 <<setw(16)<<(adr=avgDgrFric());
		contactUpdate();
		
		avgsigma=getAverageRigidPressure();
		wallctl[0].tran=0;
		wallctl[0].fixpt=getApt(1);
		if(IsPassive)
			wallctl[0].rote=TIMESTEP*vec(0,-1.0,0);
		else
			wallctl[0].rote=TIMESTEP*vec(0,1.0,0);
		if(1){
			if(af<0.05 &&g_iteration-pre_it>=20
				  ||g_iteration-pre_it>=500){
			pre_it=g_iteration;
		        updateRB(wall,wallctl,1);
			nbdry=1;

			if(g_iteration-pre_snap>=20000){
			    snapnum++;
			    sprintf(snapfile,"%s%d","snap",snapnum);
			    pre_snap=g_iteration;
			    snapshot(snapfile);
			}
			fprintf(fprslt,"%10d%6d%10.3lf%15.6lf%15.6lf%15.6lf%15.6lf",g_iteration,ActualCntctNum,adr,af,am,av,ao);
			for(rt=RBList.begin();rt!=RBList.end();++rt,++nbdry){
				if(nbdry==1)
					disp=(*rt)->getDirc();
				else
					disp=(*rt)->getApt();
				tmp=(*rt)->getNormalForce()/1000/getArea(nbdry);
				fprintf(fprslt,"%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf",
					disp.getx(),disp.gety(),disp.getz(),tmp.getx(),tmp.gety(),tmp.getz());
			}
			fprintf(fprslt,"\n");
			}
		}
	updateParticle();
	pre_af=af;
	}while(++g_iteration<1000000);
	//displayBoundary();
	fclose(fp);
	fclose(fprslt);
	printParticle(resultfile);
}
*/
