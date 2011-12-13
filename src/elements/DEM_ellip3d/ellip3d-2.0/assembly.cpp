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
#include <cstring>
#include <ctime>

#ifdef OPENMP
#include <omp.h>
#endif

//#define TIME_PROFILE

using std::cout;
using std::setw;
using std::endl;

static time_t timeStamp; // for file timestamping
static struct timeval timew1, timew2; // for wall-clock time record
static struct timeval timep1, timep2; // for internal wall-clock time profiling

namespace dem {
  
  std::ofstream progressinf;
  
  void assembly::printParticle(const char* str) const
  {
    std::ofstream ofs(str);
    if(!ofs) {
      cout<<"stream error in printParticle!"<<endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    ofs<<setw(OWID)<<TotalNum<<setw(OWID)<<1<<endl;
    ofs<<setw(OWID)<<container.getCenter().getx()
       <<setw(OWID)<<container.getCenter().gety()
       <<setw(OWID)<<container.getCenter().getz()
       <<setw(OWID)<<container.getDimx()
       <<setw(OWID)<<container.getDimy()
       <<setw(OWID)<<container.getDimz()<<endl;
    
    ofs<<setw(OWID)<<"ID"
       <<setw(OWID)<<"type"
       <<setw(OWID)<<"radius_a"
       <<setw(OWID)<<"radius_b"
       <<setw(OWID)<<"radius_c"
       <<setw(OWID)<<"position_x"
       <<setw(OWID)<<"position_y"
       <<setw(OWID)<<"position_z"
       <<setw(OWID)<<"axle_a_x"
       <<setw(OWID)<<"axle_a_y"
       <<setw(OWID)<<"axle_a_z"
       <<setw(OWID)<<"axle_b_x"
       <<setw(OWID)<<"axle_b_y"
       <<setw(OWID)<<"axle_b_z"
       <<setw(OWID)<<"axle_c_x"
       <<setw(OWID)<<"axle_c_y"
       <<setw(OWID)<<"axle_c_z"
       <<setw(OWID)<<"velocity_x"
       <<setw(OWID)<<"velocity_y"
       <<setw(OWID)<<"velocity_z"
       <<setw(OWID)<<"omga_x"
       <<setw(OWID)<<"omga_y"
       <<setw(OWID)<<"omga_z"
       <<setw(OWID)<<"force_x"
       <<setw(OWID)<<"force_y"
       <<setw(OWID)<<"force_z"
       <<setw(OWID)<<"moment_x"
       <<setw(OWID)<<"moment_y"
       <<setw(OWID)<<"moment_z"
       <<endl;
    
    vec tmp;
    std::list<particle*>::const_iterator  it;
    for (it=ParticleList.begin();it!=ParticleList.end();++it)
    {
	ofs<<setw(OWID)<<(*it)->getID()
	   <<setw(OWID)<<(*it)->getType()
	   <<setw(OWID)<<(*it)->getA()
	   <<setw(OWID)<<(*it)->getB()
	   <<setw(OWID)<<(*it)->getC();
	
	tmp=(*it)->getCurrPosition();
	ofs<<setw(OWID)<<tmp.getx()
	   <<setw(OWID)<<tmp.gety()
	   <<setw(OWID)<<tmp.getz();
	
	tmp=(*it)->getCurrDirecA();
	ofs<<setw(OWID)<<tmp.getx()
	   <<setw(OWID)<<tmp.gety()
	   <<setw(OWID)<<tmp.getz();
	
	tmp=(*it)->getCurrDirecB();
	ofs<<setw(OWID)<<tmp.getx()
	   <<setw(OWID)<<tmp.gety()
	   <<setw(OWID)<<tmp.getz();
	
	tmp=(*it)->getCurrDirecC();
	ofs<<setw(OWID)<<tmp.getx()
	   <<setw(OWID)<<tmp.gety()
	   <<setw(OWID)<<tmp.getz();
	
	tmp=(*it)->getCurrVelocity();
	ofs<<setw(OWID)<<tmp.getx()
	   <<setw(OWID)<<tmp.gety()
	   <<setw(OWID)<<tmp.getz();
	
	tmp=(*it)->getCurrOmga();
	ofs<<setw(OWID)<<tmp.getx()
	   <<setw(OWID)<<tmp.gety()
	   <<setw(OWID)<<tmp.getz();
	
	tmp=(*it)->getForce();
	ofs<<setw(OWID)<<tmp.getx()
	   <<setw(OWID)<<tmp.gety()
	   <<setw(OWID)<<tmp.getz();
	
	tmp=(*it)->getMoment();
	ofs<<setw(OWID)<<tmp.getx()
	   <<setw(OWID)<<tmp.gety()
	   <<setw(OWID)<<tmp.getz()<<endl;
    }

    ofs.close();
}


void assembly::printRectPile(const char* str)
{
    std::ofstream ofs(str, std::ios_base::app);
    if(!ofs) {
	cout<<"stream error in printRectPile!"<<endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);

    ofs<<setw(OWID)<<8<<setw(OWID)<<6<<endl;
    vec pos[8];
    for(std::list<RGDBDRY*>::iterator rt=RBList.begin();rt!=RBList.end();++rt){
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
	ofs<<setw(OWID)<<pos[i].getx()<<setw(OWID)<<pos[i].gety()<<setw(OWID)<<pos[i].getz()<<endl;

    ofs<<setw(OWID)<<1<<setw(OWID)<<2<<setw(OWID)<<6<<setw(OWID)<<5<<endl
       <<setw(OWID)<<2<<setw(OWID)<<3<<setw(OWID)<<7<<setw(OWID)<<6<<endl
       <<setw(OWID)<<3<<setw(OWID)<<4<<setw(OWID)<<8<<setw(OWID)<<7<<endl
       <<setw(OWID)<<4<<setw(OWID)<<1<<setw(OWID)<<5<<setw(OWID)<<8<<endl
       <<setw(OWID)<<1<<setw(OWID)<<4<<setw(OWID)<<3<<setw(OWID)<<2<<endl
       <<setw(OWID)<<5<<setw(OWID)<<6<<setw(OWID)<<7<<setw(OWID)<<8<<endl;

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
    std::ofstream ofs(str);
    if(!ofs) {
	cout<<"stream error in printContact!"<<endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    ofs<<setw(OWID)<<ActualCntctNum<<endl;
    ofs<<setw(OWID)<<"ptcl_1"
       <<setw(OWID)<<"ptcl_2"
       <<setw(OWID)<<"point1_x"
       <<setw(OWID)<<"point1_y"
       <<setw(OWID)<<"point1_z"
       <<setw(OWID)<<"point2_x"
       <<setw(OWID)<<"point2_y"
       <<setw(OWID)<<"point2_z"
       <<setw(OWID)<<"radius_1"
       <<setw(OWID)<<"radius_2"
       <<setw(OWID)<<"penetration"
       <<setw(OWID)<<"tangt_dispmt"
       <<setw(OWID)<<"contact_radius"
       <<setw(OWID)<<"R0"
       <<setw(OWID)<<"E0"
       <<setw(OWID)<<"normal_force"
       <<setw(OWID)<<"tangt_force"
       <<setw(OWID)<<"contact_x"
       <<setw(OWID)<<"contact_y"
       <<setw(OWID)<<"contact_z"
       <<setw(OWID)<<"normal_x"
       <<setw(OWID)<<"normal_y"
       <<setw(OWID)<<"normal_z"
       <<setw(OWID)<<"tangt_x"
       <<setw(OWID)<<"tangt_y"
       <<setw(OWID)<<"tangt_z"
       <<setw(OWID)<<"vibra_t_step"
       <<setw(OWID)<<"impact_t_step"
       <<endl;
    std::list<CONTACT>::const_iterator it;
    for (it=ContactList.begin();it!=ContactList.end();++it)
	ofs<<setw(OWID)<<it->getP1()->getID()
	   <<setw(OWID)<<it->getP2()->getID()
	   <<setw(OWID)<<it->getPoint1().getx()
	   <<setw(OWID)<<it->getPoint1().gety()
	   <<setw(OWID)<<it->getPoint1().getz()
	   <<setw(OWID)<<it->getPoint2().getx()
	   <<setw(OWID)<<it->getPoint2().gety()
	   <<setw(OWID)<<it->getPoint2().getz()
	   <<setw(OWID)<<it->getRadius1()
	   <<setw(OWID)<<it->getRadius2()
	   <<setw(OWID)<<it->getPenetration()
	   <<setw(OWID)<<it->getTgtDisp()
	   <<setw(OWID)<<it->getContactRadius()
	   <<setw(OWID)<<it->getR0()
	   <<setw(OWID)<<it->getE0()
	   <<setw(OWID)<<it->getNormalForce()
	   <<setw(OWID)<<it->getTgtForce()
	   <<setw(OWID)<<( it->getPoint1().getx()+it->getPoint2().getx() )/2
	   <<setw(OWID)<<( it->getPoint1().gety()+it->getPoint2().gety() )/2
	   <<setw(OWID)<<( it->getPoint1().getz()+it->getPoint2().getz() )/2
	   <<setw(OWID)<<it->NormalForceVec().getx()
	   <<setw(OWID)<<it->NormalForceVec().gety()
	   <<setw(OWID)<<it->NormalForceVec().getz()
	   <<setw(OWID)<<it->TgtForceVec().getx()
	   <<setw(OWID)<<it->TgtForceVec().gety()
	   <<setw(OWID)<<it->TgtForceVec().getz()
	   <<setw(OWID)<<it->getVibraTimeStep()
	   <<setw(OWID)<<it->getImpactTimeStep()
	   <<endl;
    ofs.close();
}

	
void assembly::createSample(const char* str){
    std::ifstream ifs(str);
    if(!ifs) {
	cout<<"stream error in createSample!"<<endl; exit(-1);
    }
    int RORC;
    ifs >> TotalNum >> RORC;

    REAL cx,cy,cz,dx,dy,dz;

    ifs >> cx >> cy >> cz >> dx >> dy >> dz;
    container.set(dx,dy,dz,vec(cx,cy,cz));
    Volume = container.getDimx() * container.getDimy() * container.getDimz();
    
    char s[20];
    ifs>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
       >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s;

    ParticleList.clear();
    int ID, type;
    REAL a, b, c, px,py,pz,dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz;
    REAL vx,vy,vz,omx,omy,omz,fx,fy,fz,mx,my,mz;
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
    std::list<particle*>::iterator ot, it, pt;
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
		if (   ( vfabs(v-u) < (*it)->getA() + (*pt)->getA())
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
    g_debuginf<<setw(OWID)<<timediffsec(); 
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
    std::list<particle*>::iterator it, pt;
    vec u,v;
    for (it=ParticleList.begin();it!=ParticleList.end();++it){
	u=(*it)->getCurrPosition();
	for (pt=it,++pt;pt!=ParticleList.end();++pt){
	    v=(*pt)->getCurrPosition();
	    if (   ( vfabs(v-u) < (*it)->getA() + (*pt)->getA())
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
    g_debuginf<<setw(OWID)<<timediffsec(); 
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
    std::list<particle*>::iterator ot, it, pt;
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
	if (vfabs(v-u) < (*it)->getA() + (*pt)->getA() ) {
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
    g_debuginf<<setw(OWID)<<timediffsec(); 
#endif
	
    ActualCntctNum = ContactList.size();
}

#else

*/

REAL assembly::getDensity() const{
    REAL dens=0;
    std::list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it)
	dens+=(*it)->getMass();
    return dens/=Volume;
}


REAL assembly::getAveragePenetration() const{
    int totalcntct = ContactList.size();
    if (totalcntct==0)
	return 0;
    else {
	REAL pene=0;
	for (std::list<CONTACT>::const_iterator it=ContactList.begin();it!=ContactList.end();++it)
	    pene += it->getPenetration(); 
	return pene/totalcntct;
    }
}


REAL assembly::getVibraTimeStep() const {
    int totalcntct = ContactList.size();
    if (totalcntct == 0)
	return 0;
    else {
	std::list<CONTACT>::const_iterator it=ContactList.begin();
        REAL minTimeStep = it->getVibraTimeStep();
	for (++it; it != ContactList.end(); ++it) {
	  REAL val = it->getVibraTimeStep(); 
	  minTimeStep =  val < minTimeStep ? val : minTimeStep;
	}
	return minTimeStep;
    }
}


REAL assembly::getImpactTimeStep() const {
    int totalcntct = ContactList.size();
    if (totalcntct == 0)
	return 0;
    else {
	std::list<CONTACT>::const_iterator it=ContactList.begin();
        REAL minTimeStep = it->getImpactTimeStep();
	for (++it; it != ContactList.end(); ++it) {
	  REAL val = it->getImpactTimeStep(); 
	  minTimeStep =  val < minTimeStep ? val : minTimeStep;
	}
	return minTimeStep;
    }
}
 

REAL assembly::getAverageVelocity() const{
    REAL avgv=0;
    int count=0;
    std::list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it)
	if ((*it)->getType()==0) {
	    avgv+=vfabs((*it)->getCurrVelocity());
	    count++;
	}
    return avgv/=count;
}


REAL assembly::getAverageOmga() const{
    REAL avgv=0;
    int count=0;
    std::list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it)
	if ((*it)->getType()==0){
	    avgv+=vfabs((*it)->getCurrOmga());
	    count++;
	}
    return avgv/=count;
}


REAL assembly::getAverageForce() const{
    REAL avgv=0;
    int count=0;
    std::list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it)
	if ((*it)->getType()==0){
	    avgv+=vfabs((*it)->getForce());
	    count++;
	}
    return avgv/count;
}


REAL assembly::getAverageMoment() const{
    REAL avgv=0;
    int count=0;
    std::list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it)
	if ((*it)->getType()==0){
	    avgv+=vfabs((*it)->getMoment());
	    count++;
	}
    return avgv/=count;
}


REAL assembly::getParticleVolume() const{
    REAL avgv=0;
    std::list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it)
	if ((*it)->getType()==0)
	    avgv+=(*it)->getVolume();
    return avgv;
}


vec assembly::getTopFreeParticlePosition() const{
    std::list<particle*>::const_iterator it,jt,kt;
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


REAL assembly::ellipPileForce() {
    REAL val=0;
    for(std::list<particle*>::iterator it=ParticleList.begin();it!=ParticleList.end();++it)
	if ((*it)->getType()==3) {
	    val = (*it)->getForce().getz();
	    break;
	}
    return val;
}


vec assembly::ellipPileDimn() {
    vec val;
    for(std::list<particle*>::iterator it=ParticleList.begin();it!=ParticleList.end();++it)
	if ((*it)->getType()==3) {
	    val = vec((*it)->getA(), (*it)->getB(), (*it)->getC());
	    break;
	}
    return val;
}


REAL assembly::ellipPileTipZ() {
    REAL val=0;
    for(std::list<particle*>::iterator it=ParticleList.begin();it!=ParticleList.end();++it)
	if ((*it)->getType()==3) {
	    val = (*it)->getCurrPosition().getz()-(*it)->getA();
	    break;
	}
    return val;
}


REAL assembly::ellipPilePeneVol() {
    REAL val=0;
    if (getTopFreeParticlePosition().getz()-ellipPileTipZ()<=0)
	val=0;
    else{
	// low: a signed number as lower limit for volumetric integration
	REAL low=ellipPileTipZ() + ellipPileDimn().getx() - getTopFreeParticlePosition().getz(); 
	REAL lowint=low-pow(low,3)/3.0/pow(ellipPileDimn().getx(),2);
	val = PI * ellipPileDimn().gety() * ellipPileDimn().getz()
	      *(2.0/3*ellipPileDimn().getx()-lowint);
    }
    return val;
}


void assembly::ellipPileUpdate(){
    for(std::list<particle*>::iterator it=ParticleList.begin();it!=ParticleList.end();++it){
	if ((*it)->getType()==3) {
	    (*it)->curr_velocity.setx(0);	
	    (*it)->curr_velocity.sety(0);
	    (*it)->curr_velocity.setz(-PILE_RATE);
	    (*it)->curr_position = (*it)->prev_position + (*it)->curr_velocity*TIMESTEP;
	}
    }
}


REAL assembly::getTransEnergy() const{
    REAL engy=0;
    std::list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it){
	if ((*it)->getType()==0)
	    engy+=(*it)->getTransEnergy();
    }
    return engy;
}


REAL assembly::getRotatEnergy() const{
    REAL engy=0;
    std::list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it){
	if ((*it)->getType()==0)
	    engy+=(*it)->getRotatEnergy();
    }
    return engy;
}


REAL assembly::getKinetEnergy() const{
    REAL engy=0;
    std::list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it){
	if ((*it)->getType()==0)
	    engy+=(*it)->getKinetEnergy();
    }
    return engy;
}


REAL assembly::getPotenEnergy(REAL ref) const{
    REAL engy=0;
    std::list<particle*>::const_iterator it;
    for(it=ParticleList.begin();it!=ParticleList.end();++it){
	if ((*it)->getType()==0)
	    engy+=(*it)->getPotenEnergy(ref);
    }
    return engy;
}


void assembly::clearForce(){
    for(std::list<particle*>::iterator it=ParticleList.begin();it!=ParticleList.end();++it){
	(*it)->clearForce();
    }
}


void assembly::flexiBoundaryForceZero(){
    for(std::list<particle*>::iterator it=ParticleList.begin();it!=ParticleList.end();++it){
	(*it)->flb_force=0;
	(*it)->flb_moment=0;
    }
}


void assembly::initFBForce(){
    for(std::list<particle*>::iterator it=ParticleList.begin();it!=ParticleList.end();++it){
	(*it)->force+=(*it)->flb_force;
	(*it)->moment+=(*it)->flb_moment;
    }
}


void assembly::internalForce(REAL& avgnm, REAL& avgsh){
    avgnm=0;
    avgsh=0;

    int totalcntct = ContactList.size();
    if(totalcntct==0){
	avgnm = 0;
	avgsh = 0;
    }
    else{
	std::list<CONTACT>::iterator it;
	for (it=ContactList.begin();it!=ContactList.end();++it)
	    it->checkinPreTgt(CntTgtVec); // checkin previous tangential force and displacment    
	
	CntTgtVec.clear(); // CntTgtVec must be cleared before filling in new values.

#ifdef TIME_PROFILE
	gettimeofday(&timep1,NULL); 
#endif 
	for (it=ContactList.begin();it!=ContactList.end();++it){
            bool exceed = false;
	    it->contactForce(exceed);           // cannot be parallelized as it may change a particle's force simultaneously.
	    it->checkoutTgt(CntTgtVec);   // checkout current tangential force and displacment
	    avgnm += it->getNormalForce();
	    avgsh += it->getTgtForce();
#ifdef DEBUG
	    if (exceed) {
	      char stepsstr[7];
	      char stepsfp[50];
	      sprintf(stepsstr, "%06d", g_iteration);
	      strcpy(stepsfp,"particle_");
	      strcat(stepsfp, stepsstr);
	      printParticle(stepsfp);
	    }
#endif
	}
	avgnm /= totalcntct;
	avgsh /= totalcntct;

#ifdef TIME_PROFILE
	gettimeofday(&timep2,NULL);
	g_debuginf<<setw(OWID)<<timediffsec()<<endl; 
#endif

    }
}


void assembly::updateParticle(){
    for(std::list<particle*>::iterator it=ParticleList.begin();it!=ParticleList.end();++it){
	(*it)->update();
    }
}


void assembly::createRigidBoundary(std::ifstream &ifs){
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


void assembly::createFlexiBoundary(std::ifstream &ifs){
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
    std::ifstream ifs(str);
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


void assembly::printBoundary(const char* str) const
{
    std::ofstream ofs(str);
    if(!ofs) {
	cout<<"stream error in printBoundary!"<<endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);

    ofs<<setw(OWID)<<BdryType
       <<setw(OWID)<<RgdBdryNum<<endl;
    
    std::list<RGDBDRY*>::const_iterator rt;
    for(rt=RBList.begin();rt!=RBList.end();++rt)
	(*rt)->disp(ofs);
    ofs<<endl;

    ofs.close();
}


void assembly::findParticleOnBoundary(){
    std::list<RGDBDRY*>::iterator rt;
    std::list<FLBBDRY*>::iterator ft;
    for(rt=RBList.begin();rt!=RBList.end();++rt)
	(*rt)->findParticleOnBoundary(ParticleList);
    for(ft=FBList.begin();ft!=FBList.end();++ft)
	(*ft)->findParticleOnBoundary(ParticleList);
}


void assembly::findParticleOnLine(){
    std::list<FLBBDRY*>::iterator ft;
    for(ft=FBList.begin();ft!=FBList.end();++ft)
	(*ft)->findParticleOnLine();
}


void assembly::createFlbNet(){
    std::list<FLBBDRY*>::iterator ft;
    for(ft=FBList.begin();ft!=FBList.end();++ft)
	(*ft)->createFlbNet();
}


void assembly::rigidBoundaryForce(){
  std::list<RGDBDRY*>::iterator rt;
  for(rt=RBList.begin();rt!=RBList.end();++rt)
    (*rt)->rigidBF(BdryTgtMap);

  /*
  vector<boundarytgt>::iterator it;
  std::list<RGDBDRY*>::iterator rt;

  for(rt=RBList.begin();rt!=RBList.end();++rt){	
    (*rt)->rigidBF(BdryTgtMap);
    for (it=BdryTgtMap[(*rt)->bdry_id].begin();it!=BdryTgtMap[(*rt)->bdry_id].end();++it){
      g_debuginf<<setw(OWID)<<g_iteration
		<<setw(OWID)<<(*rt)->bdry_id
		<<setw(OWID)<<BdryTgtMap[(*rt)->bdry_id].size()
		<<setw(OWID)<<it->TgtForce.getx()
		<<setw(OWID)<<it->TgtForce.gety()
		<<setw(OWID)<<it->TgtForce.getz()
		<<endl;
      //<<setw(OWID)<<it->TgtPeak<<endl;
    }
  }
  */
}


void assembly::rigidBoundaryForce(REAL penetr[],int cntnum[]){
  std::list<RGDBDRY*>::iterator rt;
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
    std::list<FLBBDRY*>::iterator ft;
    for(ft=FBList.begin();ft!=FBList.end();++ft)
	(*ft)->flxbBF();
}


vec assembly::getNormalForce(int bdry) const{
    std::list<RGDBDRY*>::const_iterator it;
    for(it=RBList.begin();it!=RBList.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getNormalForce();
    }
    return 0;
}


vec assembly::getShearForce(int bdry) const{
    std::list<RGDBDRY*>::const_iterator it;
    for(it=RBList.begin();it!=RBList.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getShearForce();
    }
    return 0;
}


REAL assembly::getAvgNormal(int bdry) const{
    std::list<RGDBDRY*>::const_iterator it;
    for(it=RBList.begin();it!=RBList.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getAvgNormal();
    }
    return 0;
}


vec assembly::getApt(int bdry) const{
    std::list<RGDBDRY*>::const_iterator it;
    for(it=RBList.begin();it!=RBList.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getApt();
    }
    return 0;
}


vec assembly::getDirc(int bdry) const{
    std::list<RGDBDRY*>::const_iterator it;
    for(it=RBList.begin();it!=RBList.end();++it){
	if((*it)->getBdryID()==bdry)
	    return (*it)->getDirc();
    }
    return 0;
}


REAL assembly::getArea(int n) const{
    std::list<RGDBDRY*>::const_iterator it;
    for(it=RBList.begin();it!=RBList.end();++it){
	if((*it)->getBdryID()==n)
	    return (*it)->area;
    }
    return 0;
}


void assembly::setArea(int n, REAL a){
    std::list<RGDBDRY*>::iterator it;
    for(it=RBList.begin();it!=RBList.end();++it){
	if((*it)->getBdryID()==n)
	    (*it)->area=a;
    }
}


REAL assembly::getAverageRigidPressure() const{
    std::list<RGDBDRY*>::const_iterator rt;
    REAL avgpres=0;
    for(rt=RBList.begin();rt!=RBList.end();++rt)
	avgpres+=vfabs((*rt)->getNormalForce())/(*rt)->getArea();
    return avgpres/=RgdBdryNum;
}


// only update CoefOfLimits[0] for specified boundaries
void assembly::updateRB(int bn[], UPDATECTL rbctl[], int num){
    for(int i=0;i<num;i++){
	for(std::list<RGDBDRY*>::iterator rt=RBList.begin();rt!=RBList.end();++rt){
	    if((*rt)->getBdryID()==bn[i]){
		(*rt)->update(rbctl[i]);
		break;
	    }
	}
    }
}


// update CoefOfLimits[1,2,3,4] for all 6 boundaries
void assembly::updateRB6(){
    for(std::list<RGDBDRY*>::iterator rt=RBList.begin();rt!=RBList.end();++rt){
	if((*rt)->getBdryID()==1 || (*rt)->getBdryID()==3){
	    for(std::list<RGDBDRY*>::iterator lt=RBList.begin();lt!=RBList.end();++lt){
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
	    for(std::list<RGDBDRY*>::iterator lt=RBList.begin();lt!=RBList.end();++lt){
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
	    for(std::list<RGDBDRY*>::iterator lt=RBList.begin();lt!=RBList.end();++lt){
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
    for(std::list<RGDBDRY*>::iterator rt=RBList.begin();rt!=RBList.end();++rt){
	if((*rt)->getBdryID()==7 || (*rt)->getBdryID()==9 ){
	    for(std::list<RGDBDRY*>::iterator lt=RBList.begin();lt!=RBList.end();++lt){
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
	    for(std::list<RGDBDRY*>::iterator lt=RBList.begin();lt!=RBList.end();++lt){
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
    std::list<FLBBDRY*>::iterator ft;
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
void assembly::deposit_RgdBdry(rectangle& container,
			       gradation& grad,
			       int   freetype,
			       int   total_steps,  
			       int   snapshots,
			       int   interval,
			       REAL  height,
			       const char* iniptclfile,   
			       const char* inibdryfile,
			       const char* particlefile, 
			       const char* contactfile,
			       const char* progressfile, 
			       const char* creparticle,
			       const char* creboundary,
			       const char* debugfile)
{
  this->container = container;
  generate_new(container, grad, iniptclfile, freetype, height); 
  // 3.0 for uniform size of (2.5e-3,*0.8,*0.6); if not uniform, it may be larger, say, 4.5
  // 3.0 for uniform spheres of 2.5e-3
  
  setBoundary(5, container, inibdryfile);
  
  deposit(total_steps,        // total_steps
	  snapshots,          // number of snapshots
	  interval,           // print interval
	  iniptclfile,        // input file, initial particles
	  inibdryfile,        // input file, initial boundaries
	  particlefile,       // output file, resulted particles, including snapshots 
	  contactfile,        // output file, resulted contacts, including snapshots 
	  progressfile,       // output file, statistical info
	  debugfile);         // output file, debug info
  
  setBoundary(6,                  
	      container,      // specimen dimension
	      "trm_boundary");// output file, containing boundaries info
  
  trim("dep_particle_end",    // input file, particles to be trimmed
       "trm_boundary",  
       creparticle,
       creboundary); 
}


// freetype:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
// ht- how many times of size would be the floating height
void assembly::generate_new(rectangle& container,
			    gradation& grad,
			    const char* particlefile,
			    int freetype,
			    REAL vht)
{
  REAL x,y,z;
  particle* newptcl;
  TotalNum = 0;
  REAL dimx     = container.getDimx();
  REAL dimy     = container.getDimy();
  REAL dimz     = container.getDimz();
  REAL diameter = grad.getMaxPtclRadius()*2.0;

  REAL x1 = container.getV1().getx() + diameter/2.0;
  REAL y1 = container.getV1().gety() + diameter/2.0;
  REAL z1 = container.getV1().getz() + diameter/2.0;
  REAL x2 = container.getV2().getx() - diameter/2.0;
  REAL y2 = container.getV2().gety() - diameter/2.0;
  REAL z2 = container.getV2().getz() - diameter/2.0;
  REAL x0 = container.getCenter().getx();
  REAL y0 = container.getCenter().gety();
  REAL z0 = container.getCenter().getz();

  if (freetype == 0) {      // just one free particle
    newptcl = new particle(TotalNum+1, 0, vec(x0,y0,z0), grad);
    ParticleList.push_back(newptcl);
    TotalNum++;
  }
  else if (freetype == 1) { // a horizontal layer of free particles
    z = container.getV2().getz();
    for (x = x1; x - x2 < NUMZERO; x += diameter)
      for (y = y1; y - y2 < NUMZERO; y += diameter) {
	newptcl = new particle(TotalNum+1, 0, vec(x,y,z), grad);
	ParticleList.push_back(newptcl);
	TotalNum++;
      }
  }
  else if (freetype == 2) { // multiple layers of free particles
    REAL offset = 0;
    if (grad.getPtclRatioBA() == 1.0 && grad.getPtclRatioCA() == 1.0)
      offset = diameter*0.05;
    
    for (z = z1; z < z1 + dimz*vht; z += diameter) {
      for (x = x1 + offset; x - x2 < NUMZERO; x += diameter)
	for (y = y1 + offset; y - y2 < NUMZERO; y += diameter) {
	  newptcl = new particle(TotalNum+1, 0, vec(x,y,z), grad);
	  ParticleList.push_back(newptcl);
	  TotalNum++;
	}	
      offset *= -1;
    }
  }
  
  printParticle(particlefile);
  
}


/*
// freetype:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
// ht- how many times of size would be the floating height
void assembly::generate(gradation& grad,
			const char* particlefile,
			int freetype,
			REAL ht)
{
    REAL x,y,z;
    particle* newptcl;
    TotalNum = 0;
    REAL est =1.02;
    int grid=9;  
    // grid: dimension of free particle array.
    // 7 - small dimn container
    // 9 - medium dimn container 
    // 11- large dimn container 

    REAL dimn=grad.dimn;
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
	REAL offset=0; // 0 for ellipsoids; dimn/2/5/5 for spheres
	if (grad.ratio_ba==1.0 && grad.ratio_ca==1.0)
	    offset = dimn/2/5/5;
	REAL z0 = -dimn/2*9/10 ;// dimn/2;
	for (z=z0; z<z0 + dimn*ht; z+=dimn/2/5) {
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
*/
 /*
// create a specimen from discreate particles through floating and then gravitation,
// boundaries are composed of fixed particles.
void assembly::deposit_PtclBdry(gradation& grad,
				int   freetype,
				REAL rsize,
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
	container.setCenter(vec(0,0,0));
	container.setDimx(grad.dimn);
	container.setDimy(grad.dimn);
	container.setDimz(grad.dimn);
	
	generate_p(grad, iniptclfile, freetype, rsize, 4.0);
	deposit_p(total_steps,        // total_steps
		  snapshots,          // number of snapshots
		  interval,           // print interval
		  grad.dimn,          // dimension of particle-composed-boundary
		  rsize,              // relative container size
		  iniptclfile,        // input file, initial particles
		  particlefile,       // output file, resulted particles, including snapshots 
		  contactfile,        // output file, resulted contacts, including snapshots 
		  progressfile,       // output file, statistical info
		  debugfile);         // output file, debug info
    }
}
 */
  /*
// freetype:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
// ht- how many times of size would be the floating height
void assembly::generate_p(gradation&  grad,
			 const char* particlefile,
			 int freetype,
			 REAL rsize,
			 REAL ht)
{
    REAL x,y,z;
    particle* newptcl;
    TotalNum = 0;
    REAL wall=2.2; // wall - wall height; ht - free particle height
    REAL est =1.02;
    int grid=static_cast<int> (nearbyint(rsize*10)-1);  

    // grid: dimension of free particle array.
    // 7 - small dimn container
    // 9 - medium dimn container 
    // 11- large dimn container 

    REAL dimn=grad.dimn;
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
*/

void assembly::scale_PtclBdry(int   total_steps,  
			      int   snapshots,
			      int   interval,
			      REAL dimn,
			      REAL rsize,
			      const char* iniptclfile,   
			      const char* particlefile, 
			      const char* contactfile,
			      const char* progressfile, 
			      const char* debugfile)
{
    deposit_p(total_steps,        // total_steps
	      snapshots,          // number of snapshots
	      interval,           // print interval
	      dimn,               // dimension of particle-composed-boundary
	      rsize,              // relative container size
	      iniptclfile,        // input file, initial particles
	      particlefile,       // output file, resulted particles, including snapshots 
	      contactfile,        // output file, resulted contacts, including snapshots 
	      progressfile,       // output file, statistical info
	      debugfile);         // output file, debug info
}

/*
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
	    interval,           // print interval
	    iniptclfile,        // input file, initial particles
	    initboundary,       // input file, boundaries
	    particlefile,       // output file, resulted particles, including snapshots 
	    contactfile,        // output file, resulted contacts, including snapshots 
	    progressfile,       // output file, statistical info
	    debugfile);         // output file, debug info
}
*/
  
void assembly::setBoundary(int   bdrynum,
			   rectangle& container,
			   const char* boundaryfile)
{
  std::ofstream ofs(boundaryfile);
  if(!ofs) { cout<<"stream error!"<<endl; exit(-1);}
  REAL x1,x2,y1,y2,z1,z2,x0,y0,z0;
  x1 = container.getV1().getx();
  y1 = container.getV1().gety();
  z1 = container.getV1().getz();
  x2 = container.getV2().getx();
  y2 = container.getV2().gety();
  z2 = container.getV2().getz();
  x0 = container.getCenter().getx();
  y0 = container.getCenter().gety();
  z0 = container.getCenter().getz();

  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs<<setw(OWID)<<0
     <<setw(OWID)<<bdrynum<<endl<<endl;
  
  if (bdrynum == 1){   // only a bottom boundary
    ofs<<setw(OWID)<<1<<endl
       <<setw(OWID)<<6
       <<setw(OWID)<<1
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<0
       <<setw(OWID)<<-1
       <<setw(OWID)<<0
       <<setw(OWID)<<0
       <<setw(OWID)<<z1
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl;
    
  }
  else if (bdrynum == 5){ // no top boundary
    // boundary 1
    ofs<<setw(OWID)<<1<<endl
       <<setw(OWID)<<1
       <<setw(OWID)<<4
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<0 
       <<setw(OWID)<<x2
       <<setw(OWID)<<y0
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<-1
       <<setw(OWID)<<0 
       <<setw(OWID)<<x0
       <<setw(OWID)<<y1
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<x0
       <<setw(OWID)<<y2
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<0
       <<setw(OWID)<<-1
       <<setw(OWID)<<x0     
       <<setw(OWID)<<y0
       <<setw(OWID)<<z1
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl<<endl
      
      // boundary 2
       <<setw(OWID)<<1<<endl
       <<setw(OWID)<<2
       <<setw(OWID)<<4
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<x0     
       <<setw(OWID)<<y2
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<0 
       <<setw(OWID)<<x2
       <<setw(OWID)<<y0      
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<-1
       <<setw(OWID)<<0
       <<setw(OWID)<<0 
       <<setw(OWID)<<x1
       <<setw(OWID)<<y0     
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<0
       <<setw(OWID)<<-1
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y0     
       <<setw(OWID)<<z1
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl<<endl
      
      // boundary 3
       <<setw(OWID)<<1<<endl
       <<setw(OWID)<<3
       <<setw(OWID)<<4
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<-1
       <<setw(OWID)<<0
       <<setw(OWID)<<0 
       <<setw(OWID)<<x1
       <<setw(OWID)<<y0
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<-1
       <<setw(OWID)<<0 
       <<setw(OWID)<<x0     
       <<setw(OWID)<<y1
       <<setw(OWID)<<y0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y2
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<0
       <<setw(OWID)<<-1
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y0     
       <<setw(OWID)<<z1
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl<<endl
      
      // boundary 4
       <<setw(OWID)<<1<<endl
       <<setw(OWID)<<4
       <<setw(OWID)<<4
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<-1
       <<setw(OWID)<<0 
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y1
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<0 
       <<setw(OWID)<<x2
       <<setw(OWID)<<y0
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<-1
       <<setw(OWID)<<0
       <<setw(OWID)<<0 
       <<setw(OWID)<<x1
       <<setw(OWID)<<y0
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<0
       <<setw(OWID)<<-1
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y0     
       <<setw(OWID)<<z1
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl<<endl
      
      // boundary 6
       <<setw(OWID)<<1<<endl
       <<setw(OWID)<<6
       <<setw(OWID)<<5
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<0
       <<setw(OWID)<<-1
       <<setw(OWID)<<x0
       <<setw(OWID)<<y0
       <<setw(OWID)<<z1
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<0
       <<setw(OWID)<<x2
       <<setw(OWID)<<y0
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<-1
       <<setw(OWID)<<0
       <<setw(OWID)<<0
       <<setw(OWID)<<x1 
       <<setw(OWID)<<y0
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y2
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<-1
       <<setw(OWID)<<0
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y1
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl;
  }
  else if (bdrynum == 6){ // all 6 boundaries
    // boundary 1
    ofs<<setw(OWID)<<1<<endl
       <<setw(OWID)<<1
       <<setw(OWID)<<5
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<0     
       <<setw(OWID)<<x2
       <<setw(OWID)<<y0
       <<setw(OWID)<<z0     
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<-1
       <<setw(OWID)<<0
       <<setw(OWID)<<x0    
       <<setw(OWID)<<y1
       <<setw(OWID)<<z0     
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y2
       <<setw(OWID)<<z0     
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<0
       <<setw(OWID)<<1
       <<setw(OWID)<<x0     
       <<setw(OWID)<<y0    
       <<setw(OWID)<<z2 
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<0 
       <<setw(OWID)<<-1
       <<setw(OWID)<<x0     
       <<setw(OWID)<<y0     
       <<setw(OWID)<<z1
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl<<endl
      
      // boundary 2
       <<setw(OWID)<<1<<endl
       <<setw(OWID)<<2
       <<setw(OWID)<<5
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<1
       <<setw(OWID)<<0     
       <<setw(OWID)<<x0    
       <<setw(OWID)<<y2
       <<setw(OWID)<<z0     
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<0
       <<setw(OWID)<<x2
       <<setw(OWID)<<y0
       <<setw(OWID)<<z0     
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<-1
       <<setw(OWID)<<0
       <<setw(OWID)<<0
       <<setw(OWID)<<x1 
       <<setw(OWID)<<y0
       <<setw(OWID)<<z0     
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<0
       <<setw(OWID)<<1
       <<setw(OWID)<<x0     
       <<setw(OWID)<<y0    
       <<setw(OWID)<<z2 
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<0 
       <<setw(OWID)<<-1
       <<setw(OWID)<<x0     
       <<setw(OWID)<<y0      
       <<setw(OWID)<<z1
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl<<endl
      
      // boundary 3
       <<setw(OWID)<<1<<endl
       <<setw(OWID)<<3
       <<setw(OWID)<<5
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<-1
       <<setw(OWID)<<0
       <<setw(OWID)<<0     
       <<setw(OWID)<<x1
       <<setw(OWID)<<y0
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<-1
       <<setw(OWID)<<0
       <<setw(OWID)<<x0     
       <<setw(OWID)<<y1
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0  
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<x0       
       <<setw(OWID)<<y2
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<0
       <<setw(OWID)<<1
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y0     
       <<setw(OWID)<<z2 
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<0 
       <<setw(OWID)<<-1
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y0      
       <<setw(OWID)<<z1
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl<<endl
      
      // boundary 4
       <<setw(OWID)<<1<<endl
       <<setw(OWID)<<4
       <<setw(OWID)<<5
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<-1
       <<setw(OWID)<<0     
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y1
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<0
       <<setw(OWID)<<x2
       <<setw(OWID)<<y0
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<-1 
       <<setw(OWID)<<0
       <<setw(OWID)<<0
       <<setw(OWID)<<x1 
       <<setw(OWID)<<y0
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<0
       <<setw(OWID)<<1
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y0     
       <<setw(OWID)<<z2 
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<0 
       <<setw(OWID)<<-1
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y0      
       <<setw(OWID)<<z1
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl<<endl
      
      // boundary 5
       <<setw(OWID)<<1<<endl
       <<setw(OWID)<<5
       <<setw(OWID)<<5
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<0
       <<setw(OWID)<<1     
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y0
       <<setw(OWID)<<z2 
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<0
       <<setw(OWID)<<x2
       <<setw(OWID)<<y0
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<-1 
       <<setw(OWID)<<0
       <<setw(OWID)<<0
       <<setw(OWID)<<x1 
       <<setw(OWID)<<y0
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y2
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<-1
       <<setw(OWID)<<0
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y1
       <<setw(OWID)<<z0
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl<<endl
      
      // boundary 6
       <<setw(OWID)<<1<<endl
       <<setw(OWID)<<6
       <<setw(OWID)<<5
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<0
       <<setw(OWID)<<-1    
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y0
       <<setw(OWID)<<z1
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<1
       <<setw(OWID)<<0 
       <<setw(OWID)<<0
       <<setw(OWID)<<x2
       <<setw(OWID)<<y0
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<-1 
       <<setw(OWID)<<0
       <<setw(OWID)<<0
       <<setw(OWID)<<x1 
       <<setw(OWID)<<y0
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y2
       <<setw(OWID)<<z0      
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl
      
       <<setw(OWID)<<1
       <<setw(OWID)<<0
       <<setw(OWID)<<-1
       <<setw(OWID)<<0 
       <<setw(OWID)<<x0      
       <<setw(OWID)<<y1
       <<setw(OWID)<<z0
       <<setw(OWID)<<0
       <<setw(OWID)<<0<<endl<<endl;
  }
  
  ofs.close();
}


void assembly::trim(const char* iniptclfile,
		    const char* inibdryfile,
		    const char* particlefile,
		    const char* boundaryfile)
{
  createSample(iniptclfile);
  createBoundary(inibdryfile);
  
  std::list<particle*>::iterator itr,itp;
  vec center;
  REAL mass = 0;
  
  if (RgdBdryNum == 1) {
  }
  else if (RgdBdryNum == 4) {
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    container.setDimx(W0); 
    container.setDimy(L0); 
    
    for(itr=ParticleList.begin();itr!=ParticleList.end();++itr){
      center=(*itr)->getCurrPosition();
      if(fabs(center.getx()) >= L0/2 ||
	 fabs(center.gety()) >= W0/2 )
	{
	  itp = itr;
	  --itr;
	  delete (*itp); // release memory
	  ParticleList.erase(itp); 
	}
    }
    
  }
  else if (RgdBdryNum == 5) {
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = -getApt(6).getz()*4;
    container.setDimx(W0); 
    container.setDimy(L0); 
    container.setDimz(H0);
    container.setCenter(0);
    Volume = W0*L0*H0;
    
    for(itr=ParticleList.begin();itr!=ParticleList.end();++itr){
      center=(*itr)->getCurrPosition();
      if(fabs(center.getx()) >= L0/2 ||
	 fabs(center.gety()) >= W0/2 ||
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
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    container.setDimx(W0); 
    container.setDimy(L0); 
    container.setDimz(H0);
    container.setCenter(0);
    Volume = W0*L0*H0;
    
    for(itr=ParticleList.begin();itr!=ParticleList.end();++itr){
      center=(*itr)->getCurrPosition();
      if(fabs(center.getx()) >= L0/2 ||
	 fabs(center.gety()) >= W0/2 ||
	 fabs(center.getz()) >= H0/2 )
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
  
  printParticle(particlefile);
  printBoundary(boundaryfile);
}


void assembly::TrimPtclBdryByHeight(REAL height,
			    const char* iniptclfile,
			    const char* particlefile)
{
    createSample(iniptclfile);

    std::list<particle*>::iterator itr, itp;
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf.precision(OPREC);
    progressinf<<setw(OWID)<<"iteration"
	       <<setw(OWID)<<"possible"
	       <<setw(OWID)<<"actual"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"translational"
	       <<setw(OWID)<<"rotational"
	       <<setw(OWID)<<"kinetic"
	       <<setw(OWID)<<"potential"
	       <<setw(OWID)<<"total"
	       <<setw(OWID)<<"void"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"coordination"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"vibra"
	       <<setw(OWID)<<"impact"
	       <<setw(OWID)<<"wall-clock" << endl
	       <<setw(OWID)<<"number"
	       <<setw(OWID)<<"contacts"
	       <<setw(OWID)<<"contacts"
	       <<setw(OWID)<<"penetration"
	       <<setw(OWID)<<"contact_normal"
	       <<setw(OWID)<<"contact_tangt"
	       <<setw(OWID)<<"velocity"
	       <<setw(OWID)<<"omga"
	       <<setw(OWID)<<"force"
	       <<setw(OWID)<<"moment"
	       <<setw(OWID)<<"energy"
	       <<setw(OWID)<<"energy"
	       <<setw(OWID)<<"energy"
	       <<setw(OWID)<<"energy"
	       <<setw(OWID)<<"energy"
	       <<setw(OWID)<<"ratio"
	       <<setw(OWID)<<"porosity"
	       <<setw(OWID)<<"number"
	       <<setw(OWID)<<"density"
	       <<setw(OWID)<<"sigma1_1"
	       <<setw(OWID)<<"sigma1_2"
	       <<setw(OWID)<<"sigma2_1"
	       <<setw(OWID)<<"sigma2_2"
	       <<setw(OWID)<<"sigma3_1"
	       <<setw(OWID)<<"sigma3_2"
	       <<setw(OWID)<<"mean_stress"
	       <<setw(OWID)<<"width"
	       <<setw(OWID)<<"length"
	       <<setw(OWID)<<"height"
	       <<setw(OWID)<<"volume"
	       <<setw(OWID)<<"epsilon_w"
	       <<setw(OWID)<<"epsilon_l"
	       <<setw(OWID)<<"epsilon_h"
	       <<setw(OWID)<<"epsilon_v"
	       <<setw(OWID)<<"t_step"
	       <<setw(OWID)<<"t_step"
	       <<setw(OWID)<<"time" << endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1); }
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from existing files.
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries.

    // pre_3: define variables used in iterations.
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
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
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<t1
		       <<setw(OWID)<<t2
		       <<setw(OWID)<<(t1+t2)
		       <<setw(OWID)<<t3
		       <<setw(OWID)<<(t1+t2+t3)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/TotalNum
		       <<setw(OWID)<<"0"
		       <<setw(OWID)<<"0"
		       <<setw(OWID)<<"0"
		       <<setw(OWID)<<"0"
		       <<setw(OWID)<<"0"
		       <<setw(OWID)<<"0"
		       <<setw(OWID)<<"0"
		       <<setw(OWID)<<"0"
		       <<setw(OWID)<<"0"
		       <<setw(OWID)<<"0"
		       <<setw(OWID)<<"0"
		       <<setw(OWID)<<"0"
		       <<setw(OWID)<<"0"
		       <<setw(OWID)<<"0"
		       <<setw(OWID)<<"0"
		       <<setw(OWID)<<"0"
	               <<setw(OWID)<<getVibraTimeStep()
	               <<setw(OWID)<<getImpactTimeStep()
		       <<setw(OWID)<<timediffsec(timew1,timew2)
		       <<endl;

	    /*
	    g_debuginf<<setw(OWID)<<g_iteration
		      <<setw(OWID)<<bdry_penetr[1]
		      <<setw(OWID)<<bdry_penetr[2]
		      <<setw(OWID)<<bdry_penetr[3]
		      <<setw(OWID)<<bdry_penetr[4]
		      <<setw(OWID)<<bdry_penetr[6]
		      <<setw(OWID)<<bdry_cntnum[1]
		      <<setw(OWID)<<bdry_cntnum[2]
		      <<setw(OWID)<<bdry_cntnum[3]
		      <<setw(OWID)<<bdry_cntnum[4]
		      <<setw(OWID)<<bdry_cntnum[6]
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
			 REAL dimn,
			 REAL rsize,
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"epsilon_v"<<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1); }
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from existing files.
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 

    // pre_3: define variables used in iterations.
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    REAL void_ratio=0;

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
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<t1
		       <<setw(OWID)<<t2
		       <<setw(OWID)<<(t1+t2)
		       <<setw(OWID)<<t3
		       <<setw(OWID)<<(t1+t2+t3)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*getActualCntctNum()/TotalNum
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"epsilon_v"<<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1); }
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from existing files.
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries.

    // pre_3: define variables used in iterations.
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];

    int         mid[2]={1,3};    // boundary 1 and 3
    UPDATECTL   midctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
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
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<t1
		       <<setw(OWID)<<t2
		       <<setw(OWID)<<(t1+t2)
		       <<setw(OWID)<<t3
		       <<setw(OWID)<<(t1+t2+t3)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    g_debuginf<<setw(OWID)<<g_iteration
		      <<setw(OWID)<<bdry_penetr[1]
		      <<setw(OWID)<<bdry_penetr[2]
		      <<setw(OWID)<<bdry_penetr[3]
		      <<setw(OWID)<<bdry_penetr[4]
		      <<setw(OWID)<<bdry_penetr[6]
		      <<setw(OWID)<<bdry_cntnum[1]
		      <<setw(OWID)<<bdry_cntnum[2]
		      <<setw(OWID)<<bdry_cntnum[3]
		      <<setw(OWID)<<bdry_cntnum[4]
		      <<setw(OWID)<<bdry_cntnum[6]
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
			 REAL sigma,			  
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3: define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         mid[2]={1,3};    // boundary 1 and 3
    int         max[2]={2,4};    // boundary 2 and 4
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   midctl[2];
    UPDATECTL   maxctl[2];
    UPDATECTL   minctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
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
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
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
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    g_debuginf<<setw(OWID)<<g_iteration
		      <<setw(OWID)<<bdry_penetr[1]
		      <<setw(OWID)<<bdry_penetr[2]
		      <<setw(OWID)<<bdry_penetr[3]
		      <<setw(OWID)<<bdry_penetr[4]
		      <<setw(OWID)<<bdry_penetr[5]
		      <<setw(OWID)<<bdry_penetr[6]
		      <<setw(OWID)<<bdry_cntnum[1]
		      <<setw(OWID)<<bdry_cntnum[2]
		      <<setw(OWID)<<bdry_cntnum[3]
		      <<setw(OWID)<<bdry_cntnum[4]
		      <<setw(OWID)<<bdry_cntnum[5]
		      <<setw(OWID)<<bdry_cntnum[6]
		      <<endl;
	}

	// 8. loop break condition
	if (   fabs(sigma1_1-sigma)/sigma < STRESS_ERROR && fabs(sigma1_2-sigma)/sigma < STRESS_ERROR
	    && fabs(sigma2_1-sigma)/sigma < STRESS_ERROR && fabs(sigma2_2-sigma)/sigma < STRESS_ERROR
	    && fabs(sigma3_1-sigma)/sigma < STRESS_ERROR && fabs(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()    
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()    
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
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
			 REAL sigma_a,
			 REAL sigma_b,
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3: define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         mid[2]={1,3};    // boundary 1 and 3
    int         max[2]={2,4};    // boundary 2 and 4
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   midctl[2];
    UPDATECTL   maxctl[2];
    UPDATECTL   minctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    REAL sigma=sigma_a;
    REAL sigma_inc=(sigma_b-sigma_a)/sigma_division;

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
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
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
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    g_debuginf<<setw(OWID)<<g_iteration
		      <<setw(OWID)<<bdry_penetr[1]
		      <<setw(OWID)<<bdry_penetr[2]
		      <<setw(OWID)<<bdry_penetr[3]
		      <<setw(OWID)<<bdry_penetr[4]
		      <<setw(OWID)<<bdry_penetr[5]
		      <<setw(OWID)<<bdry_penetr[6]
		      <<setw(OWID)<<bdry_cntnum[1]
		      <<setw(OWID)<<bdry_cntnum[2]
		      <<setw(OWID)<<bdry_cntnum[3]
		      <<setw(OWID)<<bdry_cntnum[4]
		      <<setw(OWID)<<bdry_cntnum[5]
		      <<setw(OWID)<<bdry_cntnum[6]
		      <<endl;
	}

	// 8. find the balanced status and increase confining pressure
	if (   fabs(sigma1_1-sigma)/sigma < STRESS_ERROR && fabs(sigma1_2-sigma)/sigma < STRESS_ERROR
	    && fabs(sigma2_1-sigma)/sigma < STRESS_ERROR && fabs(sigma2_2-sigma)/sigma < STRESS_ERROR
	    && fabs(sigma3_1-sigma)/sigma < STRESS_ERROR && fabs(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()    
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    sigma += sigma_inc;
	}

	// 9. loop break condition
	if (   fabs(sigma1_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma1_2-sigma_b)/sigma_b < STRESS_ERROR
	    && fabs(sigma2_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma2_2-sigma_b)/sigma_b < STRESS_ERROR
	    && fabs(sigma3_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma3_2-sigma_b)/sigma_b < STRESS_ERROR ) {
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()    
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
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
			 REAL sigma_values[],  
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3: define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         mid[2]={1,3};    // boundary 1 and 3
    int         max[2]={2,4};    // boundary 2 and 4
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   midctl[2];
    UPDATECTL   maxctl[2];
    UPDATECTL   minctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    int  i=0;
    REAL sigma=sigma_values[i];
    REAL sigma_inc=(sigma_values[i+1]-sigma_values[i])/sigma_division;
    REAL sigma_b=sigma_values[sigma_points-1];

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
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
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
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    g_debuginf<<setw(OWID)<<g_iteration
		      <<setw(OWID)<<bdry_penetr[1]
		      <<setw(OWID)<<bdry_penetr[2]
		      <<setw(OWID)<<bdry_penetr[3]
		      <<setw(OWID)<<bdry_penetr[4]
		      <<setw(OWID)<<bdry_penetr[5]
		      <<setw(OWID)<<bdry_penetr[6]
		      <<setw(OWID)<<bdry_cntnum[1]
		      <<setw(OWID)<<bdry_cntnum[2]
		      <<setw(OWID)<<bdry_cntnum[3]
		      <<setw(OWID)<<bdry_cntnum[4]
		      <<setw(OWID)<<bdry_cntnum[5]
		      <<setw(OWID)<<bdry_cntnum[6]
		      <<endl;
	}

	// 8. find the balanced status and increase confining pressure
	if (   fabs(sigma1_1-sigma)/sigma < STRESS_ERROR && fabs(sigma1_2-sigma)/sigma < STRESS_ERROR
	    && fabs(sigma2_1-sigma)/sigma < STRESS_ERROR && fabs(sigma2_2-sigma)/sigma < STRESS_ERROR
	    && fabs(sigma3_1-sigma)/sigma < STRESS_ERROR && fabs(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()    
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
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
	if (   fabs(sigma1_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma1_2-sigma_b)/sigma_b < STRESS_ERROR
	    && fabs(sigma2_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma2_2-sigma_b)/sigma_b < STRESS_ERROR
	    && fabs(sigma3_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma3_2-sigma_b)/sigma_b < STRESS_ERROR ) {
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()    
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
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
			REAL sigma_3,     
			REAL sigma_1,    
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries
 
    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];

    int min[2]={5,6};    // minimum stress acting on boundary 5 and 6
    UPDATECTL minctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    REAL sigma=sigma_3;
    REAL sigma_inc=(sigma_1-sigma_3)/sigma_division;

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
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
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
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	}

	// 8. find balanced status of odometer compression
	if (fabs(sigma3_1-sigma)/sigma < STRESS_ERROR && fabs(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()    
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    sigma += sigma_inc;
	}

	// 9. loop break condition
	if (fabs(sigma3_1-sigma_1)/sigma_1 < STRESS_ERROR && fabs(sigma3_2-sigma_1)/sigma_1 < STRESS_ERROR) {
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()    
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
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
			REAL sigma_values[],  
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries
 
    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];

    int min[2]={5,6};    // minimum stress acting on boundary 5 and 6
    UPDATECTL minctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }


    int  i=0;
    REAL sigma=sigma_values[i];
    REAL sigma_inc=(sigma_values[i+1]-sigma_values[i])/sigma_division;
    REAL sigma_b=sigma_values[sigma_points-1];

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
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
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
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	}

	// 8. find balanced status of odometer compression
	if (fabs(sigma3_1-sigma)/sigma < STRESS_ERROR && fabs(sigma3_2-sigma)/sigma < STRESS_ERROR ) {
	    balancedinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()    
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
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
	if (fabs(sigma3_1-sigma_b)/sigma_b < STRESS_ERROR && fabs(sigma3_2-sigma_b)/sigma_b < STRESS_ERROR) {
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()    
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"epsilon_v"<<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries
 
    // pre_3. define variables used in iterations
    REAL sigma3_1, sigma3_2;
    int    stepsnum=0;
    char   stepsstr[4];
    char   stepsfp[50];
    REAL avgNormal=0;
    REAL avgTangt=0;
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
	sigma3_1=vfabs(getNormalForce(5))/getArea(5); sigma3_2=vfabs(getNormalForce(6))/getArea(6);
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
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<0
		       <<setw(OWID)<<0<<setw(OWID)<<0
		       <<setw(OWID)<<0<<setw(OWID)<<0
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<endl;
/*
	// 8. loop break condition
	if (getAverageForce() < 1.0) {
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()    
		       <<setw(OWID)<<getAverageMoment()<<endl;
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
				   REAL  sigma,
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l56= 0;
    REAL sigma3_1, sigma3_2;
    REAL epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
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
	sigma3_1=vfabs(getNormalForce(5))/2.5e-3; sigma3_2=vfabs(getNormalForce(6))/2.5e-3;

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
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<0<<setw(OWID)<<0
		       <<setw(OWID)<<0<<setw(OWID)<<0
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()
		       <<setw(OWID)<<0<<setw(OWID)<<0<<setw(OWID)<<l56
		       <<setw(OWID)<<0
		       <<setw(OWID)<<0
		       <<setw(OWID)<<0
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<0
		       <<setw(OWID)<<0
		       <<setw(OWID)<<2.0*getActualCntctNum()/TotalNum
		       <<endl;

	}

	// 9. loop break condition: through displacement control mechanism
	if (   fabs(sigma3_1-sigma)/sigma < STRESS_ERROR && fabs(sigma3_2-sigma)/sigma < STRESS_ERROR )
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l56= 0;
    REAL sigma3_1, sigma3_2;
    REAL epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
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
	sigma3_1=vfabs(getNormalForce(5))/2.5e-3; sigma3_2=vfabs(getNormalForce(6))/2.5e-3;

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
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<0<<setw(OWID)<<0
		       <<setw(OWID)<<0<<setw(OWID)<<0
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()
		       <<setw(OWID)<<0<<setw(OWID)<<0<<setw(OWID)<<l56
		       <<setw(OWID)<<0
		       <<setw(OWID)<<0
		       <<setw(OWID)<<0
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<0
		       <<setw(OWID)<<0
		       <<setw(OWID)<<2.0*getActualCntctNum()/TotalNum
		       <<endl;

	}

/* Most time it is balanced, so use progressinf instead.
	// 8. find the balanced status and increase confining pressure
	if (   fabs(sigma1_1-sigma_a)/sigma_a < STRESS_ERROR && fabs(sigma1_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabs(sigma2_1-sigma_a)/sigma_a < STRESS_ERROR && fabs(sigma2_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabs(sigma3_1-sigma3_2)/(sigma3_1+sigma3_2)*2<=0.05) {
	    balancedinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()    
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)<<endl;
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
			REAL sigma_a,	  
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf.precision(OPREC);
    progressinf<<setw(OWID)<<"iteration"
	       <<setw(OWID)<<"possible"
	       <<setw(OWID)<<"actual"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"void"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"coordinate"
	       <<setw(OWID)<<"vibra"
	       <<setw(OWID)<<"impact"
	       <<setw(OWID)<<"wall-clock" << endl
	       <<setw(OWID)<<"number"
	       <<setw(OWID)<<"contacts"
	       <<setw(OWID)<<"contacts"
	       <<setw(OWID)<<"penetration"
	       <<setw(OWID)<<"contact_normal"
	       <<setw(OWID)<<"contact_tangt"
	       <<setw(OWID)<<"velocity"
	       <<setw(OWID)<<"omga"
	       <<setw(OWID)<<"force"
	       <<setw(OWID)<<"moment"
	       <<setw(OWID)<<"density"
	       <<setw(OWID)<<"sigma1_1"
	       <<setw(OWID)<<"sigma1_2"
	       <<setw(OWID)<<"sigma2_1"
	       <<setw(OWID)<<"sigma2_2"
	       <<setw(OWID)<<"sigma3_1"
	       <<setw(OWID)<<"sigma3_2"
	       <<setw(OWID)<<"mean_stress"
	       <<setw(OWID)<<"width"
	       <<setw(OWID)<<"length"
	       <<setw(OWID)<<"height"
	       <<setw(OWID)<<"volume"
	       <<setw(OWID)<<"epsilon_w"
	       <<setw(OWID)<<"epsilon_l"
	       <<setw(OWID)<<"epsilon_h"
	       <<setw(OWID)<<"epsilon_v"
	       <<setw(OWID)<<"ratio"
	       <<setw(OWID)<<"porosity"
	       <<setw(OWID)<<"number"
	       <<setw(OWID)<<"t_step"
	       <<setw(OWID)<<"t_step"
	       <<setw(OWID)<<"time" << endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf.precision(OPREC);
    balancedinf<<setw(OWID)<<"iteration"
	       <<setw(OWID)<<"possible"
	       <<setw(OWID)<<"actual"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"average"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"void"
	       <<setw(OWID)<<"sample"
	       <<setw(OWID)<<"coordinate"
	       <<setw(OWID)<<"vibra"
	       <<setw(OWID)<<"impact"
	       <<setw(OWID)<<"wall-clock" << endl
	       <<setw(OWID)<<"number"
	       <<setw(OWID)<<"contacts"
	       <<setw(OWID)<<"contacts"
	       <<setw(OWID)<<"penetration"
	       <<setw(OWID)<<"contact_normal"
	       <<setw(OWID)<<"contact_tangt"
	       <<setw(OWID)<<"velocity"
	       <<setw(OWID)<<"omga"
	       <<setw(OWID)<<"force"
	       <<setw(OWID)<<"moment"
	       <<setw(OWID)<<"density"
	       <<setw(OWID)<<"sigma1_1"
	       <<setw(OWID)<<"sigma1_2"
	       <<setw(OWID)<<"sigma2_1"
	       <<setw(OWID)<<"sigma2_2"
	       <<setw(OWID)<<"sigma3_1"
	       <<setw(OWID)<<"sigma3_2"
	       <<setw(OWID)<<"mean_stress"
	       <<setw(OWID)<<"width"
	       <<setw(OWID)<<"length"
	       <<setw(OWID)<<"height"
	       <<setw(OWID)<<"volume"
	       <<setw(OWID)<<"epsilon_w"
	       <<setw(OWID)<<"epsilon_l"
	       <<setw(OWID)<<"epsilon_h"
	       <<setw(OWID)<<"epsilon_v"
	       <<setw(OWID)<<"ratio"
	       <<setw(OWID)<<"porosity"
	       <<setw(OWID)<<"number"
	       <<setw(OWID)<<"t_step"
	       <<setw(OWID)<<"t_step"
	       <<setw(OWID)<<"time" << endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         mid[2]={1,3};    // boundary 1 and 3
    int         max[2]={2,4};    // boundary 2 and 4
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   midctl[2];
    UPDATECTL   maxctl[2];
    UPDATECTL   minctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
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
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
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
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
	               <<setw(OWID)<<getVibraTimeStep()
	               <<setw(OWID)<<getImpactTimeStep()
		       <<setw(OWID)<<timediffsec(timew1,timew2)
		       <<endl;
	    g_debuginf<<setw(OWID)<<g_iteration
		      <<setw(OWID)<<getTransEnergy()
		      <<setw(OWID)<<getRotatEnergy()
		      <<setw(OWID)<<bdry_penetr[1]
		      <<setw(OWID)<<bdry_penetr[2]
		      <<setw(OWID)<<bdry_penetr[3]
		      <<setw(OWID)<<bdry_penetr[4]
		      <<setw(OWID)<<bdry_penetr[5]
		      <<setw(OWID)<<bdry_penetr[6]
		      <<setw(OWID)<<bdry_cntnum[1]
		      <<setw(OWID)<<bdry_cntnum[2]
		      <<setw(OWID)<<bdry_cntnum[3]
		      <<setw(OWID)<<bdry_cntnum[4]
		      <<setw(OWID)<<bdry_cntnum[5]
		      <<setw(OWID)<<bdry_cntnum[6]
		      <<endl;
	}

	// Most time it is balanced, so use progressinf instead.
	// 8. find the balanced status and increase confining pressure
	if (   fabs(sigma1_1-sigma_a)/sigma_a < STRESS_ERROR && fabs(sigma1_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabs(sigma2_1-sigma_a)/sigma_a < STRESS_ERROR && fabs(sigma2_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabs(sigma3_1-sigma3_2)/(sigma3_1+sigma3_2)*2<=0.05) {
	    balancedinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()    
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
	               <<setw(OWID)<<getVibraTimeStep()
	               <<setw(OWID)<<getImpactTimeStep()
		       <<setw(OWID)<<timediffsec(timew1,timew2)
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
			REAL sigma_a,	  
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int         mid[2]={1,3};    // boundary 1 and 3
    int         max[2]={2,4};    // boundary 2 and 4
    int         min[2]={5,6};    // boundary 5 and 6
    UPDATECTL   midctl[2];
    UPDATECTL   maxctl[2];
    UPDATECTL   minctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
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
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
	void_ratio=Volume/getParticleVolume()-1;

	// displacement control
	if (g_iteration <= unload_step){ //loading
	    minctl[0].tran=vec(0,0,-TIMESTEP*COMPRESS_RATE);
	    minctl[1].tran=vec(0,0, TIMESTEP*COMPRESS_RATE);
	}
	else { 
	    if (reload==false) { // unloading
		if (fabs(sigma3_1-sigma_a)/sigma_a > STRESS_ERROR && 
		    fabs(sigma3_2-sigma_a)/sigma_a > STRESS_ERROR){
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
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    g_debuginf<<setw(OWID)<<g_iteration
		      <<setw(OWID)<<bdry_penetr[1]
		      <<setw(OWID)<<bdry_penetr[2]
		      <<setw(OWID)<<bdry_penetr[3]
		      <<setw(OWID)<<bdry_penetr[4]
		      <<setw(OWID)<<bdry_penetr[5]
		      <<setw(OWID)<<bdry_penetr[6]
		      <<setw(OWID)<<bdry_cntnum[1]
		      <<setw(OWID)<<bdry_cntnum[2]
		      <<setw(OWID)<<bdry_cntnum[3]
		      <<setw(OWID)<<bdry_cntnum[4]
		      <<setw(OWID)<<bdry_cntnum[5]
		      <<setw(OWID)<<bdry_cntnum[6]
		      <<endl;
	}

/*
	// 8. find the balanced status and increase confining pressure
	if (   fabs(sigma1_1-sigma_a)/sigma_a < STRESS_ERROR && fabs(sigma1_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabs(sigma2_1-sigma_a)/sigma_a < STRESS_ERROR && fabs(sigma2_2-sigma_a)/sigma_a < STRESS_ERROR
	    && fabs(sigma3_1-sigma3_2)/(sigma3_1+sigma3_2)*2<=0.05) {
	    balancedinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()    
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)<<endl;
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"epsilon_v"<<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);
    g_debuginf<<" iteration    end_bearing     side_friction   total_force"<<endl;

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    int    stepsnum=0;
    char   stepsstr[4];
    char   stepsfp[50];
    REAL avgNormal=0;
    REAL avgTangt=0;
    
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
	    REAL  f7=getShearForce( 7).getz();
	    REAL  f8=getShearForce( 8).getz();
	    REAL  f9=getShearForce( 9).getz();
	    REAL f10=getShearForce(10).getz();
	    REAL  fn=getNormalForce(12).getz();
	    g_debuginf<<setw(OWID)<<g_iteration
		      <<setw(OWID)<<fn
		      <<setw(OWID)<<(f7+f8+f9+f10)
		      <<setw(OWID)<<(fn+f7+f8+f9+f10)
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
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<t1
		       <<setw(OWID)<<t2
		       <<setw(OWID)<<(t1+t2)
		       <<setw(OWID)<<t3
		       <<setw(OWID)<<(t1+t2+t3)<<endl;
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
			      REAL dimn,
			      REAL rsize,
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"epsilon_v"<<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 

    // pre_3. define variables used in iterations
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    REAL void_ratio=0;
    
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
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<t1
		       <<setw(OWID)<<t2
		       <<setw(OWID)<<(t1+t2)
		       <<setw(OWID)<<t3
		       <<setw(OWID)<<(t1+t2+t3)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*getActualCntctNum()/TotalNum
		       <<endl;
	    g_debuginf<<setw(OWID)<<g_iteration
		      <<setw(OWID)<<getTopFreeParticlePosition().getz()
		      <<setw(OWID)<<ellipPileTipZ()
		      <<setw(OWID)<<getTopFreeParticlePosition().getz()-ellipPileTipZ()
		      <<setw(OWID)<<l13*l24*l56
		      <<setw(OWID)<<ellipPilePeneVol()
		      <<setw(OWID)<<Volume
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
				REAL dimn,
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"epsilon_v"<<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles
    createBoundary(inibdryfile);   // create boundaries.

    // pre_3. define variables used in iterations
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
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
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<t1
		       <<setw(OWID)<<t2
		       <<setw(OWID)<<(t1+t2)
		       <<setw(OWID)<<t3
		       <<setw(OWID)<<(t1+t2+t3)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    g_debuginf<<setw(OWID)<<g_iteration
		      <<setw(OWID)<<bdry_penetr[1]
		      <<setw(OWID)<<bdry_penetr[2]
		      <<setw(OWID)<<bdry_penetr[3]
		      <<setw(OWID)<<bdry_penetr[4]
		      <<setw(OWID)<<bdry_penetr[6]
		      <<setw(OWID)<<bdry_cntnum[1]
		      <<setw(OWID)<<bdry_cntnum[2]
		      <<setw(OWID)<<bdry_cntnum[3]
		      <<setw(OWID)<<bdry_cntnum[4]
		      <<setw(OWID)<<bdry_cntnum[6]
		      <<endl;
	    /*
	    g_debuginf<<setw(OWID)<<g_iteration
		      <<setw(OWID)<<getTopFreeParticlePosition().getz()
		      <<setw(OWID)<<ellipPileTipZ()
		      <<setw(OWID)<<getTopFreeParticlePosition().getz()-ellipPileTipZ()
		      <<setw(OWID)<<l13*l24*l56
		      <<setw(OWID)<<ellipPilePeneVol()
		      <<setw(OWID)<<Volume
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
				  REAL dimn,
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"epsilon_v"<<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles

    // pre_3. define variables used in iterations
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    REAL void_ratio=0;
    
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
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<t1
		       <<setw(OWID)<<t2
		       <<setw(OWID)<<(t1+t2)
		       <<setw(OWID)<<t3
		       <<setw(OWID)<<(t1+t2+t3)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*getActualCntctNum()/TotalNum
		       <<endl;
	    g_debuginf<<setw(OWID)<<g_iteration
		      <<setw(OWID)<<getTopFreeParticlePosition().getz()
		      <<setw(OWID)<<ellipPileTipZ()
		      <<setw(OWID)<<getTopFreeParticlePosition().getz()-ellipPileTipZ()
		      <<setw(OWID)<<l13*l24*l56
		      <<setw(OWID)<<ellipPilePeneVol()
		      <<setw(OWID)<<Volume
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
			       REAL dimn,
			       REAL force,
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"epsilon_v"<<endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
    balancedinf<<"pile penetrate..."<<endl
	       <<"   iteration   apply_force    pile_tip_pos     pile_force"<<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 

    // pre_3. define variables used in iterations
    REAL l13, l24, l56;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    REAL void_ratio=0;

    REAL zforce_inc=force/division;
    REAL zforce=zforce_inc;

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

	if(fabs(ellipPileForce()-zforce)/zforce < STRESS_ERROR ){
	    balancedinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<zforce
		       <<setw(OWID)<<getTopFreeParticlePosition().getz()-ellipPileTipZ()
		       <<setw(OWID)<<ellipPileForce()
		       <<endl;
	    zforce += zforce_inc;
	}

	if( g_iteration % interval == 0){
	    g_debuginf<<setw(OWID)<<g_iteration
		      <<setw(OWID)<<zforce
		      <<setw(OWID)<<getTopFreeParticlePosition().getz()-ellipPileTipZ()
		      <<setw(OWID)<<ellipPileForce()
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
	    REAL t1=getTransEnergy();
	    REAL t2=getRotatEnergy();
	    REAL t3=getPotenEnergy(-0.025);
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<t1
		       <<setw(OWID)<<t2
		       <<setw(OWID)<<(t1+t2)
		       <<setw(OWID)<<t3
		       <<setw(OWID)<<(t1+t2+t3)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*getActualCntctNum()/TotalNum
		       <<endl;
	}

	// 8. loop break condition
	if (fabs((zforce-force)/force)<0.001)
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
			    REAL sigma_a,     
			    REAL sigma_w,
			    REAL sigma_l,     
			    REAL sigma_h,   
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
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    std::ofstream balancedinf(balancedfile);
    if(!balancedinf) { cout<<"stream error!"<<endl; exit(-1);}
    balancedinf.setf(std::ios::scientific, std::ios::floatfield);
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
	       <<"height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
	       <<"        ratio          porosity         number"
	       <<endl;

    g_debuginf.open(debugfile);
    if(!g_debuginf) { cout<<"stream error!"<<endl; exit(-1);}
    g_debuginf.setf(std::ios::scientific, std::ios::floatfield);

    // pre_2. create particles and boundaries from files
    createSample(iniptclfile); // create container and particles, velocity and omga are set zero. 
    createBoundary(inibdryfile);   // create boundaries

    // pre_3. define variables used in iterations
    REAL W0 = getApt(2).gety()-getApt(4).gety();
    REAL L0 = getApt(1).getx()-getApt(3).getx();
    REAL H0 = getApt(5).getz()-getApt(6).getz();
    REAL l13, l24, l56, min_area, mid_area, max_area;
    REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
    REAL epsilon_w, epsilon_l, epsilon_h;
    REAL avgNormal=0;
    REAL avgTangt=0;
    int         stepsnum=0;
    char        stepsstr[4];
    char        stepsfp[50];
    
    int mid[2]={1,3};    // boundary 1 and 3
    int max[2]={2,4};    // boundary 2 and 4
    int min[2]={5,6};    // boundary 5 and 6
    UPDATECTL midctl[2];
    UPDATECTL maxctl[2];
    UPDATECTL minctl[2];
    REAL void_ratio=0;
    REAL bdry_penetr[7];
    int         bdry_cntnum[7];
    for (int i=0;i<7;++i){
	bdry_penetr[i]=0; bdry_cntnum[i]=0;
    }

    REAL sigma_w1=sigma_a;
    REAL sigma_l1=sigma_a;
    REAL sigma_h1=sigma_a;
    REAL sigma_w_inc=(sigma_w-sigma_a)/sigma_division;
    REAL sigma_l_inc=(sigma_l-sigma_a)/sigma_division;
    REAL sigma_h_inc=(sigma_h-sigma_a)/sigma_division;

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
	sigma1_1=vfabs(getNormalForce(2))/max_area; sigma1_2=vfabs(getNormalForce(4))/max_area;
	sigma2_1=vfabs(getNormalForce(1))/mid_area; sigma2_2=vfabs(getNormalForce(3))/mid_area;
	sigma3_1=vfabs(getNormalForce(5))/min_area; sigma3_2=vfabs(getNormalForce(6))/min_area;
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
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()   
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    g_debuginf<<setw(OWID)<<g_iteration
		      <<setw(OWID)<<bdry_penetr[1]
		      <<setw(OWID)<<bdry_penetr[2]
		      <<setw(OWID)<<bdry_penetr[3]
		      <<setw(OWID)<<bdry_penetr[4]
		      <<setw(OWID)<<bdry_penetr[5]
		      <<setw(OWID)<<bdry_penetr[6]
		      <<setw(OWID)<<bdry_cntnum[1]
		      <<setw(OWID)<<bdry_cntnum[2]
		      <<setw(OWID)<<bdry_cntnum[3]
		      <<setw(OWID)<<bdry_cntnum[4]
		      <<setw(OWID)<<bdry_cntnum[5]
		      <<setw(OWID)<<bdry_cntnum[6]
		      <<endl;
	}

	// 8. find the balanced status and increase confining pressure
	if (   fabs(sigma1_1-sigma_w1)/sigma_w1 < STRESS_ERROR && fabs(sigma1_2-sigma_w1)/sigma_w1 < STRESS_ERROR
	    && fabs(sigma2_1-sigma_l1)/sigma_l1 < STRESS_ERROR && fabs(sigma2_2-sigma_l1)/sigma_l1 < STRESS_ERROR
	    && fabs(sigma3_1-sigma_h1)/sigma_h1 < STRESS_ERROR && fabs(sigma3_2-sigma_h1)/sigma_h1 < STRESS_ERROR ) {
	    balancedinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()    
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
					+bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
					+bdry_cntnum[4]+bdry_cntnum[5]+bdry_cntnum[6])/TotalNum
		       <<endl;
	    sigma_w1 += sigma_w_inc;
	    sigma_l1 += sigma_l_inc;
	    sigma_h1 += sigma_h_inc;
	}

	// 9. loop break condition
	if (   fabs(sigma1_1-sigma_w)/sigma_w < STRESS_ERROR && fabs(sigma1_2-sigma_w)/sigma_w < STRESS_ERROR
	    && fabs(sigma2_1-sigma_l)/sigma_l < STRESS_ERROR && fabs(sigma2_2-sigma_l)/sigma_l < STRESS_ERROR
	    && fabs(sigma3_1-sigma_h)/sigma_h < STRESS_ERROR && fabs(sigma3_2-sigma_h)/sigma_h < STRESS_ERROR ) {
	    progressinf<<setw(OWID)<<g_iteration
		       <<setw(OWID)<<getPossCntctNum()
		       <<setw(OWID)<<getActualCntctNum()
		       <<setw(OWID)<<getAveragePenetration()
		       <<setw(OWID)<<avgNormal
		       <<setw(OWID)<<avgTangt
		       <<setw(OWID)<<getAverageVelocity() 
		       <<setw(OWID)<<getAverageOmga()
		       <<setw(OWID)<<getAverageForce()    
		       <<setw(OWID)<<getAverageMoment()
		       <<setw(OWID)<<getDensity()
		       <<setw(OWID)<<sigma1_1<<setw(OWID)<<sigma1_2
		       <<setw(OWID)<<sigma2_1<<setw(OWID)<<sigma2_2
		       <<setw(OWID)<<sigma3_1<<setw(OWID)<<sigma3_2
		       <<setw(OWID)<<getAverageRigidPressure()  // just the mean stress p
		       <<setw(OWID)<<l24<<setw(OWID)<<l13<<setw(OWID)<<l56
		       <<setw(OWID)<<Volume
		       <<setw(OWID)<<epsilon_w
		       <<setw(OWID)<<epsilon_l
		       <<setw(OWID)<<epsilon_h
		       <<setw(OWID)<<(epsilon_w+epsilon_l+epsilon_h)
		       <<setw(OWID)<<void_ratio
		       <<setw(OWID)<<void_ratio/(1+void_ratio)
		       <<setw(OWID)<<2.0*(getActualCntctNum()
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
void assembly::dircShear(REAL rate, REAL roterate,REAL stress,const char* iniptclfile,
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

	std::list<RGDBDRY*>::iterator rt;
	fprintf(fprslt,"bdry_1_norm_x  bdry_1_norm_y  bdry_1_norm_z  bdry_1_shar_x  bdry_1_shar_y  bdry_1_shar_z  \
bdry_2_norm_x  bdry_2_norm_y  bdry_2_norm_z  bdry_2_shar_x  bdry_2_shar_y  bdry_2_shar_z  \
bdry_3_norm_x  bdry_3_norm_y  bdry_3_norm_z  bdry_3_shar_x  bdry_3_shar_y  bdry_3_shar_z  \
bdry_4_norm_x  bdry_4_norm_y  bdry_4_norm_z  bdry_4_shar_x  bdry_4_shar_y  bdry_4_shar_z  \
bdry_5_norm_x  bdry_5_norm_y  bdry_5_norm_z  bdry_5_shar_x  bdry_5_shar_y  bdry_5_shar_z  \
bdry_6_norm_x  bdry_6_norm_y  bdry_6_norm_z  bdry_6_shar_x  bdry_6_shar_y  bdry_6_shar_z\n");

	REAL avgsigma;
	REAL l13, l24, l56, min_area, mid_area, max_area, lead;
	REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2;
	vec tmpnorm, tmpshar;
	REAL av=0;
	REAL ao=0;
	REAL af=0;
	REAL am=0;
	REAL avgNormal=0;
	REAL avgTangt=0;

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
		progressinf<<setw(OWID)<<g_iteration;

		findParticleOnBoundary();
		findContact();

		internalForce(avgNormal, avgTangt);
		rigidBoundaryForce();
		//track(fp,5);
		progressinf<<setw(OWID)<<getAverageVelocity()
		         <<setw(OWID)<<getAverageOmga()
		         <<setw(OWID)<<getAverageForce()
			 <<setw(OWID)<<getAverageMoment()<<endl;
		contactUpdate();
		updateParticle();

		l56=getApt(5).getz()-getApt(6).getz();
		l24=getApt(2).gety()-getApt(4).gety();
		l13=getApt(1).getx()-getApt(3).getx();
		min_area=l13*l24;
		mid_area=l56*l24;
		lead=fabs(normalize(getDirc(2))%vec(0,1,0));
		max_area=l56*l13;
		setArea(5,min_area);
		setArea(6,min_area);
		setArea(1,mid_area);
		setArea(3,mid_area);
		setArea(2,max_area);
		setArea(4,max_area);
		avgsigma=getAverageRigidPressure();
		printf("avgsigma=%15.3lf\n",avgsigma);
		sigma1_1=fabs(getNormalForce(2))/max_area;
		sigma1_2=fabs(getNormalForce(4))/max_area;
		sigma2_1=fabs(getNormalForce(1))/mid_area;
		sigma2_2=fabs(getNormalForce(3))/mid_area;
		sigma3_1=fabs(getNormalForce(5))/min_area;
		sigma3_2=fabs(getNormalForce(6))/min_area;
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
	fclose(fp);
	fclose(fprslt);
	printParticle(resultfile);
}
*/

/* 
void assembly::soft_tric(REAL _sigma3,REAL _b,const char* iniptclfile,
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

	std::list<RGDBDRY*>::iterator rt;

	int max[2]={1,2};//maximum stress acting on boundary 5 and 6
	UPDATECTL maxctl[2];
	REAL loading_rate=0.01;

	REAL avgsigma;
	REAL af, av, am, ao, adr, pre_af;
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
		progressinf<<setw(OWID)<<g_iteration;

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
		progressinf<<setw(OWID)<<(av=getAverageVelocity())
		         <<setw(OWID)<<(ao=getAverageOmga())
		         <<setw(OWID)<<(af=getAverageForce())
		         <<setw(OWID)<<(am=getAverageMoment())
			 <<setw(OWID)<<(adr=avgDgrFric());
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

	std::list<RGDBDRY*>::iterator rt;

//	int mid[2]={2,4};//intermediate stress acting on boundary 2 and 4
//	UPDATECTL midctl[2];
	int max[2]={5,6};//maximum stress acting on boundary 5 and 6
	UPDATECTL maxctl[2];
//	int min[2]={1,3};//minimum stress acting on boundary 1 and 3
//	UPDATECTL minctl[2];
	REAL loading_rate=0.01;

	REAL avgsigma;
	REAL af, av, am, ao, adr, pre_af;
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
		progressinf<<setw(OWID)<<g_iteration;

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
		progressinf<<setw(OWID)<<(av=getAverageVelocity())
		         <<setw(OWID)<<(ao=getAverageOmga())
		         <<setw(OWID)<<(af=getAverageForce())
		         <<setw(OWID)<<(am=getAverageMoment())
			 <<setw(OWID)<<(adr=avgDgrFric());

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
	fclose(fp);
	fclose(fprslt);
	printParticle(resultfile);
}
*/

/* 
void assembly::simpleShear(REAL _sigma3,REAL _b,
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

	std::list<RGDBDRY*>::iterator rt;

	int mid[2]={2,4};//intermediate stress acting on boundary 2 and 4
	UPDATECTL midctl[2];
	int max[2]={5,6};//maximum stress acting on boundary 5 and 6
	UPDATECTL maxctl[2];
	int min[2]={1,3};//minimum stress acting on boundary 1 and 3
	UPDATECTL minctl[2];
//	REAL loading_rate=0.01;
	REAL increment=0.0001;
	REAL angular_velocity=0.1;
	vec increment_velocity_x(increment,0,0);
	vec increment_velocity_y(0,increment,0);
	vec increment_velocity_z(0,0,increment);
	vec xbdry_velocity_0,xbdry_velocity_1;
	vec ybdry_velocity_0,ybdry_velocity_1;
	vec zbdry_velocity_0,zbdry_velocity_1;

	REAL avgsigma;
	REAL af, av, am, ao, adr, pre_af;
	REAL sigma1_1, sigma1_2, sigma2_1, sigma2_2, sigma3_1, sigma3_2, sigma2;
	REAL ita1_1, ita1_2, ita2_1, ita2_2, ita3_1, ita3_2;
	REAL ar;
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
		progressinf<<setw(OWID)<<g_iteration;

		findParticleOnBoundary();
		findContact();

		internalForce();
		rigidBoundaryForce();
		flexiBoundaryForce();
		//track(fp,5);
		progressinf<<setw(OWID)<<(av=getAverageVelocity())
		         <<setw(OWID)<<(ao=getAverageOmga())
		         <<setw(OWID)<<(af=getAverageForce())
		         <<setw(OWID)<<(am=getAverageMoment())
			 <<setw(OWID)<<(adr=avgDgrFric());
		contactUpdate();
		
		avgsigma=getAverageRigidPressure();
		sigma1_1=fabs(getNormalForce(5))/getArea(5);
		ita1_1=fabs(getShearForce(5))/getArea(5);
		sigma1_2=fabs(getNormalForce(6))/getArea(6);
		ita1_2=fabs(getShearForce(6))/getArea(6);
		sigma2_1=fabs(getNormalForce(2))/getArea(2);
		ita2_1=fabs(getShearForce(2))/getArea(2);
		sigma2_2=fabs(getNormalForce(4))/getArea(4);
		ita2_2=fabs(getShearForce(4))/getArea(4);
		sigma3_1=fabs(getNormalForce(1))/getArea(1);
		ita3_1=fabs(getShearForce(1))/getArea(1);
		sigma3_2=fabs(getNormalForce(3))/getArea(3);
		ita3_2=fabs(getShearForce(3))/getArea(1);
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
			
			if(af<0.02&&fabs(sigma3_1-_sigma3)<0.02*_sigma3
				  &&fabs(sigma3_2-_sigma3)<0.02*_sigma3
				  &&fabs(sigma2_1-sigma2)<0.02*_sigma3
				  &&fabs(sigma2_2-sigma2)<0.02*_sigma3
				  &&fabs(sigma1_1-_sigma3)<0.02*_sigma3
				  &&fabs(sigma1_2-_sigma3)<0.02*_sigma3
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
	fclose(fp);
	fclose(fprslt);
	printParticle(resultfile);
}
*/

/* 
void assembly::earthPressure(REAL pressure,bool IsPassive, 
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

	std::list<RGDBDRY*>::iterator rt;

	int wall[1]={1};
	UPDATECTL wallctl[1];
//	REAL loading_rate=0.001;

	REAL avgsigma;
	REAL af, av, am, ao, adr, pre_af;
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
		progressinf<<setw(OWID)<<g_iteration;

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

		progressinf<<setw(OWID)<<(av=getAverageVelocity())
		         <<setw(OWID)<<(ao=getAverageOmga())
		         <<setw(OWID)<<(af=getAverageForce())
		         <<setw(OWID)<<(am=getAverageMoment())
			 <<setw(OWID)<<(adr=avgDgrFric());
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
	fclose(fp);
	fclose(fprslt);
	printParticle(resultfile);
}
*/
