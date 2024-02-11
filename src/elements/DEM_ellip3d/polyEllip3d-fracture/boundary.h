//     1.This is a template class, for which we have to include the implementation in the header file.
//       As we cannot put using statement in a header file, we have to use std::xxx wherever we need to
//       refer anything from standard namespace.
//
//     2.A base class needs a virtual destructor, otherwise it may cause undetermined errors.
//
//     3.When inheritating a template class, it is important to know the idea "Name lookup, templates,
//       and accessing members of base classes". 
//       Reference: http://gcc.gnu.org/onlinedocs/gcc-4.0.2/gcc/Name-lookup.html#Name-lookup
//
//     all cylinder boudaries are those cylinders with vertical motherlines right now
//     but it is very convient to extend to include cylinders with non-vertical motherlines
//     all plane flexible boundaries are considered constructed by two segments of straight line

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "realtypes.h"
#include "vec.h"
#include "parameter.h"
#include "cylinder.h"
#include "boundarytgt.h"
#include <map>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

using std::cout;
using std::setw;
using std::endl;

namespace dem {

// BdryCoef is used for rigid boundary conditions
typedef struct bdryfunc{
	int order;  // x1-linear; 2-quadratic
	vec dirc;   // normal vector if plane, mother line vector if cylinder,it points out of the particles		
	vec apt;    // a point on the plane or a point on the axis of the cylinder
	REAL rad;   //zero if plane
	int side;   //zero if plane; side=1, particles are outside the cylinder; =-1, inside the cylinder
        void disp() const{
	  cout << "order: " << order << endl;
	  cout << "dirc: " << dirc.getx() << " " << dirc.gety() << " " << dirc.getz() << endl;
	  cout << "apt: " << apt.getx() << " " << apt.gety() << " " << apt.getz() << endl;
	  cout << "radius: " << rad << " side: " << side << endl;
	}
	void disp(std::ofstream &ofs) const{
	    ofs<<setw(OWID)<<order
	       <<setw(OWID)<<dirc.getx()
	       <<setw(OWID)<<dirc.gety()
	       <<setw(OWID)<<dirc.getz()
	       <<setw(OWID)<<apt.getx()
	       <<setw(OWID)<<apt.gety()
	       <<setw(OWID)<<apt.getz()
	       <<setw(OWID)<<rad
	       <<setw(OWID)<<side<<std::endl;
	}
} BdryCoef;

// LINE and CIRC are structs used in flexible boundary classes
typedef struct updatectl{
	vec tran;    // tranlation second
	vec rote;    // rotate first
	vec fixpt;   // before any update is made
	REAL expnd;// expand last
	updatectl(){tran=0;rote=0;fixpt=0;expnd=1;}
	void disp() const{
	  cout << "tran: " << tran.getx() << " " << tran.gety() << " " << tran.getz() << endl;
	  cout << "rote: " << rote.getx() << " " << rote.gety() << " " << rote.getz() << endl;
	  cout << "fixpt: " << fixpt.getx() << " " << fixpt.gety() << " " << fixpt.getz() << endl;
	  cout << "expand:" << expnd << endl;
	};
}UPDATECTL;

typedef struct line{
	vec pt1; // begining point
	vec pt2; // ending point
	void disp() const{
	  cout  << "pt1: " << pt1.getx() << " " << pt1.gety() << " " << pt1.getz() << endl;
	  cout  << "pt2: " << pt2.getx() << " " << pt2.gety() << " " << pt2.getz() << endl;
	}
	void update(UPDATECTL& ctl){
	   vec tmp1=ctl.tran+ctl.fixpt+rotateVec((pt1-ctl.fixpt),ctl.rote);   
	   vec tmp2=ctl.tran+ctl.fixpt+rotateVec((pt2-ctl.fixpt),ctl.rote); 
	   pt1=(tmp1+tmp2)/2+ctl.expnd*(tmp1-tmp2)/2;
	   pt2=(tmp1+tmp2)/2+ctl.expnd*(tmp2-tmp1)/2;
	} 
	//update pt1 and pt2 by both translating about the center and rotating about fixpt and expand with expnd;
}LINE;

typedef struct circ{
	vec center; // center of the circle
	vec norm;   // normal dirction of the circular plane, pointing out of the assembly
	int turn;   // turn=1, right-hand rule from pt1 to pt2 same direction as norm
		    // turn=-1, right-hand rule from pt1 to pt2 opposite direction as norm
	REAL radius; //the radius of the circle
	vec pt1;    // the begining point of a part circle
	vec pt2;    // the end point of a part circle;it pt1==pt2, a closure circle
	void disp() const {
	  cout << "center: " << center.getx() << " " << center.gety() << " " <<center.getz() << endl;
	  cout << "norm: " << norm.getx() << " " << norm.gety() << " " << norm.getz() << endl;
	  cout << "pt1: " << pt1.getx() << " " << pt1.gety() << " " << pt1.getz() << endl;
	  cout << "pt2: " << pt2.getx() << " " << pt2.gety() << " " << pt2.getz() << endl;
	  cout << "turn: " << turn << " radius: " << radius << endl;
	}
	void update(UPDATECTL& ctl){
		center=ctl.tran+ctl.fixpt+rotateVec(center-ctl.fixpt,ctl.rote);
		norm=rotateVec(norm,ctl.rote);
		pt1=ctl.tran+ctl.fixpt+rotateVec(pt1-ctl.fixpt,ctl.rote);
		pt2=ctl.tran+ctl.fixpt+rotateVec(pt2-ctl.fixpt,ctl.rote);
	}
}CIRC;

template<class T> class rgd_bdry{
public:
	int  bdry_id;    // the first record defines the bdry itself, the other 
	std::vector<BdryCoef> CoefOfLimits; // limitnum records define the other lines on the bdry 
	REAL avg_normal; // that give limits to the first boundary.
	REAL avg_penetr; // average penetration by particles onto this boundary
	int  cntnum;     // contact numbers by particles onto this boundary
	REAL area;       // the bounary's area
	int  limitnum;   // how many lines the boundary have
public:
	rgd_bdry(std::ifstream &ifs);
	int getBdryID() {return bdry_id;}
	virtual ~rgd_bdry() {} // base class needs a virtual destructor.
	virtual void disp() const{
		std::vector<BdryCoef>::const_iterator it;
		for(it=CoefOfLimits.begin();it!=CoefOfLimits.end();++it)
			(*it).disp();
		cout << "area: " << area << " limitnum: " << limitnum << endl;
	}
	virtual void disp(std::ofstream &ofs) const{
		std::vector<BdryCoef>::const_iterator it;
		ofs <<std::endl
		    <<setw(OWID)<<(*CoefOfLimits.begin()).order<<std::endl;
		ofs <<setw(OWID)<<bdry_id
		    <<setw(OWID)<<limitnum
		    <<setw(OWID)<<area <<std::endl;
		for(it=CoefOfLimits.begin();it!=CoefOfLimits.end();++it)
			(*it).disp(ofs);
	}
	virtual void findParticleOnBoundary(std::vector<T*>& ptcls){};
	virtual void rigidBF(std::map<int,std::vector<boundarytgt> >& BdryTgtMap)
	    {std::cout<<"parent"<<std::endl;} // calculate for each boundary particles the rigid boundary force
	virtual vec getNormalForce() const{return 0;}
	virtual REAL getAvgNormal() const{return 0;}
	virtual REAL getAvgPenetr() const{return 0;}
	virtual int         getCntnum() const{return 0;}
	virtual vec getShearForce() const{return 0;}
	virtual vec getApt() const{return 0;}
	virtual vec getDirc() const{return 0;}
	virtual void setArea(REAL a){area=a;}
	virtual REAL getArea(){return area;}
	virtual void update(UPDATECTL& ctl); //the boundary is translating with tran and rotating with rote around axis
};

template <class T>
rgd_bdry<T>::rgd_bdry(std::ifstream &ifs){
	BdryCoef tmp;
	REAL x,y,z;
	CoefOfLimits.clear();
	ifs >> bdry_id >> limitnum >> area;
	for (int k=0;k<limitnum;k++){
	    ifs >> tmp.order >> x >> y >> z;
	    tmp.dirc=vec(x,y,z);
	    ifs >> x >> y >> z;
	    tmp.apt=vec(x,y,z);
	    ifs >> tmp.rad >> tmp.side;
	    CoefOfLimits.push_back(tmp);
	}
};

template<class T>
void rgd_bdry<T>::update(UPDATECTL& ctl){
	BdryCoef tmp;
	std::vector<BdryCoef>::iterator it;
	vec nv, napt;
	for (it=CoefOfLimits.begin();it!=CoefOfLimits.end();++it){
		tmp=*it;
		nv=rotateVec(tmp.dirc,ctl.rote);
		napt=ctl.tran+ctl.fixpt+rotateVec(tmp.apt-ctl.fixpt,ctl.rote);
		(*it).dirc=nv;
		(*it).apt=napt;
	}
};

template<class T> class flb_bdry{
public:
	int bdry_id;
	virtual void disp() const{};
	virtual void findParticleOnBoundary(std::vector<T*>& ptcls){};
	virtual void update(UPDATECTL ctl[], unsigned int len){};
	virtual void findParticleOnLine(){}; // create possible particles per line
	virtual void createFlbNet(){};
	virtual void flxbBF(){};
	virtual vec triangleDstr(REAL pressure,vec norm, vec p[], T* e[]); //norm is the direction of pressure
	virtual ~flb_bdry() {};     // base class needs a virtual destructor.
};

template<class T>
vec flb_bdry<T>::triangleDstr(REAL pressure, vec norm, vec p[], T* e[]){
        //norm indicates the pressure dirction
	vec cent=(p[0]+p[1]+p[2])/3;
	REAL l1=vfabs(p[1]-p[0]);
	REAL l2=vfabs(p[2]-p[1]);
	REAL l3=vfabs(p[0]-p[2]);
	REAL hp=(l1+l2+l3)/2;
	REAL area=sqrt(hp)*sqrt(hp-l1)*sqrt(hp-l2)*sqrt(hp-l3);
	vec nm=normalize((p[0]-p[1])*(p[2]-p[1]));
	if(nm%norm<0)
		nm*=-1;
	vec nf=pressure*area*nm;
	int nonull=0;
 	int i;
	for (i=0;i<3;i++){
		if (e[i]!=NULL)
			nonull++;
	}
	for (i=0;i<3;i++){
	    if (e[i]!=NULL){
		    //e[i]->addForce(nf/nonull);
		    e[i]->flb_force+=nf/nonull;
		    //vec moment=(cent-e[i]->getCurrPosition())*nf/nonull;
		    //e[i]->addMoment(e[i]->localVec(moment));
		    //e[i]->flb_moment+=moment;
	    }
	}
	return nf;
};

template<class T> class plnrgd_bdry:public rgd_bdry<T>{
public:
	vec normal;  // normal force acting on the boundary by all contacting particles 
	vec tangt;   // tangential force acting on the boundary
	vec moment;  // moment on the boundary
	std::vector<T*> PBVec; // possible boundary particles of this specific boundary
public:
	plnrgd_bdry(std::ifstream &ifs):rgd_bdry<T>(ifs){
	    normal=0;
	    tangt=0;
	    moment=0;
	    this->avg_normal=0;
	    this->avg_penetr=0;
	    this->cntnum=0;
	};
	int getBdryID() {return this->bdry_id;}
	void disp() const;
	REAL distToBdry(vec posi) const;
	void findParticleOnBoundary(std::vector<T*>& ptcls);
	vec getApt() const;
	vec getDirc() const;
	plnrgd_bdry<T>* getBdry(int bdryid) const{
		return this;
	}
	void rigidBF(std::map<int,std::vector<boundarytgt> >& BdryTgtMap);
	vec getNormalForce() const{return normal;}
	REAL getAvgNormal() const{return this->avg_normal;}
	REAL getAvgPenetr() const{return this->avg_penetr;}
	int         getCntnum() const{return this->cntnum;}

	vec getShearForce() const{return tangt;}
};

template<class T>
vec plnrgd_bdry<T>::getApt() const{
	return (*this->CoefOfLimits.begin()).apt;
};

template<class T>
vec plnrgd_bdry<T>::getDirc() const{
	return (*this->CoefOfLimits.begin()).dirc;
};

template<class T>
void plnrgd_bdry<T>::disp() const{
	rgd_bdry<T>::disp();
	cout << "normal: " << normal.getx() << " " << normal.gety() << " " <<normal.getz() << endl;
	typename std::vector<T*>::const_iterator it;
	int i=0;
	for(it=PBVec.begin();it!=PBVec.end();++it){
		if(i++<10)
			cout << (*it)->getID();
		else{
			i=0;
			cout << (*it)->getID() << endl;
		}
	}
};

template<class T>
REAL plnrgd_bdry<T>::distToBdry(vec posi) const{
	vec dv=(*this->CoefOfLimits.begin()).dirc;	// CoefOfLimits.begin() is the boundary itself, 
	vec pt=(*this->CoefOfLimits.begin()).apt;	// others are crossing boundaries
	vec ndv=normalize(dv);
	return (posi-pt)%ndv;
};

template<class T>
void plnrgd_bdry<T>::findParticleOnBoundary(std::vector<T*>& ptcls){	// August 22, 2013
    typename std::vector<T*>::iterator it;
    std::vector<BdryCoef>::iterator bt;
    bool next;
    PBVec.clear();
    REAL dist, r;
    vec posi, ndirc;

    // for boundary predetection as in notes page 106
    // at present it can just deal with the boundaries that are normal to the axes
    vec dv=(*this->CoefOfLimits.begin()).dirc;
    vec pt=(*this->CoefOfLimits.begin()).apt;
    vec ndv=normalize(dv);

    for (it=ptcls.begin();it!=ptcls.end();++it){
	if ( (*it)->getType() == 0 ) { // only process free particles, excluding type 5
	    posi=(*it)->getCurrPosition();
	    dist=distToBdry(posi);
	    if(dist>=0 || fabs(dist) > (*it)->getMaxRadius()) // outside to CoefOfLimits[0] or inside too much
		continue;


	    // predetection, using the same method as in contact predetection, notes page 106
	    // April 14, 2014
	    REAL aplus  = (*it)->getAplus();
	    REAL aminus = (*it)->getAminus();
	    REAL bplus  = (*it)->getBplus();
	    REAL bminus = (*it)->getBminus();
	    REAL cplus  = (*it)->getCplus();
	    REAL cminus = (*it)->getCminus(); 

    	    vec local_p1 = vec(aplus, bplus, cplus);
    	    vec local_p2 = vec(-aminus, bplus, cplus);
    	    vec local_p3 = vec(-aminus, -bminus, cplus);
    	    vec local_p4 = vec(aplus, -bminus, cplus);
    	    vec local_p5 = vec(aplus, bplus, -cminus);
    	    vec local_p6 = vec(-aminus, bplus, -cminus);
    	    vec local_p7 = vec(-aminus, -bminus, -cminus);
    	    vec local_p8 = vec(aplus, -bminus, -cminus);

	    vec global_p1 = posi + (*it)->globalVec(local_p1);
	    vec global_p2 = posi + (*it)->globalVec(local_p2);
	    vec global_p3 = posi + (*it)->globalVec(local_p3);
	    vec global_p4 = posi + (*it)->globalVec(local_p4);
	    vec global_p5 = posi + (*it)->globalVec(local_p5);
	    vec global_p6 = posi + (*it)->globalVec(local_p6);
	    vec global_p7 = posi + (*it)->globalVec(local_p7);
	    vec global_p8 = posi + (*it)->globalVec(local_p8);

	    // judge which boundary this is
	    // now it can deal with boundaries that is not normal to axles
	    // check if the eight corner points are all inside the boundary or not
	    if(distToBdry(global_p1)<=0 && distToBdry(global_p2)<=0 && distToBdry(global_p3)<=0 && distToBdry(global_p4)<=0
	    && distToBdry(global_p5)<=0 && distToBdry(global_p6)<=0 && distToBdry(global_p7)<=0 && distToBdry(global_p8)<=0 )
		continue;

/*	    if(ndv.getx()>0){	// x+
		if(global_p1.getx() < pt.getx() && global_p2.getx() < pt.getx() && global_p3.getx() < pt.getx() && global_p4.getx() < pt.getx()
		&& global_p5.getx() < pt.getx() && global_p6.getx() < pt.getx() && global_p7.getx() < pt.getx() && global_p8.getx() < pt.getx() ){
		    // now is judge if the particle is inside too much, so we need to use <
		    continue;
		}
	    }
	    if(ndv.getx()<0){	// x-
		if(global_p1.getx() > pt.getx() && global_p2.getx() > pt.getx() && global_p3.getx() > pt.getx() && global_p4.getx() > pt.getx()
		&& global_p5.getx() > pt.getx() && global_p6.getx() > pt.getx() && global_p7.getx() > pt.getx() && global_p8.getx() > pt.getx() ){
		    // now is judge if the particle is inside too much, so we need to use >
		    continue;
		}
	    }
	    if(ndv.gety()>0){	// y+
		if(global_p1.gety() < pt.gety() && global_p2.gety() < pt.gety() && global_p3.gety() < pt.gety() && global_p4.gety() < pt.gety()
		&& global_p5.gety() < pt.gety() && global_p6.gety() < pt.gety() && global_p7.gety() < pt.gety() && global_p8.gety() < pt.gety() ){
		    // now is judge if the particle is inside too much, so we need to use <
		    continue;
		}
	    }
	    if(ndv.gety()<0){	// y-
		if(global_p1.gety() > pt.gety() && global_p2.gety() > pt.gety() && global_p3.gety() > pt.gety() && global_p4.gety() > pt.gety()
		&& global_p5.gety() > pt.gety() && global_p6.gety() > pt.gety() && global_p7.gety() > pt.gety() && global_p8.gety() > pt.gety() ){
		    // now is judge if the particle is inside too much, so we need to use >
		    continue;
		}
	    }
	    if(ndv.getz()>0){	// z+
		if(global_p1.getz() < pt.getz() && global_p2.getz() < pt.getz() && global_p3.getz() < pt.getz() && global_p4.getz() < pt.getz()
		&& global_p5.getz() < pt.getz() && global_p6.getz() < pt.getz() && global_p7.getz() < pt.getz() && global_p8.getz() < pt.getz() ){
		    // now is judge if the particle is inside too much, so we need to use <
		    continue;
		}
	    }
	    if(ndv.getz()<0){	// z-
		if(global_p1.getz() > pt.getz() && global_p2.getz() > pt.getz() && global_p3.getz() > pt.getz() && global_p4.getz() > pt.getz()
		&& global_p5.getz() > pt.getz() && global_p6.getz() > pt.getz() && global_p7.getz() > pt.getz() && global_p8.getz() > pt.getz() ){
		    // now is judge if the particle is inside too much, so we need to use >
		    continue;
		}
	    }
*/

	    /*
	    g_debuginf<<"boundary.h: g_iter="<<g_iteration
		      <<" bdryId="<<getBdryID()
		      <<" ptclId="<<(*it)->getID()<<std::endl;
	    */
	    next=true;
	    for (bt=++this->CoefOfLimits.begin();bt!=this->CoefOfLimits.end();++bt){ // CoefOfLimits[1,2,...]
		ndirc=normalize((*bt).dirc);
		r=vfabs((posi-(*bt).apt)-(posi-(*bt).apt)%ndirc*ndirc);
		if( ( (*bt).order==1 && (posi-(*bt).apt)%(*bt).dirc >= 0 ) ||
		    ( (*bt).order==2 && (r-(*bt).rad)*(*bt).side<0 ) ){
		    next=false; // the particle is out of boundary, process next particle
		    break;
		}
	    }
	    if(next)
		PBVec.push_back(*it);
	}
    }

    /*
    if (getBdryID() > 6) { // cavity boundaries
      g_debuginf<<"boundary.h: g_iter="<<g_iteration
		<<" bdryId="<<getBdryID()
		<<" ptcl#="<<PBVec.size();
      for(it = PBVec.begin(); it != PBVec.end(); ++it)
	g_debuginf << " " << (*it)->getID();
      g_debuginf << std::endl;
    }
    */
    
 };

/*
template<class T>
void plnrgd_bdry<T>::findParticleOnBoundary(std::vector<T*>& ptcls){
	typename std::vector<T*>::iterator it;
	std::vector<BdryCoef>::iterator bt;
	bool next;
	PBVec.clear();
 	for (it=ptcls.begin();it!=ptcls.end();++it){
		vec posi=(*it)->getCurrPosition();
		REAL dist=distToBdry(posi);
		//if (fabs(dist)>ROOM*(*it)->getA()&&dist<0)
		if(fabs(dist)>ROOM*(*it)->getA()||dist>0)
			continue;
		next=false;
		for (bt=++CoefOfLimits.begin();bt!=CoefOfLimits.end();++bt){
			vec ndirc=normalize((*bt).dirc);
			REAL r=abs((posi-(*bt).apt)-(posi-(*bt).apt)%ndirc*ndirc);
			if((*bt).order==1&&(posi-(*bt).apt)%(*bt).dirc>0||
				(*bt).order==2&&(r-(*bt).rad)*(*bt).side<0){
				next=true;//the particle is outof boundary, process next particle
				break;
			}
		}
		if(!next)
			PBVec.push_back(*it);
	}
};
*/


template<class T>
void plnrgd_bdry<T>::rigidBF(std::map<int,std::vector<boundarytgt> >& BdryTgtMap){
    typename std::vector<T*>::iterator it;
    this->avg_normal=0;
    this->avg_penetr=0;
    this->cntnum=0;
    normal=0;
    tangt=0;
    moment=0;

    // for each plane boundary, define a temparory variable vtmp to use,
    // better than define a member variable which needs to be cleared.
    // and vtmp is initialized as empty in each iteration.
    std::vector<boundarytgt> vtmp;

    // for each possible boundary particle
    REAL penetr=0;
    int count=0;
    for (it=PBVec.begin();it!=PBVec.end();++it){
	penetr=0;
	(*it)->planeRBForce(this,BdryTgtMap,vtmp,penetr);
	this->avg_penetr += penetr;
	count++;
    }
    if (count>0) this->avg_penetr /= count;
    this->cntnum=count;

    // checkout tangential forces and displacements after each particle is processed
    BdryTgtMap[this->bdry_id]=vtmp;
};

template<class T> class cylrgd_bdry:public rgd_bdry<T>{
public:
	vec normal; 
	std::vector<T*> PBVec;
public:
	cylrgd_bdry(std::ifstream &ifs):rgd_bdry<T>(ifs){normal=0;}
	void disp() const;
	REAL distToBdry(vec posi) const;
	void findParticleOnBoundary(std::vector<T*>& ptcls);
	void rigidBF();
	vec getNormalForce() const{return normal;};
};

template<class T>
void cylrgd_bdry<T>::disp() const{
	rgd_bdry<T>::disp();
	cout << "normal: " << normal.getx() << " " << normal.gety() << " " <<normal.getz() << endl;
	typename std::vector<T*>::const_iterator it;
	int i=0;
	for(it=PBVec.begin();it!=PBVec.end();++it){
		if(i++<10)
		  cout << (*it)->getID();
		else{
		  i=0;
		  cout << (*it)->getID() << endl;
		}
	}
};

template<class T>
REAL cylrgd_bdry<T>::distToBdry(vec posi) const{
	vec ndc=(*this->CoefOfLimits.begin()).dirc;
	vec napt=(*this->CoefOfLimits.begin()).apt;
	REAL r=(*this->CoefOfLimits.begin()).rad;
	vec norm=normalize(ndc);
	return fabs(r-vfabs((posi-napt)-(posi-napt)%norm*norm));
};

template<class T>
void cylrgd_bdry<T>::findParticleOnBoundary(std::vector<T*> &ptcls){	// August 22, 2013
	typename std::vector<T*>::iterator it;
	std::vector<BdryCoef>::iterator bt;
	bool next;
	PBVec.clear();
	REAL dist,r;
	vec posi, ndirc;
 	for (it=ptcls.begin();it!=ptcls.end();++it){
		posi=(*it)->getCurrPosition();
		dist=distToBdry(posi);
		if (dist>(*it)->getMaxRadius())
			continue;
		next=false;
		for (bt=++this->CoefOfLimits.begin();bt!=this->CoefOfLimits.end();++bt){
			ndirc=normalize((*bt).dirc);
			r=vfabs((posi-(*bt).apt)-(posi-(*bt).apt)%ndirc*ndirc);
			if( ( (*bt).order==1&&(posi-(*bt).apt)%(*bt).dirc>(*it)->getMaxRadius() )||
				( (*bt).order==2&&(r-(*bt).rad)*(*bt).side<0 ) ){
				next=true;//the particle is outof boundary, process next particle
				break;
			}
		}
		if(!next)
			PBVec.push_back(*it);
	}
};

template<class T>
void cylrgd_bdry<T>::rigidBF(){
	// I am temporially saitisfied with the cylinder with vertical mother line
	cylinder cyl;
	typename std::vector<T*>::iterator it;
	BdryCoef tmp;
	tmp=*this->CoefOfLimits.begin();
	cyl.setRadius(tmp.rad);
	cyl.setCenter(tmp.apt);
	normal=0;
	for (it=PBVec.begin();it!=PBVec.end();++it){
		normal-=(*it)->cylinderRBForce(this->bdry_id,cyl,tmp.side);
	}
};

template<class T> class plnflb_bdry:public flb_bdry<T>{
public:
	vec sumpressure; // sum of total water pressure on particles
	REAL confining;// confining pressure by surrounding liquid
	std::vector<LINE> framelist;//store rigid lines it
	vec norm;        // normal direction pointing outward the assembly
	int framenum;    // how many rigid lines there are, can only be 2
	std::vector<T*> PBVec;
	vec FlxbNet[100][50];
	T* RelatedP[100][50];
	std::vector<T*> PCFB[100][50];
	int np, nz;
public:
	plnflb_bdry(std::ifstream &ifs);
	virtual ~plnflb_bdry() {}; // base class needs a virtual destructor.
	void disp() const;
	void findParticleOnBoundary(std::vector<T*>& ptcls);
	void findParticleOnLine(); // create possible particles per line
	void createFlbNet();
	void flxbBF();    // FlxbNet[nz][np], RelatedP[nz][np]; if side=1, particle is in the side of >0, side=-1, <0
	void update(UPDATECTL ctl[], unsigned int len);
	void delNull();
};

template<class T>
void plnflb_bdry<T>::disp() const{
	int iz,ip,i;
	cout << "sumpressure: " << sumpressure.getx() << " " << sumpressure.gety() << " " << sumpressure.getz() << endl;
        cout << "norm: " << norm.getx() << " " << norm.gety() << " " << norm.getz() << endl;
        cout << "np:" << np << " nz: " << nz << " framenum: " << framenum << endl;
        cout << "confining pressure: " << confining << endl;
	std::vector<LINE>::const_iterator it;
	for(it=framelist.begin();it!=framelist.end();++it)
		(*it).disp();
	typename std::vector<T*>::const_iterator jt;
	i=0;
	for(jt=PBVec.begin();jt!=PBVec.end();++jt){
		if(i++<10)
			cout << (*jt)->getID();
		else{
			i=0;
			cout << (*jt)->getID() << endl;
		}
	}
	cout << "point of the net" << endl;
	for(iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
		  cout << FlxbNet[iz][ip].getx() << " " <<FlxbNet[iz][ip].gety() << " " << FlxbNet[iz][ip].getz();
		}
		cout << endl;
	}
	i=0;
	cout << "particles of the net" << endl;
	for(iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			if(i++<10) {
				if(RelatedP[iz][ip]!=NULL) cout << RelatedP[iz][ip]->getID();
			}
			else{
				i=0;
				if(RelatedP[iz][ip]!=NULL) cout << RelatedP[iz][ip]->getID() << endl;
			}
		}
	}
	cout << "possible particles around a line" << endl;
	for(iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			for (jt=PCFB[iz][ip].begin();jt!=PCFB[iz][ip].end();++jt)
				if((*jt)!=NULL) cout << (*jt)->getID();
			cout << endl;
		}
	}
};

template<class T>
plnflb_bdry<T>::plnflb_bdry(std::ifstream &ifs){
	vec tmp;
	REAL x,y,z;
	int i,j;

	framelist.clear();
	ifs >> framenum >> nz >> np >> confining;
	for (i=0;i<framenum;i++){
	    ifs >> x >> y >> z;
	    tmp=vec(x,y,z);
	    LINE tmpln;
	    tmpln.pt1=tmp;
	    ifs >> x >> y >> z;
	    tmp=vec(x,y,z);
	    tmpln.pt2=tmp;
	    framelist.push_back(tmpln);
	}
	ifs >> x >> y >> z;
	norm=vec(x,y,z);
	for (i=0;i<nz;i++){
	    for(j=0;j<np;j++){
		RelatedP[i][j]=NULL;
		PCFB[i][j].clear();
	    }
	}
};

template<class T>
void plnflb_bdry<T>::findParticleOnBoundary(std::vector<T*>& ptcls){
	/*
		1----2
		|    |
		|    |
		4----3
	*/
	typename std::vector<T*>::iterator it;
	vec pt1=(*framelist.begin()).pt1;
	vec pt2=(*framelist.begin()).pt2;
	vec pt3=(*++framelist.begin()).pt2;
	vec pt4=(*++framelist.begin()).pt1;
	PBVec.clear();
	for (it=ptcls.begin();it!=ptcls.end();++it){
		vec posi=(*it)->getCurrPosition();
		vec vdt=(posi-pt1)%normalize(norm)*normalize(norm);
		vec proj=posi-vdt;
		REAL dist=vfabs(vdt);
		if (dist>(*it)->getMaxRadius()&&vdt%norm<0)
			continue;
		vec v1=pt1-proj;
		vec v2=pt2-proj;
		vec v3=pt3-proj;
		vec v4=pt4-proj;
		if(((v1*v2)%norm)*((v2*v3)%norm)*((v3*v4)%norm)*((v4*v1)%norm)<0){
			continue;
		}
		PBVec.push_back(*it);
	}
};

template<class T>
void plnflb_bdry<T>::update(UPDATECTL ctl[], unsigned int len){
	std::vector<LINE>::iterator it;
	int i=0;
	if (framelist.size()!=len){
		perror("in plnflb_bdry::update: not enough information for update");
		exit(-1);
	}
	for (it=framelist.begin();it!=framelist.end();++it,i++){
		(*it).update(ctl[i]);
	}
	vec bch1=(*framelist.begin()).pt1-(*framelist.begin()).pt2;
	vec bch2=(*++framelist.begin()).pt1-(*framelist.begin()).pt2;
	vec nm=normalize(bch1*bch2);
	if(nm%norm<0)
		nm*=-1;
	norm=nm;
};

template<class T>
void plnflb_bdry<T>::findParticleOnLine(){
	//in z dircetion, the net is nz-1 grid and nz node
	vec pt1=(*framelist.begin()).pt1;
	vec pt2=(*framelist.begin()).pt2;
	vec pt3=(*++framelist.begin()).pt1;
	vec pt4=(*++framelist.begin()).pt2;
	int iz, ip;
	typename std::vector<T*>::iterator it;
	for (iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			PCFB[iz][ip].clear();
			vec ip_top=pt1+ip*(pt2-pt1)/(np-1);
			vec ip_bot=pt3+ip*(pt4-pt3)/(np-1);
			vec thept=ip_bot+(ip_top-ip_bot)/(nz-1)*iz;
			for (it=PBVec.begin();it!=PBVec.end();++it){
				vec v0=(*it)->getCurrPosition();
				REAL dist=vfabs((v0-thept)-(v0-thept)%normalize(norm)*normalize(norm));
				if(dist<(*it)->getMaxRadius())
					PCFB[iz][ip].push_back(*it);
			}
		}
	}
};

template<class T>
void plnflb_bdry<T>::delNull(){
	int ip=0,iz=0;
	int kp,kz;
	int dobrk=0;
	for (iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			if(RelatedP[iz][ip]!=NULL){
				dobrk=1;
				break;
			}
		}
		if(dobrk)
			break;
	}
	if(RelatedP[0][0]==NULL)
		RelatedP[0][0]=RelatedP[iz][ip];
	for(kz=0;kz<nz;kz++){
		for(kp=0;kp<np;kp++){
			if(RelatedP[kz][kp]==NULL){
				if(kp!=0)
					RelatedP[kz][kp]=RelatedP[kz][kp-1];
				else
					RelatedP[kz][kp]=RelatedP[kz-1][kp];
			}
		}
	}
};

template<class T>
void plnflb_bdry<T>::createFlbNet(){
	int iz, ip;
	typename std::vector<T*>::iterator it;
	vec pt1=(*framelist.begin()).pt1;
	vec pt2=(*framelist.begin()).pt2;
	vec pt3=(*++framelist.begin()).pt1;
	vec pt4=(*++framelist.begin()).pt2;
	/*  the order of pt1 pt2, pt3, pt4 must conform to the normal of the plane, that is
		pt1-----pt2
		 |\      |
		 |  \    |
		 |    \  |
		pt3-----pt4
		if the norm is shooting out of the plane;
	*/
	for (iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			vec ip_top=pt1+(pt2-pt1)/(np-1)*ip;
			vec ip_bot=pt3+(pt4-pt3)/(np-1)*ip;
			vec thept=ip_bot+(ip_top-ip_bot)/(nz-1)*iz;
			vec outmost=thept;
			T *op=NULL;
			T *bp=NULL;
			REAL otmst=-1.0e16;
			for (it=PCFB[iz][ip].begin();it!=PCFB[iz][ip].end();++it){
				if(norm%((*it)->getCurrPosition()-thept)>otmst){
					otmst=norm%((*it)->getCurrPosition()-thept);
					bp=*it;
				}
				vec rt[2];
				if((*it)->intersectWithLine(thept,normalize(norm),rt,1)){
					//v1.print();rt[0].print();rt[1].print();getchar();
					vec tp=(rt[0]-rt[1])%norm>0?rt[0]:rt[1];
					if ((tp-outmost)%norm>0){
						outmost=tp;
						op=*it;
					}
				}
			}

			FlxbNet[iz][ip]=outmost;
			if(bp!=NULL){
				if(op!=NULL)RelatedP[iz][ip]=op;
				else RelatedP[iz][ip]=bp;
			}else
				RelatedP[iz][ip]=NULL;
			if(op!=NULL)
				op->IsFBP=true;
		}
	}
	delNull();
};

template<class T>
void plnflb_bdry<T>::flxbBF(){
	int iz, ip;
	vec p[3];
	T* e[3];
	sumpressure=0;
	for (iz=1;iz<nz;iz++){
		for (ip=0;ip<np-1;ip++){
			vec p1=FlxbNet[iz-1][ip];
			vec p2=FlxbNet[iz-1][ip+1];
			vec p3=FlxbNet[iz][ip+1];
			vec p4=FlxbNet[iz][ip];
			T* e1=RelatedP[iz-1][ip];
			T* e2=RelatedP[iz-1][ip+1];
			T* e3=RelatedP[iz][ip+1];
			T* e4=RelatedP[iz][ip];
			p[0]=p1;p[1]=p2;p[2]=p4;
			e[0]=e1;e[1]=e2;e[2]=e4;
			sumpressure+=this->triangleDstr(confining,-norm,p,e);
			p[0]=p2;p[1]=p3;p[2]=p4;
			e[0]=e2;e[1]=e3;e[2]=e4;
			sumpressure+=this->triangleDstr(confining,-norm,p,e);
		}
	}
};

template <class T> class cylflb_bdry:public flb_bdry<T>{
public:
	vec sumpressure;
	REAL confining;
	std::vector<CIRC> framelist;// top and bottom frame, can only have two elements
	REAL alf;               // the expand from pt1 ro pt2; for a complete cylinder alf=2Pi
	int framenum;             // =2
	int side;                 // 1, the particles are outside cylinder; -1, inside
	std::vector<T*> PBVec;
	vec FlxbNet[100][50];
	T* RelatedP[100][50];
	std::vector<T*> PCFB[100][50];
	int np, nz;
public:
	cylflb_bdry(std::ifstream &ifs);
	virtual ~cylflb_bdry() {}; // base class needs a virtual destructor.
	void disp() const;
	void findParticleOnBoundary(std::vector<T*>& ptcls);
	void delNull();
	void findParticleOnLine();          // create possible particles per line
	void createFlbNet();
	void flxbBF();
	void update(UPDATECTL ctl[], unsigned int len);
};

template<class T>
void cylflb_bdry<T>::disp() const{
	int iz,ip,i;
	cout << "sumpressure: " << sumpressure.getx() << " " << sumpressure.gety() << " " << sumpressure.getz() << endl;
	cout << "confining pressure: " << confining << endl;
	cout << "side: " << side << " framenum: " << framenum << " np: " << np << " nz: " << nz << " alf: " << alf << endl;
	std::vector<CIRC>::const_iterator it;
	for(it=framelist.begin();it!=framelist.end();++it)
		(*it).disp();
	typename std::vector<T*>::const_iterator jt;
	i=0;
	for(jt=PBVec.begin();jt!=PBVec.end();++jt){
		if(i++<10)
			cout << (*jt)->getID();
		else{
			i=0;
			cout << (*jt)->getID() << endl;
		}
	}
	cout << "point of the net" << endl;
	for(iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
		  cout << FlxbNet[iz][ip].getx() << " " << FlxbNet[iz][ip].gety() << " " << FlxbNet[iz][ip].getz();
		}
		cout << endl;
	}
	cout << "particles of the net" << endl;
	for(iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			if(i++<10) {
				if(RelatedP[iz][ip]!=NULL) cout << RelatedP[iz][ip]->getID();
			}
			else{
				i=0;
				if(RelatedP[iz][ip]!=NULL) cout << RelatedP[iz][ip]->getID() << endl;
			}
		}
	}
	cout << "possible particles around a line" << endl;
	for(iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			for (jt=PCFB[iz][ip].begin();jt!=PCFB[iz][ip].end();++jt)
				if((*jt)!=NULL) cout << (*jt)->getID();
			cout << endl;
		}
	}
};

template<class T>
cylflb_bdry<T>::cylflb_bdry(std::ifstream &ifs){
	CIRC tmp;
	REAL x,y,z;
	int i,j;

	framelist.clear();
	ifs >> framenum >> nz >> np >> side >> confining;
	for(i=0;i<framenum;i++){
	    ifs >> x >> y >> z;
	    vec ipt=vec(x,y,z);
	    tmp.center=ipt;
	    ifs >> x >> y >> z;
	    ipt=vec(x,y,z);
	    tmp.norm=ipt;
	    ifs >> x >> y >> z;
	    ipt=vec(x,y,z);
	    tmp.pt1=ipt;
	    ifs >> x >> y >> z;
	    ipt=vec(x,y,z);
	    tmp.pt2=ipt;
	    ifs >> tmp.turn >> tmp.radius;
	    framelist.push_back(tmp);
	}
	for (i=0;i<nz;i++){
	    for(j=0;j<np;j++){
		RelatedP[i][j]=NULL;
		PCFB[i][j].clear();
	    }
	}
};

template<class T>
void cylflb_bdry<T>::findParticleOnBoundary(std::vector<T*>&ptcls){
	typename std::vector<T*>::iterator it;
	vec ct1=framelist.begin()->center;
	vec ct2=(++framelist.begin())->center;
	vec ml1=framelist.begin()->norm;
	vec ml2=(++framelist.begin())->norm;
	REAL r1=framelist.begin()->radius;
	REAL r2=(++framelist.begin())->radius;
	vec pt1=framelist.begin()->pt1;
	vec pt2=framelist.begin()->pt2;
	vec pt3=(++framelist.begin())->pt1;
	vec pt4=(++framelist.begin())->pt2;
	int turn=framelist.begin()->turn;
	
	if (vfabs(ml1*ml2)>1.0e-5*vfabs(ml1)||fabs(r1-r2)>1.0e-5*r1){
		perror("in cylflb_bdry::findParticleOnBoundary: the two CIRC do not build a cylinder");
		exit(-1);
	}
	if (vfabs((pt1-pt3)*ml1)>1.0e-8||
		vfabs((pt2-pt4)*ml1)>1.0e-8){
		perror("in cylflb_bdry::findParticleOnBoundary: the end points are not good");
		exit(-1);
	}

	alf=angle(pt1-ct1,pt2-ct1,turn*ml1);
	
	PBVec.clear();
	vec ct=(ct1+ct2)/2;
	vec norm=ml1;
	REAL rad=r1;
	for (it=ptcls.begin();it!=ptcls.end();++it){
		vec posi=(*it)->getCurrPosition();
		vec proj=posi-ct-(posi-ct)%normalize(norm)*normalize(norm);
		REAL dist=vfabs(proj);
		if ((dist<0.8*rad&&side==-1)||(dist>1.2*rad&&side==1))
			continue;
		if(vfabs(pt1-pt2)>1.0e-8){// a uncomplete circle
			REAL bta=angle(pt1-ct1,proj,turn*ml1);
			if (bta>alf)
				continue;
		}
		PBVec.push_back(*it);
	}
};

template<class T>
void cylflb_bdry<T>::findParticleOnLine(){
	//in z dircetion, the net is nz-1 grid and nz node
	vec pt1=(*framelist.begin()).pt1;
	vec ct1=framelist.begin()->center;
	vec ct2=(++framelist.begin())->center;
	vec nm=framelist.begin()->norm;
	int turn=framelist.begin()->turn;
	vec rote=turn*alf/(nz-1)*normalize(nm);

	int iz, ip;
	typename std::vector<T*>::iterator it;
	for (iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			PCFB[iz][ip].clear();
			vec ll=rotateVec(pt1-ct1,rote*ip);
			vec ip_top=ct1+ll;
			vec ip_bot=ct2+ll;
			vec thept=ip_bot+iz*(ip_top-ip_bot)/(nz-1);
			for (it=PBVec.begin();it!=PBVec.end();++it){
				vec v0=(*it)->getCurrPosition();
				REAL dist=vfabs((v0-thept)-(v0-thept)%normalize(ll)*normalize(ll));
				if(dist<(*it)->getMaxRadius())
					PCFB[iz][ip].push_back(*it);
			}
		}
	}
};

template<class T>
void cylflb_bdry<T>::delNull(){
	int ip=0,iz=0;
	int kp,kz;
	int dobrk=0;
	for (iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			if(RelatedP[iz][ip]!=NULL){
				dobrk=1;
				break;
			}
		}
		if(dobrk)
			break;
	}
	if(RelatedP[0][0]==NULL)
		RelatedP[0][0]=RelatedP[iz][ip];
	for(kz=0;kz<nz;kz++){
		for(kp=0;kp<np;kp++){
			if(RelatedP[kz][kp]==NULL){
				if(kp!=0)
					RelatedP[kz][kp]=RelatedP[kz][kp-1];
				else
					RelatedP[kz][kp]=RelatedP[kz-1][kp];
			}
		}
	}
};

template<class T>
void cylflb_bdry<T>::createFlbNet(){
	int iz, ip;
	typename std::vector<T*>::iterator it;

	vec pt1=(*framelist.begin()).pt1;
	vec ct1=framelist.begin()->center;
	vec ct2=(++framelist.begin())->center;
	vec nm=framelist.begin()->norm;
	int turn=framelist.begin()->turn;
	vec rote=alf*turn*normalize(nm)/(nz-1);

	for (iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			vec ll=rotateVec(pt1-ct1,rote*ip);
			vec ip_top=ct1+ll;
			vec ip_bot=ct2+ll;
			vec thept=ip_bot+iz*(ip_top-ip_bot)/(nz-1);
			ll*=-side;//ll point out of assembly
			vec outmost=thept;
			T *op=NULL;
			for (it=(PCFB[iz][ip]).begin();it!=(PCFB[iz][ip]).end();++it){
				vec rt[2];
				if((*it)->intersectWithLine(thept,normalize(ll),rt,1)){
					//v1.print();rt[0].print();rt[1].print();getchar();
					vec tp=(rt[0]-rt[1])%ll>0?rt[0]:rt[1];
					if ((tp-outmost)%ll>0){
						outmost=tp;
						op=*it;
					}
				}
			}
			FlxbNet[iz][ip]=outmost;
			RelatedP[iz][ip]=op;
			if(op!=NULL)
 				op->IsFBP=true;
		}
	}
	delNull();
};

template<class T>
void cylflb_bdry<T>::flxbBF(){
	int iz, ip;
	vec p[3];
	T* e[3];
	vec ct1=framelist.begin()->center;
	vec nm=normalize(framelist.begin()->norm);

	sumpressure=0;
	for (iz=1;iz<nz;iz++){
		for (ip=0;ip<np-1;ip++){
			
			vec p1=FlxbNet[iz-1][ip];
			vec p2=FlxbNet[iz-1][ip+1];
			vec p3=FlxbNet[iz][ip+1];
			vec p4=FlxbNet[iz][ip];
			T* e1=RelatedP[iz-1][ip];
			T* e2=RelatedP[iz-1][ip+1];
			T* e3=RelatedP[iz][ip+1];
			T* e4=RelatedP[iz][ip];

			p[0]=p1;p[1]=p2;p[2]=p4;
			e[0]=e1;e[1]=e2;e[2]=e4;
			vec tricnt=(p[0]+p[1]+p[2])/3;
			vec trinm=(tricnt-ct1)-(tricnt-ct1)%nm*nm;
			trinm*=side;
			sumpressure+=this->triangleDstr(confining,trinm,p,e);
			p[0]=p2;p[1]=p3;p[2]=p4;
			e[0]=e2;e[1]=e3;e[2]=e4;
			sumpressure+=this->triangleDstr(confining,trinm,p,e);
		}
	}
};

template<class T>
void cylflb_bdry<T>::update(UPDATECTL ctl[], unsigned int len){
	std::vector<CIRC>::iterator it;
	int i=0;
	if (framelist.size()!=len){
		perror("in plnflb_bdry::update: not enough information for update");
		exit(-1);
	}
	for (it=framelist.begin();it!=framelist.end();++it,i++){
		(*it).update(ctl[i]);
	}
};

} // namespace dem ends

#endif
