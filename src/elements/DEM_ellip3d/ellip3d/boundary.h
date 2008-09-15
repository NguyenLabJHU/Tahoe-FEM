//     1.This is a template class, for which we have to include the implementation in the header file.
//       As we can't put using statement in a header file, we have to use std::xxx wherever we need to
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

#include <map>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

#include "vec.h"
#include "cylinder.h"
#include "boundarytgt.h"

namespace dem {

// BdryCoef is used for rigid boundary conditions
typedef struct bdryfunc{
	int order;  // 1-linear; 2-quadratic
	vec dirc;   // normal vector if plane, mother line vector if cylinder,it points out of the particles		
	vec apt;    // a point on the plane or a point on the axis of the cylinder
	long double rad; //zero if plane
	int side;   //zero if plane; side=1, particles are outside the cylinder; =-1, inside the cylinder
        void disp() const{
		printf("the elements of a boundary are\n");
		printf("order:%5d\n",order);
		printf("dirc: %10.6Lf%10.6Lf%10.6Lf\n",dirc.getx(),dirc.gety(),dirc.getz());
		printf("apt : %10.6Lf%10.6Lf%10.6Lf\n",apt.getx(),apt.gety(),apt.getz());
		printf("radius %10.6Lf and side %5d:\n",rad,side);
	}
	void disp(std::ofstream &ofs) const{
	    ofs<<std::setw(10)<<order
	       <<std::setw(16)<<dirc.getx()
	       <<std::setw(16)<<dirc.gety()
	       <<std::setw(16)<<dirc.getz()
	       <<std::setw(16)<<apt.getx()
	       <<std::setw(16)<<apt.gety()
	       <<std::setw(16)<<apt.getz()
	       <<std::setw(16)<<rad
	       <<std::setw(10)<<side<<std::endl;
	}
} BdryCoef;

// LINE and CIRC are structs used in flexible boundary classes
typedef struct updatectl{
	vec tran;    // tranlation second
	vec rote;    // rotate first
	vec fixpt;   // before any update is made
	long double expnd;// expand last
	updatectl(){tran=0;rote=0;fixpt=0;expnd=1;}
	void disp() const{
		printf("tran: %10.6Lf%10.6Lf%10.6Lf\n",tran.getx(),tran.gety(),tran.getz());
		printf("rote: %10.6Lf%10.6Lf%10.6Lf\n",rote.getx(),rote.gety(),rote.getz());
		printf("fixpt: %10.6Lf%10.6Lf%10.6Lf\n",fixpt.getx(),fixpt.gety(),fixpt.getz());
		printf("expand: %Lf\n", expnd);
	};
}UPDATECTL;

typedef struct line{
	vec pt1; // begining point
	vec pt2; // ending point
	void disp() const{
		printf("pt1: %10.6Lf%10.6Lf%10.6Lf\n",pt1.getx(),pt1.gety(),pt1.getz());
		printf("pt2: %10.6Lf%10.6Lf%10.6Lf\n",pt2.getx(),pt2.gety(),pt2.getz());
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
	long double radius; //the radius of the circle
	vec pt1;    // the begining point of a part circle
	vec pt2;    // the end point of a part circle;it pt1==pt2, a closure circle
	void disp() const {
		printf("center: %10.6Lf%10.6Lf%10.6Lf\n",center.getx(),center.gety(),center.getz());
		printf("norm: %10.6Lf%10.6Lf%10.6Lf\n",norm.getx(),norm.gety(),norm.getz());
		printf("pt1: %10.6Lf%10.6Lf%10.6Lf\n",pt1.getx(),pt1.gety(),pt1.getz());
		printf("pt2: %10.6Lf%10.6Lf%10.6Lf\n",pt2.getx(),pt2.gety(),pt2.getz());
		printf("turn:%5d radius:%10.6Lf\n",turn,radius);
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
	int         bdry_id;    // the first record defines the bdry itself, the other 
	std::vector<BdryCoef> CoefOfLimits; // limitnum records define the other lines on the bdry 
	long double avg_normal;             // that give limits to the first boundary.
	long double avg_penetr; // average penetration by particles onto this boundary
	int         cntnum;     // contact numbers by particles onto this boundary
	long double area;       // the bounary's area
	int         limitnum;   // how many lines the boundary have
public:
	rgd_bdry(std::ifstream &ifs);
	int getBdryID() {return bdry_id;}
	virtual ~rgd_bdry() {} // base class needs a virtual destructor.
	virtual void disp() const{
		std::vector<BdryCoef>::const_iterator it;
		for(it=CoefOfLimits.begin();it!=CoefOfLimits.end();++it)
			(*it).disp();
		printf("area: %Lf limitnum: %d\n",area,limitnum);
	}
	virtual void disp(std::ofstream &ofs) const{
		std::vector<BdryCoef>::const_iterator it;
		ofs <<std::endl
		    <<std::setw(10)<<(*CoefOfLimits.begin()).order<<std::endl;
		ofs <<std::setw(10)<<bdry_id
		    <<std::setw(10)<<limitnum
		    <<std::setw(16)<<area <<std::endl;
		for(it=CoefOfLimits.begin();it!=CoefOfLimits.end();++it)
			(*it).disp(ofs);
	}
	virtual void createPBL(std::list<T*>& ptcls){};
	virtual void rigidBF(std::map<int,std::vector<boundarytgt> >& BdryTgtMap)
	    {std::cout<<"parent"<<std::endl;} // calculate for each boundary particles the rigid boundary force
	virtual vec getNormal() const{return 0;}
	virtual long double getAvgNormal() const{return 0;}
	virtual long double getAvgPenetr() const{return 0;}
	virtual int         getCntnum() const{return 0;}
	virtual vec getTangt() const{return 0;}
	virtual vec getApt() const{return 0;}
	virtual vec getDirc() const{return 0;}
	virtual void setArea(long double a){area=a;}
	virtual long double getArea(){return area;}
	virtual void update(UPDATECTL& ctl); //the boundary is translating with tran and rotating with rote around axis
};

template <class T>
rgd_bdry<T>::rgd_bdry(std::ifstream &ifs){
	BdryCoef tmp;
	long double x,y,z;
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
	virtual void createPBL(std::list<T*>& ptcls){};
	virtual void update(UPDATECTL ctl[], unsigned int len){};
	virtual void createPLL(){}; // create possible particles per line
	virtual void createFlbNet(){};
	virtual void flxbBF(){};
	virtual vec triangleDstr(long double pressure,vec norm, vec p[], T* e[]); //norm is the direction of pressure
	virtual ~flb_bdry() {};     // base class needs a virtual destructor.
};

template<class T>
vec flb_bdry<T>::triangleDstr(long double pressure, vec norm, vec p[], T* e[]){
        //norm indicates the pressure dirction
	vec cent=(p[0]+p[1]+p[2])/3;
	long double l1=vfabsl(p[1]-p[0]);
	long double l2=vfabsl(p[2]-p[1]);
	long double l3=vfabsl(p[0]-p[2]);
	long double hp=(l1+l2+l3)/2;
	long double area=sqrtl(hp)*sqrtl(hp-l1)*sqrtl(hp-l2)*sqrtl(hp-l3);
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
	std::list<T*> PBList; // possible boundary particles of this specific boundary
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
	long double distToBdry(vec posi) const;
	void createPBL(std::list<T*>& ptcls);
	vec getApt() const;
	vec getDirc() const;
	plnrgd_bdry<T>* getBdry(int bdryid) const{
		return this;
	}
	void rigidBF(std::map<int,std::vector<boundarytgt> >& BdryTgtMap);
	vec getNormal() const{return normal;}
	long double getAvgNormal() const{return this->avg_normal;}
	long double getAvgPenetr() const{return this->avg_penetr;}
	int         getCntnum() const{return this->cntnum;}

	vec getTangt() const{return tangt;}
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
	printf("normal: %20.6Lf%20.6Lf%20.6Lf\n", 
		normal.getx(),normal.gety(),normal.getz());
	typename std::list<T*>::const_iterator it;
	int i=0;
	for(it=PBList.begin();it!=PBList.end();++it){
		if(i++<10)
			printf("%5d",(*it)->getID());
		else{
			i=0;
			printf("%5d\n",(*it)->getID());
		}
	}
};

template<class T>
long double plnrgd_bdry<T>::distToBdry(vec posi) const{
	vec dv=(*this->CoefOfLimits.begin()).dirc;
	vec pt=(*this->CoefOfLimits.begin()).apt;
	vec ndv=normalize(dv);
	return (posi-pt)%ndv;
};

template<class T>
void plnrgd_bdry<T>::createPBL(std::list<T*>& ptcls){
    typename std::list<T*>::iterator it;
    std::vector<BdryCoef>::iterator bt;
    bool next;
    PBList.clear();
    long double dist, r;
    vec posi, ndirc;
    for (it=ptcls.begin();it!=ptcls.end();++it){
	if ( (*it)->getType() == 0 ) { // only process free particles, excluding type 5
	    posi=(*it)->getCurrPosition();
	    dist=distToBdry(posi);
	    if(dist>=0 || fabsl(dist) > (*it)->getA()) // outside to CoefOfLimits[0] or inside too much
		continue;
	    next=true;
	    //g_exceptioninf<<std::setw(10)<<g_iteration <<std::setw(10)<<getBdryID()<<std::setw(10)<<(*it)->getID();
	    for (bt=++this->CoefOfLimits.begin();bt!=this->CoefOfLimits.end();++bt){ // CoefOfLimits[1,2,...]
		ndirc=normalize((*bt).dirc);
		r=vfabsl((posi-(*bt).apt)-(posi-(*bt).apt)%ndirc*ndirc);
		if((*bt).order==1 && (posi-(*bt).apt)%(*bt).dirc >= 0 ||
		   (*bt).order==2 && (r-(*bt).rad)*(*bt).side<0){
		    next=false; // the particle is out of boundary, process next particle
		    break;
		}
	    }
	    if(next)
		PBList.push_back(*it);
	}
    }

};

/*
template<class T>
void plnrgd_bdry<T>::createPBL(std::list<T*>& ptcls){
	typename std::list<T*>::iterator it;
	std::vector<BdryCoef>::iterator bt;
	bool next;
	PBList.clear();
 	for (it=ptcls.begin();it!=ptcls.end();++it){
		vec posi=(*it)->getCurrPosition();
		long double dist=distToBdry(posi);
		//if (fabsl(dist)>ROOM*(*it)->getA()&&dist<0)
		if(fabsl(dist)>ROOM*(*it)->getA()||dist>0)
			continue;
		next=false;
		for (bt=++CoefOfLimits.begin();bt!=CoefOfLimits.end();++bt){
			vec ndirc=normalize((*bt).dirc);
			long double r=abs((posi-(*bt).apt)-(posi-(*bt).apt)%ndirc*ndirc);
			if((*bt).order==1&&(posi-(*bt).apt)%(*bt).dirc>0||
				(*bt).order==2&&(r-(*bt).rad)*(*bt).side<0){
				next=true;//the particle is outof boundary, process next particle
				break;
			}
		}
		if(!next)
			PBList.push_back(*it);
	}
};
*/


template<class T>
void plnrgd_bdry<T>::rigidBF(std::map<int,std::vector<boundarytgt> >& BdryTgtMap){
    typename std::list<T*>::iterator it;
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
    long double penetr=0;
    int count=0;
    for (it=PBList.begin();it!=PBList.end();++it){
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
	std::list<T*> PBList;
public:
	cylrgd_bdry(std::ifstream &ifs):rgd_bdry<T>(ifs){normal=0;}
	void disp() const;
	long double distToBdry(vec posi) const;
	void createPBL(std::list<T*>& ptcls);
	void rigidBF();
	vec getNormal() const{return normal;};
};

template<class T>
void cylrgd_bdry<T>::disp() const{
	rgd_bdry<T>::disp();
	printf("normal: %Lf %Lf %Lf\n", normal.getx(), normal.gety(), normal.getz());
	typename std::list<T*>::const_iterator it;
	int i=0;
	for(it=PBList.begin();it!=PBList.end();++it){
		if(i++<10)
			printf("%5d",(*it)->getID());
		else{
			i=0;
			printf("%5d\n",(*it)->getID());
		}
	}
};

template<class T>
long double cylrgd_bdry<T>::distToBdry(vec posi) const{
	vec ndc=(*this->CoefOfLimits.begin()).dirc;
	vec napt=(*this->CoefOfLimits.begin()).apt;
	long double r=(*this->CoefOfLimits.begin()).rad;
	vec norm=normalize(ndc);
	return fabsl(r-vfabsl((posi-napt)-(posi-napt)%norm*norm));
};

template<class T>
void cylrgd_bdry<T>::createPBL(std::list<T*> &ptcls){
	typename std::list<T*>::iterator it;
	std::vector<BdryCoef>::iterator bt;
	bool next;
	PBList.clear();
	long double dist,r;
	vec posi, ndirc;
 	for (it=ptcls.begin();it!=ptcls.end();++it){
		posi=(*it)->getCurrPosition();
		dist=distToBdry(posi);
		if (dist>(*it)->getA())
			continue;
		next=false;
		for (bt=++this->CoefOfLimits.begin();bt!=this->CoefOfLimits.end();++bt){
			ndirc=normalize((*bt).dirc);
			r=vfabsl((posi-(*bt).apt)-(posi-(*bt).apt)%ndirc*ndirc);
			if((*bt).order==1&&(posi-(*bt).apt)%(*bt).dirc>(*it)->getA()||
				(*bt).order==2&&(r-(*bt).rad)*(*bt).side<0){
				next=true;//the particle is outof boundary, process next particle
				break;
			}
		}
		if(!next)
			PBList.push_back(*it);
	}
};

template<class T>
void cylrgd_bdry<T>::rigidBF(){
	// I am temporially saitisfied with the cylinder with vertical mother line
	cylinder cyl;
	typename std::list<T*>::iterator it;
	BdryCoef tmp;
	tmp=*this->CoefOfLimits.begin();
	cyl.set_radius(tmp.rad);
	cyl.set_center(tmp.apt);
	normal=0;
	for (it=PBList.begin();it!=PBList.end();++it){
		normal-=(*it)->cylinderRBForce(this->bdry_id,cyl,tmp.side);
	}
};

template<class T> class plnflb_bdry:public flb_bdry<T>{
public:
	vec sumpressure; // sum of total water pressure on particles
	long double confining;// confining pressure by surrounding liquid
	std::list<LINE> framelist;//store rigid lines it
	vec norm;        // normal direction pointing outward the assembly
	int framenum;    // how many rigid lines there are, can only be 2
	std::list<T*> PBList;
	vec FlxbNet[100][50];
	T* RelatedP[100][50];
	std::list<T*> PCFB[100][50];
	int np, nz;
public:
	plnflb_bdry(std::ifstream &ifs);
	virtual ~plnflb_bdry() {}; // base class needs a virtual destructor.
	void disp() const;
	void createPBL(std::list<T*>& ptcls);
	void createPLL(); // create possible particles per line
	void createFlbNet();
	void flxbBF();    // FlxbNet[nz][np], RelatedP[nz][np]; if side=1, particle is in the side of >0, side=-1, <0
	void update(UPDATECTL ctl[], unsigned int len);
	void delNull();
};

template<class T>
void plnflb_bdry<T>::disp() const{
	int iz,ip,i;
	printf("sumpressure:%20.6Lf%20.6Lf%20.6Lf\n",
		sumpressure.getx(),sumpressure.gety(),sumpressure.getz());
	printf("norm:%10.6Lf%10.6Lf%10.6Lf\n",
		norm.getx(),norm.gety(),norm.getz());
	printf("np:%5dnz:%5dframenum:%5d\n",np,nz,framenum);
	printf("confining pressure:%20.6Lf\n",confining);
	std::list<LINE>::const_iterator it;
	for(it=framelist.begin();it!=framelist.end();++it)
		(*it).disp();
	typename std::list<T*>::const_iterator jt;
	i=0;
	for(jt=PBList.begin();jt!=PBList.end();++jt){
		if(i++<10)
			printf("%5d",(*jt)->getID());
		else{
			i=0;
			printf("%5d\n",(*jt)->getID());
		}
	}
	printf("\npoint of the net\n");
	for(iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			printf("%10.6Lf%10.6Lf%10.6Lf",
				FlxbNet[iz][ip].getx(),FlxbNet[iz][ip].gety(),FlxbNet[iz][ip].getz());
		}
		printf("\n\n");
	}
	i=0;
	printf("\nparticles of the net\n");
	for(iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			if(i++<10)
				if(RelatedP[iz][ip]!=NULL)printf("%5d",RelatedP[iz][ip]->getID());
			else{
				i=0;
				if(RelatedP[iz][ip]!=NULL)printf("%5d\n",RelatedP[iz][ip]->getID());
			}
		}
	}
	printf("\npossible particles around a line\n");
	for(iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			for (jt=PCFB[iz][ip].begin();jt!=PCFB[iz][ip].end();++jt)
				if((*jt)!=NULL)printf("%5d",(*jt)->getID());
			printf("\n");
		}
	}
};

template<class T>
plnflb_bdry<T>::plnflb_bdry(std::ifstream &ifs){
	vec tmp;
	long double x,y,z;
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
void plnflb_bdry<T>::createPBL(std::list<T*>& ptcls){
	/*
		1----2
		|    |
		|	 |
		4----3
	*/
	typename std::list<T*>::iterator it;
	vec pt1=(*framelist.begin()).pt1;
	vec pt2=(*framelist.begin()).pt2;
	vec pt3=(*++framelist.begin()).pt2;
	vec pt4=(*++framelist.begin()).pt1;
	PBList.clear();
	for (it=ptcls.begin();it!=ptcls.end();++it){
		vec posi=(*it)->getCurrPosition();
		vec vdt=(posi-pt1)%normalize(norm)*normalize(norm);
		vec proj=posi-vdt;
		long double dist=vfabsl(vdt);
		if (dist>(*it)->getA()&&vdt%norm<0)
			continue;
		vec v1=pt1-proj;
		vec v2=pt2-proj;
		vec v3=pt3-proj;
		vec v4=pt4-proj;
		if(((v1*v2)%norm)*((v2*v3)%norm)*((v3*v4)%norm)*((v4*v1)%norm)<0){
			continue;
		}
		PBList.push_back(*it);
	}
};

template<class T>
void plnflb_bdry<T>::update(UPDATECTL ctl[], unsigned int len){
	std::list<LINE>::iterator it;
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
void plnflb_bdry<T>::createPLL(){
	//in z dircetion, the net is nz-1 grid and nz node
	vec pt1=(*framelist.begin()).pt1;
	vec pt2=(*framelist.begin()).pt2;
	vec pt3=(*++framelist.begin()).pt1;
	vec pt4=(*++framelist.begin()).pt2;
	int iz, ip;
	typename std::list<T*>::iterator it;
	for (iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			PCFB[iz][ip].clear();
			vec ip_top=pt1+ip*(pt2-pt1)/(np-1);
			vec ip_bot=pt3+ip*(pt4-pt3)/(np-1);
			vec thept=ip_bot+(ip_top-ip_bot)/(nz-1)*iz;
			for (it=PBList.begin();it!=PBList.end();++it){
				vec v0=(*it)->getCurrPosition();
				long double dist=vfabsl((v0-thept)-(v0-thept)%normalize(norm)*normalize(norm));
				if(dist<(*it)->getA())
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
	typename std::list<T*>::iterator it;
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
			long double otmst=-1.0e16;
			for (it=PCFB[iz][ip].begin();it!=PCFB[iz][ip].end();++it){
				if(norm%((*it)->getCurrPosition()-thept)>otmst){
					otmst=norm%((*it)->getCurrPosition()-thept);
					bp=*it;
				}
				vec rt[2];
				if((*it)->intersectWithLine(thept,normalize(norm),rt)){
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
			//fprintf(net,"%15.8lf%15.8lf%15.8lf\n",FlxbNet[iz][ip].x, FlxbNet[iz][ip].y,FlxbNet[iz][ip].z);
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
			sumpressure+=triangleDstr(confining,-norm,p,e);
			p[0]=p2;p[1]=p3;p[2]=p4;
			e[0]=e2;e[1]=e3;e[2]=e4;
			sumpressure+=triangleDstr(confining,-norm,p,e);
		}
	}
};

template <class T> class cylflb_bdry:public flb_bdry<T>{
public:
	vec sumpressure;
	long double confining;
	std::list<CIRC> framelist;// top and bottom frame, can only have two elements
	long double alf;               // the expand from pt1 ro pt2; for a complete cylinder alf=2Pi
	int framenum;             // =2
	int side;                 // 1, the particles are outside cylinder; -1, inside
	std::list<T*> PBList;
	vec FlxbNet[100][50];
	T* RelatedP[100][50];
	std::list<T*> PCFB[100][50];
	int np, nz;
public:
	cylflb_bdry(std::ifstream &ifs);
	virtual ~cylflb_bdry() {}; // base class needs a virtual destructor.
	void disp() const;
	void createPBL(std::list<T*>& ptcls);
	void delNull();
	void createPLL();          // create possible particles per line
	void createFlbNet();
	void flxbBF();
	void update(UPDATECTL ctl[], unsigned int len);
};

template<class T>
void cylflb_bdry<T>::disp() const{
	int iz,ip,i;
	printf("sumpressure:%10.6Lf%10.6Lf%10.6Lf\n",
		sumpressure.getx(),sumpressure.gety(),sumpressure.getz());
	printf("confining pressure:%20.6Lf\n",confining);
	printf("side:%5d, framenum: %5d, np:%5d, nz:%5d,alf: %10.6Lf\n",
		side,framenum,np,nz,alf);
	std::list<CIRC>::const_iterator it;
	for(it=framelist.begin();it!=framelist.end();++it)
		(*it).disp();
	typename std::list<T*>::const_iterator jt;
	i=0;
	for(jt=PBList.begin();jt!=PBList.end();++jt){
		if(i++<10)
			printf("%5d",(*jt)->getID());
		else{
			i=0;
			printf("%5d\n",(*jt)->getID());
		}
	}
	printf("\npoint of the net\n");
	for(iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			printf("%10.6Lf%10.6Lf%10.6Lf",
				FlxbNet[iz][ip].getx(),FlxbNet[iz][ip].gety(),FlxbNet[iz][ip].getz());
		}
		printf("\n\n");
	}
	printf("\nparticles of the net\n");
	for(iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			if(i++<10)
				if(RelatedP[iz][ip]!=NULL)printf("%5d",RelatedP[iz][ip]->getID());
			else{
				i=0;
				if(RelatedP[iz][ip]!=NULL)printf("%5d\n",RelatedP[iz][ip]->getID());
			}
		}
	}
	printf("\npossible particles around a line\n");
	for(iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			for (jt=PCFB[iz][ip].begin();jt!=PCFB[iz][ip].end();++jt)
				if((*jt)!=NULL)printf("%5d",(*jt)->getID());
			printf("\n");
		}
	}
};

template<class T>
cylflb_bdry<T>::cylflb_bdry(std::ifstream &ifs){
	CIRC tmp;
	long double x,y,z;
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
void cylflb_bdry<T>::createPBL(std::list<T*>&ptcls){
	typename std::list<T*>::iterator it;
	vec ct1=framelist.begin()->center;
	vec ct2=(++framelist.begin())->center;
	vec ml1=framelist.begin()->norm;
	vec ml2=(++framelist.begin())->norm;
	long double r1=framelist.begin()->radius;
	long double r2=(++framelist.begin())->radius;
	vec pt1=framelist.begin()->pt1;
	vec pt2=framelist.begin()->pt2;
	vec pt3=(++framelist.begin())->pt1;
	vec pt4=(++framelist.begin())->pt2;
	int turn=framelist.begin()->turn;
	
	if (vfabsl(ml1*ml2)>1.0e-5*vfabsl(ml1)||fabsl(r1-r2)>1.0e-5*r1){
		perror("in cylflb_bdry::createPBL: the two CIRC do not build a cylinder");
		exit(-1);
	}
	if (vfabsl((pt1-pt3)*ml1)>1.0e-8||
		vfabsl((pt2-pt4)*ml1)>1.0e-8){
		perror("in cylflb_bdry::createPBL: the end points are not good");
		exit(-1);
	}

	alf=angle(pt1-ct1,pt2-ct1,turn*ml1);
	
	PBList.clear();
	vec ct=(ct1+ct2)/2;
	vec norm=ml1;
	long double rad=r1;
	for (it=ptcls.begin();it!=ptcls.end();++it){
		vec posi=(*it)->getCurrPosition();
		vec proj=posi-ct-(posi-ct)%normalize(norm)*normalize(norm);
		long double dist=vfabsl(proj);
		if ((dist<0.8*rad&&side==-1)||(dist>1.2*rad&&side==1))
			continue;
		if(vfabsl(pt1-pt2)>1.0e-8){// a uncomplete circle
			long double bta=angle(pt1-ct1,proj,turn*ml1);
			if (bta>alf)
				continue;
		}
		PBList.push_back(*it);
	}
};

template<class T>
void cylflb_bdry<T>::createPLL(){
	//in z dircetion, the net is nz-1 grid and nz node
	vec pt1=(*framelist.begin()).pt1;
	vec ct1=framelist.begin()->center;
	vec ct2=(++framelist.begin())->center;
	vec nm=framelist.begin()->norm;
	int turn=framelist.begin()->turn;
	vec rote=turn*alf/(nz-1)*normalize(nm);

	int iz, ip;
	typename std::list<T*>::iterator it;
	for (iz=0;iz<nz;iz++){
		for(ip=0;ip<np;ip++){
			PCFB[iz][ip].clear();
			vec ll=rotateVec(pt1-ct1,rote*ip);
			vec ip_top=ct1+ll;
			vec ip_bot=ct2+ll;
			vec thept=ip_bot+iz*(ip_top-ip_bot)/(nz-1);
			for (it=PBList.begin();it!=PBList.end();++it){
				vec v0=(*it)->getCurrPosition();
				long double dist=vfabsl((v0-thept)-(v0-thept)%normalize(ll)*normalize(ll));
				if(dist<(*it)->getA())
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
	typename std::list<T*>::iterator it;

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
				if((*it)->intersectWithLine(thept,normalize(ll),rt)){
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
			sumpressure+=triangleDstr(confining,trinm,p,e);
			p[0]=p2;p[1]=p3;p[2]=p4;
			e[0]=e2;e[1]=e3;e[2]=e4;
			sumpressure+=triangleDstr(confining,trinm,p,e);
		}
	}
};

template<class T>
void cylflb_bdry<T>::update(UPDATECTL ctl[], unsigned int len){
	std::list<CIRC>::iterator it;
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
