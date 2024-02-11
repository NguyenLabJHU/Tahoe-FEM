#ifndef NODE_H
#define NODE_H

#include "matrix.h"
#include "globfuncs.h"
#include "realtypes.h"
#include "parameters.h"


namespace memFluid{

class node{
friend class ele2dcoupled;
friend class activepoint;


private:
    REAL nodex, nodey;	// coords
    int num_node9;	// the number in node9, initially -1
    int num_node4;	// the number in node4, initially -1

    bool topNodes;
    bool botNodes;
    bool rightNodes;
    bool leftNodes;
    bool border_nodes;

    bool act_point;	// if active point

    REAL Phi91;
    REAL Phi41;
    REAL Phi92;
    REAL Phi42;
//    REAL Phi92s; // for ellipse
    int pos92;	// the number of this node in enrich node 92
    int pos42;  // initially, -1
//    int pos91;
//    int pos41; 	// never used in matlab
    bool enrich_node92;
    bool enrich_node42;
    bool enrich_node41;
    bool enrich_node91;

    REAL Grad_phi_x;
    REAL Grad_phi_y;

    void init();

public:
    node(REAL x, REAL y);

    void setTopNodes(bool is) {topNodes = is;}
    void setBotNodes(bool is) {botNodes = is;}
    void setRightNodes(bool is) {rightNodes = is;}
    void setLeftNodes(bool is) {leftNodes = is;}
    bool getTopNodes() const {return topNodes;}
    bool getBotNodes() const {return botNodes;}
    bool getRightNodes() const {return rightNodes;}
    bool getLeftNodes() const {return leftNodes;}
    bool getAct_point() const {return act_point;}

    int getNum_node9() const {return num_node9;}
    int getNum_node4() const {return num_node4;}
    void setNum_node9(const int num) {num_node9 = num;}
    void setNum_node4(const int num) {num_node4 = num;}
    REAL getPhi91() const {return Phi91;}
    REAL getPhi92() const {return Phi92;}
    REAL getPhi42() const {return Phi42;}
    REAL getPhi41() const {return Phi41;}
    REAL getGrad_phi_x() const {return Grad_phi_x;}
    REAL getGrad_phi_y() const {return Grad_phi_y;}
    int  getPos92() const {return pos92;}
    int  getPos42() const {return pos42;}
    void setEnrich_node42(bool is) {enrich_node42 = is;}
    void setEnrich_node92(bool is) {enrich_node92 = is;}
    void setEnrich_node41(bool is) {enrich_node41 = is;}
    void setEnrich_node91(bool is) {enrich_node91 = is;}
    void initPos42() {pos42 = 0;}
    void initPos92() {pos92 = 0;}
    bool getEnrich_node42() const {return enrich_node42;}
    bool getEnrich_node92() const {return enrich_node92;}
    void setPos92(int i) {pos92 = i;}
    void setPos42(int i) {pos42 = i;}
    REAL getX() const {return nodex;}
    REAL getY() const {return nodey;}
    void setPhi92(REAL val) {Phi92 = val;}
    void setPhi42(REAL val) {Phi42 = val;}
    void setPhi91(REAL val) {Phi91 = val;}
    void setPhi41(REAL val) {Phi41 = val;}

    void normalR2(std::vector<node*>::iterator , std::vector<node*>::iterator ,
		  const REAL, const REAL, const REAL);	

    // non member functions
//    friend /*inline */bool operator < (const node*, const node*);

//    friend bool less_than(node* , node*);	// sort by y & x

};

bool less_than(node * , node *);	// sort by y & x

} // end of memFluid



#endif
