#include "node.h"
#include <algorithm>
#include <list>
#include <vector>


namespace memFluid{

void node::init(){

    num_node9 = -1;	
    num_node4 = -1;	

    topNodes = false;
    botNodes = false;
    rightNodes = false;
    leftNodes = false;
    border_nodes = false;

    act_point = false;

    Phi91 = 0;
    Phi41 = 0;
    Phi92 = 0;
    Phi42 = 0;
//    REAL Phi92s; // for ellipse

    pos92 = -1;	// the number of this node in enrich node 92
    pos42 = -1;  // initially, -1
//    int pos91;
//    int pos41; 	// never used in matlab

    enrich_node92 = false;
    enrich_node42 = false;
    enrich_node41 = false;
    enrich_node91 = false;

    Grad_phi_x = 0;
    Grad_phi_y = 0;


} // init()


node::node(REAL x, REAL y){

    nodex = x; nodey = y;
    init();

} // node()


void node::normalR2(std::vector<node*>::iterator connect2D_beg, 
		    std::vector<node*>::iterator connect2D_end, const REAL Ly, const REAL Lx, const REAL element_length){

    REAL h = element_length*0.5;

    int ny = round(Ly/h+1);
    int nx = round(Lx/h+1);

    int a = 1;
    h = element_length*a*0.5;

    node* self = (this);
    std::vector<node*>::iterator i0j0 = std::find(connect2D_beg, connect2D_end, self); 
    // the same address means the same node

    REAL dx = element_length*0.5;
    REAL dy = element_length*0.5;

    if(rightNodes && botNodes){

	std::vector<node*>::iterator i0_plus = i0j0+nx;
	std::vector<node*>::iterator j0_minus = i0j0-1;
	Grad_phi_y = ( (*i0_plus)->getPhi92()-(*i0j0)->getPhi92() )/dy;
	Grad_phi_x = ( (*i0j0)->getPhi92()-(*j0_minus)->getPhi92() )/dx;

    }
    else if(rightNodes && topNodes){

	std::vector<node*>::iterator i0_minus = i0j0-nx;
	std::vector<node*>::iterator j0_minus = i0j0-1;
	Grad_phi_y = ( (*i0j0)->getPhi92()-(*i0_minus)->getPhi92() )/dy;
	Grad_phi_x = ( (*i0j0)->getPhi92()-(*j0_minus)->getPhi92() )/dx;

    }
    else if(topNodes && leftNodes){

	std::vector<node*>::iterator i0_minus = i0j0-nx;
	std::vector<node*>::iterator j0_plus = i0j0+1;
	Grad_phi_y = ( (*i0j0)->getPhi92()-(*i0_minus)->getPhi92() )/dy;
	Grad_phi_x = ( (*j0_plus)->getPhi92()-(*i0j0)->getPhi92() )/dx;

    }
    else if(botNodes && leftNodes){

	std::vector<node*>::iterator i0_plus = i0j0+nx;
	std::vector<node*>::iterator j0_plus = i0j0+1;
	Grad_phi_y = ( (*i0_plus)->getPhi92()-(*i0j0)->getPhi92() )/dy;
	Grad_phi_x = ( (*j0_plus)->getPhi92()-(*i0j0)->getPhi92() )/dx;

    }
    else if(botNodes){

	std::vector<node*>::iterator i0_plus = i0j0+nx;
	std::vector<node*>::iterator j0_plus = i0j0+1;
	std::vector<node*>::iterator j0_minus = i0j0-1;
	Grad_phi_y = ( (*i0_plus)->getPhi92()-(*i0j0)->getPhi92() )/dy;
	Grad_phi_x = ( (*j0_plus)->getPhi92()-(*j0_minus)->getPhi92() )/dx*0.5;

    }
    else if(rightNodes){

	std::vector<node*>::iterator i0_plus = i0j0+nx;
	std::vector<node*>::iterator i0_minus = i0j0-nx;
	std::vector<node*>::iterator j0_minus = i0j0-1;
	Grad_phi_y = ( (*i0_plus)->getPhi92()-(*i0_minus)->getPhi92() )/dy*0.5;
	Grad_phi_x = ( (*i0j0)->getPhi92()-(*j0_minus)->getPhi92() )/dx;

    }
    else if(topNodes){

	std::vector<node*>::iterator i0_minus = i0j0-nx;
	std::vector<node*>::iterator j0_plus = i0j0+1;
	std::vector<node*>::iterator j0_minus = i0j0-1;
	Grad_phi_y = ( (*i0j0)->getPhi92()-(*i0_minus)->getPhi92() )/dy;
	Grad_phi_x = ( (*j0_plus)->getPhi92()-(*j0_minus)->getPhi92() )/dx*0.5;

    }
    else if(leftNodes){

	std::vector<node*>::iterator i0_plus = i0j0+nx;
	std::vector<node*>::iterator i0_minus = i0j0-nx;
	Grad_phi_y = ( (*i0_plus)->getPhi92()-(*i0_minus)->getPhi92() )/dy*0.5;
	Grad_phi_x = 0;

    }
    else {

	std::vector<node*>::iterator i0_plus  = i0j0+nx;
	std::vector<node*>::iterator i0_minus = i0j0-nx;
	std::vector<node*>::iterator j0_plus  = i0j0+1;
	std::vector<node*>::iterator j0_minus = i0j0-1;
	Grad_phi_y = ( (*i0_plus)->getPhi92()-(*i0_minus)->getPhi92() )/dy*0.5;
	Grad_phi_x = ( (*j0_plus)->getPhi92()-(*j0_minus)->getPhi92() )/dx*0.5;

    }


} // normalR2()

/*
bool operator < (const node* A, const node* B){

    if(A->getY() < B->getY()){
	return true;
    }
    else if(A->getY()==B->getY()) {
	if(A->getX() < B->getX()){
	    return true;
	}
	else if(A->getX()==B->getX()) {
	    return false;
	}	
	else {
	    return false;
	}
    }
    else {
	return false;
    }

} // <
*/

// non member function, for sorting
bool less_than(node* A, node* B){

    if(A->getY() < B->getY()){
	return true;
    }
    else if(A->getY()==B->getY()) {
	if(A->getX() < B->getX()){
	    return true;
	}
	else if(A->getX()==B->getX()) {
	    return false;
	}	
	else {
	    return false;
	}
    }
    else {
	return false;
    }

} // less_than()




} // end of memFluid
