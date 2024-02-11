#ifndef SURFACETENSION_H
#define SURFACETENSION_H

#include "realtypes.h"
#include "parameters.h"
#include "matrix.h"
#include "globfuncs.h"
#include "node.h"
#include "activepoint.h"
#include "ele2dcoupled.h"
#include <vector>
#include <list>
#include <string>


namespace memFluid {


class surfacetension{

private:
    // global variables
    REAL Lx, Ly, element_length;
    std::vector<node*> node9;
    std::vector<node*> node4;
    std::vector<ele2dcoupled*> element;
    std::list<activepoint*> Active_points; 
    std::vector<node*> topNodes;
    std::vector<node*> botNodes;
    std::vector<node*> rightNodes;
    std::vector<node*> leftNodes;
    int numnode9, numnode4;
    REAL dx9, dy9;

    std::string plain_axi;

    int element_x, element_y;
    REAL R2;
    REAL c2[2];
    int close_DBC2;
    REAL a, b;

    int source;
    int quadorder;
    REAL tol;

    // construct mesh
    REAL pt1[2], pt2[2], pt3[2], pt4[2];	// four corner points
    
    // computer number of nodes, of elements
    int numelem9, numelem4;


    // boundary nodes
//    std::vector<int> topNodes4, botNodes4, rightNodes4, leftNodes4;
//    std::vector<int> middletop4, middlebot4, middleleft4, middleright4;
    std::vector<node*> markers;	// for plotting

    // get enriched nodes
    std::vector<ele2dcoupled*> split_elem2;
    int numsnode92, numsnode42;

    std::list<ele2dcoupled*> Split_ordered;	// list is convenient for unique
    std::vector<node*> connect2D;		// this can be done by sorting node9 & renumbering the num_node9
    std::vector<REAL> nodex;	// x coordinate of Ordered_points
    std::vector<REAL> nodey;	// y coordinate of Ordered_points
    REAL A0;

    // get winner node no which to discretize lagrange multipliers
    matrix Wmoes2, Moes_Assembly2;

    // Initialize C & K matrix, force vector
    int total_unknown;
    matrix v, s;
    matrix F, Fm2x, xFm2, Fm2y, yFm2, Fm1x, Fm1y; // global
    matrix K;	// global matrix
    matrix dsol;

    // initialize K F

    // I, J, X are used to construct sparse matrix,
    // I, J are for numbering, X is the value
    std::vector<int> I;
    std::vector<int> J;
    std::vector<REAL> X;	// global matrix

public:
    surfacetension(){
	Lx = 0; Ly = 0; element_length = 0;
	numnode9 = 9; numnode4 = 0;
	dx9 = 0; dy9 = 0;
	element_x = 0; element_y = 0;
    	R2 = 0;
    	c2[0] = 0; c2[1] = 0;
    	close_DBC2 = 0;
    	a = 0; b = 0;

    	source = 10;
    	quadorder = 2;
    	tol = 1e-05;
	pt1[0] = 0; pt1[1] = 0; pt2[0] = 0; pt2[1] = 0; 
	pt3[0] = 0; pt3[1] = 0; pt4[0] = 0; pt4[1] = 0;

 	numelem9 = 0; numelem4 = 0;
	numsnode92 = 0; numsnode42 = 0;
	A0 = 0;
	total_unknown = 0;
    }

    ~surfacetension(){

     	// it is important to release memory pointed to by pointers in the container,
    	// otherwise memory leaking occurs
    	for(std::vector<node*>::iterator it=node9.begin(); it<node9.end(); ++it)
	    delete (*it);
    	for(std::vector<node*>::iterator it=node4.begin(); it<node4.end(); ++it)
	    delete (*it);
    	for(std::vector<ele2dcoupled*>::iterator it2d=element.begin(); it2d<element.end(); ++it2d)
	    delete (*it2d);
    	for(std::list<activepoint*>::iterator it=Active_points.begin(); it!=Active_points.end(); ++it)
	    delete (*it);

    	for(std::vector<node*>::iterator it=topNodes.begin(); it<topNodes.end(); ++it)
	    delete (*it);
    	for(std::vector<node*>::iterator it=botNodes.begin(); it<botNodes.end(); ++it)
	    delete (*it);
    	for(std::vector<node*>::iterator it=rightNodes.begin(); it<rightNodes.end(); ++it)
	    delete (*it);
    	for(std::vector<node*>::iterator it=leftNodes.begin(); it<leftNodes.end(); ++it)
	    delete (*it);
    	for(std::vector<node*>::iterator it=markers.begin(); it<markers.end(); ++it)
	    delete (*it);

    	for(std::vector<ele2dcoupled*>::iterator it2d=split_elem2.begin(); it2d<split_elem2.end(); ++it2d)
	    delete (*it2d);
    	for(std::list<ele2dcoupled*>::iterator it2d=Split_ordered.begin(); it2d!=Split_ordered.end(); ++it2d)
	    delete (*it2d);

    	for(std::vector<node*>::iterator it=connect2D.begin(); it!=connect2D.end(); ++it)
	    delete (*it);

    	// in case of consecutive simulations
	node9.clear(); node4.clear(); element.clear(); Active_points.clear();
	topNodes.clear(); botNodes.clear(); rightNodes.clear(); leftNodes.clear(); markers.clear();
	split_elem2.clear(); Split_ordered.clear();
	connect2D.clear();
    }

    void solve(const std::string &plain_axi, const char* filename, int interval);	
    // solve the problem, interval is for print out
    
    void constructMesh();	// construct mesh
    void levelSet();		// find the level set
    void boundaryNodes();	// find the boundary nodes
    void initGridParticle();	// initialize grid based particles
    void initLevelSet();	// Init level-set/elements intersection points
    void gaussPoints(const int it);		// Gauss points & weights
    void applyBdryConditions(); // apply boundary conditions


    // this function will generate the mesh based on the pt1, pt2, pt3, pt4 and the arguments
    // finally, it will create elements in node9 and element9 or node4 and element4 based on the string
    // if string is "Q9", then 9-node element, if "Q4", then 4-node element
    // and get numelem9 and numelem4	
    REAL distancePointToEllipse(REAL, REAL, REAL, REAL, REAL, REAL);

    matrix velocity_pt_fluid(matrix &p_proj, const std::string in_out);
    REAL velocity_pt_membrane(matrix &p_proj);

    // update particles position
    void find_rotation();	// will call activepoints::calcOmega_e()

    // RESAMPLE FOOTPOINTS AND CALCULATE GEOMETRICAL QUANTITIES
    void Foot_point_ressampling_v2();
    void find_laplace_beltrami_coeff();
    void E_F_J_T_coeff_updating();
    void sort_points();	// create Split_ordered, 
    void find_proj_normal_deriv_pt_parallel(matrix &, matrix &, matrix &, matrix &);
    void find_velocity_gradient_d(matrix &d,matrix &d_cont,matrix &E,matrix &Fd,REAL &Ja, matrix &p);

    // get level-set/elements
    bool isInSplit_ordered(node* in) const;
    bool isInW_update(node* in) const;

    void meshRectangularRegion(int, int);
    void square_node_array(int nnx, int nny);
    void make_elem(matrix &node_pattern, int numx, int numy, int inc_u, int inc_v);

    void sort2D_xfm();	// need Ly, Lx, element_length_mat
    void Init_Foot_points_func();	// create Active_points
    void setFoot_points_n_Foot_points_actpnt(activepoint* pt_act);

    void changeIsnanFoot_points();

    void printFoot_points_T_plot(const char* str) const;
    void printOrdered_point(const char* str) const;
    void printAratio(const char* stepsfp, const int nAratio, 
		     std::vector<REAL>::const_iterator Aratio_beg, 
		     std::vector<REAL>::const_iterator Aratio_end ) const;
    void printFigureVariables(const char* file_name, const int it) const;
    void printNode9(const char* str) const;
    void printElement9(const char* str) const;
    void printMarkers(const char* str) const;
    void printDsol(const char* str) const;
    void printOrdered_xFm(const char* str) const;
    void printFoot_proj(const char* str) const;

    void printF_debug(const char* str) const;
    void printK_debug(const char* str) const;
    void printDsol_debug(const char* str) const;
    void readDsol_debug();

};






} // end of memFluid


#endif
