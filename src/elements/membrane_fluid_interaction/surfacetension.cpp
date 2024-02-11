#include "surfacetension.h"
#include "array.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <complex>
#include <math.h>
#include <map>
#include <algorithm>
#include <iomanip>


#if defined(_OPENMP)
#include <omp.h>
#else
inline int omp_get_thread_num() { return 0;}
inline int omp_get_num_threads() { return 1;}
#endif

// id is the id of thread, n is problem size, p is the number of threads
#define BLOCK_LOW(id, p, n)    ((id)*(n)/(p))
#define BLOCK_HIGH(id, p, n)   (BLOCK_LOW((id)+1, p, n)-1)
#define BLOCK_SIZE(id, p, n)   (BLOCK_LOW((id)+1, p, n)-BLOCK_LOW(id, p, n))


namespace memFluid{

void surfacetension::solve(const std::string& axi, const char* file_name, int interval){

    plain_axi = axi;

    if(plain_axi=="axi1") {
	    // Dimensin of the domain (it is simply a rectangular region Lx x Ly)
	    element_y = 60;
	    element_x = 40;

	    element_length = 0.1;

   	    Ly = element_y*element_length*0.5;	// y
	    Lx = element_x*element_length*0.5;	// x

	    R2 = 1.01;
	    c2[0] = -0.5*Lx;
	    c2[1] = 0;

	    close_DBC2 = 0;

    } // if "axi1"

    if(plain_axi=="axi2"){
	    // Dimension of the domain (it is simply a rectangular region Lx x Ly)
	    element_y = 60;
	    element_x = 30;
	
	    element_length = 0.1;

   	    Ly = element_y*element_length*0.5;	// y
	    Lx = element_x*element_length*0.5;	// x

	    a = Lx-0.221+10;
	    b = Ly*0.5-0.121;
	    c2[0] = -0.5*Lx;
	    c2[1] = -0.4;

    } // if "axi2"

    if(plain_axi=="plain"){
	    // Dimension of the domain (it is simply a rectangular region Lx x Ly)
	    element_y = 80;
	    element_x = 180;
	
	    element_length = 0.1;

	    Ly = element_y*element_length*0.5; // y
	    Lx = element_x*element_length*0.5; // x

	    R2 = Ly*0.5-0.521;
	    c2[0] = 0;
	    c2[1] = 0;

//	    a = 0.5*Lx-0.221;
//	    b = 0.5*Ly-0.321;

	    close_DBC2 = 1;
    } // if "plain"

    source = 10;
    quadorder = 2;
    tol = 1e-05;

    // Material properties    
    L_PLUS = 0.00000001*R2;
    L_MINUS = 0.00000001*R2;

    MU_PLUS = MU_F_PLUS/L_PLUS;
    MU_PLUS = MU_F_MINUS/L_MINUS;


    // floor position
    YP = c2[1] - R2 - 0.1;

    constructMesh();

    // LEVEL SET
    levelSet();

    // BOUNDARY NODES            
    boundaryNodes();

    //  initialize Grid Based Particle 
    initGridParticle();

    // Init level-set/elements intersection points 
    initLevelSet();

    // start time loop  
    int nit = 0;
    int nit2 = 0;
    int nit3 = 0;
    int nAratio = 0;
    REAL dt = 1e-6;
    int it = 1;
    // print information
    int  stepsnum=0;
    char stepsstr[4];
    char stepsfp[50];

    REAL errorFJ = 1.0e10;
    REAL errorV = 1.0e10;

    // start nonlinear iterations
    while (errorFJ>1e-09 || errorV>0.001*element_length ){

	nit = nit+1;
	nit2 = nit2+1;
	nit3 = nit3+1;
	it = it+1;

   	// get winner nodes on which to discretize lagrange multipliers (MOES)2
	Wmoes2.clear();
	Moes_Assembly2.clear();

	// GAUSS POINTS & WEIGHTS  
 	gaussPoints(it);

	// Initialize C & K matrix, force vector
	total_unknown = numnode9*2 + numnode4 + numsnode92 + numsnode42 + 4*Split_ordered.size()+4;


	//initialize external forces (tension (neuman), velocity (dirichlet), sources)
	matrix t(2,1);
	t(1,1) = 0; t(2,1) = 0;
	
	v = zeros(total_unknown,1);
	
        s = zeros(total_unknown ,1);

	//initialize K,F
	F = zeros(total_unknown,1);
   	Fm2x = zeros(numnode9,1); xFm2 = zeros(numnode9,1);
    	Fm2y = zeros(numnode9,1); yFm2 = zeros(numnode9,1);
    	Fm1x = zeros(numnode9,1);
    	Fm1y = zeros(numnode9,1);
	int ng = 0;
	int total_triplet = 0;

//#pragma omp parallel for
	for(std::vector<ele2dcoupled*>::iterator it2d=element.begin(); it2d<element.end(); ++it2d){
	    if( (*it2d)->getSplit_ordered()>0 ){
		(*it2d)->setEleUnknowns(27+8+2+2+2+2);
		(*it2d)->setNtriplets(total_triplet);
		total_triplet = total_triplet + (27+8+2+2+2+2)*(27+8+2+2+2+2);
	    }
	    else {
		int partial_enriched = (*it2d)->getPartial_enriched();
		int partial_enriched4 = (*it2d)->getPartial_enriched4();

		(*it2d)->setEleUnknowns(18 + 4 + partial_enriched + partial_enriched4);
		(*it2d)->setNtriplets(total_triplet);
		total_triplet = total_triplet + (18+4+partial_enriched+partial_enriched4)*(18+4+partial_enriched+partial_enriched4);
	    }
	}

	for(int i_I=0; i_I<total_triplet; ++i_I){
	    I.push_back(0);
	}

	for(int i_J=0; i_J<total_triplet; ++i_J){
	    J.push_back(0);
	}

	for(int i_X=0; i_X<total_triplet; ++i_X){
	    X.push_back(0);
	}


	std::cout << it << " STIFFNESS MATRIX COMPUTATION " << std::endl;

	//-----------loop on elements ----------------------

	int element_size = element.size();

#pragma omp parallel
  {
	int num_threads = omp_get_num_threads();
	int thread_num  = omp_get_thread_num();
    #pragma omp single
	std::cout << "number of threads: " << num_threads << std::endl;
//	std::cout << "thread number: " << thread_num << std::endl;

	int start_num = BLOCK_LOW(thread_num, num_threads, element_size);
	int end_num   = BLOCK_HIGH(thread_num, num_threads, element_size);

//	for(std::vector<ele2dcoupled*>::iterator it2d=element.begin(); it2d<element.end(); ++it2d){
	for(int it_num=start_num; it_num<end_num+1; ++it_num){

	    (element[it_num])->initIJXFelemConnect();	// initialize Ii, Jj, Xx, Felem, connect

	    // sctrBv, sctrBl_plus, sctrBl_minus, sctrBv_m 
	    // sctrBp, sctrBlp are in ele2dcoupled
	    (element[it_num])->assembly_fluid_XFEM_ndiscont_v2_parallel(numnode9, numsnode92, numnode4,
							 numsnode42, Split_ordered.size(), close_DBC2);


	    // sctrBTOT is variable in ele2dcoupled
	    (element[it_num])->calcSctrBTOT();
	
	    (element[it_num])->calcNn9Nn4();	// nn9 = length(sctrBv)


	    (element[it_num])->setAllUnknowns();
	    (element[it_num])->setAllIndex();

	    //initialize Ke, Fe
	    (element[it_num])->initKeFe();
	    (element[it_num])->initSeVe(s, v);

	    (element[it_num])->initAllKs();
	    (element[it_num])->initAllFs();

	    (element[it_num])->calcT_elementDtdx();

	    (element[it_num])->calculateKsFs(Lx, Ly, plain_axi, nodex.begin(), nodex.end(), nodey.begin(), nodey.end());	   



	    if( (element[it_num])->getSplit_ordered()>0 ){ // split element

		(element[it_num])->calcSplitKsFs(Lx, Ly, plain_axi, Active_points.begin(), Active_points.end());
	    } 


	    //------------------
	    // assembly, element
	    //------------------

	    (element[it_num])->assembly();	// this is element matrices assembly

	    // assembly, global
	    int element_unknowns = (element[it_num])->getEleUnknowns();
	    int ntriplets = (element[it_num])->getNtriplets();
	    matrix Ii_temp = (element[it_num])->getIi();
	    matrix Jj_temp = (element[it_num])->getJj();
	    matrix Xx_temp = (element[it_num])->getXx();
	    for(int i=ntriplets+1; i<ntriplets+element_unknowns*element_unknowns+1; ++i){

		I[i-1] = Ii_temp(1,i-ntriplets);
		J[i-1] = Jj_temp(1,i-ntriplets);
		X[i-1] = Xx_temp(1,i-ntriplets);		
	    }

//	    matrix connect_temp = (*it2d)->getConnect();
//	    matrix Felem_temp = (*it2d)->getFelem();
//	    for(int i=1; i<element_unknowns+1; ++i){
//
//		F(connect_temp(1,i), 1) = F(connect_temp(1,i), 1) + Felem_temp(1, i);
//	    }
	
     	} // end of element for loop

  } // end omp parallel

	// this cannot be parallelized
	for(std::vector<ele2dcoupled*>::iterator it2d=element.begin(); it2d<element.end(); ++it2d){

	    int element_unknowns = (*it2d)->getEleUnknowns();
	    matrix connect_temp = (*it2d)->getConnect();
	    matrix Felem_temp = (*it2d)->getFelem();
	    for(int i=1; i<element_unknowns+1; ++i){

		F(connect_temp(1,i), 1) = F(connect_temp(1,i), 1) + Felem_temp(1, i);
	    }
  	}


	K = zeros(total_unknown, total_unknown);
#pragma omp parallel for
	for(int i=0; i<total_triplet; ++i){	// total_triplet = length(X)
	    K(I[i], J[i]) = X[i];
	}



	//-----apply boundary condition-----

	applyBdryConditions();


	//-----------------------------------//
	//  solve for velocity and pressure  //
	//-----------------------------------//
	std::cout << it << " solve for velocity and pressure" << std::endl;

/*	
	// get parameters for sparseMatrixSolver
	std::vector<REAL> val_K;
	std::vector<int> row_index;	// 0-based
	std::vector<int> col_start;	// 0-based
	for(int jj_index=1; jj_index<total_unknown+1; ++jj_index){
	    bool is_first = true;
	    for(int ii_index=jj_index; ii_index<total_unknown+1; ++ii_index){
		int val_ij = K(ii_index,jj_index);
		if(fabs(val_ij)>1e-12){	// non-zeros
		    val_K.push_back(val_ij);
		    row_index.push_back(ii_index-1);
		    if(is_first){ // first non-zero, starts a column
			col_start.push_back(val_K.size()-1);
		    }
		    is_first = false;
		    break;
		}
	    }
	}
	col_start.push_back(val_K.size());
	
*/


	matrixEqnSolver(dsol, K, F);	// extern function
//	sparseMatrixSolver(dsol, F, total_unknown, val_K.begin(), val_K.end(), 
//			   row_index.begin(), row_index.end(), col_start.begin(), col_start.end() );


// get dsol from matlab
//readDsol_debug();






	REAL max_v = 0;
	for(int i=1; i!=2*numnode9+1; ++i){
	    REAL abs_v = fabs(dsol(i,1));
	    if(max_v < abs_v)
		max_v = abs_v;
	}
	errorV = max_v*dt;
	std::cout << "errorV: " << errorV << std::endl;

//	REAL max_x, max_y;
	for(std::list<activepoint*>::iterator itpt=Active_points.begin(); itpt!=Active_points.end(); ++itpt){
 	    
	    (*itpt)->initVel_FootsVel_membranesX_projsV();

	    matrix pt = (*itpt)->getFoot_points();

	    matrix p_proj, t, dtdx;
	    find_proj_normal_deriv_pt_parallel(p_proj, t, dtdx, pt);

	    matrix v = velocity_pt_fluid(p_proj, "in");
	    (*itpt)->setVel_Foot_points_in(v);

	    v = velocity_pt_fluid(p_proj, "out");
	    (*itpt)->setVel_Foot_points_out(v);

	    REAL v_m = velocity_pt_membrane(p_proj);
 	    (*itpt)->calcVel_Foot_membrane_normal_tangent(v_m, t);

	    (*itpt)->setX_proj_m(p_proj);

	} // end of Active_points loop


//	REAL errorM = pow(max_x*max_x + max_y*max_y, 0.5); // never used



/*
	//---------------//
	//   PLOTTING    //
	//---------------//
	if(it==2 || nit2==interval){
	    nit = 0;

	    sprintf(stepsstr, "%03d", it); 
	    strcpy(stepsfp, file_name); strcat(stepsfp, "_figure_1_Ordered_point_"); strcat(stepsfp, stepsstr);
	    printOrdered_point(stepsfp);
	    
	}

	if(it==2 || nit3==interval){
	    nit3 = 0;
	
	    nAratio = nAratio+1;
	    std::vector<REAL> Aratio;
	    if(plain_axi=="axi1"){
	
		    REAL A = 0;
		    for(std::list<ele2dcoupled*>::const_iterator it2d=Split_ordered.begin(); 
		                                                 it2d!=Split_ordered.end(); ++it2d){
			matrix node_mat = (*it2d)->getOrdered_point();			
			REAL Rzm = (node_mat(1,1)+0.5*Lx+node_mat(2,1)+0.5*Lx)*0.5;
			REAL dz = node_mat(2,2)-node_mat(1,2);
			A = A + PI*(Rzm*Rzm*dz);
		    }
		
		    Aratio.push_back(A0/A);
  	    } // if "axi1"

	    if(plain_axi=="plain"){
		    REAL Aratio_val = polyarea(nodex.begin(), nodex.end(), nodey.begin(), nodey.end());
		    Aratio_val = A0/Aratio_val;
		    Aratio.push_back(Aratio_val);

	    } // if "plain"

	    sprintf(stepsstr, "%03d", it); 
	    strcpy(stepsfp, file_name); strcat(stepsfp, "_figure_2_Aratio_"); strcat(stepsfp, stepsstr);
	    printAratio(stepsfp, nAratio, Aratio.begin(), Aratio.end());

   	}

	if(it==2 || nit==interval){
	
	    nit = 0;

	    printFigureVariables(file_name, it);

	} // end if

*/


	//----------------------------//
	// UPDATE PARATICLES POSITION //
	//----------------------------//
	std::cout << it << " update particles" << std::endl;

	for(std::list<activepoint*>::iterator itpt=Active_points.begin();
	    itpt!=Active_points.end(); ++itpt){
	
	    (*itpt)->initAllDts();
	    (*itpt)->calcFoot_points_dt2(dt);
	    (*itpt)->calcFoot_points_dt2_tangent(dt);
	}

	find_rotation();

	for(std::list<activepoint*>::iterator itpt=Active_points.begin();
	    itpt!=Active_points.end(); ++itpt){

	    matrix pt = (*itpt)->getFoot_points_dt2();
	    matrix pt_tang = (*itpt)->getFoot_points_dt2_tangent();

	    matrix p_proj2, t, dtdx;
	    find_proj_normal_deriv_pt_parallel(p_proj2, t, dtdx, pt_tang);

	    int in_on;
	    //in_on = InPolygon2(pt(1,1), pt(1,2), nodex.begin(), nodex.end(), nodey.begin(), nodey.end());
	    in_on = inpoly(pt(1,1), pt(1,2), nodex.begin(), nodex.end(), nodey.begin(), nodey.end());	    

	    std::string in_out;
	    if(in_on==1)
		in_out = "in";
	    else
		in_out = "out";

	    matrix v_n = velocity_pt_fluid(pt, in_out);	// v_n is 2x1

	    REAL v_m = velocity_pt_membrane(p_proj2);

	    (*itpt)->updateFoot_points(v_n, v_m, dt);

	} // end of activepoints for

	//----------------------------------------------------------//
	// RESAMPLE FOOTPOINTS AND CALCULATE GEOMETRICAL QUANTITIES //	
	//----------------------------------------------------------//
	Foot_point_ressampling_v2();
	find_laplace_beltrami_coeff();

	for(std::list<activepoint*>::iterator itpt=Active_points.begin(); itpt!=Active_points.end(); ++itpt){

	    (*itpt)->initT_Foot_points();

	    matrix p = (*itpt)->getFoot_points();
	    matrix d, d_cont, E, Fd;
	    REAL Ja;
	    find_velocity_gradient_d(d, d_cont, E, Fd, Ja, p);	// d, d_cont, E, Fd are all 1x2
	    
	    (*itpt)->updateFEJT_Foot_points(d, d_cont, E, Fd, Ja, dt);

	    (*itpt)->updateT_plot();	    
	}

/*	// print out
	if(it%interval==0){
	    sprintf(stepsstr, "%03d", stepsnum); 
	    strcpy(stepsfp, file_name); strcat(stepsfp, "_Foot_points_T_plot_"); strcat(stepsfp, stepsstr);
	    printFoot_points_T_plot(stepsfp);

	    ++stepsnum;
	}
*/


	E_F_J_T_coeff_updating();



 	//--------------------//
	// get enriched nodes //
	//--------------------//
	split_elem2.clear();	// split_elem2 is a vector

	// inialize enrich_node92 and enrich_node42 for all nodes
	for(std::vector<node*>::iterator itn=node9.begin(); itn!=node9.end(); ++itn){
	    (*itn)->setEnrich_node92(false);
	    (*itn)->setEnrich_node42(false);	// node9 > node4
	}

	for(std::vector<ele2dcoupled*>::iterator it2d=element.begin(); it2d!=element.end(); ++it2d){
	    matrix phi2(1,4);
	
	    (*it2d)->setSplit_elem2(false);	// initialize

	    matrix phi2_9 = (*it2d)->getPhi92();
	    phi2(1,1) = phi2_9(1,1); phi2(1,2) = phi2_9(1,2);
	    phi2(1,3) = phi2_9(1,3); phi2(1,4) = phi2_9(1,4);
	    // no need to update center_node since, in ele2dcoupled, node* center_node;
	    if( max(phi2)*min(phi2) < 0 ){

		(*it2d)->setSplit_elem2(true);
		(*it2d)->setEnrich_node42(true);
		(*it2d)->setEnrich_node92(true);
		split_elem2.push_back(*it2d);
	    }
	}

	// corresponding94 not used in matlab code
	
	numsnode92 = 0;
	int nsnode92 = 0;
     	for(std::vector<node*>::iterator itn=node9.begin(); itn<node9.end(); ++itn){
	    //initialize
	    (*itn)->initPos42(); (*itn)->initPos92();
	    
	    if( (*itn)->getEnrich_node92() ){
		(*itn)->setPos92(nsnode92+1);
		nsnode92 = nsnode92+1;
		numsnode92 = numsnode92+1;
	    }
	}

	//int nsnode42 = 0;
	numsnode42 = 0;
	int nsnode42 = 0;
 	for(std::vector<node*>::iterator itn=node4.begin(); itn<node4.end(); ++itn){
	    if( (*itn)->getEnrich_node42() ){
		(*itn)->setPos42(nsnode42+1);
		nsnode42 = nsnode42+1;
		numsnode42 = numsnode42+1;
	    }
	}


	//--------------------------------------------//
	// get level-set/elements intersection points //
	//--------------------------------------------//

	for(std::vector<ele2dcoupled*>::iterator it2d=element.begin(); it2d<element.end(); ++it){
	    // inialize xc2, yc2, xcp2, ycp2
	    (*it2d)->initXc2Yc2Xcp2Ycp2();
	    
	    if( (*it2d)->getSplit_elem2() ){
		(*it2d)->intersectls9();
	    }
	}

	sort_points();

	// node = Ordered_point;
	nodex.clear(); nodey.clear();
    	std::list<ele2dcoupled*>::const_iterator it_or;
    	for(it_or=Split_ordered.begin(); it_or!=Split_ordered.end(); ++it_or){
  	    matrix node_mat = (*it_or)->getOrdered_point();
	    nodex.push_back(node_mat(1,1));	// not include the x for the last point
	    nodey.push_back(node_mat(1,2));	// not include the y for the last point
   	}
    	it_or--;	// the last split element
    	matrix node_mat = (*it_or)->getOrdered_point();
    	nodex.push_back(node_mat(2,1));
    	nodey.push_back(node_mat(2,2));
    	it_or = Split_ordered.begin();
    	node_mat = (*it_or)->getOrdered_point();
    	nodex.push_back(node_mat(1,1));
    	nodey.push_back(node_mat(1,2));

	for(std::vector<node*>::iterator itn=node9.begin(); itn<node9.end(); ++itn){
	    if(isInSplit_ordered(*itn) == false){
		REAL inpoly_val = inpoly((*itn)->getX(), (*itn)->getY(),
					 nodex.begin(), nodex.end(), nodey.begin(), nodey.end());
		REAL temp = sign(-1.0*inpoly_val+0.5);
		(*itn)->setPhi92(temp);
	    }
	}

	for(std::list<ele2dcoupled*>::iterator it2d=Split_ordered.begin(); it2d!=Split_ordered.end(); ++it2d){
	    (*it2d)->calcPhi92();
	}
	
	// cannot be combined with above for loop
	for(std::vector<ele2dcoupled*>::iterator it2d=element.begin(); it2d<element.end(); ++it2d){
	    (*it2d)->setPhi42EqualPhi92();

	    // initialize t_element_total, dtdx_total and t_element
	    (*it2d)->initTEleTotalDtdx();
	}

	for(std::list<ele2dcoupled*>::iterator it2d=Split_ordered.begin(); it2d!=Split_ordered.end(); ++it2d){
	    for(int kn=1; kn<1+9; ++kn){
		matrix pt(1,2);
 		pt(1,1) = (*it2d)->getNode9x(kn-1);
		pt(1,2) = (*it2d)->getNode9y(kn-1);

		matrix pt_proj, t_element_pt, dtdx_pt;

	 	// find_proj_normal_deriv_pt.m is the same as find_proj_normal_deriv_pt_parallel.m 
		find_proj_normal_deriv_pt_parallel(pt_proj, t_element_pt, dtdx_pt, pt);

	    	(*it2d)->calcTEleTotalDtdx(t_element_pt, dtdx_pt, kn);

	    }
	}

    } // end of while


} // end of solve()


// construct mesh 
void surfacetension::constructMesh(){
    // Number of nodes along two directions
    int nnx9 = element_x+1;
    int nny9 = element_y+1;

    int nnx4 = nnx9 - (nnx9-1)*0.5;
    int nny4 = nny9 - (nny9-1)*0.5;

    // Four corner points
    pt1[0] = -Lx*0.5;
    pt1[1] = -Ly*0.5;

    pt2[0] = Lx*0.5;
    pt2[1] = -Ly*0.5;

    pt3[0] = Lx*0.5;
    pt3[1] = Ly*0.5;

    pt4[0] = -Lx*0.5;
    pt4[1] = Ly*0.5;

    // 9-nodes element mesh for velocity
    meshRectangularRegion(nnx9, nny9);

    dx9 = node9[1]->getX()-node9[0]->getX();	// node9(2,1)-node9(1,1);
    dy9 = node9[nnx9]->getY()-node9[0]->getY();	// node9(1+nnx9,2)-node9(1,2);
    REAL dx4 = 2*dx9;
    REAL dy4 = dx4;

    // compute number of nodes, of elements
    numnode9 = node9.size();
    numelem9 = element.size();
    numnode4 = node4.size();
    numelem4 = numelem9;	// the same number of 4-node element and 9-node element

    //std::vector<node*> numnode_S;	// not used in matlab code
//    for(std::vector<node*>::const_iterator it=node4.begin(); it<node4.end(); ++it){
//	if( (*it)->getX()+0.645<=0 && fabs((*it)->getY()-c2[1])<0.11 )
//	    numnode_S.push_back(*it);
//    }

    //std::vector<node*> numnode_v;	// not used in matlab code
//    for(std::vector<node*>::const_iterator it=node9.begin(); it<node9.end(); ++it){
//	if( (*it)->getX()+0.645<=0 && fabs((*it)->getY()-c2[1])<0.11 )
//	    numnode_v.push_back(i);
//    }

} // end of constructMesh()


// level set
void surfacetension::levelSet(){

    //--------------------------------------------------//
    //compute signed distance to vesicle value at nodes //
    //--------------------------------------------------//
    int node_size = node9.size();
    int element_size = element.size();

#pragma omp parallel
  {

    int num_threads = omp_get_num_threads();
    int thread_num  = omp_get_thread_num();

    int start_num = BLOCK_LOW(thread_num, num_threads, node_size);
    int end_num   = BLOCK_HIGH(thread_num, num_threads, node_size);

    for(int itn=start_num; itn<end_num+1; ++itn){
	node9[itn]->setPhi91(1);	// Phi91 = ones(1,numnode9);
	node9[itn]->setPhi41(1);	// Phi41 = ones(1,numnode4);
	node9[itn]->setPhi92(0);	// Phi92 = zeros(1,numelem9);
	node9[itn]->setPhi42(0);	// Phi42 = zeros(1,numnode4);
    }

    //----------------------------------//
    //compute level set2 value at nodes //
    //----------------------------------//

    //Phi92s = zeros(1,numelem9); // not used in matlab code   

    //ellipse
//    for(int n=1; n!=numnode9+1; ++n){
//	Phi92s(1,n) = ( (1.0/a*(node9(n,1)-c2[0]))*(1.0/a*(node9(n,1)-c2[0]))
//		      (1.0/b*(node9(n,2)-c2[1]))*(1.0/b*(node9(n,2)-c2[1])) - 1 );
//    }
//
//    for(int n=1; n!=numnode9+1; ++n){
//	Phi92(1,n) = sign(Phi92s(1,n))*distancePointToEllipse(node9(n,1),node9(n,2),a,b,c2[0],c2[1]);
//    }

    // circle
    for(int itn=start_num; itn<end_num+1; ++itn){
	REAL Phi92 = pow( (node9[itn]->getX()-c2[0])*(node9[itn]->getX()-c2[0]) 
                         +(node9[itn]->getY()-c2[1])*(node9[itn]->getY()-c2[1]), 0.5) - R2;

	node9[itn]->setPhi92(Phi92);
    }

    start_num = BLOCK_LOW(thread_num, num_threads, element_size);
    end_num   = BLOCK_HIGH(thread_num, num_threads, element_size);

    for(int it2d=start_num; it2d<end_num+1; ++it2d){
	element[it2d]->setPhi42EqualPhi92();	// Phi42(sctr4) = Phi92(sctr9(1:4));
	element[it2d]->setCenter_node();	// center_node(iel,:) = node9(sctr9(9),:);
    }

  } // end omp parallel

} // end of levelSet()


// boundaryNodes()
void surfacetension::boundaryNodes(){

    // initialize
    topNodes.clear(); botNodes.clear(); 
    leftNodes.clear(); rightNodes.clear();
    for(std::vector<node*>::const_iterator itn=node9.begin(); itn<node9.end(); ++itn){

	if(fabs((*itn)->getY()-0.5*Ly) < 1e-10){
	    (*itn)->setTopNodes(true);
	    topNodes.push_back(*itn);	// boundary nodes vector has no numbering
	}
	if(fabs((*itn)->getY()+0.5*Ly) < 1e-10){
	    (*itn)->setBotNodes(true);
	    botNodes.push_back(*itn);
	}
	if(fabs((*itn)->getX()+0.5*Lx) < 1e-10){
	    (*itn)->setLeftNodes(true);
	    leftNodes.push_back(*itn);
	}
	if(fabs((*itn)->getX()-0.5*Lx) < 1e-10){
	    (*itn)->setRightNodes(true);
	    rightNodes.push_back(*itn);
	}
    }

/*
    // initialize
    topNodes4.clear(); botNodes4.clear(); 
    leftNodes4.clear(); rightNodes4.clear();
    middletop4.clear(); middlebot4.clear();
    middleleft4.clear(); middleright4.clear(); 
    for(std::vector<node*>::const_iterator it=node4.begin(); it<node4.begin(); ++it){

	if(fabs((*it)->getY()-0.5*Ly) < 1e-10)
	    topNodes4.push_back((*it)->getNum());
	if(fabs((*it)->getY()+0.5*Ly) < 1e-10)
	    botNodes4.push_back((*it)->getNum());
	if(fabs((*it)->getX()+0.5*Lx) < 1e-10)
	    leftNodes4.push_back((*it)->getNum());
	if(fabs((*it)->getX()-0.5*Lx) < 1e-10)
	    rightNodes4.push_back((*it)->getNum());

	if(fabs((*it)->getX()) < 1.1*element_length && fabs((*it)->getY()-0.5*Ly) < 1e-10)
	    middletop4.push_back((*it)->getNum());
     	if(fabs((*it)->getX()) < 1.1*element_length && fabs((*it)->getY()+0.5*Ly) < 1e-10)
            middlebot4.push_back((*it)->getNum());
        if(fabs((*it)->getY()) < 1.1*element_length && fabs((*it)->getX()-0.5*Lx) < 1e-10)
            middleleft4.push_back((*it)->getNum());
        if(fabs((*it)->getY()) < 1.1*element_length && fabs((*it)->getX()+0.5*Lx) < 1e-10)
            middleright4.push_back((*it)->getNum());

    }
*/

    markers.clear();	// initialize
    for(std::vector<node*>::const_iterator itn=leftNodes.begin(); itn<leftNodes.end(); ++itn){
	if((*itn)->getY() > 0)
	    markers.push_back(*itn);	// in matlab, markers are coords for those leftnodes whose y>0
    }


    //--------------------//
    // get enriched nodes //
    //--------------------//
    split_elem2.clear();

    int node_size = node9.size();
#pragma omp parallel
  {
    int num_threads = omp_get_num_threads();
    int thread_num  = omp_get_thread_num();

    int start_num = BLOCK_LOW(thread_num, num_threads, node_size);
    int end_num   = BLOCK_HIGH(thread_num, num_threads, node_size);

    for(int itn=start_num; itn<end_num+1; ++itn){
	node9[itn]->setEnrich_node92(false);	// enrich_node92 = zeros(numnode9,1);
	node9[itn]->setEnrich_node42(false);	// enrich_node42 = zeros(numnode4,1); node9 > node4
    }
  } // end omp parallel

    for(std::vector<ele2dcoupled*>::iterator it2d=element.begin(); it2d<element.end(); ++it2d){
	
	(*it2d)->setSplit_elem2(false);	// initialize
	(*it2d)->setEnrich_element2(false);
	(*it2d)->setElement_length(0);
 	(*it2d)->initCenter_node();

	matrix phi2 = (*it2d)->getPhi42();
	(*it2d)->setCenter_node();	// in matlab, center_node is coords of each element9
	// no need to update center_node since, in ele2dcoupled, node* center_node;
	(*it2d)->calcElement_length();

	if( max(phi2)*min(phi2) < 0 ){
	    //count1 = count1+1;
	    (*it2d)->setSplit_elem2(true);
	    (*it2d)->setEnrich_node42(true);
	    (*it2d)->setEnrich_node92(true);
	    (*it2d)->setEnrich_element2(true);
	    split_elem2.push_back(*it2d);	// this function cannot be parallelized
	}
    }

    element_length = element[0]->getElement_length();	// element_length(1)
    
    numsnode92 = 0;
    int nsnode92 = 0;
    for(std::vector<node*>::iterator itn=node9.begin(); itn<node9.end(); ++itn){
	//initialize
	(*itn)->initPos42(); (*itn)->initPos92();
	    
	if( (*itn)->getEnrich_node92() ){
	    
 	    (*itn)->setPos92(nsnode92+1);
	    nsnode92++;
	    numsnode92++;
	}
    }

    numsnode42 = 0;
    int nsnode42 = 0;
    for(std::vector<node*>::iterator itn=node4.begin(); itn<node4.end(); ++itn){
 	if( (*itn)->getEnrich_node42() ){
		
	    (*itn)->setPos42(nsnode42+1);
	    nsnode42++;
	    numsnode42++;
	}
    }

#pragma omp parallel
  {
    int num_threads = omp_get_num_threads();
    int thread_num  = omp_get_thread_num();

    int start_num = BLOCK_LOW(thread_num, num_threads, node_size);
    int end_num   = BLOCK_HIGH(thread_num, num_threads, node_size);

    for(int itn=start_num; itn<end_num+1; ++itn){
	
	node9[itn]->setEnrich_node91(false);	// enrich_node91 = zeros(numnode9,1);
	node9[itn]->setEnrich_node41(false);	// enrich_node41 = zeros(numnode4,1);

    }
  } // end omp parallel


} // end of boundaryNodes()


// initialize grid based particle
void surfacetension::initGridParticle(){

/*    REAL element_length_min = element[0]->getElement_length();
    for(std::vector<ele2dcoupled*>::iterator it2d=element.begin(); it2d<element.end(); ++it2d){
	if(element_length_min > (*it2d)->getElement_length())
	    element_length_min = (*it2d)->getElement_length();
    }
*/

    int element_size = element.size();
#pragma omp parallel
  {

    int num_threads = omp_get_num_threads();
    int thread_num  = omp_get_thread_num();

    int start_num = BLOCK_LOW(thread_num, num_threads, element_size);
    int end_num   = BLOCK_HIGH(thread_num, num_threads, element_size);

    for(int it2d=start_num; it2d<end_num+1; ++it2d){

	element[it2d]->setW_update(false);	// initialize
	matrix phi;
	phi = element[it2d]->getPhi92();
	matrix phi_abs = abs(phi);
	if(min(phi_abs)<=1.1*element_length){
	    element[it2d]->setW_update(true);
     	}
    }
  } // end omp parallel
  
    sort2D_xfm();	

    connect2D = node9;

/*    REAL element_length_min_temp;
    for(std::vector<ele2dcoupled*>::const_iterator it2d=element.begin(); it2d<element.end(); ++it2d){
	if(it2d==element.begin()){
	    element_length_min_temp = (*it2d)->getElement_length();
	}
	else if(element_length_min_temp>(*it2d)->getElement_length())
	    element_length_min_temp = (*it2d)->getElement_length();

    }
//element length are the same, so element_length_min_temp = elment_length(1)
*/ 
    int node_size = node9.size();
#pragma omp parallel
  {
    int num_threads = omp_get_num_threads();
    int thread_num  = omp_get_thread_num();

    int start_num = BLOCK_LOW(thread_num, num_threads, node_size);
    int end_num   = BLOCK_HIGH(thread_num, num_threads, node_size);

    for(int itn=start_num; itn<end_num+1; ++itn){
	if(isInW_update(node9[itn])){
	    node9[itn]->normalR2(connect2D.begin(), connect2D.end(), Ly, Lx, element_length);
	}
    }
  } // end omp parallel

    element_length = element[0]->getElement_length();	// element_length(1)

    Init_Foot_points_func();	// create Active_points

    for(std::list<activepoint*>::iterator itpt=Active_points.begin(); itpt!=Active_points.end(); ++itpt){
	
	(*itpt)->initAllFoots();

	matrix temp_mat(1,2);
	temp_mat(1,1) = 0; temp_mat(1,2) = 0;

	(*itpt)->setVel_Foot_points_membrane(temp_mat);
	(*itpt)->setE_Foot_points(temp_mat);
	(*itpt)->setT_Foot_points(temp_mat);
	(*itpt)->setJ_Foot_points(1);

	temp_mat(1,1) = 1; temp_mat(1,2) = 1;
	(*itpt)->setF_Foot_points(temp_mat);

    }

    Foot_point_ressampling_v2();
    find_laplace_beltrami_coeff();

    for(std::list<activepoint*>::iterator itpt=Active_points.begin(); itpt!=Active_points.end(); ++itpt){
	
	(*itpt)->initAllFoots();

	matrix temp_mat(1,2);
	temp_mat(1,1) = 0; temp_mat(1,2) = 0;

	(*itpt)->setVel_membrane(temp_mat);
	(*itpt)->setE_Foot_points(temp_mat);
	(*itpt)->setT_Foot_points(temp_mat);
	(*itpt)->setJ_Foot_points(1);

	temp_mat(1,1) = 1; temp_mat(1,2) = 1;
	(*itpt)->setF_Foot_points(temp_mat);

    }

    E_F_J_T_coeff_updating();

} // end of initGridParticle()


// initLevelSet()
void surfacetension::initLevelSet(){

    int element_size = element.size();
#pragma omp parallel
  {
    int num_threads = omp_get_num_threads();
    int thread_num  = omp_get_thread_num();

    int start_num = BLOCK_LOW(thread_num, num_threads, element_size);
    int end_num   = BLOCK_HIGH(thread_num, num_threads, element_size);

    for(int it2d=start_num; it2d<end_num+1; ++it2d){
	element[it2d]->initXc2Yc2Xcp2Ycp2();
	element[it2d]->initLe2Le2_old();

	if( element[it2d]->getSplit_elem2() ){
	    element[it2d]->intersectls9();
	}
    }
  } // end omp parallel

    sort_points();

    // node = Ordered_point;
    nodex.clear(); nodey.clear();
    std::list<ele2dcoupled*>::const_iterator it_or;
    for(it_or=Split_ordered.begin(); it_or!=Split_ordered.end(); ++it_or){
  	matrix node_mat = (*it_or)->getOrdered_point();
	nodex.push_back(node_mat(1,1));	// not include the x for the last point
	nodey.push_back(node_mat(1,2));	// not include the y for the last point
    }
    it_or--;	// the last split element
    matrix node_mat = (*it_or)->getOrdered_point();
    nodex.push_back(node_mat(2,1));
    nodey.push_back(node_mat(2,2));
    it_or = Split_ordered.begin();
    node_mat = (*it_or)->getOrdered_point();
    nodex.push_back(node_mat(1,1));
    nodey.push_back(node_mat(1,2));

    if(plain_axi=="axi1") {
	    A0=0;
	    for(std::list<ele2dcoupled*>::const_iterator it2d=Split_ordered.begin(); it2d!=Split_ordered.end(); ++it2d){
		matrix node_mat = (*it2d)->getOrdered_point();
		REAL Rzm = (node_mat(1,1)+0.5*Lx+node_mat(2,1)+0.5*Lx)*0.5;
		REAL dz = node_mat(2,2)-node_mat(1,2);
		A0 = A0 + PI*(Rzm*Rzm*dz);
	    }
    } // if "axi1"

    if(plain_axi=="axi2"){
	    A0=0;
	    for(std::list<ele2dcoupled*>::const_iterator it2d=Split_ordered.begin(); it2d!=Split_ordered.end(); ++it2d){
		matrix node_mat = (*it2d)->getOrdered_point();
		REAL Rzm = (node_mat(1,1)+0.5*Lx+node_mat(2,1)+0.5*Lx)*0.5;
		REAL dz = node_mat(2,2)-node_mat(1,2);
		A0 = A0 + PI*(Rzm*Rzm*dz);
	    }
    } // if "axi2"

    if(plain_axi=="plain"){
	    nodex.clear(); nodey.clear();
	    std::list<ele2dcoupled*>::const_iterator it2d;
	    for(it2d=Split_ordered.begin(); it2d!=Split_ordered.end(); ++it2d){
		matrix node_mat = (*it2d)->getOrdered_point();
		nodex.push_back(node_mat(1,1));	// not include the x for the last point
		nodey.push_back(node_mat(1,2));	// not include the y for the last point
	    }
	    it2d--;	// the last split element
	    matrix node_mat = (*it2d)->getOrdered_point();
	    nodex.push_back(node_mat(2,1));
	    nodey.push_back(node_mat(2,2));

	    it2d = Split_ordered.begin();
	    node_mat = (*it2d)->getOrdered_point();
	    nodex.push_back(node_mat(1,1));
	    nodey.push_back(node_mat(1,2));	    

	    A0 = polyarea(nodex.begin(), nodex.end(), nodey.begin(), nodey.end());

    } // if "plain"

    int node_size = node9.size();
#pragma omp parallel
  {
    int num_threads = omp_get_num_threads();
    int thread_num  = omp_get_thread_num();

    int start_num = BLOCK_LOW(thread_num, num_threads, node_size);
    int end_num   = BLOCK_HIGH(thread_num, num_threads, node_size);

    for(int itn=start_num; itn<end_num+1; ++itn){
	if(isInSplit_ordered(node9[itn]) == false){
	    REAL temp = sign(-inpoly(node9[itn]->getX(), node9[itn]->getY(),
					 nodex.begin(), nodex.end(), nodey.begin(), nodey.end())+0.5);
	    node9[itn]->setPhi92(temp);
	}
    }

  } // end omp parallel

    for(std::list<ele2dcoupled*>::iterator it=Split_ordered.begin(); it!=Split_ordered.end(); ++it){
	(*it)->calcPhi92();
    }

    // cannot be combined with above for loop
    element_size = element.size();
#pragma omp parallel
  {
    int num_threads = omp_get_num_threads();
    int thread_num  = omp_get_thread_num();

    int start_num = BLOCK_LOW(thread_num, num_threads, element_size);
    int end_num   = BLOCK_HIGH(thread_num, num_threads, element_size);

    for(int it2d=start_num; it2d<end_num+1; ++it2d){
     #pragma omp critical
	element[it2d]->setPhi42EqualPhi92();

	// initialize t_element_total, dtdx_total and t_element
	element[it2d]->initTEleTotalDtdx();
    }
  } // end omp parallel

    for(std::list<ele2dcoupled*>::iterator it2d=Split_ordered.begin(); it2d!=Split_ordered.end(); ++it2d){

    	for(int kn=1; kn<1+9; ++kn){
	    matrix pt(1,2);
 	    pt(1,1) = (*it2d)->getNode9x(kn-1);
	    pt(1,2) = (*it2d)->getNode9y(kn-1);

	    matrix pt_proj, t_element_pt, dtdx_pt;

	    find_proj_normal_deriv_pt_parallel(pt_proj, t_element_pt, dtdx_pt, pt);

	    (*it2d)->calcTEleTotalDtdx(t_element_pt, dtdx_pt, kn);
    	}	
    }

} // end of initLevelSet()


// gauss points and weights
void surfacetension::gaussPoints(const int it){


int test_num2=0;
    for(std::vector<ele2dcoupled*>::iterator it2d=element.begin(); it2d<element.end(); ++it2d){
	if( (*it2d)->getSplit_elem2() ){
	    // find the gauss points coord and weight in 2D split parent element9
	    (*it2d)->discontQ4quad();	//(order) is always 3
	}
	else{
	    // find the gauss points coord and weight in 2D regular parent element9
//	    int order = 3;
	    (*it2d)->quadrature();	//(order, "GAUSS", 2); order is always be 3, always "GAUSS", always 2
	}
    }
    // Qx Qy have already been changed

    for(std::vector<ele2dcoupled*>::iterator it2d=element.begin(); it2d<element.end(); ++it2d){

	(*it2d)->initNogp2Wl2Ql2Paired_connect2();
	(*it2d)->setLe2_oldEqualLe2();
	(*it2d)->setLe2(0);

	if( (*it2d)->getSplit_elem2() ){ // split element
	   
	    (*it2d)->calcLe2();
	    if(Moes_Assembly2.num_row!=0 || Moes_Assembly2.num_col!=0){ 
		// ???????????? Moes_Assembly2 is always empty
		//if(length())
	    }
	}
    }
    
    if(it==2){
	for(std::vector<ele2dcoupled*>::iterator it2d=element.begin(); it2d<element.end(); ++it2d){
	    (*it2d)->setLe2_oldEqualLe2();
	}
    }

    for(std::list<ele2dcoupled*>::iterator it2d=Split_ordered.begin(); it2d!=Split_ordered.end(); ++it2d){

	(*it2d)->calcNogp2Wl2Ql_2D(dx9, dy9);

    }
    // Ql_2Dx, Ql_2Dy have already been changed

}// end of gaussPoints()


//applyBdryConditions()
void surfacetension::applyBdryConditions(){

    std::vector<int> bcdofy, bcdofx;
    for(std::vector<node*>::const_iterator itbot=botNodes.begin(); itbot<botNodes.end(); ++itbot){
		
       	bcdofy.push_back((*itbot)->getNum_node9()*2);
    }
    for(std::vector<node*>::const_iterator itbot=leftNodes.begin(); itbot<leftNodes.end(); ++itbot){
		
	bcdofx.push_back((*itbot)->getNum_node9()*2-1);
    }
	
    for(std::vector<int>::const_iterator itx=bcdofx.begin(); itx<bcdofx.end(); ++itx){
#pragma omp parallel for
	for(int i=1; i<total_unknown+1; ++i){	// K is total_unknown x total_unknown
	    K((*itx), i) = 0;
	    K(i, (*itx)) = 0;
	}
    }
    for(std::vector<int>::const_iterator ity=bcdofy.begin(); ity!=bcdofy.end(); ++ity){
#pragma omp parallel for
    	for(int i=1; i<total_unknown+1; ++i){	// K is total_unknown x total_unknown
	    K((*ity), i) = 0;
	    K(i, (*ity)) = 0;
	}
    }

    // K(bcdofx,bcdofx) = speye(nbcx)
    int ii=1;
    for(std::vector<int>::const_iterator iti=bcdofx.begin();
	iti<bcdofx.end(); ++iti, ++ii){
	int jj=1;
	for(std::vector<int>::const_iterator itj=bcdofx.begin();
	    itj<bcdofx.end(); ++itj, ++jj){

	    if(ii==jj)
	 	K((*iti), (*itj)) = 1;
	    else
		K((*iti), (*itj)) = 0;
	}
    }
    ii=1;
    for(std::vector<int>::const_iterator iti=bcdofy.begin();
	iti<bcdofy.end(); ++iti, ++ii){
	int jj=1;
	for(std::vector<int>::const_iterator itj=bcdofy.begin();
	    itj<bcdofy.end(); ++itj, ++jj){

	    if(ii==jj)
	 	K((*iti), (*itj)) = 1;
            else
		K((*iti), (*itj)) = 0;
	}
    }
	
    // F(bcdofx,1)=0
    for(std::vector<int>::const_iterator itx=bcdofx.begin(); itx<bcdofx.end(); ++itx){
	F((*itx), 1) = 0;
    }
    for(std::vector<int>::const_iterator ity=bcdofy.begin(); ity!=bcdofy.end(); ++ity){
	F((*ity), 1) = 0;
    }

} // end of applyBdryConditions()


bool surfacetension::isInSplit_ordered(node* in) const{

    for(std::list<ele2dcoupled*>::const_iterator it=Split_ordered.begin(); it!=Split_ordered.end(); ++it){
	if( (*it)->isNode9In(in) )
	    return true;
    }

    return false;

} // isInSplit_ordered()


bool surfacetension::isInW_update(node* in) const{

    for(std::vector<ele2dcoupled*>::const_iterator it2d=element.begin(); it2d<element.end(); ++it2d){
	if( (*it2d)-> getW_update() )
	    if( (*it2d)->isNode9In(in) )
	    	return true;
    }

    return false;

} // isInW_update()


void surfacetension::meshRectangularRegion(int nnx9, int nny9){
// obtain number of elements from number of nodes

    int numx9 = (nnx9-1)*0.5;
    int numy9 = (nny9-1)*0.5;
    square_node_array(nnx9, nny9);
    int inc_u9 = 2;
    int inc_v9 = 2*nnx9;
    matrix node_pattern9(1,9);
    node_pattern9(1,1) = 1; node_pattern9(1,2) = 3; 
    node_pattern9(1,3) = 2*nnx9+3; node_pattern9(1,4) = 2*nnx9+1;
    node_pattern9(1,5) = 2; node_pattern9(1,6) = nnx9+3; node_pattern9(1,7) = 2*nnx9+2;
    node_pattern9(1,8) = nnx9+1; node_pattern9(1,9) = nnx9+2;
    make_elem(node_pattern9, numx9, numy9, inc_u9, inc_v9);

} // meshRectangularRegion()


void surfacetension::square_node_array(int numnode_u, int numnode_v){
//create node9 & node4

    node9.clear(); node4.clear();

    matrix xi_pts = linspace(-1.0, 1.0, numnode_u);
    matrix eta_pts = linspace(-1.0, 1.0, numnode_v);

    matrix x_pts(1,4);
    matrix y_pts(1,4);
    x_pts(1,1) = pt1[0]; x_pts(1,2) = pt2[0]; x_pts(1,3) = pt3[0]; x_pts(1,4) = pt4[0];
    y_pts(1,1) = pt1[1]; y_pts(1,2) = pt2[1]; y_pts(1,3) = pt3[1]; y_pts(1,4) = pt4[1];
    int num_node9 = 1;
    int num_node4 = 1;
    for(int r=1; r<numnode_v+1; ++r){
	REAL eta = eta_pts(1,r);
	for(int c=1; c<numnode_u+1; ++c){
	    REAL xi = xi_pts(1,c);
	    matrix N(4,1);
	    N(1,1) = 0.25*(1-xi)*(1-eta);
	    N(2,1) = 0.25*(1+xi)*(1-eta);
	    N(3,1) = 0.25*(1+xi)*(1+eta);
	    N(4,1) = 0.25*(1-xi)*(1+eta);
	
	    REAL coords_x, coords_y;
	    coords_x = (x_pts*N)(1,1);
	    coords_y = (y_pts*N)(1,1);

	    node* pt = new node(coords_x, coords_y);
	    pt->setNum_node9(num_node9);
	    node9.push_back(pt);	// this function cannot be parallelized
	    num_node9++;
	
	    if(r%2!=0 && c%2!=0){	// odd row and col is node4
		pt->setNum_node4(num_node4);
		node4.push_back(pt);
		num_node4++;    
	    }
	
	}
    }

} // square_node_array


void surfacetension::make_elem(matrix &node_pattern, int num_u, int num_v, int inc_u, int inc_v){
//create element

    element.clear();

    int inc = 0;;

    std::vector<node*> temp_node;
    matrix num_nodes(1,9);
    int num_ele = 1;
    for(int row=1; row<num_v+1; ++row){
	for(int col=1; col<num_u+1; ++col){
	    num_nodes = node_pattern+inc;	// the numbers of node for element9
	    temp_node.clear();
	    for(int i=1; i<node_pattern.num_col+1; ++i){
		int num_node = num_nodes(1,i);
		temp_node.push_back(node9[num_node-1]);	// num_node starts from 1
							// this function cannot be parallelized
	    }
	    ele2dcoupled* pt = new ele2dcoupled(temp_node.begin(), temp_node.end(), num_ele);
	    element.push_back(pt);
	    num_ele++;
	    inc = inc+inc_u;
	}
	inc = row*inc_v;
    }

} // make_elem


void surfacetension::sort2D_xfm(){

    connect2D = node9;
    std::sort(connect2D.begin(), connect2D.end(), less_than);

} // sort2D_xfm()


void surfacetension::Init_Foot_points_func(){

    for(std::vector<node*>::iterator it=node9.begin(); it<node9.end(); ++it){
	REAL Phi92 = (*it)->getPhi92();	// check if Phi92 in this function is REAL
	if(fabs(Phi92)<=0.6*element_length){
	    activepoint* pt_act = new activepoint(*it);
	    setFoot_points_n_Foot_points_actpnt(pt_act);

	    // add to Active_points
	    Active_points.push_back(pt_act);	// this function cannot be parallelized
   	}
    }

} // Init_Foot_points_func()


void surfacetension::setFoot_points_n_Foot_points_actpnt(activepoint* pt_act){
// calculate Foot_points & n_Foot_points when created as in "Init_Foot_func.m"	

    matrix pt(1,2);
    pt = pt_act->getCoords();


    matrix ne0(2,1);
    ne0(1,1) = pt_act->getGrad_phi_x();
    ne0(2,1) = pt_act->getGrad_phi_y();        
 
    REAL phi0 = pt_act->getPhi92();         

    matrix Foot_points = pt-phi0*ne0.getTrans();
    pt_act->setFoot_points(Foot_points);

    // k0 = dsearchn(center_node,[pt_proj(1) pt_proj(2)]);
    REAL dist_min;
    std::vector<ele2dcoupled*>::const_iterator itn = element.begin();
    for(std::vector<ele2dcoupled*>::const_iterator it2d=element.begin(); it2d<element.end(); ++it2d){
	matrix center_node = ((*it2d)->getCenter_node());
	REAL dist = pow((Foot_points(1,1)-center_node(1,1))*(Foot_points(1,1)-center_node(1,1))
		      + (Foot_points(1,2)-center_node(1,2))*(Foot_points(1,2)-center_node(1,2)), 0.5);

	if(it2d==element.begin()){
	    dist_min = dist;
	    itn = it2d;
	}
	else{
	    if(dist<dist_min){
		dist_min = dist;
		itn = it2d;
	    }
	}

    }

//    int k0 = (*itn)->getEle_num();
    matrix center_node = (*itn)->getCenter_node();

    REAL xi, eta;
    xi = (Foot_points(1,1)-center_node(1,1))/dx9;
    eta = (Foot_points(1,2)-center_node(1,2))/dy9;
        
    // [N,~]=lagrange_basis('Q9',pt_par);
    matrix N(9,1);
    N(1,1) = 0.25*xi*eta*(xi-1)*(eta-1);
    N(2,1) = 0.25*xi*eta*(xi+1)*(eta-1);
    N(3,1) = 0.25*xi*eta*(xi+1)*(eta+1);
    N(4,1) = 0.25*xi*eta*(xi-1)*(eta+1);
    N(5,1) = -0.25*2*eta*(xi+1)*(xi-1)*(eta-1);
    N(6,1) = -0.25*2*xi*(xi+1)*(eta+1)*(eta-1);
    N(7,1) = -0.25*2*eta*(xi+1)*(xi-1)*(eta+1);
    N(8,1) = -0.25*2*xi*(xi-1)*(eta+1)*(eta-1);
    N(9,1) = 0.25*4*(xi+1)*(xi-1)*(eta+1)*(eta-1);  
    
//    int BCjump=1;
        
//    matrix N1 = xfem_Nmu2(1,BCjump,N);
        
    matrix Nfem = N;
        
    matrix n_Foot_points(1,2);
    n_Foot_points(1,1) = ( Nfem.getTrans()*(*itn)->getGrad_phi_x() )(1,1) ; 
    n_Foot_points(1,2) = ( Nfem.getTrans()*(*itn)->getGrad_phi_y() )(1,1); 

    pt_act->setN_Foot_points(n_Foot_points);

} // setFoot_points_n_Foot_points_actpnt()



//matrix surfacetension::xfem_Nmu2(int iel, int BCjump, matrix & N){
//
//    int nn = 9;
//    matrix phi1 = element[iel-1]->getPhi91();	// vector starts from 0
//    matrix phi2 = element[iel-1]->getPhi92();
//
//    matrix Nfem = N;
//
//} // xfem_Nmu2()


void surfacetension::Foot_point_ressampling_v2(){

    matrix Foot_points, n_Foot_points, vel_Foot_points_membrane, E_Foot_points, F_Foot_points;
    matrix J_Foot_points(Active_points.size(), 1);

    int ii=1;
    for(std::list<activepoint*>::const_iterator it=Active_points.begin(); it!=Active_points.end(); ++it, ++ii){
	Foot_points.appendRow( ((*it)->getFoot_points()).getRow(1) );	// this cannot be parallelized
	n_Foot_points.appendRow( ((*it)->getN_Foot_points()).getRow(1) );
	vel_Foot_points_membrane.appendRow( ((*it)->getVel_Foot_points_membrane()).getRow(1) );
 	E_Foot_points.appendRow( ((*it)->getE_Foot_points()).getRow(1) );
	F_Foot_points.appendRow( ((*it)->getF_Foot_points()).getRow(1) );
	J_Foot_points(ii,1) = (*it)->getJ_Foot_points();
    }

    if(plain_axi=="axi1"){

	    matrix temp_Foot_points = Foot_points;
	    Foot_points = zeros(2*temp_Foot_points.num_row,2);
#pragma omp parallel for
	    for(int i=1; i<temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = temp_Foot_points(i,1);
		Foot_points(i,2) = temp_Foot_points(i,2);
	    }
#pragma omp parallel for
	    for(int i=1+temp_Foot_points.num_row; i<2*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row,1))-0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-temp_Foot_points.num_row,2);
	    }

	    // n_Foot_points
	    matrix n_Foot_points_col1;
	    n_Foot_points_col1.appendCol(n_Foot_points.getCol(1));
	    n_Foot_points_col1 = -1*n_Foot_points_col1;
	    matrix sym_points2;
	    sym_points2.appendCol(n_Foot_points_col1.getCol(1));
 	    sym_points2.appendCol(n_Foot_points.getCol(2));
	    for(int i=1; i<sym_points2.num_row+1; ++i){
		n_Foot_points.appendRow(sym_points2.getRow(i));
	    }

	    // vel_Foot_points_membrane
	    matrix vel_Foot_points_col1;
	    vel_Foot_points_col1.appendCol(vel_Foot_points_membrane.getCol(1));
	    vel_Foot_points_col1 = -1*vel_Foot_points_col1;
	    matrix sym_points3;
	    sym_points3.appendCol(vel_Foot_points_col1.getCol(1));
 	    sym_points3.appendCol(vel_Foot_points_membrane.getCol(2));
	    for(int i=1; i<sym_points3.num_row+1; ++i){
		vel_Foot_points_membrane.appendRow(sym_points3.getRow(i));
	    }

	    // E_Foot_points
	    matrix temp_mat = E_Foot_points;
	    for(int i=1; i<temp_mat.num_row+1; ++i){
		E_Foot_points.appendRow(temp_mat.getRow(i));
	    }

	    // F_Foot_points
	    matrix temp_mat2 = F_Foot_points;
	    for(int i=1; i<temp_mat2.num_row+1; ++i){
		F_Foot_points.appendRow(temp_mat2.getRow(i));
	    }

	    // J_Foot_points
	    matrix temp_mat3 = J_Foot_points;
	    for(int i=1; i<temp_mat3.num_row+1; ++i){
		J_Foot_points.appendRow(temp_mat3.getRow(i));
	    }
	    
    } // if "axi1"

    if(plain_axi=="axi2"){
	    matrix temp_Foot_points = Foot_points;
	    Foot_points = zeros(3*temp_Foot_points.num_row,2);
	    for(int i=1; i<temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = temp_Foot_points(i,1);
		Foot_points(i,2) = temp_Foot_points(i,2);
	    }
	    for(int i=1+temp_Foot_points.num_row; i<2*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row,1))-0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-temp_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_Foot_points.num_row; i<3*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(-Lx*0.5+temp_Foot_points(i-2*temp_Foot_points.num_row,1))+0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-2*temp_Foot_points.num_row,2);
	    }

	    // n_Foot_points
	    matrix temp_n_Foot_points = n_Foot_points;
	    n_Foot_points = zeros(3*temp_n_Foot_points.num_row,2);
	    for(int i=1; i<temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = temp_n_Foot_points(i,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i,2);
	    }
	    for(int i=1+temp_n_Foot_points.num_row; i<2*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = -temp_n_Foot_points(i-temp_n_Foot_points.num_row,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i-temp_n_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_n_Foot_points.num_row; i<3*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = -temp_n_Foot_points(i-2*temp_n_Foot_points.num_row,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i-2*temp_n_Foot_points.num_row,2);
	    }

	    // vel_Foot_points_membrane
	    matrix temp_vel_Foot_points = vel_Foot_points_membrane;
	    vel_Foot_points_membrane = zeros(3*temp_vel_Foot_points.num_row,2);
	    for(int i=1; i<temp_vel_Foot_points.num_row+1; ++i){
		vel_Foot_points_membrane(i,1) = temp_vel_Foot_points(i,1);
		vel_Foot_points_membrane(i,2) = temp_vel_Foot_points(i,2);
	    }
	    for(int i=1+temp_vel_Foot_points.num_row; i<2*temp_vel_Foot_points.num_row+1; ++i){
		vel_Foot_points_membrane(i,1) = -temp_vel_Foot_points(i-temp_vel_Foot_points.num_row,1);
		vel_Foot_points_membrane(i,2) = temp_vel_Foot_points(i-temp_vel_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_vel_Foot_points.num_row; i<3*temp_vel_Foot_points.num_row+1; ++i){
		vel_Foot_points_membrane(i,1) = -temp_vel_Foot_points(i-2*temp_vel_Foot_points.num_row,1);
		vel_Foot_points_membrane(i,2) = temp_vel_Foot_points(i-2*temp_vel_Foot_points.num_row,2);
	    }

	    // E_Foot_points
	    matrix temp_E_Foot_points = E_Foot_points;
	    E_Foot_points = zeros(3*temp_E_Foot_points.num_row,2);
	    for(int i=1; i<temp_E_Foot_points.num_row+1; ++i){
		E_Foot_points(i,1) = temp_E_Foot_points(i,1);
		E_Foot_points(i,2) = temp_E_Foot_points(i,2);
	    }
	    for(int i=1+temp_E_Foot_points.num_row; i<2*temp_E_Foot_points.num_row+1; ++i){
		E_Foot_points(i,1) = temp_E_Foot_points(i-temp_E_Foot_points.num_row,1);
		E_Foot_points(i,2) = temp_E_Foot_points(i-temp_E_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_E_Foot_points.num_row; i<3*temp_E_Foot_points.num_row+1; ++i){
		E_Foot_points(i,1) = temp_E_Foot_points(i-2*temp_E_Foot_points.num_row,1);
		E_Foot_points(i,2) = temp_E_Foot_points(i-2*temp_E_Foot_points.num_row,2);
	    }

	    // F_Foot_points
	    matrix temp_F_Foot_points = F_Foot_points;
	    F_Foot_points = zeros(3*temp_F_Foot_points.num_row,2);
	    for(int i=1; i<temp_F_Foot_points.num_row+1; ++i){
		F_Foot_points(i,1) = temp_F_Foot_points(i,1);
		F_Foot_points(i,2) = temp_F_Foot_points(i,2);
	    }
	    for(int i=1+temp_F_Foot_points.num_row; i<2*temp_F_Foot_points.num_row+1; ++i){
		F_Foot_points(i,1) = temp_F_Foot_points(i-temp_F_Foot_points.num_row,1);
		F_Foot_points(i,2) = temp_F_Foot_points(i-temp_F_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_F_Foot_points.num_row; i<3*temp_F_Foot_points.num_row+1; ++i){
		F_Foot_points(i,1) = temp_F_Foot_points(i-2*temp_F_Foot_points.num_row,1);
		F_Foot_points(i,2) = temp_F_Foot_points(i-2*temp_F_Foot_points.num_row,2);
	    }

	    // J_Foot_points
	    matrix temp_J_Foot_points = J_Foot_points;
	    J_Foot_points = zeros(3*temp_J_Foot_points.num_row,1);
	    for(int i=1; i<temp_J_Foot_points.num_row+1; ++i){
		J_Foot_points(i,1) = temp_J_Foot_points(i,1);
	    }
	    for(int i=1+temp_J_Foot_points.num_row; i<2*temp_J_Foot_points.num_row+1; ++i){
		J_Foot_points(i,1) = temp_J_Foot_points(i-temp_J_Foot_points.num_row,1);
	    }
	    for(int i=1+2*temp_J_Foot_points.num_row; i<3*temp_J_Foot_points.num_row+1; ++i){
		J_Foot_points(i,1) = temp_J_Foot_points(i-2*temp_J_Foot_points.num_row,1);
	    }

    } // if "axi2"

    if(plain_axi=="plain"){
	
    } // if "plain"

    for(std::list<activepoint*>::iterator it=Active_points.begin(); it!=Active_points.end(); ++it){

	(*it)->initAllCoef_AllFoots();

	matrix p = (*it)->getCoords();

	std::list<REAL> D2;
	for(int i=1; i<Foot_points.num_row+1; ++i){
	    REAL D2_val = pow( (Foot_points(i,1)-p(1,1))*(Foot_points(i,1)-p(1,1))
			     + (Foot_points(i,2)-p(1,2))*(Foot_points(i,2)-p(1,2)), 0.5);
	    D2.push_back(D2_val);
	}

	std::vector<int> IDX2 = sort30(D2.begin(), D2.end());

	matrix y_sorted = zeros(30,2);
   	matrix vel_sorted = zeros(30,2);
	matrix E_sorted = zeros(30,2);
	matrix F_sorted = zeros(30,2);
	matrix J_sorted = zeros(30,1);
	matrix n_sorted = zeros(30,2);

	for(int i=0; i<30; ++i){
	    y_sorted(i+1, 1) = Foot_points(IDX2[i], 1);
	    y_sorted(i+1, 2) = Foot_points(IDX2[i], 2);

	    vel_sorted(i+1, 1) = vel_Foot_points_membrane(IDX2[i], 1);
	    vel_sorted(i+1, 2) = vel_Foot_points_membrane(IDX2[i], 2);

	    E_sorted(i+1, 1) = E_Foot_points(IDX2[i], 1);
	    E_sorted(i+1, 2) = E_Foot_points(IDX2[i], 2);

	    F_sorted(i+1, 1) = F_Foot_points(IDX2[i], 1);
	    F_sorted(i+1, 2) = F_Foot_points(IDX2[i], 2);

	    J_sorted(i+1, 1) = J_Foot_points(IDX2[i], 1);

	    n_sorted(i+1, 1) = n_Foot_points(IDX2[i], 1);
	    n_sorted(i+1, 2) = n_Foot_points(IDX2[i], 2);
	}

	//n0_old and y0 are the origin and orientation of the local coordinate system

	matrix n0_old(1,2);
	n0_old(1,1) = n_Foot_points(IDX2[0], 1);
	n0_old(1,2) = n_Foot_points(IDX2[0], 2);

	matrix y0(1,2);
	y0(1,1) = y_sorted(1,1);
	y0(1,2) = y_sorted(1,2);

	matrix v0(1,2);
	v0(1,1) = vel_sorted(1,1);
	v0(1,2) = vel_sorted(1,2);

	matrix E0(1,2);
	E0(1,1) = E_sorted(1,1);
	E0(1,2) = E_sorted(1,2);

	matrix F0(1,2);
	F0(1,1) = F_sorted(1,1);
	F0(1,2) = F_sorted(1,2);

	REAL J0 = J_sorted(1,1);

	matrix n0(1,2);
	n0(1,1) = n_sorted(1,1);
	n0(1,2) = n_sorted(1,2);

	//we now select a few neighboring Foot_points to construct the 
	//polynomials that interpolate the interface

      	matrix y_picked = y0;
    	matrix vel_picked = v0;
    	matrix E_picked = E0;
    	matrix F_picked = F0;
    	matrix J_picked(1,1); J_picked(1,1) = J0;
    	matrix n_picked = n0;

	for(int i=2; i<length(y_sorted)+1; ++i){
	    REAL d;
	    for(int j=1; j<y_picked.num_row+1; ++j){
	    	REAL d_temp = pow( (y_picked(j,1)-y_sorted(i,1))*(y_picked(j,1)-y_sorted(i,1))
			         + (y_picked(j,2)-y_sorted(i,2))*(y_picked(j,2)-y_sorted(i,2)), 0.5);

		if(j==1){
		    d = d_temp;
		}
		else{
		    if(d>d_temp)
			d = d_temp;
		}
	    }

	    if(d >= 0.2*element_length){
		y_picked.appendRow(y_sorted.getRow(i));
		vel_picked.appendRow(vel_sorted.getRow(i));
		E_picked.appendRow(E_sorted.getRow(i));
		F_picked.appendRow(F_sorted.getRow(i));
		J_picked.appendRow(J_sorted.getRow(i));
		n_picked.appendRow(n_sorted.getRow(i));
	    }

	}

	std::list<REAL> d_sort;
	for(int i=1; i<y_picked.num_row+1; ++i){
	    REAL d_sort_val = pow( (y_picked(i,1)-y0(1,1))*(y_picked(i,1)-y0(1,1))
			        + (y_picked(i,2)-y0(1,2))*(y_picked(i,2)-y0(1,2)), 0.5);

	    d_sort.push_back(d_sort_val);
	}

	std::vector<int> i_sort = sort30(d_sort.begin(), d_sort.end());

	matrix y_picked_temp = y_picked;
	matrix vel_picked_temp = vel_picked;
	matrix E_picked_temp = E_picked;
	matrix F_picked_temp = F_picked;
	matrix J_picked_temp = J_picked;
	matrix n_picked_temp = n_picked;
	for(int i=0; i<y_picked_temp.num_row; ++i){
	    y_picked(i+1, 1) = y_picked_temp(i_sort[i], 1);
	    y_picked(i+1, 2) = y_picked_temp(i_sort[i], 2);

	    vel_picked(i+1, 1) = vel_picked_temp(i_sort[i], 1);
	    vel_picked(i+1, 2) = vel_picked_temp(i_sort[i], 2);

	    E_picked(i+1, 1) = E_picked_temp(i_sort[i], 1);
	    E_picked(i+1, 2) = E_picked_temp(i_sort[i], 2);

	    F_picked(i+1, 1) = F_picked_temp(i_sort[i], 1);
	    F_picked(i+1, 2) = F_picked_temp(i_sort[i], 2);

	    J_picked(i+1, 1) = J_picked_temp(i_sort[i], 1);

	    n_picked(i+1, 1) = n_picked_temp(i_sort[i], 1);
	    n_picked(i+1, 2) = n_picked_temp(i_sort[i], 2);
	}

	int m = length(y_picked);


	if(m<4){
	    std::cout << "not enough sampled points to build order 3 polynomial" << std::endl;
//	    exit(-1);
	    continue;
	}

	matrix y_picked3 = y0;
	for(int i=2; i<length(y_sorted)+1; ++i){
	    REAL d;
	    for(int j=1; j<y_picked3.num_row+1; ++j){
	    	REAL d_temp = pow( (y_picked3(j,1)-y_sorted(i,1))*(y_picked3(j,1)-y_sorted(i,1))
			         + (y_picked3(j,2)-y_sorted(i,2))*(y_picked3(j,2)-y_sorted(i,2)), 0.5);

		if(j==1){
		    d = d_temp;
		}
		else{
		    if(d>d_temp)
			d = d_temp;
		}
	    }
	
	    if(d >= 0.2*element_length){
		y_picked3.appendRow(y_sorted.getRow(i));
	    }

	}

	d_sort.clear();
//	d_sort = y_picked3;
	for(int i=1; i<y_picked3.num_row+1; ++i){
	    REAL d_sort_val = pow( (y_picked3(i,1)-y0(1,1))*(y_picked3(i,1)-y0(1,1))
			         + (y_picked3(i,2)-y0(1,2))*(y_picked3(i,2)-y0(1,2)), 0.5);
	    d_sort.push_back(d_sort_val);
	}

	i_sort.clear();
	i_sort = sort30(d_sort.begin(), d_sort.end());

	matrix y_picked3_temp = y_picked3;
	for(int i=0; i<y_picked3_temp.num_row; ++i){
	    y_picked3(i+1, 1) = y_picked3_temp(i_sort[i], 1);
	}

	int m3 = length(y_picked3);

	if(m>=20){
	    m3 = 20;

	    matrix y_picked3_temp = y_picked3;
	    y_picked3 = zeros(20,2);
	    for(int i=1; i<21; ++i){
	    	y_picked3(i,1) = y_picked3_temp(i,1);
		y_picked3(i,2) = y_picked3_temp(i,2);
	    }
	}

	if(m>=9){
	    m = 9;

	    matrix y_picked_temp = y_picked;
	    matrix vel_picked_temp = vel_picked;
	    matrix E_picked_temp = E_picked;
	    matrix F_picked_temp = F_picked;
	    matrix J_picked_temp = J_picked;
	    matrix n_picked_temp = n_picked;
	    y_picked = zeros(9,2);
	    vel_picked = zeros(9,2);
	    E_picked = zeros(9,2);
	    F_picked = zeros(9,2);
	    J_picked = zeros(9,1);
	    n_picked = zeros(9,2);

	    for(int i=1; i<10; ++i){
	    	y_picked(i,1) = y_picked_temp(i,1);
		y_picked(i,2) = y_picked_temp(i,2);

	    	vel_picked(i,1) = vel_picked_temp(i,1);
		vel_picked(i,2) = vel_picked_temp(i,2);

	    	E_picked(i,1) = E_picked_temp(i,1);
		E_picked(i,2) = E_picked_temp(i,2);

	    	F_picked(i,1) = F_picked_temp(i,1);
		F_picked(i,2) = F_picked_temp(i,2);

	    	J_picked(i,1) = J_picked_temp(i,1);

	    	n_picked(i,1) = n_picked_temp(i,1);
		n_picked(i,2) = n_picked_temp(i,2);
	    }
	}

	//find coordinates of picked points in the local sytem defined by n0_old
	matrix y_local;
	for(int i=1; i<length(y_picked)+1; ++i){
	    matrix temp_local;
	    matrix temp_n0(2,2);
	    temp_n0(1,1) = n0_old(1,2); temp_n0(1,2) = -n0_old(1,1);
	    temp_n0(2,1) = n0_old(1,1); temp_n0(2,2) = n0_old(1,2);
	    matrix temp_y(2,1);
	    temp_y(1,1) = y_picked(i,1)-y0(1,1); temp_y(2,1) = y_picked(i,2)-y0(1,2);
	    temp_local = temp_n0*temp_y;
	    
	    y_local.appendRow(temp_local.getCol(1));
  	}

	matrix y_local3;
	for(int i=1; i<length(y_picked3)+1; ++i){
	    matrix temp_local;
	    matrix temp_n0(2,2);
	    temp_n0(1,1) = n0_old(1,2); temp_n0(1,2) = -n0_old(1,1);
	    temp_n0(2,1) = n0_old(1,1); temp_n0(2,2) = n0_old(1,2);
	    matrix temp_y(2,1);
	    temp_y(1,1) = y_picked3(i,1)-y0(1,1); temp_y(2,1) = y_picked3(i,2)-y0(1,2);
	    temp_local = temp_n0*temp_y;
	    
	    y_local3.appendRow(temp_local.getCol(1));
  	}
	int num_i = length(y_picked3);

	REAL xmin = y_local(1,1);
	REAL xmax = y_local(1,1);
	for(int i=2; i<y_local.num_row+1; ++i){
	    if(xmin > y_local(i,1))
		xmin = y_local(i,1);
	    if(xmax < y_local(i,1))
		xmax = y_local(i,1);
	}

	// build the polynomial that locally parameterizes the membrane (order 2)
	REAL sum_x = 0; REAL sum_x2 = 0; REAL sum_x3 = 0; REAL sum_x4 = 0;
	REAL sum_y = 0; REAL sum_xy = 0; REAL sum_x2y = 0;
	REAL sum_vx = 0; REAL sum_xvx = 0; REAL sum_x2vx = 0;
	REAL sum_vy = 0; REAL sum_xvy = 0; REAL sum_x2vy = 0;
	REAL sum_E1 = 0; REAL sum_xE1 = 0; REAL sum_x2E1 = 0;
	REAL sum_E2 = 0; REAL sum_xE2 = 0; REAL sum_x2E2 = 0;
	REAL sum_F1 = 0; REAL sum_xF1 = 0; REAL sum_x2F1 = 0;
	REAL sum_F2 = 0; REAL sum_xF2 = 0; REAL sum_x2F2 = 0; 
	REAL sum_J = 0; REAL sum_xJ = 0; REAL sum_x2J= 0;
	REAL sum_n1 = 0; REAL sum_xn1 = 0; REAL sum_x2n1 = 0;
	REAL sum_n2 = 0; REAL sum_xn2 = 0; REAL sum_x2n2 = 0;
	for(int i=1; i<y_local.num_row+1; ++i){
	    sum_x     += y_local(i,1);
	    sum_x2    += y_local(i,1)*y_local(i,1);
	    sum_x3    += y_local(i,1)*y_local(i,1)*y_local(i,1);
	    sum_x4    += y_local(i,1)*y_local(i,1)*y_local(i,1)*y_local(i,1);

	    sum_y     += y_local(i,2);
	    sum_xy    += y_local(i,1)*y_local(i,2);
	    sum_x2y   += y_local(i,1)*y_local(i,1)*y_local(i,2);
	
	    sum_vx    += vel_picked(i,1);
	    sum_xvx   += y_local(i,1)*vel_picked(i,1);
	    sum_x2vx  += y_local(i,1)*y_local(i,1)*vel_picked(i,1);
	
	    sum_vy    += vel_picked(i,2);
	    sum_xvy   += y_local(i,1)*vel_picked(i,2);
	    sum_x2vy  += y_local(i,1)*y_local(i,1)*vel_picked(i,2);

	    sum_E1    += E_picked(i,1);
	    sum_xE1   += y_local(i,1)*E_picked(i,1);
	    sum_x2E1  += y_local(i,1)*y_local(i,1)*E_picked(i,1);

	    sum_E2    += E_picked(i,2);
	    sum_xE2   += y_local(i,1)*E_picked(i,2);
	    sum_x2E2  += y_local(i,1)*y_local(i,1)*E_picked(i,2);

	    sum_F1    += F_picked(i,1);
	    sum_xF1   += y_local(i,1)*F_picked(i,1);
	    sum_x2F1  += y_local(i,1)*y_local(i,1)*F_picked(i,1);

	    sum_F2    += F_picked(i,2);
	    sum_xF2   += y_local(i,1)*F_picked(i,2);
	    sum_x2F2  += y_local(i,1)*y_local(i,1)*F_picked(i,2);

	    sum_J     += J_picked(i,1);
	    sum_xJ    += y_local(i,1)*J_picked(i,1);
	    sum_x2J   += y_local(i,1)*y_local(i,1)*J_picked(i,1);

	    sum_n1    += n_picked(i,1);
	    sum_xn1   += y_local(i,1)*n_picked(i,1);
	    sum_x2n1  += y_local(i,1)*y_local(i,1)*n_picked(i,1);

	    sum_n2    += n_picked(i,2);
	    sum_xn2   += y_local(i,1)*n_picked(i,2);
	    sum_x2n2  += y_local(i,1)*y_local(i,1)*n_picked(i,2);
	}
	
	matrix temp_sum_xs(3,3);
	temp_sum_xs(1,1) = m;      temp_sum_xs(1,2) = sum_x;  temp_sum_xs(1,3) = sum_x2;
	temp_sum_xs(2,1) = sum_x;  temp_sum_xs(2,2) = sum_x2; temp_sum_xs(2,3) = sum_x3;
	temp_sum_xs(3,1) = sum_x2; temp_sum_xs(3,2) = sum_x3; temp_sum_xs(3,3) = sum_x4;

	// a_coeff
	matrix temp_xy(3,1);
	temp_xy(1,1) = sum_y; temp_xy(2,1) = sum_xy; temp_xy(3,1) = sum_x2y;
	matrix a_coeff = (temp_sum_xs%temp_xy).getTrans();	// % is left division "\"
	(*it)->setA_coeff(a_coeff);

	// a_coeff_vel
	matrix temp_vx(3,1);
	temp_vx(1,1) = sum_vx; temp_vx(2,1) = sum_xvx; temp_vx(3,1) = sum_x2vx;
	matrix a_coeff_vel_1 = (temp_sum_xs%temp_vx).getTrans();

	matrix temp_vy(3,1);
	temp_vy(1,1) = sum_vy; temp_vy(2,1) = sum_xvy; temp_vy(3,1) = sum_x2vy;
	matrix a_coeff_vel_2 = (temp_sum_xs%temp_vy).getTrans();

	matrix a_coeff_vel(1,6);
	a_coeff_vel(1,1) = a_coeff_vel_1(1,1); a_coeff_vel(1,3) = a_coeff_vel_1(1,2); a_coeff_vel(1,5) = a_coeff_vel_1(1,3);
	a_coeff_vel(1,2) = a_coeff_vel_2(1,1); a_coeff_vel(1,4) = a_coeff_vel_2(1,2); a_coeff_vel(1,6) = a_coeff_vel_2(1,3);
	(*it)->setA_coeff_vel(a_coeff_vel);

	// a_coeff_E
	matrix temp_xE1(3,1);
	temp_xE1(1,1) = sum_E1; temp_xE1(2,1) = sum_xE1; temp_xE1(3,1) = sum_x2E1;
	matrix a_coeff_E_1 = (temp_sum_xs%temp_xE1).getTrans();

	matrix temp_xE2(3,1);
	temp_xE2(1,1) = sum_E2; temp_xE2(2,1) = sum_xE2; temp_xE2(3,1) = sum_x2E2;
	matrix a_coeff_E_2 = (temp_sum_xs%temp_xE2).getTrans();

	matrix a_coeff_E(1,6);
	a_coeff_E(1,1) = a_coeff_E_1(1,1); a_coeff_E(1,3) = a_coeff_E_1(1,2); a_coeff_E(1,5) = a_coeff_E_1(1,3);
	a_coeff_E(1,2) = a_coeff_E_2(1,1); a_coeff_E(1,4) = a_coeff_E_2(1,2); a_coeff_E(1,6) = a_coeff_E_2(1,3);
	(*it)->setA_coeff_E(a_coeff_E);

	// a_coeff_F
	matrix temp_xF1(3,1);
	temp_xF1(1,1) = sum_F1; temp_xF1(2,1) = sum_xF1; temp_xF1(3,1) = sum_x2F1;
	matrix a_coeff_F_1 = (temp_sum_xs%temp_xF1).getTrans();

	matrix temp_xF2(3,1);
	temp_xF2(1,1) = sum_F2; temp_xF2(2,1) = sum_xF2; temp_xF2(3,1) = sum_x2F2;
	matrix a_coeff_F_2 = (temp_sum_xs%temp_xF2).getTrans();

	matrix a_coeff_F(1,6);
	a_coeff_F(1,1) = a_coeff_F_1(1,1); a_coeff_F(1,3) = a_coeff_F_1(1,2); a_coeff_F(1,5) = a_coeff_F_1(1,3);
	a_coeff_F(1,2) = a_coeff_F_2(1,1); a_coeff_F(1,4) = a_coeff_F_2(1,2); a_coeff_F(1,6) = a_coeff_F_2(1,3);
	(*it)->setA_coeff_F(a_coeff_F);

	// a_coeff_J
	matrix temp_xJ(3,1);
	temp_xJ(1,1) = sum_J; temp_xJ(2,1) = sum_xJ; temp_xJ(3,1) = sum_x2J;
	matrix a_coeff_J = (temp_sum_xs%temp_xJ).getTrans();
	(*it)->setA_coeff_J(a_coeff_J);

	// a_coeff_n
	matrix temp_xn1(3,1);
	temp_xn1(1,1) = sum_n1; temp_xn1(2,1) = sum_xn1; temp_xn1(3,1) = sum_x2n1;
	matrix a_coeff_n_1 = (temp_sum_xs%temp_xn1).getTrans();

	matrix temp_xn2(3,1);
	temp_xn2(1,1) = sum_n2; temp_xn2(2,1) = sum_xn2; temp_xn2(3,1) = sum_x2n2;
	matrix a_coeff_n_2 = (temp_sum_xs%temp_xn2).getTrans();

	matrix a_coeff_n(1,6);
	a_coeff_n(1,1) = a_coeff_n_1(1,1); a_coeff_n(1,3) = a_coeff_n_1(1,2); a_coeff_n(1,5) = a_coeff_n_1(1,3);
	a_coeff_n(1,2) = a_coeff_n_2(1,1); a_coeff_n(1,4) = a_coeff_n_2(1,2); a_coeff_n(1,6) = a_coeff_n_2(1,3);
	(*it)->setA_coeff_n(a_coeff_n);

	// build the polynomial that locally parameterizes the membrane (order 3)

	REAL sum3_x = 0;  REAL sum3_x2 = 0; REAL sum3_x3 = 0; 
	REAL sum3_x4 = 0; REAL sum3_x5 = 0; REAL sum3_x6 = 0;
	REAL sum3_y = 0;  REAL sum3_xy = 0; REAL sum3_x2y = 0; REAL sum3_x3y = 0;
	for(int i=1; i<y_local3.num_row+1; ++i){
	    sum3_x     += y_local3(i,1);
	    sum3_x2    += y_local3(i,1)*y_local3(i,1);
	    sum3_x3    += y_local3(i,1)*y_local3(i,1)*y_local3(i,1);
	    sum3_x4    += y_local3(i,1)*y_local3(i,1)*y_local3(i,1)*y_local3(i,1);
	    sum3_x5    += y_local3(i,1)*y_local3(i,1)*y_local3(i,1)*y_local3(i,1)*y_local3(i,1);
	    sum3_x6    += y_local3(i,1)*y_local3(i,1)*y_local3(i,1)*y_local3(i,1)*y_local3(i,1)*y_local3(i,1);

	    sum3_y     += y_local3(i,2);
	    sum3_xy    += y_local3(i,1)*y_local3(i,2);
	    sum3_x2y   += y_local3(i,1)*y_local3(i,1)*y_local3(i,2);
	    sum3_x3y   += y_local3(i,1)*y_local3(i,1)*y_local3(i,1)*y_local3(i,2);
	}

	matrix temp_sum3_x(4,4);
	temp_sum3_x(1,1) = m3;       temp_sum3_x(1,2) = sum3_x;  temp_sum3_x(1,3) = sum3_x2; temp_sum3_x(1,4) = sum3_x3;
	temp_sum3_x(2,1) = sum3_x;   temp_sum3_x(2,2) = sum3_x2; temp_sum3_x(2,3) = sum3_x3; temp_sum3_x(2,4) = sum3_x4;
	temp_sum3_x(3,1) = sum3_x2;  temp_sum3_x(3,2) = sum3_x3; temp_sum3_x(3,3) = sum3_x4; temp_sum3_x(3,4) = sum3_x5;
	temp_sum3_x(4,1) = sum3_x3;  temp_sum3_x(4,2) = sum3_x4; temp_sum3_x(4,3) = sum3_x5; temp_sum3_x(4,4) = sum3_x6;

	matrix temp_3xy(4,1);
	temp_3xy(1,1) = sum3_y; temp_3xy(2,1) = sum3_xy; temp_3xy(3,1) = sum3_x2y; temp_3xy(4,1) = sum3_x3y;

	matrix a_coeff3 = (temp_sum3_x%temp_3xy).getTrans();
	(*it)->setA_coeff3(a_coeff3);

	//Compute the new foot point (closest to P, so equivalent to rebuilding the distance function)
	//first, find P in local coordinates

	matrix temp_n0(2,2);
	temp_n0(1,1) = n0_old(1,2); temp_n0(1,2) = -n0_old(1,1);
	temp_n0(2,1) = n0_old(1,1); temp_n0(2,2) = n0_old(1,2);
	matrix p_local = temp_n0*(p-y0).getTrans();	// p_local is 2x1

	a_coeff = (*it)->getA_coeff();
	std::complex<double> a = 2.0*a_coeff(1,3)*a_coeff(1,3);
	std::complex<double> b = 3.0*a_coeff(1,2)*a_coeff(1,3);
	std::complex<double> c = 1.0+a_coeff(1,2)*a_coeff(1,2)+2.0*a_coeff(1,3)*a_coeff(1,1)-2.0*a_coeff(1,3)*p_local(2,1);
	std::complex<double> d = a_coeff(1,2)*a_coeff(1,1)-a_coeff(1,2)*p_local(2,1)-p_local(1,1);

	std::complex<double> i(0, 1);

   	std::complex<double> x1 = -b/(3.0*a) - (pow(2.0,1.0/3.0)*(-b*b + 3.0*a*c)) /(3.0*a*pow( -2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt( 4.0*pow(-b*b + 3.0*a*c, 3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d, 2.0) ), 1.0/3.0) ) 
+ pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt(4.0*pow(-b*b + 3.0*a*c,3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d,2.0) ), 1.0/3.0) /  (3.0*pow(2.0,1.0/3.0)*a);

	std::complex<double> x2 = -b/(3.0*a) + ((1.0+ i*3.0)*(-b*b + 3.0*a*c)) / ( 3.0*pow(2.0,2.0/3.0)*a*pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt(4.0*pow(-b*b + 3.0*a*c,3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d,2.0) ), 1.0/3.0) ) 
- (1.0 - i*sqrt(3.0))*pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt(4.0*pow(-b*b + 3.0*a*c, 3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d,2.0) ), 1.0/3.0) / (6.0*pow(2.0,1.0/3.0)*a);

    	std::complex<double> x3 = -b/(3.0*a) + ( (1.0 - i*sqrt(3.0))*(-b*b + 3.0*a*c)) / (3.0*pow(2.0,2.0/3.0)*a*pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt(4.0*pow(-b*b + 3.0*a*c,3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d,2.0) ), 1.0/3.0) )
 - (1.0 + i*sqrt(3.0))*pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt(4.0*pow(-b*b + 3.0*a*c,3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d,2.0) ), 1.0/3.0)  / (6.0*pow(2.0,1.0/3.0)*a);


	//eliminate imaginary part if small enough

	if( imag(x1)<1e-8 ){
	    x1 = real(x1);
	}

	if( imag(x1)<1e-8 ){
	    x2 = real(x2);
	}

	if( imag(x1)<1e-8){
	    x3 = real(x3);
	}

	//find the roots (minimum distance)

	matrix roots;	// row vector
	if(imag(x1)==0){	// real
	    matrix temp(1,1);
	    temp(1,1) = real(x1);
	    roots.appendCol(temp.getCol(1));
	}

	if(imag(x2)==0){	// real
	    matrix temp(1,1);
	    temp(1,1) = real(x2);
	    roots.appendCol(temp.getCol(1));
	}

	if(imag(x3)==0){	// real
	    matrix temp(1,1);
	    temp(1,1) = real(x3);
	    roots.appendCol(temp.getCol(1));
	}

	matrix distance;	// row vector
	for(int k=1; k<length(roots)+1; ++k){
	    matrix temp(1,1);
	    temp(1,1) = sqrt( (roots(1,k)-p_local(1,1))*(roots(1,k)-p_local(1,1))
			   + pow( a_coeff(1,1)+a_coeff(1,2)*roots(1,k)+a_coeff(1,3)*roots(1,k)*roots(1,k)
			       -  p_local(2,1), 2) );
	    distance.appendCol(temp.getCol(1));
	}

	REAL distanceMin = distance(1,1);
	int j=1;
	for(int k=2; k<distance.num_col+1; ++k){
	    if(distanceMin > distance(1,k)){
		distanceMin = distance(1,k);
		j = k;
	    }
	}



	if(  !(roots(1,j)>2*xmin && roots(1,j)<2*xmax)  ){
	    std::cout << "p: " << p(1,1) << ", " << p(1,2) << std::endl;
	    std::cout << "foot point outside of bound [xmin xmax]" << std::endl;
	    exit(-1);
	}

	//find foot point in global coordinates
	temp_n0 = zeros(2,2);	
	temp_n0(1,1) = n0_old(1,2); temp_n0(1,2) = -n0_old(1,1);
	temp_n0(2,1) = n0_old(1,1); temp_n0(2,2) = n0_old(1,2);

	matrix temp_roots(2,1); 
	temp_roots(1,1) = roots(1,j);
	matrix temp_vec(3,1);
	temp_vec(1,1) = 1; temp_vec(2,1) = roots(1,j); temp_vec(3,1) = roots(1,j)*roots(1,j);
	matrix temp_mat = a_coeff*temp_vec;
	temp_roots(2,1) = temp_mat(1,1);
	matrix Foot_points_new = (temp_n0%temp_roots).getTrans() + y0;
	(*it)->setFoot_points(Foot_points_new);

//	a_coeff_vel = (*it)->getA_coeff_vel();
	matrix vel_Foot_points_new(1,2);
	vel_Foot_points_new(1,1) = a_coeff_vel(1,1)+a_coeff_vel(1,3)*roots(1,j)+a_coeff_vel(1,5)*roots(1,j)*roots(1,j);
	vel_Foot_points_new(1,2) = a_coeff_vel(1,2)+a_coeff_vel(1,4)*roots(1,j)+a_coeff_vel(1,6)*roots(1,j)*roots(1,j);
	(*it)->setVel_Foot_points_membrane(vel_Foot_points_new);

//	a_coeff_E
	matrix E_Foot_points_new(1,2);
	E_Foot_points_new(1,1) = a_coeff_E(1,1)+a_coeff_E(1,3)*roots(1,j)+a_coeff_E(1,5)*roots(1,j)*roots(1,j);
	E_Foot_points_new(1,2) = a_coeff_E(1,2)+a_coeff_E(1,4)*roots(1,j)+a_coeff_E(1,6)*roots(1,j)*roots(1,j);
	(*it)->setE_Foot_points(E_Foot_points_new);

	matrix F_Foot_points_new(1,2);
	F_Foot_points_new(1,1) = a_coeff_F(1,1)+a_coeff_F(1,3)*roots(1,j)+a_coeff_F(1,5)*roots(1,j)*roots(1,j);
	F_Foot_points_new(1,2) = a_coeff_F(1,2)+a_coeff_F(1,4)*roots(1,j)+a_coeff_F(1,6)*roots(1,j)*roots(1,j);
	(*it)->setF_Foot_points(F_Foot_points_new);
	
	REAL J_Foot_points_new;
	J_Foot_points_new = a_coeff_J(1,1)+a_coeff_J(1,2)*roots(1,j)+a_coeff_J(1,3)*roots(1,j)*roots(1,j);
	(*it)->setJ_Foot_points(J_Foot_points_new);

	//find normal at foot point and sign the distance
	REAL mat_val = ((Foot_points_new-p)/(Foot_points_new-p).getNorm()*n0_old.getTrans())(1,1);
	matrix n_Foot_points_new = sign( mat_val )
				 * (Foot_points_new-p)/(Foot_points_new-p).getNorm();
	(*it)->setN_Foot_points(n_Foot_points_new);

	distanceMin = -distanceMin*sign( mat_val );
	(*it)->setDistanceMin(distanceMin);	

	//find curvature at Foot point
	if(plain_axi=="axi1"){
		if( fabs(Foot_points_new(1,1)+Lx*0.5)>=1e-4 ){
		    REAL a = pow(1+pow(a_coeff3(1,2)+2*a_coeff3(1,3)*roots(1,j) 
			             + 3*a_coeff3(1,4)*roots(1,j)*roots(1,j), 2), 0.5);

		    matrix curvature_Foot_points(1,2);
		    curvature_Foot_points(1,1) = (2*a_coeff3(1,3)+6*a_coeff3(1,4)*roots(1,j))/(a*a*a);
		    curvature_Foot_points(1,2) = -(1.0/a)*n_Foot_points_new(1,1)/(Foot_points_new(1,1)+0.5*Lx);
		    (*it)->setCurvature_Foot_points(curvature_Foot_points);

		    matrix metric_cov_Foot_points(1,2);
		    metric_cov_Foot_points(1,1) = a*a;
		    metric_cov_Foot_points(1,2) = (0.5*Lx+Foot_points_new(1,1))*(0.5*Lx+Foot_points_new(1,1));
		    (*it)->setMetric_cov_Foot_points(metric_cov_Foot_points);
		}
		else{
		    REAL a = pow(1+pow(a_coeff3(1,2)+2*a_coeff3(1,3)*roots(1,j) 
			             + 3*a_coeff3(1,4)*roots(1,j)*roots(1,j), 2), 0.5);

		    matrix curvature_Foot_points(1,2);
		    curvature_Foot_points(1,1) = (2*a_coeff3(1,3)+6*a_coeff3(1,4)*roots(1,j))/(a*a*a);
		    curvature_Foot_points(1,2) = (2*a_coeff3(1,3)+6*a_coeff3(1,4)*roots(1,j))/(a*a*a);
		    (*it)->setCurvature_Foot_points(curvature_Foot_points);

		    matrix metric_cov_Foot_points(1,2);
		    metric_cov_Foot_points(1,1) = a*a;
		    metric_cov_Foot_points(1,2) = (0.5*Lx+Foot_points_new(1,1))*(0.5*Lx+Foot_points_new(1,1));
		    (*it)->setMetric_cov_Foot_points(metric_cov_Foot_points);
		}

    	} // if "axi1"

	if(plain_axi=="axi2"){
		if( fabs(Foot_points_new(1,1)+Lx*0.5)>=1e-3 ){
		    REAL a = pow(1+pow(a_coeff3(1,2)+2*a_coeff3(1,3)*roots(1,j) 
			             + 3*a_coeff3(1,4)*roots(1,j)*roots(1,j), 2), 0.5);

		    matrix curvature_Foot_points(1,2);
		    curvature_Foot_points(1,1) = (2*a_coeff3(1,3)+6*a_coeff3(1,4)*roots(1,j))/(a*a*a);
		    curvature_Foot_points(1,2) = -(1.0/a)*n_Foot_points_new(1,1)/(Foot_points_new(1,1)+0.5*Lx);
		    (*it)->setCurvature_Foot_points(curvature_Foot_points);

		    matrix metric_cov_Foot_points(1,2);
		    metric_cov_Foot_points(1,1) = a*a;
		    metric_cov_Foot_points(1,2) = (0.5*Lx+Foot_points_new(1,1))*(0.5*Lx+Foot_points_new(1,1));
		    (*it)->setMetric_cov_Foot_points(metric_cov_Foot_points);
		}
		else{
		    REAL a = pow(1+pow(a_coeff3(1,2)+2*a_coeff3(1,3)*roots(1,j) 
			             + 3*a_coeff3(1,4)*roots(1,j)*roots(1,j), 2), 0.5);

		    matrix curvature_Foot_points(1,2);
		    curvature_Foot_points(1,1) = (2*a_coeff3(1,3)+6*a_coeff3(1,4)*roots(1,j))/(a*a*a);
		    curvature_Foot_points(1,2) = 0;
		    (*it)->setCurvature_Foot_points(curvature_Foot_points);

		    matrix metric_cov_Foot_points(1,2);
		    metric_cov_Foot_points(1,1) = a*a;
		    metric_cov_Foot_points(1,2) = (0.5*Lx+Foot_points_new(1,1))*(0.5*Lx+Foot_points_new(1,1));
		    (*it)->setMetric_cov_Foot_points(metric_cov_Foot_points);
		}

	} // if "axi2"

	if(plain_axi=="plain"){
		    REAL a = pow(1+pow(a_coeff3(1,2)+2*a_coeff3(1,3)*roots(1,j) 
			             + 3*a_coeff3(1,4)*roots(1,j)*roots(1,j), 2), 0.5);

		    matrix curvature_Foot_points(1,2);
		    curvature_Foot_points(1,1) = (2*a_coeff3(1,3)+6*a_coeff3(1,4)*roots(1,j))/(a*a*a);
		    curvature_Foot_points(1,2) = 0;
		    (*it)->setCurvature_Foot_points(curvature_Foot_points);

		    matrix metric_cov_Foot_points(1,2);
		    metric_cov_Foot_points(1,1) = a*a;
		    metric_cov_Foot_points(1,2) = 1;
		    (*it)->setMetric_cov_Foot_points(metric_cov_Foot_points);	
	} // if "plain"

	matrix Foot_points_test = (*it)->getFoot_points();
    	if( isnan_mat(Foot_points_test) ){
	    std::list<activepoint*>::iterator it_temp = it; it_temp--;
	    matrix Foot_points_set = (*it_temp)->getFoot_points();
	    (*it)->setFoot_points(Foot_points_set); 
	}


    } // Active_points for loop


    matrix Foot_points_new, n_Foot_points_new, vel_Foot_points_new, E_Foot_points_new, F_Foot_points_new;
    matrix J_Foot_points_new(Active_points.size(), 1);

    matrix Foot_points_sym, n_Foot_points_sym, vel_Foot_points_sym, E_Foot_points_sym, F_Foot_points_sym;
    matrix J_Foot_points_sym;

    ii=1;
    for(std::list<activepoint*>::const_iterator it=Active_points.begin(); it!=Active_points.end(); ++it, ++ii){
	Foot_points_new.appendRow( ((*it)->getFoot_points()).getRow(1) );
	n_Foot_points_new.appendRow( ((*it)->getN_Foot_points()).getRow(1) );
	vel_Foot_points_new.appendRow( ((*it)->getVel_Foot_points_membrane()).getRow(1) );
 	E_Foot_points_new.appendRow( ((*it)->getE_Foot_points()).getRow(1) );
	F_Foot_points_new.appendRow( ((*it)->getF_Foot_points()).getRow(1) );
	J_Foot_points_new(ii,1) = (*it)->getJ_Foot_points();
    }

    if(plain_axi=="axi1"){

	    matrix temp_Foot_points = Foot_points_new;
	    Foot_points_sym = zeros(2*temp_Foot_points.num_row,2);
	    for(int i=1; i<temp_Foot_points.num_row+1; ++i){
		Foot_points_sym(i,1) = temp_Foot_points(i,1);
		Foot_points_sym(i,2) = temp_Foot_points(i,2);
	    }
	    for(int i=1+temp_Foot_points.num_row; i<2*temp_Foot_points.num_row+1; ++i){
		Foot_points_sym(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row,1))-0.5*Lx;
		Foot_points_sym(i,2) = temp_Foot_points(i-temp_Foot_points.num_row,2);
	    }

	    // n_Foot_points
	    matrix n_Foot_points_col1;
	    n_Foot_points_col1.appendCol(n_Foot_points_new.getCol(1));
	    n_Foot_points_col1 = -1*n_Foot_points_col1;
	    matrix sym_points2;
	    sym_points2.appendCol(n_Foot_points_col1.getCol(1));
 	    sym_points2.appendCol(n_Foot_points_new.getCol(2));
	    n_Foot_points_sym = n_Foot_points_new;
	    for(int i=1; i<sym_points2.num_row+1; ++i){
		n_Foot_points_sym.appendRow(sym_points2.getRow(i));
	    }

	    // vel_Foot_points_membrane
	    matrix vel_Foot_points_col1;
	    vel_Foot_points_col1.appendCol(vel_Foot_points_new.getCol(1));
	    vel_Foot_points_col1 = -1*vel_Foot_points_col1;
	    matrix sym_points3;
	    sym_points3.appendCol(vel_Foot_points_col1.getCol(1));
 	    sym_points3.appendCol(vel_Foot_points_new.getCol(2));
	    vel_Foot_points_sym = vel_Foot_points_new;
	    for(int i=1; i<sym_points3.num_row+1; ++i){
		vel_Foot_points_sym.appendRow(sym_points3.getRow(i));
	    }

	    // E_Foot_points
	    matrix temp_mat = E_Foot_points_new;
	    E_Foot_points_sym = E_Foot_points_new;
	    for(int i=1; i<temp_mat.num_row+1; ++i){
		E_Foot_points_sym.appendRow(temp_mat.getRow(i));
	    }

	    // F_Foot_points
	    matrix temp_mat2 = F_Foot_points_new;
	    F_Foot_points_sym = F_Foot_points_new;
	    for(int i=1; i<temp_mat2.num_row+1; ++i){
		F_Foot_points_sym.appendRow(temp_mat2.getRow(i));
	    }

	    // J_Foot_points
	    matrix temp_mat3 = J_Foot_points_new;
	    J_Foot_points_sym = J_Foot_points_new;
	    for(int i=1; i<temp_mat3.num_row+1; ++i){
		J_Foot_points_sym.appendRow(temp_mat3.getRow(i));
	    }
	    
    } // if "axi1"

    if(plain_axi=="axi2"){
	    matrix temp_Foot_points = Foot_points_new;
	    Foot_points_sym = zeros(3*temp_Foot_points.num_row,2);
	    for(int i=1; i<temp_Foot_points.num_row+1; ++i){
		Foot_points_sym(i,1) = temp_Foot_points(i,1);
		Foot_points_sym(i,2) = temp_Foot_points(i,2);
	    }
	    for(int i=1+temp_Foot_points.num_row; i<2*temp_Foot_points.num_row+1; ++i){
		Foot_points_sym(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row,1))-0.5*Lx;
		Foot_points_sym(i,2) = temp_Foot_points(i-temp_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_Foot_points.num_row; i<3*temp_Foot_points.num_row+1; ++i){
		Foot_points_sym(i,1) = -(-Lx*0.5+temp_Foot_points(i-2*temp_Foot_points.num_row,1))+0.5*Lx;
		Foot_points_sym(i,2) = temp_Foot_points(i-2*temp_Foot_points.num_row,2);
	    }

	    // n_Foot_points
	    matrix temp_n_Foot_points = n_Foot_points_new;
	    n_Foot_points_sym = zeros(3*temp_n_Foot_points.num_row,2);
	    for(int i=1; i<temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points_sym(i,1) = temp_n_Foot_points(i,1);
		n_Foot_points_sym(i,2) = temp_n_Foot_points(i,2);
	    }
	    for(int i=1+temp_n_Foot_points.num_row; i<2*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points_sym(i,1) = -temp_n_Foot_points(i-temp_n_Foot_points.num_row,1);
		n_Foot_points_sym(i,2) = temp_n_Foot_points(i-temp_n_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_n_Foot_points.num_row; i<3*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points_sym(i,1) = -temp_n_Foot_points(i-2*temp_n_Foot_points.num_row,1);
		n_Foot_points_sym(i,2) = temp_n_Foot_points(i-2*temp_n_Foot_points.num_row,2);
	    }

	    // vel_Foot_points_membrane
	    matrix temp_vel_Foot_points = vel_Foot_points_new;
	    vel_Foot_points_sym = zeros(3*temp_vel_Foot_points.num_row,2);
	    for(int i=1; i<temp_vel_Foot_points.num_row+1; ++i){
		vel_Foot_points_sym(i,1) = temp_vel_Foot_points(i,1);
		vel_Foot_points_sym(i,2) = temp_vel_Foot_points(i,2);
	    }
	    for(int i=1+temp_vel_Foot_points.num_row; i<2*temp_vel_Foot_points.num_row+1; ++i){
		vel_Foot_points_sym(i,1) = -temp_vel_Foot_points(i-temp_vel_Foot_points.num_row,1);
		vel_Foot_points_sym(i,2) = temp_vel_Foot_points(i-temp_vel_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_vel_Foot_points.num_row; i<3*temp_vel_Foot_points.num_row+1; ++i){
		vel_Foot_points_sym(i,1) = -temp_vel_Foot_points(i-2*temp_vel_Foot_points.num_row,1);
		vel_Foot_points_sym(i,2) = temp_vel_Foot_points(i-2*temp_vel_Foot_points.num_row,2);
	    }

	    // E_Foot_points
	    matrix temp_E_Foot_points = E_Foot_points_new;
	    E_Foot_points_sym = zeros(3*temp_E_Foot_points.num_row,2);
	    for(int i=1; i<temp_E_Foot_points.num_row+1; ++i){
		E_Foot_points_sym(i,1) = temp_E_Foot_points(i,1);
		E_Foot_points_sym(i,2) = temp_E_Foot_points(i,2);
	    }
	    for(int i=1+temp_E_Foot_points.num_row; i<2*temp_E_Foot_points.num_row+1; ++i){
		E_Foot_points_sym(i,1) = temp_E_Foot_points(i-temp_E_Foot_points.num_row,1);
		E_Foot_points_sym(i,2) = temp_E_Foot_points(i-temp_E_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_E_Foot_points.num_row; i<3*temp_E_Foot_points.num_row+1; ++i){
		E_Foot_points_sym(i,1) = temp_E_Foot_points(i-2*temp_E_Foot_points.num_row,1);
		E_Foot_points_sym(i,2) = temp_E_Foot_points(i-2*temp_E_Foot_points.num_row,2);
	    }

	    // F_Foot_points
	    matrix temp_F_Foot_points = F_Foot_points_new;
	    F_Foot_points_sym = zeros(3*temp_F_Foot_points.num_row,2);
	    for(int i=1; i<temp_F_Foot_points.num_row+1; ++i){
		F_Foot_points_sym(i,1) = temp_F_Foot_points(i,1);
		F_Foot_points_sym(i,2) = temp_F_Foot_points(i,2);
	    }
	    for(int i=1+temp_F_Foot_points.num_row; i<2*temp_F_Foot_points.num_row+1; ++i){
		F_Foot_points_sym(i,1) = temp_F_Foot_points(i-temp_F_Foot_points.num_row,1);
		F_Foot_points_sym(i,2) = temp_F_Foot_points(i-temp_F_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_F_Foot_points.num_row; i<3*temp_F_Foot_points.num_row+1; ++i){
		F_Foot_points_sym(i,1) = temp_F_Foot_points(i-2*temp_F_Foot_points.num_row,1);
		F_Foot_points_sym(i,2) = temp_F_Foot_points(i-2*temp_F_Foot_points.num_row,2);
	    }

	    // J_Foot_points
	    matrix temp_J_Foot_points = J_Foot_points_new;
	    J_Foot_points_sym = zeros(3*temp_J_Foot_points.num_row,1);
	    for(int i=1; i<temp_J_Foot_points.num_row+1; ++i){
		J_Foot_points_sym(i,1) = temp_J_Foot_points(i,1);
	    }
	    for(int i=1+temp_J_Foot_points.num_row; i<2*temp_J_Foot_points.num_row+1; ++i){
		J_Foot_points_sym(i,1) = temp_J_Foot_points(i-temp_J_Foot_points.num_row,1);
	    }
	    for(int i=1+2*temp_J_Foot_points.num_row; i<3*temp_J_Foot_points.num_row+1; ++i){
		J_Foot_points_sym(i,1) = temp_J_Foot_points(i-2*temp_J_Foot_points.num_row,1);
	    }

    } // if "axi2"

    if(plain_axi=="plain"){
            Foot_points_sym = Foot_points_new;
            n_Foot_points_sym = n_Foot_points_new;
            vel_Foot_points_sym = vel_Foot_points_new;
            E_Foot_points_sym = E_Foot_points_new;
            F_Foot_points_sym = F_Foot_points_new;
            J_Foot_points_sym = J_Foot_points_new;

    } // if "plain"


    // Phi92 = zeros(1,numnode9)
    for(std::vector<node*>::iterator it=node9.begin(); it<node9.end(); ++it){
	(*it)->setPhi92(0);
	(*it)->setPhi42(0);
    }

    for(std::list<activepoint*>::iterator it=Active_points.begin(); it!=Active_points.end(); ++it){

	matrix p0 = (*it)->getCoords();

	matrix d = zeros(numnode9, 1);
	std::vector<node*> neighbors; int ii=1;
	for(std::vector<node*>::const_iterator in=node9.begin(); in<node9.end(); ++in, ++ii){
	    d(ii,1) = pow( ((*in)->getX()-p0(1,1))*((*in)->getX()-p0(1,1))
			+ ((*in)->getY()-p0(1,2))*((*in)->getY()-p0(1,2)), 0.5);
	    if( fabs(d(ii,1))<1.5*element_length )
		neighbors.push_back(*in);
	}



	for(std::vector<node*>::const_iterator ik=neighbors.begin(); ik<neighbors.end(); ++ik){

	    if( (*ik)->getAct_point()==false ){

		matrix p(1,2);
		p(1,1) = (*ik)->getX();
		p(1,2) = (*ik)->getY();

//		matrix D2 = zeros(Foot_points_sym.num_row, 1);
		std::list<REAL> D2;
		for(int i=1; i<Foot_points_sym.num_row+1; ++i){
	    	    REAL D2_val = pow( (Foot_points_sym(i,1)-p(1,1))*(Foot_points_sym(i,1)-p(1,1))
			             + (Foot_points_sym(i,2)-p(1,2))*(Foot_points_sym(i,2)-p(1,2)), 0.5);
		    D2.push_back(D2_val);
		}

		std::vector<int> IDX2 = sort30(D2.begin(), D2.end());

		matrix y_sorted = zeros(30,2);
   		matrix vel_sorted = zeros(30,2);
		matrix E_sorted = zeros(30,2);
		matrix F_sorted = zeros(30,2);
		matrix J_sorted = zeros(30,1);
		matrix n_sorted = zeros(30,2);


		for(int i=0; i<30; ++i){
		    y_sorted(i+1, 1) = Foot_points_sym(IDX2[i], 1);
		    y_sorted(i+1, 2) = Foot_points_sym(IDX2[i], 2);
	
		    vel_sorted(i+1, 1) = vel_Foot_points_sym(IDX2[i], 1);
		    vel_sorted(i+1, 2) = vel_Foot_points_sym(IDX2[i], 2);
	
		    E_sorted(i+1, 1) = E_Foot_points_sym(IDX2[i], 1);
		    E_sorted(i+1, 2) = E_Foot_points_sym(IDX2[i], 2);
	
		    F_sorted(i+1, 1) = F_Foot_points_sym(IDX2[i], 1);
		    F_sorted(i+1, 2) = F_Foot_points_sym(IDX2[i], 2);
	
		    J_sorted(i+1, 1) = J_Foot_points_sym(IDX2[i], 1);

		    n_sorted(i+1, 1) = n_Foot_points_sym(IDX2[i], 1);
		    n_sorted(i+1, 2) = n_Foot_points_sym(IDX2[i], 2);
		}

		//n0_old and y0 are the origin and orientation of the local coordinate system

		matrix n0_old(1,2);
		n0_old(1,1) = n_Foot_points_sym(IDX2[0], 1);
		n0_old(1,2) = n_Foot_points_sym(IDX2[0], 2);
	
		matrix y0(1,2);
		y0(1,1) = y_sorted(1,1);
		y0(1,2) = y_sorted(1,2);
	
		matrix v0(1,2);
		v0(1,1) = vel_sorted(1,1);
		v0(1,2) = vel_sorted(1,2);
	
		matrix E0(1,2);
		E0(1,1) = E_sorted(1,1);
		E0(1,2) = E_sorted(1,2);

		matrix F0(1,2);
		F0(1,1) = F_sorted(1,1);
		F0(1,2) = F_sorted(1,2);
	
		REAL J0 = J_sorted(1,1);
	
		matrix n0(1,2);
		n0(1,1) = n_sorted(1,1);
		n0(1,2) = n_sorted(1,2);
	
		//we now select a few neighboring Foot_points to construct the 
		//polynomials that interpolate the interface

	      	matrix y_picked = y0;
	    	matrix vel_picked = v0;
	    	matrix E_picked = E0;
	    	matrix F_picked = F0;
	    	matrix J_picked(1,1); J_picked(1,1) = J0;
	    	matrix n_picked = n0;

		for(int i=2; i<length(y_sorted)+1; ++i){
		    REAL d;
		    for(int j=1; j<y_picked.num_row+1; ++j){
		    	REAL d_temp = pow( (y_picked(j,1)-y_sorted(i,1))*(y_picked(j,1)-y_sorted(i,1))
      				         + (y_picked(j,2)-y_sorted(i,2))*(y_picked(j,2)-y_sorted(i,2)), 0.5);

			if(j==1){
		       	    d = d_temp;
			}
			else{
		    	    if(d>d_temp)
				d = d_temp;
			}
	   	    }

	    	    if(d >= 0.2*element_length){
		    	y_picked.appendRow(y_sorted.getRow(i));
		   	vel_picked.appendRow(vel_sorted.getRow(i));
		    	E_picked.appendRow(E_sorted.getRow(i));
		    	F_picked.appendRow(F_sorted.getRow(i));
		    	J_picked.appendRow(J_sorted.getRow(i));
		    	n_picked.appendRow(n_sorted.getRow(i));
	    	    }

		}

//		matrix d_sort = y_picked;
		std::list<REAL> d_sort;
		for(int i=1; i<y_picked.num_row+1; ++i){
	    	    REAL d_sort_val = pow( (y_picked(i,1)-y0(1,1))*(y_picked(i,1)-y0(1,1))
			                 + (y_picked(i,2)-y0(1,2))*(y_picked(i,2)-y0(1,2)), 0.5);

		    d_sort.push_back(d_sort_val);
		}

		std::vector<int> i_sort = sort30(d_sort.begin(), d_sort.end());

		matrix y_picked_temp = y_picked;
		matrix vel_picked_temp = vel_picked;
		matrix E_picked_temp = E_picked;
		matrix F_picked_temp = F_picked;
		matrix J_picked_temp = J_picked;
		matrix n_picked_temp = n_picked;
		for(int i=0; i<y_picked_temp.num_row; ++i){
		    y_picked(i+1, 1) = y_picked_temp(i_sort[i], 1);
		    y_picked(i+1, 2) = y_picked_temp(i_sort[i], 2);

		    vel_picked(i+1, 1) = vel_picked_temp(i_sort[i], 1);
		    vel_picked(i+1, 2) = vel_picked_temp(i_sort[i], 2);
	
		    E_picked(i+1, 1) = E_picked_temp(i_sort[i], 1);
		    E_picked(i+1, 2) = E_picked_temp(i_sort[i], 2);
	
		    F_picked(i+1, 1) = F_picked_temp(i_sort[i], 1);
		    F_picked(i+1, 2) = F_picked_temp(i_sort[i], 2);
	
		    J_picked(i+1, 1) = J_picked_temp(i_sort[i], 1);
	
		    n_picked(i+1, 1) = n_picked_temp(i_sort[i], 1);
		    n_picked(i+1, 2) = n_picked_temp(i_sort[i], 2);
		}

		int m = length(y_picked);

		matrix y_picked3 = y0;
		for(int i=2; i<length(y_sorted)+1; ++i){
		    REAL d;
		    for(int j=1; j<y_picked3.num_row+1; ++j){
		    	REAL d_temp = pow( (y_picked3(j,1)-y_sorted(i,1))*(y_picked3(j,1)-y_sorted(i,1))
				         + (y_picked3(j,2)-y_sorted(i,2))*(y_picked3(j,2)-y_sorted(i,2)), 0.5);
	
			if(j==1){
			    d = d_temp;
			}
			else{
			    if(d>d_temp)
				d = d_temp;
			}
		    }
	
		    if(d >= 0.2*element_length){
			y_picked3.appendRow(y_sorted.getRow(i));
		    }

		}

		d_sort.clear();
//		d_sort = y_picked3;
		for(int i=1; i<y_picked3.num_row+1; ++i){
		    REAL d_sort_val = pow( (y_picked3(i,1)-y0(1,1))*(y_picked3(i,1)-y0(1,1))
				         + (y_picked3(i,2)-y0(1,2))*(y_picked3(i,2)-y0(1,2)), 0.5);
		    d_sort.push_back(d_sort_val);
		}

		i_sort.clear();
		i_sort = sort30(d_sort.begin(), d_sort.end());
		convert_sort(i_sort.begin(), i_sort.end());

		matrix y_picked3_temp = y_picked3;
		for(int i=0; i<y_picked3_temp.num_row; ++i){
		    y_picked3(i+1, 1) = y_picked3_temp(i_sort[i], 1);
		}

		int m3 = length(y_picked3);

		if(m<4){
		    std::cout << "not enough sampled points to build order 3 polynomial" << std::endl;
//		    exit(-1);
		    continue;
		}

		y_picked3 = y_picked;
		if(m>=20){
	    	    m3 = 20;

	   	    matrix y_picked3_temp = y_picked3;
	    	    y_picked3 = zeros(20,2);
	    	    for(int i=1; i<21; ++i){
	    		y_picked3(i,1) = y_picked3_temp(i,1);
			y_picked3(i,2) = y_picked3_temp(i,2);
	    	    }
		}

		if(m>=9){
	    	    m = 9;

	    	    matrix y_picked_temp = y_picked;
	    	    matrix vel_picked_temp = vel_picked;
	    	    matrix E_picked_temp = E_picked;
	    	    matrix F_picked_temp = F_picked;
	    	    matrix J_picked_temp = J_picked;
	    	    matrix n_picked_temp = n_picked;
	    	    y_picked = zeros(9,2);
	    	    vel_picked = zeros(9,2);
	    	    E_picked = zeros(9,2);
	    	    F_picked = zeros(9,2);
	    	    J_picked = zeros(9,1);
	    	    n_picked = zeros(9,2);

	    	    for(int i=1; i<10; ++i){
	    	    	y_picked(i,1) = y_picked_temp(i,1);
		    	y_picked(i,2) = y_picked_temp(i,2);
	
	    	    	vel_picked(i,1) = vel_picked_temp(i,1);
		    	vel_picked(i,2) = vel_picked_temp(i,2);

	    	    	E_picked(i,1) = E_picked_temp(i,1);
		    	E_picked(i,2) = E_picked_temp(i,2);

	    	    	F_picked(i,1) = F_picked_temp(i,1);
		    	F_picked(i,2) = F_picked_temp(i,2);
	
	    	    	J_picked(i,1) = J_picked_temp(i,1);

	    	    	n_picked(i,1) = n_picked_temp(i,1);
		    	n_picked(i,2) = n_picked_temp(i,2);
	       	    }
	    	}

		//find coordinates of picked points in the local sytem defined by n0_old
		matrix y_local;
		for(int i=1; i<length(y_picked)+1; ++i){
		    matrix temp_local;
		    matrix temp_n0(2,2);
		    temp_n0(1,1) = n0_old(1,2); temp_n0(1,2) = -n0_old(1,1);
		    temp_n0(2,1) = n0_old(1,1); temp_n0(2,2) = n0_old(1,2);
		    matrix temp_y(2,1);
		    temp_y(1,1) = y_picked(i,1)-y0(1,1); temp_y(2,1) = y_picked(i,2)-y0(1,2);
		    temp_local = temp_n0*temp_y;
		    
		    y_local.appendRow(temp_local.getCol(1));
  		}

		matrix y_local3;
		for(int i=1; i<length(y_picked3)+1; ++i){
		    matrix temp_local;
		    matrix temp_n0(2,2);
		    temp_n0(1,1) = n0_old(1,2); temp_n0(1,2) = -n0_old(1,1);
		    temp_n0(2,1) = n0_old(1,1); temp_n0(2,2) = n0_old(1,2);
		    matrix temp_y(2,1);
		    temp_y(1,1) = y_picked3(i,1)-y0(1,1); temp_y(2,1) = y_picked3(i,2)-y0(1,2);
		    temp_local = temp_n0*temp_y;
		    
		    y_local3.appendRow(temp_local.getCol(1));
  		}
		int num_i = length(y_picked3);


		REAL xmin = y_local(1,1);
		REAL xmax = y_local(1,1);
		for(int i=2; i<y_local.num_row+1; ++i){
		    if(xmin > y_local(i,1))
			xmin = y_local(i,1);
		    if(xmax < y_local(i,1))
			xmax = y_local(i,1);
		}

		// build the polynomial that locally parameterizes the membrane (order 2)
		REAL sum_x = 0; REAL sum_x2 = 0; REAL sum_x3 = 0; REAL sum_x4 = 0;
		REAL sum_y = 0; REAL sum_xy = 0; REAL sum_x2y = 0;
		REAL sum_vx = 0; REAL sum_xvx = 0; REAL sum_x2vx = 0;
		REAL sum_vy = 0; REAL sum_xvy = 0; REAL sum_x2vy = 0;
		REAL sum_E1 = 0; REAL sum_xE1 = 0; REAL sum_x2E1 = 0;
		REAL sum_E2 = 0; REAL sum_xE2 = 0; REAL sum_x2E2 = 0;
		REAL sum_F1 = 0; REAL sum_xF1 = 0; REAL sum_x2F1 = 0;
		REAL sum_F2 = 0; REAL sum_xF2 = 0; REAL sum_x2F2 = 0; 
		REAL sum_J = 0; REAL sum_xJ = 0; REAL sum_x2J= 0;
		REAL sum_n1 = 0; REAL sum_xn1 = 0; REAL sum_x2n1 = 0;
		REAL sum_n2 = 0; REAL sum_xn2 = 0; REAL sum_x2n2 = 0;
		for(int i=1; i<y_local.num_row+1; ++i){
		    sum_x     += y_local(i,1);
		    sum_x2    += y_local(i,1)*y_local(i,1);
		    sum_x3    += y_local(i,1)*y_local(i,1)*y_local(i,1);
		    sum_x4    += y_local(i,1)*y_local(i,1)*y_local(i,1)*y_local(i,1);
	
		    sum_y     += y_local(i,2);
		    sum_xy    += y_local(i,1)*y_local(i,2);
		    sum_x2y   += y_local(i,1)*y_local(i,1)*y_local(i,2);
	
		    sum_vx    += vel_picked(i,1);
		    sum_xvx   += y_local(i,1)*vel_picked(i,1);
		    sum_x2vx  += y_local(i,1)*y_local(i,1)*vel_picked(i,1);
		
		    sum_vy    += vel_picked(i,2);
		    sum_xvy   += y_local(i,1)*vel_picked(i,2);
		    sum_x2vy  += y_local(i,1)*y_local(i,1)*vel_picked(i,2);

		    sum_E1    += E_picked(i,1);
		    sum_xE1   += y_local(i,1)*E_picked(i,1);
		    sum_x2E1  += y_local(i,1)*y_local(i,1)*E_picked(i,1);
	
		    sum_E2    += E_picked(i,2);
		    sum_xE2   += y_local(i,1)*E_picked(i,2);
		    sum_x2E2  += y_local(i,1)*y_local(i,1)*E_picked(i,2);
	
		    sum_F1    += F_picked(i,1);
		    sum_xF1   += y_local(i,1)*F_picked(i,1);
		    sum_x2F1  += y_local(i,1)*y_local(i,1)*F_picked(i,1);
	
		    sum_F2    += F_picked(i,2);
		    sum_xF2   += y_local(i,1)*F_picked(i,2);
		    sum_x2F2  += y_local(i,1)*y_local(i,1)*F_picked(i,2);
	
		    sum_J     += J_picked(i,1);
		    sum_xJ    += y_local(i,1)*J_picked(i,1);
		    sum_x2J   += y_local(i,1)*y_local(i,1)*J_picked(i,1);
	
		    sum_n1    += n_picked(i,1);
		    sum_xn1   += y_local(i,1)*n_picked(i,1);
		    sum_x2n1  += y_local(i,1)*y_local(i,1)*n_picked(i,1);
	
		    sum_n2    += n_picked(i,2);
		    sum_xn2   += y_local(i,1)*n_picked(i,2);
		    sum_x2n2  += y_local(i,1)*y_local(i,1)*n_picked(i,2);
		}
		
		matrix temp_sum_xs(3,3);
		temp_sum_xs(1,1) = m;      temp_sum_xs(1,2) = sum_x;  temp_sum_xs(1,3) = sum_x2;
		temp_sum_xs(2,1) = sum_x;  temp_sum_xs(2,2) = sum_x2; temp_sum_xs(2,3) = sum_x3;
		temp_sum_xs(3,1) = sum_x2; temp_sum_xs(3,2) = sum_x3; temp_sum_xs(3,3) = sum_x4;
	
		// a_coeff
		matrix temp_xy(3,1);
		temp_xy(1,1) = sum_y; temp_xy(2,1) = sum_xy; temp_xy(3,1) = sum_x2y;
		matrix a_coeff_temp = (temp_sum_xs%temp_xy).getTrans();
	
		// a_coeff_vel
		matrix temp_vx(3,1);
		temp_vx(1,1) = sum_vx; temp_vx(2,1) = sum_xvx; temp_vx(3,1) = sum_x2vx;
		matrix a_coeff_vel_1 = (temp_sum_xs%temp_vx).getTrans();
	
		matrix temp_vy(3,1);
		temp_vy(1,1) = sum_vy; temp_vy(2,1) = sum_xvy; temp_vy(3,1) = sum_x2vy;
		matrix a_coeff_vel_2 = (temp_sum_xs%temp_vy).getTrans();
	
		matrix a_coeff_vel_temp(1,6);
		a_coeff_vel_temp(1,1) = a_coeff_vel_1(1,1); 
		a_coeff_vel_temp(1,3) = a_coeff_vel_1(1,2);
	 	a_coeff_vel_temp(1,5) = a_coeff_vel_1(1,3);
		a_coeff_vel_temp(1,2) = a_coeff_vel_2(1,1);
		a_coeff_vel_temp(1,4) = a_coeff_vel_2(1,2); 
		a_coeff_vel_temp(1,6) = a_coeff_vel_2(1,3);
	
		// a_coeff_E
		matrix temp_xE1(3,1);
		temp_xE1(1,1) = sum_E1; temp_xE1(2,1) = sum_xE1; temp_xE1(3,1) = sum_x2E1;
		matrix a_coeff_E_1 = (temp_sum_xs%temp_xE1).getTrans();
	
		matrix temp_xE2(3,1);
		temp_xE2(1,1) = sum_E2; temp_xE2(2,1) = sum_xE2; temp_xE2(3,1) = sum_x2E2;
		matrix a_coeff_E_2 = (temp_sum_xs%temp_xE2).getTrans();
	
		matrix a_coeff_E_temp(1,6);
		a_coeff_E_temp(1,1) = a_coeff_E_1(1,1);
	 	a_coeff_E_temp(1,3) = a_coeff_E_1(1,2);
 		a_coeff_E_temp(1,5) = a_coeff_E_1(1,3);
		a_coeff_E_temp(1,2) = a_coeff_E_2(1,1); 
		a_coeff_E_temp(1,4) = a_coeff_E_2(1,2); 
		a_coeff_E_temp(1,6) = a_coeff_E_2(1,3);

		// a_coeff_F
		matrix temp_xF1(3,1);
		temp_xF1(1,1) = sum_F1; temp_xF1(2,1) = sum_xF1; temp_xF1(3,1) = sum_x2F1;
		matrix a_coeff_F_1 = (temp_sum_xs%temp_xF1).getTrans();
	
		matrix temp_xF2(3,1);
		temp_xF2(1,1) = sum_F2; temp_xF2(2,1) = sum_xF2; temp_xF2(3,1) = sum_x2F2;
		matrix a_coeff_F_2 = (temp_sum_xs%temp_xF2).getTrans();
	
		matrix a_coeff_F_temp(1,6);
		a_coeff_F_temp(1,1) = a_coeff_F_1(1,1); 
		a_coeff_F_temp(1,3) = a_coeff_F_1(1,2); 
		a_coeff_F_temp(1,5) = a_coeff_F_1(1,3);
		a_coeff_F_temp(1,2) = a_coeff_F_2(1,1); 
		a_coeff_F_temp(1,4) = a_coeff_F_2(1,2); 
		a_coeff_F_temp(1,6) = a_coeff_F_2(1,3);

		// a_coeff_J
		matrix temp_xJ(3,1);
		temp_xJ(1,1) = sum_J; temp_xJ(2,1) = sum_xJ; temp_xJ(3,1) = sum_x2J;
		matrix a_coeff_J_temp = (temp_sum_xs%temp_xJ).getTrans();
	
		// a_coeff_n
		matrix temp_xn1(3,1);
		temp_xn1(1,1) = sum_n1; temp_xn1(2,1) = sum_xn1; temp_xn1(3,1) = sum_x2n1;
		matrix a_coeff_n_1 = (temp_sum_xs%temp_xn1).getTrans();
	
		matrix temp_xn2(3,1);
		temp_xn2(1,1) = sum_n2; temp_xn2(2,1) = sum_xn2; temp_xn2(3,1) = sum_x2n2;
		matrix a_coeff_n_2 = (temp_sum_xs%temp_xn2).getTrans();
	
		matrix a_coeff_n_temp(1,6);
		a_coeff_n_temp(1,1) = a_coeff_n_1(1,1); 
		a_coeff_n_temp(1,3) = a_coeff_n_1(1,2); 
		a_coeff_n_temp(1,5) = a_coeff_n_1(1,3);
		a_coeff_n_temp(1,2) = a_coeff_n_2(1,1); 
		a_coeff_n_temp(1,4) = a_coeff_n_2(1,2); 
		a_coeff_n_temp(1,6) = a_coeff_n_2(1,3);

		// build the polynomial that locally parameterizes the membrane (order 3)
	
		REAL sum3_x = 0;  REAL sum3_x2 = 0; REAL sum3_x3 = 0; 
		REAL sum3_x4 = 0; REAL sum3_x5 = 0; REAL sum3_x6 = 0;
		REAL sum3_y = 0;  REAL sum3_xy = 0; REAL sum3_x2y = 0; REAL sum3_x3y = 0;
		for(int i=1; i<y_local3.num_row+1; ++i){
		    sum3_x     += y_local3(i,1);
		    sum3_x2    += y_local3(i,1)*y_local3(i,1);
		    sum3_x3    += y_local3(i,1)*y_local3(i,1)*y_local3(i,1);
		    sum3_x4    += y_local3(i,1)*y_local3(i,1)*y_local3(i,1)*y_local3(i,1);
		    sum3_x5    += y_local3(i,1)*y_local3(i,1)*y_local3(i,1)*y_local3(i,1)*y_local3(i,1);
		    sum3_x5    += y_local3(i,1)*y_local3(i,1)*y_local3(i,1)*y_local3(i,1)*y_local3(i,1)*y_local3(i,1);

		    sum3_y     += y_local3(i,2);
		    sum3_xy    += y_local3(i,1)*y_local3(i,2);
		    sum3_x2y   += y_local3(i,1)*y_local3(i,1)*y_local3(i,2);
		    sum3_x3y   += y_local3(i,1)*y_local3(i,1)*y_local3(i,1)*y_local3(i,2);
		}

		matrix temp_sum3_x(4,4);
		temp_sum3_x(1,1) = m3;       temp_sum3_x(1,2) = sum3_x;  
		temp_sum3_x(1,3) = sum3_x2;  temp_sum3_x(1,4) = sum3_x3;
		temp_sum3_x(2,1) = sum3_x;   temp_sum3_x(2,2) = sum3_x2; 
		temp_sum3_x(2,3) = sum3_x3;  temp_sum3_x(2,4) = sum3_x4;
		temp_sum3_x(3,1) = sum3_x2;  temp_sum3_x(3,2) = sum3_x3; 
		temp_sum3_x(3,3) = sum3_x4;  temp_sum3_x(3,4) = sum3_x5;
		temp_sum3_x(4,1) = sum3_x3;  temp_sum3_x(4,2) = sum3_x4; 
		temp_sum3_x(4,3) = sum3_x5;  temp_sum3_x(4,4) = sum3_x6;

		matrix temp_3xy(4,1);	
		temp_3xy(1,1) = sum3_y; temp_3xy(2,1) = sum3_xy; 
		temp_3xy(3,1) = sum3_x2y; temp_3xy(4,1) = sum3_x3y;

		matrix a_coeff3_temp = (temp_sum3_x%temp_3xy).getTrans();
	
		//Compute the new foot point (closest to P, so equivalent to rebuilding the distance function)
		//first, find P in local coordinates
	
		matrix temp_n0(2,2);
		temp_n0(1,1) = n0_old(1,2); temp_n0(1,2) = -n0_old(1,1);
		temp_n0(2,1) = n0_old(1,1); temp_n0(2,2) = n0_old(1,2);
		matrix p_local = temp_n0*(p-y0).getTrans();	// p_local is 2x1
	
		std::complex<double> a = 2.0*a_coeff_temp(1,3)*a_coeff_temp(1,3);
		std::complex<double> b = 3.0*a_coeff_temp(1,2)*a_coeff_temp(1,3);
		std::complex<double> c = 1.0+a_coeff_temp(1,2)*a_coeff_temp(1,2)
					 + 2.0*a_coeff_temp(1,3)*a_coeff_temp(1,1)
					 - 2.0*a_coeff_temp(1,3)*p_local(2,1);
		std::complex<double> d = a_coeff_temp(1,2)*a_coeff_temp(1,1)
				       - a_coeff_temp(1,2)*p_local(2,1)-p_local(1,1);

		std::complex<double> i(0, 1);

   		std::complex<double> x1 = -b/(3.0*a) - (pow(2.0,1.0/3.0)*(-b*b + 3.0*a*c)) /(3.0*a*pow( -2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt( 4.0*pow(-b*b + 3.0*a*c, 3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d, 2.0) ), 1.0/3.0) ) 
+ pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt(4.0*pow(-b*b + 3.0*a*c,3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d,2.0) ), 1.0/3.0) /  (3.0*pow(2.0,1.0/3.0)*a);

		std::complex<double> x2 = -b/(3.0*a) + ((1.0 + i*3.0)*(-b*b + 3.0*a*c)) / ( 3.0*pow(2.0,2.0/3.0)*a*pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt(4.0*pow(-b*b + 3.0*a*c,3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d,2.0) ), 1.0/3.0) ) 
- (1.0 - i*sqrt(3.0))*pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt(4.0*pow(-b*b + 3.0*a*c, 3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d,2.0) ), 1.0/3.0) / (6.0*pow(2.0,1.0/3.0)*a);

    		std::complex<double> x3 = -b/(3.0*a) + ( (1.0 - i*sqrt(3.0))*(-b*b + 3.0*a*c)) / (3.0*pow(2.0,2.0/3.0)*a*pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt(4.0*pow(-b*b + 3.0*a*c,3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d,2.0) ), 1.0/3.0) )
 - (1.0 + i*sqrt(3.0))*pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt(4.0*pow(-b*b + 3.0*a*c,3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d,2.0) ), 1.0/3.0)  / (6.0*pow(2.0,1.0/3.0)*a);


		//eliminate imaginary part if small enough
	
		if( imag(x1)<1e-8 ){
	    	x1 = real(x1);
		}
	
		if( imag(x1)<1e-8 ){
		    x2 = real(x2);
		}
	
		if( imag(x1)<1e-8){
		    x3 = real(x3);
		}

		//find the roots (minimum distance)
	
		matrix roots;	// row vector
		if(imag(x1)==0){	// real
		    matrix temp(1,1);
		    temp(1,1) = real(x1);
		    roots.appendCol(temp.getCol(1));
		}

		if(imag(x2)==0){	// real
		    matrix temp(1,1);
		    temp(1,1) = real(x2);
		    roots.appendCol(temp.getCol(1));
		}
	
		if(imag(x3)==0){	// real
		    matrix temp(1,1);
		    temp(1,1) = real(x3);
		    roots.appendCol(temp.getCol(1));
		}

		matrix distance;	// row vector
		for(int k=1; k<length(roots)+1; ++k){
		    matrix temp(1,1);
		    temp(1,1) = sqrt( (roots(1,k)-p_local(1,1))*(roots(1,k)-p_local(1,1))
			   + pow( a_coeff_temp(1,1)+a_coeff_temp(1,2)*roots(1,k)+a_coeff_temp(1,3)*roots(1,k)*roots(1,k)
			       -  p_local(2,1), 2) );
		    distance.appendCol(temp.getCol(1));
		}
	
		if(min(distance)<=1.1*element_length){

		    activepoint* act_pt = new activepoint(*ik);

		    REAL distanceMin = distance(1,1);
		    int j=1;
		    for(int k=2; k<distance.num_col+1; ++k){
	    	    	if(distanceMin > distance(1,k)){
			    distanceMin = distance(1,k);
			    j = k;
	    	    	}
		    }
//		    (*it)->setDistanceMin(distanceMin);

//		    if(  !(roots(1,j)>2*xmin && roots(1,j)<2*xmax)  ){
//	    		std::cout << "p: " << p(1,1) << ", " << p(1,2) << std::endl;
//	    		std::cout << "foot point outside of bound [xmin xmax]" << std::endl;
//	   	 	exit(-1);
//		    }

		    //find foot point in global coordinates
		    temp_n0 = zeros(2,2);	
		    temp_n0(1,1) = n0_old(1,2); temp_n0(1,2) = -n0_old(1,1);
		    temp_n0(2,1) = n0_old(1,1); temp_n0(2,2) = n0_old(1,2);

		    matrix temp_roots(2,1); 
		    temp_roots(1,1) = roots(1,j);
		    matrix temp_vec(3,1);
		    temp_vec(1,1) = 1; temp_vec(2,1) = roots(1,j); temp_vec(3,1) = roots(1,j)*roots(1,j);
		    matrix temp_mat = a_coeff_temp*temp_vec;
		    temp_roots(2,1) = temp_mat(1,1);
		    matrix Foot_points_new = (temp_n0%temp_roots).getTrans() + y0;
//		    if(!isnan_mat(Foot_points_new)){
		        (act_pt)->setFoot_points(Foot_points_new);
//		    }
		    //find normal at foot point and sign the distance
		    REAL mat_val = ((Foot_points_new-p)/(Foot_points_new-p).getNorm()*n0_old.getTrans())(1,1);
		    matrix n_Foot_points_new = sign( mat_val ) * (Foot_points_new-p)/(Foot_points_new-p).getNorm();
		    (act_pt)->setN_Foot_points(n_Foot_points_new);

		    distanceMin = -distanceMin*sign( mat_val );
		    (act_pt)->setDistanceMin(distanceMin);	

		    (act_pt)->setA_coeff(a_coeff_temp);
		    (act_pt)->setA_coeff_vel(a_coeff_vel_temp);
		    (act_pt)->setA_coeff_E(a_coeff_E_temp);	
		    (act_pt)->setA_coeff_F(a_coeff_F_temp);	  
		    (act_pt)->setA_coeff_J(a_coeff_J_temp);  
		    (act_pt)->setA_coeff_n(a_coeff_n_temp);
		    (act_pt)->setA_coeff3(a_coeff3_temp);

		    //find curvature at Foot point
		    if(plain_axi=="axi1"){
			    if( fabs(Foot_points_new(1,1)+Lx*0.5)>=1e-4 ){
		   		REAL a = pow(1+pow(a_coeff3_temp(1,2)+2*a_coeff3_temp(1,3)*roots(1,j) 
			             	     + 3*a_coeff3_temp(1,4)*roots(1,j)*roots(1,j), 2), 0.5);

		    		matrix curvature_Foot_points(1,2);
		    		curvature_Foot_points(1,1) = (2*a_coeff3_temp(1,3)+6*a_coeff3_temp(1,4)*roots(1,j))/(a*a*a);
		    		curvature_Foot_points(1,2) = -(1.0/a)*n_Foot_points_new(1,1)/(Foot_points_new(1,1)+0.5*Lx);
		    		(act_pt)->setCurvature_Foot_points(curvature_Foot_points);

		    		matrix metric_cov_Foot_points(1,2);
		    		metric_cov_Foot_points(1,1) = a*a;
		    		metric_cov_Foot_points(1,2) = (0.5*Lx+Foot_points_new(1,1))*(0.5*Lx+Foot_points_new(1,1));
		    		(act_pt)->setMetric_cov_Foot_points(metric_cov_Foot_points);

			    }
			    else{
		   		REAL a = pow(1+pow(a_coeff3_temp(1,2)+2*a_coeff3_temp(1,3)*roots(1,j) 
			             	         + 3*a_coeff3_temp(1,4)*roots(1,j)*roots(1,j), 2), 0.5);

		    		matrix curvature_Foot_points(1,2);
		    		curvature_Foot_points(1,1) = (2*a_coeff3_temp(1,3)+6*a_coeff3_temp(1,4)*roots(1,j))/(a*a*a);
		    		curvature_Foot_points(1,2) = (2*a_coeff3_temp(1,3)+6*a_coeff3_temp(1,4)*roots(1,j))/(a*a*a);
		    		(act_pt)->setCurvature_Foot_points(curvature_Foot_points);

		    		matrix metric_cov_Foot_points(1,2);
		    		metric_cov_Foot_points(1,1) = a*a;
		    		metric_cov_Foot_points(1,2) = (0.5*Lx+Foot_points_new(1,1))*(0.5*Lx+Foot_points_new(1,1));
		    		(act_pt)->setMetric_cov_Foot_points(metric_cov_Foot_points);

			}

 		    } // if "axi1"

	    	    if(plain_axi=="axi2"){
			    if( fabs(Foot_points_new(1,1)+Lx*0.5)>=1e-3 ){
		    		REAL a = pow(1+pow(a_coeff3_temp(1,2)+2*a_coeff3_temp(1,3)*roots(1,j) 
			        		 + 3*a_coeff3_temp(1,4)*roots(1,j)*roots(1,j), 2), 0.5);

		    		matrix curvature_Foot_points(1,2);
		    		curvature_Foot_points(1,1) = (2*a_coeff3_temp(1,3)+6*a_coeff3_temp(1,4)*roots(1,j))/(a*a*a);
		    		curvature_Foot_points(1,2) = -(1.0/a)*n_Foot_points_new(1,1)/(Foot_points_new(1,1)+0.5*Lx);
		    		(act_pt)->setCurvature_Foot_points(curvature_Foot_points);

		    		matrix metric_cov_Foot_points(1,2);
		    		metric_cov_Foot_points(1,1) = a*a;
		    		metric_cov_Foot_points(1,2) = (0.5*Lx+Foot_points_new(1,1))*(0.5*Lx+Foot_points_new(1,1));
		    		(act_pt)->setMetric_cov_Foot_points(metric_cov_Foot_points);
			    }
			    else{
		    		REAL a = pow(1+pow(a_coeff3_temp(1,2)+2*a_coeff3_temp(1,3)*roots(1,j) 
			        	         + 3*a_coeff3_temp(1,4)*roots(1,j)*roots(1,j), 2), 0.5);

		    		matrix curvature_Foot_points(1,2);
			    	curvature_Foot_points(1,1) = (2*a_coeff3_temp(1,3)+6*a_coeff3_temp(1,4)*roots(1,j))/(a*a*a);
		    		curvature_Foot_points(1,2) = 0;
		    		(act_pt)->setCurvature_Foot_points(curvature_Foot_points);

		    		matrix metric_cov_Foot_points(1,2);
		   	 	metric_cov_Foot_points(1,1) = a*a;
		   	 	metric_cov_Foot_points(1,2) = (0.5*Lx+Foot_points_new(1,1))*(0.5*Lx+Foot_points_new(1,1));
		    		(act_pt)->setMetric_cov_Foot_points(metric_cov_Foot_points);
			    }

		    } // if "axi2"

	    	    if(plain_axi=="plain"){
		    	    REAL a = pow(1+pow(a_coeff3_temp(1,2)+2*a_coeff3_temp(1,3)*roots(1,j) 
			             	     + 3*a_coeff3_temp(1,4)*roots(1,j)*roots(1,j), 2), 0.5);

		    	    matrix curvature_Foot_points(1,2);
		    	    curvature_Foot_points(1,1) = (2*a_coeff3_temp(1,3)+6*a_coeff3_temp(1,4)*roots(1,j))/(a*a*a);
		    	    curvature_Foot_points(1,2) = 0;
		    	    (act_pt)->setCurvature_Foot_points(curvature_Foot_points);

		    	    matrix metric_cov_Foot_points(1,2);
		    	    metric_cov_Foot_points(1,1) = a*a;
		    	    metric_cov_Foot_points(1,2) = 1;
		    	    (act_pt)->setMetric_cov_Foot_points(metric_cov_Foot_points);	
		    } // if "plain"

		    Active_points.push_back(act_pt);
		} // if(min(distance)


	    } // not active point

	} // neigbors for loop

    } // Active_points for loop


    changeIsnanFoot_points();


    // Phi92(1,Active_points) = distanceMin;
    for(std::list<activepoint*>::iterator it=Active_points.begin(); it!=Active_points.end(); ++it){
	(*it)->setPhi92EqualDistanceMin();
    }

    for(std::vector<ele2dcoupled*>::iterator it2d=element.begin(); it2d<element.end(); ++it2d){
	(*it2d)->setPhi42EqualPhi92();
    }

    for(std::list<activepoint*>::iterator it=Active_points.begin(); it!=Active_points.end(); /*nothing*/){
	if( (*it)->getDistanceMin()>1.1*element_length ){ 
	    it = Active_points.erase(it);
	}
	else
	    ++it;
    }


} // Foot_point_ressampling_v2()


void surfacetension::find_laplace_beltrami_coeff(){

    matrix Coord_Active_points, Foot_points, n_Foot_points, curvature_Foot_points, metric_cov_Foot_points;
 
    for(std::list<activepoint*>::const_iterator it=Active_points.begin(); it!=Active_points.end(); ++it){

	Coord_Active_points.appendRow( ((*it)->getCoords()).getRow(1) );
	Foot_points.appendRow( ((*it)->getFoot_points()).getRow(1) );
	n_Foot_points.appendRow( ((*it)->getN_Foot_points()).getRow(1) );
	curvature_Foot_points.appendRow( ((*it)->getCurvature_Foot_points()).getRow(1) );
 	metric_cov_Foot_points.appendRow( ((*it)->getMetric_cov_Foot_points()).getRow(1) );
    }

    if(plain_axi=="axi1"){

	    matrix temp_Foot_points = Foot_points;
	    Foot_points = zeros(2*temp_Foot_points.num_row,2);
	    for(int i=1; i<temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = temp_Foot_points(i,1);
		Foot_points(i,2) = temp_Foot_points(i,2);
	    }
	    for(int i=1+temp_Foot_points.num_row; i<2*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row,1))-0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-temp_Foot_points.num_row,2);
	    }

	    // n_Foot_points
	    matrix n_Foot_points_col1;
	    n_Foot_points_col1.appendCol(n_Foot_points.getCol(1));
	    n_Foot_points_col1 = -1*n_Foot_points_col1;
	    matrix sym_points2;
	    sym_points2.appendCol(n_Foot_points_col1.getCol(1));
 	    sym_points2.appendCol(n_Foot_points.getCol(2));
	    for(int i=1; i<sym_points2.num_row+1; ++i){
		n_Foot_points.appendRow(sym_points2.getRow(i));
	    }

	    // curvature_Foot_points
	    matrix temp_mat = curvature_Foot_points;
	    for(int i=1; i<temp_mat.num_row+1; ++i){
		curvature_Foot_points.appendRow(temp_mat.getRow(i));
	    }

	    // metric_cov_Foot_points
	    matrix temp_mat2 = metric_cov_Foot_points;
	    for(int i=1; i<temp_mat2.num_row+1; ++i){
		metric_cov_Foot_points.appendRow(temp_mat2.getRow(i));
	    }

	    // Coord_Active_points
	    matrix temp_coord = Coord_Active_points;
	    Coord_Active_points = zeros(2*temp_coord.num_row,2);
	    for(int i=1; i<temp_coord.num_row+1; ++i){
		Coord_Active_points(i,1) = temp_coord(i,1);
		Coord_Active_points(i,2) = temp_coord(i,2);
	    }
	    for(int i=1+temp_coord.num_row; i<2*temp_coord.num_row+1; ++i){
		Coord_Active_points(i,1) = -(Lx*0.5+temp_coord(i-temp_coord.num_row,1))-0.5*Lx;
		Coord_Active_points(i,2) = temp_coord(i-temp_coord.num_row,2);
	    }
	    
    } // if "axi1"

    if(plain_axi=="axi2"){
	    matrix temp_Foot_points = Foot_points;
	    Foot_points = zeros(3*temp_Foot_points.num_row,2);
	    for(int i=1; i<temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = temp_Foot_points(i,1);
		Foot_points(i,2) = temp_Foot_points(i,2);
	    }
	    for(int i=1+temp_Foot_points.num_row; i<2*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row,1))-0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-temp_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_Foot_points.num_row; i<3*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(-Lx*0.5+temp_Foot_points(i-2*temp_Foot_points.num_row,1))+0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-2*temp_Foot_points.num_row,2);
	    }

	    // Coord_Active_points
	    matrix temp_coord = Coord_Active_points;
	    Coord_Active_points = zeros(3*temp_coord.num_row,2);
	    for(int i=1; i<temp_coord.num_row+1; ++i){
		Coord_Active_points(i,1) = temp_coord(i,1);
		Coord_Active_points(i,2) = temp_coord(i,2);
	    }
	    for(int i=1+temp_coord.num_row; i<2*temp_coord.num_row+1; ++i){
		Coord_Active_points(i,1) = -(Lx*0.5+temp_coord(i-temp_coord.num_row,1))-0.5*Lx;
		Coord_Active_points(i,2) = temp_coord(i-temp_coord.num_row,2);
	    }
	    for(int i=1+2*temp_coord.num_row; i<3*temp_coord.num_row+1; ++i){
		Coord_Active_points(i,1) = -(-Lx*0.5+temp_coord(i-2*temp_coord.num_row,1))+0.5*Lx;
		Coord_Active_points(i,2) = temp_coord(i-2*temp_coord.num_row,2);
	    }

	    // n_Foot_points
	    matrix temp_n_Foot_points = n_Foot_points;
	    n_Foot_points = zeros(3*temp_n_Foot_points.num_row,2);
	    for(int i=1; i<temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = temp_n_Foot_points(i,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i,2);
	    }
	    for(int i=1+temp_n_Foot_points.num_row; i<2*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = -temp_n_Foot_points(i-temp_n_Foot_points.num_row,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i-temp_n_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_n_Foot_points.num_row; i<3*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = -temp_n_Foot_points(i-2*temp_n_Foot_points.num_row,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i-2*temp_n_Foot_points.num_row,2);
	    }

	    // curvature_Foot_points
	    matrix temp_curv_Foot_points = curvature_Foot_points;
	    curvature_Foot_points = zeros(3*temp_curv_Foot_points.num_row,2);
	    for(int i=1; i<temp_curv_Foot_points.num_row+1; ++i){
		curvature_Foot_points(i,1) = temp_curv_Foot_points(i,1);
		curvature_Foot_points(i,2) = temp_curv_Foot_points(i,2);
	    }
	    for(int i=1+temp_curv_Foot_points.num_row; i<2*temp_curv_Foot_points.num_row+1; ++i){
		curvature_Foot_points(i,1) = temp_curv_Foot_points(i-temp_curv_Foot_points.num_row,1);
		curvature_Foot_points(i,2) = temp_curv_Foot_points(i-temp_curv_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_curv_Foot_points.num_row; i<3*temp_curv_Foot_points.num_row+1; ++i){
		curvature_Foot_points(i,1) = temp_curv_Foot_points(i-2*temp_curv_Foot_points.num_row,1);
		curvature_Foot_points(i,2) = temp_curv_Foot_points(i-2*temp_curv_Foot_points.num_row,2);
	    }

	    // metric_cov_Foot_points
	    matrix temp_metric_Foot_points = metric_cov_Foot_points;
	    metric_cov_Foot_points = zeros(3*temp_metric_Foot_points.num_row,2);
	    for(int i=1; i<temp_metric_Foot_points.num_row+1; ++i){
		metric_cov_Foot_points(i,1) = temp_metric_Foot_points(i,1);
		metric_cov_Foot_points(i,2) = temp_metric_Foot_points(i,2);
	    }
	    for(int i=1+temp_metric_Foot_points.num_row; i<2*temp_metric_Foot_points.num_row+1; ++i){
		metric_cov_Foot_points(i,1) = temp_metric_Foot_points(i-temp_metric_Foot_points.num_row,1);
		metric_cov_Foot_points(i,2) = temp_metric_Foot_points(i-temp_metric_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_metric_Foot_points.num_row; i<3*temp_metric_Foot_points.num_row+1; ++i){
		metric_cov_Foot_points(i,1) = temp_metric_Foot_points(i-2*temp_metric_Foot_points.num_row,1);
		metric_cov_Foot_points(i,2) = temp_metric_Foot_points(i-2*temp_metric_Foot_points.num_row,2);
	    }

    } // if "axi2"

    if(plain_axi=="plain"){
	
    } // if "plain"


    for(std::list<activepoint*>::iterator it=Active_points.begin(); it!=Active_points.end(); ++it){

	(*it)->initA_coeff_curvature_metric_cov();

	matrix p = (*it)->getCoords();

//	matrix D2 = zeros(Foot_points.num_row, 1);
	std::list<REAL> D2;
	for(int i=1; i<Foot_points.num_row+1; ++i){
	    REAL D2_val = pow( (Foot_points(i,1)-p(1,1))*(Foot_points(i,1)-p(1,1))
			     + (Foot_points(i,2)-p(1,2))*(Foot_points(i,2)-p(1,2)), 0.5);
	    D2.push_back(D2_val);
	}

	std::vector<int> IDX2 = sort30(D2.begin(), D2.end());

	matrix y_sorted = zeros(30,2);
   	matrix curvature_sorted = zeros(30,2);
	matrix metric_cov_sorted = zeros(30,2);
	matrix n_sorted = zeros(30,2);

	for(int i=0; i<30; ++i){
	    y_sorted(i+1, 1) = Foot_points(IDX2[i], 1);
	    y_sorted(i+1, 2) = Foot_points(IDX2[i], 2);

	    curvature_sorted(i+1, 1) = curvature_Foot_points(IDX2[i], 1);
	    curvature_sorted(i+1, 2) = curvature_Foot_points(IDX2[i], 2);

	    metric_cov_sorted(i+1, 1) = metric_cov_Foot_points(IDX2[i], 1);
	    metric_cov_sorted(i+1, 2) = metric_cov_Foot_points(IDX2[i], 2);

	    n_sorted(i+1, 1) = n_Foot_points(IDX2[i], 1);
	    n_sorted(i+1, 2) = n_Foot_points(IDX2[i], 2);
	}

	//n0_old and y0 are the origin and orientation of the local coordinate system

	matrix n0_old(1,2);
	n0_old(1,1) = n_Foot_points(IDX2[0], 1);
	n0_old(1,2) = n_Foot_points(IDX2[0], 2);

	matrix y0(1,2);
	y0(1,1) = y_sorted(1,1);
	y0(1,2) = y_sorted(1,2);

	matrix curvature0(1,2);
	curvature0(1,1) = curvature_sorted(1,1);
	curvature0(1,2) = curvature_sorted(1,2);

	matrix metric_cov0(1,2);
	metric_cov0(1,1) = metric_cov_sorted(1,1);
	metric_cov0(1,2) = metric_cov_sorted(1,2);

	matrix n0(1,2);
	n0(1,1) = n_sorted(1,1);
	n0(1,2) = n_sorted(1,2);

	//we now select a few neighboring Foot_points to construct the 
	//polynomials that interpolate the interface

      	matrix y_picked = y0;
    	matrix curvature_picked = curvature0;
    	matrix metric_cov_picked = metric_cov0;

	for(int i=2; i<length(y_sorted)+1; ++i){
	    REAL d;
	    for(int j=1; j<y_picked.num_row+1; ++j){
	    	REAL d_temp = pow( (y_picked(j,1)-y_sorted(i,1))*(y_picked(j,1)-y_sorted(i,1))
			         + (y_picked(j,2)-y_sorted(i,2))*(y_picked(j,2)-y_sorted(i,2)), 0.5);

		if(j==1){
		    d = d_temp;
		}
		else{
		    if(d>d_temp)
			d = d_temp;
		}
	    }

	    if(d >= 0.25*element_length){
		y_picked.appendRow(y_sorted.getRow(i));
		curvature_picked.appendRow(curvature_sorted.getRow(i));
		metric_cov_picked.appendRow(metric_cov_sorted.getRow(i));
	    }

	}

//	matrix d_sort = y_picked;
	std::list<REAL> d_sort;
	for(int i=1; i<y_picked.num_row+1; ++i){
	    REAL d_sort_val = pow( (y_picked(i,1)-y0(1,1))*(y_picked(i,1)-y0(1,1))
			         + (y_picked(i,2)-y0(1,2))*(y_picked(i,2)-y0(1,2)), 0.5);
	    d_sort.push_back(d_sort_val);
	}

	std::vector<int> i_sort = sort30(d_sort.begin(), d_sort.end());

	matrix y_picked_temp = y_picked;
	matrix curvature_picked_temp = curvature_picked;
	matrix metric_cov_picked_temp = metric_cov_picked;
	for(int i=0; i<y_picked_temp.num_row; ++i){
	    y_picked(i+1, 1) = y_picked_temp(i_sort[i], 1);
	    y_picked(i+1, 2) = y_picked_temp(i_sort[i], 2);

	    curvature_picked(i+1, 1) = curvature_picked_temp(i_sort[i], 1);
	    curvature_picked(i+1, 2) = curvature_picked_temp(i_sort[i], 2);

	    metric_cov_picked(i+1, 1) = metric_cov_picked_temp(i_sort[i], 1);
	}

	int m = length(y_picked);

	if(m<4){
	    std::cout << "not enough sampled points to build order 3 polynomial" << std::endl;
//	    exit(-1);
	    continue;
	}

	if(m>=20){
	    m = 20;

	    matrix y_picked_temp = y_picked;
	    matrix curvature_picked_temp = curvature_picked;
	    matrix metric_cov_picked_temp = metric_cov_picked;
	    y_picked = zeros(20,2);
	    curvature_picked = zeros(20,2);
	    metric_cov_picked = zeros(20,2);
	    for(int i=1; i<21; ++i){
	    	y_picked(i,1) = y_picked_temp(i,1);
		y_picked(i,2) = y_picked_temp(i,2);

	    	curvature_picked(i,1) = curvature_picked_temp(i,1);
		curvature_picked(i,2) = curvature_picked_temp(i,2);

	    	metric_cov_picked(i,1) = metric_cov_picked_temp(i,1);
		metric_cov_picked(i,2) = metric_cov_picked_temp(i,2);
	    }
	}

	//find coordinates of picked points in the local sytem defined by n0_old
	matrix y_local;
	for(int i=1; i<length(y_picked)+1; ++i){
	    matrix temp_local;
	    matrix temp_n0(2,2);
	    temp_n0(1,1) = n0_old(1,2); temp_n0(1,2) = -n0_old(1,1);
	    temp_n0(2,1) = n0_old(1,1); temp_n0(2,2) = n0_old(1,2);
	    matrix temp_y(2,1);
	    temp_y(1,1) = y_picked(i,1)-y0(1,1); temp_y(2,1) = y_picked(i,2)-y0(1,2);
	    temp_local = temp_n0*temp_y;
	    
	    y_local.appendRow(temp_local.getCol(1));
  	}

	REAL xmin = y_local(1,1);
	REAL xmax = y_local(1,1);
	for(int i=2; i<y_local.num_row+1; ++i){
	    if(xmin > y_local(i,1))
		xmin = y_local(i,1);
	    if(xmax < y_local(i,1))
		xmax = y_local(i,1);
	}

	// build the polynomial that locally parameterizes the membrane (order 2)
	REAL sum_x = 0;  REAL sum_x2 = 0; REAL sum_x3 = 0; 
	REAL sum_x4 = 0; REAL sum_x5 = 0; REAL sum_x6 = 0;

	REAL sum_curvature1 = 0; REAL sum_xcurvature1 = 0; REAL sum_x2curvature1 = 0; REAL sum_x3curvature1 = 0;
	REAL sum_curvature2 = 0; REAL sum_xcurvature2 = 0; REAL sum_x2curvature2 = 0; REAL sum_x3curvature2 = 0;

	REAL sum_metric_cov1 = 0; REAL sum_xmetric_cov1 = 0; REAL sum_x2metric_cov1 = 0;
	REAL sum_metric_cov2 = 0; REAL sum_xmetric_cov2 = 0; REAL sum_x2metric_cov2 = 0;

	for(int i=1; i<y_local.num_row+1; ++i){
	    sum_x     += y_local(i,1);
	    sum_x2    += y_local(i,1)*y_local(i,1);
	    sum_x3    += y_local(i,1)*y_local(i,1)*y_local(i,1);
	    sum_x4    += y_local(i,1)*y_local(i,1)*y_local(i,1)*y_local(i,1);
	    sum_x5    += y_local(i,1)*y_local(i,1)*y_local(i,1)*y_local(i,1)*y_local(i,1);
	    sum_x6    += y_local(i,1)*y_local(i,1)*y_local(i,1)*y_local(i,1)*y_local(i,1)*y_local(i,1);

	    sum_curvature1     += curvature_picked(i,1);
	    sum_xcurvature1    += y_local(i,1)*curvature_picked(i,1);
	    sum_x2curvature1   += y_local(i,1)*y_local(i,1)*curvature_picked(i,1);
	    sum_x3curvature1   += y_local(i,1)*y_local(i,1)*y_local(i,1)*curvature_picked(i,1);

	    sum_curvature2     += curvature_picked(i,2);
	    sum_xcurvature2    += y_local(i,1)*curvature_picked(i,2);
	    sum_x2curvature2   += y_local(i,1)*y_local(i,1)*curvature_picked(i,2);
	    sum_x3curvature2   += y_local(i,1)*y_local(i,1)*y_local(i,1)*curvature_picked(i,2);

	    sum_metric_cov1    += metric_cov_picked(i,1);
	    sum_xmetric_cov1   += y_local(i,1)*metric_cov_picked(i,1);
	    sum_x2metric_cov1  += y_local(i,1)*y_local(i,1)*metric_cov_picked(i,1);

	    sum_metric_cov2    += metric_cov_picked(i,2);
	    sum_xmetric_cov2   += y_local(i,1)*metric_cov_picked(i,2);
	    sum_x2metric_cov2  += y_local(i,1)*y_local(i,1)*metric_cov_picked(i,2);

	}
	
	matrix temp_sum_x(4,4);
	temp_sum_x(1,1) = m;       temp_sum_x(1,2) = sum_x;  temp_sum_x(1,3) = sum_x2; temp_sum_x(1,4) = sum_x3;
	temp_sum_x(2,1) = sum_x;   temp_sum_x(2,2) = sum_x2; temp_sum_x(2,3) = sum_x3; temp_sum_x(2,4) = sum_x4;
	temp_sum_x(3,1) = sum_x2;  temp_sum_x(3,2) = sum_x3; temp_sum_x(3,3) = sum_x4; temp_sum_x(3,4) = sum_x5;
	temp_sum_x(4,1) = sum_x3;  temp_sum_x(4,2) = sum_x4; temp_sum_x(4,3) = sum_x5; temp_sum_x(4,4) = sum_x6;

	// a_coeff_curvature
	matrix temp_curv1(4,1);
	temp_curv1(1,1) = sum_curvature1;   temp_curv1(2,1) = sum_xcurvature1; 
	temp_curv1(3,1) = sum_x2curvature1; temp_curv1(4,1) = sum_x3curvature1;
	matrix a_coeff_curvature_1 = (temp_sum_x%temp_curv1).getTrans();

	matrix temp_curv2(4,1);
	temp_curv2(1,1) = sum_curvature2;   temp_curv2(2,1) = sum_xcurvature2; 
	temp_curv2(3,1) = sum_x2curvature2; temp_curv2(4,1) = sum_x3curvature2;
	matrix a_coeff_curvature_2 = (temp_sum_x%temp_curv2).getTrans();

	matrix a_coeff_curvature(1,8);
	a_coeff_curvature(1,1) = a_coeff_curvature_1(1,1); a_coeff_curvature(1,3) = a_coeff_curvature_1(1,2);
	a_coeff_curvature(1,5) = a_coeff_curvature_1(1,3); a_coeff_curvature(1,7) = a_coeff_curvature_1(1,4);

	a_coeff_curvature(1,2) = a_coeff_curvature_2(1,1); a_coeff_curvature(1,4) = a_coeff_curvature_2(1,2); 
	a_coeff_curvature(1,6) = a_coeff_curvature_2(1,3); a_coeff_curvature(1,8) = a_coeff_curvature_2(1,4);
	(*it)->setA_coeff_curvature(a_coeff_curvature);


	matrix temp_sum_x3(3,3);
	temp_sum_x3(1,1) = m;       temp_sum_x3(1,2) = sum_x;  temp_sum_x3(1,3) = sum_x2; 
	temp_sum_x3(2,1) = sum_x;   temp_sum_x3(2,2) = sum_x2; temp_sum_x3(2,3) = sum_x3; 
	temp_sum_x3(3,1) = sum_x2;  temp_sum_x3(3,2) = sum_x3; temp_sum_x3(3,3) = sum_x4; 


	// a_coeff_F
	matrix temp_metric1(3,1);
	temp_metric1(1,1) = sum_metric_cov1; temp_metric1(2,1) = sum_xmetric_cov1; temp_metric1(3,1) = sum_x2metric_cov1;
	matrix a_coeff_metric_cov_1 = (temp_sum_x3%temp_metric1).getTrans();

	matrix temp_metric2(3,1);
	temp_metric2(1,1) = sum_metric_cov2; temp_metric2(2,1) = sum_xmetric_cov2; temp_metric2(3,1) = sum_x2metric_cov2;
	matrix a_coeff_metric_cov_2 = (temp_sum_x3%temp_metric2).getTrans();

	matrix a_coeff_metric_cov(1,6);
	a_coeff_metric_cov(1,1) = a_coeff_metric_cov_1(1,1); a_coeff_metric_cov(1,3) = a_coeff_metric_cov_1(1,2);
	a_coeff_metric_cov(1,5) = a_coeff_metric_cov_1(1,3);
	a_coeff_metric_cov(1,2) = a_coeff_metric_cov_2(1,1); a_coeff_metric_cov(1,4) = a_coeff_metric_cov_2(1,2);
 	a_coeff_metric_cov(1,6) = a_coeff_metric_cov_2(1,3);
	(*it)->setA_coeff_metric_cov(a_coeff_metric_cov);

    } // Active_points for loop


} // find_laplace_beltrami_coeff()


void surfacetension::E_F_J_T_coeff_updating(){

    matrix Foot_points, n_Foot_points, E_Foot_points, F_Foot_points, T_Foot_points;
    matrix J_Foot_points(Active_points.size(), 1);

    int ii=1;
    for(std::list<activepoint*>::const_iterator it=Active_points.begin(); it!=Active_points.end(); ++it, ++ii){
	Foot_points.appendRow( ((*it)->getFoot_points()).getRow(1) );
	n_Foot_points.appendRow( ((*it)->getN_Foot_points()).getRow(1) );
	T_Foot_points.appendRow( ((*it)->getT_Foot_points()).getRow(1) );
 	E_Foot_points.appendRow( ((*it)->getE_Foot_points()).getRow(1) );
	F_Foot_points.appendRow( ((*it)->getF_Foot_points()).getRow(1) );
	J_Foot_points(ii,1) = (*it)->getJ_Foot_points();
    }


    if(plain_axi=="axi1"){

	    matrix temp_Foot_points = Foot_points;
	    Foot_points = zeros(2*temp_Foot_points.num_row,2);
	    for(int i=1; i<temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = temp_Foot_points(i,1);
		Foot_points(i,2) = temp_Foot_points(i,2);
	    }
	    for(int i=1+temp_Foot_points.num_row; i<2*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row,1))-0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-temp_Foot_points.num_row,2);
	    }

	    // n_Foot_points
	    matrix n_Foot_points_col1;
	    n_Foot_points_col1.appendCol(n_Foot_points.getCol(1));
	    n_Foot_points_col1 = -1*n_Foot_points_col1;
	    matrix sym_points2;
	    sym_points2.appendCol(n_Foot_points_col1.getCol(1));
 	    sym_points2.appendCol(n_Foot_points.getCol(2));
	    for(int i=1; i<sym_points2.num_row+1; ++i){
		n_Foot_points.appendRow(sym_points2.getRow(i));
	    }

	    // E_Foot_points
	    matrix temp_mat = E_Foot_points;
	    for(int i=1; i<temp_mat.num_row+1; ++i){
		E_Foot_points.appendRow(temp_mat.getRow(i));
	    }

	    // F_Foot_points
	    matrix temp_mat2 = F_Foot_points;
	    for(int i=1; i<temp_mat2.num_row+1; ++i){
		F_Foot_points.appendRow(temp_mat2.getRow(i));
	    }

	    // J_Foot_points
	    matrix temp_mat3 = J_Foot_points;
	    for(int i=1; i<temp_mat3.num_row+1; ++i){
		J_Foot_points.appendRow(temp_mat3.getRow(i));
	    }

	    // T_Foot_points
	    matrix temp_matT = T_Foot_points;
	    for(int i=1; i<temp_matT.num_row+1; ++i){
		T_Foot_points.appendRow(temp_matT.getRow(i));
	    }
	    
    } // if "axi1"

    if(plain_axi=="axi2"){
	    matrix temp_Foot_points = Foot_points;
	    Foot_points = zeros(3*temp_Foot_points.num_row,2);
	    for(int i=1; i<temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = temp_Foot_points(i,1);
		Foot_points(i,2) = temp_Foot_points(i,2);
	    }
	    for(int i=1+temp_Foot_points.num_row; i<2*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row,1))-0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-temp_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_Foot_points.num_row; i<3*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(-Lx*0.5+temp_Foot_points(i-2*temp_Foot_points.num_row,1))+0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-2*temp_Foot_points.num_row,2);
	    }

	    // n_Foot_points
	    matrix temp_n_Foot_points = n_Foot_points;
	    n_Foot_points = zeros(3*temp_n_Foot_points.num_row,2);
	    for(int i=1; i<temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = temp_n_Foot_points(i,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i,2);
	    }
	    for(int i=1+temp_n_Foot_points.num_row; i<2*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = -temp_n_Foot_points(i-temp_n_Foot_points.num_row,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i-temp_n_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_n_Foot_points.num_row; i<3*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = -temp_n_Foot_points(i-2*temp_n_Foot_points.num_row,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i-2*temp_n_Foot_points.num_row,2);
	    }

	    // E_Foot_points
	    matrix temp_E_Foot_points = E_Foot_points;
	    E_Foot_points = zeros(3*temp_E_Foot_points.num_row,2);
	    for(int i=1; i<temp_E_Foot_points.num_row+1; ++i){
		E_Foot_points(i,1) = temp_E_Foot_points(i,1);
		E_Foot_points(i,2) = temp_E_Foot_points(i,2);
	    }
	    for(int i=1+temp_E_Foot_points.num_row; i<2*temp_E_Foot_points.num_row+1; ++i){
		E_Foot_points(i,1) = temp_E_Foot_points(i-temp_E_Foot_points.num_row,1);
		E_Foot_points(i,2) = temp_E_Foot_points(i-temp_E_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_E_Foot_points.num_row; i<3*temp_E_Foot_points.num_row+1; ++i){
		E_Foot_points(i,1) = temp_E_Foot_points(i-2*temp_E_Foot_points.num_row,1);
		E_Foot_points(i,2) = temp_E_Foot_points(i-2*temp_E_Foot_points.num_row,2);
	    }

	    // F_Foot_points
	    matrix temp_F_Foot_points = F_Foot_points;
	    F_Foot_points = zeros(3*temp_F_Foot_points.num_row,2);
	    for(int i=1; i<temp_F_Foot_points.num_row+1; ++i){
		F_Foot_points(i,1) = temp_F_Foot_points(i,1);
		F_Foot_points(i,2) = temp_F_Foot_points(i,2);
	    }
	    for(int i=1+temp_F_Foot_points.num_row; i<2*temp_F_Foot_points.num_row+1; ++i){
		F_Foot_points(i,1) = temp_F_Foot_points(i-temp_F_Foot_points.num_row,1);
		F_Foot_points(i,2) = temp_F_Foot_points(i-temp_F_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_F_Foot_points.num_row; i<3*temp_F_Foot_points.num_row+1; ++i){
		F_Foot_points(i,1) = temp_F_Foot_points(i-2*temp_F_Foot_points.num_row,1);
		F_Foot_points(i,2) = temp_F_Foot_points(i-2*temp_F_Foot_points.num_row,2);
	    }

	    // J_Foot_points
	    matrix temp_J_Foot_points = J_Foot_points;
	    J_Foot_points = zeros(3*temp_J_Foot_points.num_row,1);
	    for(int i=1; i<temp_J_Foot_points.num_row+1; ++i){
		J_Foot_points(i,1) = temp_J_Foot_points(i,1);
	    }
	    for(int i=1+temp_J_Foot_points.num_row; i<2*temp_J_Foot_points.num_row+1; ++i){
		J_Foot_points(i,1) = temp_J_Foot_points(i-temp_J_Foot_points.num_row,1);
	    }
	    for(int i=1+2*temp_J_Foot_points.num_row; i<3*temp_J_Foot_points.num_row+1; ++i){
		J_Foot_points(i,1) = temp_J_Foot_points(i-2*temp_J_Foot_points.num_row,1);
	    }

	    // T_Foot_points
	    matrix temp_T_Foot_points = T_Foot_points;
	    T_Foot_points = zeros(3*temp_T_Foot_points.num_row,2);
	    for(int i=1; i<temp_T_Foot_points.num_row+1; ++i){
		T_Foot_points(i,1) = temp_T_Foot_points(i,1);
		T_Foot_points(i,2) = temp_T_Foot_points(i,2);
	    }
	    for(int i=1+temp_T_Foot_points.num_row; i<2*temp_T_Foot_points.num_row+1; ++i){
		T_Foot_points(i,1) = temp_T_Foot_points(i-temp_T_Foot_points.num_row,1);
		T_Foot_points(i,2) = temp_T_Foot_points(i-temp_T_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_T_Foot_points.num_row; i<3*temp_T_Foot_points.num_row+1; ++i){
		T_Foot_points(i,1) = temp_T_Foot_points(i-2*temp_T_Foot_points.num_row,1);
		T_Foot_points(i,2) = temp_T_Foot_points(i-2*temp_T_Foot_points.num_row,2);
	    }

    } // if "axi2"

    if(plain_axi=="plain"){
	
    } // if "plain"


    for(std::list<activepoint*>::iterator it=Active_points.begin(); it!=Active_points.end(); ++it){

	(*it)->initCoeff_E_F_T_J();

	matrix p = (*it)->getFoot_points();

//	matrix D = zeros(Foot_points.num_row, 1);
	std::list<REAL> D;
	for(int i=1; i<Foot_points.num_row+1; ++i){
	    REAL D_val = pow( (Foot_points(i,1)-p(1,1))*(Foot_points(i,1)-p(1,1))
			    + (Foot_points(i,2)-p(1,2))*(Foot_points(i,2)-p(1,2)), 0.5);
	    D.push_back(D_val);
	}

	std::vector<int> IDX = sort30(D.begin(), D.end());

	matrix y_sorted = zeros(30,2);
   	matrix T_sorted = zeros(30,2);
	matrix E_sorted = zeros(30,2);
	matrix F_sorted = zeros(30,2);
	matrix J_sorted = zeros(30,1);
	matrix n_sorted = zeros(30,2);

	for(int i=0; i<30; ++i){
	    y_sorted(i+1, 1) = Foot_points(IDX[i], 1);
	    y_sorted(i+1, 2) = Foot_points(IDX[i], 2);

	    T_sorted(i+1, 1) = T_Foot_points(IDX[i], 1);
	    T_sorted(i+1, 2) = T_Foot_points(IDX[i], 2);

	    E_sorted(i+1, 1) = E_Foot_points(IDX[i], 1);
	    E_sorted(i+1, 2) = E_Foot_points(IDX[i], 2);

	    F_sorted(i+1, 1) = F_Foot_points(IDX[i], 1);
	    F_sorted(i+1, 2) = F_Foot_points(IDX[i], 2);

	    J_sorted(i+1, 1) = J_Foot_points(IDX[i], 1);

	    n_sorted(i+1, 1) = n_Foot_points(IDX[i], 1);
	    n_sorted(i+1, 2) = n_Foot_points(IDX[i], 2);
	}

	//n0_old and y0 are the origin and orientation of the local coordinate system

	matrix n0_old(1,2);
	n0_old(1,1) = n_Foot_points(IDX[0], 1);
	n0_old(1,2) = n_Foot_points(IDX[0], 2);

	matrix y0(1,2);
	y0(1,1) = y_sorted(1,1);
	y0(1,2) = y_sorted(1,2);

	matrix E0(1,2);
	E0(1,1) = E_sorted(1,1);
	E0(1,2) = E_sorted(1,2);

	matrix F0(1,2);
	F0(1,1) = F_sorted(1,1);
	F0(1,2) = F_sorted(1,2);

	matrix T0(1,2);
	T0(1,1) = T_sorted(1,1);
	T0(1,2) = T_sorted(1,2);

	REAL J0 = J_sorted(1,1);

	matrix n0(1,2);
	n0(1,1) = n_sorted(1,1);
	n0(1,2) = n_sorted(1,2);

	//we now select a few neighboring Foot_points to construct the 
	//polynomials that interpolate the interface

      	matrix y_picked = y0;
    	matrix E_picked = E0;
    	matrix F_picked = F0;
    	matrix T_picked = T0;
    	matrix J_picked(1,1); J_picked(1,1) = J0;
    	matrix n_picked = n0;

	for(int i=2; i<length(y_sorted)+1; ++i){
	    REAL d;
	    for(int j=1; j<y_picked.num_row+1; ++j){
	    	REAL d_temp = pow( (y_picked(j,1)-y_sorted(i,1))*(y_picked(j,1)-y_sorted(i,1))
			         + (y_picked(j,2)-y_sorted(i,2))*(y_picked(j,2)-y_sorted(i,2)), 0.5);

		if(j==1){
		    d = d_temp;
		}
		else{
		    if(d>d_temp)
			d = d_temp;
		}
	    }

	    if(d >= 0.2*element_length){
		y_picked.appendRow(y_sorted.getRow(i));
		E_picked.appendRow(E_sorted.getRow(i));
		F_picked.appendRow(F_sorted.getRow(i));
		T_picked.appendRow(T_sorted.getRow(i));
		J_picked.appendRow(J_sorted.getRow(i));
	    }

	}

//	matrix d_sort = y_picked;
	std::list<REAL> d_sort;
	for(int i=1; i<y_picked.num_row+1; ++i){
	    REAL d_sort_val = pow( (y_picked(i,1)-y0(1,1))*(y_picked(i,1)-y0(1,1))
			         + (y_picked(i,2)-y0(1,2))*(y_picked(i,2)-y0(1,2)), 0.5);
	    d_sort.push_back(d_sort_val);
	}

	std::vector<int> i_sort = sort30(d_sort.begin(), d_sort.end());

	matrix y_picked_temp = y_picked;
	matrix E_picked_temp = E_picked;
	matrix F_picked_temp = F_picked;
	matrix T_picked_temp = T_picked;
	matrix J_picked_temp = J_picked;
	for(int i=0; i<y_picked_temp.num_row; ++i){
	    y_picked(i+1, 1) = y_picked_temp(i_sort[i], 1);
	    y_picked(i+1, 2) = y_picked_temp(i_sort[i], 2);

	    E_picked(i+1, 1) = E_picked_temp(i_sort[i], 1);
	    E_picked(i+1, 2) = E_picked_temp(i_sort[i], 2);

	    F_picked(i+1, 1) = F_picked_temp(i_sort[i], 1);
	    F_picked(i+1, 2) = F_picked_temp(i_sort[i], 2);

	    T_picked(i+1, 1) = T_picked_temp(i_sort[i], 1);
	    T_picked(i+1, 2) = T_picked_temp(i_sort[i], 2);

	    J_picked(i+1, 1) = J_picked_temp(i_sort[i], 1);

	}

	int m = length(y_picked);

	if(m<4){
	    std::cout << "not enough sampled points to build order 3 polynomial" << std::endl;
//	    exit(-1);
	    continue;
	}

	if(m>=9){
	    m = 9;

	    matrix y_picked_temp = y_picked;
	    matrix E_picked_temp = E_picked;
	    matrix F_picked_temp = F_picked;
	    matrix T_picked_temp = T_picked;
	    matrix J_picked_temp = J_picked;
	    y_picked = zeros(9,2);
	    T_picked = zeros(9,2);
	    E_picked = zeros(9,2);
	    F_picked = zeros(9,2);
	    J_picked = zeros(9,1);

	    for(int i=1; i<10; ++i){
	    	y_picked(i,1) = y_picked_temp(i,1);
		y_picked(i,2) = y_picked_temp(i,2);

	    	E_picked(i,1) = E_picked_temp(i,1);
		E_picked(i,2) = E_picked_temp(i,2);

	    	F_picked(i,1) = F_picked_temp(i,1);
		F_picked(i,2) = F_picked_temp(i,2);

	    	T_picked(i,1) = T_picked_temp(i,1);
		T_picked(i,2) = T_picked_temp(i,2);

	    	J_picked(i,1) = J_picked_temp(i,1);

	    }
	}

	//find coordinates of picked points in the local sytem defined by n0_old
	matrix y_local;
	for(int i=1; i<length(y_picked)+1; ++i){
	    matrix temp_local;
	    matrix temp_n0(2,2);
	    temp_n0(1,1) = n0_old(1,2); temp_n0(1,2) = -n0_old(1,1);
	    temp_n0(2,1) = n0_old(1,1); temp_n0(2,2) = n0_old(1,2);
	    matrix temp_y(2,1);
	    temp_y(1,1) = y_picked(i,1)-y0(1,1); temp_y(2,1) = y_picked(i,2)-y0(1,2);
	    temp_local = temp_n0*temp_y;
	    
	    y_local.appendRow(temp_local.getCol(1));
  	}

	REAL xmin = y_local(1,1);
	REAL xmax = y_local(1,1);
	for(int i=2; i<y_local.num_row+1; ++i){
	    if(xmin > y_local(i,1))
		xmin = y_local(i,1);
	    if(xmax < y_local(i,1))
		xmax = y_local(i,1);
	}

	// build the polynomial that locally parameterizes the membrane (order 2)
	REAL sum_x = 0; REAL sum_x2 = 0; REAL sum_x3 = 0; REAL sum_x4 = 0;
	REAL sum_y = 0; REAL sum_xy = 0; REAL sum_x2y = 0;
	REAL sum_E1 = 0; REAL sum_xE1 = 0; REAL sum_x2E1 = 0;
	REAL sum_E2 = 0; REAL sum_xE2 = 0; REAL sum_x2E2 = 0;
	REAL sum_F1 = 0; REAL sum_xF1 = 0; REAL sum_x2F1 = 0;
	REAL sum_F2 = 0; REAL sum_xF2 = 0; REAL sum_x2F2 = 0; 
	REAL sum_T1 = 0; REAL sum_xT1 = 0; REAL sum_x2T1 = 0;
	REAL sum_T2 = 0; REAL sum_xT2 = 0; REAL sum_x2T2 = 0; 
	REAL sum_J = 0; REAL sum_xJ = 0; REAL sum_x2J= 0;
	for(int i=1; i<y_local.num_row+1; ++i){
	    sum_x     += y_local(i,1);
	    sum_x2    += y_local(i,1)*y_local(i,1);
	    sum_x3    += y_local(i,1)*y_local(i,1)*y_local(i,1);
	    sum_x4    += y_local(i,1)*y_local(i,1)*y_local(i,1)*y_local(i,1);

	    sum_y     += y_local(i,2);
	    sum_xy    += y_local(i,1)*y_local(i,2);
	    sum_x2y   += y_local(i,1)*y_local(i,1)*y_local(i,2);

	    sum_E1    += E_picked(i,1);
	    sum_xE1   += y_local(i,1)*E_picked(i,1);
	    sum_x2E1  += y_local(i,1)*y_local(i,1)*E_picked(i,1);

	    sum_E2    += E_picked(i,2);
	    sum_xE2   += y_local(i,1)*E_picked(i,2);
	    sum_x2E2  += y_local(i,1)*y_local(i,1)*E_picked(i,2);

	    sum_F1    += F_picked(i,1);
	    sum_xF1   += y_local(i,1)*F_picked(i,1);
	    sum_x2F1  += y_local(i,1)*y_local(i,1)*F_picked(i,1);

	    sum_F2    += F_picked(i,2);
	    sum_xF2   += y_local(i,1)*F_picked(i,2);
	    sum_x2F2  += y_local(i,1)*y_local(i,1)*F_picked(i,2);

	    sum_T1    += T_picked(i,1);
	    sum_xT1   += y_local(i,1)*T_picked(i,1);
	    sum_x2T1  += y_local(i,1)*y_local(i,1)*T_picked(i,1);

	    sum_T2    += T_picked(i,2);
	    sum_xT2   += y_local(i,1)*T_picked(i,2);
	    sum_x2T2  += y_local(i,1)*y_local(i,1)*T_picked(i,2);

	    sum_J     += J_picked(i,1);
	    sum_xJ    += y_local(i,1)*J_picked(i,1);
	    sum_x2J   += y_local(i,1)*y_local(i,1)*J_picked(i,1);

	}
	
	matrix temp_sum_xs(3,3);
	temp_sum_xs(1,1) = m;      temp_sum_xs(1,2) = sum_x;  temp_sum_xs(1,3) = sum_x2;
	temp_sum_xs(2,1) = sum_x;  temp_sum_xs(2,2) = sum_x2; temp_sum_xs(2,3) = sum_x3;
	temp_sum_xs(3,1) = sum_x2; temp_sum_xs(3,2) = sum_x3; temp_sum_xs(3,3) = sum_x4;

	// a_coeff
//	matrix temp_xy(3,1);
//	temp_xy(1,1) = sum_y; temp_xy(2,1) = sum_xy; temp_xy(3,1) = sum_x2y;
//	matrix a_coeff = (temp_sum_xs%temp_xy).getTrans();
//	(*it)->setA_coeff(a_coeff);

	// a_coeff_E
	matrix temp_xE1(3,1);
	temp_xE1(1,1) = sum_E1; temp_xE1(2,1) = sum_xE1; temp_xE1(3,1) = sum_x2E1;
	matrix a_coeff_E_1 = (temp_sum_xs%temp_xE1).getTrans();

	matrix temp_xE2(3,1);
	temp_xE2(1,1) = sum_E2; temp_xE2(2,1) = sum_xE2; temp_xE2(3,1) = sum_x2E2;
	matrix a_coeff_E_2 = (temp_sum_xs%temp_xE2).getTrans();

	matrix a_coeff_E(1,6);
	a_coeff_E(1,1) = a_coeff_E_1(1,1); a_coeff_E(1,3) = a_coeff_E_1(1,2); a_coeff_E(1,5) = a_coeff_E_1(1,3);
	a_coeff_E(1,2) = a_coeff_E_2(1,1); a_coeff_E(1,4) = a_coeff_E_2(1,2); a_coeff_E(1,6) = a_coeff_E_2(1,3);
	(*it)->setA_coeff_E(a_coeff_E);

	// a_coeff_F
	matrix temp_xF1(3,1);
	temp_xF1(1,1) = sum_F1; temp_xF1(2,1) = sum_xF1; temp_xF1(3,1) = sum_x2F1;
	matrix a_coeff_F_1 = (temp_sum_xs%temp_xF1).getTrans();

	matrix temp_xF2(3,1);
	temp_xF2(1,1) = sum_F2; temp_xF2(2,1) = sum_xF2; temp_xF2(3,1) = sum_x2F2;
	matrix a_coeff_F_2 = (temp_sum_xs%temp_xF2).getTrans();

	matrix a_coeff_F(1,6);
	a_coeff_F(1,1) = a_coeff_F_1(1,1); a_coeff_F(1,3) = a_coeff_F_1(1,2); a_coeff_F(1,5) = a_coeff_F_1(1,3);
	a_coeff_F(1,2) = a_coeff_F_2(1,1); a_coeff_F(1,4) = a_coeff_F_2(1,2); a_coeff_F(1,6) = a_coeff_F_2(1,3);
	(*it)->setA_coeff_F(a_coeff_F);

	// a_coeff_T
	matrix temp_xT1(3,1);
	temp_xT1(1,1) = sum_T1; temp_xT1(2,1) = sum_xT1; temp_xT1(3,1) = sum_x2T1;
	matrix a_coeff_T_1 = (temp_sum_xs%temp_xT1).getTrans();

	matrix temp_xT2(3,1);
	temp_xT2(1,1) = sum_T2; temp_xT2(2,1) = sum_xT2; temp_xT2(3,1) = sum_x2T2;
	matrix a_coeff_T_2 = (temp_sum_xs%temp_xT2).getTrans();

	matrix a_coeff_T(1,6);
	a_coeff_T(1,1) = a_coeff_T_1(1,1); a_coeff_T(1,3) = a_coeff_T_1(1,2); a_coeff_T(1,5) = a_coeff_T_1(1,3);
	a_coeff_T(1,2) = a_coeff_T_2(1,1); a_coeff_T(1,4) = a_coeff_T_2(1,2); a_coeff_T(1,6) = a_coeff_T_2(1,3);
	(*it)->setA_coeff_T(a_coeff_T);

	// a_coeff_J
	matrix temp_xJ(3,1);
	temp_xJ(1,1) = sum_J; temp_xJ(2,1) = sum_xJ; temp_xJ(3,1) = sum_x2J;
	matrix a_coeff_J = (temp_sum_xs%temp_xJ).getTrans();
	(*it)->setA_coeff_J(a_coeff_J);

    } // Active_points for loop


} // E_F_J_T_coeff_updating()


void surfacetension::find_rotation(){

    // vel_membrane_update_normal is vn_Foot_points         1x2
    // vel_membrane_update_tangent is vt_Foot_points        REAL
    matrix Foot_points, n_Foot_points, vn_Foot_points;

    for(std::list<activepoint*>::const_iterator it=Active_points.begin(); it!=Active_points.end(); ++it){
	Foot_points.appendRow( ((*it)->getFoot_points()).getRow(1) );
	n_Foot_points.appendRow( ((*it)->getN_Foot_points()).getRow(1) );
	vn_Foot_points.appendRow( ((*it)->getVel_membrane_update_normal()).getRow(1) );
// 	vt_Foot_points.appendRow( ((*it)->getVel_membrane_update_tangent()).getRow(1) );
    }


    if(plain_axi=="axi1"){
	    matrix temp_Foot_points = Foot_points;
	    Foot_points = zeros(2*temp_Foot_points.num_row,2);
	    for(int i=1; i<temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = temp_Foot_points(i,1);
		Foot_points(i,2) = temp_Foot_points(i,2);
	    }
	    for(int i=1+temp_Foot_points.num_row; i<2*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row,1))-0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-temp_Foot_points.num_row,2);
	    }

	    // n_Foot_points
	    matrix n_Foot_points_col1;
	    n_Foot_points_col1.appendCol(n_Foot_points.getCol(1));
	    n_Foot_points_col1 = -1*n_Foot_points_col1;
	    matrix sym_points2;
	    sym_points2.appendCol(n_Foot_points_col1.getCol(1));
 	    sym_points2.appendCol(n_Foot_points.getCol(2));
	    for(int i=1; i<sym_points2.num_row+1; ++i){
		n_Foot_points.appendRow(sym_points2.getRow(i));
	    }

	    // vn_Foot_points
	    matrix temp_mat = vn_Foot_points;
	    for(int i=1; i<temp_mat.num_row+1; ++i){
		vn_Foot_points.appendRow(temp_mat.getRow(i));
	    }
	    
    } // if "axi1"

    if(plain_axi=="axi2"){
	    matrix temp_Foot_points = Foot_points;
	    Foot_points = zeros(3*temp_Foot_points.num_row,2);
	    for(int i=1; i<temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = temp_Foot_points(i,1);
		Foot_points(i,2) = temp_Foot_points(i,2);
	    }
	    for(int i=1+temp_Foot_points.num_row; i<2*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row,1))-0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-temp_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_Foot_points.num_row; i<3*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(-Lx*0.5+temp_Foot_points(i-2*temp_Foot_points.num_row,1))+0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-2*temp_Foot_points.num_row,2);
	    }

	    // n_Foot_points
	    matrix temp_n_Foot_points = n_Foot_points;
	    n_Foot_points = zeros(3*temp_n_Foot_points.num_row,2);
	    for(int i=1; i<temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = temp_n_Foot_points(i,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i,2);
	    }
	    for(int i=1+temp_n_Foot_points.num_row; i<2*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = -temp_n_Foot_points(i-temp_n_Foot_points.num_row,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i-temp_n_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_n_Foot_points.num_row; i<3*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = -temp_n_Foot_points(i-2*temp_n_Foot_points.num_row,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i-2*temp_n_Foot_points.num_row,2);
	    }

	    // vn_Foot_points
	    matrix temp_vn_Foot_points = vn_Foot_points;
	    vn_Foot_points = zeros(3*temp_vn_Foot_points.num_row,2);
	    for(int i=1; i<temp_vn_Foot_points.num_row+1; ++i){
		vn_Foot_points(i,1) = temp_vn_Foot_points(i,1);
		vn_Foot_points(i,2) = temp_vn_Foot_points(i,2);
	    }
	    for(int i=1+temp_vn_Foot_points.num_row; i<2*temp_vn_Foot_points.num_row+1; ++i){
		vn_Foot_points(i,1) = temp_vn_Foot_points(i-temp_vn_Foot_points.num_row,1);
		vn_Foot_points(i,2) = temp_vn_Foot_points(i-temp_vn_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_vn_Foot_points.num_row; i<3*temp_vn_Foot_points.num_row+1; ++i){
		vn_Foot_points(i,1) = temp_vn_Foot_points(i-2*temp_vn_Foot_points.num_row,1);
		vn_Foot_points(i,2) = temp_vn_Foot_points(i-2*temp_vn_Foot_points.num_row,2);
	    }

    } // if "axi2"

    if(plain_axi=="plain"){
	
    } // if "plain"


    for(std::list<activepoint*>::iterator it=Active_points.begin(); it!=Active_points.end(); ++it){

	matrix a_coeff_vn(1,3);

	matrix p = (*it)->getFoot_points();

//	matrix D2 = zeros(Foot_points.num_row, 1);
	std::list<REAL> D2;
	for(int i=1; i<Foot_points.num_row+1; ++i){
	    REAL D2_val = pow( (Foot_points(i,1)-p(1,1))*(Foot_points(i,1)-p(1,1))
		 	     + (Foot_points(i,2)-p(1,2))*(Foot_points(i,2)-p(1,2)), 0.5);
	    D2.push_back(D2_val);
	}

	std::vector<int> IDX = sort30(D2.begin(), D2.end());

	matrix y_sorted = zeros(30,2);
   	matrix vn_sorted = zeros(30,2);

	for(int i=0; i<30; ++i){
	    y_sorted(i+1, 1) = Foot_points(IDX[i], 1);
	    y_sorted(i+1, 2) = Foot_points(IDX[i], 2);

	    vn_sorted(i+1, 1) = vn_Foot_points(IDX[i], 1);
	    vn_sorted(i+1, 2) = vn_Foot_points(IDX[i], 2);
	}

	//n0_old and y0 are the origin and orientation of the local coordinate system

	matrix n0_old(1,2);
	n0_old(1,1) = n_Foot_points(IDX[0], 1);
	n0_old(1,2) = n_Foot_points(IDX[0], 2);

	matrix y0(1,2);
	y0(1,1) = y_sorted(1,1);
	y0(1,2) = y_sorted(1,2);

	matrix vn0(1,2);
	vn0(1,1) = vn_sorted(1,1);
	vn0(1,2) = vn_sorted(1,2);

	//we now select a few neighboring Foot_points to construct the 
	//polynomials that interpolate the interface

      	matrix y_picked = y0;
    	matrix vn_picked = vn0;

	for(int i=2; i<length(y_sorted)+1; ++i){
	    REAL d;
	    for(int j=1; j<y_picked.num_row+1; ++j){
	    	REAL d_temp = pow( (y_picked(j,1)-y_sorted(i,1))*(y_picked(j,1)-y_sorted(i,1))
			         + (y_picked(j,2)-y_sorted(i,2))*(y_picked(j,2)-y_sorted(i,2)), 0.5);

		if(j==1){
		    d = d_temp;
		}
		else{
		    if(d>d_temp)
			d = d_temp;
		}
	    }

	    if(d >= 0.2*element_length){
		y_picked.appendRow(y_sorted.getRow(i));
		vn_picked.appendRow(vn_sorted.getRow(i));
	    }

	}

	int m = length(y_picked);

	if(m<4){
	    std::cout << "not enough sampled points to build order 3 polynomial" << std::endl;
//	    exit(-1);
	    continue;
	}

	if(m>=9){
	    m = 9;

	    matrix y_picked_temp = y_picked;
	    matrix vn_picked_temp = vn_picked;
	    y_picked = zeros(9,2);
	    vn_picked = zeros(9,2);

	    for(int i=1; i<10; ++i){
	    	y_picked(i,1) = y_picked_temp(i,1);
		y_picked(i,2) = y_picked_temp(i,2);

	    	vn_picked(i,1) = vn_picked_temp(i,1);
		vn_picked(i,2) = vn_picked_temp(i,2);
	    }
	}

	//find coordinates of picked points in the local sytem defined by n0_old
	matrix y_local;
	for(int i=1; i<length(y_picked)+1; ++i){
	    matrix temp_local;
	    matrix temp_n0(2,2);
	    temp_n0(1,1) = n0_old(1,2); temp_n0(1,2) = -n0_old(1,1);
	    temp_n0(2,1) = n0_old(1,1); temp_n0(2,2) = n0_old(1,2);
	    matrix temp_y(2,1);
	    temp_y(1,1) = y_picked(i,1)-y0(1,1); temp_y(2,1) = y_picked(i,2)-y0(1,2);
	    temp_local = temp_n0*temp_y;
	    
	    y_local.appendRow(temp_local.getCol(1));
  	}

	REAL xmin = y_local(1,1);
	REAL xmax = y_local(1,1);
	for(int i=2; i<y_local.num_row+1; ++i){
	    if(xmin > y_local(i,1))
		xmin = y_local(i,1);
	    if(xmax < y_local(i,1))
		xmax = y_local(i,1);
	}

	// build the polynomial that locally parameterizes the membrane (order 2)
	REAL sum_x = 0; REAL sum_x2 = 0; REAL sum_x3 = 0; REAL sum_x4 = 0;
	REAL sum_vn1 = 0; REAL sum_xvn1 = 0; REAL sum_x2vn1 = 0;

	for(int i=1; i<y_local.num_row+1; ++i){
	    sum_x      += y_local(i,1);
	    sum_x2     += y_local(i,1)*y_local(i,1);
	    sum_x3     += y_local(i,1)*y_local(i,1)*y_local(i,1);
	    sum_x4     += y_local(i,1)*y_local(i,1)*y_local(i,1)*y_local(i,1);

	    sum_vn1    += vn_picked(i,1);
	    sum_xvn1   += y_local(i,1)*vn_picked(i,1);
	    sum_x2vn1  += y_local(i,1)*y_local(i,1)*vn_picked(i,1);

	}
		
	REAL a = sqrt( 1+((*it)->getA_coeff3())(1,2) );

 	REAL curvature_cov = (2*((*it)->getA_coeff3())(1,3))/a;

	matrix temp_sum_xs(3,3);
	temp_sum_xs(1,1) = m;      temp_sum_xs(1,2) = sum_x;  temp_sum_xs(1,3) = sum_x2;
	temp_sum_xs(2,1) = sum_x;  temp_sum_xs(2,2) = sum_x2; temp_sum_xs(2,3) = sum_x3;
	temp_sum_xs(3,1) = sum_x2; temp_sum_xs(3,2) = sum_x3; temp_sum_xs(3,3) = sum_x4;

	// a_coeff_vn
	matrix temp_xvn1(3,1);
	temp_xvn1(1,1) = sum_vn1; temp_xvn1(2,1) = sum_xvn1; temp_xvn1(3,1) = sum_x2vn1;
	a_coeff_vn = (temp_sum_xs%temp_xvn1).getTrans();

	REAL Omega_e = a_coeff_vn(1,2) + ((*it)->getVel_membrane_update_tangent())*curvature_cov;

 	(*it)->setOmega_e(Omega_e);
    } // Active_points for loop

} // find_rotation()


void surfacetension::sort_points(){

    std::map<array, int> points;
    std::vector<ele2dcoupled*>::const_iterator it = split_elem2.begin();
    for(int i=1; i<2*split_elem2.size()+1; i+=2, ++it){
 	REAL x_tot1 = round( ((*it)->getXc2())(1,1)*1e10);
	REAL x_tot2 = round( ((*it)->getYc2())(1,1)*1e10);
	array x_tot(x_tot1, x_tot2);
	points[x_tot]++;

    }

    it = split_elem2.begin();
    for(int i=2; i<2*split_elem2.size()+1; i+=2, ++it){
 	REAL x_tot1 = round( ((*it)->getXc2())(1,2)*1e10) ;
	REAL x_tot2 = round( ((*it)->getYc2())(1,2)*1e10);
	array x_tot(x_tot1, x_tot2);
	points[x_tot]++;

    }

    // implemented: before points = 1e-10*unique(x_tot, 'rows');

    std::list<REAL> points_x;
    std::list<REAL> points_y;

    for(std::map<array, int>::iterator it_pt=points.begin(); it_pt!=points.end(); ++it_pt){
	points_x.push_back( (it_pt->first.getX())*1e-10 );
	points_y.push_back( (it_pt->first.getY())*1e-10 );
    }


    std::vector<REAL> Ordered_pointx;
    std::vector<REAL> Ordered_pointy;

    if(plain_axi=="axi1"){
	    std::list<REAL>::iterator ity = points_y.begin();
	    int i = 0;
	    REAL xmin;
	    std::list<REAL>::iterator ity_min, itx_min;
	    for(std::list<REAL>::iterator itx=points_x.begin(); itx!=points_x.end(); ++itx, ++ity){
		if((*itx)==-0.5*Lx){	// find(points(:,1)==-Lx/2)
		    if(i==0){
			xmin = (*itx);
			itx_min = itx;
			ity_min = ity;
		    }
		    else{
			if(xmin>(*itx)){
			    xmin = (*itx);
			    itx_min = itx;
			    ity_min = ity;
			}
		    }

		    i++;
		} // end if -0.5*Lx
 	    } // end for

	    Ordered_pointx.push_back((*itx_min));
	    Ordered_pointy.push_back((*ity_min));

	    points_x.erase(itx_min);
	    points_y.erase(ity_min);

	    while(!points_x.empty()){
		std::vector<REAL>::iterator px = Ordered_pointx.end();
	    	px--;
		std::vector<REAL>::iterator py = Ordered_pointy.end();
	 	py--;
		
		std::list<REAL>::iterator IDX2x, IDX2y;
		ity = points_y.begin();
		REAL dmin;
		for(std::list<REAL>::iterator itx=points_x.begin(); itx!=points_x.end(); ++itx, ++ity){
		    REAL d = sqrt( ((*itx)-(*px))*((*itx)-(*px))+((*ity)-(*py))*((*ity)-(*py)) );
		    if(itx==points_x.begin()){
			dmin = d;
			IDX2x = itx;
			IDX2y = ity;
		    }
		    else{
			if(dmin>d){
			    dmin = d;
			    IDX2x = itx;
			    IDX2y = ity;
			}
		    }
		} // end for

		Ordered_pointx.push_back((*IDX2x));
		Ordered_pointy.push_back((*IDX2y));

		points_x.erase(IDX2x);
		points_y.erase(IDX2y);

	    } // end while

    } // if "axi1"

    if(plain_axi=="axi2"){
	    std::list<REAL>::iterator ity = points_y.begin();
	    int i = 0;
	    REAL xmin;
	    std::list<REAL>::iterator ity_min, itx_min;
	    for(std::list<REAL>::iterator itx=points_x.begin(); itx!=points_x.end(); ++itx, ++ity){
		if((*itx)==-0.5*Lx){	// find(points(:,1)==-Lx/2)
		    if(i==0){
			xmin = (*itx);
			itx_min = itx;
			ity_min = ity;
		    }
		    else{
			if(xmin>(*itx)){
			    xmin = (*itx);
			    itx_min = itx;
			    ity_min = ity;
			}
		    }

		    i++;
		} // end if -0.5*Lx
 	    } // end for

	    Ordered_pointx.push_back((*itx_min));
	    Ordered_pointy.push_back((*ity_min));

	    points_x.erase(itx_min);
	    points_y.erase(ity_min);

	    while(!points_x.empty()){
		std::vector<REAL>::iterator px = Ordered_pointx.end();
	    	px--;
		std::vector<REAL>::iterator py = Ordered_pointy.end();
	 	py--;
		
		std::list<REAL>::iterator IDX2x, IDX2y;
		ity = points_y.begin();
		REAL dmin;
		for(std::list<REAL>::iterator itx=points_x.begin(); itx!=points_x.end(); ++itx, ++ity){
		    REAL d = sqrt( ((*itx)-(*px))*((*itx)-(*px))+((*ity)-(*py))*((*ity)-(*py)) );
		    if(itx==points_x.begin()){
			dmin = d;
			IDX2x = itx;
			IDX2y = ity;
		    }
		    else{
			if(dmin>d){
			    dmin = d;
			    IDX2x = itx;
			    IDX2y = ity;
			}
		    }
		} // end for

		Ordered_pointx.push_back((*IDX2x));
		Ordered_pointy.push_back((*IDX2y));

		points_x.erase(IDX2x);
		points_y.erase(IDX2y);

	    } // end while

    } // if "axi2"
	
    if(plain_axi=="plain"){
	    std::list<REAL>::iterator itx = points_x.begin();
	    std::list<REAL>::iterator ity = points_y.begin();

	    Ordered_pointx.push_back(*itx);
	    Ordered_pointy.push_back(*ity);

	    points_x.erase(itx);
	    points_y.erase(ity);

	    int np = 0;
	    while(!points_x.empty()){
		np = np+1;
		std::vector<REAL>::iterator px = Ordered_pointx.end();
	    	px--;
		std::vector<REAL>::iterator py = Ordered_pointy.end();
	 	py--;
		std::vector<REAL> D2;
		ity = points_y.begin();
		for(itx=points_x.begin(); itx!=points_x.end(); ++itx, ++ity){
		    REAL d = sqrt( ((*itx)-(*px))*((*itx)-(*px))+((*ity)-(*py))*((*ity)-(*py)) );
		    D2.push_back(d);
		}
		std::sort(D2.begin(), D2.end());
		std::vector<REAL>::iterator itD2_1=D2.begin();
		std::vector<REAL>::iterator itD2_2=itD2_1+1;
		std::vector<REAL>::iterator itD2_3=itD2_1+2;

		std::list<REAL>::iterator IDX2x1, IDX2x2, IDX2x3;	// [~, IDX2] = sort(D2);
		std::list<REAL>::iterator IDX2y1, IDX2y2, IDX2y3;	// the first three smallest D2
		
		int num=0;	// if num==3, means we find three smallest D2
		ity=points_y.begin();
		for(itx=points_x.begin(); itx!=points_x.end(); ++itx, ++ity){
		    REAL d = sqrt( ((*itx)-(*px))*((*itx)-(*px))+((*ity)-(*py))*((*ity)-(*py)) );
		    if(fabs(d-(*itD2_1)) < 1e-12){
			IDX2x1 = itx;
			IDX2y1 = ity;
			num++;
		    }
		    if(fabs(d-(*itD2_2)) < 1e-12){
			IDX2x2 = itx;
			IDX2y2 = ity;
			num++;
		    }
		    if(fabs(d-(*itD2_3)) < 1e-12){
			IDX2x3 = itx;
			IDX2y3 = ity;
			num++;
		    }
		}
		
		if(num!=3){
		    std::cout << "number in sort_points() in surfacetension.cpp is not 3: " << num << std::endl;
		    exit(-1);
		}

		matrix direction=zeros(1,2);	// initialized in else
		std::vector<REAL>::iterator npx = Ordered_pointx.end(); npx--;
		std::vector<REAL>::iterator npy = Ordered_pointy.end(); npy--;
		if(np>1 && np<5){

		    matrix direction1(1,2);
		    direction1(1,1) = (*IDX2x1)-(*npx);
		    direction1(1,2) = (*IDX2y1)-(*npy);
		    direction1 = direction1/sqrt( ((*IDX2x1)-(*npx))*((*IDX2x1)-(*npx))
						 +((*IDX2y1)-(*npy))*((*IDX2y1)-(*npy)) );
		    matrix direction2(1,2);
		    direction2(1,1) = (*IDX2x2)-(*npx);
		    direction2(1,2) = (*IDX2y2)-(*npy);
		    direction2 = direction2/sqrt( ((*IDX2x2)-(*npx))*((*IDX2x2)-(*npx))
						 +((*IDX2y2)-(*npy))*((*IDX2y2)-(*npy)) );
		    matrix direction3(1,2);
		    direction3(1,1) = (*IDX2x3)-(*npx);
		    direction3(1,2) = (*IDX2y3)-(*npy);
		    direction3 = direction3/sqrt( ((*IDX2x3)-(*npx))*((*IDX2x3)-(*npx))
						 +((*IDX2y3)-(*npy))*((*IDX2y3)-(*npy)) );

		    if( (direction1*direction.getTrans())(1,1)>0 && (direction2*direction.getTrans())(1,1)>0){
			Ordered_pointx.push_back(*IDX2x1);
			Ordered_pointy.push_back(*IDX2y1);
			points_x.erase(IDX2x1);
			points_y.erase(IDX2y1);
			direction = direction1;
		    }
		    else if( (direction1*direction.getTrans())(1,1)<0 && (direction2*direction.getTrans())(1,1)>0){
			Ordered_pointx.push_back(*IDX2x2);
			Ordered_pointy.push_back(*IDX2y2);
			points_x.erase(IDX2x2);
			points_y.erase(IDX2y2);
			direction = direction2;
		    }
		    else if( (direction2*direction.getTrans())(1,1)<0 && (direction1*direction.getTrans())(1,1)>0){
			Ordered_pointx.push_back(*IDX2x1);
			Ordered_pointy.push_back(*IDX2y1);
			points_x.erase(IDX2x1);
			points_y.erase(IDX2y1);
			direction = direction1;
		    }
		    else{
			Ordered_pointx.push_back(*IDX2x3);
			Ordered_pointy.push_back(*IDX2y3);
			points_x.erase(IDX2x3);
			points_y.erase(IDX2y3);
			direction = direction3;
		    }
		} // end if
		else{
		    Ordered_pointx.push_back(*IDX2x1);
		    Ordered_pointy.push_back(*IDX2y1);

		    direction(1,1) = (*IDX2x1)-(*npx);
		    direction(1,2) = (*IDX2y1)-(*npy);
		    direction = direction/sqrt( ((*IDX2x1)-(*npx))*((*IDX2x1)-(*npx))
					       +((*IDX2y1)-(*npy))*((*IDX2y1)-(*npy)) );
		    points_x.erase(IDX2x1);
		    points_y.erase(IDX2y1);
		} 

	    } // end while
    } // if "plain"

    if(plain_axi=="axi1"){

	    std::vector<REAL>::iterator itx = Ordered_pointx.begin();
	    std::vector<REAL>::iterator ity = Ordered_pointy.begin();
	    Ordered_pointx.push_back(*itx);
	    Ordered_pointy.push_back(*ity);

	    Split_ordered.clear();
	
	    ity = Ordered_pointy.begin();
	    for(itx=Ordered_pointx.begin(); itx<Ordered_pointx.end()-1; ++itx, ++ity){
		matrix pt_elem1(1,2);
	  	matrix pt_elem2(1,2);
		pt_elem1(1,1) = (*itx);	pt_elem1(1,2) = (*ity);
		pt_elem2(1,1) = (*(itx+1)); pt_elem2(1,2) = (*(ity+1));

		matrix pt_middle = (pt_elem1+pt_elem2)*0.5;

		if(pt_middle(1,1) > -0.5*Lx+1e-9){
		    REAL Dmin;
		    std::vector<ele2dcoupled*>::iterator it_min;
		    for(std::vector<ele2dcoupled*>::iterator it2d=element.begin(); it2d<element.end(); ++it2d){
			matrix center_node = (*it2d)->getCenter_node();
			REAL D = sqrt( (center_node(1,1)-pt_middle(1,1))*(center_node(1,1)-pt_middle(1,1))
				  +(center_node(1,2)-pt_middle(1,2))*(center_node(1,2)-pt_middle(1,2)) );

			if(it2d==element.begin()){
			    Dmin = D;
			    it_min = it2d;
			}
			else{
			    if(Dmin>D){
				Dmin = D;
				it_min = it2d;
			    }
			}
		    } // end for
		    matrix pt_elem(2,2);
		    pt_elem(1,1) = pt_elem1(1,1); pt_elem(1,2) = pt_elem1(1,2);
		    pt_elem(2,1) = pt_elem2(1,1); pt_elem(2,2) = pt_elem2(1,2);
		    (*it_min)->setOrdered_point(pt_elem);
		    Split_ordered.push_back(*it_min);
		} // end if

	    } // end for

	    Split_ordered.unique();

    } // if "axi1"

    if(plain_axi=="axi2"){
	    Split_ordered.clear();
	
	    std::vector<REAL>::iterator itx = Ordered_pointx.begin();
	    std::vector<REAL>::iterator ity = Ordered_pointy.begin();
	    for(itx=Ordered_pointx.begin(); itx<Ordered_pointx.end()-1; ++itx, ++ity){
		matrix pt_elem1(1,2);
	  	matrix pt_elem2(1,2);
		pt_elem1(1,1) = (*itx);	pt_elem1(1,2) = (*ity);
		pt_elem2(1,1) = (*(itx+1)); pt_elem2(1,2) = (*(ity+1));

		matrix pt_middle = (pt_elem1+pt_elem2)*0.5;

		if(pt_middle(1,1) > -0.5*Lx+1e-9){
		    REAL Dmin;
		    std::vector<ele2dcoupled*>::iterator it_min;
		    for(std::vector<ele2dcoupled*>::iterator it2d=element.begin(); it2d<element.end(); ++it2d){
			matrix center_node = (*it2d)->getCenter_node();
			REAL D = sqrt( (center_node(1,1)-pt_middle(1,1))*(center_node(1,1)-pt_middle(1,1))
				  +(center_node(1,2)-pt_middle(1,2))*(center_node(1,2)-pt_middle(1,2)) );

			if(it2d==element.begin()){
			    Dmin = D;
			    it_min = it2d;
			}
			else{
			    if(Dmin>D){
				Dmin = D;
				it_min = it2d;
			    }
			}
		    } // end for
		    matrix pt_elem(2,2);
		    pt_elem(1,1) = pt_elem1(1,1); pt_elem(1,2) = pt_elem1(1,2);
		    pt_elem(2,1) = pt_elem2(1,1); pt_elem(2,2) = pt_elem2(1,2);
		    (*it_min)->setOrdered_point(pt_elem);
		    Split_ordered.push_back(*it_min);
		} // end if

	    } // end for

	    Split_ordered.unique();

    } // if "axi2"

    if(plain_axi=="plain"){
	    Split_ordered.clear();
	
	    std::vector<REAL>::iterator itx = Ordered_pointx.begin();
	    std::vector<REAL>::iterator ity = Ordered_pointy.begin();
	    for(itx=Ordered_pointx.begin(); itx<Ordered_pointx.end()-1; ++itx, ++ity){
		matrix pt_elem1(1,2);
	  	matrix pt_elem2(1,2);
		pt_elem1(1,1) = (*itx);	pt_elem1(1,2) = (*ity);
		pt_elem2(1,1) = (*(itx+1)); pt_elem2(1,2) = (*(ity+1));

		matrix pt_middle = (pt_elem1+pt_elem2)*0.5;
	
		REAL Dmin;
		std::vector<ele2dcoupled*>::iterator it_min;
		for(std::vector<ele2dcoupled*>::iterator it2d=element.begin(); it2d<element.end(); ++it2d){
		    matrix center_node = (*it2d)->getCenter_node();
		    REAL D = sqrt( (center_node(1,1)-pt_middle(1,1))*(center_node(1,1)-pt_middle(1,1))
				  +(center_node(1,2)-pt_middle(1,2))*(center_node(1,2)-pt_middle(1,2)) );

		    if(it2d==element.begin()){
			Dmin = D;
			it_min = it2d;
		    }
		    else{
			if(Dmin>D){
			    Dmin = D;
			    it_min = it2d;
			}
		    }
		} // end for
		matrix pt_elem(2,2);
		pt_elem(1,1) = pt_elem1(1,1); pt_elem(1,2) = pt_elem1(1,2);
		pt_elem(2,1) = pt_elem2(1,1); pt_elem(2,2) = pt_elem2(1,2);
		(*it_min)->setOrdered_point(pt_elem);
		Split_ordered.push_back(*it_min);

	    } // end for

	    Split_ordered.unique();

    } // if "plain"

    // set int Split_ordered in ele2dcoupled
    int Split_num = 0;
    for(std::list<ele2dcoupled*>::iterator it2d=Split_ordered.begin(); it2d!=Split_ordered.end(); ++it2d){
	Split_num++;
	(*it2d)->setSplit_ordered(Split_num);
    }

} // sort_points()


void surfacetension::find_proj_normal_deriv_pt_parallel(matrix & p_proj, matrix & t, matrix & dtdx, matrix & p){

    matrix Coord_Active_points, Foot_points, n_Foot_points, a_coeff, a_coeff_n;

    for(std::list<activepoint*>::const_iterator it=Active_points.begin(); it!=Active_points.end(); ++it){

	Coord_Active_points.appendRow( ((*it)->getCoords()).getRow(1) );
	Foot_points.appendRow( ((*it)->getFoot_points()).getRow(1) );
	n_Foot_points.appendRow( ((*it)->getN_Foot_points()).getRow(1) );
	a_coeff.appendRow( ((*it)->getA_coeff()).getRow(1) );
 	a_coeff_n.appendRow( ((*it)->getA_coeff_n()).getRow(1) );
    }

    matrix a_coeff_sym, a_coeff_n_sym;
    if(plain_axi=="axi1"){
	    matrix temp_Foot_points = Foot_points;
	    Foot_points = zeros(2*temp_Foot_points.num_row,2);
#pragma omp parallel for
	    for(int i=1; i<temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = temp_Foot_points(i,1);
		Foot_points(i,2) = temp_Foot_points(i,2);
	    }
#pragma omp parallel for
	    for(int i=1+temp_Foot_points.num_row; i<2*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row,1))-0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-temp_Foot_points.num_row,2);
	    }

	    // n_Foot_points
	    matrix n_Foot_points_col1;
	    n_Foot_points_col1.appendCol(n_Foot_points.getCol(1));
	    n_Foot_points_col1 = -1*n_Foot_points_col1;
	    matrix sym_points2;
	    sym_points2.appendCol(n_Foot_points_col1.getCol(1));
 	    sym_points2.appendCol(n_Foot_points.getCol(2));
	    for(int i=1; i<sym_points2.num_row+1; ++i){
		n_Foot_points.appendRow(sym_points2.getRow(i));
	    }

	    // a_coeff
	    a_coeff_sym = a_coeff;
	    for(int i=1; i<a_coeff.num_row+1; ++i){
		a_coeff_sym.appendRow(a_coeff.getRow(i));
	    }

	    // a_coeff_n
	    a_coeff_n_sym = a_coeff_n;
	    for(int i=1; i<a_coeff_n.num_row+1; ++i){
		a_coeff_n_sym.appendRow(a_coeff_n.getRow(i));
	    }

	    // Coord_Active_points
	    matrix temp_coord = Coord_Active_points;
	    Coord_Active_points = zeros(2*temp_coord.num_row,2);
	    for(int i=1; i<temp_coord.num_row+1; ++i){
		Coord_Active_points(i,1) = temp_coord(i,1);
		Coord_Active_points(i,2) = temp_coord(i,2);
	    }
	    for(int i=1+temp_coord.num_row; i<2*temp_coord.num_row+1; ++i){
		Coord_Active_points(i,1) = -(Lx*0.5+temp_coord(i-temp_coord.num_row,1))-0.5*Lx;
		Coord_Active_points(i,2) = temp_coord(i-temp_coord.num_row,2);
	    }
	    
    } // if "axi1"

    if(plain_axi=="axi2"){
	    matrix temp_Foot_points = Foot_points;
	    Foot_points = zeros(3*temp_Foot_points.num_row,2);
	    for(int i=1; i<temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = temp_Foot_points(i,1);
		Foot_points(i,2) = temp_Foot_points(i,2);
	    }
	    for(int i=1+temp_Foot_points.num_row; i<2*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row,1))-0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-temp_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_Foot_points.num_row; i<3*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(-Lx*0.5+temp_Foot_points(i-2*temp_Foot_points.num_row,1))+0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-2*temp_Foot_points.num_row,2);
	    }

	    // Coord_Active_points
	    matrix temp_coord = Coord_Active_points;
	    Coord_Active_points = zeros(3*temp_coord.num_row,2);
	    for(int i=1; i<temp_coord.num_row+1; ++i){
		Coord_Active_points(i,1) = temp_coord(i,1);
		Coord_Active_points(i,2) = temp_coord(i,2);
	    }
	    for(int i=1+temp_coord.num_row; i<2*temp_coord.num_row+1; ++i){
		Coord_Active_points(i,1) = -(Lx*0.5+temp_coord(i-temp_coord.num_row,1))-0.5*Lx;
		Coord_Active_points(i,2) = temp_coord(i-temp_coord.num_row,2);
	    }
	    for(int i=1+2*temp_coord.num_row; i<3*temp_coord.num_row+1; ++i){
		Coord_Active_points(i,1) = -(-Lx*0.5+temp_coord(i-2*temp_coord.num_row,1))+0.5*Lx;
		Coord_Active_points(i,2) = temp_coord(i-2*temp_coord.num_row,2);
	    }

	    // n_Foot_points
	    matrix temp_n_Foot_points = n_Foot_points;
	    n_Foot_points = zeros(3*temp_n_Foot_points.num_row,2);
	    for(int i=1; i<temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = temp_n_Foot_points(i,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i,2);
	    }
	    for(int i=1+temp_n_Foot_points.num_row; i<2*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = -temp_n_Foot_points(i-temp_n_Foot_points.num_row,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i-temp_n_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_n_Foot_points.num_row; i<3*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = -temp_n_Foot_points(i-2*temp_n_Foot_points.num_row,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i-2*temp_n_Foot_points.num_row,2);
	    }

	    // a_coeff
	    a_coeff_sym = zeros(3*a_coeff.num_row,2);
	    for(int i=1; i<a_coeff.num_row+1; ++i){
		a_coeff_sym(i,1) = a_coeff(i,1);
		a_coeff_sym(i,2) = a_coeff(i,2);
	    }
	    for(int i=1+a_coeff.num_row; i<2*a_coeff.num_row+1; ++i){
		a_coeff_sym(i,1) = a_coeff(i-a_coeff.num_row,1);
		a_coeff_sym(i,2) = a_coeff(i-a_coeff.num_row,2);
	    }
	    for(int i=1+2*a_coeff.num_row; i<3*a_coeff.num_row+1; ++i){
		a_coeff_sym(i,1) = a_coeff(i-2*a_coeff.num_row,1);
		a_coeff_sym(i,2) = a_coeff(i-2*a_coeff.num_row,2);
	    }

	    // a_coeff_n
	    a_coeff_n_sym = zeros(3*a_coeff_n.num_row,2);
	    for(int i=1; i<a_coeff_n.num_row+1; ++i){
		a_coeff_n_sym(i,1) = a_coeff_n(i,1);
		a_coeff_n_sym(i,2) = a_coeff_n(i,2);
	    }
	    for(int i=1+a_coeff_n.num_row; i<2*a_coeff_n.num_row+1; ++i){
		a_coeff_n_sym(i,1) = a_coeff_n(i-a_coeff_n.num_row,1);
		a_coeff_n_sym(i,2) = a_coeff_n(i-a_coeff_n.num_row,2);
	    }
	    for(int i=1+2*a_coeff_n.num_row; i<3*a_coeff_n.num_row+1; ++i){
		a_coeff_n_sym(i,1) = a_coeff_n(i-2*a_coeff_n.num_row,1);
		a_coeff_n_sym(i,2) = a_coeff_n(i-2*a_coeff_n.num_row,2);
	    }

    } // if "axi2"

    if(plain_axi=="plain"){
	    a_coeff_sym = a_coeff;
	    a_coeff_n_sym = a_coeff_n;

    } // if "plain"

//    matrix D = zeros(Coord_Active_points.num_row, 1);
    std::list<REAL> D;
    for(int i=1; i<Coord_Active_points.num_row+1; ++i){
	REAL D_val = pow( (Coord_Active_points(i,1)-p(1,1))*(Coord_Active_points(i,1)-p(1,1))
		        + (Coord_Active_points(i,2)-p(1,2))*(Coord_Active_points(i,2)-p(1,2)), 0.5);
	D.push_back(D_val);
    }

    std::vector<int> IDX = sort30(D.begin(), D.end());

//    matrix D2 = zeros(30,1);
    std::list<REAL> D2;
    for(int i=0; i<30; ++i){
    	REAL D2_val = pow( (Foot_points(IDX[i],1)-p(1,1))*(Foot_points(IDX[i],1)-p(1,1))
		     	 + (Foot_points(IDX[i],2)-p(1,2))*(Foot_points(IDX[i],2)-p(1,2)), 0.5);
	D2.push_back(D2_val);
    }

    std::vector<int> IDX2 = sort30(D2.begin(), D2.end());	// numbering start from 1

    matrix y_sorted(30,2);
	
    for(int i=0; i<30; ++i){
	y_sorted(i+1, 1) = Foot_points(IDX[IDX2[i]-1], 1);
	y_sorted(i+1, 2) = Foot_points(IDX[IDX2[i]-1], 2);
    }

    matrix n0_old(1,2);
    n0_old(1,1) = n_Foot_points(IDX[IDX2[0]-1], 1);
    n0_old(1,2) = n_Foot_points(IDX[IDX2[0]-1], 2);

    matrix y0(1,2);
    y0(1,1) = y_sorted(1,1);
    y0(1,2) = y_sorted(1,2);

    int in2 = IDX[IDX2[0]-1];

    //Compute the new foot point (closest to P, so equivalent to rebuilding the distance function)
    //first, find P in local coordinates

    matrix temp_n0(2,2);
    temp_n0(1,1) = n0_old(1,2); temp_n0(1,2) = -n0_old(1,1);
    temp_n0(2,1) = n0_old(1,1); temp_n0(2,2) = n0_old(1,2);
    matrix p_local = temp_n0*(p-y0).getTrans();	// p_local is 2x1

    std::complex<double> a = 2.0*a_coeff_sym(in2,3)*a_coeff_sym(in2,3);
    std::complex<double> b = 3.0*a_coeff_sym(in2,2)*a_coeff_sym(in2,3);
    std::complex<double> c = 1.0+a_coeff_sym(in2,2)*a_coeff_sym(in2,2)
			     + 2.0*a_coeff_sym(in2,3)*a_coeff_sym(in2,1)-2*a_coeff_sym(in2,3)*p_local(2,1);
    std::complex<double> d = a_coeff_sym(in2,2)*a_coeff_sym(in2,1)-a_coeff_sym(in2,2)*p_local(2,1)-p_local(1,1);

    std::complex<double> i(0, 1);

    std::complex<double> x1 = -b/(3.0*a) - (pow(2.0,1.0/3.0)*(-b*b + 3.0*a*c)) /(3.0*a*pow( -2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt( 4.0*pow(-b*b + 3.0*a*c, 3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d, 2.0) ), 1.0/3.0) ) 
+ pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt(4.0*pow(-b*b + 3.0*a*c,3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d,2.0) ), 1.0/3.0) /  (3.0*pow(2.0,1.0/3.0)*a);

    std::complex<double> x2 = -b/(3.0*a) + ((1.0 + i*3.0)*(-b*b + 3.0*a*c)) / ( 3.0*pow(2.0,2.0/3.0)*a*pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt(4.0*pow(-b*b + 3.0*a*c,3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d,2.0) ), 1.0/3.0) ) 
- (1.0 - i*sqrt(3.0))*pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt(4.0*pow(-b*b + 3.0*a*c, 3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d,2.0) ), 1.0/3.0) / (6.0*pow(2.0,1.0/3.0)*a);

    std::complex<double> x3 = -b/(3.0*a) + ( (1.0 - i*sqrt(3.0))*(-b*b + 3.0*a*c)) / (3.0*pow(2.0,2.0/3.0)*a*pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt(4.0*pow(-b*b + 3.0*a*c,3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d,2.0) ), 1.0/3.0) )
 - (1.0 + i*sqrt(3.0))*pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d + sqrt(4.0*pow(-b*b + 3.0*a*c,3.0) + pow(-2.0*b*b*b + 9.0*a*b*c - 27.0*a*a*d,2.0) ), 1.0/3.0)  / (6.0*pow(2.0,1.0/3.0)*a);

    //eliminate imaginary part if small enough

    if( imag(x1)<1e-8 ){
	x1 = real(x1);
    }

    if( imag(x1)<1e-8 ){
	x2 = real(x2);
    }

    if( imag(x1)<1e-8){
	x3 = real(x3);
    }

    //find the roots (minimum distance)

    matrix roots;	// row vector
    if(imag(x1)==0){	// real
 	matrix temp(1,1);
	temp(1,1) = real(x1);
	roots.appendCol(temp.getCol(1));
    }

    if(imag(x2)==0){	// real
	matrix temp(1,1);
	temp(1,1) = real(x2);
	roots.appendCol(temp.getCol(1));
    }

    if(imag(x3)==0){	// real
	matrix temp(1,1);
	temp(1,1) = real(x3);
	roots.appendCol(temp.getCol(1));
    }

    matrix distance;	// row vector
    for(int k=1; k<length(roots)+1; ++k){
	matrix temp(1,1);
	temp(1,1) = sqrt( (roots(1,k)-p_local(1,1))*(roots(1,k)-p_local(1,1))
			+ pow( a_coeff_sym(in2,1)+a_coeff_sym(in2,2)*roots(1,k)
			+ a_coeff_sym(in2,3)*roots(1,k)*roots(1,k)-p_local(2,1), 2) );
	distance.appendCol(temp.getCol(1));
    }

    REAL distanceMin = distance(1,1);
    int j=1;
    for(int k=2; k<distance.num_col+1; ++k){
	if(distanceMin > distance(1,k)){
	    distanceMin = distance(1,k);
	    j = k;
	}
    }

    //find foot point in global coordinates
    temp_n0 = zeros(2,2);	
    temp_n0(1,1) = n0_old(1,2); temp_n0(1,2) = -n0_old(1,1);
    temp_n0(2,1) = n0_old(1,1); temp_n0(2,2) = n0_old(1,2);

    matrix temp_roots(2,1); 
    temp_roots(1,1) = roots(1,j);
    matrix temp_vec(3,1);
    temp_vec(1,1) = 1; temp_vec(2,1) = roots(1,j); temp_vec(3,1) = roots(1,j)*roots(1,j);
    matrix a_coeff_sym_in2(1,3);
    a_coeff_sym_in2(1,1) = a_coeff_sym(in2,1);
    a_coeff_sym_in2(1,2) = a_coeff_sym(in2,2);
    a_coeff_sym_in2(1,3) = a_coeff_sym(in2,3);
    matrix temp_mat = a_coeff_sym_in2*temp_vec;
    temp_roots(2,1) = temp_mat(1,1);
    p_proj = (temp_n0%temp_roots).getTrans() + y0;

    t = zeros(1,2);
    t(1,1) = (a_coeff_n_sym(in2,2)+a_coeff_n_sym(in2,4)*roots(1,j)+a_coeff_n_sym(in2,6)*roots(1,j)*roots(1,j));
    t(1,2) = -(a_coeff_n_sym(in2,1)+a_coeff_n_sym(in2,3)*roots(1,j)+a_coeff_n_sym(in2,5)*roots(1,j)*roots(1,j));

    matrix temp_coeff(1,2);
    matrix temp_root(2,1);
    temp_coeff(1,1) = a_coeff_n_sym(in2,4); temp_coeff(1,2) = 2*a_coeff_n_sym(in2,6);
    temp_root(1,1) = 1; temp_root(2,1) = roots(1,j);

    REAL dt1dx1 = ( temp_coeff*temp_root*t(1,1) )(1,1);
    REAL dt1dx2 = ( temp_coeff*temp_root*t(1,2) )(1,1);

    REAL dt2dx1 = -( temp_coeff*temp_root*t(1,1) )(1,1);

    temp_coeff = zeros(1,2);
    temp_coeff(1,1) = a_coeff_n_sym(in2,3); temp_coeff(1,2) = 2*a_coeff_n_sym(in2,5);
    REAL dt2dx2 = -( temp_coeff*temp_root*t(1,2) )(1,1);

    dtdx = zeros(2,2);
    dtdx(1,1) = dt1dx1; dtdx(1,2) = dt1dx2;
    dtdx(2,1) = dt2dx1; dtdx(2,2) = dt2dx2;

} // find_proj_normal_deriv_pt_parallel()


void surfacetension::find_velocity_gradient_d(matrix &d,matrix &d_cont,matrix &E,matrix &Fd,REAL &Ja,matrix &p){

    matrix Foot_points, n_Foot_points, a_coeff3, a_coeff_vel, a_coeff_n, a_coeff_E, a_coeff_F, a_coeff_J;

    for(std::list<activepoint*>::const_iterator it=Active_points.begin(); it!=Active_points.end(); ++it){

	Foot_points.appendRow( ((*it)->getFoot_points()).getRow(1) );
	n_Foot_points.appendRow( ((*it)->getN_Foot_points()).getRow(1) );
	a_coeff3.appendRow( ((*it)->getA_coeff3()).getRow(1) );
	a_coeff_vel.appendRow( ((*it)->getA_coeff_vel()).getRow(1) );
 	a_coeff_n.appendRow( ((*it)->getA_coeff_n()).getRow(1) );
	a_coeff_E.appendRow( ((*it)->getA_coeff_E()).getRow(1) );
	a_coeff_F.appendRow( ((*it)->getA_coeff_F()).getRow(1) );
	a_coeff_J.appendRow( ((*it)->getA_coeff_J()).getRow(1) );
    }

    matrix a_coeff3_sym, a_coeff_vel_sym, a_coeff_n_sym, a_coeff_E_sym, a_coeff_F_sym, a_coeff_J_sym;
    if(plain_axi=="axi1"){

	    matrix temp_Foot_points = Foot_points;
	    Foot_points = zeros(2*temp_Foot_points.num_row,2);
	    for(int i=1; i<temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = temp_Foot_points(i,1);
		Foot_points(i,2) = temp_Foot_points(i,2);
	    }
	    for(int i=1+temp_Foot_points.num_row; i<2*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row,1))-0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-temp_Foot_points.num_row,2);
	    }

	    // n_Foot_points
	    matrix n_Foot_points_col1;
	    n_Foot_points_col1.appendCol(n_Foot_points.getCol(1));
	    n_Foot_points_col1 = -1*n_Foot_points_col1;
	    matrix sym_points2;
	    sym_points2.appendCol(n_Foot_points_col1.getCol(1));
 	    sym_points2.appendCol(n_Foot_points.getCol(2));
	    for(int i=1; i<sym_points2.num_row+1; ++i){
		n_Foot_points.appendRow(sym_points2.getRow(i));
	    }

	    // a_coeff3
	    a_coeff3_sym = a_coeff3;
	    for(int i=1; i<a_coeff3.num_row+1; ++i){
		a_coeff3_sym.appendRow(a_coeff3.getRow(i));
	    }

	    // a_coeff_vel
	    a_coeff_vel_sym = a_coeff_vel;
	    for(int i=1; i<a_coeff_vel.num_row+1; ++i){
		a_coeff_vel_sym.appendRow(a_coeff_vel.getRow(i));
	    }

	    // a_coeff_n
	    a_coeff_n_sym = a_coeff_n;
	    for(int i=1; i<a_coeff_n.num_row+1; ++i){
		a_coeff_n_sym.appendRow(a_coeff_n.getRow(i));
	    }

	    // a_coeff_E
	    a_coeff_E_sym = a_coeff_E;
	    for(int i=1; i<a_coeff_E.num_row+1; ++i){
		a_coeff_E_sym.appendRow(a_coeff_E.getRow(i));
	    }

	    // a_coeff_F
	    a_coeff_F_sym = a_coeff_F;
	    for(int i=1; i<a_coeff_F.num_row+1; ++i){
		a_coeff_F_sym.appendRow(a_coeff_F.getRow(i));
	    }

	    // a_coeff_J
	    a_coeff_J_sym = a_coeff_J;
	    for(int i=1; i<a_coeff_J.num_row+1; ++i){
		a_coeff_J_sym.appendRow(a_coeff_J.getRow(i));
	    }

	    
    } // if "axi1"

    if(plain_axi=="axi2"){
	    matrix temp_Foot_points = Foot_points;
	    Foot_points = zeros(3*temp_Foot_points.num_row,2);
	    for(int i=1; i<temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = temp_Foot_points(i,1);
		Foot_points(i,2) = temp_Foot_points(i,2);
	    }
	    for(int i=1+temp_Foot_points.num_row; i<2*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row,1))-0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-temp_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_Foot_points.num_row; i<3*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(-Lx*0.5+temp_Foot_points(i-2*temp_Foot_points.num_row,1))+0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-2*temp_Foot_points.num_row,2);
	    }

	    // n_Foot_points
	    matrix temp_n_Foot_points = n_Foot_points;
	    n_Foot_points = zeros(3*temp_n_Foot_points.num_row,2);
	    for(int i=1; i<temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = temp_n_Foot_points(i,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i,2);
	    }
	    for(int i=1+temp_n_Foot_points.num_row; i<2*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = -temp_n_Foot_points(i-temp_n_Foot_points.num_row,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i-temp_n_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_n_Foot_points.num_row; i<3*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = -temp_n_Foot_points(i-2*temp_n_Foot_points.num_row,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i-2*temp_n_Foot_points.num_row,2);
	    }

	    // a_coeff3
	    a_coeff3_sym = a_coeff3;
	    for(int i=1; i<a_coeff3.num_row+1; ++i){
		a_coeff3_sym.appendRow(a_coeff3.getRow(i));
	    }
	    for(int i=1; i<a_coeff3.num_row+1; ++i){
		a_coeff3_sym.appendRow(a_coeff3.getRow(i));
	    }

	    // a_coeff_vel
	    a_coeff_vel_sym = a_coeff_vel;
	    for(int i=1; i<a_coeff_vel.num_row+1; ++i){
		a_coeff_vel_sym.appendRow(a_coeff_vel.getRow(i));
	    }
	    for(int i=1; i<a_coeff3.num_row+1; ++i){
		a_coeff3_sym.appendRow(a_coeff3.getRow(i));
	    }


	    // a_coeff_n
	    a_coeff_n_sym = a_coeff_n;
	    for(int i=1; i<a_coeff_n.num_row+1; ++i){
		a_coeff_n_sym.appendRow(a_coeff_n.getRow(i));
	    }
	    for(int i=1; i<a_coeff_n.num_row+1; ++i){
		a_coeff_n_sym.appendRow(a_coeff_n.getRow(i));
	    }

	    // a_coeff_E
	    a_coeff_E_sym = a_coeff_E;
	    for(int i=1; i<a_coeff_E.num_row+1; ++i){
		a_coeff_E_sym.appendRow(a_coeff_E.getRow(i));
	    }
	    for(int i=1; i<a_coeff_E.num_row+1; ++i){
		a_coeff_E_sym.appendRow(a_coeff_E.getRow(i));
	    }

	    // a_coeff_F
	    a_coeff_F_sym = a_coeff_F;
	    for(int i=1; i<a_coeff_F.num_row+1; ++i){
		a_coeff_F_sym.appendRow(a_coeff_F.getRow(i));
	    }
	    for(int i=1; i<a_coeff_F.num_row+1; ++i){
		a_coeff_F_sym.appendRow(a_coeff_F.getRow(i));
	    }

	    // a_coeff_J
	    a_coeff_J_sym = a_coeff_J;
	    for(int i=1; i<a_coeff_J.num_row+1; ++i){
		a_coeff_J_sym.appendRow(a_coeff_J.getRow(i));
	    }
	    for(int i=1; i<a_coeff_J.num_row+1; ++i){
		a_coeff_J_sym.appendRow(a_coeff_J.getRow(i));
	    }

    } // if "axi2"

    if(plain_axi=="plain"){
            a_coeff3_sym = a_coeff3;
            a_coeff_vel_sym = a_coeff_vel;
            a_coeff_n_sym = a_coeff_n;
            a_coeff_E_sym = a_coeff_E;
            a_coeff_F_sym = a_coeff_F;
            a_coeff_J_sym = a_coeff_J;

    } // if "plain"

//    matrix D2 = zeros(Foot_points.num_row, 1);
    REAL D2_min;
    int IDX2;
    for(int i=1; i<Foot_points.num_row+1; ++i){
	REAL D2 = pow( (Foot_points(i,1)-p(1,1))*(Foot_points(i,1)-p(1,1))
		     + (Foot_points(i,2)-p(1,2))*(Foot_points(i,2)-p(1,2)), 0.5);
 	if(i==1){
	    D2_min = D2;
	    IDX2 = i;
    	}
	else if(D2_min>D2){
	    D2_min = D2;
	    IDX2 = i;
	}
    }

    matrix n0_old(1,2);
    n0_old(1,1) = n_Foot_points(IDX2, 1);
    n0_old(1,2) = n_Foot_points(IDX2, 2);

    matrix y0(1,2);
    y0(1,1) = Foot_points(IDX2, 1);
    y0(1,2) = Foot_points(IDX2, 2);

    int in2 = IDX2;

    //Compute the new foot point (closed to P, so equivalent to rebuilding the distance function)
    //first, find P in local coordinates
    matrix p_local;
    matrix temp_n0(2,2);
    temp_n0(1,1) = n0_old(1,2); temp_n0(1,2) = -n0_old(1,1);
    temp_n0(2,1) = n0_old(1,1); temp_n0(2,2) = n0_old(1,2);

    p_local = temp_n0*(p-y0).getTrans();	// p_local is 2x1

    E(1,1) = a_coeff_E_sym(in2,1)+a_coeff_E_sym(in2,3)*p_local(1,1)+a_coeff_E_sym(in2,5)*p_local(1,1)*p_local(1,1);
    E(1,2) = a_coeff_E_sym(in2,2)+a_coeff_E_sym(in2,4)*p_local(1,1)+a_coeff_E_sym(in2,6)*p_local(1,1)*p_local(1,1);

    Fd(1,1) = a_coeff_F_sym(in2,1)+a_coeff_F_sym(in2,3)*p_local(1,1)+a_coeff_F_sym(in2,5)*p_local(1,1)*p_local(1,1);
    Fd(1,2) = a_coeff_F_sym(in2,2)+a_coeff_F_sym(in2,4)*p_local(1,1)+a_coeff_F_sym(in2,6)*p_local(1,1)*p_local(1,1);

    Ja = a_coeff_J_sym(in2,1)+a_coeff_J_sym(in2,2)*p_local(1,1)+a_coeff_J_sym(in2,3)*p_local(1,1)*p_local(1,1);

    matrix v(1,2);
    v(1,1) = a_coeff_vel_sym(in2,1)+a_coeff_vel_sym(in2,3)*p_local(1,1)+a_coeff_vel_sym(in2,5)*p_local(1,1)*p_local(1,1);
    v(1,2) = a_coeff_vel_sym(in2,2)+a_coeff_vel_sym(in2,4)*p_local(1,1)+a_coeff_vel_sym(in2,6)*p_local(1,1)*p_local(1,1);

    matrix n(1,2);
    n(1,1) = a_coeff_n_sym(in2,1)+a_coeff_n_sym(in2,3)*p_local(1,1)+a_coeff_n_sym(in2,5)*p_local(1,1)*p_local(1,1);
    n(1,2) = a_coeff_n_sym(in2,2)+a_coeff_n_sym(in2,4)*p_local(1,1)+a_coeff_n_sym(in2,6)*p_local(1,1)*p_local(1,1);

    REAL v_curvi = v(1,1)*n(1,2)- v(1,2)*n(1,1);

    REAL dvds1 = (a_coeff_vel_sym(in2,3)+2*a_coeff_vel_sym(in2,5)*p_local(1,1))*n(1,2) 
	       + v(1,1)*(a_coeff_n_sym(in2,4)+2*a_coeff_n_sym(in2,6)*p_local(1,1))
    	       - (a_coeff_vel_sym(in2,4)+2*a_coeff_vel_sym(in2,6)*p_local(1,1))*n(1,1) 
	       - v(1,2)*(a_coeff_n_sym(in2,3)+2*a_coeff_n_sym(in2,5)*p_local(1,1));

    REAL vn = v(1,1)*n(1,1)+ v(1,2)*n(1,2);


    if(plain_axi=="axi1"){
    	
            if(p(1,1)+0.5*Lx>1e-4){
            	REAL a = pow( 1+pow(a_coeff3_sym(in2,2)+2*a_coeff3_sym(in2,3)*p_local(1,1) 
		                    + 3*a_coeff3_sym(in2,4)*p_local(1,1)*p_local(1,1), 2),  0.5);
            
            	REAL a_prime = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))
			         * ( a_coeff3_sym(in2,2)+2*a_coeff3_sym(in2,3)*p_local(1,1)
				   + 3*a_coeff3_sym(in2,4)*p_local(1,1)*p_local(1,1) )/a;
            
	    	matrix curvature_cov(1,2);          
            	curvature_cov(1,1) = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))/a;
            	curvature_cov(1,2) =-(1/a)*n(1,1)*(p(1,1)+Lx*0.5);
            
            	d(1,1) = 2*(a*a*(dvds1 + v_curvi*a_prime/a)) - 2*vn*curvature_cov(1,1);
            	d(1,2) = 2*(((p(1,1)+Lx*0.5)*n(1,2))*v_curvi) - 2*vn*curvature_cov(1,2);
            
            	d_cont(1,1) = 2*((dvds1 + v_curvi*a_prime/a)) 
			    - 2*vn*curvature_cov(1,1)/(a*a);
            	d_cont(1,2) = 2*(n(1,2)*v_curvi/(p(1,1)+Lx*0.5)) 
			    - 2*vn*curvature_cov(1,2)/((p(1,1)+Lx*0.5)*(p(1,1)+Lx*0.5));
	    }
	    else{
            	REAL a = pow( 1+pow(a_coeff3_sym(in2,2)+2*a_coeff3_sym(in2,3)*p_local(1,1) 
		                    + 3*a_coeff3_sym(in2,4)*p_local(1,1)*p_local(1,1), 2),  0.5);
            
            	REAL a_prime = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))
			         * ( a_coeff3_sym(in2,2)+2*a_coeff3_sym(in2,3)*p_local(1,1)
				   + 3*a_coeff3_sym(in2,4)*p_local(1,1)*p_local(1,1) )/a;
            
	    	matrix curvature_cov(1,2);          
            	curvature_cov(1,1) = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))/a;
            	curvature_cov(1,2) =-(1/a)*n(1,1)*(p(1,1)+Lx*0.5);
            
            	d(1,1) = 2*(a*a*(dvds1 + v_curvi*a_prime/a)) - 2*vn*curvature_cov(1,1);
            	d(1,2) = 2*(((p(1,1)+Lx*0.5)*n(1,2))*v_curvi) - 2*vn*curvature_cov(1,2);
            
            	d_cont(1,1) = 2*((dvds1 + v_curvi*a_prime/a)) 
			    - 2*vn*curvature_cov(1,1)/(a*a);
            	d_cont(1,2) = 2*((dvds1 + v_curvi*a_prime/a)) 
			    - 2*vn*curvature_cov(1,1)/(a*a);

   	    }

    } // if "axi1"
        
    if(plain_axi=="axi2"){
            if(p(1,1)+0.5*Lx>1e-4){
            	REAL a = pow( 1+pow(a_coeff3_sym(in2,2)+2*a_coeff3_sym(in2,3)*p_local(1,1) 
		                    + 3*a_coeff3_sym(in2,4)*p_local(1,1)*p_local(1,1), 2),  0.5);
            
            	REAL a_prime = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))
			         * ( a_coeff3_sym(in2,2)+2*a_coeff3_sym(in2,3)*p_local(1,1)
				   + 3*a_coeff3_sym(in2,4)*p_local(1,1)*p_local(1,1) )/a;
            
	    	matrix curvature_cov(1,2);          
            	curvature_cov(1,1) = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))/a;
            	curvature_cov(1,2) =-(1/a)*n(1,1)*(p(1,1)+Lx*0.5);
            
            	d(1,1) = 2*(a*a*(dvds1 + v_curvi*a_prime/a)) - 2*vn*curvature_cov(1,1);
            	d(1,2) = 2*(((p(1,1)+Lx*0.5)*n(1,2))*v_curvi) - 2*vn*curvature_cov(1,2);
            
            	d_cont(1,1) = 2*(dvds1 + v_curvi*a_prime/a) 
			    - 2*vn*curvature_cov(1,1)/(a*a);
            	d_cont(1,2) = 2*(n(1,2)*v_curvi/(p(1,1)+Lx*0.5)) 
			    - 2*vn*curvature_cov(1,2)/((p(1,1)+Lx*0.5)*(p(1,1)+Lx*0.5));
	    }
	    else{
            	REAL a = pow( 1+pow(a_coeff3_sym(in2,2)+2*a_coeff3_sym(in2,3)*p_local(1,1) 
		                    + 3*a_coeff3_sym(in2,4)*p_local(1,1)*p_local(1,1), 2),  0.5);
            
            	REAL a_prime = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))
			         * ( a_coeff3_sym(in2,2)+2*a_coeff3_sym(in2,3)*p_local(1,1)
				   + 3*a_coeff3_sym(in2,4)*p_local(1,1)*p_local(1,1) )/a;
            
	    	matrix curvature_cov(1,2);          
            	curvature_cov(1,1) = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))/a;
            	curvature_cov(1,2) =-(1/a)*n(1,1)*(p(1,1)+Lx*0.5);
            
            	d(1,1) = 2*(a*a*(dvds1 + v_curvi*a_prime/a)) - 2*vn*curvature_cov(1,1);
            	d(1,2) = 2*(((p(1,1)+Lx*0.5)*n(1,2))*v_curvi) - 2*vn*curvature_cov(1,2);
            
            	d_cont(1,1) = 2*(dvds1 + v_curvi*a_prime/a) 
			    - 2*vn*curvature_cov(1,1)/(a*a);
            	d_cont(1,2) = 2*((dvds1 + v_curvi*a_prime/a)) 
			    - 2*vn*curvature_cov(1,1)/(a*a);

   	    }

    } // if "axi2"      
	
    if(plain_axi=="plain"){
	    REAL a = pow( 1+pow(a_coeff3_sym(in2,2)+2*a_coeff3_sym(in2,3)*p_local(1,1) 
		                + 3*a_coeff3_sym(in2,4)*p_local(1,1)*p_local(1,1), 2),  0.5);
            
            REAL a_prime = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))
			     * ( a_coeff3_sym(in2,2)+2*a_coeff3_sym(in2,3)*p_local(1,1)
				+ 3*a_coeff3_sym(in2,4)*p_local(1,1)*p_local(1,1) )/a;
            
	    matrix curvature_cov(1,2);          
            curvature_cov(1,1) = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))/a;
            curvature_cov(1,2) = 0;
            
            d(1,1) = 2*(a*a*(dvds1 + v_curvi*a_prime/a)) - 2*vn*curvature_cov(1,1);
            d(1,2) = 0;
            
            d_cont(1,1) = 2*(dvds1 + v_curvi*a_prime/a) 
		        - 2*vn*curvature_cov(1,1)/(a*a);
            d_cont(1,2) = 0;   
       
    } // if "plain"

} // find_velocity_gradient_d()


matrix surfacetension::velocity_pt_fluid(matrix & pt, const std::string in_out){

    int BCjump;
    if(in_out=="in"){
	BCjump = 1;
    }
    if(in_out=="out"){
	BCjump = 2;
    }

    // [~, k0] = min(D);
    std::vector<ele2dcoupled*>::const_iterator itk0;
    REAL Dmin;
    for(std::vector<ele2dcoupled*>::const_iterator it2d=element.begin(); it2d<element.end(); ++it2d){
	REAL D =  sqrt(( ((*it2d)->getCenter_node())(1,1)-pt(1,1))*( ((*it2d)->getCenter_node())(1,1)-pt(1,1))
		     + ( ((*it2d)->getCenter_node())(1,2)-pt(1,2))*( ((*it2d)->getCenter_node())(1,2)-pt(1,2)));

	if(it2d==element.begin()){
	    Dmin = D;
	    itk0 = it2d;
	}
	else if(Dmin > D){
	    Dmin = D;
	    itk0 = it2d;
	}
    }

    matrix t_element = zeros(9,2);
    if((*itk0)->getSplit_ordered()>0){
	matrix t_element_total = (*itk0)->getTEleTotal();
	for(int i=1; i<10; ++i){
	    t_element(i,1) = t_element_total(1,2*i-1);
	    t_element(i,2) = t_element_total(1,2*i);
	}
    }
//    else
//	t_element = zeros(9,2);

    int num_Split_ordered = Split_ordered.size();
    (*itk0)->assembly_fluid_XFEM_ndiscont_parallel(numnode9, numsnode92, numnode4, numsnode42,
					           num_Split_ordered, close_DBC2);

    matrix sctrBv = (*itk0)->getSctrBv();

    int nn9 = sctrBv.num_col;

    matrix Ve(nn9,1);

    for(int i=1; i<nn9+1; ++i){
	Ve(i,1) = dsol(sctrBv(1,i), 1);
    }

    matrix pt_par(1,2);
    pt_par(1,1) = ( pt(1,1)-((*itk0)->getCenter_node())(1,1) )/dx9;
    pt_par(1,2) = ( pt(1,2)-((*itk0)->getCenter_node())(1,2) )/dy9;

    (*itk0)->calcN_9(pt_par(1,1), pt_par(1,2));
    (*itk0)->xfem_Nmu9_parallel(BCjump);
    matrix N9 = (*itk0)->getN9();

    matrix Nv(2,nn9);

    for(int i=1; i<10; ++i){
	Nv(1,2*i-1) = N9(i,1);
	Nv(2,2*i)   = N9(i,1);
    }
    for(int i=19; i<nn9+1; ++i){
	if(i-9>length(N9) || i-18>9)
	    break;
	Nv(1,i) = N9(i-9,1)*t_element(i-18,1);
	Nv(2,i) = N9(i-9,1)*t_element(i-18,2);
    }

    matrix v = Nv*Ve;	// v is 2x1
   
    return v.getTrans();

} // velocity_pt_fluid()


REAL surfacetension::velocity_pt_membrane(matrix & pt){

    // [~, k0] = min(D);
    std::list<ele2dcoupled*>::const_iterator itk;
    REAL Dmin;
    for(std::list<ele2dcoupled*>::const_iterator it2d=Split_ordered.begin(); 
				   		 it2d!=Split_ordered.end(); ++it2d){

	REAL D =  sqrt(( ((*it2d)->getCenter_node())(1,1)-pt(1,1))*( ((*it2d)->getCenter_node())(1,1)-pt(1,1))
		     + ( ((*it2d)->getCenter_node())(1,2)-pt(1,2))*( ((*it2d)->getCenter_node())(1,2)-pt(1,2)));

	if(it2d==Split_ordered.begin()){
	    Dmin = D;
	    itk = it2d;
	}
	else if(Dmin > D){
	    Dmin = D;
	    itk = it2d;
	}
    }

    matrix pt1;
    pt1.appendRow( ((*itk)->getOrdered_point()).getRow(1) );
    matrix pt2;
    pt2.appendRow( ((*itk)->getOrdered_point()).getRow(2) );

    matrix p = project_point_to_line_segment(pt1, pt2, pt);

    REAL p_par = -1 + 2*sqrt( (p(1,1)-pt1(1,1))*(p(1,1)-pt1(1,1))
			     +(p(1,2)-pt1(1,2))*(p(1,2)-pt1(1,2)) )
		      / sqrt( (pt2(1,1)-pt1(1,1))*(pt2(1,1)-pt1(1,1))
			     +(pt2(1,2)-pt1(1,2))*(pt2(1,2)-pt1(1,2)) );


    matrix Nl(1,2);
    REAL xi = p_par;
    Nl(1,1) = (1-xi)*0.5;
    Nl(1,2) = (1+xi)*0.5;

    int num_Split_ordered = Split_ordered.size();
    (*itk)->assembly_fluid_XFEM_ndiscont_v2_parallel(numnode9, numsnode92, numnode4, numsnode42,
						     num_Split_ordered, close_DBC2);

    matrix sctrBv_m = (*itk)->getSctrBv_m();

    matrix dsol_part(2,1);	// sctrBv_m is 1x2

    dsol_part(1,1) = dsol(sctrBv_m(1,1), 1);
    dsol_part(2,1) = dsol(sctrBv_m(1,2), 1);

    REAL v = (Nl*dsol_part)(1,1);

    return v;

} // velocity_pt_membrane()


void surfacetension::changeIsnanFoot_points(){

    for(std::list<activepoint*>::iterator it=Active_points.begin(); it!=Active_points.end(); ++it){
	matrix Foot_points = (*it)->getFoot_points();
    	if( isnan_mat(Foot_points) ){
	    std::list<activepoint*>::iterator it_temp = it; it_temp--;
	    matrix Foot_points_temp = (*it_temp)->getFoot_points();
	    (*it)->setFoot_points(Foot_points_temp);
	}
    }

} // changeIsnanFoot_points()


void surfacetension::printFoot_points_T_plot(const char* str) const{

    std::ofstream ofs(str);
    if(!ofs) {
    	std::cout << "stream error!" << std::endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    ofs << std::setw(OWID) << "Foot_points_x"
	<< std::setw(OWID) << "Foot_points_y"
	<< std::setw(OWID) << "T_plot_x"
	<< std::setw(OWID) << "T_plot_y" << std::endl;

    matrix Foot_points;
    matrix T_plot;
    for(std::list<activepoint*>::const_iterator it=Active_points.begin(); it!=Active_points.end(); ++it){
	Foot_points = (*it)->getFoot_points();
	T_plot = (*it)->getT_plot();
	
	ofs << std::setw(OWID) << Foot_points(1,1)
	    << std::setw(OWID) << Foot_points(1,2)
	    << std::setw(OWID) << T_plot(1,1)
	    << std::setw(OWID) << T_plot(1,2) << std::endl;
    }

    ofs.close();

} // printFoot_points_T_plot()


void surfacetension::printOrdered_point(const char* str) const{

    std::ofstream ofs(str);
    if(!ofs) {
    	std::cout << "stream error!" << std::endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    ofs << std::setw(OWID) << "Ordered_point_x"
	<< std::setw(OWID) << "Ordered_point_y" << std::endl;

    matrix Ordered_point;
    std::list<ele2dcoupled*>::const_iterator it2d;
    for(it2d=Split_ordered.begin(); it2d!=Split_ordered.end(); ++it2d){
	
	Ordered_point = (*it2d)->getOrdered_point();
	ofs << std::setw(OWID) << Ordered_point(1,1)
	    << std::setw(OWID) << Ordered_point(1,2) << std::endl;
    }

    it2d--;
    Ordered_point = (*it2d)->getOrdered_point();
    ofs << std::setw(OWID) << Ordered_point(1,1)
        << std::setw(OWID) << Ordered_point(1,2) << std::endl;

    ofs.close();

} // printOrdered_point()


void surfacetension::printAratio(const char* str, const int nAratio, 
		     std::vector<REAL>::const_iterator Aratio_beg, 
		     std::vector<REAL>::const_iterator Aratio_end ) const{

    std::ofstream ofs(str);
    if(!ofs) {
    	std::cout << "stream error!" << std::endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    ofs << std::setw(OWID) << nAratio << std::endl;

    for(std::vector<REAL>::const_iterator it=Aratio_beg; it<Aratio_end; ++it){

	ofs << std::setw(OWID) << (*it) << std::endl;
    }

    ofs.close();

} // printAratio()


void surfacetension::printFigureVariables(const char* file_name, const int it) const{

    char stepsstr[4];
    char stepsfp[50];

    // print node9
    sprintf(stepsstr, "%03d", it); 
    strcpy(stepsfp, file_name); strcat(stepsfp, "_figure_node9_"); strcat(stepsfp, stepsstr);
    printNode9(stepsfp);

    // print element9
    sprintf(stepsstr, "%03d", it); 
    strcpy(stepsfp, file_name); strcat(stepsfp, "_figure_element9_"); strcat(stepsfp, stepsstr);
    printElement9(stepsfp);

    // print markers
    sprintf(stepsstr, "%03d", it); 
    strcpy(stepsfp, file_name); strcat(stepsfp, "_figure_markers_"); strcat(stepsfp, stepsstr);
    printMarkers(stepsfp);

    // print dsol
    sprintf(stepsstr, "%03d", it); 
    strcpy(stepsfp, file_name); strcat(stepsfp, "_figure_dsol_"); strcat(stepsfp, stepsstr);
    printDsol(stepsfp);

    // print Ordered_ponit, xFm2, yFm2, Fm2x, Fm2y, loop over Split_ordered
    sprintf(stepsstr, "%03d", it); 
    strcpy(stepsfp, file_name); strcat(stepsfp, "_figure_Ordered_xFm_"); strcat(stepsfp, stepsstr);
    printOrdered_xFm(stepsfp);

    // print Foot_points, x_proj_m, vel_Foot_points_membrane
    sprintf(stepsstr, "%03d", it); 
    strcpy(stepsfp, file_name); strcat(stepsfp, "_figure_Foot_proj_"); strcat(stepsfp, stepsstr);
    printFoot_proj(stepsfp);

} // printFigureVariables()


void surfacetension::printNode9(const char* str) const {

    std::ofstream ofs(str);
    if(!ofs) {
    	std::cout << "stream error!" << std::endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    ofs << std::setw(OWID) << "node9_x"
	<< std::setw(OWID) << "node9_y" << std::endl;

    for(std::vector<node*>::const_iterator it=node9.begin(); it<node9.end(); ++it){
	ofs << std::setw(OWID) << (*it)->getX()
	    << std::setw(OWID) << (*it)->getY() << std::endl;
    }

    ofs.close();

} // printNode9()


void surfacetension::printElement9(const char* str) const {

    std::ofstream ofs(str);
    if(!ofs) {
    	std::cout << "stream error!" << std::endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    ofs << std::setw(OWID) << "element9_1"
	<< std::setw(OWID) << "element9_2"
	<< std::setw(OWID) << "element9_3" 
	<< std::setw(OWID) << "element9_4"
	<< std::setw(OWID) << "element9_5"
	<< std::setw(OWID) << "element9_6" 
	<< std::setw(OWID) << "element9_7"
	<< std::setw(OWID) << "element9_8"
	<< std::setw(OWID) << "element9_9" << std::endl;

    matrix sctr9;
    for(std::vector<ele2dcoupled*>::const_iterator it2d=element.begin(); it2d<element.end(); ++it2d){
	sctr9 = (*it2d)->getSctr9();
	ofs << std::setw(OWID) << sctr9(1,1)
	    << std::setw(OWID) << sctr9(1,2) 
	    << std::setw(OWID) << sctr9(1,3)
	    << std::setw(OWID) << sctr9(1,4)
	    << std::setw(OWID) << sctr9(1,5) 
	    << std::setw(OWID) << sctr9(1,6)
	    << std::setw(OWID) << sctr9(1,7) 
	    << std::setw(OWID) << sctr9(1,8) 
	    << std::setw(OWID) << sctr9(1,9)  << std::endl;
    }

    ofs.close();

} // printElement9()


void surfacetension::printMarkers(const char* str) const {

    std::ofstream ofs(str);
    if(!ofs) {
    	std::cout << "stream error!" << std::endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    ofs << std::setw(OWID) << "markers_x"
	<< std::setw(OWID) << "markers_y" << std::endl;

    for(std::vector<node*>::const_iterator it=markers.begin(); it<markers.end(); ++it){
	ofs << std::setw(OWID) << (*it)->getX()
	    << std::setw(OWID) << (*it)->getY() << std::endl;
    }

    ofs.close();

} // printMarkers()


void surfacetension::printDsol(const char* str) const {

    std::ofstream ofs(str);
    if(!ofs) {
    	std::cout << "stream error!" << std::endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);

    for(int i=1; i<dsol.num_row+1; ++i){
	ofs << std::setw(OWID) << dsol.getVal(i,1) << std::endl;
    }

    ofs.close();

} // printDsol()


void surfacetension::printOrdered_xFm(const char* str) const{

    std::ofstream ofs(str);
    if(!ofs) {
    	std::cout << "stream error!" << std::endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    ofs << std::setw(OWID) << "Ordered_point_x"
	<< std::setw(OWID) << "Ordered_point_y" 
	<< std::setw(OWID) << "xFm2" 
	<< std::setw(OWID) << "yFm2" 
	<< std::setw(OWID) << "Fm2x" 
	<< std::setw(OWID) << "Fm2y" << std::endl;

    matrix Ordered_point;
    std::list<ele2dcoupled*>::const_iterator it2d;
    for(it2d=Split_ordered.begin(); it2d!=Split_ordered.end(); ++it2d){
	
	Ordered_point = (*it2d)->getOrdered_point();
	ofs << std::setw(OWID) << Ordered_point(1,1)
	    << std::setw(OWID) << Ordered_point(1,2)
	    << std::setw(OWID) << (*it2d)->getxFm2()
	    << std::setw(OWID) << (*it2d)->getyFm2()
	    << std::setw(OWID) << (*it2d)->getFm2x()
	    << std::setw(OWID) << (*it2d)->getFm2y() << std::endl;
    }

    it2d--;
    Ordered_point = (*it2d)->getOrdered_point();
    ofs << std::setw(OWID) << Ordered_point(1,1)
        << std::setw(OWID) << Ordered_point(1,2) << std::endl;

    ofs.close();

} // printOrdered_xFm()


void surfacetension::printFoot_proj(const char* str) const{

    std::ofstream ofs(str);
    if(!ofs) {
    	std::cout << "stream error!" << std::endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    ofs << std::setw(OWID) << "Foot_points_x"
	<< std::setw(OWID) << "Foot_points_y" 
	<< std::setw(OWID) << "vel_F_p_mem_x" 
	<< std::setw(OWID) << "vel_F_p_mem_y" 
	<< std::setw(OWID) << "x_proj_m_x" 
	<< std::setw(OWID) << "x_proj_m_y" << std::endl;

    matrix Foot_points;
    matrix vel_F_p_mem;
    matrix x_proj_m;
    for(std::list<activepoint*>::const_iterator it=Active_points.begin(); it!=Active_points.end(); ++it){
	Foot_points = (*it)->getFoot_points();
	vel_F_p_mem = (*it)->getVel_Foot_points_membrane();
	x_proj_m    = (*it)->getX_proj_m();
	
	ofs << std::setw(OWID) << Foot_points(1,1)
	    << std::setw(OWID) << Foot_points(1,2)
	    << std::setw(OWID) << vel_F_p_mem(1,1)
	    << std::setw(OWID) << vel_F_p_mem(1,2)
	    << std::setw(OWID) << x_proj_m(1,1)
	    << std::setw(OWID) << x_proj_m(1,2) << std::endl;
    }

    ofs.close();

} // printFoot_proj()


void surfacetension::printF_debug(const char* str) const{

    std::ofstream ofs(str);
    if(!ofs) {
    	std::cout << "stream error!" << std::endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);

    ofs << F.print();

    ofs.close();

} // printF_debug()


void surfacetension::printK_debug(const char* str) const{

    std::ofstream ofs(str);
    if(!ofs) {
    	std::cout << "stream error!" << std::endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);

    ofs << K.print();

    ofs.close();

} // printK_debug()


void surfacetension::printDsol_debug(const char* str) const{

    std::ofstream ofs(str);
    if(!ofs) {
    	std::cout << "stream error!" << std::endl; exit(-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);

    ofs << dsol.print();

    ofs.close();

} // printDsol_debug()


void surfacetension::readDsol_debug(){

    std::ifstream ifs("dsol_it2_mat");
    if(!ifs) {
	std::cout << "stream error!" << std::endl; exit(-1);
    }
    dsol = zeros(total_unknown,1);
    for(int i=1; i<total_unknown+1; ++i){
	ifs >> dsol(i,1);
    }


} // readDsol_debug()




} // end of memFluid
