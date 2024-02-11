#include "ele2dcoupled.h"
#include <complex>


namespace memFluid{

void ele2dcoupled::init(){

    element_length = 1;

    center_node = ele9_node[8];
    enrich_element2 = false;
    split_elem2 = false;
    W_update = false;
    border_elements = false;

    Split_ordered = -1;	// initialized when vector<ele2dcoupled*> Split_ordered is created

    // get level-set/elements intersection points
    xc2 = zeros(1,2);	// 1x2
    yc2 = zeros(1,2);	// 1x2
    xcp2 = zeros(1,2);	// 1x2
    ycp2 = zeros(1,2);	// 1x2


    Le2 = 0;
    Le2_old = 0;

    Ordered_point = zeros(2,2);	// coords of two ordered points in the split element

    // Gauss ponits & weights
    W = zeros(1,28);
    Qx = zeros(1,28);
    Qy = zeros(1,28);
    ngp = 0;
    nogp2 = 0;
    Wl2 = zeros(1,3);
    Ql_1D = zeros(1,3);
    Ql_1D2=zeros(1,3);
    Ql_2Dx = zeros(1,3);
    Ql_2Dy = zeros(1,3);
    Paired_connect2 = zeros(1,3);


    t_element_total = zeros(1,18);
    dtdx_total = zeros(1,36);
    t_element = zeros(9,2);
    dtdx = zeros(9,4);

    element_unknowns = 0;
    ntriplets = 0;
    sctrBv = zeros(1,18);
    sctrBl_plus = zeros(1,9);
    sctrBl_minus = zeros(1,9);
    sctrBv_m = zeros(1,9);
    sctrBTOT = zeros(1,43);
    sctrBp = zeros(1,4);
    sctrBlp = zeros(1,4);
    nn9 = 0;
    nn4 = 0;

    Ii = zeros(1,43*43);
    Jj = zeros(1,43*43);
    Xx = zeros(1,43*43);
    Felem = zeros(1,43);
    connect = zeros(1,43);

    v_unknown = 0;
    p_unknown = 0;
    lp_unknown = 0;
    l_plus_unknown = 0;
    l_minus_unknown = 0;
    v_m_unknown = 0;

/*
    matrix indexv;
    matrix indexp;
    matrix indexlp;
    matrix indexl_plus;
    matrix indexl_minus;
    matrix indexv_m;

    matrix Ke;
    matrix Fe;
    matrix se;
    matrix ve;

    // all Ks
    matrix K_vv;
    matrix K_vp;
    matrix K_pv;
    matrix K_lpp;
    matrix K_plp;
    matrix K_l_plusv;
    matrix K_vl_plus;
    matrix K_l_minusv;
    matrix K_vl_minus;
    matrix K_v_mv;
    matrix K_v_mv_m;
    matrix K_vv_m;
    matrix K_l_plusv_m;
    matrix K_v_ml_plus;
    matrix K_l_minusv_m;
    matrix K_v_ml_minus;
    
    // all Fs
    matrix Fp;
    matrix Fv;
    matrix Flp;
    matrix Fl_plus;
    matrix Fl_minus;
    matrix Fv_m;
*/

    // loop over gauss points
    N_9 = zeros(9,1);
    dNdxi_9 = zeros(9,2);
    N_4 = zeros(4,1);
    dNdxi_4 = zeros(4,2);

    // compute the jabocian, old values of v,F and J
    J0 = zeros(2,2);
    J09 = zeros(2,2);

/*
    matrix N9;
    matrix N4;
    matrix B9;

    matrix Nv;
    matrix Bv;
    matrix B_dot_v;
    matrix Np;
    matrix Nv_tangent;
    matrix Bv_tangent;

    // split element
    matrix Nl;
    matrix normal_proj_matrix;
    matrix Gpt_proj;
*/

    N_1D = zeros(2,1);
    t = zeros(1,2);
    force_repulsion = zeros(1,2);
    n_Gpt = zeros(1,2);
    force_membrane = zeros(1,2);
    Fm2x = 0;
    Fm2y = 0;
    xFm2 = 0;
    yFm2 = 0;

} // init()


ele2dcoupled::ele2dcoupled(std::vector<node*>::const_iterator node_beg, 
		 	   std::vector<node*>::const_iterator node_end, int num){

    num_ele = num;
    int i=0;
    for(std::vector<node*>::const_iterator it=node_beg; it<node_end; ++it, ++i){
	ele9_node[i] = (*it);
	if(i<4){
	    ele4_node[i] = (*it);
	}
    }

    init();

} // ele2dcoupled()


matrix ele2dcoupled::getCoord9() const{

    matrix coord(9,2);
    for(int i=0; i<9; ++i){
	coord(i+1, 1) = ele9_node[i]->nodex;
	coord(i+1, 2) = ele9_node[i]->nodey;
    }
    
    return coord;

} // getCoord9()


matrix ele2dcoupled::getSctr9() const{

    matrix sctr9(1,9);
    for(int i=0; i<9; i++){
	sctr9(1,i+1) = ele9_node[i]->num_node9;
    }

    return sctr9;

} // getSctr9()


matrix ele2dcoupled::getSctr4() const{

    matrix sctr4(1,4);
    for(int i=0; i<4; i++){
	sctr4(1,i+1) = ele4_node[i]->num_node4;
    }

    return sctr4;

} // getSctr4()


matrix ele2dcoupled::getPhi91() const{
    matrix Phi91(1,9);
    for(int i=0; i<9; ++i)
   	Phi91(1,i+1) = ele9_node[i]->getPhi91();

    return Phi91;

} // getPhi91()


matrix ele2dcoupled::getPhi92() const{
    matrix Phi92(1,9);
    for(int i=0; i<9; ++i)
   	Phi92(1,i+1) = ele9_node[i]->getPhi92();

    return Phi92;

} // getPhi92()


matrix ele2dcoupled::getPhi41() const{
    matrix Phi41(1,4);
    for(int i=0; i<4; ++i)
   	Phi41(1,i+1) = ele4_node[i]->getPhi41();

    return Phi41;

} // getPhi41()


matrix ele2dcoupled::getPhi42() const{
    matrix Phi42(1,4);
    for(int i=0; i<4; ++i)
   	Phi42(1,i+1) = ele4_node[i]->getPhi42();

    return Phi42;

} // getPhi42()


matrix ele2dcoupled::getGrad_phi_x() const{
    matrix Grad_phi_x(9,1);
    for(int i=0; i<9; ++i)
   	Grad_phi_x(i+1,1) = ele9_node[i]->getGrad_phi_x();

    return Grad_phi_x;

} // getGrad_phi_x()


matrix ele2dcoupled::getGrad_phi_y() const{
    matrix Grad_phi_y(9,1);
    for(int i=0; i<9; ++i)
   	Grad_phi_y(i+1,1) = ele9_node[i]->getGrad_phi_y();

    return Grad_phi_y;

} // getGrad_phi_y()


matrix ele2dcoupled::getNodeCoord9() const{

    matrix coord(9,2);
    for(int i=0; i<9; ++i){
	coord(i,1) = ele9_node[i]->nodex;
	coord(i,2) = ele9_node[i]->nodey;
    }

    return coord;

} // getNodeCoord9()


matrix ele2dcoupled::getNodeCoord4() const{

    matrix coord(4,2);
    for(int i=0; i<4; ++i){
	coord(i,1) = ele4_node[i]->nodex;
	coord(i,2) = ele4_node[i]->nodey;
    }

    return coord;

} // getNodeCoord4()


matrix ele2dcoupled::getCenter_node() const{

    matrix coord(1,2);
    coord(1,1) = center_node->nodex;
    coord(1,2) = center_node->nodey;

    return coord;

} // getCenter_node()


void ele2dcoupled::initNogp2Wl2Ql2Paired_connect2(){

    nogp2=0;
    Wl2=zeros(1,3);
    Ql_1D2=zeros(1,3);
    Ql_1D=zeros(1,3);
    Ql_2Dx=zeros(1,3);
    Ql_2Dy=zeros(1,3);
    Paired_connect2 = zeros(1,3);

} // initNogp2Wl2Ql2Paired_connect2()


void ele2dcoupled::calcLe2(){

    REAL x1_1, x1_2;
    REAL x2_1, x2_2;
    x1_1 = xc2(1,2); 
    x1_2 = yc2(1,2);
    x2_1 = xc2(1,1);
    x2_2 = yc2(1,1);
    // level set length in element
    Le2 = pow( ((x2_1 -x1_1)*(x2_1-x1_1) + (x2_2 -x1_2)*(x2_2-x1_2)), 0.5);

} // calcLe2()


void ele2dcoupled::calcNogp2Wl2Ql_2D(const REAL &dx9, const REAL &dy9){

//    int k = Split_ordered;

    int order = 3;
//    nogp2 = order;
    // find the level-set gauss points coord and weight in 1D parent element
    Find_Gauss_Euler1D();	// order is always 3
    nogp2 = order;

    // find the level-set gauss points coord in 2D parent element
    
    matrix n1(1,2);
    matrix n2(1,2);
    n1(1,1) = (Ordered_point(1,1)-center_node->nodex)/dx9;
    n1(1,2) = (Ordered_point(1,2)-center_node->nodey)/dy9;
    n2(1,1) = (Ordered_point(2,1)-center_node->nodex)/dx9;
    n2(1,2) = (Ordered_point(2,2)-center_node->nodey)/dy9;

    matrix nA(1,2);
    matrix nD(1,2);
    nA(1,1) = 0.5*(n1(1,1)+n2(1,1));
    nA(1,2) = 0.5*(n1(1,2)+n2(1,2));
    nD(1,1) = 0.5*(n1(1,1)-n2(1,1));
    nD(1,2) = 0.5*(n1(1,2)-n2(1,2));

    Ql_2Dx(1,1) = nA(1,1) + Ql_1D(1,1)*nD(1,1);
    Ql_2Dx(1,2) = nA(1,1) + Ql_1D(1,2)*nD(1,1);
    Ql_2Dx(1,3) = nA(1,1) + Ql_1D(1,3)*nD(1,1);
    Ql_2Dy(1,1) = nA(1,2) + Ql_1D(1,1)*nD(1,2);
    Ql_2Dy(1,2) = nA(1,2) + Ql_1D(1,2)*nD(1,2);
    Ql_2Dy(1,3) = nA(1,2) + Ql_1D(1,3)*nD(1,2);

} // calcNogp2Wl2Ql_2D()


int ele2dcoupled::getPartial_enriched() const {
	
    int sum = 0;
    for(int i=0; i<9; ++i){
	if( ele9_node[i]->enrich_node92 )
	    sum++;
    }

    return sum;

} // getPartial_enriched()


int ele2dcoupled::getPartial_enriched4() const {
	
    int sum = 0;
    for(int i=0; i<4; ++i){
	if( ele9_node[i]->enrich_node92 )
	    sum++;
    }

    return sum;

} // getPartial_enriched4()


void ele2dcoupled::calcSctrBTOT(){

    sctrBTOT = zeros(1,43);
    if(element_unknowns > sctrBTOT.num_col)
	sctrBTOT = zeros(1, element_unknowns);
    for(int i=1; i<sctrBv.num_col+1; ++i){
	if(i>element_unknowns) break;
	sctrBTOT(1, i) = sctrBv(1,i);
    }
    for(int i=1; i<sctrBp.num_col+1; ++i){
	if(i+sctrBv.num_col>element_unknowns) break;
	sctrBTOT(1, sctrBv.num_col+i) = sctrBp(1,i);
    }
    for(int i=1; i<sctrBlp.num_col+1; ++i){
	if(sctrBv.num_col+sctrBp.num_col+i > element_unknowns) break;
	sctrBTOT(1, sctrBv.num_col+sctrBp.num_col+i) = sctrBlp(1,i);
    }
    for(int i=1; i<sctrBl_plus.num_col+1; ++i){
	if(sctrBv.num_col+sctrBp.num_col+sctrBlp.num_col+i>element_unknowns) break;
	sctrBTOT(1, sctrBv.num_col+sctrBp.num_col+sctrBlp.num_col+i) = sctrBl_plus(1,i);
    }
    for(int i=1; i<sctrBl_minus.num_col+1; ++i){
	if(sctrBv.num_col+sctrBp.num_col+sctrBlp.num_col+sctrBl_plus.num_col+i>element_unknowns) break;
	sctrBTOT(1, sctrBv.num_col+sctrBp.num_col+sctrBlp.num_col+sctrBl_plus.num_col+i) = sctrBl_minus(1,i);
    }
    for(int i=1; i<sctrBv_m.num_col+1; ++i){
	if(sctrBv.num_col+sctrBp.num_col+sctrBlp.num_col+sctrBl_plus.num_col+sctrBl_minus.num_col+i>element_unknowns) break;
	sctrBTOT(1, sctrBv.num_col+sctrBp.num_col+sctrBlp.num_col+sctrBl_plus.num_col+sctrBl_minus.num_col+i) = sctrBv_m(1,i);
    }

} // calcSctrBTOT()


void ele2dcoupled::setAllUnknowns(){

    v_unknown = size(sctrBv,2);
    p_unknown = size(sctrBp,2);
    lp_unknown = size(sctrBlp,2);
    l_plus_unknown = size(sctrBl_plus,2);
    l_minus_unknown = size(sctrBl_minus,2);
    v_m_unknown = size(sctrBv_m,2);

} // end of setAllUnknowns()


void ele2dcoupled::setAllIndex(){
    // indexv
    indexv.clear();
    for(int i=1; i<v_unknown+1; ++i){
	std::vector<REAL> temp;
	temp.push_back(i);	
	indexv.appendCol(temp);	// row vector
	temp.clear();
    }
    // indexp
    indexp.clear();
    for(int i=v_unknown+1; i<p_unknown+v_unknown+1; ++i){
	std::vector<REAL> temp;
	temp.push_back(i);	
	indexp.appendCol(temp);	// row vector
	temp.clear();
    }
    // indexlp
    indexlp.clear();
    for(int i=p_unknown+v_unknown + 1; i<p_unknown+v_unknown + lp_unknown+1; ++i){
	std::vector<REAL> temp;
	temp.push_back(i);	
	indexlp.appendCol(temp);	// row vector
	temp.clear();
    }
    // indexl_plus
    indexl_plus.clear();
    for(int i=p_unknown+v_unknown+ lp_unknown + 1; i<p_unknown+v_unknown+ lp_unknown + l_plus_unknown+1; ++i){
	std::vector<REAL> temp;
	temp.push_back(i);	
	indexl_plus.appendCol(temp);	// row vector
	temp.clear();
    }
    // indexl_minus
    indexl_minus.clear();
    for(int i=p_unknown+v_unknown+ lp_unknown + l_plus_unknown + 1; i<p_unknown+v_unknown+ lp_unknown + l_plus_unknown + l_minus_unknown+1; ++i){
	std::vector<REAL> temp;
	temp.push_back(i);	
	indexl_minus.appendCol(temp);	// row vector
	temp.clear();
    }
    // indexv_m
    indexv_m.clear();
    for(int i=p_unknown+v_unknown+ lp_unknown + l_plus_unknown + l_minus_unknown + 1; i<p_unknown+v_unknown+ lp_unknown + l_plus_unknown + l_minus_unknown + v_m_unknown+1; ++i){
	std::vector<REAL> temp;
	temp.push_back(i);	
	indexv_m.appendCol(temp);	// row vector
	temp.clear();
    }

} // end of setAllIndex()


void ele2dcoupled::initSeVe(matrix&s, matrix&v){

    se = zeros(sctrBp.num_col, 1);
    for(int i=1; i!=sctrBp.num_col+1; ++i){
	se(i,1) = s(sctrBp(1,i), 1);
    }

    ve = zeros(sctrBv.num_col, 1);
    for(int i=1; i!=sctrBv.num_col+1; ++i){
 	ve(i,1) = v(sctrBv(1,i), 1);
    }

} // initSeVe()


void ele2dcoupled::initAllKs(){

    if(length(indexv)==0){
	K_vv.clear();
    }
    else {
    	K_vv = zeros(length(indexv),length(indexv));
    }

    if(length(indexv)==0 || length(indexp)==0){
	K_vp.clear();
    }
    else{
    	K_vp = zeros(length(indexv),length(indexp));
    }

    if(length(indexp)==0 || length(indexv)==0){
	K_pv.clear();
    }
    else{
    	K_pv = zeros(length(indexp),length(indexv));
    }

    if(length(indexlp)==0 || length(indexp)==0){
	K_lpp.clear();
    }
    else{
    	K_lpp = zeros(length(indexlp),length(indexp));
    }

    if(length(indexp)==0 || length(indexlp)==0){
	K_plp.clear();
    }
    else {
    	K_plp = zeros(length(indexp),length(indexlp));
    }

    if(length(indexl_plus)==0 || length(indexv)==0){
   	K_l_plusv.clear();
    }
    else{
    	K_l_plusv = zeros(length(indexl_plus),length(indexv));
    }

    if(length(indexv)==0 || length(indexl_plus)==0){
	K_vl_plus.clear();
    }
    else{
    	K_vl_plus = zeros(length(indexv),length(indexl_plus));
    }

    if(length(indexl_minus)==0 || length(indexv)==0){
	K_l_minusv.clear();
    }
    else
    	K_l_minusv = zeros(length(indexl_minus),length(indexv));
  
    if(length(indexv)==0 || length(indexl_minus)==0)
	K_vl_minus.clear();
    else
    	K_vl_minus = zeros(length(indexv),length(indexl_minus));
    
    if(length(indexv_m)==0 || length(indexv)==0)
	K_v_mv.clear();
    else
	K_v_mv = zeros(length(indexv_m),length(indexv));
    
    if(length(indexv_m)==0 || length(indexv_m)==0)
	K_v_mv_m.clear();
    else
	K_v_mv_m = zeros(length(indexv_m),length(indexv_m));
    
    if(length(indexv)==0 || length(indexv_m)==0)
	K_vv_m.clear();
    else
	K_vv_m = zeros(length(indexv),length(indexv_m));
    
    if(length(indexl_plus)==0 || length(indexv_m)==0)
	K_l_plusv_m.clear();
    else
	K_l_plusv_m = zeros(length(indexl_plus),length(indexv_m));

    if(length(indexv_m)==0 || length(indexl_plus)==0)
   	K_v_ml_plus.clear();
    else
    	K_v_ml_plus = zeros(length(indexv_m),length(indexl_plus));

    if(length(indexl_minus)==0 || length(indexv_m)==0)
	K_l_minusv_m.clear();
    else
    	K_l_minusv_m = zeros(length(indexl_minus),length(indexv_m));

    if(length(indexv_m)==0 || length(indexl_minus)==0)
 	K_v_ml_minus.clear();
    else
    	K_v_ml_minus = zeros(length(indexv_m),length(indexl_minus));

} // end of initAllKs()


void ele2dcoupled::initIJXFelemConnect(){

    Ii = zeros(1,43*43);
    Jj = zeros(1,43*43);
    Xx = zeros(1,43*43);
    Felem = zeros(1,43);
    connect = zeros(1,43);

} // end of initIJXFelemConnect()


void ele2dcoupled::initAllFs(){

    if(length(indexp)==0)
	Fp.clear();
    else
    	Fp = zeros(length(indexp),1);
    
    if(length(indexv)==0)
	Fv.clear();
    else
	Fv = zeros(length(indexv),1);
    
    if(length(indexlp)==0)
	Flp.clear();
    else
	Flp = zeros(length(indexlp),1);

    if(length(indexl_plus)==0)
	Fl_plus.clear();
    else
    	Fl_plus = zeros(length(indexl_plus),1);
    
    if(length(indexl_minus)==0)
	Fl_minus.clear();
    else
	Fl_minus = zeros(length(indexl_minus),1);
    
    if(length(indexv_m)==0)
	Fv_m.clear();
    else
     	Fv_m = zeros(length(indexv_m),1);

} // end of initAllFs()



void ele2dcoupled::calcN_9(const REAL x, const REAL y){

    REAL xi = x; REAL eta = y;

    N_9(1,1) = 0.25*xi*eta*(xi-1)*(eta-1);
    N_9(2,1) = 0.25*xi*eta*(xi+1)*(eta-1);
    N_9(3,1) = 0.25*xi*eta*(xi+1)*(eta+1);
    N_9(4,1) = 0.25*xi*eta*(xi-1)*(eta+1);
    N_9(5,1) = -2*0.25*eta*(xi+1)*(xi-1)*(eta-1);
    N_9(6,1) = -2*0.25*xi*(xi+1)*(eta+1)*(eta-1);
    N_9(7,1) = -2*0.25*eta*(xi+1)*(xi-1)*(eta+1);
    N_9(8,1) = -2*0.25*xi*(xi-1)*(eta+1)*(eta-1);
    N_9(9,1) = 4*0.25*(xi+1)*(xi-1)*(eta+1)*(eta-1);

} // calcN_9()


void ele2dcoupled::calcDNdxi_9(const REAL x, const REAL y){
      
    REAL xi = x; REAL eta = y;

    dNdxi_9(1,1) = 0.25*eta*(2*xi-1)*(eta-1);
    dNdxi_9(1,2) = 0.25*xi*(xi-1)*(2*eta-1);

    dNdxi_9(2,1) = 0.25*eta*(2*xi+1)*(eta-1);
    dNdxi_9(2,2) = 0.25*xi*(xi+1)*(2*eta-1);
                 
    dNdxi_9(3,1) = 0.25*eta*(2*xi+1)*(eta+1);
    dNdxi_9(3,2) = 0.25*xi*(xi+1)*(2*eta+1);
                 
    dNdxi_9(4,1) = 0.25*eta*(2*xi-1)*(eta+1);
    dNdxi_9(4,2) = 0.25*xi*(xi-1)*(2*eta+1);
                
    dNdxi_9(5,1) = -4*0.25*xi*eta*(eta-1);
    dNdxi_9(5,2) = -2*0.25*(xi+1)*(xi-1)*(2*eta-1);
         
    dNdxi_9(6,1) = -2*0.25*(2*xi+1)*(eta+1)*(eta-1);
    dNdxi_9(6,2) = -4*0.25*xi*eta*(xi+1);
                
    dNdxi_9(7,1) = -4*0.25*xi*eta*(eta+1);
    dNdxi_9(7,2) = -2*0.25*(xi+1)*(xi-1)*(2*eta+1);
         
    dNdxi_9(8,1) = -2*0.25*(2*xi-1)*(eta+1)*(eta-1);
    dNdxi_9(8,2) = -4*0.25*xi*eta*(xi-1);
                 
    dNdxi_9(9,1) = 8*0.25*xi*(eta*eta-1);
    dNdxi_9(9,2) = 8*0.25*eta*(xi*xi-1);

} // calcDNdxi_9()


void ele2dcoupled::calcN_4(const REAL x, const REAL y){

    REAL xi = x; REAL eta = y;

    N_4(1,1) = 0.25*(1-xi)*(1-eta);
    N_4(2,1) = 0.25*(1+xi)*(1-eta);
    N_4(3,1) = 0.25*(1+xi)*(1+eta);
    N_4(4,1) = 0.25*(1-xi)*(1+eta);

} // calcN_4()


void ele2dcoupled::calcDNdxi_4(const REAL x, const REAL y){

    REAL xi = x; REAL eta = y;      

    dNdxi_4(1,1) = 0.25*(-(1-eta));
    dNdxi_4(1,2) = 0.25*(-(1-xi));

    dNdxi_4(2,1) = 0.25*(1-eta);
    dNdxi_4(2,2) = 0.25*(-(1+xi));
                 
    dNdxi_4(3,1) = 0.25*(1+eta);
    dNdxi_4(3,2) = 0.25*(1+xi);
                 
    dNdxi_4(4,1) = 0.25*(-(1+eta));
    dNdxi_4(4,2) = 0.25*(1-xi);
                
} // calcDNdxi_4()


void ele2dcoupled::calcN_1D(const REAL pt1D){

    N_1D = zeros(2,1);
    REAL xi = pt1D;
    N_1D(1,1) = (1-xi)*0.5;
    N_1D(2,1) = (1+xi)*0.5;

} // calcN_1D()


void ele2dcoupled::calculateKsFs(const REAL Lx, const REAL Ly, const std::string &plain_axi, 
std::vector<REAL>::const_iterator nodex_beg, std::vector<REAL>::const_iterator nodex_end, 	      std::vector<REAL>::const_iterator nodey_beg, std::vector<REAL>::const_iterator nodey_end){

    for(int ix=1; ix!=ngp+1; ++ix){
	REAL ptx = Qx(1,ix);
	REAL pty = Qy(1,ix);

	//lagrange_basis("Q9", ptx, pty);	// calculate N_9 and dNdxi_9
	calcN_9(ptx, pty);
	calcDNdxi_9(ptx, pty);

        matrix node_coord9 = getCoord9();
	matrix Gpt = N_9.getTrans()*node_coord9;
		
	REAL r = Gpt(1,1)+Lx*0.5;	// Lx is extern

	int in_on;
	//in_on = InPolygon2(Gpt(1,1), Gpt(1,2), nodex_beg, nodex_end, nodey_beg, nodey_end);
        in_on = inpoly(Gpt(1,1), Gpt(1,2), nodex_beg, nodex_end, nodey_beg, nodey_end);	

	int phi, BCjump, g;
	REAL mu;
	if(in_on==1){
	    phi = -1;
 	    BCjump = 1;
	    g = 10;
	    mu = MU_F_MINUS;
	}
	else{
            phi=1;
            BCjump = 2;
            g= 0;
            mu = MU_F_PLUS;
	}
		
	//-----------------------------------
	// compute the jacobian, old values of v,F and J at guass point
	// and get Nf, Bf, Nv, Bv, B.v, NJ, BJ
	xfem_Nmu9_parallel(BCjump);	// calculate ele2dcoupled.N9

	//lagrange_basis("Q4", ptx, pty);	// calculate N_4, dNdxi_4
	calcN_4(ptx, pty);
	calcDNdxi_4(ptx, pty);

	xfem_Nmu4_parallel(BCjump);	// calculate ele2dcoupled.N4
	
	xfem_Bmu9_parallel(BCjump);	// calculate ele2dcoupled.B9 , J09

	//setJ0EqualJ09();			// set J0 = J09;
	J0 = J09;

	// Nv Bv B_dot_v Np Nv_tangent Bv_tangent in ele2dcoupled
	get_NvBvNp_fluid_tdisc(r, plain_axi);	

	//--------------------------------------//
	//   compute Cinte, Kinite, Finite      //
 	//--------------------------------------//
	if(plain_axi=="plain"){

		K_vv = K_vv + Bv.getTrans()*Bv*mu*W(1,ix)*det(J0);
		
		matrix bigNp;
		bigNp.appendRow(Np.getRow(1));
		bigNp.appendRow(Np.getRow(1));
		std::vector<REAL> zeros_vec;
		for(int i=1; i!=length(Np)+1; ++i)
		    zeros_vec.push_back(0);
		bigNp.appendRow(zeros_vec);
		bigNp.appendRow(zeros_vec);
                K_vp = K_vp - Bv.getTrans()*bigNp*W(1,ix)*det(J0);
                K_pv = K_pv - Np.getTrans()*B_dot_v*W(1,ix)*det(J0);
                Fp = Fp - Np.getTrans()*Np*se*W(1,ix)*det(J0);

		matrix zeros_mat(2,1);
		zeros_mat(1,1) = 0; zeros_mat(2,1) = 0;
                Fv = Fv + Nv.getTrans()*zeros_mat*W(1,ix)*det(J0);

    	} // if "plain"
	
	if(plain_axi=="axi1"){

		K_vv = K_vv + Bv.getTrans()*Bv*mu*W(1,ix)*det(J0)*r;
		matrix bigNp;
		bigNp.appendRow(Np.getRow(1));
		bigNp.appendRow(Np.getRow(1));
		std::vector<REAL> zeros_vec;
		for(int i=1; i!=length(Np)+1; ++i)
		    zeros_vec.push_back(0);
		bigNp.appendRow(zeros_vec);
		bigNp.appendRow(Np.getRow(1));
                K_vp = K_vp - Bv.getTrans()*bigNp*W(1,ix)*det(J0)*r;

                K_pv = K_pv - Np.getTrans()*B_dot_v*W(1,ix)*det(J0)*r;
                Fp = Fp - Np.getTrans()*Np*se*W(1,ix)*det(J0)*r;

		matrix grav(2,1);
		grav(1,1) = 0; grav(2,1) = -RHO*g;
                Fv = Fv + Nv.getTrans()*grav*W(1,ix)*det(J0)*r;

	} // if "axi1"

	if(plain_axi=="axi2"){
		K_vv = K_vv + Bv.getTrans()*Bv*mu*W(1,ix)*det(J0)*r;

		matrix bigNp;
		bigNp.appendRow(Np.getRow(1));
		bigNp.appendRow(Np.getRow(1));
		std::vector<REAL> zeros_vec;
		for(int i=1; i!=length(Np)+1; ++i)
		    zeros_vec.push_back(0);
		bigNp.appendRow(zeros_vec);
		bigNp.appendRow(Np.getRow(1));
                K_vp = K_vp - Bv.getTrans()*bigNp*W(1,ix)*det(J0)*r;
                K_pv = K_pv - Np.getTrans()*B_dot_v*W(1,ix)*det(J0)*r;
                Fp = Fp - Np.getTrans()*Np*se*W(1,ix)*det(J0)*r;

		matrix zeros_mat(2,1);
		zeros_mat(1,1) = 0; zeros_mat(2,1) = 0;
                Fv = Fv + Nv.getTrans()*zeros_mat*W(1,ix)*det(J0)*r;

	} // if "axi2"
		
    }// end of gauss points for loop

} // endo of calculateKsFs()


void ele2dcoupled::xfem_Nmu9_parallel(const int BCjump){

    int nn = 9;
    matrix phi91 = getPhi91();
    matrix phi92 = getPhi92();

    matrix Nfem = N_9;
    matrix Nxfem1, Nxfem2;

    // switch between non-enriched and enriched elements
    if( ele9_node[0]->enrich_node91==false && ele9_node[1]->enrich_node91==false && ele9_node[2]->enrich_node91==false
     && ele9_node[3]->enrich_node91==false && ele9_node[4]->enrich_node91==false && ele9_node[5]->enrich_node91==false 
     && ele9_node[6]->enrich_node91==false && ele9_node[7]->enrich_node91==false && ele9_node[8]->enrich_node91==false ){ 
     // Non-enriched elements

	Nxfem1.clear();
    }
    else{
	// loop on nodes, check node9 is enriched
	for(int in=0; in<nn; ++in){
	    if(ele9_node[in]->enrich_node91 == true){ 	// H(x) enriched node9
		REAL Hgp;
		// Enrichment function, H(x) at global Gauss point
		if(BCjump==1){ // condition inside levelset
		    Hgp = 0;
		}
		else if(BCjump==2){	// condition outside levelset
		    Hgp = 1;
		}
		else{
		    REAL phig = (phi91*N_9)(1,1);
		    REAL dist = phig;
		    Hgp = heaviside(dist);	// heaviside implemented in globfuns.h
		}
		// Enrichment function, H(x) at node9 "in"
		REAL dist = phi91(1,in+1);
		REAL Hi = heaviside(dist);
		// Bxfem at node9 "in"
		std::vector<REAL> NI_enr;
		NI_enr.push_back( N_9(in+1, 1)*(Hgp-Hi) );	// N_9 is 9x1
		
		Nxfem1.appendRow(NI_enr);
		NI_enr.clear();

	    } // end if
	} // end for
    } // end else

    if( ele9_node[0]->enrich_node92==false && ele9_node[1]->enrich_node92==false && ele9_node[2]->enrich_node92==false
     && ele9_node[3]->enrich_node92==false && ele9_node[4]->enrich_node92==false && ele9_node[5]->enrich_node92==false 
     && ele9_node[6]->enrich_node92==false && ele9_node[7]->enrich_node92==false && ele9_node[8]->enrich_node92==false ){ 
     // Non-enriched elements

	Nxfem2.clear();
    }
    else{
	// loop on nodes, check node9 is enriched
	for(int in=0; in<nn; ++in){
	    if(ele9_node[in]->enrich_node92 == true){ 	// H(x) enriched node9
		REAL Hgp;
		// Enrichment function, H(x) at global Gauss point
		if(BCjump==1){ // condition inside levelset
		    Hgp = 0;
		}
		else if(BCjump==2){	// condition outside levelset
		    Hgp = 1;
		}
		else{
		    REAL phig = (phi92*N_9)(1,1);
		    REAL dist = phig;
		    Hgp = heaviside(dist);	// heaviside implemented in globfuns.h
		}
		// Enrichment function, H(x) at node9 "in"
		REAL dist = phi92(1,in+1);
		REAL Hi = heaviside(dist);
		// Bxfem at node9 "in"
		std::vector<REAL> NI_enr;
		NI_enr.push_back( N_9(in+1, 1)*(Hgp-Hi) );	// N_9 is 9x1
		
		Nxfem1.appendRow(NI_enr);
		NI_enr.clear();

	    } // end if
	} // end for
    } // end else

    N9 = Nfem;
    for(int i=1; i<Nxfem1.num_row+1; ++i)
	N9.appendRow(Nxfem1.getRow(i));
    for(int i=1; i<Nxfem2.num_row+1; ++i)
	N9.appendRow(Nxfem2.getRow(i));


} // xfem_Nmu9_parallel()


void ele2dcoupled::xfem_Nmu4_parallel(const int BCjump){

    int nn = 4;
    matrix phi41 = getPhi41();
    matrix phi42 = getPhi42();

    matrix Nfem = N_4;
    matrix Nxfem1, Nxfem2;

    // switch between non-enriched and enriched elements
    if( ele4_node[0]->enrich_node41==false && ele4_node[1]->enrich_node41==false
     && ele4_node[2]->enrich_node41==false && ele4_node[3]->enrich_node41==false ){ 
     // Non-enriched elements

	Nxfem1.clear();
    }
    else{
	// loop on nodes, check node4 is enriched
	for(int in=0; in<nn; ++in){
	    if(ele4_node[in]->enrich_node41 == true){ 	// H(x) enriched node4
		REAL Hgp;
		// Enrichment function, H(x) at global Gauss point
		if(BCjump==1){ // condition inside levelset
		    Hgp = 0;
		}
		else if(BCjump==2){	// condition outside levelset
		    Hgp = 1;
		}
		else{
		    REAL phig = (phi41*N_4)(1,1);
		    REAL dist = phig;
		    Hgp = heaviside(dist);	// heaviside implemented in globfuns.h
		}
		// Enrichment function, H(x) at node4 "in"
		REAL dist = phi41(1,in+1);
		REAL Hi = heaviside(dist);
		// Bxfem at node4 "in"
		std::vector<REAL> NI_enr;
		NI_enr.push_back( N_4(in+1, 1)*(Hgp-Hi) );	// N_4 is 4x1
		
		Nxfem1.appendRow(NI_enr);
		NI_enr.clear();

	    } // end if
	} // end for
    } // end else

    // switch between non-enriched and enriched elements
    if( ele4_node[0]->enrich_node42==false && ele4_node[1]->enrich_node42==false 
     && ele4_node[2]->enrich_node42==false && ele4_node[3]->enrich_node42==false ){ 
     // Non-enriched elements

	Nxfem2.clear();
    }
    else{
	// loop on nodes, check node4 is enriched
	for(int in=0; in<nn; ++in){
	    if(ele4_node[in]->enrich_node42 == true){ 	// H(x) enriched node4
		REAL Hgp;
		// Enrichment function, H(x) at global Gauss point
		if(BCjump==1){ // condition inside levelset
		    Hgp = 0;
		}
		else if(BCjump==2){	// condition outside levelset
		    Hgp = 1;
		}
		else{
		    REAL phig = (phi42*N_4)(1,1);
		    REAL dist = phig;
		    Hgp = heaviside(dist);	// heaviside implemented in globfuns.h
		}
		// Enrichment function, H(x) at node4 "in"
		REAL dist = phi42(1,in+1);
		REAL Hi = heaviside(dist);
		// Bxfem at node4 "in"
		std::vector<REAL> NI_enr;
		NI_enr.push_back( N_4(in+1, 1)*(Hgp-Hi) );	// N_4 is 4x1
		
		Nxfem1.appendRow(NI_enr);
		NI_enr.clear();

	    } // end if
	} // end for
    } // end else

    N4 = Nfem;
    for(int i=1; i<Nxfem1.num_row+1; ++i)
	N4.appendRow(Nxfem1.getRow(i));
    for(int i=1; i<Nxfem2.num_row+1; ++i)
	N4.appendRow(Nxfem2.getRow(i));


} // xfem_Nmu4_parallel()


void ele2dcoupled::xfem_Bmu9_parallel(const int BCjump){

    int nn = 9;
    matrix phi91 = getPhi91();
    matrix phi92 = getPhi92();

    J09 = getCoord9().getTrans()*dNdxi_9;
    matrix invJ09 = J09.getInvs();	// J09 is 2x2

    matrix dNdx = dNdxi_9*invJ09;

    // Bfem is always computed
    matrix Bfem;	// Bfem is 2x9
    Bfem.appendRow(dNdx.getCol(1));
    Bfem.appendRow(dNdx.getCol(2));

    matrix Bxfem1, Bxfem2;
   
    // switch between non-enriched and enriched elements
    if( ele9_node[0]->enrich_node91==false && ele9_node[1]->enrich_node91==false && ele9_node[2]->enrich_node91==false
     && ele9_node[3]->enrich_node91==false && ele9_node[4]->enrich_node91==false && ele9_node[5]->enrich_node91==false 
     && ele9_node[6]->enrich_node91==false && ele9_node[7]->enrich_node91==false && ele9_node[8]->enrich_node91==false ){ 
     // Non-enriched elements

	Bxfem1.clear();
    }
    else{
//	int k = 0;
	// loop on nodes, check node9 is enriched
	for(int in=0; in<nn; ++in){
	    if(ele9_node[in]->enrich_node91 == true){ 	// H(x) enriched node9
		REAL Hgp;
		// Enrichment function, H(x) at global Gauss point
		if(BCjump==1){ // condition inside levelset
		    Hgp = 0;
		}
		else if(BCjump==2){	// condition outside levelset
		    Hgp = 1;
		}
		else{
		    REAL phig = (phi91*N_9)(1,1);
		    REAL dist = phig;
		    Hgp = heaviside(dist);	// heaviside implemented in globfuns.h
		}
		// Enrichment function, H(x) at node9 "in"
		REAL dist = phi91(1,in+1);
		REAL Hi = heaviside(dist);
		// Bxfem at node9 "in"
		std::vector<REAL> BI_enr;
		BI_enr.push_back( dNdx(in+1,1)*(Hgp-Hi) );
		BI_enr.push_back( dNdx(in+1,2)*(Hgp-Hi) );
		
		Bxfem1.appendCol(BI_enr);
		BI_enr.clear();

	    } // end if
	} // end for
    } // end else

    if( ele9_node[0]->enrich_node92==false && ele9_node[1]->enrich_node92==false && ele9_node[2]->enrich_node92==false
     && ele9_node[3]->enrich_node92==false && ele9_node[4]->enrich_node92==false && ele9_node[5]->enrich_node92==false
     && ele9_node[6]->enrich_node92==false && ele9_node[7]->enrich_node92==false && ele9_node[8]->enrich_node92==false ){ 
     // Non-enriched elements

	Bxfem2.clear();
    }
    else{
	// loop on nodes, check node9 is enriched
	for(int in=0; in<nn; ++in){
	    if(ele9_node[in]->enrich_node92 == true){ 	// H(x) enriched node9
		REAL Hgp;
		// Enrichment function, H(x) at global Gauss point
		if(BCjump==1){ // condition inside levelset
		    Hgp = 0;
		}
		else if(BCjump==2){	// condition outside levelset
		    Hgp = 1;
		}
		else{
		    REAL phig = (phi92*N_9)(1,1);
		    REAL dist = phig;
		    Hgp = heaviside(dist);	// heaviside implemented in globfuns.h
		}
		// Enrichment function, H(x) at node9 "in"
		REAL dist = phi92(1,in+1);
		REAL Hi = heaviside(dist);
		// Bxfem at node9 "in"
		std::vector<REAL> BI_enr;
		BI_enr.push_back( dNdx(in+1,1)*(Hgp-Hi) );
		BI_enr.push_back( dNdx(in+1,2)*(Hgp-Hi) );
		
		Bxfem2.appendCol(BI_enr);
		BI_enr.clear();

	    } // end if
	} // end for
    } // end else

    B9 = Bfem;
    for(int i=1; i<Bxfem1.num_col+1; ++i)
	B9.appendCol(Bxfem1.getCol(i));
    for(int i=1; i<Bxfem2.num_col+1; ++i)
	B9.appendCol(Bxfem2.getCol(i));

} // xfem_Bmu9_parallel()


void ele2dcoupled::get_NvBvNp_fluid_tdisc(const REAL r, const std::string &plain_axi){

    if(plain_axi=="plain"){
	    // Nv
	    Nv = zeros(2,nn9);
	    for(int i=1; i<10; ++i){
		Nv(1, 2*i-1) = N9(i,1);
		Nv(2, 2*i)   = N9(i,1);
	    }
	    for(int i=19; i<nn9+1; ++i){
		Nv(1,i) = N9(i-9, 1)*t_element(i-18, 1);
		Nv(2,i) = N9(i-9, 1)*t_element(i-18, 2);
	    }

	    // Nv_tangent
	    Nv_tangent = zeros(2,9);
	    for(int i=1; i<10; ++i){
		Nv_tangent(1,i) = N9(i,1)*t_element(i,1);
		Nv_tangent(2,i) = N9(i,1)*t_element(i,2);
	    }

	    // Bv
	    Bv = zeros(4,nn9);
	    for(int i=1; i<10; ++i){
		Bv(1, 2*i-1) = B9(1,i);
		Bv(2, 2*i)   = B9(2,i);
		Bv(3, 2*i-1) = B9(2,i);
		Bv(4, 2*i)   = B9(1,i);
	    }
	    for(int i=19; i<nn9+1; ++i){
		Bv(1,i) = B9(1,i-9)*t_element(i-18,1)+N9(i-9,1)*dtdx(i-18,1);
		Bv(2,i) = B9(2,i-9)*t_element(i-18,2)+N9(i-9,1)*dtdx(i-18,2);
		Bv(3,i) = B9(2,i-9)*t_element(i-18,1)+N9(i-9,1)*dtdx(i-18,3);
		Bv(4,i) = B9(1,i-9)*t_element(i-18,2)+N9(i-9,1)*dtdx(i-18,4);
	    }

	    // Bv_tangent
	    Bv_tangent = zeros(4,9);
	    for(int i=1; i<10; ++i){
		Bv_tangent(1,i) = B9(1,i)*t_element(i,1)+N9(i,1)*dtdx(i,1);
		Bv_tangent(2,i) = B9(2,i)*t_element(i,2)+N9(i,1)*dtdx(i,2);
		Bv_tangent(3,i) = B9(2,i)*t_element(i,1)+N9(i,1)*dtdx(i,3);
		Bv_tangent(4,i) = B9(1,i)*t_element(i,2)+N9(i,1)*dtdx(i,4);
	    }

	    // B_dot_v
	    B_dot_v = zeros(1,nn9);
	    for(int i=1; i<10; ++i){
		B_dot_v(1,2*i-1) = B9(1,i);
	 	B_dot_v(1,2*i)   = B9(2,i);
	    }
	    for(int i=19; i<nn9+1; ++i){
		B_dot_v(1,i) = B9(1,i-9)*t_element(i-18,1)+B9(2,i-9)*t_element(i-18,2)
			     + N9(i-9,1)*(dtdx(i-18,1)+dtdx(i-18,2));
	    }

    } // if "plain"

    if(plain_axi=="axi1"){
	    // Nv
	    Nv = zeros(2,nn9);
	    for(int i=1; i<10; ++i){
		Nv(1, 2*i-1) = N9(i,1);
		Nv(2, 2*i)   = N9(i,1);
	    }
	    for(int i=19; i<nn9+1; ++i){
		Nv(1,i) = N9(i-9, 1)*t_element(i-18, 1);
		Nv(2,i) = N9(i-9, 1)*t_element(i-18, 2);
	    }

	    // Nv_tangent
	    Nv_tangent = zeros(2,9);
	    for(int i=1; i<10; ++i){
		Nv_tangent(1,i) = N9(i,1)*t_element(i,1);
		Nv_tangent(2,i) = N9(i,1)*t_element(i,2);
	    }

	    // Bv
	    Bv = zeros(4,nn9);
	    for(int i=1; i<10; ++i){
		Bv(1, 2*i-1) = B9(1,i);
		Bv(2, 2*i)   = B9(2,i);
		Bv(3, 2*i-1) = 1.0/(pow(2,0.5))*B9(2,i);
		Bv(3, 2*i)   = 1.0/(pow(2,0.5))*B9(1,i);
		Bv(4, 2*i-1) = N9(i,1)/r;
	    }

	    for(int i=19; i<nn9+1; ++i){
		Bv(1,i) = B9(1,i-9)*t_element(i-18,1)+N9(i-9,1)*dtdx(i-18,1);
		Bv(2,i) = B9(2,i-9)*t_element(i-18,2)+N9(i-9,1)*dtdx(i-18,2);
		Bv(3,i) = 1.0/pow(2,0.5)* (B9(2,i-9)*t_element(i-18,1)+N9(i-9,1)*dtdx(i-18,3))
			+ 1.0/pow(2,0.5)* (B9(1,i-9)*t_element(i-18,2)+N9(i-9,1)*dtdx(i-18,4));
		Bv(4,i) = N9(i-9,1)*t_element(i-18,1)/r;
	    }

	    // Bv_tangent
	    Bv_tangent = zeros(4,9);
	    for(int i=1; i<10; ++i){
		Bv_tangent(1,i) = B9(1,i)*t_element(i,1)+N9(i,1)*dtdx(i,1);
		Bv_tangent(2,i) = B9(2,i)*t_element(i,2)+N9(i,1)*dtdx(i,2);
		Bv_tangent(3,i) = 1.0/pow(2,0.5)*(B9(2,i)*t_element(i,1)+N9(i,1)*dtdx(i,3))
				+ 1.0/pow(2,0.5)*(B9(1,i)*t_element(i,2)+N9(i,1)*dtdx(i,4));
		Bv_tangent(4,i) = N9(i,1)*t_element(i,1)/r;
	    }

	    // B_dot_v
	    B_dot_v = zeros(1,nn9);
	    for(int i=1; i<10; ++i){
		B_dot_v(1,2*i-1) = B9(1,i)+N9(i,1)/r;
	 	B_dot_v(1,2*i)   = B9(2,i);
	    }

	    for(int i=19; i<nn9+1; ++i){
		B_dot_v(1,i) = B9(1,i-9)*t_element(i-18,1)+B9(2,i-9)*t_element(i-18,2)
			     + N9(i-9,1)*(dtdx(i-18,1)+dtdx(i-18,2))
			     + N9(i-9,1)*t_element(i-18,1)/r;
	    }

    } // if "axi1"

    if(plain_axi=="axi2"){
	    // Nv
	    Nv = zeros(2,nn9);
	    for(int i=1; i<10; ++i){
		Nv(1, 2*i-1) = N9(i,1);
		Nv(2, 2*i)   = N9(i,1);
	    }
	    for(int i=19; i<nn9+1; ++i){
		Nv(1,i) = N9(i-9, 1)*t_element(i-18, 1);
		Nv(2,i) = N9(i-9, 1)*t_element(i-18, 2);
	    }

	    // Nv_tangent
	    Nv_tangent = zeros(2,9);
	    for(int i=1; i<10; ++i){
		Nv_tangent(1,i) = N9(i,1)*t_element(i,1);
		Nv_tangent(2,i) = N9(i,1)*t_element(i,2);
	    }

	    // Bv
	    Bv = zeros(4,nn9);
	    for(int i=1; i<10; ++i){
		Bv(1, 2*i-1) = B9(1,i);
		Bv(2, 2*i)   = B9(2,i);
		Bv(3, 2*i-1) = 1.0/(pow(2,0.5))*B9(2,i);
		Bv(3, 2*i)   = 1.0/(pow(2,0.5))*B9(1,i);
		Bv(4, 2*i-1) = N9(i,1)/r;
	    }
	    for(int i=19; i<nn9+1; ++i){
		Bv(1,i) = B9(1,i-9)*t_element(i-18,1)+N9(i-9,1)*dtdx(i-18,1);
		Bv(2,i) = B9(2,i-9)*t_element(i-18,2)+N9(i-9,1)*dtdx(i-18,2);
		Bv(3,i) = 1.0/pow(2,0.5)* (B9(2,i-9)*t_element(i-18,1)+N9(i-9,1)*dtdx(i-18,3))
			+ 1.0/pow(2,0.5)* (B9(1,i-9)*t_element(i-18,2)+N9(i-9,1)*dtdx(i-18,4));
		Bv(4,i) = (N9(i-9,1)*t_element(i-18,1)+N9(i-9,1)*t_element(i-18,2))/r;
	    }

	    // Bv_tangent
	    Bv_tangent = zeros(4,9);
	    for(int i=1; i<10; ++i){
		Bv_tangent(1,i) = B9(1,i)*t_element(i,1)+N9(i,1)*dtdx(i,1);
		Bv_tangent(2,i) = B9(2,i)*t_element(i,2)+N9(i,1)*dtdx(i,2);
		Bv_tangent(3,i) = 1.0/pow(2,0.5)*(B9(1,i)*t_element(i,1)+N9(i,1)*dtdx(i,3))
				+ 1.0/pow(2,0.5)*(B9(1,i)*t_element(i,2)+N9(i,1)*dtdx(i,4));
		Bv_tangent(4,i) = (N9(i,1)*t_element(i,1)+N9(i,1)*t_element(i,2))/r;
	    }

	    // B_dot_v
	    B_dot_v = zeros(1,nn9);
	    for(int i=1; i<10; ++i){
		B_dot_v(1,2*i-1) = B9(1,i)+N9(i,1)/r;
	 	B_dot_v(1,2*i)   = B9(2,i);
	    }
	    for(int i=19; i<nn9+1; ++i){
		B_dot_v(1,i) = B9(1,i-9)*t_element(i-18,1)+B9(2,i-9)*t_element(i-18,2)
			     + N9(i-9,1)*(dtdx(i-18,1)+dtdx(i-18,2))
			     + (N9(i-9,1)*t_element(i-18,1)+N9(i-9,1)*t_element(i-18,2))/r;
	    }

    } // if "axi2"

//    Np = zeros(1,nn4);
    Np.clear();
    Np.appendRow(N4.getCol(1));

} // get_NvBvNp_fluid_tdisc()


void ele2dcoupled::calcSplitKsFs(const REAL Lx, const REAL Ly, const std::string &plain_axi,
				 std::list<activepoint*>::const_iterator Active_points_beg,
		       		 std::list<activepoint*>::const_iterator Active_points_end ){
    
    for(int ix=1; ix!=nogp2+1; ++ix){

	REAL pt2Dx = Ql_2Dx(1,ix);
	REAL pt2Dy = Ql_2Dy(1,ix);

	REAL pt1D = -Ql_1D(1,ix);

	calcN_9(pt2Dx, pt2Dy);
	calcDNdxi_9(pt2Dx, pt2Dy);

	matrix node_coord9 = getCoord9();
	matrix Gpt = N_9.getTrans()*node_coord9;	// Gpt is 1x2 
		
	REAL r = Gpt(1,1)+Lx*0.5;	// Lx is extern

	//-------------------------------------
	// and get Nf, Bf, Nv, Bv, B.v, NJ, BJ on the minus side
	int BCjump = 1;

	calcN_1D(pt1D);	// lagrange_basis('L2', pt1D)

	xfem_Nmu9_parallel(BCjump);	// calculate ele2dcoupled.N9

	//lagrange_basis("Q4", ptx, pty);	// calculate N_4, dNdxi_4
	calcN_4(pt2Dx, pt2Dy);
	calcDNdxi_4(pt2Dx, pt2Dy);

	xfem_Nmu4_parallel(BCjump);	// calculate ele2dcoupled.N4
	
	xfem_Bmu9_parallel(BCjump);	// calculate ele2dcoupled.B9 , J09

	// Nv Bv B_dot_v Np Nv_tangent Bv_tangent in ele2dcoupled
	get_NvBvNp_fluid_tdisc(r, plain_axi);



	//--------------------------------
	// get normal tangent, strain, membrane force
	find_proj_normal_driv_pt_parallel(Lx, Ly, Gpt, plain_axi, 
					  Active_points_beg, Active_points_end); // calculate Gpt_proj, t


	REAL distance = (Gpt_proj(1,2)-YP)-0.02;	// YP is extern

	REAL force_repulsion_n;
	if(distance<0){
	    force_repulsion_n = 5e-6*REPULSION1*(distance)*exp(-0.05*REPULSION2*distance); 
	}
	else{
	    force_repulsion_n = 0;
	}
	
	n_Gpt(1,1) = -t(1,2);
    	n_Gpt(1,2) = t(1,1);

	force_repulsion = force_repulsion_n*n_Gpt;

	find_membrane_force_pt_parallel23(Lx, Ly, plain_axi, Active_points_beg, Active_points_end);

	REAL J0_val = Le2*0.5;

	Nl = N_1D.getTrans();

	normal_proj_matrix.clear();
	if(plain_axi=="plain"){
	   
		r=1;
		normal_proj_matrix = zeros(2,4);
		normal_proj_matrix(1,1) = n_Gpt(1,1);
//		normal_proj_matrix(1,2) = 0;
		normal_proj_matrix(1,3) = n_Gpt(1,2);
//		normal_proj_matrix(1,4) = 0;
//		normal_proj_matrix(2,1) = 0;
		normal_proj_matrix(2,2) = n_Gpt(1,2);
//		normal_proj_matrix(2,3) = 0;
		normal_proj_matrix(2,4) = n_Gpt(1,1);

	} // if "plain"

	if(plain_axi=="axi1"){
		normal_proj_matrix = zeros(2,4);
		normal_proj_matrix(1,1) = n_Gpt(1,1);
//		normal_proj_matrix(1,2) = 0;
		normal_proj_matrix(1,3) = n_Gpt(1,2);
//		normal_proj_matrix(1,4) = 0;
//		normal_proj_matrix(2,1) = 0;
		normal_proj_matrix(2,2) = n_Gpt(1,2);
		normal_proj_matrix(2,3) = n_Gpt(1,1);
//		normal_proj_matrix(2,4) = 0;

	} // if "axi1"

	if(plain_axi=="axi2"){
		normal_proj_matrix = zeros(2,4);
		normal_proj_matrix(1,1) = n_Gpt(1,1);
//		normal_proj_matrix(1,2) = 0;
		normal_proj_matrix(1,3) = n_Gpt(1,2);
//		normal_proj_matrix(1,4) = 0;
//		normal_proj_matrix(2,1) = 0;
		normal_proj_matrix(2,2) = n_Gpt(1,2);
		normal_proj_matrix(2,3) = n_Gpt(1,1);
//		normal_proj_matrix(2,4) = 0;	
	} // if "axi2"

	matrix Np_part(1,4);
	Np_part(1,1) = Np(1,1); Np_part(1,2) = Np(1,2);
	Np_part(1,3) = Np(1,3); Np_part(1,4) = Np(1,4);

	matrix temp_mat;
	temp_mat.clear();
	temp_mat.appendCol(K_lpp.getCol(5)); temp_mat.appendCol(K_lpp.getCol(6));
	temp_mat.appendCol(K_lpp.getCol(7)); temp_mat.appendCol(K_lpp.getCol(8));
	temp_mat = temp_mat + Nl.getTrans()*Np_part*Wl2(1,ix)*J0_val*r;
	for(int i=1; i!=K_lpp.num_row+1; ++i){
	    K_lpp(i,5) = temp_mat(i,1);
	    K_lpp(i,6) = temp_mat(i,2);
	    K_lpp(i,7) = temp_mat(i,3);
	    K_lpp(i,8) = temp_mat(i,4);
	}

	temp_mat.clear();
	temp_mat.appendRow(K_plp.getRow(5)); temp_mat.appendRow(K_plp.getRow(6));
	temp_mat.appendRow(K_plp.getRow(7)); temp_mat.appendRow(K_plp.getRow(8));
	temp_mat = temp_mat + Np_part.getTrans()*Nl*Wl2(1,ix)*J0_val*r;
	for(int i=1; i!=K_plp.num_col+1; ++i){
	    K_plp(5,i) = temp_mat(1,i);
	    K_plp(6,i) = temp_mat(2,i);
	    K_plp(7,i) = temp_mat(3,i);
	    K_plp(8,i) = temp_mat(4,i);
	}

	temp_mat.clear();
	temp_mat.appendCol(K_l_plusv.getCol(19)); temp_mat.appendCol(K_l_plusv.getCol(20));
	temp_mat.appendCol(K_l_plusv.getCol(21)); temp_mat.appendCol(K_l_plusv.getCol(22));
	temp_mat.appendCol(K_l_plusv.getCol(23)); temp_mat.appendCol(K_l_plusv.getCol(24));
	temp_mat.appendCol(K_l_plusv.getCol(25)); temp_mat.appendCol(K_l_plusv.getCol(26));
	temp_mat.appendCol(K_l_plusv.getCol(27));
	temp_mat = temp_mat + Nl.getTrans()*(t*MU_F_PLUS*normal_proj_matrix*Bv_tangent
		 + MU_PLUS*t*Nv_tangent)*Wl2(1,ix)*J0_val*r;
	for(int i=1; i!=K_l_plusv.num_row+1; ++i){
	    K_l_plusv(i,19) = temp_mat(i,1);
	    K_l_plusv(i,20) = temp_mat(i,2);
	    K_l_plusv(i,21) = temp_mat(i,3);
	    K_l_plusv(i,22) = temp_mat(i,4);
	    K_l_plusv(i,23) = temp_mat(i,5);
	    K_l_plusv(i,24) = temp_mat(i,6);
	    K_l_plusv(i,25) = temp_mat(i,7);
	    K_l_plusv(i,26) = temp_mat(i,8);	
	    K_l_plusv(i,27) = temp_mat(i,9);
	}

	temp_mat.clear();
	temp_mat.appendRow(K_vl_plus.getRow(19)); temp_mat.appendRow(K_vl_plus.getRow(20));
	temp_mat.appendRow(K_vl_plus.getRow(21)); temp_mat.appendRow(K_vl_plus.getRow(22));
	temp_mat.appendRow(K_vl_plus.getRow(23)); temp_mat.appendRow(K_vl_plus.getRow(24));
	temp_mat.appendRow(K_vl_plus.getRow(25)); temp_mat.appendRow(K_vl_plus.getRow(26));
	temp_mat.appendRow(K_vl_plus.getRow(27));
	temp_mat = temp_mat + (Bv_tangent.getTrans()*MU_F_PLUS*(normal_proj_matrix.getTrans()
                 * t.getTrans()*Nl) 
                 + MU_PLUS*Nv_tangent.getTrans()*(t.getTrans()*Nl))*Wl2(1,ix)*J0_val*r;
	for(int i=1; i!=K_vl_plus.num_col+1; ++i){
	    K_vl_plus(19,i) = temp_mat(1,i);
	    K_vl_plus(20,i) = temp_mat(2,i);
	    K_vl_plus(21,i) = temp_mat(3,i);
	    K_vl_plus(22,i) = temp_mat(4,i);
	    K_vl_plus(23,i) = temp_mat(5,i);
	    K_vl_plus(24,i) = temp_mat(6,i);
	    K_vl_plus(25,i) = temp_mat(7,i);
	    K_vl_plus(26,i) = temp_mat(8,i);	
	    K_vl_plus(27,i) = temp_mat(9,i);
	}

        K_l_plusv = K_l_plusv + Nl.getTrans()*(t*Nv*(-MU_MINUS + MU_PLUS))*Wl2(1,ix)*J0_val*r;
                
        K_vl_plus = K_vl_plus 
                  + Nv.getTrans()*(t.getTrans()*Nl*(-MU_MINUS + MU_PLUS))*Wl2(1,ix)*J0_val*r;
                
        K_l_plusv_m = K_l_plusv_m + Nl.getTrans()*((MU_MINUS - MU_PLUS)*Nl)*Wl2(1,ix)*J0_val*r;
               
        K_v_ml_plus = K_v_ml_plus + Nl.getTrans()*((MU_MINUS - MU_PLUS)*Nl)*Wl2(1,ix)*J0_val*r;

        K_l_minusv = K_l_minusv + Nl.getTrans()*(t  *MU_F_MINUS* normal_proj_matrix*Bv
                   + MU_MINUS*t*Nv )*Wl2(1,ix)*J0_val*r;
                
        K_vl_minus = K_vl_minus + (MU_F_MINUS*Bv.getTrans()
                                 *(normal_proj_matrix.getTrans()*t.getTrans()*Nl)
                   + MU_MINUS*Nv.getTrans()*(t.getTrans()*Nl))*Wl2(1,ix)*J0_val*r;

        K_l_minusv_m = K_l_minusv_m - Nl.getTrans()*(MU_MINUS*Nl)*Wl2(1,ix)*J0_val*r;
                
        K_v_ml_minus = K_v_ml_minus - Nl.getTrans()*(MU_MINUS*Nl)*Wl2(1,ix)*J0_val*r;
                
        K_v_mv_m = K_v_mv_m + Nl.getTrans()*(Nl*(MU_MINUS + MU_PLUS))*Wl2(1,ix)*J0_val*r;
                
        K_v_mv = K_v_mv + Nl.getTrans()*(-1.0*t*Nv*(MU_MINUS + MU_PLUS))*Wl2(1,ix)*J0_val*r;


	temp_mat.clear();
	temp_mat.appendCol(K_v_mv.getCol(19)); temp_mat.appendCol(K_v_mv.getCol(20));
	temp_mat.appendCol(K_v_mv.getCol(21)); temp_mat.appendCol(K_v_mv.getCol(22));
	temp_mat.appendCol(K_v_mv.getCol(23)); temp_mat.appendCol(K_v_mv.getCol(24));
	temp_mat.appendCol(K_v_mv.getCol(25)); temp_mat.appendCol(K_v_mv.getCol(26));
	temp_mat.appendCol(K_v_mv.getCol(27));
	temp_mat = temp_mat + Nl.getTrans()*(-MU_PLUS*t*Nv_tangent)*Wl2(1,ix)*J0_val*r;
	for(int i=1; i!=K_v_mv.num_row+1; ++i){
	    K_v_mv(i,19) = temp_mat(i,1);
	    K_v_mv(i,20) = temp_mat(i,2);
	    K_v_mv(i,21) = temp_mat(i,3);
	    K_v_mv(i,22) = temp_mat(i,4);
	    K_v_mv(i,23) = temp_mat(i,5);
	    K_v_mv(i,24) = temp_mat(i,6);
	    K_v_mv(i,25) = temp_mat(i,7);
	    K_v_mv(i,26) = temp_mat(i,8);	
	    K_v_mv(i,27) = temp_mat(i,9);
	}

        Fv_m = Fv_m + Nl.getTrans()*(t*force_membrane.getTrans())*Wl2(1,ix)*J0_val*r;

                
        Flp  = Flp  + ( (n_Gpt*((force_membrane + force_repulsion).getTrans()))(1,1)*(Nl.getTrans()))
                    * Wl2(1,ix)*J0_val*r;

                
        Fv  = Fv    + ( (n_Gpt*((force_membrane + force_repulsion).getTrans()))(1,1)
                    * (Nv.getTrans())*(n_Gpt.getTrans()))*Wl2(1,ix)*J0_val*r;

	if(ix==3){
	    Fm2x =  force_membrane(1,1);
            Fm2y =  force_membrane(1,2);
            xFm2 = Gpt(1,1);
            yFm2 = Gpt(1,2);
	}

    } // end of gauss point for loop

} // end of calcSplitKsFs()


void ele2dcoupled::assembly(){

    for(int i=1; i!=length(indexv)+1; ++i){
	Fe(indexv(1,i), 1) = Fe(indexv(1,i), 1) + Fv(i,1);	// Fe, Fv are col vectors
    }

    for(int i=1; i!=length(indexp)+1; ++i){
	Fe(indexp(1,i), 1) = Fe(indexp(1,i), 1) + Fp(i,1);	// index are row vectors
    }

    for(int i=1; i!=length(indexlp)+1; ++i){
	Fe(indexlp(1,i), 1) = Fe(indexlp(1,i), 1) + Flp(i,1);	// index are row vectors
    }

    for(int i=1; i!=length(indexv_m)+1; ++i){
	Fe(indexv_m(1,i), 1) = Fe(indexv_m(1,i), 1) + Fv_m(i,1);	// index are row vectors
    }

    // not implement Fe(isnan(Fe))

    Felem = Felem + Fe.getTrans();
    connect = sctrBTOT;

    for(int i=1; i!=length(indexv)+1; ++i){
	for(int j=1; j!=length(indexv)+1; ++j){
	    Ke(indexv(1,i), indexv(1,j)) = K_vv(i,j);
	}
    }

    for(int i=1; i!=length(indexv)+1; ++i){
	for(int j=1; j!=length(indexp)+1; ++j){
	    Ke(indexv(1,i), indexp(1,j)) = K_vp(i,j);
	}
    }

    for(int i=1; i!=length(indexp)+1; ++i){
	for(int j=1; j!=length(indexv)+1; ++j){
	    Ke(indexp(1,i), indexv(1,j)) = K_pv(i,j);
	}
    }

    for(int i=1; i!=length(indexlp)+1; ++i){
	for(int j=1; j!=length(indexp)+1; ++j){
	    Ke(indexlp(1,i), indexp(1,j)) = K_lpp(i,j);
	}
    }

    for(int i=1; i!=length(indexp)+1; ++i){
	for(int j=1; j!=length(indexlp)+1; ++j){
	    Ke(indexp(1,i), indexlp(1,j)) = K_plp(i,j);
	}
    }

    for(int i=1; i!=length(indexl_plus)+1; ++i){
	for(int j=1; j!=length(indexv)+1; ++j){
	    Ke(indexl_plus(1,i), indexv(1,j)) = K_l_plusv(i,j);
	}
    }

    for(int i=1; i!=length(indexv)+1; ++i){
	for(int j=1; j!=length(indexl_plus)+1; ++j){
	    Ke(indexv(1,i), indexl_plus(1,j)) = K_vl_plus(i,j);
	}
    }

    for(int i=1; i!=length(indexl_plus)+1; ++i){
	for(int j=1; j!=length(indexv_m)+1; ++j){
	    Ke(indexl_plus(1,i), indexv_m(1,j)) = K_l_plusv_m(i,j);
	}
    }

    for(int i=1; i!=length(indexv_m)+1; ++i){
	for(int j=1; j!=length(indexl_plus)+1; ++j){
	    Ke(indexv_m(1,i), indexl_plus(1,j)) = K_v_ml_plus(i,j);
	}
    }

    for(int i=1; i!=length(indexl_minus)+1; ++i){
	for(int j=1; j!=length(indexv)+1; ++j){
	    Ke(indexl_minus(1,i), indexv(1,j)) = K_l_minusv(i,j);
	}
    }

    for(int i=1; i!=length(indexv)+1; ++i){
	for(int j=1; j!=length(indexl_minus)+1; ++j){
	    Ke(indexv(1,i), indexl_minus(1,j)) = K_vl_minus(i,j);
	}
    }

    for(int i=1; i!=length(indexl_minus)+1; ++i){
	for(int j=1; j!=length(indexv_m)+1; ++j){
	    Ke(indexl_minus(1,i), indexv_m(1,j)) = K_l_minusv_m(i,j);
	}
    }

    for(int i=1; i!=length(indexv_m)+1; ++i){
	for(int j=1; j!=length(indexl_minus)+1; ++j){
	    Ke(indexv_m(1,i), indexl_minus(1,j)) = K_v_ml_minus(i,j);
	}
    }

    for(int i=1; i!=length(indexv_m)+1; ++i){
	for(int j=1; j!=length(indexv_m)+1; ++j){
	    Ke(indexv_m(1,i), indexv_m(1,j)) = K_v_mv_m(i,j);
	}
    }

    for(int i=1; i!=length(indexv_m)+1; ++i){
	for(int j=1; j!=length(indexv)+1; ++j){
	    Ke(indexv_m(1,i), indexv(1,j)) = K_v_mv(i,j);
	}
    }

    Ii = zeros(1,43*43);
    Jj = zeros(1,43*43);
    Xx = zeros(1,43*43);

    int nt = 0;
    for(int krow=1; krow!=element_unknowns+1; ++krow){
	for(int kcol=1; kcol!=element_unknowns+1; ++kcol){
	    nt = nt+1;
            Ii(1,nt) = sctrBTOT(1,krow);
            Jj(1,nt) = sctrBTOT(1,kcol);
            Xx(1,nt) = Ke(krow,kcol);
	}
    }

} // end of assembly


void ele2dcoupled::setEnrich_node42(bool is){

    for(int i=0; i<4; ++i){
    	ele4_node[i]->enrich_node42 = is;
    }

} // setEnrich_node42()


void ele2dcoupled::setEnrich_node92(bool is){

    for(int i=0; i<9; ++i){
    	ele9_node[i]->enrich_node92 = is;
    }

} // setEnrich_node42()


bool ele2dcoupled::isNode9In(const node* in) const{

    for(int i=0; i<9; ++i){
	if(in == ele9_node[i])
	    return true;
    }

    return false;

} // isNode9In()


void ele2dcoupled::calcPhi92(){

    if(ele9_node[0]->Phi92*ele9_node[1]->Phi92 > 0)
	ele9_node[4]->Phi92 = (ele9_node[0]->Phi92 + ele9_node[1]->Phi92)*0.5;

    if(ele9_node[1]->Phi92*ele9_node[2]->Phi92 > 0)
	ele9_node[5]->Phi92 = (ele9_node[1]->Phi92 + ele9_node[2]->Phi92)*0.5;

    if(ele9_node[2]->Phi92*ele9_node[3]->Phi92 > 0)
	ele9_node[6]->Phi92 = (ele9_node[2]->Phi92 + ele9_node[3]->Phi92)*0.5;

    if(ele9_node[3]->Phi92*ele9_node[0]->Phi92 > 0)
	ele9_node[7]->Phi92 = (ele9_node[3]->Phi92 + ele9_node[0]->Phi92)*0.5;

} // calcPhi92()


void ele2dcoupled::setPhi42EqualPhi92(){

    for(int i=0; i<4; ++i){
    	ele4_node[i]->Phi42 = ele9_node[i]->Phi92;
    }

} // setPhi42EqualPhi92()


void ele2dcoupled::initTEleTotalDtdx(){

    t_element_total = zeros(1,18);
    dtdx_total = zeros(1,36);
    t_element=zeros(9,2);

} // initTEleTotalDtdx()


void ele2dcoupled::calcTEleTotalDtdx(matrix &t_element_pt, matrix & dtdx_pt, const int kn){

    REAL norm_t_ele = norm(t_element_pt);	// t_element_pt is 1x2
    t_element(kn,1) = t_element_pt(1,1)/norm_t_ele;
    t_element(kn,2) = t_element_pt(1,2)/norm_t_ele;
        
    t_element_total(1,2*kn-1) = t_element(kn,1);
    t_element_total(1,2*kn) = t_element(kn,2);

    dtdx_total(1,4*kn-3) = dtdx_pt(1,1);
    dtdx_total(1,4*kn-2) = dtdx_pt(2,2);
    dtdx_total(1,4*kn-1) = dtdx_pt(1,2);
    dtdx_total(1,4*kn) = dtdx_pt(2,1);          

} // calcTEleTotalDtdx()


void ele2dcoupled::calcT_elementDtdx(){

    if(Split_ordered>0){
	t_element = zeros(9,2);
	dtdx = zeros(9,4);
	for(int i=1; i!=9+1; ++i){
	    t_element(i,1) = t_element_total(1, 1+(i-1)*2);
	    t_element(i,2) = t_element_total(1, 2+(i-1)*2);
	    dtdx(i,1) = dtdx_total(1, 1+(i-1)*4);
	    dtdx(i,2) = dtdx_total(1, 2+(i-1)*4);
	    dtdx(i,3) = dtdx_total(1, 3+(i-1)*4);
	    dtdx(i,4) = dtdx_total(1, 4+(i-1)*4);
	}
    }
    else{
	    t_element = zeros(9,2);
	    dtdx = zeros(9,4);
    }

} // calcT_elementDtdx()


void ele2dcoupled::intersectls9(){

    node* corner[5];
    corner[0] = ele9_node[0]; corner[1] = ele9_node[1];
    corner[2] = ele9_node[2]; corner[3] = ele9_node[3];
    corner[4] = ele9_node[0];

    int n = 0;

    //loop on element 9 edges
    for(int i=1; i<5; ++i){
	node* n1 = corner[i-1];
	node* n2 = corner[i];
if(n>=2) break;
	if( (n1->Phi92)*(n2->Phi92)<0 ){
	    if(i==1){
	 	n = n+1;
	
	 	xc2(1,n) = (n1->nodex)-(n1->Phi92)*((n2->nodex)-(n1->nodex))/((n2->Phi92)-(n1->Phi92));
		yc2(1,n) = n1->nodey;

		xcp2(1,n) = 2*(xc2(1,n)-((n1->nodex)+(n2->nodex))*.05)/fabs((n1->nodex)-(n2->nodex));
		ycp2(1,n) = -1;
	    }
	    else if(i==2){
		n = n+1;

		yc2(1,n) = (n1->nodey)-(n1->Phi92)*((n2->nodey)-(n1->nodey))/((n2->Phi92)-(n1->Phi92));
		xc2(1,n) = n1->nodex;

		xcp2(1,n) = 1;
		ycp2(1,n) = 2*(yc2(1,n)-((n1->nodey)+(n2->nodey))*0.5)/fabs((n1->nodey)-(n2->nodey));
	    }
	    else if(i==3){
		n = n+1;
	
		xc2(1,n) = (n1->nodex)-(n1->Phi92)*((n2->nodex)-(n1->nodex))/((n2->Phi92)-(n1->Phi92));
		yc2(1,n) = n1->nodey;

		xcp2(1,n) = 2*(xc2(1,n)-((n1->nodex)+(n2->nodex))*0.5)/fabs((n1->nodex)-(n2->nodex));
		ycp2(1,n) = 1;
	    }
	    else if(i==4){
		n = n+1;
		
		yc2(1,n) = (n1->nodey)-(n1->Phi92)*((n2->nodey)-(n1->nodey))/((n2->Phi92)-(n1->Phi92));
		xc2(1,n) = n1->nodex;

		xcp2(1,n) = -1;
		ycp2(1,n) = 2*(yc2(1,n)-((n1->nodey)+(n2->nodey))*0.5)/fabs((n1->nodey)-(n2->nodey));

	    }

	}

    } // end of for


} // intersectls9()


void ele2dcoupled::discontQ4quad(){
// order is always 3
/*
    matrix node(6,2);
    node(1,1) = -1; node(1,2) = -1; node(2,1) = 1;  node(2,2) = -1;
    node(3,1) = 1;  node(3,2) = 1;  node(4,1) = -1; node(4,2) = 1;
    node(5,1) = xcp2(1,1); node(5,2) = ycp2(1,1);
    node(6,1) = xcp2(1,2); node(6,2) = ycp2(1,2);

    // get decomposed triagnles
    matrix nodex; nodex.appendCol(node.getCol(1));
    matrix nodey; nodey.appendCol(node.getCol(2));
    matrix tri = delaunay(nodex, nodey);
*/
    std::vector<REAL> nodex;
    nodex.push_back(-1); nodex.push_back(1);
    nodex.push_back(1);  nodex.push_back(-1);
    nodex.push_back(xcp2(1,1)); nodex.push_back(xcp2(1,2));

    std::vector<REAL> nodey;
    nodey.push_back(-1); nodey.push_back(-1);
    nodey.push_back(1);  nodey.push_back(1);
    nodey.push_back(ycp2(1,1)); nodey.push_back(ycp2(1,2));

    matrix tri = delaunay(nodex.begin(), nodex.end(), nodey.begin(), nodey.end());

    tri = tricheck(nodex.begin(), nodey.begin(), tri);	// tri is nx3


    // loop over subtriangles to get quadrature points and weights
    int pt = 1;
    for(int e=1; e<tri.num_row+1; ++e){
	// [w,q]=quadrature(order,'TRIANGULAR',2)
	matrix q(7,2);
	matrix w(7,1);

        q(1,1) = 0.1012865073235; q(1,2) = 0.1012865073235;
        q(2,1) = 0.7974269853531; q(2,2) = 0.1012865073235;
        q(3,1) = 0.1012865073235; q(3,2) = 0.7974269853531; 
        q(4,1) = 0.4701420641051; q(4,2) = 0.0597158717898;
        q(5,1) = 0.4701420641051; q(5,2) = 0.4701420641051;
        q(6,1) = 0.0597158717898; q(6,2) = 0.4701420641051; 
        q(7,1) = 0.3333333333333; q(7,2) = 0.3333333333333;
        
        w(1,1) = 0.1259391805448; 
        w(2,1) = 0.1259391805448; 
        w(3,1) = 0.1259391805448; 
        w(4,1) = 0.1323941527885;
        w(5,1) = 0.1323941527885;
        w(6,1) = 0.1323941527885;
        w(7,1) = 0.2250000000000; 

	// transform quadrature points into the parent element
	matrix coord(3,2);	// coord is 3x2
	coord(1,1) = nodex[tri(e,1)-1]; coord(1,2) = nodey[tri(e,1)-1];
	coord(2,1) = nodex[tri(e,2)-1]; coord(2,2) = nodey[tri(e,2)-1];
	coord(3,1) = nodex[tri(e,3)-1]; coord(3,2) = nodey[tri(e,3)-1];

	std::vector<REAL> temp_col;
	temp_col.push_back(1); temp_col.push_back(1); temp_col.push_back(1);
	matrix coord_3x3;
	coord_3x3.appendCol(coord.getCol(1)); coord_3x3.appendCol(coord.getCol(2));
	coord_3x3.appendCol(temp_col);
	REAL a = det(coord_3x3)*0.5;

	if(a<0){	// need to swap connectivity
	    matrix coord_temp = coord;
	    coord.clear();
	    coord.appendRow(coord_temp.getRow(2));
	    coord.appendRow(coord_temp.getRow(1));
	    coord.appendRow(coord_temp.getRow(3));
	    
	    coord_3x3.clear();
	    coord_3x3.appendCol(coord.getCol(1)); coord_3x3.appendCol(coord.getCol(2));
	    coord_3x3.appendCol(temp_col);
	    a = det(coord_3x3)*0.5;
	}

	if(a!=0){
	    for(int n=1; n<8; ++n){
		matrix N(3,1);
		REAL xi = q(n,1);
		REAL eta = q(n,2);
		N(1,1) = 1-xi-eta; N(2,1) = xi; N(3,1) = eta;
		matrix Q = N.getTrans()*coord;	// Q is 1x2
		Qx(1,pt) = Q(1,1);
		Qy(1,pt) = Q(1,2);
		W(1,pt)  = 2*w(n,1)*a;
		pt = pt+1;
	    }
	}
        if(pt>28) break;

    } // end for


} // discontQ4quad()


void ele2dcoupled::quadrature(){
// order is always be 3, always "GAUSS", always 2

    matrix r1pt(1,3);
    r1pt(1,1) = 0.774596669241483;
    r1pt(1,2) =-0.774596669241483;
    r1pt(1,3) = 0.000000000000000;

    matrix r1wt(1,3);
    r1wt(1,1) = 0.555555555555556;
    r1wt(1,2) = 0.555555555555556; 
    r1wt(1,3) = 0.888888888888889;   

    matrix quadpoint(9,2);
    matrix quadweight(9,1);
    int n=1;
    for(int i=1; i<4; ++i){
	for(int j=1; j<4; ++j){
	    Qx(1,n) = r1pt(1,i);
	    Qy(1,n) = r1pt(1,j);
	    W(1,n)  = r1wt(1,i)*r1wt(1,j);
	   
	    n = n+1;
	}
    }

    ngp = 9;

} // quadrature()


void ele2dcoupled::assembly_fluid_XFEM_ndiscont_parallel(const int numnode9, const int numsnode92, 
					      const int numnode4, const int numsnode42, int num_Split_ordered, int close_DBC2){
// only change sctrBv sctrBp


    int nn9 = 9;
    int nn4 = 4;

    for(int k=1; k<nn9+1; ++k){
	sctrBv(1, 2*k-1) = 2*(ele9_node[k-1]->num_node9)-1;
	sctrBv(1, 2*k)   = 2*(ele9_node[k-1]->num_node9);
    }

    for(int k=1; k<nn4+1; ++k){
	sctrBp(1,k) = (ele4_node[k-1]->num_node4)+numnode9*2;
    }


    if( ele4_node[0]->enrich_node42==false && ele4_node[1]->enrich_node42==false
     && ele4_node[2]->enrich_node42==false && ele4_node[3]->enrich_node42==false ){

//	sctrBl_plus.clear();
//	sctrBl_minus.clear();
//	sctrBv_m.clear();
//	sctrBlp.clear();
//	matrx sctrBv_n;
	
	/* // calculate sctrB_lag
	if(border_element){
	    for(int k=1; k<4; ++k){
		int sct[3] = {2, 6, 3};
		// posLag = find(border_nodes==unique(element9(e,sct(k))));
	    }

	}
	*/
    }
    else{
	matrix sctrBXFEMv, sctrBXFEMp;	// row vector
	
	//int snt = 0;
	for(int k=0; k<nn9; ++k){
	    if(ele9_node[k]->enrich_node92){
		//snt = snt+1;
		std::vector<REAL> temp_vec;
		temp_vec.push_back(2*numnode9 + numnode4 + (ele9_node[k]->pos92));
		sctrBXFEMv.appendCol(temp_vec);
	    }
	}

	for(int k=0; k<nn4; ++k){
	    if(ele4_node[k]->enrich_node42){
		std::vector<REAL> temp_vec;
		temp_vec.push_back(2*numnode9 + numnode4 + numsnode92 + (ele4_node[k]->pos42));
		sctrBXFEMp.appendCol(temp_vec);
	    }
	}

/* don't need to update sctrBlv & sctrBlp
	// dirichlet_elem2 = Split_ordered
	if(Split_ordered>0){
	    int k = Split_ordered;
	    matrix Moes_A2(1,2);
	    Moes_A2(1,1) = k; Moes_A2(1,2) = k+1;

	    if(k==num_Split_ordered && close_DBC2==1){
		Moes_A2(1,2) = 1;
	    }

//	    extra = 0;
//	    if(close_DBC2==0){
//		extra = 1;
//	    }

	    for(k=1; k<length(Moes_A2)+1; ++k){
		
           	sctrBlp(1,k) = 2*numnode9 + numsnode92 + numnode4 + numsnode42  + Moes_A2(1,k) ;
            	sctrBlv(1,k) = 2*numnode9 + numsnode92 + numnode4 
				  + numsnode42 + num_Split_ordered + Moes_A2(1,k) ;
	    }
	}
	else{
            sctrBlp.clear();
	    sctrBlv.clear();
	}
*/

	for(int i=1; i<sctrBXFEMv.num_col+1; ++i){
	    sctrBv.appendCol(sctrBXFEMv.getCol(i));
	}
	for(int i=1; i<sctrBXFEMp.num_col+1; ++i){
	    sctrBp.appendCol(sctrBXFEMp.getCol(i));
	}

    } // end else


} // assembly_fluid_XFEM_ndiscont_parallel()


void ele2dcoupled::assembly_fluid_XFEM_ndiscont_v2_parallel(const int numnode9, const int numsnode92, 
					      const int numnode4, const int numsnode42, int num_Split_ordered, int close_DBC2){

    int nn9 = 9;
    int nn4 = 4;

    for(int k=1; k<nn9+1; ++k){
	sctrBv(1, 2*k-1) = 2*(ele9_node[k-1]->num_node9)-1;
	sctrBv(1, 2*k)   = 2*(ele9_node[k-1]->num_node9);
    }

    for(int k=1; k<nn4+1; ++k){
	sctrBp(1,k) = (ele4_node[k-1]->num_node4)+numnode9*2;
    }

    int extra = 0;
    if(close_DBC2==0){
	extra = 1;
    }

    if( ele4_node[0]->enrich_node42==false && ele4_node[1]->enrich_node42==false
     && ele4_node[2]->enrich_node42==false && ele4_node[3]->enrich_node42==false ){

	sctrBl_plus.clear();
	sctrBl_minus.clear();
	sctrBv_m.clear();
	sctrBlp.clear();
//	matrx sctrBv_n;
	
	/* // calculate sctrB_lag
	if(border_element){
	    for(int k=1; k<4; ++k){
		int sct[3] = {2, 6, 3};
		// posLag = find(border_nodes==unique(element9(e,sct(k))));
	    }

	}
	*/
    }
    else{
	matrix sctrBXFEMv, sctrBXFEMp;	// row vector
	
	//int snt = 0;
	for(int k=0; k<nn9; ++k){
	    if(ele9_node[k]->enrich_node92){
		//snt = snt+1;
		std::vector<REAL> temp_vec;
		temp_vec.push_back(2*numnode9 + numnode4 + (ele9_node[k]->pos92));
		sctrBXFEMv.appendCol(temp_vec);
	    }
	}

	for(int k=0; k<nn4; ++k){
	    if(ele4_node[k]->enrich_node42){
		std::vector<REAL> temp_vec;
		temp_vec.push_back(2*numnode9 + numnode4 + numsnode92 + (ele4_node[k]->pos42));
		sctrBXFEMp.appendCol(temp_vec);
	    }
	}

	// dirichlet_elem2 = Split_ordered
	if(Split_ordered>0){
	    int k = Split_ordered;
	    matrix Moes_A2(1,2);
	    Moes_A2(1,1) = k; Moes_A2(1,2) = k+1;
	    // this is important for index vectors
	    sctrBlp = zeros(1,2); sctrBl_plus = zeros(1,2);
	    sctrBl_minus = zeros(1,2); sctrBv_m = zeros(1,2);	

	    if(k==num_Split_ordered && close_DBC2==1){
		Moes_A2(1,2) = 1;
	    }
	    extra = 0;
	    if(close_DBC2==0){
		extra = 1;
	    }

	    for(k=1; k<length(Moes_A2)+1; ++k){
		
           	sctrBlp(1,k) = 2*numnode9 + numsnode92 + numnode4 + numsnode42  + Moes_A2(1,k) ;
            	sctrBl_plus(1,k) = 2*numnode9 + numsnode92 + numnode4 
				 + numsnode42 + num_Split_ordered +extra + (Moes_A2(1,k)) ;
            	sctrBl_minus(1,k) = 2*numnode9 + numsnode92 + numnode4 
				  + numsnode42 + 2*num_Split_ordered +2*extra + Moes_A2(1,k) ;
            	sctrBv_m(1,k) = 2*numnode9 + numsnode92 + numnode4
			      + numsnode42 + 3*num_Split_ordered +3*extra + Moes_A2(1,k) ;
//           	sctrBv_n(1,k) = 2*numnode9 + numsnode92 + numnode4
//			      + numsnode42 + 4*num_Split_ordered +4*extra + Moes_A2(1,k) ;
	    }
	}
	else{
            sctrBl_plus.clear();
            sctrBl_minus.clear();
            sctrBv_m.clear();
            sctrBlp.clear();
	}

	for(int i=1; i<sctrBXFEMv.num_col+1; ++i){
	    sctrBv.appendCol(sctrBXFEMv.getCol(i));
	}
	for(int i=1; i<sctrBXFEMp.num_col+1; ++i){
	    sctrBp.appendCol(sctrBXFEMp.getCol(i));
	}

    } // end else

} // assembly_fluid_XFEM_ndiscont_v2_parallel()


void ele2dcoupled::Find_Gauss_Euler1D(){

    Wl2 = zeros(1,3); Ql_1D = zeros(1,3);

    Ql_1D(1,1) = 0.774596669241483;
    Ql_1D(1,2) =-0.774596669241483;
    Ql_1D(1,3) = 0.000000000000000;

    Wl2(1,1) = 0.555555555555556;
    Wl2(1,2) = 0.555555555555556; 
    Wl2(1,3) = 0.888888888888889; 

} // Find_Gauss_Euler1D()



void ele2dcoupled::find_proj_normal_driv_pt_parallel(const REAL Lx, const REAL Ly, matrix&p, const std::string &plain_axi,
					   	     std::list<activepoint*>::const_iterator Active_points_beg,
					   	     std::list<activepoint*>::const_iterator Active_points_end ){

// calculate Gpt_proj, t

//    std::list<activepoint*> Active_points(Active_points_beg, Active_points_end);

    matrix Coord_Active_points, Foot_points, n_Foot_points, a_coeff, a_coeff_n;

    for(std::list<activepoint*>::const_iterator it=Active_points_beg; it!=Active_points_end; ++it){

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
	    for(int i=1; i<temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = temp_Foot_points(i,1);
		Foot_points(i,2) = temp_Foot_points(i,2);
	    }
	    for(int i=1+temp_Foot_points.num_row; i<2*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row, 1))-0.5*Lx;
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
		Coord_Active_points(i,1) = -(Lx*0.5+temp_coord(i-temp_coord.num_row, 1))-0.5*Lx;
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
		Foot_points(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row, 1))-0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-temp_Foot_points.num_row,2);
	    }
	    for(int i=1+2*temp_Foot_points.num_row; i<3*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(-Lx*0.5+temp_Foot_points(i-2*temp_Foot_points.num_row, 1))+0.5*Lx;
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
		Coord_Active_points(i,1) = -(Lx*0.5+temp_coord(i-temp_coord.num_row, 1))-0.5*Lx;
		Coord_Active_points(i,2) = temp_coord(i-temp_coord.num_row,2);
	    }
	    for(int i=1+2*temp_coord.num_row; i<3*temp_coord.num_row+1; ++i){
		Coord_Active_points(i,1) = -(-Lx*0.5+temp_coord(i-2*temp_coord.num_row, 1))+0.5*Lx;
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
		n_Foot_points(i,1) = -temp_n_Foot_points(i-temp_n_Foot_points.num_row, 1);
		n_Foot_points(i,2) = temp_n_Foot_points(i-temp_n_Foot_points.num_row, 2);
	    }
	    for(int i=1+2*temp_n_Foot_points.num_row; i<3*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = -temp_n_Foot_points(i-2*temp_n_Foot_points.num_row, 1);
		n_Foot_points(i,2) = temp_n_Foot_points(i-2*temp_n_Foot_points.num_row, 2);
	    }

	    // a_coeff
	    a_coeff_sym = zeros(3*a_coeff.num_row,2);
	    for(int i=1; i<a_coeff.num_row+1; ++i){
		a_coeff_sym(i,1) = a_coeff(i,1);
		a_coeff_sym(i,2) = a_coeff(i,2);
	    }
	    for(int i=1+a_coeff.num_row; i<2*a_coeff.num_row+1; ++i){
		a_coeff_sym(i,1) = a_coeff(i-a_coeff.num_row, 1);
		a_coeff_sym(i,2) = a_coeff(i-a_coeff.num_row, 2);
	    }
	    for(int i=1+2*a_coeff.num_row; i<3*a_coeff.num_row+1; ++i){
		a_coeff_sym(i,1) = a_coeff(i-2*a_coeff.num_row, 1);
		a_coeff_sym(i,2) = a_coeff(i-2*a_coeff.num_row, 2);
	    }

	    // a_coeff_n
	    a_coeff_n_sym = zeros(3*a_coeff_n.num_row,2);
	    for(int i=1; i<a_coeff_n.num_row+1; ++i){
		a_coeff_n_sym(i,1) = a_coeff_n(i, 1);
		a_coeff_n_sym(i,2) = a_coeff_n(i, 2);
	    }
	    for(int i=1+a_coeff_n.num_row; i<2*a_coeff_n.num_row+1; ++i){
		a_coeff_n_sym(i,1) = a_coeff_n(i-a_coeff_n.num_row, 1);
		a_coeff_n_sym(i,2) = a_coeff_n(i-a_coeff_n.num_row, 2);
	    }
	    for(int i=1+2*a_coeff_n.num_row; i<3*a_coeff_n.num_row+1; ++i){
		a_coeff_n_sym(i,1) = a_coeff_n(i-2*a_coeff_n.num_row,1);
		a_coeff_n_sym(i,2) = a_coeff_n(i-2*a_coeff_n.num_row, 2);
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
			     + 2.0*a_coeff_sym(in2,3)*a_coeff_sym(in2,1)-2.0*a_coeff_sym(in2,3)*p_local(2,1);
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
    Gpt_proj = (temp_n0%temp_roots).getTrans() + y0;	// % is left division "\" for matrix


    t(1,1) = (a_coeff_n_sym(in2,2)+a_coeff_n_sym(in2,4)*roots(1,j)+a_coeff_n_sym(in2,6)*roots(1,j)*roots(1,j));
    t(1,2) = -(a_coeff_n_sym(in2,1)+a_coeff_n_sym(in2,3)*roots(1,j)+a_coeff_n_sym(in2,5)*roots(1,j)*roots(1,j));

/*
    matrix temp_coeff(1,2);
    matrix temp_root(2,1);
    temp_coeff(1,1) = a_coeff_n_sym(in2,4); temp_coeff(1,2) = 2*a_coeff_n_sym(in2,6);
    temp_root(1,1) = 1; temp_root(2,1) = roots(1,j);

    REAL dt1dx1 = temp_coeff*temp_root*t(1,1);
    REAL dt1dx2 = temp_coeff*temp_root*t(1,2);

    REAL dt2dx1 = -temp_coeff*temp_root*t(1,1);

    temp_coeff = zeros(1,2);
    temp_coeff(1,1) = a_coeff_n_sym(in2,3); temp_coeff(1,2) = 2*a_coeff_n_sym(in2,5);
    REAL dt2dx2 = -temp_coeff*temp_root*t(1,2);

    dtdx = zeros(2,2);
    dtdx(1,1) = dt1dx1; dtdx(1,2) = dt1dx2;
    dtdx(2,1) = dt2dx1; dtdx(2,2) = dt2dx2;
*/

} // find_proj_normal_driv_pt_parallel()


void ele2dcoupled::find_membrane_force_pt_parallel23( const REAL Lx, const REAL Ly, const std::string &plain_axi,
					   	     std::list<activepoint*>::const_iterator Active_points_beg,
					   	     std::list<activepoint*>::const_iterator Active_points_end ){


//    std::list<activepoint*> Active_points(Active_points_beg, Active_points_end);

    matrix Foot_points, n_Foot_points, a_coeff3, a_coeff_n, a_coeff_T, a_coeff_curvature, a_coeff_metric_cov;


    for(std::list<activepoint*>::const_iterator it=Active_points_beg; it!=Active_points_end; ++it){

	Foot_points.appendRow( ((*it)->getFoot_points()).getRow(1) );
	n_Foot_points.appendRow( ((*it)->getN_Foot_points()).getRow(1) );
	a_coeff3.appendRow( ((*it)->getA_coeff3()).getRow(1) );
 	a_coeff_n.appendRow( ((*it)->getA_coeff_n()).getRow(1) );
	a_coeff_T.appendRow( ((*it)->getA_coeff_T()).getRow(1) );
	a_coeff_curvature.appendRow( ((*it)->getA_coeff_curvature()).getRow(1) );
	a_coeff_metric_cov.appendRow( ((*it)->getA_coeff_metric_cov()).getRow(1) );
    }

    matrix a_coeff3_sym, a_coeff_n_sym, a_coeff_T_sym, a_coeff_curvature_sym, a_coeff_metric_cov_sym;
    if(plain_axi=="axi1"){
	    matrix temp_Foot_points = Foot_points;
	    Foot_points = zeros(2*temp_Foot_points.num_row,2);
	    for(int i=1; i<temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = temp_Foot_points(i,1);
		Foot_points(i,2) = temp_Foot_points(i,2);
	    }
	    for(int i=1+temp_Foot_points.num_row; i<2*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row, 1))-0.5*Lx;
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

	    // a_coeff_n
	    a_coeff_n_sym = a_coeff_n;
	    for(int i=1; i<a_coeff_n.num_row+1; ++i){
		a_coeff_n_sym.appendRow(a_coeff_n.getRow(i));
	    }

	    // a_coeff_T
	    a_coeff_T_sym = a_coeff_T;
	    for(int i=1; i<a_coeff_T.num_row+1; ++i){
		a_coeff_T_sym.appendRow(a_coeff_T.getRow(i));
	    }

	    // a_coeff_curvature
	    a_coeff_curvature_sym = a_coeff_curvature;
	    for(int i=1; i<a_coeff_curvature.num_row+1; ++i){
		a_coeff_curvature_sym.appendRow(a_coeff_curvature.getRow(i));
	    }

	    // a_coeff_metric_cov
	    a_coeff_metric_cov_sym = a_coeff_metric_cov;
	    for(int i=1; i<a_coeff_metric_cov.num_row+1; ++i){
		a_coeff_metric_cov_sym.appendRow(a_coeff_metric_cov.getRow(i));
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
		Foot_points(i,1) = -(Lx*0.5+temp_Foot_points(i-temp_Foot_points.num_row, 1))-0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-temp_Foot_points.num_row, 2);
	    }
	    for(int i=1+2*temp_Foot_points.num_row; i<3*temp_Foot_points.num_row+1; ++i){
		Foot_points(i,1) = -(-Lx*0.5+temp_Foot_points(i-2*temp_Foot_points.num_row, 1))+0.5*Lx;
		Foot_points(i,2) = temp_Foot_points(i-2*temp_Foot_points.num_row, 2);
	    }

	    // n_Foot_points
	    matrix temp_n_Foot_points = n_Foot_points;
	    n_Foot_points = zeros(3*temp_n_Foot_points.num_row,2);
	    for(int i=1; i<temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = temp_n_Foot_points(i,1);
		n_Foot_points(i,2) = temp_n_Foot_points(i,2);
	    }
	    for(int i=1+temp_n_Foot_points.num_row; i<2*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = -temp_n_Foot_points(i-temp_n_Foot_points.num_row, 1);
		n_Foot_points(i,2) = temp_n_Foot_points(i-temp_n_Foot_points.num_row, 2);
	    }
	    for(int i=1+2*temp_n_Foot_points.num_row; i<3*temp_n_Foot_points.num_row+1; ++i){
		n_Foot_points(i,1) = -temp_n_Foot_points(i-2*temp_n_Foot_points.num_row, 1);
		n_Foot_points(i,2) = temp_n_Foot_points(i-2*temp_n_Foot_points.num_row, 2);
	    }

	    // a_coeff3
	    a_coeff3_sym = a_coeff3;
	    for(int i=1; i<a_coeff3.num_row+1; ++i){
		a_coeff3_sym.appendRow(a_coeff3.getRow(i));
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

	    // a_coeff_T
	    a_coeff_T_sym = a_coeff_T;
	    for(int i=1; i<a_coeff_T.num_row+1; ++i){
		a_coeff_T_sym.appendRow(a_coeff_T.getRow(i));
	    }
	    for(int i=1; i<a_coeff_T.num_row+1; ++i){
		a_coeff_T_sym.appendRow(a_coeff_T.getRow(i));
	    }

	    // a_coeff_curvature
	    a_coeff_curvature_sym = a_coeff_curvature;
	    for(int i=1; i<a_coeff_curvature.num_row+1; ++i){
		a_coeff_curvature_sym.appendRow(a_coeff_curvature.getRow(i));
	    }
	    for(int i=1; i<a_coeff_curvature.num_row+1; ++i){
		a_coeff_curvature_sym.appendRow(a_coeff_curvature.getRow(i));
	    }

	    // a_coeff_metric_cov
	    a_coeff_metric_cov_sym = a_coeff_metric_cov;
	    for(int i=1; i<a_coeff_metric_cov.num_row+1; ++i){
		a_coeff_metric_cov_sym.appendRow(a_coeff_metric_cov.getRow(i));
	    }
	    for(int i=1; i<a_coeff_metric_cov.num_row+1; ++i){
		a_coeff_metric_cov_sym.appendRow(a_coeff_metric_cov.getRow(i));
	    }

    } // if "axi2"

    if(plain_axi=="plain"){
            a_coeff3_sym = a_coeff3;
            a_coeff_n_sym = a_coeff_n;
            a_coeff_T_sym = a_coeff_T;
            a_coeff_curvature_sym = a_coeff_curvature;
            a_coeff_metric_cov_sym = a_coeff_metric_cov;
    } // if "plain"

//    matrix D2 = zeros(Foot_points.num_row, 1);
    REAL D2_min;
    int IDX2;
    for(int i=1; i<Foot_points.num_row+1; ++i){
	REAL D2 = pow( (Foot_points(i,1)-Gpt_proj(1,1))*(Foot_points(i,1)-Gpt_proj(1,1))
		     + (Foot_points(i,2)-Gpt_proj(1,2))*(Foot_points(i,2)-Gpt_proj(1,2)), 0.5);
 	if(i==1){
	    D2_min = D2;
	    IDX2 = i;
    	}
	else if(D2_min>D2){
	    D2_min = D2;
	    IDX2 = i;
	}
    }

//    std::vector<int> IDX2 = sort30(D2.begin(), D2.end());

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

    p_local = temp_n0*(Gpt_proj-y0).getTrans();	// p_local is 2x1

    matrix n(1,2);
    n(1,1) = a_coeff_n_sym(in2,1)+a_coeff_n_sym(in2,3)*p_local(1,1)+a_coeff_n_sym(in2,5)*p_local(1,1)*p_local(1,1);
    n(1,2) = a_coeff_n_sym(in2,2)+a_coeff_n_sym(in2,4)*p_local(1,1)+a_coeff_n_sym(in2,6)*p_local(1,1)*p_local(1,1);

    matrix t(1,2);
    t(1,1) = -n(1,2);
    t(1,2) = n(1,1);

    matrix curvature2(1,2);
    curvature2(1,1) = a_coeff_curvature_sym(in2,1)+a_coeff_curvature_sym(in2,3)*p_local(1,1)
  		    + a_coeff_curvature_sym(in2,5)*p_local(1,1)*p_local(1,1)
		    + a_coeff_curvature_sym(in2,7)*p_local(1,1)*p_local(1,1)*p_local(1,1);
    curvature2(1,2) = a_coeff_curvature_sym(in2,2)+a_coeff_curvature_sym(in2,4)*p_local(1,1)
		    + a_coeff_curvature_sym(in2,6)*p_local(1,1)*p_local(1,1)
		    + a_coeff_curvature_sym(in2,8)*p_local(1,1)*p_local(1,1)*p_local(1,1);

    matrix stress(1,2);
    stress(1,1) = a_coeff_T_sym(in2,1)+a_coeff_T_sym(in2,3)*p_local(1,1)+a_coeff_T_sym(in2,5)*p_local(1,1)*p_local(1,1);
    stress(1,2) = a_coeff_T_sym(in2,2)+a_coeff_T_sym(in2,4)*p_local(1,1)+a_coeff_T_sym(in2,6)*p_local(1,1)*p_local(1,1);

    REAL dT11ds1 = a_coeff_T_sym(in2,3) + 2*a_coeff_T_sym(in2,5)*p_local(1,1);

    REAL dcurvature1ds1 = a_coeff_curvature_sym(in2,3)+2*a_coeff_curvature_sym(in2,5)*p_local(1,1)
			+ 3*a_coeff_curvature_sym(in2,7)*p_local(1,1)*p_local(1,1);
    REAL dcurvature2ds1 = a_coeff_curvature_sym(in2,4)+2*a_coeff_curvature_sym(in2,6)*p_local(1,1)
			+ 3*a_coeff_curvature_sym(in2,8)*p_local(1,1)*p_local(1,1);

    REAL dcurvature12ds12 = 2*a_coeff_curvature_sym(in2,5)+6*a_coeff_curvature_sym(in2,7)*p_local(1,1);
    REAL dcurvature22ds12 = 2*a_coeff_curvature_sym(in2,6)+6*a_coeff_curvature_sym(in2,8)*p_local(1,1);

    matrix metric_cov(1,2);
    metric_cov(1,1) = a_coeff_metric_cov_sym(in2,1)+a_coeff_metric_cov_sym(in2,3)*p_local(1,1)
		    + a_coeff_metric_cov_sym(in2,5)*p_local(1,1)*p_local(1,1);
    metric_cov(1,2) = a_coeff_metric_cov_sym(in2,2)+a_coeff_metric_cov_sym(in2,4)*p_local(1,1)
		    + a_coeff_metric_cov_sym(in2,6)*p_local(1,1)*p_local(1,1);

    REAL dmetric_cov1ds1 = a_coeff_metric_cov_sym(in2,3) + 2*a_coeff_metric_cov_sym(in2,5)*p_local(1,1);
    REAL dmetric_cov2ds1 = a_coeff_metric_cov_sym(in2,4) + 2*a_coeff_metric_cov_sym(in2,6)*p_local(1,1);


    matrix curvature_cont(1,2);
    matrix curvature_cov(1,2);
    matrix Christoffel_1(1,2);
    REAL Christoffel_2;
    if(plain_axi=="axi1"){
   
            REAL a = pow( 1+pow(a_coeff3_sym(in2,2)+2*a_coeff3_sym(in2,3)*p_local(1,1) 
		              + 3*a_coeff3_sym(in2,4)*p_local(1,1)*p_local(1,1), 2),  0.5);
            
            REAL a_prime = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))
		         * ( a_coeff3_sym(in2,2)+2*a_coeff3_sym(in2,3)*p_local(1,1)
			   + 3*a_coeff3_sym(in2,4)*p_local(1,1)*p_local(1,1) )/a;
            
            curvature_cont(1,1) = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))/(a*a*a);
            curvature_cont(1,2) =-(1/a)*n(1,1)/(Gpt_proj(1,1)+Lx*0.5);
            
            
            curvature_cov(1,1) = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))/a;
            curvature_cov(1,2) =-(1/a)*n(1,1)*(Gpt_proj(1,1)+Lx*0.5);
            
            Christoffel_1(1,1) = a_prime/a;
            Christoffel_1(1,2) = -(1/(a*a))*n(1,2)*(Gpt_proj(1,1)+Lx*0.5);
            
            Christoffel_2 = n(1,2)/(Gpt_proj(1,1)+Lx*0.5);

    } // if "axi1"
        
    if(plain_axi=="axi2"){
        
            REAL a = pow( 1+pow(a_coeff3_sym(in2,2)+2*a_coeff3_sym(in2,3)*p_local(1,1) 
		              + 3*a_coeff3_sym(in2,4)*p_local(1,1)*p_local(1,1), 2),  0.5);
            
            REAL a_prime = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))
			 * ( a_coeff3_sym(in2,2)+2*a_coeff3_sym(in2,3)*p_local(1,1) 
			   + 3*a_coeff3_sym(in2,4)*p_local(1,1)*p_local(1,1) )/a;
            
            curvature_cont(1,1) = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))/(a*a*a);
            curvature_cont(1,2) =-(1/a)*n(1,1)/(Gpt_proj(1,1)+Lx*0.5);
            
            
            curvature_cov(1,1) = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))/a;
            curvature_cov(1,2) =-(1/a)*n(1,1)*(Gpt_proj(1,1)+Lx*0.5);
            
            Christoffel_1(1,1) = a_prime/a;
            Christoffel_1(1,2) = -(1/(a*a))*n(1,2)*(Gpt_proj(1,1)+Lx*0.5);
            
            Christoffel_2 = n(1,2)/(Gpt_proj(1,1)+Lx*0.5);

    } // if "axi2"       
	
    if(plain_axi=="plain"){
        
            REAL a = pow( 1+pow(a_coeff3_sym(in2,2)+2*a_coeff3_sym(in2,3)*p_local(1,1) 
			      + 3*a_coeff3_sym(in2,4)*p_local(1,1)*p_local(1,1), 2),  0.5);
        
            REAL a_prime = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))
			 * ( a_coeff3_sym(in2,2)+2*a_coeff3_sym(in2,3)*p_local(1,1) 
			   + 3*a_coeff3_sym(in2,4)*p_local(1,1)*p_local(1,1) )/a;
        
            curvature_cont(1,1) = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))/(a*a*a);
            
            curvature_cont(1,2) = 0;
        
            curvature_cov(1,1) = (2*a_coeff3_sym(in2,3)+6*a_coeff3_sym(in2,4)*p_local(1,1))/a;
            curvature_cov(1,2) = 0;
        
            Christoffel_1(1,1) = a_prime/a;
            Christoffel_1(1,2) = 0;
        
            Christoffel_2 = 0; 
       
    } // if "plain"

    REAL delta_H = 1/(sqrt(metric_cov(1,1)*metric_cov(1,2)))
	         * ( ((dmetric_cov1ds1*metric_cov(1,2)+dmetric_cov2ds1*metric_cov(1,2))*1/(2*sqrt(metric_cov(1,1)*metric_cov(1,2))) 
	             + sqrt(metric_cov(1,1)*metric_cov(1,2))*dmetric_cov1ds1)*(dcurvature1ds1+dcurvature2ds1) 
	             + sqrt(metric_cov(1,1)*metric_cov(1,2))*metric_cov(1,1)*(dcurvature12ds12+dcurvature22ds12) );

    REAL force_membrane_bending = -KC*( delta_H+0.5*(curvature_cont(1,1)+curvature_cont(1,2))
					   	   *(  (curvature_cont(1,1)+curvature_cont(1,2))
						      *(curvature_cont(1,1)+curvature_cont(1,2))
						     -4*(curvature_cont(1,1)*curvature_cont(1,2)) 
						     ) 
				       ) 
				  + GAMMA*(curvature_cont(1,1)+curvature_cont(1,2));

    matrix force_membrane_stretch = (curvature_cont(1,1)*stress(1,1)+curvature_cont(1,2)*stress(1,2))*n
				  - ((dT11ds1+(2*Christoffel_1(1,1))*stress(1,1)))*t;

    force_membrane  =  force_membrane_bending*n + force_membrane_stretch;

    //mean_curvature = (curvature_cont(1)+curvature_cont(2));


} // find_membrane_force_pt_parallel23()



}// end of memFluid
