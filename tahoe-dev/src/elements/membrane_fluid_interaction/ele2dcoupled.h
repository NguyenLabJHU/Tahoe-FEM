#ifndef ELE2DCOUPLED_H
#define ELE2DCOUPLED_H

#include "node.h"
#include "activepoint.h"
#include "matrix.h"
#include "globfuncs.h"
#include "realtypes.h"
#include "parameters.h"
#include <list>
#include <vector>

namespace memFluid{


class ele2dcoupled {

private:
    int num_ele;

    REAL element_length;

    node* ele9_node[9];
    node* ele4_node[4];

    node* center_node;
    bool enrich_element2;
    bool split_elem2;
    bool W_update;
    bool border_elements;

    // dirichlet_elem2 is vector<>::Split_ordered
    int Split_ordered;	// initialized when vector<ele2dcoupled*> Split_ordered is created
			// number in the Split_ordered, if not Split_ordered -1
    // get level-set/elements intersection points
    matrix xc2;		// 1x2
    matrix yc2;		// 1x2
    matrix xcp2;	// 1x2
    matrix ycp2;	// 1x2


    REAL Le2;
    REAL Le2_old;

    matrix Ordered_point;	// coords of two ordered points in the split element

    // Gauss ponits & weights
    matrix W;
    matrix Qx;
    matrix Qy;
    int ngp, nogp2;
    matrix Wl2;
    matrix Ql_1D;
    matrix Ql_1D2;
    matrix Ql_2Dx;
    matrix Ql_2Dy;
    matrix Paired_connect2;


    matrix t_element_total;
    matrix dtdx_total;
    matrix t_element;
    matrix dtdx;

    int element_unknowns;
    int ntriplets;
    matrix sctrBv;
    matrix sctrBl_plus;
    matrix sctrBl_minus;
    matrix sctrBv_m;
    matrix sctrBTOT;
    matrix sctrBp;
    matrix sctrBlp;
    int nn9;
    int nn4;

    matrix Ii;
    matrix Jj;
    matrix Xx;
    matrix Felem;
    matrix connect;

    int v_unknown;
    int p_unknown;
    int lp_unknown;
    int l_plus_unknown;
    int l_minus_unknown;
    int v_m_unknown;
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

    // loop over gauss points
    matrix N_9;
    matrix dNdxi_9;
    matrix N_4;
    matrix dNdxi_4;

    // compute the jabocian, old values of v,F and J
    matrix N9;	// nn9x1
    matrix N4;	// nn4x1
    matrix B9;	// 2x(nn9-9)
    matrix J0;	// 2x2
    matrix J09;	// 2x2
    matrix Nv;	// 2xnn9
    matrix Bv;	// 4xnn9
    matrix B_dot_v;	// 1xnn9
    matrix Np;		// 1xnn4
    matrix Nv_tangent;	// 2x9
    matrix Bv_tangent;	// 4x9

    // split element
    matrix N_1D;
    matrix Gpt_proj;
    matrix t;
    matrix force_repulsion;
    matrix n_Gpt;
    matrix force_membrane;
    matrix Nl;
    matrix normal_proj_matrix;
    REAL Fm2x;
    REAL Fm2y;
    REAL xFm2;
    REAL yFm2;

    void init();

public:
    ele2dcoupled(std::vector<node*>::const_iterator beg, std::vector<node*>::const_iterator end, int);
/*    ~ele2dcoupled(){

        // it is important to release memory pointed to by pointers in the container,
        // otherwise memory leaking occurs
	for(int i=0; i<9; ++i)
	    delete ele9_node[i];    	
	for(int i=0; i<4; ++i)
	    delete ele4_node[i]; 

    	delete center_node;
    }
*/
    int getEle_num() const {return num_ele;}
    matrix getCoord9() const;
    matrix getSctr9() const;
    matrix getSctr4() const;
    matrix getPhi91() const;
    matrix getPhi92() const;
    matrix getPhi41() const;
    matrix getPhi42() const;
    matrix getGrad_phi_x() const; // for 9-node ele
    matrix getGrad_phi_y() const; // for 9-node ele
    matrix getSctrBTOT() const {return sctrBTOT;}
    matrix getOrdered_point() const {return Ordered_point;}
    matrix getN9() const {return N9;}
    matrix getXc2() const {return xc2;}
    matrix getYc2() const {return yc2;}
    matrix getXcp2() const {return xcp2;}
    matrix getYcp2() const {return ycp2;}
    REAL getFm2x() const {return Fm2x;}
    REAL getFm2y() const {return Fm2y;}
    REAL getxFm2() const {return xFm2;}
    REAL getyFm2() const {return yFm2;}
    REAL getNode9x(int i) {return ele9_node[i]->nodex;}	// i is from 0
    REAL getNode9y(int i) {return ele9_node[i]->nodey;}
    REAL getLe2_old() const {return Le2_old;}
   
    void discontQ4quad();	// order is always 3
    void quadrature();		// order is always be 3, always "GAUSS", always 2
    void initNogp2Wl2Ql2Paired_connect2();
    void calcNogp2Wl2Ql_2D(const REAL &dx9, const REAL &dy9);
    void Find_Gauss_Euler1D();
    void setLe2_oldEqualLe2() {Le2_old = Le2;}
    void setLe2(REAL val) {Le2 = val;}
    void calcLe2();
    void setTEleTotal(matrix &t_ele) {t_element_total = t_ele;}
    void setDtdxTotal(matrix &dtdx) {dtdx_total = dtdx;}
    void setTElement(matrix &t_ele) {t_element = t_ele;}
//    void setDtdx(matrix &dtdx_mat) {dtdx = dtdx_mat;}

    void setEleUnknowns(const int unknown) {element_unknowns = unknown;}
    void setNtriplets(const int ntrip) {ntriplets = ntrip;}
    void calcSctrBTOT();
    void calcNn9Nn4() {nn9 = length(sctrBv); nn4 = length(sctrBp);}
    matrix getSctrBv() const {return sctrBv;}
    matrix getSctrBp() const {return sctrBp;}
    matrix getSctrBl_plus() const {return sctrBl_plus;}
    matrix getSctrBl_minus() const {return sctrBl_minus;}
    matrix getSctrBv_m() const {return sctrBv_m;}
    matrix getSctrBlp() const {return sctrBlp;}
    int getEleUnknowns() const {return element_unknowns;}
    void setVUnknown() {v_unknown = size(sctrBv,2);}
    void setPUnknown() {p_unknown = size(sctrBp,2);}
    void setLpUnknown() {lp_unknown = size(sctrBlp,2);}
    void setLPlusUnknown() {l_plus_unknown = size(sctrBl_plus,2);}
    void setLMinusUnknown() {l_minus_unknown = size(sctrBl_minus,2);}
    void setVMUnknown() {v_m_unknown = size(sctrBv_m,2);}
    void setAllUnknowns();
    void setAllIndex();
  
    void initIJXFelemConnect();	// initialize Ii, Jj, Xx, Felem, connect
    void initKe() {Ke = zeros(element_unknowns, element_unknowns);}
    void initFe() {Fe = zeros(43,1);}
    void initKeFe() {Ke = zeros(element_unknowns, element_unknowns); Fe = zeros(43,1);}
    void initSeVe(matrix &s, matrix &v);
    void initAllKs();
    void initAllFs();

    matrix getTEleTotal() const {return t_element_total;}
    matrix getDtdxTotal() const {return dtdx_total;}
    int getNtriplets() const {return ntriplets;}
    int getNn9() const {return nn9;}
    int getNn4() const {return nn4;}
    matrix getIi() const {return Ii;}
    matrix getJj() const {return Jj;}
    matrix getXx() const {return Xx;}
    matrix getKe() const {return Ke;}
    matrix getFe() const {return Fe;}
    matrix getFelem() const {return Felem;}
    matrix getConnect() const {return connect;}
    int getPartial_enriched() const;
    int getPartial_enriched4() const;

    int getV_unknown() const {return v_unknown;}
    int getP_unknown() const {return p_unknown;}
    int getLp_unknown() const {return lp_unknown;}
    int getL_plus_unknown() const {return l_plus_unknown;}
    int getL_minus_unknown() const {return l_minus_unknown;}
    int getV_m_unknown() const {return v_m_unknown;}
    matrix getIndexv() const {return indexv;}
    matrix getIndexp() const {return indexp;}
    matrix getIndexlp() const {return indexlp;}
    matrix getIndexl_plus() const {return indexl_plus;}
    matrix getIndexl_minus() const {return indexl_minus;}
    matrix getIndexv_m() const {return indexv_m;}
    matrix getK_vv() const {return K_vv;}
    matrix getK_vp() const {return K_vp;}
    matrix getK_pv() const {return K_pv;}
    matrix getFp() const {return Fp;}
    matrix getFv() const {return Fv;}

//    int getNgp() const {return ngp;}
    matrix getW()  const {return W;}
    matrix getQx() const {return Qx;}
    matrix getQy() const {return Qy;}
    matrix getQl_1D() const {return Ql_1D;}
    matrix getQl_2Dx() const {return Ql_2Dx;}
    matrix getQl_2Dy() const {return Ql_2Dy;}
    matrix getWl2() const {return Wl2;}
    matrix getNodeCoord9() const;
    matrix getNodeCoord4() const;
    matrix getCenter_node() const;
    
    // guass integrate to get Ks Fs
    void calculateKsFs(const REAL, const REAL, const std::string &plain_axi, 
                       std::vector<REAL>::const_iterator nodex_beg, std::vector<REAL>::const_iterator nodex_end,
		       std::vector<REAL>::const_iterator nodey_beg, std::vector<REAL>::const_iterator nodey_end);
    void calcSplitKsFs(const REAL, const REAL, const std::string &plain_axi, 
		       std::list<activepoint*>::const_iterator Active_points_beg,
		       std::list<activepoint*>::const_iterator Active_points_end );	// integrate split Ks Fs

    void calcN_9(const REAL x, const REAL y);	// calculate N_9 at (x,y)
    void calcDNdxi_9(const REAL x, const REAL y); // calculate dNdxi at (x,y)
    void calcN_4(const REAL x, const REAL y);	// calculate N_4 at (x,y)
    void calcDNdxi_4(const REAL x, const REAL y); // calculate dNdxi at (x,y)
    void calcN_1D(const REAL pt1D);
    matrix getN_9() const {return N_9;}
    matrix getDNdxi_9() const {return dNdxi_9;}
    void xfem_Nmu9_parallel(const int BCjump);
    void xfem_Nmu4_parallel(const int BCjump);
    void xfem_Bmu9_parallel(const int BCjump);	// calculate B9, J09
    void get_NvBvNp_fluid_tdisc(const REAL r, const std::string &plain_axi);
    void setJ0EqualJ09() {J0 = J09;}		// set J0 = J09

    // same as that in surfacetension.h
    void find_proj_normal_driv_pt_parallel(const REAL, const REAL, matrix &Gpt, const std::string &plain_axi,
					   std::list<activepoint*>::const_iterator Active_points_beg,
					   std::list<activepoint*>::const_iterator Active_points_end ); 
    void find_membrane_force_pt_parallel23(const REAL, const REAL, const std::string &plain_axi,
					   std::list<activepoint*>::const_iterator Active_points_beg,
					   std::list<activepoint*>::const_iterator Active_points_end );

    // assembly
    void assembly();

    // get enriched node
    void setSplit_ordered(int val) {Split_ordered = val;}
    void setSplit_elem2(bool is) {split_elem2 = is;}
    void setEnrich_element2(bool is) {enrich_element2 = is;}
    void setW_update(bool is) {W_update = is;}
    void setEnrich_node42(bool);
    void setEnrich_node92(bool);
    void initCenter_node() {center_node = NULL;}
    void setElement_length(REAL val) {element_length = val;}
    void setOrdered_point(matrix &mat) {Ordered_point = mat;}
    void calcElement_length() {element_length = ele9_node[1]->nodex-ele9_node[0]->nodex;}
    REAL getElement_length() const {return element_length;}

    // get level-set/elements
    void initXc2Yc2Xcp2Ycp2() {xc2=zeros(1,2); yc2=zeros(1,2); xcp2=zeros(1,2); ycp2=zeros(1,2);}
    void initLe2Le2_old() {Le2 = 0; Le2_old = 0;}
    int  getSplit_ordered() const {return Split_ordered;}
    bool getSplit_elem2() const {return split_elem2;}
    bool getW_update() const {return W_update;}
    void intersectls9();
    bool isNode9In(const node*) const;
    void normalR2(matrix &connect2D);
    void calcPhi92();
    void setPhi42EqualPhi92();
    void setCenter_node() {center_node = ele9_node[8];}
    void initTEleTotalDtdx();
    void calcTEleTotalDtdx(matrix &t_element_pt, matrix & dtdx_pt, const int kn);
    void calcT_elementDtdx();

    // only change sctrBv sctrBp
    void assembly_fluid_XFEM_ndiscon_parallel(const int, const int, const int, const int, int, int);
    void assembly_fluid_XFEM_ndiscont_parallel(const int, const int, const int, const int, int, int);
    void assembly_fluid_XFEM_ndiscont_v2_parallel(const int, const int, const int, const int, int, int);

};




} // end of memFluid



#endif
