#ifndef ACTIVEPOINT_H
#define ACTIVEPOINT_H

#include "realtypes.h"
#include "matrix.h"
#include "globfuncs.h"
#include "parameters.h"
#include "node.h"

namespace memFluid{


class activepoint {	// length(Foot_points) == length(Active_points) in surface_tension.m

private:
    node* act_node;

    matrix Foot_points;	// the coord of foot points
    matrix n_Foot_points;

    // initialize Grid Based Particles
    REAL J_Foot_points;
    REAL distanceMin;
    matrix E_Foot_points;
    matrix T_Foot_points;
    matrix F_Foot_points;
    matrix vel_Foot_points_membrane;
    matrix vel_membrane;
    matrix a_coeff;
    matrix a_coeff3;
    matrix a_coeff_vel;
    matrix a_coeff_n;

    matrix a_coeff_E;
    matrix a_coeff_F;
    matrix a_coeff_T;
    matrix a_coeff_J;

    matrix curvature_Foot_points;
    matrix metric_cov_Foot_points;

    matrix a_coeff_curvature;
    matrix a_coeff_metric_cov;

    // solve for velocity and pressure
    matrix vel_Foot_points_in;
    matrix vel_Foot_points_out;
    matrix vel_membrane_update_normal;
    REAL vel_membrane_update_tangent;
    REAL vel_membrane_gradient;
    matrix x_proj;
    matrix x_proj_m;
    matrix v_div_pt;

    // update particles position
    matrix Foot_points_dt2;
    matrix n_Foot_points_dt2;
    matrix n_Foot_points_dt20;
    matrix Foot_points_dt2_tangent;
    matrix rotation_magnitude;
    REAL Omega_e;

    // ressampling footpoints and calculate geometry
    matrix T_plot;	

    void init();


public:
    activepoint(node*);
/*    ~activepoint(){
        // it is important to release memory pointed to by pointers in the container,
        // otherwise memory leaking occurs
	delete act_node;
    }	
*/
    // initialize Grid Based Particles
    void initAllFoots();
    void initAllCoef_AllFoots();
    void initVel_FootsVel_membranesX_projsV();
    void initA_coeff_curvature_metric_cov();
    void initCoeff_E_F_T_J();
    void setVel_membrane(matrix &mat) {vel_membrane = mat;}
    void setVel_Foot_points_membrane(matrix &mat) {vel_Foot_points_membrane = mat;}
    void setFoot_points(matrix &mat) {Foot_points = mat;}
    void setN_Foot_points(matrix &mat) {n_Foot_points = mat;}
    void setE_Foot_points(matrix &mat) {E_Foot_points = mat;}
    void setF_Foot_points(matrix &mat) {F_Foot_points = mat;}
    void setT_Foot_points(matrix &mat) {T_Foot_points = mat;}
    void setJ_Foot_points(REAL val) {J_Foot_points = val;}
    void setCurvature_Foot_points(matrix &mat) {curvature_Foot_points = mat;}
    void setMetric_cov_Foot_points(matrix &mat) {metric_cov_Foot_points = mat;}
    void setVel_Foot_points_in(matrix &mat) {vel_Foot_points_in = mat;}
    void setVel_Foot_points_out(matrix &mat) {vel_Foot_points_out = mat;}
    void setA_coeff(matrix &mat) {a_coeff = mat;}
    void setA_coeff_vel(matrix &mat) {a_coeff_vel = mat;}
    void setA_coeff_E(matrix &mat) {a_coeff_E = mat;}
    void setA_coeff_F(matrix &mat) {a_coeff_F = mat;}
    void setA_coeff_J(matrix &mat) {a_coeff_J = mat;}
    void setA_coeff_T(matrix &mat) {a_coeff_T = mat;}
    void setA_coeff_n(matrix &mat) {a_coeff_n = mat;}
    void setA_coeff3(matrix &mat)  {a_coeff3 = mat;}
    void setA_coeff_curvature(matrix &mat) {a_coeff_curvature = mat;}
    void setA_coeff_metric_cov(matrix &mat) {a_coeff_metric_cov = mat;}
    void setOmega_e(const REAL &val) {Omega_e = val;}
    void setDistanceMin(REAL val) {distanceMin = val;}
    void setPhi92EqualDistanceMin() {act_node->Phi92 = distanceMin;}

    // update particles position
    void initAllDts();	// initialize variables in "update particles position" part
    void calcFoot_points_dt2(const REAL&);
    void calcFoot_points_dt2_tangent(const REAL&);
    void calcVel_Foot_membrane_normal_tangent(const REAL &v_m, matrix & t);

    void calcN_Foot_points_dt2_Rotation_magnitude(const REAL &dt);
    void updateFoot_points(matrix &v_n, const REAL &v_m, const REAL &dt);
    void setX_projEqualFoot_points() {x_proj = Foot_points;}
    void setX_proj_m(matrix& mat) {x_proj_m = mat;}


    // ressampling and calculate geometry
    void initT_Foot_points() {T_Foot_points = zeros(1,2);}
    matrix getCoords() const;
    matrix getFoot_points() const {return Foot_points;}
    matrix getN_Foot_points() const {return n_Foot_points;}
    matrix getE_Foot_points() const {return E_Foot_points;}
    matrix getF_Foot_points() const {return F_Foot_points;}
    matrix getT_Foot_points() const {return T_Foot_points;}
    REAL   getJ_Foot_points() const {return J_Foot_points;}
    matrix getFoot_points_dt2() const {return Foot_points_dt2;}
    matrix getFoot_points_dt2_tangent() const {return Foot_points_dt2_tangent;}
    matrix getVel_membrane_update_normal() const {return vel_membrane_update_normal;}
    REAL   getVel_membrane_update_tangent() const {return vel_membrane_update_tangent;}
    matrix getVel_Foot_points_membrane() const {return vel_Foot_points_membrane;}
    matrix getCurvature_Foot_points() const {return curvature_Foot_points;}
    matrix getMetric_cov_Foot_points() const {return metric_cov_Foot_points;}
    matrix getA_coeff() const {return a_coeff;}
    matrix getA_coeff3() const {return a_coeff3;}
    matrix getA_coeff_E() const {return a_coeff_E;}
    matrix getA_coeff_F() const {return a_coeff_F;}
    matrix getA_coeff_J() const {return a_coeff_J;}
    matrix getA_coeff_n() const {return a_coeff_n;}
    matrix getA_coeff_T() const {return a_coeff_T;}
    matrix getA_coeff_vel() const {return a_coeff_vel;}
    matrix getA_coeff_curvature() const {return a_coeff_curvature;}
    matrix getA_coeff_metric_cov() const {return a_coeff_metric_cov;}
    matrix getT_plot() const {return T_plot;}
    matrix getX_proj_m() const {return x_proj_m;}
    REAL   getDistanceMin() const {return distanceMin;}
    REAL   getPhi92() const {return act_node->Phi92;}
    REAL   getGrad_phi_x() const {return act_node->Grad_phi_x;}
    REAL   getGrad_phi_y() const {return act_node->Grad_phi_y;}
    void updateFEJT_Foot_points(matrix &d, matrix &d_cont, matrix &E, matrix &Fd, const REAL &Ja, const REAL &dt);
    void updateT_plot();


}; // end of activepoints


} // end of memFluid

#endif 
