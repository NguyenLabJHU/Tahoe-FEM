#include "activepoint.h"


namespace memFluid{

void activepoint::init(){

    act_node->act_point = true;

// initialize variables
    Foot_points = zeros(1,2);	// the coord of foot points
    n_Foot_points = zeros(1,2);

    // initialize Grid Based Particles
    J_Foot_points = 0;
    distanceMin = 0;
    E_Foot_points = zeros(1,2);
    T_Foot_points = zeros(1,2);
    F_Foot_points = zeros(1,2);
    vel_Foot_points_membrane = zeros(1,2);
    vel_membrane = zeros(1,2);
    a_coeff = zeros(1,3);
    a_coeff3 = zeros(1,4);
    a_coeff_vel = zeros(1,6);
    a_coeff_n = zeros(1,6);

    a_coeff_E = zeros(1,6);
    a_coeff_F = zeros(1,6);
    a_coeff_T = zeros(1,6);
    a_coeff_J = zeros(1,3);

    curvature_Foot_points = zeros(1,2);
    metric_cov_Foot_points = zeros(1,2);

    a_coeff_curvature = zeros(1,8);
    a_coeff_metric_cov = zeros(1,6);

    // solve for velocity and pressure
    vel_Foot_points_in = zeros(1,2);
    vel_Foot_points_out = zeros(1,2);
    vel_membrane_update_normal = zeros(1,2);
    vel_membrane_update_tangent = 0;
    vel_membrane_gradient = 0;
    x_proj = zeros(1,2);
    x_proj_m = zeros(1,2);
    v_div_pt = zeros(1,2);

    // update particles position
    Foot_points_dt2 = zeros(1,2);
    n_Foot_points_dt2 = zeros(1,2);
    n_Foot_points_dt20 = zeros(1,2);
    Foot_points_dt2_tangent = zeros(1,2);
    rotation_magnitude = zeros(1,2);
    Omega_e = 0;

    // ressampling footpoints and calculate geometry
    T_plot = zeros(1,2);

} // init()


activepoint::activepoint(node* pt){

    act_node = pt;
    init();

} // activepoint()


void activepoint::initAllFoots(){

    vel_Foot_points_membrane = zeros(1,2);
    E_Foot_points = zeros(1,2);
    T_Foot_points = zeros(1,2);
    F_Foot_points = zeros(1,2);
    vel_membrane = zeros(1,2);
    J_Foot_points = 0;

} // initAllFoots()


void activepoint::initCoeff_E_F_T_J(){

//    a_coeff = zeros(1,3);
    a_coeff_E = zeros(1,6);
    a_coeff_F = zeros(1,6);
    a_coeff_T = zeros(1,6);
    a_coeff_J = zeros(1,3);

} // initCoeff_E_F_T_J()


void activepoint::initAllDts(){

    Foot_points_dt2 = zeros(1,2);
    n_Foot_points_dt2 = zeros(1,2);
    n_Foot_points_dt20 = zeros(1,2);
    Foot_points_dt2_tangent = zeros(1,2);
    rotation_magnitude = zeros(1,2);

} // end of initAllDts()


void activepoint::initAllCoef_AllFoots(){

    a_coeff =zeros(1,3);
    a_coeff3 =zeros(1,4);
    a_coeff_n =zeros(1,6);
    a_coeff_vel =zeros(1,6);
    a_coeff_E =zeros(1,6);
    a_coeff_F =zeros(1,6);
    a_coeff_J =zeros(1,3);
    distanceMin = 0;
    Foot_points = zeros(1,2);	// for some nan results
    n_Foot_points = zeros(1,2);
    vel_Foot_points_membrane = zeros(1,2);
    E_Foot_points = zeros(1,2);
    F_Foot_points = zeros(1,2);
    J_Foot_points = 0;
    curvature_Foot_points = zeros(1,2);
    metric_cov_Foot_points = zeros(1,2);

} // initAllCoef_AllFoots()


void activepoint::initA_coeff_curvature_metric_cov(){
    
    a_coeff_curvature = zeros(1,8);
    a_coeff_metric_cov = zeros(1,6);

} // initA_coeff_curvature_metric_cov()


void activepoint::initVel_FootsVel_membranesX_projsV(){

    vel_Foot_points_in = zeros(1,2);
    vel_Foot_points_out = zeros(1,2);
    vel_Foot_points_membrane = zeros(1,2);
    vel_membrane_update_normal = zeros(1,2);
    vel_membrane_update_tangent = 0;
    vel_membrane_gradient = 0;
    x_proj = zeros(1,2);
    x_proj_m = zeros(1,2);
    v_div_pt = zeros(1,2);


} // initVel_FootsVel_membranesX_projsV()


void activepoint::calcFoot_points_dt2(const REAL &dt){

    matrix pts_vec(1,2);
    pts_vec(1,1) = n_Foot_points(1,2);
    pts_vec(1,2) = -n_Foot_points(1,1);
    Foot_points_dt2 = Foot_points + (vel_membrane_update_normal(1,1)
                                  * n_Foot_points+vel_membrane_update_tangent
                                  * pts_vec)*dt*0.5;

} // end of calcFoot_points_dt2()


void activepoint::calcFoot_points_dt2_tangent(const REAL &dt){

    Foot_points_dt2_tangent = Foot_points+ (vel_membrane_update_tangent)*dt*0.5;

} // calcFoot_points_dt2_tangent()


void activepoint::calcN_Foot_points_dt2_Rotation_magnitude(const REAL &dt){
    
    matrix n_Foot_points_dt2_temp(3,1);
    matrix temp_Omega(3,3);
    temp_Omega(1,1) = 0;
    temp_Omega(1,2) = -Omega_e;
    temp_Omega(1,3) = 0;
    temp_Omega(2,1) = Omega_e;
    temp_Omega(2,2) = 0;
    temp_Omega(2,3) = 0;
    temp_Omega(3,1) = 0;
    temp_Omega(3,2) = 0;
    temp_Omega(3,3) = 0;
    matrix n_Foot_points_temp(3,1);
    n_Foot_points_temp(1,1) = n_Foot_points(1,1);
    n_Foot_points_temp(2,1) = n_Foot_points(1,2);
    n_Foot_points_temp(3,1) = 0;

    n_Foot_points_dt2_temp = expm(0.5*dt*temp_Omega)*n_Foot_points_temp;
    n_Foot_points_dt2(1,1) = n_Foot_points_dt2_temp(1,1);
    n_Foot_points_dt2(1,2) = n_Foot_points_dt2_temp(2,1);

//    matrix t_Foot_points_dt2_temp(2,1);
//    t_Foot_points_dt2_temp(1,1) = n_Foot_points_dt2_temp(2,1);
//    t_Foot_points_dt2_temp(2,1) = -n_Foot_points_dt2_temp(1,1);

    matrix temp_mat(1,2);
    temp_mat(1,1) = n_Foot_points_dt2_temp(1,1);
    temp_mat(1,2) = n_Foot_points_dt2_temp(2,1);
    rotation_magnitude = Omega_e*temp_mat;

} // calcN_Foot_points_dt2_Rotation_magnitude()


void activepoint::calcVel_Foot_membrane_normal_tangent(const REAL &v_m, matrix &t){

    matrix t_vec(2,1);
    t_vec(1,1) = -t(1,2); t_vec(2,1) = t(1,1);

    vel_Foot_points_membrane = (vel_Foot_points_in*t_vec+vel_Foot_points_out*t_vec)*0.5*(t_vec.getTrans())+v_m*t;
    matrix mat_temp = vel_Foot_points_in*t_vec;
    vel_membrane_update_normal(1,1) = ( mat_temp(1,1)+mat_temp(1,1) )*0.5;
    vel_membrane_update_tangent = v_m;

} // calcVel_Foot_membrane_normal_tangent()


void activepoint::updateFoot_points(matrix &v_n, const REAL &v_m, const REAL &dt){

    matrix n_Foot_points_dt2_temp(3,1);
    matrix temp_Omega(3,3);
    temp_Omega(1,1) = 0;
    temp_Omega(1,2) = -Omega_e;
    temp_Omega(1,3) = 0;
    temp_Omega(2,1) = Omega_e;
    temp_Omega(2,2) = 0;
    temp_Omega(2,3) = 0;
    temp_Omega(3,1) = 0;
    temp_Omega(3,2) = 0;
    temp_Omega(3,3) = 0;
    matrix n_Foot_points_temp(3,1);
    n_Foot_points_temp(1,1) = n_Foot_points(1,1);
    n_Foot_points_temp(2,1) = n_Foot_points(1,2);
    n_Foot_points_temp(3,1) = 0;

    n_Foot_points_dt2_temp = expm(0.5*dt*temp_Omega)*n_Foot_points_temp;

    matrix n_Foot_temp_2(1,2);
    n_Foot_temp_2(1,1) = n_Foot_points_dt2_temp(1,1);
    n_Foot_temp_2(1,2) = n_Foot_points_dt2_temp(2,1);

    matrix t_Foot_points_dt2_temp(1,2);
    t_Foot_points_dt2_temp(1,1) = n_Foot_points_dt2_temp(2,1);
    t_Foot_points_dt2_temp(1,2) = -n_Foot_points_dt2_temp(1,1);

    matrix vel_Foot_points_update = (n_Foot_temp_2*v_n.getTrans())*n_Foot_temp_2+v_m*t_Foot_points_dt2_temp;
    
    Foot_points = Foot_points+ vel_Foot_points_update*dt;



} // updateFoot_points()


void activepoint::updateFEJT_Foot_points(matrix &d, matrix &d_cont, matrix &E, matrix &Fd, const REAL &Ja, const REAL &dt){
// d, d_cont, E, Fd are all 1x2

    F_Foot_points(1,1) = Fd(1,1)*(1 + d*dt)(1,1);
    F_Foot_points(1,2) = Fd(1,2)*(1 + d*dt)(1,2);

    E_Foot_points(1,1) = E(1,1) + F_Foot_points(1,1)*d(1,1)*F_Foot_points(1,1)*dt;
    E_Foot_points(1,2) = E(1,2) + F_Foot_points(1,2)*d(1,2)*F_Foot_points(1,2)*dt;

    J_Foot_points = Ja*(1 + (d_cont(1,1)+d_cont(1,2))*dt);

    T_Foot_points(1,1) = 1.0/J_Foot_points*F_Foot_points(1,1)
               * ((LAMBDA_MEMBRANE+2*MU_MEMBRANE)*E_Foot_points(1,1) + LAMBDA_MEMBRANE*E_Foot_points(1,2) )
               * F_Foot_points(1,1);
    T_Foot_points(1,2) = 1.0/J_Foot_points*F_Foot_points(1,2)
               * ((LAMBDA_MEMBRANE+2*MU_MEMBRANE)*E_Foot_points(1,2) + LAMBDA_MEMBRANE*E_Foot_points(1,1))
               * F_Foot_points(1,2);

} // updateFEJT_Foot_points()


void activepoint::updateT_plot(){

    T_plot =  T_Foot_points(1,1)*n_Foot_points;

} // updateT_plot()


matrix activepoint::getCoords() const{

    matrix coords(1,2);
    coords(1,1) = act_node->nodex;
    coords(1,2) = act_node->nodey;

    return coords;

} // getCoords()



} // end of memFluid
