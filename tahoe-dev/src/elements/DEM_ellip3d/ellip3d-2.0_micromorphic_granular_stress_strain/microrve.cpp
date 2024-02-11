#include "microrve.h"
#include <cstring>


namespace dem{

microrve::microrve(matrix xc, int i, int j, int k, REAL dx_t, REAL dy_t, REAL dz_t){

    x_center = xc;	// center of the RVE

    dx = dx_t; dy = dy_t; dz = dz_t;

    // the boundary coordinates of the RVE
    x_min = xc(1,1)-dx*0.5;
    x_max = xc(1,1)+dx*0.5;
    y_min = xc(2,1)-dy*0.5;
    y_max = xc(2,1)+dy*0.5;
    z_min = xc(3,1)-dz*0.5;
    z_max = xc(3,1)+dz*0.5;

    sigma_grain = zeros(3,3);	// granular stress in this RVE
    xi_arm_nplus = zeros(3,1);	// moment arm of this RVE, pointing from the center of Omega to the center of RVE,
				// which is based on all the particles in the RVE, refer to notes pg79
    xi_arm_n = zeros(3,1);
    xi_arm_nminus = zeros(3,1);
    
    rel_posi_init = zeros(3,1);
    rel_posi_curr = zeros(3,1);

    v = dx*dy*dz;		// volume of this RVE, dxdxd
    rho = 0;		// density of this

    num_pctl_nplus = 0;

    num_x = i;
    num_y = j;
    num_z = k;

} // microrve()


void microrve::calculateA_xi(REAL dt){

    a_xi = (xi_arm_nplus-2*xi_arm_n+xi_arm_nminus)/(dt*dt);

} // calculateA_xi()


void microrve::updateRVE_z(REAL dz_RVE, REAL zOmega_min){

    dz = dz_RVE;
    x_center(3,1) = zOmega_min+0.5*dz_RVE+(num_z-1)*dz_RVE;

    z_min = x_center(3,1)-dz*0.5;
    z_max = x_center(3,1)+dz*0.5;

    v = dx*dy*dz;

} // updateRVE_z()



void microrve::updateRVE_xyz(REAL dx_RVE, REAL dy_RVE, REAL dz_RVE, 
			     REAL xOmega_min, REAL yOmega_min, REAL zOmega_min){

    dx = dx_RVE;
    dy = dy_RVE;
    dz = dz_RVE;
    x_center(1,1) = xOmega_min+0.5*dx_RVE+(num_x-1)*dx_RVE;
    x_center(2,1) = yOmega_min+0.5*dy_RVE+(num_y-1)*dy_RVE;
    x_center(3,1) = zOmega_min+0.5*dz_RVE+(num_z-1)*dz_RVE;


    x_min = x_center(1,1)-dx*0.5;
    x_max = x_center(1,1)+dx*0.5;

    y_min = x_center(2,1)-dy*0.5;
    y_max = x_center(2,1)+dy*0.5;

    z_min = x_center(3,1)-dz*0.5;
    z_max = x_center(3,1)+dz*0.5;

    v = dx*dy*dz;

} // updateRVE_z()



} // namespace dem
