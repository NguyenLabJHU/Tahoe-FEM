#include "microrve.h"
#include <cstring>


namespace dem{

/*microrve::microrve(matrix xc, REAL d){

    x_center = xc;	// center of the RVE

    // the boundary coordinates of the RVE
    x_min = xc(1,1)-d*0.5;
    x_max = xc(1,1)+d*0.5;
    y_min = xc(1,2)-d*0.5;
    y_max = xc(1,2)+d*0.5;
    z_min = xc(1,3)-d*0.5;
    z_max = xc(1,3)+d*0.5;

    sigma_grain = zeros(3,3);	// granular stress in this RVE
    xi_arm_nplus = zeros(3,1);	// moment arm of this RVE, pointing from the center of Omega to the center of RVE,
				// which is based on all the particles in the RVE, refer to notes pg79
    xi_arm_n = zeros(3,1);
    xi_arm_nminus = zeros(3,1);
    
    v = d*d*d;		// volume of this RVE, dxdxd
    rho = 0;		// density of this

    num_pctl_nplus = 0;
    num_pctl_n = 0;

} // microrve()


microrve::calculateA_xi(){

    a_xi = 0.5*(xi_arm_nplus-2*xi_arm_n+xi_arm_nminus);

} // calculateA_xi()


*/

} // namespace dem
