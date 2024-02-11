#ifndef MICRORVE_H
#define MICRORVE_H

#include "realtypes.h"
#include "matrix.h"



namespace dem {

class microrve{

/*friend class assembly;

private:
    matrix x_center;	// center of the RVE, row vector

    // the boundary coordinates of the RVE
    REAL x_min;
    REAL x_max;
    REAL y_min;
    REAL y_max;
    REAL z_min;
    REAL z_max;

    matrix sigma_grain;	// granular stress in this RVE
    matrix xi_arm_nplus;	// moment arm of this RVE, pointing from the center of Omega to the center of RVE,
				// which is based on all the particles in the RVE, refer to notes pg79
				// row vector
    matrix xi_arm_n;
    matrix xi_arm_nminus;  

    matrix a_xi;		// acceleration of moment arm, row vector


    REAL v;		// volume of this RVE, dxdxd
    REAL rho;		// density of this RVE

    int num_pctl_nplus;	// number of particles in this RVE
    int num_pctl_n;	

public:
    microrve(matrix xc, REAL d);	// xc is the coordinates of center, d is the dimension of this cuboidal RVE
					// xc should be row vector (1,3)
   
    void setSigma_grain(matrix sigma){sigma_grain = sigma;}
    void setXi_arm(matrix xi){xi_arm = xi;}

    void setRho(REAL val){rho = val;}
    void calculateA_xi();		// calculate the acceleration of moment arm
*/
}; // class microrve



} // namespace dem



#endif
