#ifndef MICRORVE_H
#define MICRORVE_H

#include "realtypes.h"
#include "matrix.h"



namespace dem {

class microrve{

friend class assembly;

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

    matrix rel_posi_init;	// relative position vector in reference configuration
    matrix rel_posi_curr;	// relative position vector in current configuration

    matrix a_xi;		// acceleration of moment arm, row vector


    REAL v;		// volume of this RVE, dxdxd
    REAL rho;		// density of this RVE

    int num_pctl_nplus;	// number of particles in this RVE	

    // the number coordinate of this RVE
    // in x,y and z directions
    int num_x;
    int num_y;
    int num_z;

    // dimensions
    REAL dx;
    REAL dy;
    REAL dz;

public:
    microrve(matrix xc, int i, int j, int k, REAL dx, REAL dy, REAL dz);	
    // xc is the coordinates of center, d is the dimension of this cuboidal RVE
    // xc should be row vector (1,3)
   
    REAL getDx(){return dx;}
    REAL getDy(){return dy;}
    REAL getDz(){return dz;}
    matrix getCurrRelativePosition(){return rel_posi_curr;}
    matrix getInitRelativePosition(){return rel_posi_init;}
    int getNumX() const {return num_x;}
    int getNumY() const {return num_y;}
    int getNumZ() const {return num_z;}

    void setSigma_grain(matrix sigma){sigma_grain = sigma;}
    void setXi_arm_nplus(matrix xi){xi_arm_nplus = xi;}

    void setRho(REAL val){rho = val;}
    void calculateA_xi(REAL dt);		// calculate the acceleration of moment arm

    void updateRVE_z(REAL , REAL);	// change min, max and v
    void updateRVE_xyz(REAL dx_RVE, REAL dy_RVE, REAL dz_RVE, REAL xOmega_min, REAL yOmega_min, REAL zOmega_min); // March 17, 2014

    void setInitRelativePosition() {rel_posi_init = rel_posi_curr;}	
    void setInitRelativePosition(matrix a) {rel_posi_init = a;}


}; // class microrve



} // namespace dem



#endif
