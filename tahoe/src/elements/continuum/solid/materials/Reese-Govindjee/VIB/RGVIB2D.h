/* $Id: RGVIB2D.h,v 1.3 2002-10-14 16:04:02 thao Exp $ */
/* created: TDN (01/22/2001) */

#ifndef _RG_VIB_2D_H_
#define _RG_VIB_2D_H_

/* base classes */
#include "RGBaseT.h"
#include "Material2DT.h"
#include "ViscVIB.h"

namespace Tahoe {

/* forward declarations */
class CirclePointsT;

/** 2D Isotropic ViscVIB using Ogden's spectral formulation */
class RGVIB2D: public RGBaseT, public Material2DT, public ViscVIB
{
  public:
  
	/* constructor */
	RGVIB2D(ifstreamT& in, const FiniteStrainT& element);

	/* destructor */
	~RGVIB2D(void);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/*initialize history variables*/
	virtual void PointInitialize(void);
 
	/* class specific initializations */ 
        virtual void Initialize(void); 

	virtual bool Need_F_last(void) const { return true; };

        /*compute output variables*/ 
        virtual int NumOutputVariables() const; 
        virtual void OutputLabels(ArrayT<StringT>& labels) const; 
        virtual void ComputeOutput(dArrayT& output); 
 
	/* strain energy density */
	virtual double StrainEnergyDensity(void);

        /* spatial description */ 
        virtual const dMatrixT& c_ijkl(void); 
        virtual const dSymMatrixT& s_ij(void); 
 
        /* material description */ 
        virtual const dMatrixT& C_IJKL(void); // material tangent moduli 
        virtual const dSymMatrixT& S_IJ(void); // PK2 stress 

  protected:
  
        enum EnergyType {Inelastic=0, Elastic=1}; 

	/*principal elastic stretches*/
	virtual void ComputeEigs_e(const dArrayT& eigenstretch, 
				   dArrayT& eigenstretch_e,double& corr_factor,
				   const double& dt, dArrayT& eigenstress, 
				   dSymMatrixT& eigenmodulus);
  
	/* stresses and moduli*/
  	void sigA(const dArrayT& eigenstretch, dArrayT& eigenstress, 
			  int etype);

  	void dtauAdepB(const dArrayT& eigenstretch, 
				    dArrayT& eigenstress, dSymMatrixT& eigenmodulus, 
				    int etype);

  	void Calgorithm(const dArrayT& eigenstretch, 
			    dArrayT& eigenstress, dSymMatrixT& eigenmodulus, 
			    dMatrixT& calg);

	/* return true of model is purely 2D, plain stress */
	virtual bool PurePlaneStress(void) const { return true; };

  private:

	void ComputeiKAB(double& Jv, double& Je, const double& dt, 
			 dArrayT& eigenstress, dSymMatrixT& eigenmodulus);

  	/* calculates "bond" lengths from Lagrangian stretch eigenvalues */
	void ComputeLengths(const dArrayT& eigenstretch, int etype);
  
  	/* initialize angle tables */
  	void Construct(void); 

  protected:
	
  	/* integration point generator */
  	CirclePointsT*	fCircle;  
 
  private:  
        /* work space */ 
        dSymMatrixT fb; 
        dArrayT     fEigs; 
        dArrayT     fEigs_e; 
        dArrayT     fsigA_E; 
        dSymMatrixT fdtauAdepB_E; 
        dArrayT     fsigA_I; 
        dSymMatrixT fdtauAdepB_I; 
        dMatrixT    fCalg; 
        dMatrixT    fModMat; 
         
        /* return values */ 
        dMatrixT        fModulus; 
        dSymMatrixT     fStress; 
         
        /*jacobian of F*/ 
        double fJ; 

  	/*inelastic moduli*/
  	dMatrixT fiKAB;

	/*inverse viscosities*/
	double fietaS;
	double fietaB;
	
	/*2D geometric constraint*/
	double fconst;

	/*ratio of corrector to predictor*/
	dArrayT fcorr_ratio;
};
}
#endif /* _RG_VIB_2D_H_ */
