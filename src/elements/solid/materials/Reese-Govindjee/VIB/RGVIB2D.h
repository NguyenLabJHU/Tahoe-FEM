/* $Id: RGVIB2D.h,v 1.2 2003-03-25 06:30:24 thao Exp $ */
/* created: TDN (01/22/2001) */

#ifndef _RG_VIB_2D_H_
#define _RG_VIB_2D_H_

/* base classes */
#include "RGBaseT.h"
#include "ViscVIB.h"

namespace Tahoe {

/* forward declarations */
class CirclePointsT;

/** 2D Isotropic ViscVIB using Ogden's spectral formulation */
class RGVIB2D: public RGBaseT, public ViscVIB
{
  public:
  
	/* constructor */
	RGVIB2D(ifstreamT& in, const FSMatSupportT& support);

	/* destructor */
	~RGVIB2D(void);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/* class specific initializations */ 
        virtual void Initialize(void); 

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
	void ComputeEigs_e(const dArrayT& eigenstretch,dArrayT& eigenstretch_e, 
			   dArrayT& eigenstress,dSymMatrixT& eigenmodulus);
  
	/* stresses and moduli*/
  	void dWdE(const dArrayT& eigenstretch, dArrayT& eigenstress,int etype);

  	void ddWddE(const dArrayT& eigenstretch, 
		    dArrayT& eigenstress, dSymMatrixT& eigenmodulus,int etype);

  	void Calgorithm(const dArrayT& eigenstretch, const dArrayT& eigenstretch_e,
			dArrayT& eigenstress, dSymMatrixT& eigenmodulus,dMatrixT& Calg);

	/* return true of model is purely 2D, plain stress */
	virtual bool PurePlaneStress(void) const { return true; };

  private:

	void ComputeiKAB(const double& J, const double& Je, const dArrayT& eigenstress, 
			 const dSymMatrixT& eigenmodulus);

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
        dArrayT     ftau_E; 
        dSymMatrixT fDtauDep_E; 
        dArrayT     ftau_I; 
        dSymMatrixT fDtauDep_I; 
        dMatrixT    fModMat; 
         
        /* return values */ 
        dMatrixT        fModulus; 
        dSymMatrixT     fStress; 
         
  	/*inelastic moduli*/
  	dMatrixT fiKAB;

	/*inverse viscosities*/
	double fietaS;
	double fietaB;
	
	/*2D geometric constraint*/
	const double fconst;  //fconst = 0.5 for 2D formulation.
};
}
#endif /* _RG_VIB_2D_H_ */
