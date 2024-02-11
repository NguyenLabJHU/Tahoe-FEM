/* $Id: RGVIB2D.h,v 1.2 2006-11-12 18:26:55 thao Exp $ */
/* created: TDN (01/22/2001) */

#ifndef _RG_VIB_2D_H_
#define _RG_VIB_2D_H_

#ifdef VIB_MATERIAL

/* base classes */
#include "RGBaseT.h"
#include "ViscVIB.h"
#include "FSMatSupportT.h"
#include "ofstreamT.h"
#include "StringT.h"

namespace Tahoe {

/* forward declarations */
class CirclePointsT;
class FSMatSupportT;
class ofstreamT;

/** 2D Isotropic ViscVIB using Ogden's spectral formulation */
class RGVIB2D: public RGBaseT, public ViscVIB
{
  public:
  
	/* constructor */
	RGVIB2D(ifstreamT& in, const FSMatSupportT& support);

	/* destructor */
	~RGVIB2D(void);

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
	virtual void ComputeEigs_e(const dArrayT& eigenstretch, 
				   dArrayT& eigenstretch_e, 
				   dArrayT& eigenstress, 
				   dSymMatrixT& eigenmodulus);
  
	/* stresses and moduli*/
  	void dWdE(const dArrayT& eigenstretch, dArrayT& eigenstress, 
			  int etype);

  	void ddWddE(const dArrayT& eigenstretch, 
				    dArrayT& eigenstress, dSymMatrixT& eigenmodulus, 
				    int etype);

  	void Calgorithm(const dArrayT& eigenstretch, 
			const dArrayT& eigenstretch_e, dArrayT& eigenstress, 
			dSymMatrixT& eigenmodulus, dMatrixT& calg);

	/* return true of model is purely 2D, plain stress */
	virtual bool PurePlaneStress(void) const { return true; };

	virtual bool HasDissipVar(void) const{ return false;};

  private:
	void ComputeiKAB(double& Jv, double& Je, dArrayT& eigenstress, 
			 dSymMatrixT& eigenmodulus);

  	/* calculates "bond" lengths from Lagrangian stretch eigenvalues */
	void ComputeLengths(const dArrayT& eigenstretch, int etype);
  
  	/* initialize angle tables */
  	void Construct(void); 

  protected:
	
  	/* integration point generator */
  	CirclePointsT*	fCircle;  
 
  private:  
	
        SpectralDecompT fSpectralDecompSpat;
	SpectralDecompT fSpectralDecompTrial;
	SpectralDecompT fSpectralDecompRef;

        /* work space */ 
        dSymMatrixT fb; 
	dSymMatrixT fb_tr;
        dArrayT     fEigs; 
        dArrayT     fEigs_e;
	dArrayT     fEigs_v;

        dArrayT     ftau_E; 
        dArrayT     ftau_I; 

        dSymMatrixT fDtauDep_E; 
        dSymMatrixT fDtauDep_I; 

        dMatrixT    fCalg; 
        dMatrixT    fModMat; 

  	dMatrixT fiKAB;
	dMatrixT fGAB;
	dMatrixT fDAB;

        /* return values */ 
        dMatrixT        fModulus; 
        dSymMatrixT     fStress; 
         
	/*inverse viscosities*/
	double fietaS;
	double fietaB;
	
	/*2D geometric constraint*/
	double fconst;
};
}
#endif /* _RG_VIB_2D_H_ */
#endif /*VIB_MATERIAL*/