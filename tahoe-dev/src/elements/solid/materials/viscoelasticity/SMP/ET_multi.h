/* $Id: ET_multi.h,v 1.1 2016-03-17 13:43:05 tahoe.vickynguyen Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _ET_multi_
#define _ET_multi_

/* base class */
#include "LAdMatrixT.h"
#include "RGSplitT2.h"
#include "InvLangevin.h"
#include "Gamma.h"

namespace Tahoe {

class ET_multi: public RGSplitT2
{
   public:

	/* constructor/destructor */
	ET_multi(void);

	/*Bookkeeping: initialize, update/reset internal variables and the beginning and end of each time step */
	virtual void PointInitialize(void);              
	virtual void UpdateHistory(void); // element at a time
	virtual void ResetHistory(void);  // element at a time


	/** \name Input/Output: 
		Defines the model parameters
		Reads in model parameters input file
		Overloads ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);

	/*compute output variables*/
	virtual int NumOutputVariables() const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

	/*These are the main functions*/
	/*stress response*/
	/*returns the free energy density the current deformation and internal state variables*/
	virtual double StrainEnergyDensity(void); 
	/*returns the spatial tangent modulus in reduced matrix form*/
	virtual const dMatrixT& c_ijkl(void);
	/*returns the symmetric Cauchy stress tensor*/
	virtual const dSymMatrixT& s_ij(void);
	/*returns inverse(F_T)*/
	virtual const dMatrixT& ThermalDeformation_Inverse(void);
    virtual bool Need_F_last(void) const;
	
	/*Heat generation*/
    /** incremental heat generation */
	virtual double IncrementalHeat(void);
    
	/** this model does generate heat */
	virtual bool HasIncrementalHeat(void) const { return true; };
    

    
   private:
	/*fictive temperature*/
	/*returns the fictive temperature Tf given the departure from equilibrium delta_bar_neq.*/
	/*relaxation spectrum*/
	/*returns  tau_R/tauR0, the normalized strutural relaxation time*/
	double StructuralRelaxationFunc(const double Temperature, const double Sc);
	/*returns  tau_S/tauS0, normalized stress relaxation time*/
	double StressRelaxationFunc(const double Temperature, const double smag);
	/*returns  tau_Y/tauY0, normalized  relaxation time for evolution of yield strength*/
	
	/*returns the relaxation spectrum for x=tau/taukWW*/
	double StretchedExponentialSpectrum(const double x, const double tauKWW, const double betakWW);
	
	/*calculate discrete spectrum*/
	void Compute_discrete_spectrum(const double beta, const double tauWWW, const double taumin, const double taumax, const int num,
		dArrayT& times, dArrayT& spectrum, StringT filename);
	/*Numerical integration of the evolution equations for the eigenvalues of the elastic deformation tensor be, and delta_bar_neq*/
	void Compute_Tei(const dArrayT& eigs_n, const dArrayT& eigs, const dArrayT& eigs_tr, const dArrayT& DScSe,const dArrayT& stretche,const dArrayT& stressk, const dArrayT& Tfk_n, dArrayT& Tfk);
	void Compute_le(const ArrayT<dSymMatrixT>& C_vn, ArrayT<dSymMatrixT>& C_v, const dArrayT& Tfk_n, dArrayT& Tfk, double& Heat);
	void Compute_Kneq(dMatrixT& Modulus1, dMatrixT& Modulus2);

   protected:
   /*Gamma Function*/
   Gamma fGamma;
   
	/*Reference Temperature*/
	double fT0; /*Initial temperature*/
	double fTg;			/* glass transition temperature*/
	double fC1, fC2;	/*WLF constants*/
	

    /*excess heat capacity cr-cg*/
    double fdeltac;
    double fcg;
	/*shear moduli*/
	double fmur;		/*the rubbery, high temperature shear modulus*/
	double fmug;		/*the glassy, low temperature shear modulus*/
	
    double fafrac, fb_decay,fTess, fBB,  fT2;

	/*low temp viscous flow*/
	double fQS;			/*activation volume for viscoplastic flow*/
	
	/*discrete structural relaxation spectrum*/
	int fNumR;	
	dArrayT ftimesR, fdalpha;  
	StringT fInputR;  

	/*discrete stress relaxation spectrum*/
	int fNumS;
	dArrayT ftimesS, fdmu;
	StringT fInputS;  
	
	/*accessors for internal variables*/
	dArrayT fTfk;
	dArrayT fTfk_n;
    double* fHeat;
    double* fHeat_n;
    
    /*workspaces*/	
	dArrayT	fl_tr;
	dArrayT	fle;
	dArrayT ftau;
    dArrayT fstressk;
	
	LAdMatrixT fKdel;
	dArrayT fRdel;
	
	LAdMatrixT fKAB;
	LAdMatrixT fKAB2;
	dArrayT fRes;
	
	dArrayT fGA0;
	dArrayT fGA1;
	dArrayT fGA2;
	dMatrixT fDAB;
	dMatrixT fMat;
    
    /*new*/
    dArrayT  fEigs_last, fEigs_e_last, flatenth, fGB0, fGB1, fGB2, fG10, fG11, fG12, fG20, fG21, fG22;
    dMatrixT fF3D_last, fDAB2;
    dSymMatrixT fb_last;
    LAdMatrixT fKrTe, fKTer, fInverse1, fInverse2, fInverse3, fKdel2, fK1, fK2;
};
    inline bool ET_multi::Need_F_last(void) const { return true; }
}
#endif /* _ET_multi_ */
