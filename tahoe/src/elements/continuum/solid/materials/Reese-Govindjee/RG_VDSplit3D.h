/* $Id: RG_VDSplit3D.h,v 1.3 2002-11-14 17:06:09 paklein Exp $ */
/* created: TDN (01/22/2001) */

#ifndef _RG_VDSplit_3D_
#define _RG_VDSplit_3D_

/* base classes */
#include "RGBaseT.h"

namespace Tahoe {

/* forward declarations */
class IsoPotentialT;

class RG_VDSplit3D: public RGBaseT
{
   public:
  
	/* constructor */
	RG_VDSplit3D(ifstreamT& in, const FDMatSupportT& support);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/* strain energy density */
	virtual double StrainEnergyDensity(void);
	virtual const dMatrixT& c_ijkl(void);
	virtual const dSymMatrixT& s_ij(void);
	virtual const dMatrixT& C_IJKL(void);
	virtual const dSymMatrixT& S_IJ(void);

        /*compute output variables*/ 
        virtual int NumOutputVariables() const; 
	virtual void OutputLabels(ArrayT<StringT>& labels) const; 
        virtual void ComputeOutput(dArrayT& output);

   private:
	void ComputeEigs_e(const dArrayT& eigenstretch, 
			   dArrayT& eigenstretch_e, dArrayT& eigenstress, 
			    dSymMatrixT& eigenmodulus);
	void ComputeiKAB(dSymMatrixT& eigenmodulus, double& bulkmodulus);

   protected:

	enum Spring {kEquilibrium = 0, kNonEquilibrium = 1};
 
	virtual double Phi(const dArrayT& eigenstretch_bar, const double& J,
			   const int SpringType) = 0;
	virtual void devTau(const dArrayT& eigenstretch_bar, 
			    dArrayT& eigenstress, const int SpringType) = 0;
	virtual double meanTau(const double& J, const int SpringType) = 0;
	virtual void DdevTauDepsilon(const dArrayT& eigenstretch_bar, 
				     dSymMatrixT&  eigenmodulus, 
				     const int SpringType)  = 0;
	virtual double DmeanTauDepsilon(const double& J, 
					const int SpringType) = 0;

   protected:
	/* return values */
	dMatrixT	fModulus;
	dSymMatrixT     fStress;
	
	const double fthird;
   private:  
	/*work space*/
	dSymMatrixT fb;
	dArrayT     fEigs;
	dArrayT     fEigs_e;
	dArrayT     fEigs_bar;
	dArrayT     fEigs_ebar;
	dArrayT	    ftau_EQ;
	dSymMatrixT fDtauDe_EQ;
	dArrayT     ftau_NEQ;
	dSymMatrixT fDtauDe_NEQ;
	dMatrixT    fModMat;
  	dMatrixT    fiKAB;
	
	/*jacobian of F*/
	double fJ;
	double fJe;

  	/*viscosities*/
	double fietaS;
	double fietaB;
};
}
#endif /* _RG_VDSplit_3D_ */
