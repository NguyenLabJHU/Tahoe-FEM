/* created: TDN (01/22/2001) */

#ifndef _RGSplitT_
#define _RGSplitT_

/* base classes */
#include "RGViscoelasticityT.h"

namespace Tahoe {

/* forward declarations */
class PotentialT;

class RGSplitT: public RGViscoelasticityT
{
   public:
  
	/* constructor/destructor */
	RGSplitT(ifstreamT& in, const FSMatSupportT& support);

	~RGSplitT(void);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	virtual void Initialize(void);

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
	void ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
	                   dArrayT& eigenstress, dSymMatrixT& eigenmodulus);
	void ComputeiKAB(dSymMatrixT& eigenmodulus, double& bulkmodulus);
    
   protected:
	/* return values */
	dMatrixT	fModulus;
	dSymMatrixT     fStress;
	
	/* free energy potential */
	PotentialT* fPot_EQ;
	PotentialT* fPot_NEQ;

	const double fthird;
  
    private:  
	/* spectral operations */
	SpectralDecompT fSpectralDecompSpat;
	
	/*work space*/
	dSymMatrixT fStretch;
	dSymMatrixT fInverse;
	dSymMatrixT fb3D;
	dSymMatrixT fbe;
	dSymMatrixT fb_tr;
	dMatrixT fF3D;

	dArrayT     fEigs_dev;
	dArrayT     fEigs;
	dArrayT     fEigs_e;
	dArrayT     fEigs_v;
	dArrayT		fEigs_tr;

	dArrayT	    ftau_EQ;
	dArrayT     ftau_NEQ;

	dSymMatrixT fStress3D;
	dSymMatrixT fDtauDe_EQ;
	dSymMatrixT fDtauDe_NEQ;
	dMatrixT fCalg;

	dMatrixT    fModulus3D;
	dMatrixT    fModMat;
  	dMatrixT    fiKAB;
	
  	/*viscosities*/
	double fietaS;
	double fietaB;
};
}
#endif /* _RGSplitT_ */
