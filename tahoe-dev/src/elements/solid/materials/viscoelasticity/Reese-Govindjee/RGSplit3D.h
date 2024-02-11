/* created: TDN (01/22/2001) */

#ifndef _RGSplit_3D_
#define _RGSplit_3D_

/* base classes */
#include "RGBaseT.h"

namespace Tahoe {

/* forward declarations */
class PotentialT;

class RGSplit3D: public RGBaseT
{
   public:
  
	/* constructor/destructor */
	RGSplit3D(ifstreamT& in, const FSMatSupportT& support);

	~RGSplit3D(void);

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

    /*material inelastic stress measure*/
	virtual bool HasDissipVar(void) const {return true;};

	/* return true of model is purely 2D, plain stress */
	virtual bool PurePlaneStress(void) const { return false;};
	
   private:
	void ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
	                   dArrayT& eigenstress, dSymMatrixT& eigenmodulus);
	void ComputeiKAB(dSymMatrixT& eigenmodulus, double& bulkmodulus);
    
   protected:
	/* return values */
	dMatrixT	fModulus;
	dSymMatrixT     fStress;
	
	/*internal inelastic stress measure*/
	dSymMatrixT fMatInStress;

	/* free energy potential */
	PotentialT* fPot_EQ;
	PotentialT* fPot_NEQ;

	const double fthird;
   private:  
	/* spectral operations */
	SpectralDecompT fSpectralDecompSpat;
	SpectralDecompT fSpectralDecompRef;
	SpectralDecompT fSpectralDecompTrial;

	/*work space*/
	dSymMatrixT fb;
	dSymMatrixT fb3D;
	dSymMatrixT fbe;
	dSymMatrixT fb_tr;
	dMatrixT fF3D;

	dArrayT     fEigs;
	dArrayT     fEigs_e;
	dArrayT     fEigs_bar;
	dArrayT     fEigs_ebar;
	dArrayT     fEigs_v;

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
#endif /* _RGSplit_3D_ */
