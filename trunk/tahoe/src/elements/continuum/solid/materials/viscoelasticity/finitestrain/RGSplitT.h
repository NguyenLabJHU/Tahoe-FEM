/* $Id: RGSplitT.h,v 1.2 2004-07-15 08:29:30 paklein Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _RGSplitT_
#define _RGSplitT_

/* base class */
#include "RGViscoelasticityT.h"

/* direct members */
#include "SpectralDecompT.h"

namespace Tahoe {

/* forward declarations */
class PotentialT;

class RGSplitT: public RGViscoelasticityT
{
   public:
  
	/* constructor/destructor */
	RGSplitT(void);

	~RGSplitT(void);
	
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

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

   private:
	void ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
	                   dArrayT& eigenstress, dSymMatrixT& eigenmodulus);
	void ComputeiKAB(dSymMatrixT& eigenmodulus, double& bulkmodulus);
    
   protected:

	/* return values */
	dMatrixT fModulus;
	dSymMatrixT fStress;
	
	/* free energy potential */
	PotentialT* fPot_EQ;
	PotentialT* fPot_NEQ;

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
	dSymMatrixT fC_v_2D;
	
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
#endif /* _RGSplitT_ */
