/* $Id: BergstromBoyce.h,v 1.1 2009-04-23 15:00:25 thao Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _BergstromBoyce_
#define _BergstromBoyce_

/* base class */
#include "RGSplitT2.h"
#include "InvLangevin.h"

namespace Tahoe {

class BergstromBoyce: public RGSplitT2
{
   public:
  
	/* constructor/destructor */
	BergstromBoyce(void);

	/*compute output variables*/ 
	virtual int NumOutputVariables() const; 
	virtual void OutputLabels(ArrayT<StringT>& labels) const; 
	virtual void ComputeOutput(dArrayT& output);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);

	protected:
	double ShearViscosity(const double smag, const double sy, const double lv);	
	double DVisc_Dlamv(const double smag, const double sy, const double lv);	

    private:
	virtual void Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m, 
						const double& dtau_m, dMatrixT& Calg, const int type);

	virtual void ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
	                   dArrayT& eigenstress, dSymMatrixT& eigenmodulus, const int process_num);
    
	enum ViscType {kBergstromBoyce = 0,
					kExpBergstromBoyce};
	
   protected:
	
	/*viscous flow*/
	int fViscType;

	double fetaS0;      /*initial viscosity*/
	double fsy0;		/*yield strength*/
	double fm;			/*strain sensitivity*/
	double fepsilon;	/*regularization constant*/
	
	/*Residual*/
	dArrayT fRes;
	dArrayT fDelta;
	dMatrixT fGAB;
	dMatrixT fDAB;
	dMatrixT fDABbar;
	dMatrixT fMat;
			
};
}
#endif /* _BergstromBoyce_ */
