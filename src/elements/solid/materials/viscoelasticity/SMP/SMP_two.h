/* $Id: SMP_two.h,v 1.1 2009-04-23 02:51:50 thao Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _SMP_two_
#define _SMP_two_

/* base class */
#include "SMP_simple.h"
#include "InvLangevin.h"

namespace Tahoe {

class SMP_two: public SMP_simple
{
   public:

	enum EnergyType {kMooneyRivlin=0,
					kLangevin=1};
					 
	enum ViscType {kNone = -1, 
		       kSimple=0, 
					kPower=1,
					kBergStromBoyce=2}; 
  
	/* constructor/destructor */
	SMP_two(void);

	/* initialize, update/reset internal variables */
	virtual void PointInitialize(void);              
	virtual void UpdateHistory(void); // element at a time
	virtual void ResetHistory(void);  // element at a time


	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
		
	virtual double ShearViscosityConst(const double Temperature, const double deltaneq);	
	
	protected:
	virtual void Initialize(void);

   private:
	/* set inverse of thermal transformation - return true if active */
//	virtual bool SetInverseThermalTransformation(dMatrixT& F_trans_inv);  			

	virtual void Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m, 
						const double& dtau_m, dMatrixT& Calg, const int type);

	virtual void ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
	                   dArrayT& eigenstress, dSymMatrixT& eigenmodulus, const int process_num);
    
   protected:
	
	/*low temp viscous flow*/
	double fetaS1;		/*reference shear viscosity*/		
	
};
}
#endif /* _SMP_two_ */
