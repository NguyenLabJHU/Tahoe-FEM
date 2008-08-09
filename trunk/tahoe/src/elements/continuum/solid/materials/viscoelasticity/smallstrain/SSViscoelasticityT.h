/* $Id: SSViscoelasticityT.h,v 1.4 2008-08-09 15:00:38 tdnguye Exp $ */
/* created: TDN (5/31/2001) */
#ifndef _SS_VISCO_H_
#define _SS_VISCO_H_
 
#include "SSSolidMatT.h"
#include "SS_Visc_Support.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/** small strain linear viscoelastic constitutive law */
class SSViscoelasticityT: public SSSolidMatT, public SS_Visc_Support
{
	public:

	/** constructor */
	SSViscoelasticityT(void);
			
	/* apply pre-conditions at the current time step */
	virtual void InitStep(void){SSSolidMatT::InitStep();}

	/** return true if the material has history variables */
	virtual bool HasHistory(void) const { return true; };

	/* initialize history variable */
	virtual bool NeedsPointInitialization(void) const {return true;}; // declare true
	virtual void PointInitialize(void);                               // assigns storage space
	
	/* update/reset internal variables */
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

	virtual double StrainEnergyDensity(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. \note NOT IMPLEMENTED */
	virtual double Pressure(void) const {return(1.0/3.0*fStress.Trace());	};
 
	/* spatial description */ 
	virtual const dMatrixT& c_ijkl(void); // spatial tangent moduli 
	virtual const dSymMatrixT& s_ij(void); // Cauchy stress 
 
	/* material description */ 
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli 
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress 
	
	/*compute output variables*/
	virtual int NumOutputVariables() const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

	/*@}*/

	protected: 
	/* strain */
	dSymMatrixT fStrain3D;

	/* stress/modulus */
	dMatrixT fModulus; 
	dSymMatrixT fStress; 

};

} // namespace Tahoe 
#endif /*_SS_VISCO_H_*/
