/* $Id: QLV_Nfibers.h,v 1.1 2008-09-25 13:32:32 thao Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _QLV_Nfibers_
#define _QLV_Nfibers_ 
 
/* base class */
#include "FSFiberMatT.h"
#include "C1FunctionT.h"
#include "PotentialT.h"
#include "dArray2DT.h"
#include "Array2DT.h"

namespace Tahoe {

class QLV_Nfibers: public FSFiberMatT
{
   public:
  
	/* constructor/destructor */
	QLV_Nfibers(void);
	
/* destructor */
	~QLV_Nfibers(void);
	
	/* strain energy density */
	virtual double StrainEnergyDensity(void);

	virtual GlobalT::SystemTypeT TangentType(void) const{return GlobalT::kSymmetric;};
	virtual bool Need_F_last(void) const { return true; };

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/* material description */
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress

	/*compute output variables*/
	virtual int NumOutputVariables() const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/*Initialize history variable*/
	virtual bool NeedsPointInitialization(void) const {return true;}; 
	virtual void PointInitialize(void);              

	/* update/reset internal variables */
	/** return true if the material has history variables */
	virtual bool HasHistory(void) const { return true; };
	virtual void UpdateHistory(void); // element at a time
	virtual void ResetHistory(void);  // element at a time
	
	void Load(ElementCardT& element, int ip);
	void Store(ElementCardT& element, int ip);
	
		/* Dimension internal state variables*/
	/*derived class must call RGViscoelaticity::SetStateVariables(fNumProcess)
	  to dimension internal state variable arrays if fNumProcess > 1 (default value)*/
	void SetStateVariables (const int num_mat_process,const int num_fib_process);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;
	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** describe the parameters */
	virtual void DefineParameters(ParameterListT& list) const;
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** types of analysis */
	enum FiberPot {
	         kStandard = 0,      // W(I) = k(I-1)^2
		      kExponential = 1,  // W(I) = k1/k2 (exp(k2(I-1)) - k2 I) veronda westmann
		     kExpQuad = 2        // W(I) = k1/2k2 (exp(k2(I-1)^2) - 2 k2 I) holzapfel gasser
		   };
		

	/****************************** fiber ***********************************************/
	/*computes eq fiber stress in local frame*/
	virtual void ComputeFiberStress (const dSymMatrixT& Stretch, dSymMatrixT& Stress);
	
	/*computes eq fiber moduli in local frame*/
	virtual void ComputeFiberMod (const dSymMatrixT& Stretch, dSymMatrixT& dSymMatrixT, dMatrixT& Mod);

	/*computes fiber stress in local frame*/
	/*redefined from base class*/
	virtual void ComputeFiberOverStress (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_n,  dSymMatrixT& Stress, dArrayT& h, 
		const dArrayT& h_n,	const int process_index, const int fillmode = dSymMatrixT::kAccumulate);
	
	/*computes  moduli in local frame*/
	virtual void ComputeFiberAlgMod (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_n, dSymMatrixT& Stress, dArrayT& h, 
		const dArrayT& h_n, dMatrixT& Mod, const int process_index, const int fillmode = dSymMatrixT::kAccumulate);

	/****************************** matrix ***********************************************/
	/*computes eq isotropic matrix stress*/
	virtual void ComputeMatrixStress (const dSymMatrixT& Stretch, dSymMatrixT& Stress);

	/*computes eq matrix moduli*/
	virtual void ComputeMatrixMod (const dSymMatrixT& Stretch, dSymMatrixT& Stress, dMatrixT& Mod);

	/*computes isotropic matrix stress*/
	virtual void ComputeMatrixOverStress (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_n, 
		dSymMatrixT& H, const dSymMatrixT& H_n, const int process_index, const int fillmode = dSymMatrixT::kOverwrite);

	/*computes matrix moduli*/
	virtual void ComputeMatrixAlgMod (const dSymMatrixT& Stretch, const dSymMatrixT& Stretch_n,
		dSymMatrixT& H, const dSymMatrixT& H_n, dMatrixT& Mod, const int process_index, const int fillmode = dSymMatrixT::kOverwrite);
protected:
	

	/* number of nonequilibrium processes*/
	/* must be set in derived classes before TakeParameterList is called*/
	/* default value is 1*/
	
	int fNumFibProcess;
	int fNumMatProcess;
	bool fsame;
 	/* work space needed to compute fiber viscoelastic behavior*/
	/*dimension fNumFibProcess + 1*/
//	int fPotType_f;
	Array2DT<C1FunctionT*> fPot_f;
	dArray2DT ftime_f;
	
	/*moduli for Mooney Rivlin Potential, n_process x 2 (c1, c2)*/
	ArrayT<PotentialT*> fPot_m;
	dArrayT ftime_m;
				
	/*number of state variables*/
	int fnstatev;
	
	/* Internal state variables array*/
	dArrayT fstatev;

	/*current values*/
	/*deviatoric*/
	ArrayT<dSymMatrixT> fHm;
	ArrayT<dArrayT> fhf;
	dArrayT ftau;

	/*preceding values*/		 
	/*deviatoric*/
	ArrayT<dSymMatrixT> fHm_n;
	ArrayT<dArrayT> fhf_n;
	dArrayT ftau_n;
	
	/*stretch*/
	dSymMatrixT fb;
	dArrayT fEigs;

	/*previous deformation tensor*/
	dSymMatrixT fC_n;
	dSymMatrixT fb_n;
	dSymMatrixT fFiberStretch_n;
	
	SpectralDecompT fSpectralDecompSpat;


	dMatrixT fModMat;
	dMatrixT fMod3;
	dSymMatrixT fdtau_dep;
};
	
/*equilibrium values*/
}
#endif /* _QLV_Nfibers_ */
