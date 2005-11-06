/* $Id: RodT.h,v 1.20 2005-11-06 00:37:58 paklein Exp $ */
/* created: paklein (10/22/1996) */

#ifndef _ROD_T_H_
#define _ROD_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "RodMaterialT.h"
#include "LocalArrayT.h"
#include "dSymMatrixT.h"
#include "pArrayT.h"

namespace Tahoe {

/** \note the RodT class doesn't provide complete support for the
 * different time integration schemes implemented using the
 * controller classes. need to add something like the
 * ComputeEffectiveDVA functions from the continuum element
 * classes to provide contributions to the global equations
 * which are consistent with the time integration algorithm.
 * PAK (05/30/1999) */
class RodT: public ElementBaseT
{
public:

	/** constructor */
	RodT(const ElementSupportT& support);

	/** form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* NOT implemented. Returns an zero force vector */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);
			
	/* returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void);
	
	/* writing output */
	virtual void RegisterOutput(void);
	virtual void WriteOutput(void);

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);

	/* initialize/finalize time increment */
	virtual void CloseStep(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	ParameterInterfaceT* NewSub(const StringT& name) const;
		
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected: /* for derived classes only */
	 	
	/* called by FormRHS and FormLHS */
	virtual void LHSDriver(GlobalT::SystemTypeT);
	virtual void RHSDriver(void);

	/* increment current element */
	virtual bool NextElement(void);

	/** return true if connectivities are changing */
	virtual bool ChangingGeometry(void) const { return false; };

private: /* MD related computational functions */

	void ComputeInstKE(void);
	void ComputeAvgKE(void);
	void ComputeInstPE(void);
	void ComputeAvgPE(void);
	void ComputeInstTotalE(void);
	void ComputeAvgTotalE(void);
	void ComputeInstTemperature(void);
	void ComputeAvgTemperature(void);
	void ComputeInstPressure(void);
	void ComputeAvgPressure(void);
	int PrintMDToFile(void);

	/* Hardy-related functions */
	void ComputeHardyStress(void);
	void ComputeHardyHeatFlux(void);
	void LocalizationFunction(double distance, double radius, double& gauss);
	void PairIntegral(void);

protected:

	/* material data */
	pArrayT<RodMaterialT*> fMaterialsList; 	
	RodMaterialT*	       fCurrMaterial;

	/** list of nodes used by the group both in and not in
	 * current bonding configuration */
	iArrayT fGroupNodes;

private:

	/** output diagnostic data */
	bool fOutputDiagnostic;

	/** \name work space */
	/*@{*/
	/** constant matrix needed to compute the stiffness */
	dMatrixT fOneOne;

	/** current pair vector */
	dArrayT fBond;

	/** reference pair vector */
	dArrayT fBond0;

	/** current coordinates for one pair bond */
//	dArray2DT fPairCoords;

	/** local accelerations */
	LocalArrayT fLocAcc;
	dArrayT fNEE_vec;
	/*@}*/

	/** \name MD related variables */
	/*@{*/
	double fKb;
	double fInstKE, fInstPE, fInstTotalE, fInstTemp, fInstPressure;
	double fAvgKE, fAvgPE, fAvgTotalE, fAvgTemp, fAvgPressure;
	double fSumKE, fSumPE, fSumTotalE, fSumTemp, fSumPressure;
	LocalArrayT fLocVel;
	dMatrixT fHardyStress;
	dArrayT fHardyHeatFlux;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _ROD_T_H_ */
