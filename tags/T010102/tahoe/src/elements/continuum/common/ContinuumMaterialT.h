/* $Id: ContinuumMaterialT.h,v 1.4 2001-10-24 02:11:24 paklein Exp $ */
/* created: paklein (11/20/1996) */

#ifndef _CONTINUUM_MATERIAL_T_H_
#define _CONTINUUM_MATERIAL_T_H_

#include "Environment.h"
#include "GlobalT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ElementCardT;
class dArrayT;
template <class TYPE> class ArrayT;
class StringT;
class ContinuumElementT;

/** interface for continuum materials. */
class ContinuumMaterialT
{
public:

	/** constructor.
	 * \param element reference to the host element */
	ContinuumMaterialT(const ContinuumElementT& element);

	/** destructor */
	virtual ~ContinuumMaterialT(void);

	/** form of tangent matrix. \return symmetric by default */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/** reference to the host element */
	const ContinuumElementT& ContinuumElement(void) const;

	/** number of degrees of freedom (per node) in the host
	 * element group. */
	int NumDOF(void) const;

	/** number of spatial dimensions in the host element group. */
	int NumSD(void) const;

	/** the total number of integration points per element in the
	 * host element group. */
	int NumIP(void) const;

	/** reference to the current integration point within the
	 * element of evaluation. */
	const int& CurrIP(void) const;

	/** return the total number of elements in the host element
	 * group. */
	int NumElements(void) const;

	/** return the number of the current element of evaluation. */
	int CurrElementNumber(void) const;

	/** reference to the ElementCardT for the  specified element. */
	ElementCardT& ElementCard(int card) const;

	/** reference to the ElementCardT for the current element of
	 * evaluation */
	ElementCardT& CurrentElement(void) const;

	/** initialization. Called immediately after constructor to allow
	 * class specific initializations. */
	virtual void Initialize(void);

	/** return true if model needs ContinuumMaterialT::PointInitialize
	 * to be called for every integration point of every element as
	 * part of the model initialization. \return false by default. */
	virtual bool NeedsPointInitialization(void) const;
	
	/** model initialization. Called per integration point for every
	 * element using the model. Deformation variables are available
	 * during this call. */
	virtual void PointInitialize(void);

	/** apply pre-conditions at the current time step. Called once for
	 * the model at the beginning of a time increment */
	virtual void InitStep(void);

	/** finalize the current time step. Called once for the model at 
	 * the end of a time increment */
	virtual void CloseStep(void);

	/** update internal variables. Called once per element for all
	 * elements using the model, hence no deformation variables are
	 * available during this call. */
	virtual void UpdateHistory(void);

	/** restore internal variables to their state at the beginning of
	 * the current time increment. Called once per element for all
	 * elements using the model, hence no deformation variables are
	 * available during this call. */
	virtual void ResetHistory(void);

	/** write parameters to the output stream. */
	virtual void Print(ostream& out) const = 0;

	/** write the model name to the output stream. */
	virtual void PrintName(ostream& out) const = 0;

	/** return the number of constitutive model output parameters
	 * per evaluation point. Used by the host element group in
	 * conjunction with ContinuumMaterialT::OutputLabels and
	 * ContinuumMaterialT::ComputeOutput to collect model variables
	 * for output. \return zero by default */
	virtual int NumOutputVariables(void) const;

	/** return the labels for model output parameters
	 * per evaluation point. Used by the host element group in
	 * conjunction with ContinuumMaterialT::NumOutputVariables and
	 * ContinuumMaterialT::ComputeOutput to collect model variables
	 * for output.
	 * \param labels destination for the variable labels. Returns
	 *        with length ContinuumMaterialT::NumOutputVariables */
	virtual void OutputLabels(ArrayT<StringT>& labels) const;

	/** return material output variables. Used by the host element group 
	 * in conjunction with ContinuumMaterialT::NumOutputVariables and
	 * ContinuumMaterialT::OutputLabels to collect model variables
	 * for output. Called per integration point. Deformation variables
	 * are available.
	 * \param output destination for the output. Must be passed in with
	 *        length ContinuumMaterialT::NumOutputVariables */
	virtual void ComputeOutput(dArrayT& output);

	/** returns true if two materials have compatible output variables.
	 * Used by the host element to determine whether the two material
	 * models can be used within the same host element group when
	 * requesting model-specific, materials output. */
	static bool CompatibleOutput(const ContinuumMaterialT& m1, const ContinuumMaterialT& m2);
	
protected:

	/** host element group */
	const ContinuumElementT& fContinuumElement;
	
	/** number of degrees of freedom */
	int fNumDOF;

	/** number of degrees of spatial dimensions */
	int fNumSD;
	
	/** number of integration points */
	int fNumIP;

	/** reference to the current integration point for the
	 * current element of evaluation. */
	const int& fCurrIP;
};

/* inlines */
inline int ContinuumMaterialT::NumDOF(void) const { return fNumDOF; }
inline int ContinuumMaterialT::NumSD(void) const { return fNumSD; }
inline int ContinuumMaterialT::NumIP(void) const { return fNumIP; }
inline const int& ContinuumMaterialT::CurrIP(void) const { return fCurrIP; }
inline const ContinuumElementT& ContinuumMaterialT::ContinuumElement(void) const
{ return fContinuumElement; }

#endif /* _CONTINUUM_MATERIAL_T_H_ */
