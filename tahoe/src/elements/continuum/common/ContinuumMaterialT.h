/* $Id: ContinuumMaterialT.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (11/20/1996)                                          */
/* Defines the interface for continuum materials.                         */

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

class ContinuumMaterialT
{
public:

	/* constructor */
	ContinuumMaterialT(const ContinuumElementT& element);

	/* destructor */
	virtual ~ContinuumMaterialT(void);

	/* form of tangent matrix (symmetric by default) */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* element information */
	const ContinuumElementT& ContinuumElement(void) const;

	/* number of nodal degrees of freedom */
	int NumDOF(void) const;

	/* integration information */
	int NumIP(void) const;
	int CurrIP(void) const;

	/* element card data */
	int NumElements(void) const;
	int CurrElementNumber(void) const;
	ElementCardT& ElementCard(int card) const;
	ElementCardT& CurrentElement(void) const;

	/* initialization */
	virtual void Initialize(void);

	/* storage initialization */
	virtual bool NeedsPointInitialization(void) const; // false by default
	virtual void PointInitialize(void);                // ip at a time

	/* apply pre-conditions at the current time step */
	virtual void InitStep(void);

	/* update/reset internal variables */
	virtual void UpdateHistory(void); // element at a time
	virtual void ResetHistory(void);  // element at a time

	/* print parameters */
	virtual void Print(ostream& out) const = 0;
	virtual void PrintName(ostream& out) const = 0;

	/* returns the number of variables computed for nodal extrapolation
	 * during for element output, ie. internal variables. Returns 0
	 * by default */
	virtual int NumOutputVariables(void) const; // 0 by default
	virtual void OutputLabels(ArrayT<StringT>& labels) const; // none by default
	virtual void ComputeOutput(dArrayT& output);
	
protected:

	/* element group */
	const ContinuumElementT& fContinuumElement;
	
	/* nodal degrees of freedom */
	int fNumDOF;
	
	/* integration point info */
	int fNumIP;
	const int& fCurrIP;  	
};

/* inlines */
inline int ContinuumMaterialT::NumDOF(void) const { return fNumDOF; }
inline int ContinuumMaterialT::NumIP(void) const { return fNumIP; }
inline int ContinuumMaterialT::CurrIP(void) const { return fCurrIP; }
inline const ContinuumElementT& ContinuumMaterialT::ContinuumElement(void) const
{ return fContinuumElement; }

#endif /* _CONTINUUM_MATERIAL_T_H_ */
