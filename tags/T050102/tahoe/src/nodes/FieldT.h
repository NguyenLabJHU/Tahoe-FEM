/* $Id: FieldT.h,v 1.1 2002-04-21 07:13:33 paklein Exp $ */

#ifndef _FIELD_T_H_
#define _FIELD_T_H_

/* direct members */
#include "ArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "StringT.h"
#include "IC_CardT.h"
#include "KBC_CardT.h"
#include "FBC_CardT.h"
#include "pArrayT.h"

/* forward declarations */
class nControllerT;
class KBC_ControllerT;
class FBC_ControllerT;

/** field of unknowns */
class FieldT
{
public:
	
	/** constructor */
	FieldT(const StringT& name, nControllerT& controller, int ndof);

	/** \name initialization */
	/*@{*/
	/** set field labels */
	void SetLabels(const ArrayT<StringT>& labels);

	/** set number of nodes. Allocates memory. */
	void Dimension(int nnd);
	/*@}*/
	
	/** \name accessors */
	/*@{*/
	/** field name */
	const StringT& Name(void) const { return fName; };
	
	/** the field labels */
	const ArrayT<StringT>& Labels(void) const { return fLabels; };

	/** the specified derivative of the field */ 
	const dArray2DT& Field(int order) const { return fField[order]; };
	/*@}*/

	/** \name time integration */
	/*@{*/
	/** beginning of time series */
	void InitialCondition(void);
	
	/** apply predictor to all degrees of freedom */
	void InitStep(void);

	/** update history */
	void CloseStep(void);

	/** update the active degrees of freedom */
	void Update(const dArrayT& update, int eq_start, int eq_stop);

	/** reset displacements (and configuration to the last known solution) */
	void ResetStep(void);
	/*@}*/

private:

	/** name */
	StringT fName;

	/** time integrator */
	nControllerT& fnController;
	
	/** the field [nderiv]: [nnd] x [ndof] */
	ArrayT<dArray2DT> fField;

	/** field history */
	ArrayT<dArray2DT> fField_last;

	/** equation numbers [nnd] x [ndof] */
	iArray2DT fEqnos;
	
	/** field dof labels [ndof] */
	ArrayT<StringT> fLabels;
	
	/** \name initial and boundary conditions */
	/*@{*/
	/** initial conditions */
	ArrayT<IC_CardT> fIC;
	  	
	/** kinematic boundary conditions */
	ArrayT<KBC_CardT> fKBC;

	/** special KBC objects */
	pArrayT<KBC_ControllerT*> fKBCControllers;

	/** force boundary conditions */
	ArrayT<FBC_CardT> fFBC;

	/** special FBC objects */
	pArrayT<FBC_ControllerT*> fFBCControllers;
	/*@}*/
};

#endif /* _FIELD_T_H_ */
