/* $Id: FieldT.h,v 1.1.2.1 2002-04-22 07:06:05 paklein Exp $ */

#ifndef _FIELD_T_H_
#define _FIELD_T_H_

/* direct members */
#include "iArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "StringT.h"
#include "IC_CardT.h"
#include "KBC_CardT.h"
#include "FBC_CardT.h"
#include "pArrayT.h"

/* forward declarations */
class LocalArrayT;
class nControllerT;
class KBC_ControllerT;
class FBC_ControllerT;
template <class TYPE> class RaggedArray2DT;
class ifstreamT;
class ofstreamT;

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

	/** set the number of external nodes. Dimensions work space for external nodes */
	void DimensionExternal(int nnd);

	/** register the local array with its source */
	void RegisterLocal(LocalArrayT& array) const;

	/** set the group number */
	void SetGroup(int group) { fGroup = group; };	
	/*@}*/
	
	/** \name accessors */
	/*@{*/
	/** field name */
	const StringT& Name(void) const { return fName; };

	/** set the group number */
	int Group(void) { return fGroup; };	
	
	/** the field labels */
	const ArrayT<StringT>& Labels(void) const { return fLabels; };

	/** the specified derivative of the field */ 
	const dArray2DT& Field(int order) const { return fField[order]; };
	
	/** number of nodes */
	int NumNodes(void) const { return fEqnos.MajorDim(); };
	
	/** number of degrees of freedom per node */
	int NumDOF(void) const { return fEqnos.MinorDim(); };
	
	/** number of time derivatives stored by the field */
	int Order(void) const { return fField.Length(); };
	
	/** the equation numbers array */
	iArray2DT& Equations(void) { return fEqnos; };
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

	/** \name collecting equation numbers */
	/*@{*/
	/** collection equation numbers.
	 * \param nodes element connectivities: [nel] x [nen]
	 * \param eqnos destination for equation numbers: [nel] x [nen*ndof] */
	void SetLocalEqnos(const iArray2DT& nodes, iArray2DT& eqnos) const;

	/** collection equation numbers. Connectivities are passed in a RaggedArray2DT, 
	 * which allows an arbitrary number of nodes per element.
	 * \param nodes element connectivities: [nel] x [nen_i]
	 * \param eqnos destination for equation numbers: [nel] x [nen_i*ndof] */
	void SetLocalEqnos(const RaggedArray2DT<int>& nodes, RaggedArray2DT<int>& eqnos) const;

	/** collect equation numbers */
	void SetLocalEqnos(const iArrayT& tags, iArray2DT& eqnos) const;
	/*@}*/

	/** \name restart functions
	 * The restart functions should read/write any data that overrides the 
	 * default values */
	/*@{*/ 
	void ReadRestart(ifstreamT& in);
	void WriteRestart(ofstreamT& out) const;
	/*@}*/ 

	/** \name accessors to data for external nodes */
	/*@{*/
	iArray2DT& ExternalEquations(void) { return fExEqnos; };
	dArray2DT& ExternalUpdate(void) { return fExUpdate; };
	/*@}*/

	/** write field equation numbers to the output stream */
	void WriteEquationNumbers(ostream& out, const iArrayT* node_map) const;

private:

	/** apply the IC_CardT to the field */
	void Apply_IC(const IC_CardT& card);

private:

	/** name */
	StringT fName;
	
	/** solution set number */
	int fGroup;

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
	pArrayT<KBC_ControllerT*> fKBC_Controllers;

	/** force boundary conditions */
	ArrayT<FBC_CardT> fFBC;

	/** special FBC objects */
	pArrayT<FBC_ControllerT*> fFBC_Controllers;
	/*@}*/

	/** \name external nodes */
	/*@{*/
	iArray2DT fExEqnos;
	dArray2DT fExUpdate;
	/*@}*/
};

/* inlines */
inline void FieldT::SetLocalEqnos(const iArrayT& tags,
	iArray2DT& eqnos) const
{
	eqnos.RowCollect(tags,fEqnos);
}

#endif /* _FIELD_T_H_ */
