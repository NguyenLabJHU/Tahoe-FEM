/* $Id: BasicFieldT.h,v 1.4 2002-07-05 22:28:30 paklein Exp $ */

#ifndef _BASIC_FIELD_T_H_
#define _BASIC_FIELD_T_H_

/* direct members */
#include "StringT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"

namespace Tahoe {

/* forward declarations */
class iArrayT;

/** basic container for field data */
class BasicFieldT
{
public:

	/** constructor */
	BasicFieldT(const StringT& name, int ndof, int order);

	/** \name initialization */
	/*@{*/
	/** set field labels */
	void SetLabels(const ArrayT<StringT>& labels);

	/** set number of nodes. (Re-)allocates memory. */
	void Dimension(int nnd);
	/*@}*/
	
	/** \name accessors */
	/*@{*/
	/** field name */
	const StringT& Name(void) const { return fName; };
	
	/** the field labels */
	const ArrayT<StringT>& Labels(void) const { return fLabels; };

	/** reference to the specified derivative of the field */ 
	dArray2DT& operator[](int order) { return fField[order]; };

	/** const reference to the specified derivative of the field */ 
	const dArray2DT& operator[](int order) const { return fField[order]; };
	
	/** number of nodes */
	int NumNodes(void) const { return fEqnos.MajorDim(); };
	
	/** number of degrees of freedom per node */
	int NumDOF(void) const { return fEqnos.MinorDim(); };
	
	/** number of time derivatives stored by the field */
	int Order(void) const { return fField.Length() - 1; };
	/*@}*/

	/** \name equation numbers */
	/*@{*/
	/** return the equation of the degree of freedom of the specified node */
	int EquationNumber(int node, int dof) const { return fEqnos(node, dof); };
	
	/** const access to the equation numbers */
	const iArray2DT& Equations(void) const { return fEqnos; };

	/** non-const access to the equation numbers. Modify these at your own risk */
	iArray2DT& Equations(void) { return fEqnos; };

	/** write field equation numbers to the output stream */
	void WriteEquationNumbers(ostream& out, const iArrayT* node_map) const;
	/*@}*/

protected:

	/** name */
	StringT fName;
	
	/** the field [nderiv]: [nnd] x [ndof] */
	ArrayT<dArray2DT> fField;

	/** field dof labels [ndof] */
	ArrayT<StringT> fLabels;	

	/** equation array: [nnd] x [ndof] */
	iArray2DT fEqnos;
};

} // namespace Tahoe 
#endif /* _BASIC_FIELD_T_H_ */
