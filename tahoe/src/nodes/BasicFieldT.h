/* $Id: BasicFieldT.h,v 1.1.2.1 2002-04-24 01:29:25 paklein Exp $ */

#ifndef _BASIC_FIELD_T_H_
#define _BASIC_FIELD_T_H_

/* direct members */
#include "StringT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"

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
	int Order(void) const { return fField.Length(); };
	/*@}*/

	/** \name equation numbers */
	/*@{*/
	/** determine the number of equations for the field */
	int NumEquations(void) const;
	
	/** the global equation numbers array. Resets the cached count of active equations */
	iArray2DT& Equations(void);
	
	/** const access to the equation numbers */
	const iArray2DT& Equations(void) const { return fEqnos; };

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

	/** \name equation numbers */
	/*@{*/
	/** equation array: [nnd] x [ndof] */
	iArray2DT fEqnos;

	/** (cached) count of the number of active equations */
	int fNumActiveEquations;
	/*@{*/	
};

/* the global equation numbers array. Resets the cached count of active equations */
inline iArray2DT& BasicFieldT::Equations(void)
{
	fNumActiveEquations = -1;
	return fEqnos;
}

/* determine the number of equations for the field */
int BasicFieldT::NumEquations(void) const
{
	/* recalculate */
	if (fNumActiveEquations == -1)
	{
		int neq = 0;
		int *peq = fEqnos.Pointer();
		int len  = fEqnos.Length();
		for (int i = 0; i < len; i++)
			if (*peq++ > 0) neq++;
	
		/* not so const */
		const_cast<BasicFieldT*>(this)->fNumActiveEquations = neq;
	}
	return fNumActiveEquations;
}

#endif /* _BASIC_FIELD_T_H_ */
