/* $Id: BasicFieldT.h,v 1.8 2003-08-14 05:34:21 paklein Exp $ */
#ifndef _BASIC_FIELD_T_H_
#define _BASIC_FIELD_T_H_

/* direct members */
#include "StringT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "nArray2DGroupT.h"

namespace Tahoe {

/* forward declarations */
class iArrayT;

/** basic container for field data */
class BasicFieldT
{
public:

	/** constructor */
	BasicFieldT(void);

	/** destructor */
	virtual ~BasicFieldT(void) {};

	/** \name initialization */
	/*@{*/
	/** initialize field */
	void Initialize(const StringT& name, int ndof, int order);
	
	/** set field labels */
	void SetLabels(const ArrayT<StringT>& labels);

	/** set number of nodes. Resizes all arrays to the new number of nodes. 
	 * Excess nodes at the tail of all arrays are discarded. Additional nodes 
	 * added to all arrays is not initialized.
	 * \param nnd number of nodes 
	 * \param copy_in if true, values that fit are copied in. Otherwise,
	 *        new array is unititalized */
	virtual void Dimension(int nnd, bool copy_in);
	
	/** set all field values to 0.0 */
	virtual void Clear(void);
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
	void WriteEquationNumbers(ostream& out, const ArrayT<int>* node_map) const;
	/*@}*/

protected:

	/** \name array registration
	 * Register arrays with the dynamic memory managers */
	/*@{*/
	/** register an double array */
	void RegisterArray2D(nArray2DT<double>& array2D);

	/** register an integer array */
	void RegisterArray2D(nArray2DT<int>& array2D);
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

private:

	/** \name dynamic memory managers
	 * Handle dynamic resizing of arrays with dimension: [nnd] x [ndof] */
	/*@{*/
	nArray2DGroupT<double> fdArray2DGroup;
	nArray2DGroupT<int>    fiArray2DGroup;
	/*@}*/	
};

/* register an double array */
inline void BasicFieldT::RegisterArray2D(nArray2DT<double>& array2D)
{
	fdArray2DGroup.Register(array2D);
}

/* register an integer array */
inline void BasicFieldT::RegisterArray2D(nArray2DT<int>& array2D)
{
	fiArray2DGroup.Register(array2D);
}

} // namespace Tahoe 
#endif /* _BASIC_FIELD_T_H_ */
