/* $Id: ElementSupportT.h,v 1.30 2004-07-15 08:25:44 paklein Exp $ */
#ifndef _ELEMENT_SUPPORT_T_H_
#define _ELEMENT_SUPPORT_T_H_

/* direct members */
#include "BasicSupportT.h"

/* direct members */
#include "dArray2DT.h"
#ifdef _FRACTURE_INTERFACE_LIBRARY_
#include "StringT.h"
#endif

namespace Tahoe {

/* forward declarations */
class dMatrixT;
class GroupAverageT;
class dArrayT;
class iArrayT;
class StringT;
class OutputSetT;
class LocalArrayT;

/** support for the ElementBaseT class hierarchy. A limited interface to get 
 * information in and out of an ElementBaseT */
class ElementSupportT: public BasicSupportT
{
public:

	/** constructor */
	ElementSupportT(void);

	/** destructor */
	~ElementSupportT(void);

	/** return the index of the given element group */
	int ElementGroupNumber(const ElementBaseT* group) const;

	/** \name assembly functions */
	/*@{*/
	void AssembleLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& eqnos) const;
	void AssembleRHS(int group, const nArrayT<double>& elRes, const nArrayT<int>& eqnos) const;
	void AssembleLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& row_eqnos, const nArrayT<int>& col_eqnos) const;
	void AssembleLHS(int group, const nArrayT<double>& diagonal_elMat, const nArrayT<int>& eqnos) const;
	/*@}*/

	/** \name nodal averaging */
	/*@{*/
	/** initialize work space to the number of values to be averaged */
	void ResetAverage(int n_values) const;

	/** assemble values 
	 * \param nodes list of nodes for the values being assembled: [nnd] 
	 * \param vals values to be assembled: [nnd] x [nvals] */
	void AssembleAverage(const iArrayT& nodes, const dArray2DT& vals) const;

	/** average assembled values and return the array of averages 
	 * values: [nnd] x [nvals] */
	const dArray2DT& OutputAverage(void) const;

	/** return averaged values for the nodes with assembled values. Returned
	 * nodes are ordered by increasing node number */
	void OutputUsedAverage(dArray2DT& average_values) const;
	/*@}*/

	/** \name input/output */
	/*@{*/
#ifdef _FRACTURE_INTERFACE_LIBRARY_
	/** register the output set. returns the ID that should be used with
	 * ElementSupport::WriteOutput. or SIERRA interface, we're just holding 
	 * the data till they want it */
	int RegisterOutput(ArrayT<StringT>& n_labels, ArrayT<StringT>& e_labels);
#endif

	/** write results for a single output set
	 * \param ID output set ID for the given data
	 * \param n_values nodal output values
	 * \param e_values element output values */
	void WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values) const;

	/** return true if output is going to be written for the current time step */
	bool WriteOutput(void) const { return BasicSupportT::WriteOutput(); };

	/** write a snapshot */
	void WriteOutput(const StringT& file, const dArray2DT& coords, const iArrayT& node_map,
		const dArray2DT& values, const ArrayT<StringT>& labels) const;
	/*@}*/

#ifdef _FRACTURE_INTERFACE_LIBRARY_
	/** Parameters normally read from input stream must be passed through ElementSupport */
	enum CodeT { kGeometryCode = 0, /**< Topology of surface element */
	    		    kNumIntPts = 1, /**< Number of integration points */
	             kCloseSurface = 2, /**< Initially close cohesive surfaces? */
				   kOutputArea = 3, /**< Output fracture area */
	             kMaterialCode = 4};/**< Which cohesive law to use */ 

	/** set the number of nodes in the fracture interface */
	void SetNumNodes(int nn);	

	/** set the memory used to hold the reference configuration */
	void SetInitialCoordinates(dArray2DT *InitialCoords);
	void SetInitialCoordinates(double *InitialCoords);
	
	/** set the memory used to hold the current configuration */
	void SetCurrentCoordinates(dArray2DT *CurrentCoords);
	void SetCurrentCoordinates(double *CurrentCoords);

	/** use the displacements and reference configuration to update the current one */
	void UpdateCurrentCoordinates(double *displacements);

	/** set the time step in the fracture interface */
	void SetTimeStep(double dt);

	/** set the model manager in the fracture interface */
	void SetModelManager(ModelManagerT* modelManager);

	/** set the number of elements in the fracture interface */
	void SetNumElements(int nelem);

	/** accessor for the number of elements in the fracture interface */
	int NumElements(void) const;

	/** accessor for element floating point input when streams are not available */
	dArrayT *FloatInput(void) const;
	
	/** accessor for element integer input when streams are not available */
	iArrayT *IntInput(void) const;
	
	/** generate equation numbers based on connectivity information */
	void SetEqnos(int *conn, const int& nelem, const int& nElemNodes, const int&nNodes);

	dArrayT& Residual(void) const { return *fResidual; };
	
	dMatrixT& Stiffness(void) const { return *fStiffness; };
	
	void SetMaterialInput(double *inputFloats, int length);
	
	void SetElementInput(int *inputInts, int length);
		
	int ReturnInputInt(CodeT label);
	
	double *StateVariableArray(void);
	
	void SetStateVariableArray(double *incomingArray);
	
	void SetBlockID(StringT &Id);
	
	StringT& BlockID(void);
	
	void OutputSize(int& nNodeOutputVars, int& nElemOutputVars);
	
	void SetOutputCodes(iArrayT& fNodalOutputCodes, iArrayT& fElementOutputCodes);
	
	void SetOutputPointers(double *nodalOutput, double *elemOutput);
	
	void SetResidual(double *nodalForces);

#endif /* _FRACTURE_INTERFACE_LIBRARY_ */

#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/** \name command line flags */
	/*@{*/
	const ArrayT<StringT>& Argv(void) const;
	bool CommandLineOption(const char* str) const;

	/** returns the index of the requested option or -1 if not bound */
	bool CommandLineOption(const char* str, int& index) const;
	/*@}*/
#endif

private:

#ifdef _FRACTURE_INTERFACE_LIBRARY_	 
	dArrayT *fResidual;
	dMatrixT *fStiffness;
	
	dArray2DT *fInitialCoordinates;
	dArray2DT *fCurrentCoordinates;

	int fElem;	
	GroupAverageT* fGroupAverage;

	dArrayT *fparams;
	iArrayT *iparams;
	
	iArrayT *ieqnos;
	
	double *fStateVars, *fNodalOutput, *fElemOutput;

	StringT sBlockID;
	
	ArrayT<StringT> fNodeOutputLabels;
	ArrayT<StringT> fElemOutputLabels;
#endif
	
};

#ifdef _FRACTURE_INTERFACE_LIBRARY_
inline int ElementSupportT::NumElements(void) const { return fElem; }
inline dArrayT *ElementSupportT::FloatInput(void) const { return fparams; }
inline iArrayT *ElementSupportT::IntInput(void) const { return iparams; }
#endif

inline void ElementSupportT::AssembleLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& row_eqnos,
	const nArrayT<int>& col_eqnos) const 
{ 
	BasicSupportT::AssembleLHS(group, elMat, row_eqnos, col_eqnos); 
}

inline void ElementSupportT::AssembleLHS(int group, const nArrayT<double>& diagonal_elMat, const nArrayT<int>& eqnos) const
{
	BasicSupportT::AssembleLHS(group, diagonal_elMat, eqnos); 
};

} /* namespace Tahoe */

#endif /* _ELEMENT_SUPPORT_T_H_ */
