/* $Id: MFGPElementSupportT.h,v 1.4 2005-06-08 21:40:59 paklein Exp $ */
#ifndef _MFGP_SUPPORT_T_H_
#define _MFGP_SUPPORT_T_H_

/* base classes */
#include "ParameterInterfaceT.h"

/* direct members */
#include "iArrayT.h"
#include "dMatrixT.h"
#include "LocalArrayGroupT.h"
#include "RaggedArray2DT.h"
#include "nArrayGroupT.h"
#include "nMatrixGroupT.h"
#include "dArray2DT.h"
#include "IOBaseT.h"
#include "MeshFreeT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class MeshFreeShapeFunctionT;
class MeshFreeNodalShapeFunctionT;
class MeshFreeSupportT;
class ElementCardT;
class StringT;
class ElementBaseT;
class ModelManagerT;

/** support for meshfree calculations */
class MFGPElementSupportT: public ParameterInterfaceT
{
public:

	/** constructor */
	MFGPElementSupportT(void);

	/** destructor */
	virtual ~MFGPElementSupportT(void) { };
	
	/** accessors */
	MeshFreeSupportT& MeshFreeSupport(void) const;

	/** set pointer to the shape functions */
	void SetShape(MeshFreeShapeFunctionT* mf_shapes);

	/** set pointer to the shape functions 
	  * assume that first and second parameters of the SetShape are
	  * displacement and plastic multiplier shape functions respectively
	  * only one call to this SetShape is needed */
	void SetShape(MeshFreeShapeFunctionT* mf_shapes_displ, MeshFreeShapeFunctionT* mf_shapes_plast);

	/** \name accessors */
	/*@{*/
	/** return array of interpolant points */
	const iArrayT& InterpolantNodes(void) { return fAllFENodes; };
	
	/** returnt the array of off-grid points */
	const iArrayT& OffGridNodes(void) { return fOffGridNodes; };

	/** element nodes */
	const RaggedArray2DT<int>& ElementNodes(void);

	/** element equations */
	RaggedArray2DT<int>& ElementEquations(void) { return fElemEqnosEX; };
	/*@}*/

	/** \name memory managers for variable numbers of element nodes */
	/*@{*/
	/** accessor */
	int NumElementNodes(void) const;

	/** resize number of field nodes, returning number of element nodes */
	int SetElementNodes(int element);

	/** redimension the off-diagonal, non-square matrices */
	void SetOffDiagMatrix(int element);
	
	/** register local arrays */
	void Register(LocalArrayT& u);

	/** register vectors that need to be length number of element equations */
	void Register(dArrayT& v);

	/** register matrix that need square dimensions of element equations */
	void Register(dMatrixT& m);
	
	/** register matrix that need non-square dimensions of element equations */
	// if field_1 and field_2 have same dof, use MeshFreeElementSupportT::Register
	void Register(dMatrixT& m, LocalArrayT& field_1, LocalArrayT& field_2);
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** mark "dead" cells (no active equations) returns the number of active */
	int MarkActiveCells(AutoArrayT<ElementCardT>& elem_cards);

	/** initialization of meshless information. This method must be called once after 
	 * a call to MeshFreeElementSupportT::TakeParameterList */
	virtual void InitSupport(ostream& out, AutoArrayT<ElementCardT>& elem_cards, 
		const iArrayT& surface_nodes, int numDOF, int max_node_num, ModelManagerT* model);
	
	/** initialization of meshless information. This method must be called once after 
	 * a call to MeshFreeElementSupportT::TakeParameterList */
	virtual void InitSupport(ostream& out, AutoArrayT<ElementCardT>& elem_cards_displ, 
		AutoArrayT<ElementCardT>& elem_cards_plast, const iArrayT& surface_nodes, 
		int numDOF_displ, int numDOF_plast, int max_node_num, ModelManagerT* model);
	
	/** \name construct nodal field */
	/*@{*/
	void SetNodalField(const dArray2DT& dof);
	void GetNodalField(const dArray2DT& dof, const iArrayT& nodes,
		dArray2DT& field) const;
	void FreeNodalField(void); // free memory associated with reconstruction
	/*@}*/

	/** write data for any cell containing the specified node as
	 * well as the nodes own neighborhood. (map == NULL) means no map. */
	void TraceNode(ostream& out, int node, const ElementBaseT& element_group);

	/** weight meshfree computational cost at nodes */
	void WeightNodes(iArrayT& weight) const;

private:

	/** called during MeshFreeElementSupportT::InitSupport */
	void CollectNodesData(ostream& out, int max_node_num, ModelManagerT* model);

	/* collect nodes with interpolating shape functions */
	void SetAllFENodes(const iArrayT& fe_nodes);

protected:

	/* mesh-free shape functions */
	MeshFreeShapeFunctionT* fMFShapes;
	MeshFreeShapeFunctionT* fMFShapes_displ;
	MeshFreeShapeFunctionT* fMFShapes_plast;
	MeshFreeNodalShapeFunctionT* fNodalShapes;

	/* manager for dynamic local arrays */
	LocalArrayGroupT fLocGroup;
	
	LocalArrayT fLocField_1;
	LocalArrayT fLocField_2;
	
	/* variable length element arrays */
	int fNumElemenNodes;
	int fNumElemenNodes_displ;
	int fNumElemenNodes_plast;
	const RaggedArray2DT<int>* fElemNodesEX;
	const RaggedArray2DT<int>* fElemNodesEX_displ;
	const RaggedArray2DT<int>* fElemNodesEX_plast;
	RaggedArray2DT<int>        fElemEqnosEX;
	RaggedArray2DT<int>        fElemEqnosEX_displ;
	RaggedArray2DT<int>        fElemEqnosEX_plast;
	ArrayT<iArrayT> fUNodeLists; // pointers to fElemNodesEX data
	ArrayT<iArrayT> fUNodeLists_displ; // pointers to fElemNodesEX_displ data
	ArrayT<iArrayT> fUNodeLists_plast; // pointers to fElemNodesEX_plast data
	
	/* variable length workspace managers */
	nArrayGroupT<double>  fNEEArray;
	nMatrixGroupT<double> fNEEMatrix;
	nMatrixGroupT<double> fNEEMatrixOD;

	/* mesh-free data */
	iArrayT fFENodes;   // interpolant nodes
	iArrayT fEFGNodes;  // pure EFG nodes
	iArrayT fAllFENodes;	
	iArrayT fOffGridNodes;
	
	/* nodal field reconstruction */
	bool      fFieldSet;
	dArray2DT fNodalU;
	iArrayT   fGlobalToNodesUsedMap;
	int       fMapShift;

private:

	/** \name pointers to lists of class parameters used during initialization */
	/*@{*/
	const ParameterListT* fOffGridID;
	const ParameterListT* fInterpolantID;
	const ParameterListT* fMeshlessID;
	/*@}*/
};

inline void MFGPElementSupportT::SetShape(MeshFreeShapeFunctionT* mf_shapes) { fMFShapes = mf_shapes; }
inline void MFGPElementSupportT::SetShape(MeshFreeShapeFunctionT* mf_shapes_displ, 
			MeshFreeShapeFunctionT* mf_shapes_plast) { 
	fMFShapes_displ = mf_shapes_displ;
	fMFShapes_plast = mf_shapes_plast; 
}
inline int MFGPElementSupportT::NumElementNodes(void) const { return fNumElemenNodes; }
inline void MFGPElementSupportT::Register(LocalArrayT& u) { fLocGroup.Register(u); }
inline void MFGPElementSupportT::Register(dArrayT& v) { fNEEArray.Register(v); }
inline void MFGPElementSupportT::Register(dMatrixT& m) { fNEEMatrix.Register(m); }

/* element nodes */
inline const RaggedArray2DT<int>& MFGPElementSupportT::ElementNodes(void) {
	if (!fElemNodesEX) 
		ExceptionT::GeneralFail("MFGPElementSupportT::ElementNodes", "pointer not set");
	return *fElemNodesEX;
	
	if (!fElemNodesEX_displ) 
		ExceptionT::GeneralFail("MFGPElementSupportT::ElementNodes", "pointer not set");
	return *fElemNodesEX_displ;
	
	if (!fElemNodesEX_plast) 
		ExceptionT::GeneralFail("MFGPElementSupportT::ElementNodes", "pointer not set");
	return *fElemNodesEX_plast;
}

} /* namespace Tahoe */


#endif /* _MFGP_SUPPORT_T_H_ */
