/* $Id: MeshFreeElementSupportT.h,v 1.12 2005-11-18 06:31:25 paklein Exp $ */
/* created: paklein (11/12/1999) */
#ifndef _MFREE_SUPPORT_T_H_
#define _MFREE_SUPPORT_T_H_

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
class MeshFreeElementSupportT: public ParameterInterfaceT
{
public:

	/** constructor */
	MeshFreeElementSupportT(void);

	/** destructor */
	virtual ~MeshFreeElementSupportT(void) { };
	
	/** accessors */
	MeshFreeSupportT& MeshFreeSupport(void) const;

	/** set pointer to the shape functions */
	void SetShape(MeshFreeShapeFunctionT* mf_shapes);

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

	/** register local arrays */
	void Register(LocalArrayT& u);

	/** register vectors that need to be length number of element equations */
	void Register(dArrayT& v);

	/** register matrix that need square dimensions of element equations */
	void Register(dMatrixT& m);
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
	virtual void InitSupport(const ParameterListT& params, ostream& out,
		AutoArrayT<ElementCardT>& elem_cards, const iArrayT& surface_nodes, 
		int numDOF, int max_node_num, ModelManagerT* model);

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
	MeshFreeNodalShapeFunctionT* fNodalShapes;

	/* manager for dynamic local arrays */
	LocalArrayGroupT fLocGroup;
	
	/* variable length element arrays */
	int fNumElemenNodes;
	const RaggedArray2DT<int>* fElemNodesEX;
	RaggedArray2DT<int>        fElemEqnosEX;
	ArrayT<iArrayT> fUNodeLists; // pointers to fElemNodesEX data
	
	/* variable length workspace managers */
	nArrayGroupT<double>  fNEEArray;
	nMatrixGroupT<double> fNEEMatrix;

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

inline void MeshFreeElementSupportT::SetShape(MeshFreeShapeFunctionT* mf_shapes) { fMFShapes = mf_shapes; }
inline int MeshFreeElementSupportT::NumElementNodes(void) const { return fNumElemenNodes; }
inline void MeshFreeElementSupportT::Register(LocalArrayT& u) { fLocGroup.Register(u); }
inline void MeshFreeElementSupportT::Register(dArrayT& v) { fNEEArray.Register(v); }
inline void MeshFreeElementSupportT::Register(dMatrixT& m) { fNEEMatrix.Register(m); }

/* element nodes */
inline const RaggedArray2DT<int>& MeshFreeElementSupportT::ElementNodes(void) {
	if (!fElemNodesEX) 
		ExceptionT::GeneralFail("MeshFreeElementSupportT::ElementNodes", "pointer not set");
	return *fElemNodesEX;
}

} /* namespace Tahoe */


#endif /* _MFREE_SUPPORT_T_H_ */
