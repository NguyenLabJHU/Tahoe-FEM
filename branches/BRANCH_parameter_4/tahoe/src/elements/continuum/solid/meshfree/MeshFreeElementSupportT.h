/* $Id: MeshFreeElementSupportT.h,v 1.10 2004-01-27 01:26:19 cjkimme Exp $ */
/* created: paklein (11/12/1999) */

#ifndef _MFREE_SUPPORT_T_H_
#define _MFREE_SUPPORT_T_H_

/* direct members */
#include "iArrayT.h"
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

class MeshFreeElementSupportT
{
public:

	/** constructor */
	MeshFreeElementSupportT(ifstreamT& in);

	/** destructor */
	virtual ~MeshFreeElementSupportT(void) { };
	
	/* accessors */
	MeshFreeSupportT& MeshFreeSupport(void) const;

	/** return the auto-border flag */
	int AutoBorder(void) const { return fAutoBorder; };

protected:

	 /* write parameters to stream */
	virtual void PrintControlData(ostream& out) const;

	/* initialization */
	virtual void InitSupport(ifstreamT& in, ostream& out,
		AutoArrayT<ElementCardT>& elem_cards, const iArrayT& surface_nodes,
		int numDOF, int max_node_num, ModelManagerT* model);

	/* resize number of field nodes - returns number
	 * of element nodes */
	int SetElementNodes(int element);
	int NumElementNodes(void) const;

	/* construct nodal field */
	void SetNodalField(const dArray2DT& dof);
	void GetNodalField(const dArray2DT& dof, const iArrayT& nodes,
		dArray2DT& field) const;
	void FreeNodalField(void); // free memory associated with reconstruction

	/* mark "dead" cells (no active equations) returns the number of active */
	int MarkActiveCells(AutoArrayT<ElementCardT>& elem_cards);

	/* write data for any cell containing the specified node as
	 * well as the nodes own neighborhood. (map == NULL) means no map. */
	void TraceNode(ostream& out, int node, const ElementBaseT& element_group);

	/* weight nodes */
	void WeightNodes(iArrayT& weight) const;

private:

	/* initialization steps */
	void EchoNodesData(ifstreamT& in, ostream& out, int max_node_num,
		ModelManagerT* model);

	/* collect nodes with interpolating shape functions */
	void SetAllFENodes(const iArrayT& fe_nodes);

protected:

	/* mesh-free parameters */
	int fAutoBorder; // 1 => make all "surface" nodes exact

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
};

inline int MeshFreeElementSupportT::NumElementNodes(void) const { return fNumElemenNodes; }

} // namespace Tahoe 
#endif /* _MFREE_SUPPORT_T_H_ */
