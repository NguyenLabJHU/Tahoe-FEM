/* $Id: MFGPElementSupportT.h,v 1.1 2005-04-26 22:32:28 kyonten Exp $ */
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
class D3MeshFreeShapeFunctionT;
class D3MeshFreeSupportT;
class ElementCardT;
class StringT;
class ElementBaseT;
class ModelManagerT;

/** support for meshfree calculations */
class MFGPElementSupportT//: public ParameterInterfaceT
{
public:

	/** constructor */
	MFGPElementSupportT(void);

	/** destructor */
	//virtual ~MFGPElementSupportT(void) { };
	
	/** accessors */
	D3MeshFreeSupportT& D3MeshFreeSupport(void) const;

	/** set pointer to the shape functions 
	  * assume that first and second parameters of the SetShape are
	  * displacement and plastic multiplier shape functions respectively
	  * only one call to this SetShape is needed */
	void SetShape(D3MeshFreeShapeFunctionT* mf_shapes_displ, D3MeshFreeShapeFunctionT* mf_shapes_plast);

	/** \name memory managers for variable numbers of element nodes */
	/*@{*/
	/** redimension the off-diagonal, non-square matrices */
	void MFGPElementSupportT::SetOffDiagMatrix(int element);

	/** register matrix that need non-square dimensions of element equations */
	// if field_1 and field_2 have same dof, use MeshFreeElementSupportT::Register
	void Register(dMatrixT& m, LocalArrayT& field_1, LocalArrayT& field_2);
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	//virtual void DefineSubs(SubListT& sub_list) const;

	/** accept parameter list */
	//virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** initialization of meshless information. This method must be called once after 
	 * a call to MFGPElementSupportT::TakeParameterList */
	virtual void InitSupport(ostream& out, AutoArrayT<ElementCardT>& elem_cards_displ, 
		AutoArrayT<ElementCardT>& elem_cards_plast, const iArrayT& surface_nodes, int numDOF_displ, 
		int numDOF_plast, int max_node_num, ModelManagerT* model);

protected:

	/* mesh-free shape functions */
	D3MeshFreeShapeFunctionT* fShapes_displ;
	D3MeshFreeShapeFunctionT* fShapes_plast;

	/* local arrays */
	LocalArrayT fLocField_1;
	LocalArrayT fLocField_2;
	
	/* variable length element arrays */
	int fNumElemNodesDispl;
	int fNumElemNodesPlast;
	const RaggedArray2DT<int>* fElemNodesEXDispl;
	RaggedArray2DT<int>        fElemEqnosEXDispl;
	const RaggedArray2DT<int>* fElemNodesEXPlast;
	RaggedArray2DT<int>        fElemEqnosEXPlast;
	ArrayT<iArrayT> fUNodeLists; // pointers to fElemNodesEX data
	
	/* variable length workspace managers */
	nMatrixGroupT<double> fNEEMatrix;  // Kulambda or Klambdau  
};

inline void MFGPElementSupportT::SetShape(D3MeshFreeShapeFunctionT* mf_shapes_displ, 
			D3MeshFreeShapeFunctionT* mf_shapes_plast) { 
	fShapes_displ = mf_shapes_displ;
	fShapes_plast = mf_shapes_plast; 
}


} /* namespace Tahoe */


#endif /* _MFGP_SUPPORT_T_H_ */
