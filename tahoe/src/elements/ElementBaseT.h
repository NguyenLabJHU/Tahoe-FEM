/* $Id: ElementBaseT.h,v 1.1.1.1 2001-01-29 08:20:34 paklein Exp $ */
/* created: paklein (05/24/1996)                                          */

#ifndef _ELEMENTBASE_T_H_
#define _ELEMENTBASE_T_H_

#include "GlobalT.h"

/* base class */
#include "iConsoleObjectT.h"

/* direct members */
#include "iArray2DT.h"
#include "ElementMatrixT.h"
#include "dMatrixT.h"
#include "ElementCardT.h"
#include "dArrayT.h"
#include "AutoArrayT.h"
#include "IOBaseT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class FEManagerT;
class NodeManagerT;
class LocalArrayT;
class LoadTime;
class eControllerT;
template <class TYPE> class RaggedArray2DT;
class iAutoArrayT;
class dArray2DT;
class StringT;

class ElementBaseT: public iConsoleObjectT
{
public:

	/* constructor */
	ElementBaseT(FEManagerT& fe_manager);

	/* destructor */
	virtual ~ElementBaseT(void);

	/* accessors */
	const FEManagerT& FEManager(void) const;
	const GlobalT::StateT& RunState(void) const;
	int NumSD(void) const;
	int NumDOF(void) const;

	/* set the controller */
	virtual void SetController(eControllerT* controller);

	/* allocates space and reads connectivity data */
	virtual void Initialize(void);

	/* re-initialize: signal to element group that the global
	 * equations numbers are going to be reset so that the group
	 * has the opportunity to reconnect and should reinitialize
	 * an dependencies on global equation numbers obtained from the
	 * NodeManagerT.
	 *
	 * NOTE: any memory allocated after initial construction (or since
	 * the last Reinitialize) should be "shuffled down" at this point, ie.
	 * reallocated and copied, to make room for the global stiffness
	 * matrix */
	virtual void Reinitialize(void);

	/* form of tangent matrix - symmetric by default */
	virtual GlobalT::SystemTypeT TangentType(void) const = 0;
		
	/* solution calls */
	void FormLHS(void);
	void FormRHS(void);
	virtual void AddNodalForce(int node, dArrayT& force) = 0;

	/* initialize/finalize time increment */
	virtual void InitStep(void);
	virtual void CloseStep(void);
	virtual void ResetStep(void); // restore last converged state

	/* element level reconfiguration for the current solution */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/* returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void) = 0;

	/* append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
	
	/* writing output */
	virtual void RegisterOutput(void) = 0;
	virtual void WriteOutput(IOBaseT::OutputModeT mode) = 0;

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode) = 0;
	
	/* appends group connectivities to the array (X -> geometry, U -> field) */
	virtual void ConnectsX(AutoArrayT<const iArray2DT*>& connects) const;
	virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	             AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;
	
	/* returns a pointer to the specified LoadTime function */
	LoadTime* GetLTfPtr(int num) const;
	
	/* initial condition/restart functions (per time sequence)
	 *
	 * Set to initial conditions.  The restart functions
	 * should read/write any data that overrides the default
	 * values */
	virtual void InitialCondition(void);
	virtual void ReadRestart(istream& in);
	virtual void WriteRestart(ostream& out) const;

	/* element card data */
	int NumElements(void) const;
	int CurrElementNumber(void) const;
	ElementCardT& ElementCard(int card) const;
	ElementCardT& CurrentElement(void) const;

	/* returns 1 if DOF's are interpolants of the nodal values */
	virtual int InterpolantDOFs(void) const;
	virtual void NodalDOFs(const iArrayT& nodes, dArray2DT& DOFs) const;

	/* block ID for the specified element */
	int ElementBlockID(int element) const;

	/* weight the computational effort of every node */
	virtual void WeightNodalCost(iArrayT& weight) const;

protected: /* for derived classes only */

	/* block info indexes */
	enum BlockIndexT {kID = 0,
	            kStartNum = 1,
	            kBlockDim = 2,
	            kBlockMat = 3,
	       kBlockDataSize = 4}; // number of items in the block

	/* get local element data, X for geometry, U for
	 * field variables */
	const LocalArrayT& SetLocalX(LocalArrayT& localarray); // the geometry
	const LocalArrayT& SetLocalU(LocalArrayT& localarray); // the degrees of freedom

		 	
	/* called by FormRHS and FormLHS */
	virtual void LHSDriver(void) = 0;
	virtual void RHSDriver(void) = 0;

	/* assembling the left and right hand sides */
	void AssembleRHS(void) const;
	void AssembleLHS(void) const;
	
	/* element loop operations */
	void Top(void);
	virtual bool NextElement(void);

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
	
	/* element data */
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);
	virtual void ReadConnectivity(ifstreamT& in, ostream& out);
	virtual void WriteConnectivity(ostream& out) const;

	/* generate connectivities with local numbering -
* returns the number of nodes used by the element group */
	int MakeLocalConnects(iArray2DT& localconnects);
	void NodesUsed(iArrayT& nodes_used) const;

	/* return pointer to block data given the ID */
	const int* BlockData(int block_ID) const;

	/* write all current element information to the stream */
	virtual void CurrElementInfo(ostream& out) const;

private:

	/* set element cards array */
	void SetElementCards(void);

	/* I/O formats */
	void ReadConnectivity_ASCII(ifstreamT& in, ostream& out, int num_blocks);
	void ReadConnectivity_TahoeII(ifstreamT& in, ostream& out, int num_blocks);
	void ReadConnectivity_ExodusII(ifstreamT& in, ostream& out, int num_blocks);
	void WriteConnectivity_ASCII(ostream& out) const;

	/* return the default number of element nodes */
	virtual int DefaultNumElemNodes(void) const;
	//NOTE: needed because ExodusII does not store ANY information about
	//      empty element groups, which causes trouble for parallel execution
	//      when a partition contains no element from a group.
	
protected:

	FEManagerT&   fFEManager;
	NodeManagerT* fNodes;

	/* element controller */
	eControllerT* fController;

	/* derived data */
	int	fNumSD;
	int	fNumDOF;
	int fNumElemNodes;	
	int	fNumElemEqnos;
	int	fAnalysisCode;
	
	/* element-by-element info */
	int	fNumElements;
	AutoArrayT<ElementCardT> fElementCards;
	
	/* grouped element arrays */
	iArray2DT fConnectivities;		
	iArray2DT fEqnos;			
	
	/* element tangent matrix and force vector */								
	ElementMatrixT fLHS;
	dArrayT        fRHS;
	
	/* data for multiple connectivity blocks */
	iArray2DT fBlockData; // [ID] [1st group element] [size] per block
};

/* inline functions */

/* up */
inline const FEManagerT& ElementBaseT::FEManager(void) const { return fFEManager; }
inline int ElementBaseT::NumSD(void) const { return fNumSD; }
inline int ElementBaseT::NumDOF(void) const { return fNumDOF; }

/* currElement operations */
inline void ElementBaseT::Top(void) { fElementCards.Top(); }
inline bool ElementBaseT::NextElement(void) { return fElementCards.Next(); }

/* element card */
inline int ElementBaseT::NumElements(void) const { return fNumElements; }
inline int ElementBaseT::CurrElementNumber(void) const { return fElementCards.Position(); }
inline ElementCardT& ElementBaseT::CurrentElement(void) const { return fElementCards.Current(); }
inline ElementCardT& ElementBaseT::ElementCard(int card) const { return fElementCards[card]; }

/* called by FormRHS and FormLHS */
inline void ElementBaseT::LHSDriver(void) { }
inline void ElementBaseT::RHSDriver(void) { }

#endif /* _ELEMENTBASE_T_H_ */
