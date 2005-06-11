/* $Id: CSEBaseT.h,v 1.22 2005-03-15 07:15:35 paklein Exp $ */
/* created: paklein (11/19/1997) */
#ifndef _CSE_BASE_T_H_
#define _CSE_BASE_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "iArrayT.h"
#include "LocalArrayT.h"
#include "GeometryT.h"
#include "ofstreamT.h"

namespace Tahoe {

/* forward declarations */
class SurfaceShapeT;
class StringT;

/** base class for cohesive surface elements */
class CSEBaseT: public ElementBaseT
{
public:

	/* flags for derived class support */
	enum FormulationT {Isotropic = 0,
	                 Anisotropic = 1, 
	         NoRotateAnisotropic = 2,
	         	ModeIAnisotropic = 3,
	            RigidAnisotropic = 4,
	       NodalRigidAnisotropic = 5};

	/** indicies for nodal output */
	enum NodalOutputCodeT {
	                  NodalCoord = 0, /**< reference coordinates */
                       NodalDisp = 1, /**< displacements */
                   NodalDispJump = 2, /**< opening displacements */
                   NodalTraction = 3, /**< traction */
                    MaterialData = 4, /**< output from constitutive relations */ 
                    InternalData = 5  /**< 'output' within Tahoe */
	};
	
	/** indicies for element output */
	enum ElementOutputCodeT {
	             Centroid = 0, /**< reference coordinates of centroid */
           CohesiveEnergy = 1, /**< dissipated energy */
                 Traction = 2  /**< element-averaged traction */ };

	/* constructors */
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	CSEBaseT(const ElementSupportT& support);
#else
	CSEBaseT(ElementSupportT& support);
#endif

	/* destructor */
	~CSEBaseT(void);

	/* start of new time sequence */
	virtual void InitialCondition(void);

	/* finalize time increment */
	virtual void CloseStep(void);

	/* resets to the last converged solution */
	virtual GlobalT::RelaxCodeT ResetStep(void);

#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/* solution calls */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);
#endif

	/* returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void); //not implemented

	/* writing output */
	virtual void RegisterOutput(void);
	virtual void WriteOutput(void);

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);

#ifdef _FRACTURE_INTERFACE_LIBRARY_	
	/* Initialize fields passed in from the outside */
	virtual void InitStep(void);
#endif

	/** return true if element type is axisymmetric */
	bool Axisymmetric(void) const { return fAxisymmetric; };

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
protected:

	/** \name construction of connectivities */
	/*@{*/
	/** extract element block info from parameter list to be used. Method is
	 * used in conjunction with ElementBaseT::DefineElements to initialize
	 * the element group connectivities. CSEBaseT::CollectBlockInfo first calls
	 * ElementBaseT::CollectBlockInfo and generates corrected connectivities
	 * for improper higher order elements. */
	virtual void CollectBlockInfo(const ParameterListT& list, ArrayT<StringT>& block_ID,  ArrayT<int>& mat_index) const;
	/*@}*/

	/** \name output data */
	/*@{*/
	virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
		const iArrayT& e_codes, dArray2DT& e_values) = 0;

	/** set output labels array */
	virtual void GenerateOutputLabels(const iArrayT& n_codes, ArrayT<StringT>& n_labels,
		const iArrayT& e_codes, ArrayT<StringT>& e_labels) const;
	/*@}*/

	/* write current element information to the output */
	void CurrElementInfo(ostream& out) const;
	
	/* number of facet nodes as a function of number of element nodes */
	virtual int NumFacetNodes(void) { return NumElementNodes()/2; }
	
	/* return the default number of element nodes */
	virtual int DefaultNumElemNodes(void) const;
	//NOTE: needed because ExodusII does not store ANY information about
	//      empty element groups, which causes trouble for parallel execution
	//      when a partition contains no element from a group.
	
private:

	/* close surfaces to zero gap */
	void CloseSurfaces(void) const;
	
protected:

	/* parameters */
	GeometryT::CodeT fGeometryCode;
	int fNumIntPts;
	bool fAxisymmetric;
	bool fCloseSurfaces;
	bool fOutputArea;

	/** \name output control */
	/*@{*/
	int fOutputID;
	iArrayT	fNodalOutputCodes;
	iArrayT	fElementOutputCodes;
	iArrayT fNodesUsed;
	bool fOutputGlobalTractions;
	/*@}*/

	/** \name output of fracture area */
	/*@{*/
	ofstreamT farea_out;
	double fFractureArea;
	/*@}*/

	/* shape functions */
	SurfaceShapeT* fShapes;

	/* local arrays */
	LocalArrayT	fLocInitCoords1; // ref geometry of 1st facet
	LocalArrayT	fLocCurrCoords;  // current geometry
	iArrayT		fNodes1;         // nodes on 1st facet

	/* work space */
	dArrayT  fNEEvec;	
	dMatrixT fNEEmat;

	/* parameters */
	static const int NumNodalOutputCodes;
	static const int NumElementOutputCodes;
	
	/** output connectivities. For low-order element types, these will
	 * be the same as in ElementBaseT::fBlockData. These will be
	 * different if the connectivities needed for the element calculations
	 * is not compatible with the element topologies supported by
	 * most database types or post-processors */
	ArrayT<StringT> fOutputBlockID;	 
};

} // namespace Tahoe 
#endif /* _CSE_BASE_T_H_ */
