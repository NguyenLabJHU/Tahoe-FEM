/* $Id: CSEBaseT.h,v 1.4.8.1 2002-04-28 22:26:21 paklein Exp $ */
/* created: paklein (11/19/1997) */

#ifndef _CSE_BASE_T_H_
#define _CSE_BASE_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "iArrayT.h"
#include "LocalArrayT.h"
#include "GeometryT.h"
#include "fstreamT.h"

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
	         NoRotateAnisotropic = 2};

	/** indicies for nodal output */
	enum NodalOutputCodeT {
	                  NodalCoord = 0, /**< reference coordinates */
                       NodalDisp = 1, /**< displacements */
                   NodalDispJump = 2, /**< opening displacements */
                   NodalTraction = 3, /**< traction */
                    MaterialData = 4  /**< output from constitutive relations */ };

	/** indicies for element output */
	enum ElementOutputCodeT {
	             Centroid = 0, /**< reference coordinates of centroid */
           CohesiveEnergy = 1, /**< dissipated energy */
                 Traction = 2  /**< element-averaged traction */ };

	/* constructor */
	CSEBaseT(const ElementSupportT& support, const FieldT& field);

	/* destructor */
	~CSEBaseT(void);

	/* allocates space and reads connectivity data */
	virtual void Initialize(void);

	/* start of new time sequence */
	virtual void InitialCondition(void);

	/* finalize time increment */
	virtual void CloseStep(void);

	/* resets to the last converged solution */
	virtual void ResetStep(void);

	/* solution calls */
	virtual void AddNodalForce(int node, dArrayT& force);

	/* returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void); //not implemented

	/* writing output */
	virtual void RegisterOutput(void);
	virtual void WriteOutput(IOBaseT::OutputModeT mode);

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);

protected: /* for derived classes only */

	/* element status flags */
	enum StatusT {kOFF = 0,
                   kON = 1,
               kMarked = 2};

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;

	/* nodal value calculations */
	virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
		const iArrayT& e_codes, dArray2DT& e_values) = 0;

	/* construct output labels array */
	virtual void GenerateOutputLabels(const iArrayT& n_codes, ArrayT<StringT>& n_labels,
		const iArrayT& e_codes, ArrayT<StringT>& e_labels) const;

	/* write current element information to the output */
	void CurrElementInfo(ostream& out) const;
	
private:

	/* close surfaces to zero gap */
	void CloseSurfaces(void) const;

	/* return the default number of element nodes */
	virtual int DefaultNumElemNodes(void) const;
	//NOTE: needed because ExodusII does not store ANY information about
	//      empty element groups, which causes trouble for parallel execution
	//      when a partition contains no element from a group.
	
protected:

	/* parameters */
	GeometryT::CodeT fGeometryCode;
	int fNumIntPts;
	int fCloseSurfaces;
	int fOutputArea;

	/* output control */
	int fOutputID;
	iArrayT	fNodalOutputCodes;
	iArrayT	fElementOutputCodes;
	iArrayT fNodesUsed;

	/* output stream for fracture area */
	ofstreamT farea_out;
	double fFractureArea;

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
};

#endif /* _CSE_BASE_T_H_ */
