/* $Id: AdhesionT.h,v 1.9 2004-07-15 08:25:48 paklein Exp $ */
#ifndef _ADHESION_T_H_
#define _ADHESION_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "pArrayT.h"
#include "LocalArrayT.h"
#include "dArray2DT.h"
#include "nVariArray2DT.h"
#include "GeometryT.h"
#include "RaggedArray2DT.h"
#include "nArrayGroupT.h"
#include "nVariMatrixT.h"
#include "nArray2DGroupT.h"
#include "nMatrixGroupT.h"

namespace Tahoe {

/* forward declarations */
class SurfaceShapeT;
class C1FunctionT;
class iGridManagerT;
class ScheduleT;

/** class to calculate surface adhesion forces between bodies */
class AdhesionT: public ElementBaseT
{
public:

	/** constructor */
	AdhesionT(const ElementSupportT& support);

	/** destructor */
	~AdhesionT(void);

	/** form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const { return GlobalT::kSymmetric; };

	/** element level reconfiguration for the current solution */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/** return the force exerted on the specified node */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);

	/** returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void); // not implemented
	
	/** writing output */
	virtual void RegisterOutput(void);
	virtual void WriteOutput(void);

	/** compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);  // not implemented

	/** \name connectivities */
	/*@{*/
	virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;
	virtual void ConnectsX(AutoArrayT<const iArray2DT*>& connects) const;
	/*@}*/

	/** collecting element group equation numbers. See ElementBaseT::Equations
	 * for more documentation. */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** surface specification modes */
	enum SurfaceSpecModeT {kNodesOnFacet = 0,
                               kSideSets = 1,
                           kBodyBoundary = 2};

	/** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
	/*@{*/
	/** form group contribution to the stiffness matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

	/** form group contribution to the residual */
	virtual void RHSDriver(void);
	/*@}*/

	/** \name initialization steps */
	/*@{*/
	/** construct the adhesive surfaces */
	void ExtractSurfaces(const ParameterListT& list);

	/** construct class work space */
	virtual void SetWorkSpace(void);
	/*@}*/

	/** generate face-interaction data - return true if configuration has
	 * changed since the last call */
	bool SetConfiguration(void);

	/** \name surface input methods */
	/*@{*/
	/** specify facets as side sets */
	void InputSideSets(const ParameterListT& list, GeometryT::CodeT& geom, iArray2DT& facets);

	/** specify facets automatically from body boundaries */
	void InputBodyBoundary(const ParameterListT& list, ArrayT<GeometryT::CodeT>& geom,
		ArrayT<iArray2DT>& surfaces);
	/*@}*/

	/** return the number of integration points to use for the given face geometry */
	int NumIP(GeometryT::CodeT code) const;

protected:

	/*@{*/
	/** a flag to allow the penalization of penetration */
	bool fPenalizePenetration;
	bool fAllowSameSurface;
	/*@}*/

	/** \name surface data */
	/*@{*/
	/** nodes on the surface faces */
	ArrayT<iArray2DT> fSurfaces;
	
	/** shape functions over undeformed configuration. Surface shape 
	 * functions particular to the topology of faces in each surface.
	 * Shape functions over the undeformed configuration are needed to
	 * define the integration rule.
	 * \note We could calculate and store the integration parameters
	 *       over the undeformed configuration instead of recalculating
	 *       them over the faces during integration. */
	ArrayT<SurfaceShapeT*> fShapes;

	/** shape functions over current configuration. Surface shape 
	 * functions particular to the topology of faces in each surface. */
	ArrayT<SurfaceShapeT*> fCurrShapes;

	/** initial coordinates. Initial coordinates in local ordering 
	 * particular to the topology of faces in each surface. Initial
	 * coordinates are needed to define the integration rules. */
	ArrayT<LocalArrayT> fLocInitCoords;

	/** current coordinates. Current coordinates in local ordering 
	 * particular to the topology of faces in each surface. */
	ArrayT<LocalArrayT> fLocCurrCoords;

	/** time scaling functions. Allows the interaction associated with
	 * a given surface to scale in time. This allows adaptive solution
	 * methods to introduce the effects of surface interactions gradually.
	 * The interaction of a surface pair is scaled by the average of the
	 * scaling functions for each surface */
	ArrayT<const ScheduleT*> fScaling;

	/** IDs for writing output data. Because each surface may have
	 * its own face geometry, they each need a separate output ID. */
	iArrayT fOutputID;

	/** total force on each face. For each surface, the integrated force
	 * on each face. This is calculated and stored during AdhesionT::RHSDriver. */
	ArrayT<dArray2DT> fFaceForce;

	/** current face area. For each surface, the current area
	 * of each face. This is calculated and stored during AdhesionT::RHSDriver. */
	ArrayT<dArrayT> fCurrentFaceArea;
	/*@}*/
	
	/** \name grouped facet data */
	/*@{*/
	/** centroid of faces. Used for determining interacting faces.
	 * Values are actually average coordinate of the face not the
	 * true centroid. The approximate value is sufficient for determining
	 * the approximate location of the face. */
	dArray2DT fFaceCentroids;

	/** centroid index array. For each face, contains the surface number
	 * and the local index within than surface */
	iArray2DT fFaceIndex;
	
	/** enum for the information in the AdhesionT::fFaceIndex array */
	enum FaceIndexT {
        kSurface = 0, /**< surface containing the face */
     kLocalIndex = 1, /**< local index of the face on the surface */
   kFaceIndexDim = 2  /**< minor dimension of the AdhesionT::fFaceIndex array */
		};
	/*@}*/
	
	/** \name interacting faces
	 * Faces at matching indecies of the two arrays, AdhesionT::Surface1 and
	 * AdhesionT::Surface2, are interacting. Due to the searching scheme the 
	 * index of the surface for the face in the AdhesionT::Surface1 array will 
	 * always be less than or equal the index of the matching surface in 
	 * AdhesionT::Surface2. */
	/*@{*/
	AutoArrayT<int> fSurface1;
	AutoArrayT<int> fSurface2;
	
	/** face-pair connectivies */
	RaggedArray2DT<int> fFaceConnectivities;

	/** face-pair equation numbers */
	RaggedArray2DT<int> fFaceEquations;
	/*@}*/
	
	/** search grid */
	iGridManagerT* fGrid;
	
	/** \name surface interaction */
	/*@{*/
	double fCutOff;
	C1FunctionT* fAdhesion;
	/*@}*/
	
	/** \name work space */
	/*@{*/
	/** matrix of jump gradients (see SurfaceShapeT::Grad_d): [nen] x [nsd] */
	dMatrixT fGrad_d;

	/** integration point coordinates: [nip] x [nsd] */
	dArray2DT fIPCoords2;

	/** normal vectors at the integration points: [nip] x [nsd] */
	dArray2DT fIPNorm2;
	
	/** work space vector: [element DOF] */
	dArrayT fNEE_vec; 

	/** work space matrix: [element DOF] x [element DOF] */
	dMatrixT fNEE_mat;
	/*@}*/

	/** \name dynamic resizing */
	/*@{*/
	/** manager of vectors length number element equations */
	nArrayGroupT<double> fNEE_vec_man;

	/** manager of matricies dimension number element equations */
	nMatrixGroupT<double> fNEE_mat_man;
	
	/** memory manager for AdhesionT::fGrad_d */
	nVariMatrixT<double> fGrad_d_man;
	
	/** manager of data from the face that is in the inner loop
	 * of the integration of quantities over a face pair. */
	nArray2DGroupT<double> fFace2_man;
	/*@}*/

	/** link surfaces in ConnectsU - for graph */
	iArray2DT fSurfaceLinks;
};

} /* namespace Tahoe */
#endif /* _ADHESION_T_H_ */
