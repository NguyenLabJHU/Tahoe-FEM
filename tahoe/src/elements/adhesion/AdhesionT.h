/* $Id: AdhesionT.h,v 1.1.2.1 2002-10-17 04:24:24 paklein Exp $ */
#ifndef _ADHESION_T_H_
#define _ADHESION_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "pArrayT.h"
#include "LocalArrayT.h"
#include "dArray2DT.h"
#include "nVariArray2DT.h"
#include "iGridManagerT.h"
#include "GeometryT.h"

namespace Tahoe {

/* forward declarations */
class SurfaceShapeT;

/** class to calculate surface adhesion forces between bodies */
class AdhesionT: public ElementBaseT
{
public:

	/** constructor */
	AdhesionT(const ElementSupportT& support, const FieldT& field);

	/** destructor */
	~AdhesionT(void);

	/** form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const { return GlobalT::kSymmetric; };

	/** element level reconfiguration for the current solution */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/** initialization after constructor */
	virtual void Initialize(void);

	/** return the force exerted on the specified node */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);

	/** returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void); // not implemented
	
	/** writing output */
	virtual void RegisterOutput(void);
	virtual void WriteOutput(IOBaseT::OutputModeT mode);

	/** compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);  // not implemented

	/** \name connectivities */
	/*@{*/
	virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;
	virtual void ConnectsX(AutoArrayT<const iArray2DT*>& connects) const;
	/*@}*/
	 	
protected:

	/** surface specification modes */
	enum SurfaceSpecModeT {kNodesOnFacet = 0,
                               kSideSets = 1,
                           kBodyBoundary = 2};

	/** print element group data */
	virtual void PrintControlData(ostream& out) const;
	
	/** \name initialization steps */
	/*@{*/
	/** construct the adhesive surfaces. \note This implementation is
	 * adapted from ContactT::EchoConnectivityData */
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);

	/** construct class work space */
	virtual void SetWorkSpace(void);
	/*@}*/

	/** generate face-interaction data - return true if configuration has
	 * changed since the last call */
	bool SetConfiguration(void);

	/* steps in setting configuration */
	virtual bool SetActiveInteractions(void) = 0; // "internal" data
	virtual void SetConnectivities(void) = 0; // "external" data - interface to FEManager

	/** \name surface input methods 
	 * \note All these methods have been adapted from ContactT. */
	/*@{*/
	/** specify facets as lists of nodes */
	void InputNodesOnFacet(ifstreamT& in, GeometryT::CodeT& geom, iArray2DT& facets);

	/** specify facets as side sets */
	void InputSideSets(ifstreamT& in, GeometryT::CodeT& geom, iArray2DT& facets);

	/** specify facets automatically from body boundaries */
	void InputBodyBoundary(ifstreamT& in, ArrayT<GeometryT::CodeT>& geom,
		ArrayT<iArray2DT>& surfaces);
	/*@}*/

	/** return the number of integration points to use for the given face geometry */
	int NumIP(GeometryT::CodeT code) const;

protected:

	/** \name surface data */
	/*@{*/
	/** nodes on the surface faces */
	ArrayT<iArray2DT> fSurfaces;
	
	/** shape functions. Surface shape functions particular to the
	 * topology of faces in each surface. */
	ArrayT<SurfaceShapeT*> fShapes;

	/** current coordinates. Current coordinates in local ordering 
	 * particular to the topology of faces in each surface. */
	ArrayT<LocalArrayT> fLocCurrCoords;
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
	 * Faces at matching indecies of the two arrays are interacting.
	 * Due to the searching scheme the index of the surface for the
	 * face in the AdhesionT::Surface1 array will always be less than or equal
	 * the index of the matching surface in AdhesionT::Surface2. */
	/*@{*/
	AutoArrayT<int> fSurface1;
	AutoArrayT<int> fSurface2;
	/*@}*/
	
	/** search grid */
	iGridManagerT fGrid;

	/** link surfaces in ConnectsU - for graph */
	iArray2DT fSurfaceLinks;
};

} // namespace Tahoe 
#endif /* _ADHESION_T_H_ */
