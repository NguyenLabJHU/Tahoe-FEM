/* $Id: SurfaceT.h,v 1.8 2001-04-24 00:33:22 rjones Exp $ */

#ifndef _SURFACE_T_H_
#define _SURFACE_T_H_

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class FEManagerT;
class NodeManagerT;

/* direct members */
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "ArrayT.h"
#include "RaggedArray2DT.h"
#include "FaceT.h"
#include "GeometryT.h"

class SurfaceT
{
  public:

	/* constructor */
	SurfaceT(void);

	virtual ~SurfaceT(void);

	/* print data */
	void PrintConnectivityData(ostream& out);
	void PrintKinematicData(ostream& out);

	/* allocate and input face nodes */
	void InputSideSets
		(const FEManagerT& fe_manager, ifstreamT& in, ostream& out);

	/* compute neighbors and initalize coordinates */
	void Initialize (const NodeManagerT* node_manager);

	/* update kinetimatics */
	void UpdateConfiguration();

	inline void SetTag(int tag) {fTag = tag;}

	/* access functions */
	inline const int Tag(void) const {return fTag;}
	inline const int NumNodes(void) {return fGlobalNodes.Length();}
	inline const int NumFaces(void) {return fFaces.Length();}
	inline const int NumSD(void) {return fNumSD;}
	inline iArrayT&   GlobalNodes(void) {return fGlobalNodes;}
	inline dArray2DT& Coordinates(void) {return fCoordinates;}
	inline ArrayT<FaceT*>& Faces(void) {return fFaces;}
//inline ArrayT<FaceT*>& NeighborFaces(void) {return fFaces;}
	inline const double* Position(int i) {return fCoordinates(i);}
	inline const double* Normal(int i)   {return fNormals(i);}
	inline const double* Tangent1(int i) {return fTangent1s(i);}
	inline const double* Tangent2(int i) {return fTangent2s(i);}

	/* these are predicated on the surfaces being homogeneous */
	inline int NumNodesPerFace(void)
		{return fFaces[0]->NumNodes();}
	inline GeometryT::CodeT GeometryType(void)
		{return fFaces[0]->GeometryType();}


  protected:
        /* surface specification modes */
        enum SurfaceSpecModeT {kNodesOnFacet = 0,
                               kSideSets = 1,
                           kBodyBoundary = 2};

	/* FE boundary faces, pointers to base class*/
	ArrayT<FaceT*> fFaces ;

	/* list of global node numbers i.e local->global map */
	iArrayT fGlobalNodes;

 	/* Nodal data */
	/* current surface coordinates */
	dArray2DT fCoordinates;
	/* nodal outward unit normals and tangents */
	dArray2DT fNormals; 
	dArray2DT fTangent1s; 
	dArray2DT fTangent2s; 

	/* neighbors */
	RaggedArray2DT <FaceT*>  fNodeNeighbors ; // for averaging
	RaggedArray2DT <int>     fLocalNodeInNeighbors ; // for averaging
	RaggedArray2DT <FaceT*>  fFaceNeighbors ; // for contact tracking

  private:
	int fNumSD;
	int fTag;

	void ComputeNeighbors(void);
	void ComputeNeighbors2D(void);
	void ComputeNeighbors3D(void);

	void ComputeSurfaceBasis(void);

	const NodeManagerT* kNodeManager;


};

/* inlines */


#endif /* _SURFACE_T_H_ */
