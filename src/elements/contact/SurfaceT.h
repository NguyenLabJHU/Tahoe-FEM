/* $Id: SurfaceT.h,v 1.11 2001-04-27 00:55:26 rjones Exp $ */

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
	inline const int NumNodes(void) const {return fGlobalNodes.Length();}
	inline const int NumFaces(void) const {return fFaces.Length();}
	inline const int NumSD(void) const {return fNumSD;}
	inline const iArrayT&   GlobalNodes(void) const {return fGlobalNodes;}
	inline const dArray2DT& Coordinates(void) const {return fCoordinates;}
	inline const ArrayT<FaceT*>& Faces(void)  const  {return fFaces;}
//inline ArrayT<FaceT*>& NeighborFaces(void) {return fFaces;}
	inline const double* Position(int i) const {return fCoordinates(i);}
	inline const double* Normal(int i) const   {return fNormals(i);}
	inline const double* Tangent1(int i) const {return fTangent1s(i);}
	inline const double* Tangent2(int i) const {return fTangent2s(i);}
	/* these are predicated on the surfaces being homogeneous */
	inline int NumNodesPerFace(void) const
		{return fFaces[0]->NumNodes();}
	inline GeometryT::CodeT GeometryType(void) const
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
//RaggedArray2DT <FaceT*>  fFaceNeighbors ; // for contact tracking

	int fNumSD;
	int fTag;
	const NodeManagerT* kNodeManager;

  private:
	void ComputeNeighbors(void);
	void ComputeSurfaceBasis(void);



};

/* inlines */


#endif /* _SURFACE_T_H_ */
