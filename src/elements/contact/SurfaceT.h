/* $Id: SurfaceT.h,v 1.3 2001-04-09 22:28:56 rjones Exp $ */

#ifndef _SURFACE_T_H_
#define _SURFACE_T_H_

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class FEManagerT;
class NodeManagerT;
class FaceT;

/* direct members */
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "ArrayT.h"
#include "RaggedArray2DT.h"
#include "GeometryT.h"

class SurfaceT
{
  public:

	/* constructor */
	SurfaceT(void);

	virtual ~SurfaceT(void);

	/* print data */
	void PrintData(ostream& out);

	/* allocate and input face nodes */
	void InputSideSets
		(const FEManagerT& fe_manager, ifstreamT& in, ostream& out);

	/* compute neighbors and initalize coordinates */
	void Initialize (const NodeManagerT* node_manager);

	/* update kinetimatics */
	void UpdateConfiguration();

	/* access functions */
	inline int NumNodes(void) {return fNumNodes;}
	inline int NumFaces(void) {return fNumFaces;}
	inline int NumSD(void) {return fNumSD;}
	inline dArray2DT& Coordinates(void) {return fCoordinates;}
	/* these are predicated on the surfaces being homogeneous */
	inline int NumNodesPerFace(void)
		{return 2;}
//{return fFaces[0]->NumNodes();}
//inline GeometryT::CodeT GeometryType(void)
//{return fFaces[0]->GeometryType();}


  protected:
        /* surface specification modes */
        enum SurfaceSpecModeT {kNodesOnFacet = 0,
                               kSideSets = 1,
                           kBodyBoundary = 2};

	/* FE boundary faces, pointers to base class*/
	ArrayT<FaceT*> fFaces ;

	/* list of global node numbers i.e local->global map */
	iArrayT fNodes;

	/* current surface coordinates */
	dArray2DT fCoordinates;

	/* nodal outward unit normals */
	dArray2DT fNormals; 

	/* neighbors */
	RaggedArray2DT <int>     fNodeNeighbors ; // for averaging, ordered CCW
	RaggedArray2DT <FaceT*>  fFaceNeighbors ; // for contact tracking

  private:
	void ComputeNeighbors(void);
	void ComputeSurfaceNormals(void);

	const NodeManagerT* kNodeManager;

	int fNumNodes;
	int fNumFaces;
	int fNumSD;

};

/* inlines */


#endif /* _SURFACE_T_H_ */
