/* $Id: SurfaceT.h,v 1.2 2001-03-22 22:30:00 rjones Exp $ */


#ifndef _SURFACE_T_H_
#define _SURFACE_T_H_

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class FEManagerT;

/* direct members */
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "dArray2DT.h"


class SurfaceT
{
  public:

	/* constructor */
	SurfaceT(void);

	/* print data */
	void PrintData(ostream& out);

	/* input processing and allocation */
	void InputSideSets
		(const FEManagerT& kFEManager, ifstreamT& in, ostream& out);

	/* compute neighbors */
	void Initialize (void)

	/* update kinematics */
	void UpdateConfiguration(void);

	/* access functions */
	int NumNodes(void) {return fNumNodes;}
	int NumFaces(void) {return fNumFaces;}

  protected:
        /* surface specification modes */
        enum SurfaceSpecModeT {kNodesOnFacet = 0,
                               kSideSets = 1,
                           kBodyBoundary = 2};

	/* connectivity list in local node numbers*/
	ArrayT <FaceT*> fFaces ;

	/* list of global node numbers i.e local->global map */
	iArrayT fNodes;

	/* neighbors, linked list ?? */
	// RAGGED ARRAYS?
	ArrayT <iArrayT> fNodeNeighbors ; // for averaging, ordered CCW
	ArrayT <FaceT*>  fFaceNeighbors ; // for contact tracking

	/* nodal jacobians ?? */ 
	//dArrayT fJacobians;

	/* nodal outward unit normals */
	dArray2DT fNormals; 

	/* current coordinates */
	dArray2DT fCoordinates;

	int fNumSD;
  private:
	int fNumNodes;
	int fNumFaces;

};

#endif /* _SURFACE_T_H_ */
