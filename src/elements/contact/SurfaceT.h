/* $Id: SurfaceT.h,v 1.23 2007-10-09 23:24:47 rjones Exp $ */
#ifndef _SURFACE_T_H_
#define _SURFACE_T_H_

/* direct members */
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "ArrayT.h"
#include "RaggedArray2DT.h"
#include "FaceT.h"
#include "GeometryT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class ElementSupportT;

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
		(const ElementSupportT& support, StringT ss_ID, ostream& out);
//		(const ElementSupportT& support, ArrayT<StringT>& ss_ID, ostream& out);

	/* compute neighbors and initalize coordinates */
	void Initialize (const ElementSupportT& support);

	/* update kinetimatics */
	void UpdateConfiguration();

	inline void SetTag(int tag) {fTag = tag;}

	/* access functions */
	inline const int Tag(void) const {return fTag;}
	inline const int NumNodes(void) const {return fGlobalNodes.Length();}
	inline const int NumFaces(void) const {return fFaces.Length();}
	inline const int NumSD(void) const {return fNumSD;}
	inline const int NodeNumber(int local_num) const {return fGlobalNodes[local_num];}
	inline const iArrayT&   GlobalNodes(void) const {return fGlobalNodes;}
	// NOTE this copy of fGlobalNodes
	inline const iArray2DT& GlobalNodeNumbers(void) const 
		{return fGlobalNodeNumbers;}
	inline const dArray2DT& Coordinates(void) const {return fCoordinates;}
	inline const ArrayT<FaceT*>& Faces(void)  const  {return fFaces;}
//inline ArrayT<FaceT*>& NeighborFaces(void) {return fFaces;}
	inline const double* Position(int i) const {return fCoordinates(i);}
	inline const double* Normal(int i) const   {return fNormals(i);}
	inline const double* Tangent1(int i) const {return fTangent1s(i);}
	inline const double* Tangent2(int i) const {return fTangent2s(i);}
	/* these are predicated on the surfaces being homogeneous */
	inline int NumNodesPerFace(void) const
		{return fNumNodesPerFace;}
	// NOTE if a surface has a single geometry type, why have FaceT's?
	inline GeometryT::CodeT GeometryType(void) const
		{return fGeometryType;}
//inline int NumIPs(void) const
//{return fFaces[0]->NumIPs();}


  protected:
	int fNumNodesPerFace;

	GeometryT::CodeT fGeometryType;

	/* surface specification modes */
   	enum SurfaceSpecModeT {kNodesOnFacet = 0,
                               kSideSets = 1,
                           kBodyBoundary = 2};

	/* FE boundary faces, pointers to base class*/
	ArrayT<FaceT*> fFaces ;

	/* list of global node numbers i.e local->global map */
	iArrayT fGlobalNodes;
	iArray2DT fGlobalNodeNumbers;

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
	const ElementSupportT* fElementSupport;

  private:
	void ComputeNeighbors(void);
	void ComputeSurfaceBasis(void);



};

/* inlines */


} // namespace Tahoe 
#endif /* _SURFACE_T_H_ */
