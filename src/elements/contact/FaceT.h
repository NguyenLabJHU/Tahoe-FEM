/* $Id: FaceT.h,v 1.6 2001-04-16 17:30:51 rjones Exp $ */

#ifndef _FACE_T_H_
#define _FACE_T_H_

/* direct members */
#include "iArrayT.h"
#include "dArray2DT.h"
#include "ArrayT.h"
#include "GeometryT.h"
#if 0
#include <math.h>
#endif


/* forward declarations */
class SurfaceT;
class iArrayT;
class dArrayT;
class dMatrixT;

class FaceT 
{
public:

        /* constructor */
        FaceT	(SurfaceT& surface, 	
		dArray2DT& surface_coordinates,
		int num_face_nodes,
		int* connectivity);


        /* (virtual) destructor */
        virtual ~FaceT(void);

	/* initialization after construction */
	virtual void Initialize(void)=0;

	/* geometric computation */
        virtual void ComputeCentroid(double& centroid) =0; 
	virtual double ComputeRadius(void)=0;
        virtual void ComputeNormal
		(dArrayT& local_coordinates,double& normal)=0; 
	virtual void ComputeShapeFunctions
		(dArrayT& local_coordinates, dArrayT& shape_functions)=0;
	virtual void ComputeShapeFunctions
		(dArrayT& local_coordinates, dMatrixT& shape_functions)=0;
	virtual double ComputeJacobian
		(dArrayT& local_coordinates)=0;
        virtual bool Projection
		(double& point, double& normal, 
		dArrayT& local_coordinates, double gap)=0; 

	/* access functions */
	inline const int NumNodes(void) const 
		{return fConnectivity.Length();}
	inline const GeometryT::CodeT GeometryType(void) const 
		{return fGeometryType;}
 	inline const iArrayT& Connectivity(void) const 
		{return fConnectivity;} 
	inline const int NumVertexNodes(void) 
		{return fNumVertexNodes;}
#if 0
	//2D
	inline const int Left() {return fConnectivity[0];}
	inline const int Right() {return fConnectivity[1];}
	//3D (2d too)
	inline const int Next(int i)
		{return mod(i + 1, fNumVertexNodes);}
	inline const int Prev(int i)
		{return mod(i - 1, fNumVertexNodes);}
	inline const int LocalNodeNumber(int node_num)
		{for (int i = 0; i < fConnectivity.Length(); i++) {
			if (node_num == fConnectivity[i]) return i ; } }
#endif

#if 0
// need Face2DT and Face3DT?
// these functions should know the tolerances...
	/* check functions */   // shouldnot recompute face normal
	inline bool CheckOpposition(double& nm, double& face_nm)
	inline bool CheckOpposition(double* nm)
		{return 0; }
	inline bool CheckLocalCoordinates(double* xi, double tol_xi)
		{return 0; }
	inline bool CheckGap(double gap, double tol_g)
		{ {gap < tol_g} ? {return 1} : {return 0} ;}
#endif
		

protected:
	/* geometry type */
	GeometryT::CodeT fGeometryType;

	/* number of vertex nodes */
	int fNumVertexNodes;

	/* reference to parent surface */
        const SurfaceT& fSurface;

	/* reference to nodal coordinates of surface */
        const dArray2DT& fSurfaceCoordinates;

	/* connectivity, in node numbers local to surface */
	iArrayT fConnectivity;

private:
	/* face neighbors */
//ArrayT<FaceT*> fNeighbors;

};

#endif /* _FACE_T_H_ */

