/* $Id: FaceT.h,v 1.2 2001-04-09 22:28:55 rjones Exp $ */

#ifndef _FACE_T_H_
#define _FACE_T_H_

/* direct members */
#include "iArrayT.h"
#include "dArray2DT.h"
#include "ArrayT.h"
#include "GeometryT.h"


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
	inline int NumNodes(void) const 
		{return fNumNodes;}
	inline GeometryT::CodeT GeometryType(void) const 
		{return fGeometryType;}
 	inline const iArrayT& Connectivity(void) const 
		{return fConnectivity;} 

protected:

	/* number of nodes */
	int fNumNodes;

	/* geometry type */
	GeometryT::CodeT fGeometryType;

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

