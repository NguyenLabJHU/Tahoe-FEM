/* $Id */

#ifndef _FACE_T_H_
#define _FACE_T_H_

/* direct members */
#include "iArray2DT.h"

/* forward declarations */
class VectorT;

/* derived from SurfaceT ? */
class FaceT : public SurfaceT
{
public:

        /* constructor */
        FaceT(SurfaceT& surface,iArrayT& connectivity, dArrayT& coordinates);

        /* destructor */
        virtual ~FaceT(void);

        virtual void ComputeCentroid(Vector& centroid) =0; 
	virtual double ComputeRadius(void)=0;
        virtual void ComputeNormal(Vector& normal)=0; 
	virtual void ComputeShapeFunctions
		(double& local_coordinates, double& shape_functions)=0;
	virtual void ComputeShapeFunctions
		(double& local_coordinates, MatrixT shape_functions)=0;
	virtual void ComputeJacobians
		(double& local_coordinates, double& jacobians)=0;
        virtual bool Projection
		(Vector& point, Vector& normal, 
		double& local_coordinates, double gap)=0; 
	virtual int NumNodes() {return fNumNodes;}
	/* access functions */
 	inline const iArrayT& Connectivity() { return fConnectivity:} const;
protected:

private:
	/* number of nodes */
	int fNumNodes;
	int fNumVertexNodes;
	int fNumEdgeNodes;
	int fNumInteriorNodes;

	/* reference to parent surface */
        const SurfaceT& fSurface;

	/* list of node numbers local to the surface */
	/* these are a CCW ring */
	iArrayT fConnectivity;

	/* reference to nodal coordinates of surface */
        const dArray2DT& fCoordinates;

	/* face neighbors */
        ArrayT<FaceT*> fNeighbors;

	


};

#endif /* _FACE_T_H_ */

