/* $Id: MeshFreeSupport3DT.cpp,v 1.9 2004-01-27 01:21:11 cjkimme Exp $ */
/* created: paklein (09/13/1998) */

#include "MeshFreeSupport3DT.h"

#include <math.h>
#include <string.h>

#include "ExceptionT.h"
#include "toolboxConstants.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "Vector3T.h"

/* constructor */

using namespace Tahoe;

MeshFreeSupport3DT::MeshFreeSupport3DT(const ParentDomainT* domain, const dArray2DT& coords,
	const iArray2DT& connects, const iArrayT& nongridnodes, ifstreamT& in):
	MeshFreeSupportT(domain, coords, connects, nongridnodes, in)
{

}

/* cutting facet functions */
void MeshFreeSupport3DT::SetCuttingFacets(const dArray2DT& facet_coords,
	int num_facet_nodes)
{
	/* inherited */
	MeshFreeSupportT::SetCuttingFacets(facet_coords, num_facet_nodes);

	if (fNumFacetNodes != 0 &&
	    fNumFacetNodes != 3 && 
	    fNumFacetNodes != 4)
	{
		cout << "\n MeshFreeSupport3DT::SetCuttingFacets: 3D cutting facets must\n"
		     <<   "     have 3 or 4 nodes: " << fNumFacetNodes << endl;
		throw ExceptionT::kSizeMismatch;
	}
	if (fNumFacetNodes == 0 && facet_coords.MajorDim() != 0) {
		cout << "\n MeshFreeSupport3DT::SetCuttingFacets: found facets nodes = 0 with\n" 
		     <<   "     non-zero number of facets (" << facet_coords.MajorDim()
		     << ")" << endl;
		throw ExceptionT::kSizeMismatch;	
	}
}	

/*************************************************************************
* Private
*************************************************************************/

/* process boundaries - nodes marked as "inactive" at the
* current x_node by setting nodal_params = -1.0 */
void MeshFreeSupport3DT::ProcessBoundaries(const dArray2DT& coords,
	const dArrayT& x_node, dArray2DT& nodal_params)
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (coords.MajorDim() != nodal_params.MajorDim()) throw ExceptionT::kSizeMismatch;
	if (coords.MinorDim() != x_node.Length()) throw ExceptionT::kSizeMismatch;
#endif

	/* quick exit if there are no facets */
	if (!fCutCoords || fCutCoords->MajorDim() == 0) return;

	const double* x = x_node.Pointer();
	for (int i = 0; i < coords.MajorDim(); i++)
		if (!Visible(x, coords(i)))
			nodal_params.SetRow(i, -1.0);
}	

/* returns 1 if the path x1-x2 is visible */
int MeshFreeSupport3DT::Visible(const double* x1, const double* x2)
{
	/* quick exit */
	if (!fCutCoords || fCutCoords->MajorDim() == 0) return 1;	

	if (fNumFacetNodes == 3)
	{
		/* exhaustive search for now */
		for (int j = 0; j < fCutCoords->MajorDim(); j++)
		{
			const double* v0 = (*fCutCoords)(j);
			const double* v1 = v0 + 3;
			const double* v2 = v1 + 3;
			if (Intersect(v0, v1, v2, x1, x2)) return 0;
		}
	}
	else if (fNumFacetNodes == 4)
	{
		/* exhaustive search for now */
		for (int j = 0; j < fCutCoords->MajorDim(); j++)
		{
			const double* v0 = (*fCutCoords)(j);
			const double* v1 = v0 + 3;
			const double* v2 = v1 + 3;
			const double* v3 = v2 + 3;
			if (Intersect(v0, v1, v2, x1, x2) ||
			    Intersect(v0, v2, v3, x1, x2)) return 0;
		}
	}
	else
	{
		cout << "\n MeshFreeSupport3DT::Visible: unsupported number of\n"
		     <<   "     facet nodes: " << fNumFacetNodes << endl;
		throw ExceptionT::kOutOfRange;
	}
	return 1;
}

int MeshFreeSupport3DT::Intersect(const double* x0, const double* x1,
	const double* x2, const double *vA, const double* vB) const
{
	/* facet centroid */
	Vector3T<double> xc;
	xc.Average(x0, x1, x2);
	
	/* "expand" facet to insure overlap */
//	const double eps = 0.01;
	const double eps = 0.02;
	double eps1 = 1.0 + eps;

	Vector3T<double> v0, v1, v2;
	v0.Combine(eps1, x0, -eps, xc);
	v1.Combine(eps1, x1, -eps, xc);
	v2.Combine(eps1, x2, -eps, xc);

	/* workspace vectors */
	Vector3T<double> v01, v12, vAB, N;

	/* edge vectors */
	v01.Diff(v1, v0);
	v12.Diff(v2, v1);

	/* facet normal (direction) = v01 x v12 */
	N.Cross(v01, v12);
	
	/* connecting vector */
	vAB.Diff(vB, vA);

	/* intersection distance */
	double den = Vector3T<double>::Dot(N, vAB);
	if (fabs(den) < kSmall) return 0;
	double s = (Vector3T<double>::Dot(N,v0) - Vector3T<double>::Dot(N,vA))/den;
	
	/* AB doesn't cross facet plane */
	if (s < 0.0 || s > 1.0)
		return 0;
	else
	{
		Vector3T<double> vP, Ni, viP;
		vP.Combine(1.0, vA, s, vAB);
	
		/* edge 1 */
		viP.Diff(vP, v0);
		Ni.Cross(v01, viP);
		if (Vector3T<double>::Dot(N, Ni) < 0.0)
			return 0;
		else
		{		
			/* edge 2 */
			viP.Diff(vP, v1);
			Ni.Cross(v12, viP);
			if (Vector3T<double>::Dot(N, Ni) < 0.0)
				return 0;
			else
			{
				Vector3T<double> v20;
				v20.Diff(v0, v2);

				/* edge 3 */
				viP.Diff(vP, v2);
				Ni.Cross(v20, viP);
				if (Vector3T<double>::Dot(N, Ni) < 0.0)
					return 0;
				else
					return 1;
			}
		}
	}
}

/* TEMP: fixed cutting surface shapes */
void MeshFreeSupport3DT::CutCircle(const dArray2DT& coords, const dArrayT& x_node,
	dArrayT& dmax)
{
	/* circular crack centered at {0,0,zc} with radius r
	 * and crack-face normal in the z-direction */

	/* cutting disc parameters */
	double x_c =  0.0;
	double y_c =-15.0;
	double z_c =  0.0;
	double r   = 15.0;
	
	const double* pa = x_node.Pointer();
	for (int i = 0; i < coords.MajorDim(); i++)
	{
		const double* pb = coords(i);
	
		/* must be on opposite size of z = zc */
		if ((pa[2] - z_c)*(pb[2] - z_c) < kSmall)
		{
			double v_ab[] = {pb[0] - pa[0],
			                 pb[1] - pa[1],
			                 pb[2] - pa[2]};

			double l_ab = sqrt(v_ab[0]*v_ab[0] +
			                   v_ab[1]*v_ab[1] +
			                   v_ab[2]*v_ab[2]);
			
			double z_ab = v_ab[2];
			double z_a  = z_c - pa[2];
			double s_0  = fabs(z_a/z_ab);

			/* intersecion v_ab and the z_c plane */
			double  v_a0[] = {pa[0] + s_0*v_ab[0] - x_c,
			                  pa[1] + s_0*v_ab[1] - y_c,
			                  pa[2] + s_0*v_ab[2]};
			
			/* check */
			if (fabs(z_c - v_a0[2]) > kSmall)
			{
				cout << "\n MeshFreeSupport3DT::CutCircle: project error" << endl;
				throw ExceptionT::kGeneralFail;
			}
	
			if ((v_a0[0]*v_a0[0] + v_a0[1]*v_a0[1]) < r*r) dmax[i] = -1.0;
		}
	}
}	
	
void MeshFreeSupport3DT::CutEllipse(const dArray2DT& coords, const dArrayT& x_node,
	dArrayT& dmax)
{
	/* elliptical crack centered at {x_c, y_c, z_c} with radius r
	 * and crack-face normal in the z-direction */

	/* cutting ellipse parameters */
	double x_c =  0.0;
	double y_c =-10.0;
	double z_c =  0.0;
	double a_x = 10.0; // x-direction axis
	double a_y =  7.5; // y-direction axis
	
	const double* pa = x_node.Pointer();
	for (int i = 0; i < coords.MajorDim(); i++)
	{
		const double* pb = coords(i);
	
		/* must be on opposite size of z = zc */
		if ((pa[2] - z_c)*(pb[2] - z_c) < kSmall)
		{
			double v_ab[] = {pb[0] - pa[0],
			                 pb[1] - pa[1],
			                 pb[2] - pa[2]};

			double l_ab = sqrt(v_ab[0]*v_ab[0] +
			                   v_ab[1]*v_ab[1] +
			                   v_ab[2]*v_ab[2]);
			
			double z_ab = v_ab[2];
			double z_a  = z_c - pa[2];
			double s_0  = fabs(z_a/z_ab);

			/* intersecion v_ab and the z_c plane */
			double  v_a0[] = {pa[0] + s_0*v_ab[0] - x_c,
			                  pa[1] + s_0*v_ab[1] - y_c,
			                  pa[2] + s_0*v_ab[2]};
			
			/* check */
			if (fabs(z_c - v_a0[2]) > kSmall)
			{
				cout << "\n MeshFreeSupport3DT::CutEllipse: project error" << endl;
				throw ExceptionT::kGeneralFail;
			}
	
			if ((v_a0[0]*v_a0[0]/(a_x*a_x) +
			     v_a0[1]*v_a0[1]/(a_y*a_y)) < 1.0) dmax[i] = -1.0;
		}
	}
}	
	
