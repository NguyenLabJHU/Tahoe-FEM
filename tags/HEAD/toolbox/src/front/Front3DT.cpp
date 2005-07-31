/* $Id: Front3DT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (03/18/1999)                                          */

#include "Front3DT.h"

#include "FrontSegmentT.h"
#include "FrontNodeT.h"
#include "iArrayT.h"

/* constants */
const double Pi = acos(-1.0);

/* vector functions */
static inline void Average(const double* A, const double* B, double* avg)
{
	avg[0] = 0.5*(A[0] + B[0]);
	avg[1] = 0.5*(A[1] + B[1]);
	avg[2] = 0.5*(A[2] + B[2]);
};

inline static void CrossProduct(const double* A, const double* B, double* AxB)
{   AxB[0] = A[1]*B[2] - A[2]*B[1];
	AxB[1] = A[2]*B[0] - A[0]*B[2];
	AxB[2] = A[0]*B[1] - A[1]*B[0];
};

inline static double Dot(const double* A, const double* B)
{ return A[0]*B[0] + A[1]*B[1] + A[2]*B[2]; };

/* constructor */
Front3DT::Front3DT(double cone, double da, double da_s, int num_pts):
	FrontT(3, 3, cone, da, da_s, num_pts)
{

}

/* destructor */
Front3DT::~Front3DT(void)
{
	/* free front segments data */
	for (int i = 0; i < fFrontLines.Length(); i++)
		delete fFrontLines[i];

	/* empty list */
	fFrontLines.Allocate(0);		
}

/* construct initial front */
void Front3DT::Initialize(const dArray2DT& facet_coords, const iArrayT& fr_facets,
		const iArrayT& fr_edges)
{
	/* checks */
	if (fr_facets.Length() != fr_edges.Length()) throw eSizeMismatch;

	/* facet coords */
	const double* pA;
	const double* pB;
	const double* pC;
	FrontSegmentT* last_seg = NULL;
	for (int i = 0; i < fr_facets.Length(); i++)
	{
		/* order facets coords */
		SetEdgePointers(facet_coords(fr_facets[i]), fr_edges[i], &pA, &pB, &pC);

		/* new segment */
		FrontSegmentT* new_seg = new FrontSegmentT(pA, pB, pC);
		if (!new_seg) throw eOutOfMemory;

		/* add to front */
		fFrontLines.Append(new_seg);

		/* construct front node */
		if (last_seg != NULL)
		{
			FrontNodeT* node = NewFrontNode(*last_seg, *new_seg);
			fFrontNodes.Append(node);
		}

		/* next */
		last_seg = new_seg;
	}
}  		


/* extend front at the specified points along the given direction - returns new
* cutting facets */
const dArray2DT& Front3DT::NewFacets(const ArrayT<int>& extend_pts,
	const ArrayT<int>& extend_dir)
{
#if __option(extended_errorcheck)
	if (extend_pts.Length() != extend_dir.Length()) throw eSizeMismatch;
#endif

	/* each extension generates 2 new facets */
	fNewFacetMan.SetMajorDimension(2*extend_pts.Length(), false);

	/* generate front segments and facets */
	int facet = 0;
	double facet_coords[3*3];
	double* x_n = facet_coords;
	double*   A = x_n + 3;
	double*   B = A + 3;
	for (int i = 0; i < extend_pts.Length(); i++)
	{
		/* node and sampling point */
		int node = extend_pts[i];
		int dir  = extend_dir[i];
		FrontNodeT& frontnode = *fFrontNodes[node];
	
		/* origin and direction */
		const double* x = frontnode.Coords();
		const double* n = frontnode.Direction(dir);
		
		/* new front point */
		x_n[0] = x[0] + fda*n[0];
		x_n[1] = x[1] + fda*n[1];
		x_n[2] = x[2] + fda*n[2];

		/* get surrounding front segments */
		FrontSegmentT& seg1 = *fFrontLines[node];
		FrontSegmentT& seg2 = *fFrontLines[node+1];
		
		/* reset segments/generate facets */
		seg1.x1(A);
		seg1.x2(B);
		seg1.Reset(A, x_n, B);
		fNewFacets.SetRow(facet++, facet_coords);

		seg2.x1(A);
		seg2.x2(B);
		seg2.Reset(x_n, B, A);
		fNewFacets.SetRow(facet++, facet_coords);
		
		/* reset front node */
		ResetFrontNode(seg1, seg2, frontnode);
	}
	
	/* reset ALL front nodes (probably redundant, except at ends) */
	for (int j = 0; j < fFrontNodes.Length(); j++)
	{
		/* node */
		FrontNodeT& frontnode = *fFrontNodes[j];
	
		/* get surrounding front segments */
		const FrontSegmentT& seg1 = *fFrontLines[j];
		const FrontSegmentT& seg2 = *fFrontLines[j + 1];
		
		/* reset front node */
		ResetFrontNode(seg1, seg2, frontnode);
	}

	/* check front segment length */
	CheckFrontSpacing();

	/* check segment included angles */
	RemoveKinks(fMinAngle);

	return fNewFacets;
}

/* activate kink control by specifying angle > 0 degrees */
void Front3DT::SetKinkAngle(double angle)
{
	if (angle < 0.0)
		fMinAngle = -1.0;
	else
		fMinAngle = angle;
}

/************************************************************************
* Private
************************************************************************/

const dArray2DT& Front3DT::NewFacets_insert(const ArrayT<int>& extend_pts,
	const ArrayT<int>& extend_dir)
{
#if __option(extended_errorcheck)
	if (extend_pts.Length() != extend_dir.Length()) throw eSizeMismatch;
#endif

	/* quick exit */
	if (extend_pts.Length() == 0)
	{
		fNewFacetMan.SetMajorDimension(0, false);
		return fNewFacets;
	}
	else
	{
		//TEMP - work space
		dArray2DT new_facet(3,3);
	
		/* math wrapper */
		iArrayT i_extend_pts;
		i_extend_pts.Alias(extend_pts);

		/* each extension generates 2 new facets */
		fNewFacetMan.SetMajorDimension(2*extend_pts.Length(), false);

		/* generate front segments and facets */
		int facet = 0;
		double facet_coords[3*3];
		double* x_n = facet_coords;
		double*   A = x_n + 3;
		double*   B = A + 3;
		for (int i = 0; i < extend_pts.Length(); i++)
		{
			/* node and sampling point */
			int node = extend_pts[i];
			int dir  = extend_dir[i];
			FrontNodeT& frontnode = *fFrontNodes[node];
	
			/* origin and direction */
			const double* x = frontnode.Coords();
			const double* n = frontnode.Direction(dir);
			
			/* new front point */
			x_n[0] = x[0] + fda*n[0];
			x_n[1] = x[1] + fda*n[1];
			x_n[2] = x[2] + fda*n[2];

			/* get front segments */
			FrontSegmentT& seg1 = *fFrontLines[node];
			FrontSegmentT& seg2 = *fFrontLines[node+1];

			/* first segment */
			seg1.x1(A);
			seg1.x2(B);
			double d_xA = sqrt(Dot(x_n, A));
			if (1 || d_xA < fda*fda_s)
			{
				/* move segment */
				seg1.Reset(A, x_n, B);
				fNewFacets.SetRow(facet++, facet_coords);
			}
			else
			{	
				/* segment midpoint */
				double mid[3];
				Average(A, x_n, mid);

				/* set segments */
				FrontSegmentT* new_seg = new FrontSegmentT(A, mid, B);
				if (!new_seg) throw eOutOfMemory;
				seg1.Reset(mid, x_n, B);
			
				/* new front node */
				FrontNodeT* new_node = NewFrontNode(*new_seg, seg1);
				if (!new_node) throw eOutOfMemory;

				/* insert */
				fFrontLines.InsertAt( new_seg, node);
				fFrontNodes.InsertAt(new_node, node);
			
				/* shift */
				node += 1;
				i_extend_pts += 1;

				/* need extra space */
				fNewFacetMan.SetMajorDimension(fNewFacets.MajorDim() + 1, true);
			
				/* new facets */
				new_facet.SetRow(0, A);		
				new_facet.SetRow(1, mid);		
				new_facet.SetRow(2, B);		
				fNewFacets.SetRow(facet++, new_facet.Pointer());

				new_facet.SetRow(0, mid);		
				new_facet.SetRow(1, x_n);		
				new_facet.SetRow(2, B);		
				fNewFacets.SetRow(facet++, new_facet.Pointer());
			}

			/* next front segment */
			seg2.x1(A);
			seg2.x2(B);
			double d_Bx = sqrt(Dot(B, x_n));
			if (d_Bx < fda*fda_s)
			{
				/* move segment */
				seg2.Reset(x_n, B, A);
				fNewFacets.SetRow(facet++, facet_coords);
			}
			else
			{
				/* segment midpoint */
				double mid[3];
				Average(x_n, B, mid);

				/* set segments */
				seg2.Reset(x_n, mid, A);
				FrontSegmentT* new_seg = new FrontSegmentT(mid, B, A);
				if (!new_seg) throw eOutOfMemory;
			
				/* new front node */
				FrontNodeT* new_node = NewFrontNode(seg2, *new_seg);
				if (!new_node) throw eOutOfMemory;

				/* insert */
				fFrontLines.InsertAt( new_seg, node + 2);
				fFrontNodes.InsertAt(new_node, node + 1);
			
				/* shift */
				i_extend_pts += 1;

				/* need extra space */
				fNewFacetMan.SetMajorDimension(fNewFacets.MajorDim() + 1, true);
			
				/* new facets */
				new_facet.SetRow(0, x_n);		
				new_facet.SetRow(1, mid);		
				new_facet.SetRow(2, A);		
				fNewFacets.SetRow(facet++, new_facet.Pointer());

				new_facet.SetRow(0, mid);		
				new_facet.SetRow(1, B);		
				new_facet.SetRow(2, A);		
				fNewFacets.SetRow(facet++, new_facet.Pointer());
			}
		
			/* reset front node */
			ResetFrontNode(seg1, seg2, frontnode);
		}
	
		/* reset ALL front nodes (probably redundant, except at ends) */
		for (int j = 0; j < fFrontNodes.Length(); j++)
		{
			/* node */
			FrontNodeT& frontnode = *fFrontNodes[j];
			
			/* get surrounding front segments */
			const FrontSegmentT& seg1 = *fFrontLines[j];
			const FrontSegmentT& seg2 = *fFrontLines[j + 1];
		
			/* reset front node */
			ResetFrontNode(seg1, seg2, frontnode);
		}

		/* check segment included angles */
		RemoveKinks(fMinAngle);

		return fNewFacets;
	}
}

/* check for minimum included angle - returns facets added to remove
* kinks in the front */
const dArray2DT& Front3DT::RemoveKinks(double min_angle)
{
	/* new facet data */
	dArray2DT new_facet;
	double N_t1t2[3];

	/* mininum included cosine */
	double cos_min = cos(min_angle*Pi/180.0);

	/* loop over front segments (until no additions) */
	int count = 0;
	int done = 0;
	while (!done && count++ < 10)
	{
		/* set flag */
		done = 1;
		for (int i = 0; i < fFrontLines.Length()-1; i++)
		{
			/* consecutive segments */
			FrontSegmentT& seg1 = *fFrontLines[i];
			FrontSegmentT& seg2 = *fFrontLines[i+1];
			
			/* "upward" crack plane normal */
			const double* N_nt = fFrontNodes[i]->Normal();

			/* internal corner */
			CrossProduct(seg1.N_t(), seg2.N_t(), N_t1t2);
			if (Dot(N_nt, N_t1t2) < 0.0)
			{
				/* included angle */
				double cos_12 = Dot(seg1.N_n(), seg2.N_n());
				
				/* found kink */
				if (cos_12 < cos_min)
				{
					/* to repeat */
					done = 0;
								
					/* reset segment */
					seg1.Reset(seg1.x1(), seg2.x2(), seg2.x1());
					
					/* new facet */					
					int facet_num = fNewFacets.MajorDim();
					fNewFacetMan.SetMajorDimension(facet_num + 1, true);
					new_facet.Set(3, 3, fNewFacets(facet_num));					
					new_facet.SetRow(0, seg1.x1());
					new_facet.SetRow(1, seg2.x2());
					new_facet.SetRow(2, seg2.x1());
					
					/* delete */
					fFrontLines.DeleteAt(i+1);
					fFrontNodes.DeleteAt(i);				
				}							
			}
		}
	}

	return fNewFacets;
}

/* set segment pointers (assuming triangular facets in 3D) */  		
void Front3DT::SetEdgePointers(const double* coords, int edge, const double** ppA,
	const double** ppB, const double** ppC) const
{
	switch (edge)
	{
		case 0:
			*ppA = coords;
			*ppB = *ppA + 3;
			*ppC = *ppB + 3;
			break;
		case 1:
			*ppC = coords;
			*ppA = *ppC + 3;
			*ppB = *ppA + 3;
			break;
		case 2:		
			*ppB = coords;
			*ppC = *ppB + 3;
			*ppA = *ppC + 3;
			break;
		default:
			throw eOutOfRange;
	}
}

/* initialize new front node */
FrontNodeT* Front3DT::NewFrontNode(const FrontSegmentT& seg1,
	const FrontSegmentT& seg2) const
{
	/* end points */
	const double* x12 = seg1.x2();
	const double* x21 = seg2.x1();

	/* segments must be continuous */
	if (fabs(x12[0] - x21[0]) > kSmall ||
	    fabs(x12[1] - x21[1]) > kSmall ||
	    fabs(x12[2] - x21[2]) > kSmall)
	{
		cout << "\n Front3DT::NewFrontNode: segments not continuous" << endl;
		throw eGeneralFail;
	}

	/* average normals */
	double N_n[3], N_t_avg[3];
	Average(seg1.N_n(), seg2.N_n(), N_n);
	Average(seg1.N_t(), seg2.N_t(), N_t_avg);
	
	/* construct local orthogonal basis */
	double N_nt[3], N_t[3];
	CrossProduct(N_n, N_t_avg, N_nt);
	CrossProduct(N_nt, N_n, N_t);

	/* construct node */
	FrontNodeT* node = new FrontNodeT(3, x12, N_n, N_t, fcone, fda*fda_s, fnum_pts);
	if (!node) throw eOutOfMemory;

	return node;
}

/* reset front node */
void Front3DT::ResetFrontNode(const FrontSegmentT& seg1, const FrontSegmentT& seg2,
	FrontNodeT& node) const
{
	/* end points */
	const double* x12 = seg1.x2();
	const double* x21 = seg2.x1();

	/* segments must be continuous */
	if (fabs(x12[0] - x21[0]) > kSmall ||
	    fabs(x12[1] - x21[1]) > kSmall ||
	    fabs(x12[2] - x21[2]) > kSmall)
	{
		cout << "\n Front3DT::NewFrontNode: segments not continuous" << endl;
		throw eGeneralFail;
	}

	/* average normals */
	double N_n[3], N_t_avg[3];
	Average(seg1.N_n(), seg2.N_n(), N_n);
	Average(seg1.N_t(), seg2.N_t(), N_t_avg);
	
	/* construct local orthogonal basis */
	double N_nt[3], N_t[3];
	CrossProduct(N_n, N_t_avg, N_nt);
	CrossProduct(N_nt, N_n, N_t);

	/* reset node */
	node.Reset(3, x12, N_n, N_t, fcone, fda*fda_s, fnum_pts);
}

/* check front spacing */
void Front3DT::CheckFrontSpacing(void)
{
	/* spacing limit = stress sampling distance */
	double max_length = (fda_s > 1.0) ? fda*fda_s : fda;

	/* loop over front segments (until no additions) */
	int count = 0;
	int done = 0;
	while (!done && count++ < 10)
	{
		/* set flag */
		done = 1;

		for (int i = 0; i < fFrontLines.Length(); i++)
		{
			FrontSegmentT& seg = *fFrontLines[i];
			if (seg.Length() > max_length)
			{
				/* to repeat */
				done = 0;
				
				/* get midpoint */
				double x_mid[3];
				seg.MidPoint(x_mid);
				
				/* reset segments */
				FrontSegmentT* new_seg = new FrontSegmentT(seg.x1(), x_mid, seg);
				if (!new_seg) throw eOutOfMemory;
				seg.Reset(x_mid, seg.x2(), seg);
				
				/* generate new front node */
				FrontNodeT* new_node = NewFrontNode(*new_seg, seg);

				/* insert */
				fFrontLines.InsertAt( new_seg, i);
				fFrontNodes.InsertAt(new_node, i);
			}
		}
	}
}
