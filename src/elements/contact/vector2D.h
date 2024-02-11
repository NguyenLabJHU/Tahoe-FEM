/* $Id: vector2D.h,v 1.9 2003-02-03 04:40:18 paklein Exp $ */
#ifndef _VECTOR_2D_H_
#define _VECTOR_2D_H_

/* 2D vector functions */

namespace Tahoe {

inline static void LCross(const double* v,  double* vXe3)
{
        vXe3[0] = -v[1];
        vXe3[1] =  v[0];
};

inline static void RCross(const double* v,  double* e3Xv)
{
// this follows the exoII sideset convection CCW coor
        e3Xv[0] =  v[1];
        e3Xv[1] = -v[0];
};


inline static double Dot(const double* v1, const double* v2)
{
        return v1[0]*v2[0] + v1[1]*v2[1];
};

inline static void Diff(const double* start, const double* end, double* v)
{
        v[0] = end[0] - start[0];
        v[1] = end[1] - start[1];
};

inline static void Add(const double* v1, const double* v2, double* v)
{
        v[0] = v1[0] + v2[0];
        v[1] = v1[1] + v2[1];
};

inline static void Add(double* v1, double scale, const double* v2)
{
        v1[0] +=  scale*v2[0];
        v1[1] +=  scale*v2[1];
};


inline static void Ave(const double* v1, const double* v2, double* v)
{
        v[0] = 0.5 * ( v1[0] + v2[0]);
        v[1] = 0.5 * ( v1[1] + v2[1]);
};


inline static double Mag(const double* v)
{
        return  sqrt (Dot(v,v)) ;
};

inline static void Normalize(double* v)
{
        double scale = 1.0/ Mag(v) ;
        v[0] *= scale ;
        v[1] *= scale ;
};

inline static void Scale(double scale, double* v)
{
        v[0] *= scale ;
        v[1] *= scale ;
};


inline static void Proj(double* v, double* n, double* proj_v)
{
	double dot = v[0]*n[0] + v[1]*n[1];
	proj_v[0] = v[0] - dot * n[0] ;
	proj_v[1] = v[1] - dot * n[1] ;
};

#if 0
inline static void Permutation(dMatrixT& p_mat)
{
	p_mat[0][0] = 0.0; p_mat[0][1] =-1.0;
	p_mat[1][0] = 1.0; p_mat[1][1] = 0.0;
};
#endif


} // namespace Tahoe 
#endif /* _VECTOR_2D_H_ */
