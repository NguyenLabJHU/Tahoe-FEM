/* $Id: vector3D.h,v 1.9 2003-11-21 22:54:35 paklein Exp $ */
#ifndef _VECTOR_3D_H_
#define _VECTOR_3D_H_

// need only basic ops: dot,scale,cross,add
/* 3D vector functions */

namespace Tahoe {

inline static void Cross(const double* v, const double* w, double* vXw)
{
        vXw[0] = v[1]*w[2] - v[2]*w[1];
        vXw[1] = v[2]*w[0] - v[0]*w[2];
        vXw[2] = v[0]*w[1] - v[1]*w[0];
};

inline static double Dot(const double* v, const double* w)
{       return v[0]*w[0] + v[1]*w[1] + v[2]*w[2]; };

inline static void Diff(const double* start, const double* end, double* v)
{
        v[0] = end[0] - start[0];
        v[1] = end[1] - start[1];
        v[2] = end[2] - start[2];
};

inline static void Add(const double* v1, const double* v2, double* v)
{
        v[0] = v1[0] + v2[0];
        v[1] = v1[1] + v2[1];
        v[2] = v1[2] + v2[2];
};

inline static void Add(const double a1,const double* v1, 
                       const double a2,const double* v2, double* v)
{
        v[0] = a1*v1[0] + a2*v2[0];
        v[1] = a1*v1[1] + a2*v2[1];
        v[2] = a1*v1[2] + a2*v2[2];
};

inline static void Ave(const double* v1,const double* v2, const double* v3, double* v)
{
        v[0] = 0.25*(v1[0] + v2[0] + v3[0] );
        v[1] = 0.25*(v1[1] + v2[1] + v3[1] );
        v[2] = 0.25*(v1[2] + v2[2] + v3[2] );
};


inline static void Ave(const double* v1, const double* v2, const double* v3, const double* v4, double* v)
{
        v[0] = 0.25*(v1[0] + v2[0] + v3[0] + v4[0]);
        v[1] = 0.25*(v1[1] + v2[1] + v3[1] + v4[1]);
        v[2] = 0.25*(v1[2] + v2[2] + v3[2] + v4[2]);
};

inline static double Magnitude(const double* v)
{
        return  sqrt (Dot(v,v)) ;
};


inline static void Normalize(double* v)
{
        double scale = 1.0/ Magnitude(v) ;
        v[0] *= scale ;
        v[1] *= scale ;
        v[2] *= scale ;
};

inline static double TripleProduct(const double* v1,const double* v2,const double* v3)
{
        return    v1[1]*v2[2]*v3[0] - v1[2]*v2[1]*v3[0]
               +  v1[2]*v2[0]*v3[1] - v1[0]*v2[2]*v3[1]
               +  v1[0]*v2[1]*v3[2] - v1[1]*v2[0]*v3[2];
};

inline static void Proj(const double* v,const double* n, double* proj_v)
{
        double dot = v[0]*n[0] + v[1]*n[1] + v[2]*n[2];
        proj_v[0] = proj_v[0] - dot * n[0] ;
        proj_v[1] = proj_v[1] - dot * n[1] ;
        proj_v[2] = proj_v[2] - dot * n[2] ;

};

} // namespace Tahoe 
#endif /* _VECTOR_3D_H_ */
