/* $Id: vector2D.h,v 1.1 2001-04-16 17:30:52 rjones Exp $ */

#ifndef _VECTOR_2D_H_
#define _VECTOR_2D_H_


/* 2D vector functions */
inline static void LCross(const double* v,  double* vXe3)
{
        vXe3[0] = -v[1];
        vXe3[1] =  v[0];
};

inline static void RCross(const double* v,  double* e3Xv)
{
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

inline static void Ave(const double* v1, const double* v2, double* v)
{
        v[0] = 0.5 * ( v1[0] - v2[0]);
        v[1] = 0.5 * ( v1[1] - v2[1]);
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

#endif /* _VECTOR_2D_H_ */
