/* $Id: ContactNodeT.h,v 1.4 2001-04-24 18:17:37 rjones Exp $ */


#ifndef _CONTACT_NODE_T_H_
#define _CONTACT_NODE_T_H_

/* direct members */
#include "SurfaceT.h"
#include "FaceT.h"

/* forward declarations */
class ofstreamT;

class ContactNodeT 
{
  public:

	/* constructor */
	ContactNodeT(SurfaceT& surface, int node_tag);

	/* constructor */
	~ContactNodeT(void);

	/* print data */
	void PrintData(ostream& out);

	/* clear opposing data */
	inline void ClearOpposing(void) 
		{ fOpposingSurface = NULL; fOpposingFace= NULL; fGap = 1.0e8;}

	/* assign opposing point on surface */
	bool AssignOpposing
		(const SurfaceT& opposing_surface, 
		const FaceT& opposing_face,
		double* xi, double g) ;
	void UpdateOpposing(double* xi, double g);

  protected:
        /* data */
	SurfaceT&  fSurface;
	int        fNodeTag; // need to protect the value of the tag?
	const SurfaceT*  fOpposingSurface ; 
	const FaceT*     fOpposingFace ; 
	double     fxi[2] ;
	double     fGap ;

  public:
        /* access functions */ 
	inline const double* Position(void) const
		{return fSurface.Position(fNodeTag);}
	inline const double* Normal(void) const
		{return fSurface.Normal(fNodeTag);}
	inline const double* Tangent1(void) const
		{return fSurface.Tangent1(fNodeTag);}
	inline const double* Tangent2(void) const
		{return fSurface.Tangent2(fNodeTag);}
        inline const SurfaceT* OpposingSurface(void) const
		{return fOpposingSurface;}
        inline const FaceT* OpposingFace(void) const 
		{return fOpposingFace;}
        inline const double* OpposingLocalCoordinates(void) const
		{return fxi;}
        inline const double Gap(void) const 		
		{return fGap;}


  private:

};

#endif /* _CONTACT_NODE_T_H_ */
