/* $Id: ContactNodeT.h,v 1.7 2001-08-06 20:55:12 rjones Exp $ */


#ifndef _CONTACT_NODE_T_H_
#define _CONTACT_NODE_T_H_

/* direct members */
#include "SurfaceT.h"
#include "FaceT.h"
#include "nMatrixT.h"

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

	enum ContactNodeStatusT { kNoProjection = -1,
				  kProjection,
				  kContact,
				  kSlip};

	/* clear opposing data */
	inline void ClearOpposing(void) 
		{ fStatus = kNoProjection; fOpposingSurface = NULL; 
		fOpposingFace= NULL; fGap = 1.0e8;}

	/* assign opposing point on surface */
	bool AssignOpposing
		(const SurfaceT& opposing_surface, 
		const FaceT& opposing_face,
		double* xi, double g) ;

	void UpdateOpposing(double* xi, double g);

	/* can't jump surfaces */
	void AssignOriginal(void);

	void AssignStatus(nMatrixT<dArrayT>& enforcement_parameters);
	inline void AssignOriginalStatus(void)
		{fOriginalStatus = fStatus;}

	inline void ResetStatus(void)
		{fStatus = kNoProjection; fGap = 1.0e8;}

	void ComputeSlip(double* slip);
				  
  protected:
        /* data */
	SurfaceT&  fSurface;
	int        fNodeTag; // need to protect the value of the tag?
	const SurfaceT*  fOpposingSurface ; 
	const FaceT*     fOpposingFace ; 
	double     fxi[2] ;
	double     fGap ;
	int	   fStatus;
	const FaceT*     fOriginalOpposingFace ; 
	double     fxiO[2] ;
	int	   fOriginalStatus;
	

  public:
        /* access functions */ 
	inline const int Tag(void) const
		{return fNodeTag;}
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
        inline const int Status(void) const 		
		{return fStatus;}
        inline const FaceT* OriginalOpposingFace(void) const 
		{return fOriginalOpposingFace;}
        inline const double* OriginalLocalCoordinates(void) const
		{return fxiO;}
        inline const int OriginalStatus(void) const 		
		{return fOriginalStatus;}


  private:

};

#endif /* _CONTACT_NODE_T_H_ */
