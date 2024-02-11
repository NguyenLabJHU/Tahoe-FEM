/* $Id: ContactNodeT.h,v 1.20 2003-12-20 01:22:14 rjones Exp $ */
#ifndef _CONTACT_NODE_T_H_
#define _CONTACT_NODE_T_H_

/* direct members */
#include "ContactSurfaceT.h"
#include "FaceT.h"
#include "nMatrixT.h"

namespace Tahoe {

/* forward declarations */
class ofstreamT;

class ContactNodeT 
{
  public:

	/* constructor */
	ContactNodeT(ContactSurfaceT& surface, int node_tag);

	/* constructor */
	~ContactNodeT(void);

	/* print data */
	void PrintData(ostream& out);

	enum ContactNodeStatusT { 	kNoProjection = -1,
								kProjection};

	/** initialize data */
	inline void Initialize(void) { 
			fStatus = kNoProjection; 
			fEnforcementStatus = -10;
			fOpposingSurface = NULL; 
			fOpposingFace= NULL; 
			fGap = 1.0e8;
			fOriginalOpposingFace    = NULL;
			fxi[0]           = 0.0 ;
			fxi[1]           = 0.0 ;
			fPressure        = 0.0;
			fTangTraction[0] = 0.0;
			fTangTraction[1] = 0.0;
	}

	/** assign opposing point on surface */
	bool AssignOpposing
		(const SurfaceT& opposing_surface, 
		const FaceT& opposing_face,
		double* xi, double g) ;

	/** assign status at the beginning of the time step */
	inline void AssignOriginal(void) {
		fOriginalOpposingFace = fOpposingFace;
		fxiO[0] = fxi[0];
		fxiO[1] = fxi[1];
	}

	/** reset (current) status when node loses contact */
	inline void ResetStatus(void) {
			fStatus = kNoProjection; 
			fGap = 1.0e8;
	}

	/** assign value of gap for last timestep */
	inline void AssignLastGap(void) {
			fLastGap = fGap; 
			if (fLastGap < fMinGap )  fMinGap=fLastGap  ;
	}

	void ComputeSlip(double* slip);
	double ComputeSlip(void);

	
  protected:
	/* data */
	ContactSurfaceT&  fSurface;
	/* local node number in surface */
	int        fNodeTag; // need to protect the value of the tag?
	const ContactSurfaceT*  fOpposingSurface ; 
	const FaceT*     fOpposingFace ; 
	double     fxi[2] ;
	double     fGap ;
	double     fLastGap ;
	double     fMinGap ;
	int	   fStatus;
	int	   fOriginalStatus;
	int	   fEnforcementStatus;
	const FaceT*     fOriginalOpposingFace ; 
	double     fxiO[2] ;
	double fPressure;
	double fTangTraction[2];
	

  public:
	/* access functions */ 
	bool HasProjection(void) {return fStatus > kNoProjection;}
	bool HasMultiplier(void) {return fSurface.HasMultiplier(fNodeTag);}
	double& Pressure(void) {return fSurface.Multiplier(fNodeTag,0);}
	double& nPressure(void) {return fPressure;}
	double& TangTraction(void) {return fSurface.Multiplier(fNodeTag,1);}
	double& TangTraction(int i) {return fSurface.Multiplier(fNodeTag,i);}
	double& nTangTraction(int i) {return fTangTraction[i];}
	inline const int Tag(void) const
		{return fNodeTag;}
	inline const double* Position(void) const
		{return fSurface.Position(fNodeTag);}
	inline const double* Normal(void) const
		{return fSurface.Normal(fNodeTag);}
	inline const double* Tangent(void) const
		{return fSurface.Tangent1(fNodeTag);}
	inline const double* Tangent1(void) const
		{return fSurface.Tangent1(fNodeTag);}
	inline const double* Tangent2(void) const
		{return fSurface.Tangent2(fNodeTag);}
	inline const ContactSurfaceT* OpposingSurface(void) const
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
	inline int& EnforcementStatus(void) 
		{return fEnforcementStatus;}
	inline double& LastGap(void) 
		{return fLastGap;}
	inline double& MinGap(void) 
		{return fMinGap;}

  private:

};

} // namespace Tahoe 
#endif /* _CONTACT_NODE_T_H_ */
