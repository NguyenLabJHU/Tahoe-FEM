/* $Id: ContactNodeT.h,v 1.16 2003-02-03 04:40:18 paklein Exp $ */
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

	/** clear opposing data */
	inline void ClearOpposing(void) 
		{ fStatus = kNoProjection; fOpposingSurface = NULL; 
		fOpposingFace= NULL; fGap = 1.0e8;}

	/** assign opposing point on surface */
	bool AssignOpposing
		(const SurfaceT& opposing_surface, 
		const FaceT& opposing_face,
		double* xi, double g) ;

	/** assign status at the beginning of the time step */
	void AssignOriginal(void);

	/** reset (current) status when node loses contact */
	inline void ResetStatus(void)
		{fStatus = kNoProjection; fGap = 1.0e8;}

	void ComputeSlip(double* slip);

	
  protected:
	/* data */
	ContactSurfaceT&  fSurface;
	/* local node number in surface */
	int        fNodeTag; // need to protect the value of the tag?
	const ContactSurfaceT*  fOpposingSurface ; 
	const FaceT*     fOpposingFace ; 
	double     fxi[2] ;
	double     fGap ;
	int	   fStatus;
	int	   fOriginalStatus;
	int	   fEnforcementStatus;
	const FaceT*     fOriginalOpposingFace ; 
	double     fxiO[2] ;
	double fPressure;
	

  public:
	/* access functions */ 
	bool HasProjection(void) {return fStatus > kNoProjection;}
	bool HasMultiplier(void) {return fSurface.HasMultiplier(fNodeTag);}
	double& Pressure(void) {return fSurface.Multiplier(fNodeTag,0);}
	double& nPressure(void) {return fPressure;}
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

  private:

};

} // namespace Tahoe 
#endif /* _CONTACT_NODE_T_H_ */
