/* $Id: PenaltyContactElement2DT.h,v 1.18 2003-05-12 22:01:29 dzeigle Exp $ */
// created by : rjones 2001
#ifndef _PENALTY_CONTACT_ELEMENT_2D_T_H_
#define _PENALTY_CONTACT_ELEMENT_2D_T_H_

/* base classes */
#include "ContactElementT.h"

#include "pArrayT.h"

namespace Tahoe {

class C1FunctionT;

class PenaltyContactElement2DT: public ContactElementT
{
  public:

	/* constructor */
	PenaltyContactElement2DT(const ElementSupportT& support, const FieldT& field);
	
	/* initialize */
	virtual void Initialize(void);

	/* writing output */
	virtual void WriteOutput(void);

    enum EnforcementParametersT { 
                                kConsistentTangent = 0 ,
                                kPenalty ,
                                kPenaltyType ,
			kNumEnfParameters=8}; //this has to the max of all the Penalty types

	enum PenaltyTypesT {
								kLinear = 0,
								kModSmithFerrante,
								kGreenwoodWilliamson,
								kMajumdarBhushan,
		kNumPenaltyTypes};
	
// material constants for the various penalty types
	enum LinParametersT { 
						};
								
	enum SFParametersT {
								kSmithFerranteA=3,
								kSmithFerranteB
						};
	
	enum GWParametersT {
                                kAsperityHeightMean=3,
                                kAsperityHeightStandardDeviation,
                               	kAsperityDensity,
                               	kAsperityTipRadius,
                               	kHertzianModulus  
						};
						
	enum MBParametersT {
								kSigma=3,
								kFractalDimension,
								kRoughnessScale,
								kEPrime,
								kAreaFraction
						};
	 	
  protected:
	/* look-up for symmetric matrix stored as a vector */
	inline int LookUp (int s1,int s2,int n) 
		{return (s1>s2) ? (n*s2+s1) : (n*s1+s2);} 

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
		 	
	/* construct the residual force vector, called before LHS */
	virtual void RHSDriver(void);
	
	/* construct the effective mass matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT);
	
	/* total _real_ area of contact for each surface */
	dArrayT fRealArea; 

    /* penalty models */
	pArrayT<C1FunctionT*> fPenaltyFunctions;


};

} // namespace Tahoe

#endif /* _PENALTY_CONTACT_ELEMENT_2D_T_H_ */

