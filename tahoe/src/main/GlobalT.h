/* $Id: GlobalT.h,v 1.8 2002-07-05 22:28:07 paklein Exp $ */
/* created: paklein (02/03/1999) */

#ifndef _GLOBAL_T_H_
#define _GLOBAL_T_H_

#include "Environment.h"
#include "ios_fwd_decl.h"

namespace Tahoe {

/** class to handle "global" enumerated types */
class GlobalT
{
public:

	/** types of analysis */
	enum AnalysisCodeT {
	         kNoAnalysis = 0,
		      kLinStatic = 1,
		     kLinDynamic = 2,
		       kNLStatic = 3,
		      kNLDynamic = 4,
		             kDR = 5,  /**< this will be converted to a nonlinear solver method */
		  kLinExpDynamic = 6,
		   kNLExpDynamic = 7,
		  kLinStaticHeat = 19, /**< linear static heat conduction */
		   kLinTransHeat = 20, /**< linear transient heat conduction */
		            kPML = 30, /**< perfectly matched layer formulation */
			 kMultiField = 99  /**< generalized analysis code */
		   };
		
	/** stream extraction operator */
	friend istream& operator>>(istream& in, GlobalT::AnalysisCodeT& code);
	
	/** deprecated analysis codes */
	enum OldAnalysisCodeT {
		kVarNodeNLStatic = 15, /**< variables nodes supported through ModelManagerT */
		kVarNodeNLExpDyn = 16, /**< variables nodes supported through ModelManagerT */
		       kCBStatic = 8, /**< converted to KBC controller: PAK (12/10/2000) */
		   kAugLagStatic = 17, /**< moved to general support of element DOF: PAK (08/22/2001) */
	     kNLStaticKfield = 11, /**< converted to KBC controller: PAK (09/10/2000) */
		 kNLExpDynKfield = 18  /**< converted to KBC controller: PAK (09/10/2000) */
		};

// Currently nonlinear <=> large deformation, will probably
// need to break this into:
//    (a) linear, small deformation
//    (b) nonlinear, small deformation
//    (c) (nonlinear) large deformation

	/** analysis stage */
	enum StateT {
	            kNone = 0,
        kConstruction = 1,
      kInitialization = 2,
    kInitialCondition = 3,
         kReadRestart = 4,
	        kInitStep = 5,
	         kFormRHS = 6,
	         kFormLHS = 7,
	       kResetStep = 8,
	       kCloseStep = 9,
	     kWriteOutput =10,
		kWriteRestart =11,
		   kException =12,
		 kDestruction =13};

	/** global system types, ordered so n_1 > n_2 implies that n_1 is a
	 * superset of n_2. */
	enum SystemTypeT {
	       kUndefined =-1,
		    kDiagonal = 0,
		   kSymmetric = 1,
		kNonSymmetric = 2};

	/** returns the type with higher restrictions */
	static SystemTypeT MaxPrecedence(SystemTypeT code1, SystemTypeT code2);

	/** relaxation level */
	enum RelaxCodeT {
		kNoRelax = 0, /**< do nothing */
		   kReEQ = 1, /**< reset global equation numbers, but still at force equilirbium */
		  kRelax = 2, /**< relax, ie. re-find equilibrium */
	  kReEQRelax = 3  /**< reset global equation numbers and relax */ };

	/** returns flag with precedence */
	static RelaxCodeT MaxPrecedence(RelaxCodeT code1, RelaxCodeT code2);

	/** equation numbering scope mainly for parallel solvers */
	enum EquationNumberScopeT {
		kLocal  = 0, /**< equations numbered per processor */
		kGlobal = 1  /**< equations numbered over entire system */};
};

/* inlines */
inline GlobalT::SystemTypeT GlobalT::MaxPrecedence(SystemTypeT code1, SystemTypeT code2)
{
	return (code1 > code2) ? code1 : code2;
}

} // namespace Tahoe 
#endif // _GLOBAL_T_H_
