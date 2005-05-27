/* $Id: GlobalT.h,v 1.13.8.1 2005-05-27 19:55:20 paklein Exp $ */
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
		  kLinStaticHeat = 19, /**< linear steady-state heat conduction */
		   kLinTransHeat = 20, /**< linear transient heat conduction */
		   kNLStaticHeat = 21, /**< nonlinear steady-state heat conduction */
		    kNLTransHeat = 22, /**< nonlinear transient heat conduction */
		            kPML = 30, /**< perfectly matched layer formulation */
			 kMultiField = 99  /**< generalized analysis code */
		   };
		
	/** stream extraction operator */
//	friend istream& operator>>(istream& in, GlobalT::AnalysisCodeT& code);

	/** solver codes */
	enum SolverTypeT {kNewtonSolver = 0, /**< standard Newton solver */
                   kK0_NewtonSolver = 1, /**< initial tangent, Newton solver */
                   kModNewtonSolver = 2, /**< modified Newton solver (development) */
                    kExpCD_DRSolver = 3, /**< central difference, dynamic relaxation */
                   kNewtonSolver_LS = 4, /**< Newton solver with line search */
                      kPCGSolver_LS = 5, /**< preconditioned, nonlinear conjugate gradient */
                  kiNewtonSolver_LS = 6, /**< interactive Newton solver (with line search) */
                         kNOXSolver = 7, /**< NOX library solver */
                      kLinearSolver = 8, /**< linear problems */
                          kDRSolver = 9,  /**< dynamic relaxation */                               
                  kNewtonSolver_LSX = 104 /**< temporary extention to GlobalT::kNewtonSolver_LS */
                               };
	
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
	     kRelaxSystem =10,
	     kWriteOutput =11,
		kWriteRestart =12,
		   kException =13,
		 kDestruction =14};

	/** global system types, ordered so n_1 > n_2 implies that n_1 is a
	 * superset of n_2. */
	enum SystemTypeT {
	       kUndefined =-1,
		    kDiagonal = 0,
		   kSymmetric = 1,
		kNonSymmetric = 2};

	/** returns the type with higher restrictions */
	static SystemTypeT MaxPrecedence(SystemTypeT code1, SystemTypeT code2);

	/** result flag for InitStep */
	enum InitStatusT {
		kContinue = 0,      /**< no signal needed */
		kNewInteractions,   /**< new interaction with existing equation numbers */
		kAssignEquations,   /**< (re-)assign equation numbers */
	};

	/** returns the type with higher restrictions */
	static InitStatusT MaxPrecedence(InitStatusT status_1, InitStatusT status_2);

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

	/** amount of runtime logging information */
	enum LoggingT {
		kVerbose = 0,
		kModerate = 1,
		kSilent = 2
	};
	static LoggingT int2LoggingT(int i);
};

/* inlines */
inline GlobalT::SystemTypeT GlobalT::MaxPrecedence(SystemTypeT code1, SystemTypeT code2)
{
	return (code1 > code2) ? code1 : code2;
}

} // namespace Tahoe 
#endif // _GLOBAL_T_H_
