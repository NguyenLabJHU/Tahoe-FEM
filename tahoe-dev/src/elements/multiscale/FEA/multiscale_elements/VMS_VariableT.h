// $Id: VMS_VariableT.h,v 1.2 2003-02-03 04:40:25 paklein Exp $
#ifndef _VMS_VARIABLE_T_H_ 
#define _VMS_VARIABLE_T_H_ 

namespace Tahoe {

/** VMS_VariableT: This class contains methods pertaining to kinematics of
 * a dual field formulation. These include deformation gradients Fe and Fp
 * and gradients such as Grad(ue) and Grad(up) as examples.
 * Sandia National Laboratory and the University of Michigan */
class VMS_VariableT
{
  public:

		/** constructor */
		VMS_VariableT 	(void) { }
		VMS_VariableT 	(const FEA_dMatrixT& GRAD_ua, const FEA_dMatrixT& GRAD_ub);

		/** data initialization */
		void Construct 	(const FEA_dMatrixT& GRAD_ua, const FEA_dMatrixT& GRAD_ub);

		/** delete variables */
		void Delete_Vars	( void );

		/** Print Routine */
		void Print  (void);
		void Print  (char*);

		/** Compute and store F^-1, Fa, Fb, GRAD(ua), GRAD(ub) : recursive routine */ 
 	 	void Allocate_and_Compute_Variables(VMS::VarT kVariable);

 	 	/** Retrieve either Fa, Fb, GRAD_ua, or GRAD_ub from class workspace **/
		const FEA_dMatrixT& Get(VMS::VarT variable); 

		/** Fill (*this) with a+b */
		void SumOf (VMS_VariableT &a,VMS_VariableT &b); 
		void operator  =  (const VMS_VariableT &a);
		void operator +=  (const double &a);
		void operator -=  (const double &a);
		void operator *=  (const double &a);
		void operator /=  (const double &a); 

 		//protected:

    ArrayT <FEA_dMatrixT> fVars;         /** Variables : Fa, Fb, GRAD_ua, GRAD_ub, etc stored here */

	private:

		int n_vars;
};

//---------------------------------------------------------------------
} // namespace Tahoe 
#endif /* _VMS_VARIABLE_T_H_ */



