// $Id: APS_VariableT.h,v 1.11 2004-02-04 00:40:44 raregue Exp $
#ifndef _APS_VARIABLE_T_H_ 
#define _APS_VARIABLE_T_H_ 

#include "APS_FEA.h"
#include "APS_EnumT.h"

namespace Tahoe {

/** APS_VariableT: This class contains methods pertaining to kinematics of
 * a dual field formulation. These include gradients such as grad(u).
 **/
 
class APS_VariableT
{
  public:

		/** constructor */
		APS_VariableT 	(void) { }
		APS_VariableT 	(const FEA_dMatrixT& grad_u, const FEA_dMatrixT& grad_u_surf, const FEA_dVectorT& gammap, 
						 const FEA_dVectorT& gamma_p_surf, const FEA_dMatrixT& grad_gammap, FEA_dVectorT& state);

		/** data initialization */
		void Construct 	(const FEA_dMatrixT& grad_u, const FEA_dMatrixT& grad_u_surf, const FEA_dVectorT& gammap, 
						 const FEA_dVectorT& gamma_p_surf, const FEA_dMatrixT& grad_gammap, FEA_dVectorT& state);

		/** delete variables */
		void Delete_Vars	( void );

		/** Print Routine */
		void Print  (void);
		void Print  (char*);
		
		/** Compute and store ... : recursive routine */ 
 	 	void Allocate_and_Compute_Variables(APS::VarT_vector kVariable);
 	 	void Allocate_and_Compute_Variables(APS::VarT_matrix kVariable);

 	 	/** Retrieve either grad_u, gammap from class workspace **/
		const FEA_dVectorT& Get(APS::VarT_vector variable); 
		const FEA_dMatrixT& Get(APS::VarT_matrix variable); 
		
		/** Put state variables in class workspace **/
		void Put(APS::VarT_vector variable, FEA_dVectorT&); 
		
		/** Update state variables from class workspace **/
		void Update(APS::VarT_vector variable, FEA_dVectorT&); 

		/** Fill (*this) with a+b */
		void SumOf (APS_VariableT &a, APS_VariableT &b); 
		void operator  =  (const APS_VariableT &a);
		void operator +=  (const double &a);
		void operator -=  (const double &a);
		void operator *=  (const double &a);
		void operator /=  (const double &a); 

 		//protected:

    	ArrayT <FEA_dVectorT> fVars_vector; //Variables : gammap stored here
    	ArrayT <FEA_dMatrixT> fVars_matrix; //Variables : grad_u, grad_gammap stored here

	private:

		int n_vars_vector, n_vars_matrix;
};

//---------------------------------------------------------------------
} // namespace Tahoe 
#endif /* _APS_VARIABLE_T_H_ */



