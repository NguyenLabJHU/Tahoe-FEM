// $Id: MFGP_VariableT.h,v 1.1 2004-06-18 19:18:41 kyonten Exp $
#ifndef _MFGP_VARIABLE_T_H_ 
#define _MFGP_VARIABLE_T_H_ 

#include "MFGP_MFA.h"
#include "MFGP_EnumT.h"

namespace Tahoe 
{

/** MFGP_VariableT: This class contains methods pertaining to kinematics of
 * a dual field formulation. These include gradients such as grad(u).
 **/
 
class MFGP_VariableT
{
  public:

		/** constructor */
		MFGP_VariableT 	(void) { }
		MFGP_VariableT 	(const FEA_dMatrixT& grad_u, const FEA_dVectorT& gammap, 
						 const FEA_dMatrixT& grad_gammap, FEA_dVectorT& state);

		/** data initialization */
		void Construct 	(const FEA_dMatrixT& grad_u, const FEA_dVectorT& gammap, 
						 const FEA_dMatrixT& grad_gammap, FEA_dVectorT& state);

		/** delete variables */
		void Delete_Vars	( void );

		/** Print Routine */
		void Print  (void);
		void Print  (char*);
		
		/** Compute and store ... : recursive routine */ 
 	 	void Allocate_and_Compute_Variables(MFGP::VarT_vector kVariable);
 	 	void Allocate_and_Compute_Variables(MFGP::VarT_matrix kVariable);

 	 	/** Retrieve either grad_u, gammap from class workspace **/
		const FEA_dVectorT& Get(MFGP::VarT_vector variable); 
		const FEA_dMatrixT& Get(MFGP::VarT_matrix variable); 
		
		/** Put state variables in class workspace **/
		void Put(MFGP::VarT_vector variable, FEA_dVectorT&); 
		
		/** Update state variables from class workspace **/
		void Update(MFGP::VarT_vector variable, FEA_dVectorT&); 

		/** Fill (*this) with a+b */
		void SumOf (MFGP_VariableT &a, MFGP_VariableT &b); 
		void operator  =  (const MFGP_VariableT &a);
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
#endif /* _MFGP_VARIABLE_T_H_ */



