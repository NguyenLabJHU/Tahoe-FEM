// $Id: Coupled_EqsT.h,v 1.1 2005-01-26 02:13:22 kyonten Exp $
#ifndef _COUPLED_EQS_T_H_ 
#define _COUPLED_EQS_T_H_ 

#include "StringT.h"
#include "MFGP_Data_ProcessorT.h"
#include "D3MeshFreeShapeFunctionT.h"
#include "D3MeshFreeSupportT.h"

/* material model */
#include "GRAD_MRSSKStV.h"
#include "GRAD_MRSSKStV2D.h"

namespace Tahoe 
{

/** This class contains methods which build stiffness matricies 
 *  and formulate the non-linear Newton-Raphson equations Kd = -R
 *  for a coupled approach to implementation of 
 *  the Balance of Linear Momentum and Consistency Condition in weak form 
 **/

class Coupled_EqsT 
{
public:
	
	enum Eqn_TypeT { kCoupled_Eqs };

	/* constructors */
	Coupled_EqsT(void);
		
	/* destructor */				
	//~Coupled_EqsT(void);

	void Initialize(int&, D3MeshFreeShapeFunctionT*, D3MeshFreeShapeFunctionT*, 
					int& fTime_Step, double fdelta_t = 0.0); 

  	void Form_KUU_KULam(dMatrixT& Kuu, dMatrixT& Kulam); // add delta_t for dynamics
  	void Form_KLamU_KLamLam(dMatrixT& Klamu, dMatrixT& Klamlam);
  	void Form_FU_int(dArrayT& Fu_int); 
  	void Form_FLambda_int(dArrayT& Flambda_int);
	void Form_B_List(void);  // Strain Displacement Matricies
 	void Form_C_List(void);  // Constants List

	protected:
		dMatrixT B1, B3, B4, psi_lam;  
  		dMatrixT Cuu1, Cuu2, Culam1, Culam2;
  		dMatrixT Clamu1, Clamu2, Clamlam1, Clamlam2;
  		double yield;
		dSymMatrixT stress;
		double delta_t;
		int time_step;
		int ip, n_sd, n_str, n_en_displ, n_en_plast, 
			n_sd_x_n_en_displ, n_sd_x_n_en_plast, n_state;
	
	protected:
		MFGP_Data_ProcessorT Data_Pro;
		GRAD_MRSSKStV2D Mat2D;
		GRAD_MRSSKStV Mat3D;
};


} // namespace Tahoe 
#endif /* _COUPLED_EQS_T_H_*/

