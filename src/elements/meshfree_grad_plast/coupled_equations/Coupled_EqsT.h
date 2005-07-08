// $Id: Coupled_EqsT.h,v 1.3 2005-07-08 01:13:27 kyonten Exp $
#ifndef _COUPLED_EQS_T_H_ 
#define _COUPLED_EQS_T_H_ 

#include "StringT.h"
#include "D3MeshFreeShapeFunctionT.h"
#include "D3MeshFreeSupportT.h"
#include "MFGPMaterialT.h"

namespace Tahoe 
{

/** this class contains methods which build stiffness matricies 
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

	/** initialize shape functions, material, etc */
	void Initialize(int&, D3MeshFreeShapeFunctionT*, D3MeshFreeShapeFunctionT*, 
					MFGPMaterialT* curr_mat); 

  	/** form the element stiffness matrices KUU & KULambda */
  	void Form_KUU_KULam(dMatrixT& Kuu, dMatrixT& Kulam); // add delta_t for dynamics
  	
  	/** form the element stiffness matrices KLambdaU & KLambdaLambda */
  	void Form_KLamU_KLamLam(dMatrixT& Klamu, dMatrixT& Klamlam);
  	
  	/** calculate the internal force contribution ("-k*d") */
  	void Form_FU_int(dArrayT& Fu_int);
  	 
  	/** calculate the internal force contribution ("-k*lambda") */
  	void Form_FLambda_int(dArrayT& Flambda_int);
	
	/** form B matrices */
	void Form_B_List(void);  // strain displacement matricies
 	
 	/** form C matrices */
 	void Form_C_List(void);  // C matrices
 	
 	/** set the \e B1 matrix using the given shape function derivatives
	 * \param first derivatives of shape function derivatives: [nsd] x [nen]
	 * \param B1 destination for B1: [nstr]x[nsd*nnd] */
	void Set_B1(dMatrixT& B1);
		
	/** set the \e B3 matrix using the given shape function derivatives
	 * \param third derivatives of shape function derivatives: [nsd*nsd] x [nen]
	 * \param B3 destination for B3: [nstr]x[nsd*nnd] */
	void Set_B3(dMatrixT& B3);
                        
    /** set the \e psi_lam matrix using the given shape function
	 * \param shape function: [1] x [nnd]
	 * \param psi_lam destination for psi_lam: [1]x[nnd] */
    /* psi_lam: [1]x[nnd] */
    void Set_psi_lam(dMatrixT& psi_lam);
        
    /** set the \e B4 matrix using the given shape function derivatives
	 * \param second derivatives of shape function derivatives: [nstr] x [nen]
	 * \param B4 destination for B4: [1]x[nnd] */
	void Set_B4(dMatrixT& B4);

	protected:
		dMatrixT B1, B3, B4, psi_lam;  
  		dMatrixT Cuu1, Cuu2, Culam1, Culam2;
  		dMatrixT Clamu1, Clamu2, Clamlam1, Clamlam2;
		double delta_t;
		int time_step;
		int ip, n_sd, n_str, n_en_displ, n_en_plast, 
			n_sd_x_n_en_displ, n_sd_x_n_en_plast, n_state;
		
		const double* N;
	    dArray2DT DN, DDN, DDDN;
		
		MFGPMaterialT* CurrMat;
		D3MeshFreeShapeFunctionT* ShapeDispl;
		D3MeshFreeShapeFunctionT* ShapePlast;
};


} // namespace Tahoe 
#endif /* _COUPLED_EQS_T_H_*/

