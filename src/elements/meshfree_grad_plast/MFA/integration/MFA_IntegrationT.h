// $Id: MFA_IntegrationT.h,v 1.1 2004-06-18 19:22:04 kyonten Exp $
#ifndef _MFA_INTEGRATIONT_H_ 
#define _MFA_INTEGRATIONT_H_ 

/** IMPORTANT NOTE: This integrator will transpose the first argument B, for a resulting 
 *  usual integral of ( B^T.D.B ). Do not transpose before sending it. */
 
/* LHS and RHS integrals are computed */
 
namespace Tahoe 
{

class MFA_IntegrationT 
{

	public:

		MFA_IntegrationT	();
		MFA_IntegrationT	(double Weights);
		void Construct 		(double Weights);

  	//LHS integrals
  	dMatrixT of (dMatrixT &B1, dMatrixT &C, 	dMatrixT &B2 ); 
  	dMatrixT of ( dMatrixT &B1, dArrayT &C, 	double c ); 
  	dMatrixT of ( double c, dArrayT &C, dMatrixT &B ); 
  	double of(double c, double cc, double ccc); 
  	
  	// RHS integrals
  	dArrayT of (dMatrixT &B,  dArrayT &b ); 
  	dArrayT of ( double c,  dArrayT &b ); 
  	double of(double c, double cc); 
  	
  	
				//-- The following are faster (no copying) but less elegant
		
  			//LHS integrals
  			void of (dMatrixT &B1, dMatrixT &C, dMatrixT &B2 );
  			void of (dMatrixT &B1, dArrayT &C, 	double c );
  			void of ( double c, dArrayT &C, 	dMatrixT &B );
  			void of (double c, double cc, double ccc);
  			
  			//RHS integrals
  			void of (dMatrixT &B,  dArrayT &b );
  			void of ( double c,  dArrayT &b );
  	 		void of (double c, double cc);

	protected:
		
		double W; // Gauss Quadrature Weights of ips
				  // Nodal vol??
	private:

		int n_ip,n_rows,n_cols;
};

} // namespace Tahoe 
#endif /* _MFA_INTEGRATIONT_H_ */

