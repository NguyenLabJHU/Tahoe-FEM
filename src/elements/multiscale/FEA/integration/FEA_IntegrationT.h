// $Id: FEA_IntegrationT.h,v 1.4 2003-09-16 16:41:22 raregue Exp $
#ifndef _FEA_INTEGRATIONT_H_ 
#define _FEA_INTEGRATIONT_H_ 

/** IMPORTANT NOTE: This integrator will transpose the first argument B, for a resulting 
 *  usual integral of ( B^T.D.B ). Do not transpose before sending it. */
 
namespace Tahoe {

class FEA_IntegrationT 
{

	public:

		FEA_IntegrationT	();
		FEA_IntegrationT	(FEA_dScalarT &J, FEA_dScalarT &Weights);
		void Construct 		(FEA_dScalarT &J, FEA_dScalarT &Weights);
					
  	dMatrixT of ( FEA_dMatrixT &K );
  	dMatrixT of ( double c, FEA_dMatrixT &K );
  	dMatrixT of ( FEA_dScalarT &s, FEA_dMatrixT &K );
  	dMatrixT of ( double c, FEA_dScalarT &s, FEA_dMatrixT &K );
  	
  	dMatrixT of ( FEA_dVectorT &B1,	double c, FEA_dVectorT &B2 );
  	
  	dMatrixT of ( FEA_dMatrixT &B1, FEA_dMatrixT &B2 );
  	dMatrixT of ( FEA_dMatrixT &B1,	double c, FEA_dMatrixT &B2 );
  	dMatrixT of ( FEA_dMatrixT &B1,	FEA_dScalarT &s, FEA_dMatrixT &B2 );
  	dMatrixT of ( FEA_dMatrixT &B1,	double c, FEA_dScalarT &s, FEA_dMatrixT &B2 );

  	dMatrixT of ( FEA_dMatrixT &B1, FEA_dMatrixT &C, 	FEA_dMatrixT &B2 );
  	dMatrixT of ( FEA_dMatrixT &B1, double c, 				FEA_dMatrixT &C, FEA_dMatrixT &B2 );
  	dMatrixT of ( FEA_dMatrixT &B1, FEA_dScalarT &s, 	FEA_dMatrixT &C, FEA_dMatrixT &B2 );
  	dMatrixT of ( FEA_dMatrixT &B1, double c, FEA_dScalarT &s, FEA_dMatrixT &C, FEA_dMatrixT &B2 );
  	
  	dArrayT of ( FEA_dVectorT &B1,	double c );
				
  	 dArrayT of ( FEA_dMatrixT &B,  FEA_dVectorT &b );
  	 dArrayT of ( FEA_dMatrixT &B,  double &c, FEA_dVectorT &b );
  	 dArrayT of ( FEA_dMatrixT &B,  FEA_dScalarT &s, FEA_dVectorT &b );
  	 dArrayT of ( FEA_dMatrixT &B,  double &c, FEA_dScalarT &s, FEA_dVectorT &b );

				//-- The following are faster (no copying) but less elegant
		
  			void of ( FEA_dMatrixT &K, dMatrixT &k );
  			void of ( double c, FEA_dMatrixT &K, dMatrixT &k );
  			void of ( FEA_dScalarT &s, FEA_dMatrixT &K, dMatrixT &k );
  			void of ( double c, FEA_dScalarT &s, FEA_dMatrixT &K, dMatrixT &k );

  			void of ( FEA_dMatrixT &B1, FEA_dMatrixT &B2, dMatrixT &k );
  			void of ( FEA_dMatrixT &B1,	double c, FEA_dMatrixT &B2, dMatrixT &k );
  			void of ( FEA_dMatrixT &B1,	FEA_dScalarT &s, FEA_dMatrixT &B2, dMatrixT &k );
  			void of ( FEA_dMatrixT &B1,	double c, FEA_dScalarT &s, FEA_dMatrixT &B2, dMatrixT &k );

  			void of ( FEA_dMatrixT &B1, FEA_dMatrixT &C, 	FEA_dMatrixT &B2, dMatrixT &k );
  			void of ( FEA_dMatrixT &B1, double c, 				FEA_dMatrixT &C, FEA_dMatrixT &B2, dMatrixT &k );
  			void of ( FEA_dMatrixT &B1, FEA_dScalarT &s, 	FEA_dMatrixT &C, FEA_dMatrixT &B2, dMatrixT &k );
  			void of ( FEA_dMatrixT &B1, double c, FEA_dScalarT &s, FEA_dMatrixT &C, FEA_dMatrixT &B2, dMatrixT &k );
					
  	 		void of ( FEA_dMatrixT &B,  FEA_dVectorT &b, dArrayT &f );
  	 		void of ( FEA_dMatrixT &B,  double &c, FEA_dVectorT &b, dArrayT &f );
  	 		void of ( FEA_dMatrixT &B,  FEA_dScalarT &s, FEA_dVectorT &b, dArrayT &f );
  	 		void of ( FEA_dMatrixT &B,  double &c, FEA_dScalarT &s, FEA_dVectorT &b, dArrayT &f );

	protected:
		
		FEA_dScalarT W; // Gauss Quadrature Weights of ips
		FEA_dScalarT J; // Jacobian of Isoparametic Mapping dxdeta

	private:

		int n_ip,n_rows,n_cols;
};

} // namespace Tahoe 
#endif /* _FEA_INTEGRATIONT_H_ */

