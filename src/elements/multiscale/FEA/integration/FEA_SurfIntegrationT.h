// $Id: FEA_SurfIntegrationT.h,v 1.2 2003-10-12 23:47:51 raregue Exp $
#ifndef _FEA_SURFINTEGRATIONT_H_ 
#define _FEA_SURFINTEGRATIONT_H_ 

/** IMPORTANT NOTE: This integrator will transpose the first argument B, for a resulting 
 *  usual integral of ( B^T.D.B ). Do not transpose before sending it. */
 
namespace Tahoe {

class FEA_SurfIntegrationT 
{

	public:

		FEA_SurfIntegrationT	();
		FEA_SurfIntegrationT	(FEA_dScalarT &J, FEA_dScalarT &Weights);
		void Construct 			(FEA_dScalarT &J, FEA_dScalarT &Weights);
  	
  		dMatrixT of ( FEA_dVectorT &B1,	double &c, FEA_dVectorT &B2 );
  	
  		dArrayT of ( FEA_dVectorT &B1, double &c, FEA_dScalarT &s );

	protected:
		
		FEA_dScalarT W; // Gauss Quadrature Weights of ips
		FEA_dScalarT J; // Jacobian of Isoparametic Mapping dxdeta

	private:

		int n_ip,n_rows,n_cols;
};

} // namespace Tahoe 
#endif /* _FEA_SURFINTEGRATIONT_H_ */

