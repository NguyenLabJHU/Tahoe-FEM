/* $Id: AnisoCorneaIVisco.cpp,v 1.2 2006-09-05 23:10:23 thao Exp $ */
/* created: TDN (01/22/2001) */

#include "AnisoCorneaIVisco.h"
#include <math.h>
#include "toolboxConstants.h"
#include "C1FunctionT.h"
#include "ParameterContainerT.h"

/* point generator */
#include "EvenSpacePtsT.h"

#include "FungType.h"
#include "PowerTrig.h"
#include "ScaledCsch.h"
#include "LinearExponentialT.h"

const double Pi = acos(-1.0);

using namespace Tahoe;

/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
/* constructors */
AnisoCorneaIVisco::AnisoCorneaIVisco(void):
  AnisoCorneaVisco(),
  ParameterInterfaceT("aniso_scalar_visco_cornea")
{}
	
void AnisoCorneaIVisco::TakeParameterList(const ParameterListT& list)
{
	AnisoCorneaVisco::TakeParameterList(list);
	fMod3.Dimension(fNumFibStress);
}
/***********************************************************************
 * Protected
 ***********************************************************************/

void AnisoCorneaIVisco::ComputeCalg(const dSymMatrixT& FiberStretch, 
		const dSymMatrixT& FiberStretch_v,  dMatrixT& Calg, const int pindex)
{
	/*get time step*/
	const double dt = fFSFiberMatSupport->TimeStep();

	/*get flow stress*/
	ComputeFlowStress(FiberStretch, FiberStretch_v, fFlowStress, pindex);
	double sig = sqrt(fFlowStress[0]*fFlowStress[0]+fFlowStress[1]*fFlowStress[1]
			+2.0*fFlowStress[2]*fFlowStress[2]);
	
	/*calculate inverse viscosity*/
	double eta = fViscosity[pindex]->Function(sig); 
	fieta = 1.0/eta;
//	cout << "\nfieta: "<<fieta;

	/*compute Kij = del_ij - dt*V^-1_ik 2dS_k/dCfvt_j*/
	/*Cfvt = {Cfv11, Cfv22, 2 Cfv12}*/

	/*viscosity terms St dyad S*/
	double coeff;
	if ( (sig < kSmall) )
		coeff = 0.0;
	else
		coeff = -fieta/sig*fViscosity[pindex]->DFunction(sig);

	fMod3(0,0) = coeff*fFlowStress[0]*fFlowStress[0];
	fMod3(0,1) = coeff*fFlowStress[0]*fFlowStress[1];
	fMod3(0,2) = coeff*fFlowStress[0]*fFlowStress[2];
		
	fMod3(1,0) = coeff*fFlowStress[1]*fFlowStress[0];
	fMod3(1,1) = coeff*fFlowStress[1]*fFlowStress[1];
	fMod3(1,2) = coeff*fFlowStress[1]*fFlowStress[2];

	fMod3(2,0) = 2.0*coeff*fFlowStress[2]*fFlowStress[0];
	fMod3(2,1) = 2.0*coeff*fFlowStress[2]*fFlowStress[1];
	fMod3(2,2) = 2.0*coeff*fFlowStress[2]*fFlowStress[2];

	/*stress terms*/
	dFlowdCv(FiberStretch, FiberStretch_v, fMod1, pindex);
	fMod1.ToMatrix(fMod2);
		
	fMod2(0,2) *= 0.5;
	fMod2(2,0) *= 2.0;
	fMod2(1,2) *= 0.5;
	fMod2(2,1) *= 2.0;
	
	fiK.MultAB(fMod3,fMod2);
	fiK += fMod2;
	fiK *= -dt*fieta;
	fiK(0,0) += 1.0;
	fiK(1,1) += 1.0;
	fiK(2,2) += 1.0;

	fiK.Inverse();
	

	/*compute G=dS/dCf. Note for this model, dFlowStress/dC = - dSNEQ/dCfv*/
	dFlowdC(FiberStretch, FiberStretch_v, fMod1,  pindex);
	fMod1.ToMatrix(fMod2);

	fMod2(0,2) *= 0.5;
	fMod2(2,0) *= 2.0;
	fMod2(1,2) *= 0.5;
	fMod2(2,1) *= 2.0;


	fG.MultAB(fMod3,fMod2);
	fG += fMod2;
	fG *= dt*fieta;
	
	/*compute dSNEQ_I/dCfv_M K^-1_MK */
	dFlowdC(FiberStretch, FiberStretch_v, fMod1,  pindex);

	fMod2(0,0) = -(fMod1(0,0)*fiK(0,0)+fMod1(0,1)*fiK(1,0)+fMod1(0,2)*fiK(2,0));
	fMod2(1,1) = -(fMod1(1,0)*fiK(0,1)+fMod1(1,1)*fiK(1,1)+fMod1(1,2)*fiK(2,1));
	fMod2(2,2) = -(fMod1(2,0)*fiK(0,2)+fMod1(2,1)*fiK(1,2)+fMod1(2,2)*fiK(2,2));

	fMod2(0,1) = -(fMod1(0,0)*fiK(0,1)+fMod1(0,1)*fiK(1,1)+fMod1(0,2)*fiK(2,1));
	fMod2(0,2) = -(fMod1(0,0)*fiK(0,2)+fMod1(0,1)*fiK(1,2)+fMod1(0,2)*fiK(2,2));
	fMod2(1,2) = -(fMod1(1,0)*fiK(0,2)+fMod1(1,1)*fiK(1,2)+fMod1(1,2)*fiK(2,2));

	fMod2(1,0) = -(fMod1(1,0)*fiK(0,0)+fMod1(1,1)*fiK(1,0)+fMod1(1,2)*fiK(2,0));
	fMod2(2,0) = -(fMod1(2,0)*fiK(0,0)+fMod1(2,1)*fiK(1,0)+fMod1(2,2)*fiK(2,0));
	fMod2(2,1) = -(fMod1(2,0)*fiK(0,1)+fMod1(2,1)*fiK(1,1)+fMod1(2,2)*fiK(2,1));
	
	Calg.MultAB(fMod2, fG);
}

/*local newton loop for viscous stretch tensor*/ 
void AnisoCorneaIVisco::Compute_Cv(const dSymMatrixT& C_n, const dSymMatrixT& C, 
		const dSymMatrixT& C_vn, dSymMatrixT& C_v, const int pindex)
{
	/*compute inverse viscosity  based on last values of stretch matrices, C_n and Cv_n*/
	ComputeFiberStretch(C_n, fFiberStretch_n);
	ComputeFiberStretch(C_vn, fFiberStretch_vn);
	
	/*store Cv_n in fiber frame*/
	const double Cfvn0 = fFiberStretch_vn[0];
	const double Cfvn1 = fFiberStretch_vn[1];
	const double Cfvn2 = fFiberStretch_vn[2];
	
	/*get time step*/
	const double dt = fFSFiberMatSupport->TimeStep();
	
	/*compute current stretch in fiber plane*/
	ComputeFiberStretch(C, fFiberStretch);
	const dSymMatrixT& Cf = fFiberStretch;

	/*compute flow stress based on current values of stretch matrices, C and Cv*/
	ComputeFiberStretch(C_v, fFiberStretch_v);
	dSymMatrixT& Cfv = fFiberStretch_v;
	
	dSymMatrixT& R = fResidual;
	dSymMatrixT& S = fFlowStress;
	
	/*get flow stress*/
	ComputeFlowStress(Cf, Cfv, S, pindex);
	double sig = sqrt(S[0]*S[0]+S[1]*S[1]+2.0*S[2]*S[2]);
//	cout << "\nCfv: "<<Cfv;
//	cout << "\nflow stress: "<<S;
	/*calculate inverse viscosity*/
	double eta = fViscosity[pindex]->Function(sig); 
	fieta = 1.0/eta;
	
	/*residual*/
	R[0] = Cfv[0] - 2.0*dt*fieta*S[0] - Cfvn0;
	R[1] = Cfv[1] - 2.0*dt*fieta*S[1] - Cfvn1;
	R[2] = 2.0*Cfv[2] - 4.0*dt*fieta*S[2] - 2.0*Cfvn2;

//	cout << "\nR: "<<R;
	double error;
	int iteration  = 0;
	
	do 
	{	
		/*compute Kij = del_ij - dt*V^-1_ik 2dS_k/dCfvt_j*/
		/*Cfvt = {Cfv11, Cfv22, 2 Cfv12}*/

		/*viscosity terms St dyad S*/
		double coeff;
		if ( (sig < kSmall) )
			coeff = 0.0;
		else
			coeff = -fieta/sig*fViscosity[pindex]->DFunction(sig);

//		cout << "\ncoeff:"<<coeff;
		
		fMod3(0,0) = coeff*S[0]*S[0];
		fMod3(0,1) = coeff*S[0]*S[1];
		fMod3(0,2) = coeff*S[0]*S[2];
		
		fMod3(1,0) = coeff*S[1]*S[0];
		fMod3(1,1) = coeff*S[1]*S[1];
		fMod3(1,2) = coeff*S[1]*S[2];

		fMod3(2,0) = 2.0*coeff*S[2]*S[0];
		fMod3(2,1) = 2.0*coeff*S[2]*S[1];
		fMod3(2,2) = 2.0*coeff*S[2]*S[2];
		
		/*stress terms*/
		dFlowdCv(Cf, Cfv, fMod1, pindex);
		fMod1.ToMatrix(fMod2);
		
		fMod2(0,2) *= 0.5;
		fMod2(2,0) *= 2.0;
		fMod2(1,2) *= 0.5;
		fMod2(2,1) *= 2.0;
	
//		cout << "\nMod2: "<<fMod2;
//		cout << "\nMod3: "<<fMod3;
					
		fiK.MultAB(fMod3,fMod2);
//		cout << "\niK1: "<<fiK;
		fiK += fMod2;
//		cout << "\niK2: "<<fiK;
		fiK *= -dt*fieta;
//		cout << "\niK3: "<<fiK;
		fiK(0,0) += 1.0;
		fiK(1,1) += 1.0;
		fiK(2,2) += 1.0;

//		cout << "\niK: "<<fiK;
		fiK.Inverse();
		
		/*calculate increment of Cfvt*/
		double dCfvt0 = -(fiK(0,0)*R[0] + fiK(0,1)*R[1] + fiK(0,2)*R[2]);
		double dCfvt1 = -(fiK(1,0)*R[0] + fiK(1,1)*R[1] + fiK(1,2)*R[2]);
		double dCfvt2 = -(fiK(2,0)*R[0] + fiK(2,1)*R[1] + fiK(2,2)*R[2]);
		
		/*calculate update*/
		Cfv[0] += dCfvt0;
		Cfv[1] += dCfvt1;
		Cfv[2] += 0.5*dCfvt2;
		
//		cout << "\niteration: "<<iteration;
//		cout <<setprecision(12)<< "\nCf: "<<Cf;			
//		cout << "\nCfv: "<<Cfv;
//		cout << "\nflowstress: " << S;
//		cout << "\nfieta: "<<fieta;
//		cout << "\nresidual: " << R;
//		cout << "\nS dyad S: "<<fMod3;
//		cout << "\n2dFdCfvt: "<<fMod2;
//		cout << "\nfiK: "<<fiK;
//		cout << "\ndCfv: "<<dCfv0<<"\t"<<dCfv1<<"\t"<<dCfv2;

		/*get flow stress*/
		ComputeFlowStress(Cf, Cfv, S, pindex);
		double sig = sqrt(S[0]*S[0]+S[1]*S[1]+2.0*S[2]*S[2]);
	
		/*calculate inverse viscosity*/
		double eta = fViscosity[pindex]->Function(sig); 
		fieta = 1.0/eta;

		/*residual*/
		R[0] = Cfv[0] - 2.0*dt*fieta*S[0] - Cfvn0;
		R[1] = Cfv[1] - 2.0*dt*fieta*S[1] - Cfvn1;
		R[2] = 2.0*Cfv[2] - 4.0*dt*fieta*S[2] - 2.0*Cfvn2;
		
		error = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
//		cout << "\nerror: "<<error;
		iteration++;
	}while (error > kSmall && iteration < 6);
	if (iteration >= 10) 
		ExceptionT::GeneralFail("AnisoCorneaIVisco::Compute_Cv", 
			"number of iteration exceeds maximum of 10");
}

