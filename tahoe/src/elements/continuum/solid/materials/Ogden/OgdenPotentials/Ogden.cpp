/* created:   TDN (5/31/2001) */
/* Phi(I1,J) = mu/2*(I1-3)+kappa/4*(J^2-1-2*ln(J)) */
/* I1 = trace(C); J=sqrt(det(C)) */

#include "Ogden.h"

#include "ExceptionT.h"
#include "ifstreamT.h"
#include <math.h>
#include <iostream.h>
#include "ifstreamT.h"

using namespace Tahoe;

Ogden::Ogden(ifstreamT& in): fthird(1.0/3.0)
{
	/*read in material parameters*/
	in >> fr;
	falpha_r.Dimension(fr);
	fmu_r.Dimension(fr);
  
	for (int i = 0; i < fr; i++) {
		in >> falpha_r[i];
		in >> fmu_r[i];
	}
	in >> fkappa;

	fmu = 0.0;
	for (int i = 0; i < fr; i++) {
		fmu += falpha_r[i]*fmu_r[i];
	}
	fmu *= 0.5;	
}

void Ogden::Print(ostream& out) const
{
	out<<"      Shear Modulus = "<<fmu<<'\n';
	out<<"      Bulk Modulus = "<<fkappa<<'\n';
}

void Ogden::PrintName(ostream& out) const
{
	out << "Compressible Ogden Potential\n";
}
void Ogden::Initialize(void)
{}

double Ogden::Energy(const dArrayT& lambda_bar, const double& J)
{
	const double& l0 = lambda_bar[0];
	const double& l1 = lambda_bar[1];
	const double& l2 = lambda_bar[2];
 
	double phi = 0.0;

	for (int i = 0; i < fr; i++) {
		const double& mu = fmu_r[i];
		const double& alpha = falpha_r[i];

		double I1 = pow(l0,alpha*0.5) + pow(l1,alpha*0.5) + pow(l2,alpha*0.5);

		phi += mu/alpha*(I1-3);
	}
   
	phi += fkappa*(J*J-1-2*log(J));
	return(phi);
}
void Ogden::DevStress(const dArrayT& lambda_bar,dArrayT& tau)
{  
	const double& l0 = lambda_bar[0];
	const double& l1 = lambda_bar[1];
	const double& l2 = lambda_bar[2];
  
	tau = 0.0;
 
	for (int i = 0; i < fr; i++) {	
		const double& mu = fmu_r[i];
		const double& alpha = falpha_r[i];
		
		tau[0] += mu*fthird*(2.0*pow(l0,alpha*0.5)-pow(l1,alpha*0.5)-pow(l2,alpha*0.5));
		tau[1] += mu*fthird*(2.0*pow(l1,alpha*0.5)-pow(l2,alpha*0.5)-pow(l0,alpha*0.5));
		tau[2] += mu*fthird*(2.0*pow(l2,alpha*0.5)-pow(l0,alpha*0.5)-pow(l1,alpha*0.5));
	}
}

double Ogden::MeanStress(const double& J) {return(0.5*fkappa*(J*J-1));}

void Ogden::DevMod(const dArrayT& lambda_bar, dSymMatrixT& eigenmodulus)
{
	double ninth = fthird*fthird;
  
	const double& l0 = lambda_bar[0];
	const double& l1 = lambda_bar[1];
	const double& l2 = lambda_bar[2];
  
	eigenmodulus = 0.0;
	
	for (int i = 0; i < fr; i++) {	
		const double& mu = fmu_r[i];
		const double& alpha = falpha_r[i];

		eigenmodulus[0] += alpha*mu*ninth*(4.0*pow(l0,alpha*0.5)+pow(l1,alpha*0.5)+pow(l2,alpha*0.5));
		eigenmodulus[1] += alpha*mu*ninth*(4.0*pow(l1,alpha*0.5)+pow(l2,alpha*0.5)+pow(l0,alpha*0.5));
		eigenmodulus[2] += alpha*mu*ninth*(4.0*pow(l2,alpha*0.5)+pow(l0,alpha*0.5)+pow(l1,alpha*0.5));	
		
		eigenmodulus[3] += alpha*mu*ninth*(-2.0*pow(l1,alpha*0.5)-2.0*pow(l2,alpha*0.5)+pow(l0,alpha*0.5));
		eigenmodulus[4] += alpha*mu*ninth*(-2.0*pow(l0,alpha*0.5)-2.0*pow(l2,alpha*0.5)+pow(l1,alpha*0.5));
		eigenmodulus[5] += alpha*mu*ninth*(-2.0*pow(l0,alpha*0.5)-2.0*pow(l1,alpha*0.5)+pow(l2,alpha*0.5));
	}
}

double Ogden::MeanMod(const double& J) {return(fkappa*J*J);}

