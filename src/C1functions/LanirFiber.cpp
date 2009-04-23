/* $Id: LanirFiber.cpp,v 1.1 2009-04-23 14:51:52 thao Exp $ */

#include "LanirFiber.h"
#include <iostream.h>
#include <math.h>
#include "ExceptionT.h"
#include "dArrayT.h"

#include "Gamma.h"


/* constructors */
const double ep_min = 0.001;
using namespace Tahoe;

LanirFiber::LanirFiber(double K, double alpha, double beta): 
	fK(K), 
	falpha(alpha), 
	fbeta(beta) 
{ 
	SetName("lanir_fiber");
}

LanirFiber::LanirFiber(void): 
	fK(0), 
	falpha(0), 
	fbeta(0) 
{ 
	SetName("lanir_fiber");
}
double LanirFiber::PDF(double x) const
{
	double a, b;
	a = 1.0/fbeta;
	b = pow(x*a, falpha-1.0);
	double weibull =  falpha*a * b * exp(-b * x*a);
	
	/*the same for a gamma function*/
	/*
	//gamma function
	Gamma gamma_f;

	//coefficient for gamma distribution
	double coeff = 1.0/(pow(fbeta, falpha)*gamma_f.Function(falpha));
	double gammafun = coeff*pow(x,falpha-1.0)*exp(-x/fbeta);
	*/
	return(weibull);
}

double LanirFiber::Function(double r) const
{
	const char caller[] = "LanirFiber::Function";

	if (r > 1.0)
	{
		/*integrate*/
		double ep = 0.5*(r-1.0);
		double nstep = floor(ep/ep_min)+1;	
		
		double dE = ep/nstep;
		double H_n = 0.0; 
		double H2_n = 0.0; 
		double Psi = 0.0;
		double E = 0.0;

		for (int i = 0; i<nstep; i++)
		{
			/*strain at i+1*/
			E += dE;
			/*midpoint strain*/
			double Emid = E-0.5*dE;
		
			/* Distribution evaluated at midpoint strain*/
			double Pmid = PDF(Emid);
			double temp = Pmid/((1.0+2.0*Emid)*(1.0+2.0*Emid));
		
			/*H = int_0^E 1/(beta^alpha Gamma(alpha)) x^(alpha-1)/(1+2x)^2 dx*/
			double H = H_n + temp*dE;
		
			/*H2 = int_0^E 1/(beta^alpha Gamma(alpha)) x^(alpha-1) x/(1+2x)^2 dx*/
			double H2 = H2_n + temp*dE*Emid;

			/*fiber energy*/
			Psi += fK*dE* (H_n*(E - 0.5*dE) - H2_n + 0.5*temp*(0.25*dE*dE) );

			/*Update */
			H_n = H;
			H2_n = H2;
		}
		return(Psi);
	}
	else return(0.0); 
}

double LanirFiber::DFunction(double r) const
{
	const char caller[] = "LanirFiber::DFunction";

	if (r>1.0)
	{
		/*integrate*/
		double ep = 0.5*(r-1.0);
		double nstep = floor(ep/ep_min)+1;	

		double dE = ep/nstep;
		double E = 0.0;
		double H_n = 0.0; 
		double Sf = 0.0; 
		
		for (int i = 0; i<nstep; i++)
		{
			/*strain at i+1*/
			E += dE;
			/*midpoint strain*/
			double Emid = E-0.5*dE;
		
			/* Distribution evaluated at midpoint strain*/
			double Pmid = PDF(Emid);
			double temp = Pmid/((1.0+2.0*Emid)*(1.0+2.0*Emid));
		
			/*H = int_0^E 1/(beta^alpha Gamma(alpha)) x^(alpha-1)/(1+2x)^2 dx*/
			double H = H_n + temp*dE;
		
			/*fiber dPsi/dlambda*/
			Sf += 0.5*fK*dE*(H_n + 0.5*temp*dE );

	//		cout << "\nvals: "<<Emid<<"\t"<<Pmid<<"\t"<<H<<"\t"<<Sf;
			H_n = H;
		}
		return(Sf);
	}
	else return(0.0);
}

double LanirFiber::DDFunction(double r) const
{
	const char caller[] = "LanirFiber::DFunction";

	if (r>1.0)
	{
		/*integrate*/
		double ep = 0.5*(r-1.0);
		double nstep = floor(ep/ep_min)+1;	

		double dE = ep/nstep;
		double E = 0.0;
		double Cf = 0.0;

		for (int i = 0; i<nstep; i++)
		{
			/*strain at i+1*/
			E +=dE;
	
			/*midpoint strain*/
			double Emid = E-0.5*dE;
		
			/* Distribution evaluated at midpoint strain*/
			double Pmid = PDF(Emid);
			double temp = Pmid/((1.0+2.0*Emid)*(1.0+2.0*Emid));
		
			/*H = int_0^E 1/(beta^alpha Gamma(alpha)) x^(alpha-1)/(1+2x)^2 dx*/
			Cf += 0.25*fK*temp*dE;
		}		
		return(Cf);
	}
	else return(0.0);
}

/* returning values in groups */
dArrayT& LanirFiber::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		*pU++ = Function(r);
	}
	return(out);
}

dArrayT& LanirFiber::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
//		cout <<"\nr: "<<r;
//		cout<<"\ndU: "<<(r-1.0);
		*pdU++ = DFunction(r);
	}
	return(out);
}

dArrayT& LanirFiber::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		*pddU++ = DDFunction(r);
	}
	return(out);
}

void LanirFiber::ParameterDerivatives(double r, dArrayT& data) const
{
	const char caller[] = "LanirFiber::ParameterDerivatives";
	
	if (data.Length() !=3)
		ExceptionT::GeneralFail("length of data array does not equal number of parameters");
		
	if (r > 1.0)
	{
		/*integrate*/
		double ep = 0.5*(r-1.0);
		double nstep = floor(ep/ep_min)+1;	

		double dE = ep/nstep;
		double E = 0.0;
		double H_K_n = 0.0; 
		double H_dalpha_n = 0.0; 
		double H_dbeta_n = 0.0; 

		for (int i = 0; i<nstep; i++)
		{
			/*strain at i+1*/
			E += dE;
			/*midpoint strain*/
			double Emid = E-0.5*dE;
		
			/* Distribution evaluated at midpoint strain*/
			double a, b, c;
			a = 1.0/fbeta;
			b = pow(Emid*a, falpha-1.0);
			c = pow(Emid*a, falpha);

			double Pmid = PDF(Emid);
			double dP_dalpha = Pmid/falpha *(1.0- falpha*log(Emid*a)*(c - 1.0));
			double dP_dbeta = Pmid * falpha*a*(c-1);

			double coeff = 1.0/((1.0+2.0*Emid)*(1.0+2.0*Emid));
		
			/*H = int_0^E 1/(beta^alpha Gamma(alpha)) x^(alpha-1)/(1+2x)^2 dx*/
			double H_K = H_K_n + Pmid*coeff*dE;
			double H_dalpha = H_dalpha_n + dP_dalpha*coeff*dE;
			double H_dbeta = H_dbeta_n + dP_dbeta*coeff*dE;
		
			/*derivatives*/
			data[0] += 0.5*dE*(H_K_n + 0.5*Pmid*dE );
			data[1] += 0.5*fK*dE*(H_dalpha_n + 0.5*coeff*dP_dalpha*dE );
			data[2] += 0.5*fK*dE*(H_dbeta_n + 0.5*coeff*dP_dbeta*dE );
		
			H_K_n = H_K;		
			H_dalpha_n = H_dalpha;
			H_dbeta_n = H_dbeta;
		}
	}
	else data = 0.0;
}
void LanirFiber::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);
	
	ParameterT K(ParameterT::Double, "fiber_stiffness_K");
	K.AddLimit(0,LimitT::Lower);
	
/*	ParameterT alpha(ParameterT::Double, "gamma_shape_param_alpha");
	alpha.AddLimit(1,LimitT::LowerInclusive);

	ParameterT beta(ParameterT::Double, "gamma_scale_param_beta");
	beta.AddLimit(0,LimitT::Lower);

	list.AddParameter(K);
	list.AddParameter(alpha);
	list.AddParameter(beta);
	
	list.SetDescription("f(I) = int_0^y 1/(beta^alpha Gamma(alpha)) x^(alpha-1) exp(-x/beta) fK*((y-x)/(1+2x))^2 dx");	
*/
	ParameterT alpha(ParameterT::Double, "Weibull_shape_param_alpha");
	alpha.AddLimit(1,LimitT::LowerInclusive);

	ParameterT beta(ParameterT::Double, "Weibull_scale_param_beta");
	beta.AddLimit(0,LimitT::Lower);

	list.AddParameter(K);
	list.AddParameter(alpha);
	list.AddParameter(beta);
	
	/* set the description */
	list.SetDescription("f(I) = int_0^y alpha/(beta) (x/beta)^(alpha-1) exp(-(x/beta)^alpha) fK*((y-x)/(1+2x))^2 dx");	
}

void LanirFiber::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);
	fK = list.GetParameter("fiber_stiffness_K");

/*	falpha = list.GetParameter("gamma_shape_param_alpha");
	fbeta = list.GetParameter("gamma_scale_param_beta");

	//gamma function
	Gamma gamma_f;

	//coefficient for gamma distribution
	double coeff = 1.0/(pow(fbeta, falpha)*gamma_f.Function(falpha));
*/
	falpha = list.GetParameter("Weibull_shape_param_alpha");
	fbeta = list.GetParameter("Weibull_scale_param_beta");
	
}

