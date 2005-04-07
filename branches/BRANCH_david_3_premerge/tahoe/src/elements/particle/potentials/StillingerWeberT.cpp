/* $Id: StillingerWeberT.cpp,v 1.6 2004-12-09 01:41:17 cjkimme Exp $ */
#include "StillingerWeberT.h"
#include <iostream.h>
#include <math.h>
#include "dMatrixT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* initialize static parameters */
double StillingerWeberT::s_eps = 1.0;
double StillingerWeberT::s_sigma = 1.0;
double StillingerWeberT::s_b =-1.0;
double StillingerWeberT::s_A = 0.0;
double StillingerWeberT::s_B = 0.0;
double StillingerWeberT::s_p = 0.0;
double StillingerWeberT::s_q = 0.0;
double StillingerWeberT::s_lambda = 0.0;
double StillingerWeberT::s_gamma = 0.0;
double StillingerWeberT::s_costheta_ideal = 0.0;
bool StillingerWeberT::use_pow = false;

/* constructor */
StillingerWeberT::StillingerWeberT(void):
	f_eps(0.0),
	f_sigma(0.0),
	f_b(0.0),
	f_A(0.0),
	f_B(0.0),
	f_p(0.0),
	f_q(0.0),
	f_lambda(0.0),
	f_gamma(0.0),
	f_costheta_ideal(0.0)
{
	SetName("Stillinger_Weber");
}

/* return a pointer to the energy function */
PairPropertyT::EnergyFunction StillingerWeberT::getEnergyFunction(void)
{
	/* copy my data to static */
	s_eps = f_eps;
	s_sigma = f_sigma;
	s_b = f_b;
	s_A = f_A*f_eps; // compute energy scale now
	s_B = f_B;
	s_p = f_p;
	s_q = f_q;

	/* return function pointer */
	return StillingerWeberT::TwoBodyEnergy;
}

PairPropertyT::ForceFunction StillingerWeberT::getForceFunction(void)
{
	/* copy my data to static */
	s_eps = f_eps;
	s_sigma = f_sigma;
	s_b = f_b;
	s_A = f_A*f_eps; // compute energy scale now
	s_B = f_B;
	s_p = f_p;
	s_q = f_q;

	/* return function pointer */
	return StillingerWeberT::TwoBodyForce;
}

PairPropertyT::StiffnessFunction StillingerWeberT::getStiffnessFunction(void)
{
	/* copy my data to static */
	s_eps = f_eps;
	s_sigma = f_sigma;
	s_b = f_b;
	s_A = f_A*f_eps; // compute energy scale now
	s_B = f_B;
	s_p = f_p;
	s_q = f_q;

	/* return function pointer */
	return StillingerWeberT::TwoBodyStiffness;
}

/* return a pointer to the energy function */
ThreeBodyPropertyT::EnergyFunction StillingerWeberT::getThreeBodyEnergyFunction(void)
{
	/* copy my data to static */
	s_eps = f_eps;
	s_sigma = f_sigma;
	s_b = f_b;
	s_lambda = f_lambda*f_eps; // compute energy scale now
	s_gamma = f_gamma;
	s_costheta_ideal = f_costheta_ideal;

	/* return function pointer */
	return StillingerWeberT::ThreeBodyEnergy;
}

ThreeBodyPropertyT::ForceFunction StillingerWeberT::getThreeBodyForceFunction(void)
{
	/* copy my data to static */
	s_eps = f_eps;
	s_sigma = f_sigma;
	s_b = f_b;
	s_lambda = f_lambda*f_eps; // compute energy scale now
	s_gamma = f_gamma;
	s_costheta_ideal = f_costheta_ideal;

	/* return function pointer */
	return StillingerWeberT::ThreeBodyForce;
}

ThreeBodyPropertyT::StiffnessFunction StillingerWeberT::getThreeBodyStiffnessFunction(void)
{
	/* copy my data to static */
	s_eps = f_eps;
	s_sigma = f_sigma;
	s_b = f_b;
	s_lambda = f_lambda*f_eps; // compute energy scale now
	s_gamma = f_gamma;
	s_costheta_ideal = f_costheta_ideal;

	/* return function pointer */
	return StillingerWeberT::ThreeBodyStiffness;
}

/* describe the parameters needed by the interface */
void StillingerWeberT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	ParameterT mass(fMass, "mass");
	mass.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(mass, ParameterListT::ZeroOrOnce);

	ParameterT eps(f_eps, "energy_scale");
	eps.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(eps, ParameterListT::ZeroOrOnce);

	ParameterT sigma(f_sigma, "length_scale");
	sigma.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(sigma, ParameterListT::ZeroOrOnce);

	ParameterT b(f_b, "cut_off_distance");
	b.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(b, ParameterListT::ZeroOrOnce);

	ParameterT A(f_A, "two_body_energy_scale");
	A.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(A, ParameterListT::ZeroOrOnce);

	ParameterT B(f_B, "repulsive_coefficient");
	B.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(B, ParameterListT::ZeroOrOnce);

	ParameterT p(f_p, "repulsive_exponent");
	p.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(p, ParameterListT::ZeroOrOnce);
	
	ParameterT q(f_q, "attractive_exponent");
	q.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(q, ParameterListT::ZeroOrOnce);
	
	ParameterT lambda(f_lambda, "three_body_energy_scale");
	lambda.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(lambda, ParameterListT::ZeroOrOnce);
	
	ParameterT gamma(f_gamma, "three_body_exp_param");
	gamma.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(gamma, ParameterListT::ZeroOrOnce);
	
	ParameterT costheta(f_costheta_ideal, "cos_theta");
	costheta.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(costheta, ParameterListT::ZeroOrOnce);
}

/* accept parameter list */
void StillingerWeberT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* all parameters default to Si unless specified */
	const ParameterT* opt_param = list.Parameter("mass");
	if (opt_param)
		fMass = *opt_param;
	else
		fMass = 28.09;
	opt_param = list.Parameter("energy_scale");
	if (opt_param)
		f_eps = *opt_param;
	else
		f_eps = 209.2; // eV
	opt_param = list.Parameter("length_scale");
	if (opt_param)
		f_sigma = *opt_param;
	else
		f_sigma = 2.0951; // Angstroms
	opt_param = list.Parameter("cut_off_distance");
	if (opt_param) 
		f_b = *opt_param;
	else
		f_b = 1.8; // units of sigma
	opt_param = list.Parameter("two_body_energy_scale");
	if (opt_param)
		f_A = *opt_param;
	else
		f_A = 7.049556277; // pure number
	opt_param = list.Parameter("repulsive_coefficient");
	if (opt_param)
		f_B = *opt_param;
	else
		f_B = 0.6022245584; // pure number
	opt_param = list.Parameter("repulsive_exponent");
	if (opt_param) {
		f_p = *opt_param;
		if (f_p != 4.0)
			use_pow = true; 
	} else
		f_p = 4.0; // pure number
	opt_param = list.Parameter("attractive_exponent");
	if (opt_param) {
		f_q = *opt_param;
		if (f_q != 0.0)
			use_pow = true;
	} else
		f_q = 0.0; // pure number
	opt_param = list.Parameter("three_body_energy_scale");
	if (opt_param) 
		f_lambda = *opt_param;
	else
		f_lambda = 21.0; // pure number
	opt_param = list.Parameter("three_body_exp_param");
	if (opt_param)
		f_gamma = *opt_param;
	else
		f_gamma = 1.20; // pure number
	opt_param = list.Parameter("cos_theta");
	if (opt_param)
		f_costheta_ideal = *opt_param;
	else
		f_costheta_ideal = 1./3.; // pure number	
	
	SetRange(f_sigma*f_b);
	SetNearestNeighbor(pow(2.0,1.0/6.0)*f_sigma*f_b);	

}

/***********************************************************************
 * Private
 ***********************************************************************/

double StillingerWeberT::TwoBodyEnergy(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)
	double r_c = r_ab/s_sigma;
	if (r_c > s_b)
		return 0.0;
	else {
		if (use_pow) {
			return 0.; // later
		} else {
			
			double r = 1./r_c;
			double r_4 = r*r*r*r*s_B;
			r = exp(1./(r_c - s_b));	
			
			return s_A*(r_4-1)*r;
		}
	}
}

double StillingerWeberT::TwoBodyForce(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

	double r_c = r_ab/s_sigma;
	if (r_c > s_b)
		return 0.0;
	else {
		if (use_pow) {
			return 0.; // later
		} else {
	  		
	  		double r = 1./r_c;
	  		double rms = 1./(r_c - s_b);
	  		double r_4 = r*r*r*r*s_B;
	  		double r_5 = -4.0*r*r_4;
	  		r = exp(rms);
		
			return s_A/s_sigma*r*(r_5 - (r_4 - 1.)*rms*rms);
	  	}
	}
}

double StillingerWeberT::TwoBodyStiffness(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

	double r_c = r_ab/s_sigma;
	if (r_c > s_b)
		return 0.0;
	else {
	  	if (use_pow) {
	  		return 0.0; // later
	  	} else {
	  				
	  		double r = 1./r_c;
	  		double r_4 = r*r*r*r*s_B;
	  		double r_5 = -4.0*r*r_4;
	  		double r_6 = -5.0*r*r_5;
	  		double rms = 1./(r_c - s_b);
	  		r = exp(rms);
	  		r_4 -= 1.;
	  		r_c = rms*rms;
	  			
	  		return  s_A/s_sigma/s_sigma*r*(r_6 + r_c*(r_4*r_c -  2.*(r_5 - r_4*rms)));
	  	} 
	}
}

double StillingerWeberT::ThreeBodyEnergy(const double* ri, const double* rj, const double* rk)
{
	//compute rij, rik, thetaijk;
	double r_ij[3], r_ik[3], rij, rik, costheta;
	rij = rik = costheta = 0.;
	for (int i = 0; i < 3; i++) {
		r_ij[i] = ri[i] - rj[i];
		r_ik[i] = ri[i] - rk[i];
		rij += r_ij[i]*r_ij[i];
		rik += r_ik[i]*r_ik[i];
		costheta += r_ij[i]*r_ik[i];
	}	
	rij = sqrt(rij);
	rik = sqrt(rik);
	costheta /= rij*rik;
	rij /= s_sigma;
	rik /= s_sigma;
	
	// compare with cutoff
	if (rij > s_b || rik > s_b)
		return 0.0;
	else {
		costheta += s_costheta_ideal;
		costheta *= costheta;
		
		rij = 1./(rij - s_b);
		rik = 1./(rik - s_b);
		rij = exp(s_gamma*(rij + rik));
		
		return s_lambda*rij*costheta;
	}
}

double* StillingerWeberT::ThreeBodyForce(const double* ri, const double* rj, const double* rk, 
										double *fij, double *fik)
{
	//compute rij, rik, cos_theta_jik;
	double r_ij[3], r_ik[3], rij, rik, costheta;
	rij = rik = costheta = 0.;
	for (int i = 0; i < 3; i++) {
		r_ij[i] = ri[i] - rj[i];
		r_ik[i] = ri[i] - rk[i];
		rij += r_ij[i]*r_ij[i];
		rik += r_ik[i]*r_ik[i];
		costheta += r_ij[i]*r_ik[i];
	}	
	rij = sqrt(rij);
	rik = sqrt(rik);
	costheta /= rij*rik;
	rij /= s_sigma;
	rik /= s_sigma;

	// compare with cutoff
	if (rij > s_b || rik > s_b) {
		return NULL;
	} else {	
		double exp_ij, exp_ik;
		exp_ij = 1./(rij - s_b);
		exp_ik = 1./(rik - s_b);
		double com_fact = s_lambda*exp(s_gamma*(exp_ij + exp_ik));
		double angle_stretch = costheta + s_costheta_ideal;
		com_fact *= angle_stretch;
		exp_ij *= s_gamma*exp_ij;		
		exp_ik *= s_gamma*exp_ik;
		
		rij *= s_sigma;
		rik *= s_sigma;
		
		for (int i = 0; i < 3; i++) {
			r_ij[i] *= com_fact/rij; // r_ij^hat \left(\cos \theta + \frac{1}{3}\right)\lambda\exp{}
			r_ik[i] *= com_fact/rik;
			fij[i] = -r_ij[i]*angle_stretch*exp_ij/s_sigma;
			fik[i] = -r_ik[i]*angle_stretch*exp_ik/s_sigma;
			fij[i] += 2.*(r_ik[i] - r_ij[i]*costheta)/rij;
			fik[i] += 2.*(r_ij[i] - r_ik[i]*costheta)/rik;
		}
		
		return fij;
	}
	
}

double* StillingerWeberT::ThreeBodyStiffness(const double* ri, const double* rj, 
											const double* rk, dMatrixT& K_ijk)
{
	//compute rij, rik, cos_theta_jik;
	double r_ij[3], r_ik[3], rij, rik, costheta;
	rij = rik = costheta = 0.;
	for (int i = 0; i < 3; i++) {
		r_ij[i] = ri[i] - rj[i];
		r_ik[i] = ri[i] - rk[i];
		rij += r_ij[i]*r_ij[i];
		rik += r_ik[i]*r_ik[i];
		costheta += r_ij[i]*r_ik[i];
	}	
	rij = sqrt(rij);
	rik = sqrt(rik);
	costheta /= rij*rik;
	rij /= s_sigma;
	rik /= s_sigma;
	
	// compare with cutoff
	if (rij > s_b || rik > s_b) {
		return NULL;
	} else {
		double exp_ij, exp_ik;
		exp_ij = 1./(rij - s_b);
		exp_ik = 1./(rik - s_b);
		double com_fact = s_lambda*exp(s_gamma*(exp_ij + exp_ik));	
		double angle_stretch = costheta + s_costheta_ideal;
		double angle_stretch_2 = angle_stretch*angle_stretch;
		exp_ij *= s_gamma*exp_ij/s_sigma;		
		exp_ik *= s_gamma*exp_ik/s_sigma;
		
		rij *= s_sigma;
		rik *= s_sigma;
		for (int i = 0; i < 3; i++) {
			r_ij[i] /= rij;
			r_ik[i] /= rik;
		}	
		
		/* K_jk entry proportional to r_ik x r_ij */
		/*double term1, term2, term3, term4, term5, term6;
		term1 = angle_stretch_2*exp_ij*exp_ik;
		term2 = 2.*angle_stretch*costheta*(exp_ij)/rik;
		term3 = 2.*angle_stretch*costheta*(1./rij)/rik;
		term4 = 2.*(costheta*costheta)/rij/rik;
		term5 = 2.*(1/rik)/rij;
		term6 = 2.*(exp_ik*angle_stretch*costheta)/rij;*/
		double jk_fac_kj = angle_stretch_2*exp_ij*exp_ik;
		jk_fac_kj += 2.*angle_stretch*costheta*(exp_ij/rik + exp_ik/rij + 1./(rij*rik));
		jk_fac_kj += 2.*(costheta*costheta)/rij/rik;
		
		/* K_jk entry proportional to r_ij x r_ik */
		double jk_fac_jk = 2./(rik*rij);
		
		/* K_jk entry proportional to r_ij x r_ij */
		double jk_fac_jj = -2./rik*(angle_stretch*exp_ij + angle_stretch/rij + costheta/rij);
		
		/* K_jk entry proportional to r_ik x r_ik */
		double jk_fac_kk = -2./rij*(angle_stretch*exp_ik + angle_stretch/rik + costheta/rik);
		
		/* K_jk entry proportional to 1 */
		double jk_fac_diag = 2.*angle_stretch/rij/rik;
		
		/* K_jj entry proportional to r_ij x r_ik */
		double jj_fac_jk = -2./rij*(angle_stretch*exp_ij + angle_stretch/rij + costheta/rij);

		/* K_jj entry proportional to r_ij x r_ij */
		double jj_fac_jj = angle_stretch_2*exp_ij*(1./rij+exp_ij) + angle_stretch*costheta/rij*(4.*exp_ij + 6./rij);
		jj_fac_jj += 2.*costheta*costheta/rij/rij + 2.*angle_stretch_2*exp_ij/s_sigma/(rij/s_sigma-s_b);
		
		/* K_jj entry proportional to r_ik x r_ik */
		double jj_fac_kk = 2./rij/rij;
		
		/* K_jj entry proportional to 1 */
		double jj_fac_diag = -(angle_stretch_2*exp_ij + 2.*angle_stretch*costheta/rij)/rij;
		
		/* K_kk entry proportional to r_ij x r_ik */
		double kk_fac_jk = -2./rik*(angle_stretch*exp_ik + angle_stretch/rik + costheta/rik); 
		
		/* K_kk entry proportional to r_ik x r_ik */
		double kk_fac_kk =  angle_stretch_2*exp_ik*(1./rik+exp_ik) + angle_stretch*costheta/rik*(4.*exp_ik + 6./rik);
		kk_fac_kk += 2.*costheta*costheta/rik/rik + 2.*angle_stretch_2*exp_ik/s_sigma/(rik/s_sigma-s_b);
		
		/* K_kk entry proportional to r_ij x r_ij */
		double kk_fac_jj = 2./rik/rik;
		
		/* K_kk entry proportional to 1 */
		double kk_fac_diag = -(angle_stretch_2*exp_ik + 2.*angle_stretch*costheta/rik)/rik;;
		
		dMatrixT d_rij_rij(3,3), d_rik_rik(3,3), d_rij_rik(3,3), tempSpace(3,3);
		dArrayT v1, v2;
		v1.Set(3,r_ij);
		v2.Set(3,r_ij);
		d_rij_rik.Outer(v1,v2,jk_fac_jj);
		v2.Set(3,r_ik);
		d_rij_rik.Outer(v1,v2,jk_fac_kj,dMatrixT::kAccumulate);
		d_rij_rik.Outer(v2,v1,jk_fac_jk,dMatrixT::kAccumulate);
		v1.Set(3,r_ik);
		d_rij_rik.Outer(v1,v2,jk_fac_kk,dMatrixT::kAccumulate);
		d_rij_rik.PlusIdentity(jk_fac_diag);
		
		v1.Set(3,r_ij);
		v2.Set(3,r_ij);
		d_rij_rij.Outer(v1,v2,jj_fac_jj);
		v2.Set(3,r_ik);
		d_rij_rij.Outer(v1,v2,jj_fac_jk,dMatrixT::kAccumulate);
		d_rij_rij.Outer(v2,v1,jj_fac_jk,dMatrixT::kAccumulate);
		v1.Set(3,r_ik);
		d_rij_rij.Outer(v1,v2,jj_fac_kk,dMatrixT::kAccumulate);
		d_rij_rij.PlusIdentity(jj_fac_diag);
		
		v1.Set(3,r_ik);
		v2.Set(3,r_ik);
		d_rik_rik.Outer(v1,v2,kk_fac_kk);
		v2.Set(3,r_ij);
		d_rik_rik.Outer(v1,v2,kk_fac_jk,dMatrixT::kAccumulate);
		d_rik_rik.Outer(v2,v1,kk_fac_jk,dMatrixT::kAccumulate);
		v1.Set(3,r_ij);
		d_rik_rik.Outer(v1,v2,kk_fac_jj,dMatrixT::kAccumulate);
		d_rik_rik.PlusIdentity(kk_fac_diag);
		
		/* final part of chain rule d r_kl/d r_{i,j, or k} */
		tempSpace = d_rij_rik;
		tempSpace.Transpose(d_rij_rik,dMatrixT::kAccumulate);
		tempSpace += d_rij_rij;
		tempSpace += d_rik_rik;
		K_ijk.SetBlock(0,0,tempSpace); // d^2 by d_r_i^2
		tempSpace.Transpose(d_rij_rik);
		tempSpace += d_rij_rij;
		tempSpace *= -1.;
		K_ijk.SetBlock(0,3,tempSpace); // d^2 by d_r_i d_r_j
		tempSpace.Transpose();
		K_ijk.SetBlock(3,0,tempSpace);
		tempSpace.Transpose(d_rij_rik);
		tempSpace += d_rik_rik;
		tempSpace *= -1.;
		K_ijk.SetBlock(6,0,tempSpace); // d^2 by d_r_i d_r_k
		tempSpace.Transpose();
		K_ijk.SetBlock(0,6,tempSpace);
		K_ijk.SetBlock(3,3,d_rij_rij); // d^2 by d_r_j^2
		K_ijk.SetBlock(3,6,d_rij_rik); // d^2 by d_r_j d_r_k
		K_ijk.SetBlock(6,6,d_rik_rik); // d^2 by d_r_k^2
		
		K_ijk *= com_fact;
		
		return K_ijk.Pointer();
	}	
}
