/* $Id: J2IsoVIB3DLinHardT.cpp,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (10/12/1998)                                          */
/* VIB plus principal stretch elasticity                                  */
/* Interface for a elastoplastic material that is linearly                */
/* isotropically elastic subject to the Huber-von Mises yield             */
/* condition as fYield with kinematic/isotropic hardening laws            */
/* given by:                                                              */
/* 		H(a) = (1 - ftheta) fH_bar a                                         */
/* K(a) = fYield + ftheta fH_bar a                                        */
/* 		where a is the internal hardening variable                           */
/* 	Note: all calculations are peformed in 3D.                            */

#include "J2IsoVIB3DLinHardT.h"

#include <iostream.h>
#include <math.h>

#include "ElasticT.h"
#include "Constants.h"
#include "iArrayT.h"
#include "ElementCardT.h"
#include "StringT.h"
#include "C1FunctionT.h"

/* flags */
const int kNumFlags = 2;
const int kEP   = 0;
const int kInit = 1;
	//not used
	
/* elastic/plastic flag values */
const int kE0        =-2; // elastic but intermediate state not initialized
const int kP0        =-1; // hit plastic but intermediate state not initialized
const int kNotInit   = 0;
const int kIsPlastic = 1;
const int kIsElastic = 2;
const int kReset     = 3; // indicate not to repeat update

/* init flag values */
const int kIsNotInit = 0;
const int kIsInit    = 1;
	
/* class constants */
const int	kNumInternal   = 5;
const int	kalpha         = 0; /* isotropic hardening         */
const int	kstressnorm    = 1; /* norm of the relative stress */
const int	kdgamma        = 2; /* consistency parameter       */
const int	kftrial        = 3; /* yield function value        */
const int	kDetF_tot      = 4; /* determinant of total F      */

const double sqrt23        = sqrt(2.0/3.0);
const double kYieldTol     = 1.0e-10;

const int kNSD = 3;

/* element output data */
const int kNumOutput = 4;
static const char* Labels[kNumOutput] = {"s_max", "s_min", "VM stress", "alpha"};

/* constructor */
J2IsoVIB3DLinHardT::J2IsoVIB3DLinHardT(ifstreamT& in, const ElasticT& element):
	IsoVIB3D(in, element),
	J2PrimitiveT(in),

	fLocLastDisp(element.LastDisplacements()),

//TEMP
	fEigs(kNSD),
	fBeta(kNSD),

	fb_elastic(kNSD),
	fEPModuli(kNSD),

	/* work space */
	fMatrixTemp1(kNSD),
	fMatrixTemp2(kNSD),
	fdev_beta(kNSD),
	
	/* deformation gradient stuff */
	fFtot(3),
	ffrel(3),
	fF_temp(3)
{
	/* check last displacements */
	if (!fLocLastDisp.IsRegistered() ||
		 fLocLastDisp.MinorDim() != NumDOF())
	{
		cout << "\n J2IsoVIB3DLinHardT::J2IsoVIB3DLinHardT: last local displacement vector is invalid" << endl;
		throw eGeneralFail;
	}
}

/* update internal variables */
void J2IsoVIB3DLinHardT::UpdateHistory(void)
{
	/* get flags */
	ElementCardT& element = CurrentElement();
	iArrayT& Flags = element.IntegerData();
	
	/* update plastic variables */
	int num_ip = NumIP();
	for (int ip = 0; ip < num_ip; ip++)
		if (Flags[ip] != kReset && Flags[ip] != kNotInit)
	  //if (Flags[ip] == kIsPlastic || Flags[ip] == kP0)
		{
			/* do not repeat if called again */
			Flags[ip] = kReset; //why was I worried about this????
						
			/* fetch element data */
			LoadData(element,ip);
	
			/* plastic increment */
			fInternal[kalpha] += sqrt23*fInternal[kdgamma];
			
			/* set principal values and directions */
			fSpectral.SpectralDecomp(fb_tr, true);
			
			/* log_e -> eig[b_e] */
			flog_e[0] = exp(2.0*flog_e[0]);
			flog_e[1] = exp(2.0*flog_e[1]);
			flog_e[2] = exp(2.0*flog_e[2]);
					
			/* compute intermediate config */
			fb_n = fSpectral.EigsToRank2(flog_e);
		}
}

/* reset internal variables to last converged solution */
void J2IsoVIB3DLinHardT::ResetHistory(void)
{
	/* flag not to update again */
	//(element->IntegerData()) = kReset;
	//need to check if initialized!
	
	/* get flags */
	ElementCardT& element = CurrentElement();
	iArrayT& Flags = element.IntegerData();

	for (int i = 0; i < Flags.Length(); i++)
		if (Flags[i] == kIsElastic || Flags[i] == kIsPlastic)
			Flags[i] = kReset;
}

/* print parameters */
void J2IsoVIB3DLinHardT::Print(ostream& out) const
{
	/* inherited */
	IsoVIB3D::Print(out);
	J2PrimitiveT::Print(out);
}

/* modulus */
const dMatrixT& J2IsoVIB3DLinHardT::c_ijkl(void)
{
	cout << "\n J2IsoVIB3DLinHardT::c_ijkl: no tangent implemented" << endl;
	throw eGeneralFail;
	
	return fddW; // dummy
}
	
/* stress */
const dSymMatrixT& J2IsoVIB3DLinHardT::s_ij(void)
{
	/* Compute F_total and f_relative 3D */
	ComputeGradients();

	/* get trial stretch */
	int ip = CurrIP();
	const dSymMatrixT& b_tr = TrialStretch(fFtot, ffrel, ip);

	/* spectral decomposition */
	fSpectral.SpectralDecomp(b_tr, true);

	/* principal stretches */
	const dArrayT& b_eigs = fSpectral.Eigenvalues();

	/* det F: won't change since flow is isochoric */
	double J = sqrt(b_eigs.Product());
	
	/* compute trial principal stresses */
	ComputeBeta(b_eigs, fBeta);

	/* plastic correction */
	ReturnMapping(b_tr, b_eigs, fBeta, ip);

	/* Kirchhoff -> Cauchy */	
	fBeta /= J;

	/* build from spectral */
	return fSpectral.EigsToRank2(fBeta);
}

/* material description */
const dMatrixT& J2IsoVIB3DLinHardT::C_IJKL(void)
{
	cout << "\n J2IsoVIB3DLinHardT::C_IJKL: not optimized for total Lagrangian formulation.";
	cout <<   "    use updated Lagrangian formulation." << endl;
	throw eGeneralFail;

	return fddW; // dummy
}

const dSymMatrixT& J2IsoVIB3DLinHardT::S_IJ(void)
{
	cout << "\n J2IsoVIB3DLinHardT::S_IJ: not optimized for total Lagrangian formulation.";
	cout <<   "    use updated Lagrangian formulation." << endl;
	throw eGeneralFail;

	return fb_elastic; // dummy
}

/* strain energy density */
double J2IsoVIB3DLinHardT::StrainEnergyDensity(void)
{
	/* Compute F_total and f_relative 3D */
	ComputeGradients();

	/* principal stretches */
	const dSymMatrixT& b_tr = TrialStretch(fFtot, ffrel, CurrIP());
	b_tr.PrincipalValues(fEigs);

	/* stretched bonds */
	ComputeLengths(fEigs);

	/* update potential table */
	fPotential->MapFunction(fLengths,fU);

	/* sum contributions */
	double  energy = 0.0;
	double* pU = fU.Pointer();
	double* pj = fjacobian.Pointer();	
	for (int i = 0; i < fLengths.Length(); i++)
		energy += (*pU++)*(*pj++);
	
	return energy;
}

/* required parameter flags */
bool J2IsoVIB3DLinHardT::NeedLastDisp(void) const { return true; }

/* returns number of internal variables */
int J2IsoVIB3DLinHardT::NumOutputVariables(void) const { return kNumOutput; }
void J2IsoVIB3DLinHardT::OutputLabels(ArrayT<StringT>& labels) const
{
	/* set size */
	labels.Allocate(kNumOutput);
	
	/* copy labels */
	for (int i = 0; i < kNumOutput; i++)
		labels[i] = Labels[i];
}

void J2IsoVIB3DLinHardT::ComputeOutput(dArrayT& output)
{
	/* compute Cauchy stress (sets fBeta) */
	s_ij();
	
	/* in-plane principal values (Cauchy) */
	if (fBeta[0] > fBeta[1])
	{
		output[0] = fBeta[0];
		output[1] = fBeta[1];
	}
	else
	{
		output[0] = fBeta[1];
		output[1] = fBeta[0];
	}

	/* Cauchy -> Kirchhoff */
	fBeta *= sqrt( fEigs[0]*fEigs[1]*fEigs[2] );
	
	/* Von Mises equivalent (Kirchhoff) stress */
	fBeta    -= fBeta.Average(); //deviatoric part
	output[2] = fBeta.Magnitude()/sqrt23;

	/* plastic evolution parameter */
	if (CurrentElement().IsAllocated())
		output[3] = fInternal[kalpha] + sqrt23*fInternal[kdgamma];
		//element data already loaded in sij(). alpha not set
		//until Update(), which hasn't occurred yet
	else
		output[3] = 0.0;
}

/***********************************************************************
* Protected
***********************************************************************/

/* print name */
void J2IsoVIB3DLinHardT::PrintName(ostream& out) const
{
	/* inherited */
	IsoVIB3D::PrintName(out);

	out << "    J2 plasticity, principal stretch return mapping\n";
}

/* returns the elastic stretch */
const dSymMatrixT& J2IsoVIB3DLinHardT::TrialStretch(const dMatrixT& F_total,
	const dMatrixT& f_relative, int ip)
{
	/* element pointer */
	ElementCardT& element = CurrentElement();
	
	/* compute left Cauchy-Green */
	if (element.IsAllocated()) /* element has been plastic */
	{
		/* load internal variables */
		LoadData(element, ip);
	
		/* check intermediate state */
		iArrayT& Flags = element.IntegerData();
		if (Flags[ip] <= kNotInit)
		{
			/* initialize intermediate config */
			InitIntermediate(F_total, f_relative);
			
			/* reset EP flag */
			Flags[ip] = (Flags[ip] == kP0) ? kIsPlastic : kIsElastic;
		}
	
		/* compute (and store) trial stretch */
		fb_tr.MultQBQT(f_relative, fb_n);

		return fb_tr;
	}
	else /* element is elastic */
	{
		fb_elastic.MultAAT(F_total);
	
		return fb_elastic;
	}
}

/* return the correction to principal stress vector computed by the
* mapping the stress back to the yield surface, if needed, returning
* 1 of a plastic correction is needed (even possible and valid) */
void J2IsoVIB3DLinHardT::ReturnMapping(const dSymMatrixT& b_tr, const dArrayT& b_eigs,
	dArrayT& beta, int ip)
{
	/* deviatoric part of principal stresses */
	fdev_beta = beta;
	fdev_beta -= fdev_beta.Average();

	/* element pointer */
	ElementCardT& element = CurrentElement();

	/* check consistency and initialize plastic element */
	if (PlasticLoading(fdev_beta, element, ip) &&
	    !element.IsAllocated())
	{
		/* new element */
		AllocateElement(element);
							
		/* set element data */
		PlasticLoading(fdev_beta, element, ip);
		
		/* initialize trial stretch */
		fb_tr = b_tr;		
	}

	/* for plastic elements */
	if (element.IsAllocated())
	{
		/* fetch data */
		double  ftrial = fInternal[kftrial];
		double& dgamma = fInternal[kdgamma];
				
		/* store log elastic stretches */
		flog_e[0] = 0.5*log(b_eigs[0]);
		flog_e[1] = 0.5*log(b_eigs[1]);
		flog_e[2] = 0.5*log(b_eigs[2]);
		
		/* is a plastic increment */
		if (ftrial > 0.0)
		{
			/* temp work space */
			dArrayT  res_e(kNSD);    // log stretches residual
			dArrayT  log_e_tr(kNSD); // log stretch of trial elastic state
			dArrayT  eigs_e(kNSD);   // eigenvalues of b_e
			dArrayT  d_beta(kNSD), d_log_e(kNSD);
			dMatrixT ddw_inv(kNSD), h(kNSD), dBdB(kNSD);

			/* set d_Phi/d_beta d_beta */
			dBdB.Outer(fUnitNorm,fUnitNorm);
			dBdB += (1.0/3.0);
			dBdB.PlusIdentity(-1.0);
			dBdB *= (-1.0/fInternal[kstressnorm]);

			/* initialize log elastic stretches */
			log_e_tr = flog_e;
			eigs_e   = b_eigs;

			/* initialize errors */
			dgamma = 0.0;
			res_e  = 0.0; // since dgamma = 0.0
			double mag_res_e = 0.0;

			/* return mapping - set plastic increment */
			int count = 0;
			while ((mag_res_e > kYieldTol || fabs(ftrial) > kYieldTol)
			       && ++count < 20)
			{		
				/* compute modified elastic moduli */
				Computeddw(eigs_e, ddw_inv);
				ddw_inv.Inverse();
				h.SetToCombination(1.0, ddw_inv, dgamma, dBdB);
				h.Inverse();
				
				/* consistency parameter increment */
				double ddgamma = (ftrial + h.MultmBn(fUnitNorm,res_e))/
				                     h.MultmBn(fUnitNorm,fUnitNorm);
	
				/* beta increment */
				res_e.AddScaled(-ddgamma, fUnitNorm);
				h.Multx(res_e, d_beta);

				/* log elastic stretch increment */
				ddw_inv.Multx(d_beta, d_log_e);
	
				/* compute new elastic state */
				dgamma += ddgamma;
				flog_e += d_log_e;
				
				/* elastic principal stresses */
				eigs_e[0] = exp(2.0*flog_e[0]);
				eigs_e[1] = exp(2.0*flog_e[1]);
				eigs_e[2] = exp(2.0*flog_e[2]);
				ComputeBeta(eigs_e, beta);
				
				/* deviatoric principal stresses */	
				fdev_beta  = beta;
				fdev_beta -= fdev_beta.Average();
	
				/* yield condition */		
				ftrial = YieldCondition(fdev_beta, fInternal[kalpha]);
	
				/* log stretch residual */
				res_e.SetToCombination(-1.0,flog_e,
				                        1.0,log_e_tr,
				                    -dgamma,fUnitNorm);
				mag_res_e = res_e.Magnitude();
			}
	
			/* return mapping failed */
			if (count == 20)
			{
				cout << "\n J2IsoVIB2DLinHardT::ReturnMapping: failed to converge after ";
				cout << count << " iterations\n";
				cout << " tolerance         = " << kYieldTol << '\n';
				cout << " consistency error = " << ftrial    << '\n';
				cout << " log stretch error = " << mag_res_e << endl;
				throw eGeneralFail;
			}
		}
	}
}	

/* return a pointer to a new plastic element object constructed with
* the data from element */
void J2IsoVIB3DLinHardT::AllocateElement(ElementCardT& element)
{
	/* number of integration points */
	int num_ip = NumIP();

	/* determine storage */
	int i_size = 0;
	i_size += num_ip; //flags

	int d_size = 0;
	d_size += dSymMatrixT::NumValues(kNSD)*num_ip; // fb_n
	d_size += dSymMatrixT::NumValues(kNSD)*num_ip; // fb_tr
	d_size += kNSD*num_ip;                  // fbeta_tr
	d_size += kNSD*num_ip;                  // flog_e
	d_size += kNSD*num_ip;                  // fUnitNorm - principal space
	d_size += kNumInternal*num_ip;          // fInternal

	/* construct new plastic element */
	element.Allocate(i_size, d_size);

	/* initialize values */
	element.IntegerData() = kNotInit;
	element.DoubleData()  = 0.0;
}

/***********************************************************************
* Private
***********************************************************************/

/* compute F_total and f_relative */
void J2IsoVIB3DLinHardT::ComputeGradients(void)
{
	/* compute relative displacement gradient */
	fFtot = F();
	fF_temp.Inverse(F(fLocLastDisp));
	ffrel.MultAB(fFtot,fF_temp);
}

/* initialize intermediate state from F_n (for ) */
void J2IsoVIB3DLinHardT::InitIntermediate(const dMatrixT& F_total,
	const dMatrixT& f_relative)
{
	/* compute F_n */
	fMatrixTemp1.Inverse(f_relative);
	fMatrixTemp2.MultAB(fMatrixTemp1, F_total);
	
	/* compute (and store) b */
	fb_n.MultAAT(fMatrixTemp2);
}

/* load element data for the specified integration point */
void J2IsoVIB3DLinHardT::LoadData(const ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* decode */
	int stressdim = dSymMatrixT::NumValues(kNSD);
	int blocksize = stressdim + stressdim + kNSD + kNSD + kNSD + kNumInternal;
	int dex       = ip*blocksize;
	
	     fb_n.Set(        kNSD, &d_array[dex             ]);
	    fb_tr.Set(        kNSD, &d_array[dex += stressdim]);
	 fbeta_tr.Set(        kNSD, &d_array[dex += stressdim]);
	   flog_e.Set(        kNSD, &d_array[dex += kNSD     ]);
	fUnitNorm.Set(        kNSD, &d_array[dex += kNSD     ]);
	fInternal.Set(kNumInternal, &d_array[dex += kNSD     ]);     	
}

/* returns 1 if the trial elastic strain state lies outside of the
* yield surface.
*
* NOTE: pass (element = NULL) if element is not yet plastic, ie. has
*       no stored internal variables.
*/
int J2IsoVIB3DLinHardT::PlasticLoading(const dArrayT& dev_beta,
	ElementCardT& element, int ip)
{
	/* not yet plastic */
	if (!element.IsAllocated())
		return (YieldCondition(dev_beta, 0.0) > kYieldTol);
	/* already plastic */
	else
	{
		/* get flags */
		iArrayT& Flags = element.IntegerData();

		/* load internal variables */
		LoadData(element,ip);
		
		/* set internal variables */		
		fInternal[kftrial] = YieldCondition(dev_beta, fInternal[kalpha]);
		
		/* plastic */
		if (fInternal[kftrial] > kYieldTol)
		{
			/* compute unit normal */
			double& norm = fInternal[kstressnorm];

			norm = dev_beta.Magnitude();
			fUnitNorm.SetToScaled(1.0/norm, dev_beta);
		
			/* set flag */
			Flags[ip] = (Flags[ip] == kNotInit) ? kP0 : kIsPlastic;
	
			return 1;
		}
		/* elastic */
		else
		{
			/* set flag */
			Flags[ip] = (Flags[ip] == kNotInit) ? kE0 : kIsElastic;
			
			return 0;
		}
	}
}

double J2IsoVIB3DLinHardT::YieldCondition(const dArrayT& devpstress,
	double alpha) const
{
	return devpstress.Magnitude() - sqrt23*K(alpha);
}

/* integrate Kirchhoff stress principal values */
void J2IsoVIB3DLinHardT::ComputeBeta(const dArrayT& eigs, dArrayT& beta)
{
	/* stretched bonds */
	ComputeLengths(eigs);

	/* derivatives of the potential */
	fPotential->MapDFunction(fLengths, fdU);

	/* initialize kernel pointers */
	double* pdU = fdU.Pointer();
	double* pl  = fLengths.Pointer();
	double* pj  = fjacobian.Pointer();

	double* p1  = fStressTable(0);
	double* p2  = fStressTable(1);
	double* p3  = fStressTable(2);
	
	/* PK2 principal values */	
	double s1 = 0.0;
	double s2 = 0.0;
	double s3 = 0.0;
	for (int i = 0; i < fLengths.Length(); i++)
	{
		double factor = (*pj++)*(*pdU++)/(*pl++);
		
		s1 += factor*(*p1++);
		s2 += factor*(*p2++);
		s3 += factor*(*p3++);
	}

	/* PK2 -> Kirchhoff stress */
	beta[0] = s1*eigs[0];
	beta[1] = s2*eigs[1];
	beta[2] = s3*eigs[2];
}

void J2IsoVIB3DLinHardT::Computeddw(const dArrayT& eigs, dMatrixT& ddw)
{
	/* stretched bonds */
	ComputeLengths(eigs);

	/* derivatives of the potential */
	fPotential->MapDFunction(fLengths, fdU);
	fPotential->MapDDFunction(fLengths, fddU);	

	/* initialize kernel pointers */
	double* pdU  = fdU.Pointer();
	double* pddU = fddU.Pointer();
	double* pl   = fLengths.Pointer();
	double* pj   = fjacobian.Pointer();

	/* stress */
	double *ps1  = fStressTable(0);
	double *ps2  = fStressTable(1);
	double *ps3  = fStressTable(2);

	/* modulus */
	double *pc11  = fModuliTable(0);
	double *pc22  = fModuliTable(1);
	double *pc33  = fModuliTable(2);

	double *pc23  = fModuliTable(3);
	double *pc13  = fModuliTable(4);
	double *pc12  = fModuliTable(5);
	
	/* PK2 principal values */	
	double s1 = 0.0;
	double s2 = 0.0;
	double s3 = 0.0;

	ddw = 0.0;
	double& c11 = ddw(0,0);
	double& c22 = ddw(1,1);
	double& c33 = ddw(2,2);

	double& c23 = ddw(1,2);
	double& c13 = ddw(0,2);
	double& c12 = ddw(0,1);

	for (int i = 0; i < fLengths.Length(); i++)
	{
		double sfactor = (*pj)*(*pdU)/(*pl);
		double cfactor = (*pj++)*((*pddU++)/((*pl)*(*pl)) -
		                          (*pdU++)/((*pl)*(*pl)*(*pl)));	
		pl++;
	
		s1 += sfactor*(*ps1++);
		s2 += sfactor*(*ps2++);
		s3 += sfactor*(*ps3++);

		c11 += cfactor*(*pc11++);
		c22 += cfactor*(*pc22++);
		c33 += cfactor*(*pc33++);

		c23 += cfactor*(*pc23++);
		c13 += cfactor*(*pc13++);
		c12 += cfactor*(*pc12++);
	}
	
	c11 *= eigs[0]*eigs[0];
	c22 *= eigs[1]*eigs[1];
	c33 *= eigs[2]*eigs[2];

	c23 *= eigs[1]*eigs[2];
	c13 *= eigs[0]*eigs[2];
	c12 *= eigs[0]*eigs[1];
	
	c11 += 2.0*s1*eigs[0];
	c22 += 2.0*s2*eigs[1];
	c33 += 2.0*s3*eigs[2];
	
	/* copy symmetric */
	ddw(2,1) = ddw(1,2);
	ddw(2,0) = ddw(0,2);
	ddw(1,0) = ddw(0,1);
}	
