#include "WLC.h"
#include "FSMatSupportT.h"
#include "ifstreamT.h"
#include "ExceptionT.h"
#include <math.h>
//#include <iostream.h>
#include <stdlib.h>

using namespace Tahoe;
 
static const double fk = 1.3806503e-23; 
static const int fn = 4;
/***********************************************************************
 * Public
 ***********************************************************************/

/* constructors */
WLC::WLC(void):
	ParameterInterfaceT("Bischoff-Arruda_WLC"){}

#if 0
/* print parameters */
void WLC::Print(ostream& out) const
{
    out << "\n Anisotropic Bischoff-Aruda model with WLC statistics";
	out << "\n Material Parameters: ";
	out << "\n N: " << fN << "\t k: " << fk << "\t T: "<<  fT <<"\t A: "<< fA;
			
	/*unit cell parameters*/
	out << "\n Unit Cell :";
	out << "\n l1: " << fl1
		<< "\n\t " << f_e1;
	out << "\n l2: " << fl2
		<< "\n\t " << f_e2;
	out << "\n l3: " << fl3
		<< "\n\t " << f_e3;	
		
	/*bulk response parameters*/
	out << "\n gamma: "<< fgamma
		<< "\t beta: " << fbeta;
}
#endif

double WLC::StrainEnergyDensity(void)
{
	/*calculates equilibrium part*/
	Compute_C(fC);
	if (NumSD() == 2)
	{
		fStretch[0] = fC[0];
		fStretch[1] = fC[1];
		fStretch[2] = 1.0;
       
		fStretch[3] = 0.0;
		fStretch[4] = 0.0;
		fStretch[5] = fC[2];
	}
	else fStretch = fC;
	
	double energy = 0.0;
	
	/*entropy part*/
	Compute_R(fR, fStretch);
	const double* pr = fR.Pointer();
	for (int i = 0; i < fn; i++) {
		energy += 0.5*(*pr)*(*pr)/fL - 0.25*(*pr) + 0.25*fL/(1-(*pr)/fL);
		pr++; 
	}
	energy *= 0.25*fN*fk*fT/fA;
	
	/*repulsion part*/
	Compute_eigs(fEigs, fStretch);
	double p1 = sqrt(2.0*fA/fL);
	double p2 = -0.25*fN*fk*fT/(fA*fA*sqrt(2.0*fL/fA)) * (p1 - 0.25 + 0.25/(1-p1));

	energy -= 0.5*p2*(fl1*fl1*log(fEigs[0]) + fl2*fl2*log(fEigs[1]) + fl3*fl3*log(fEigs[3]));

	/*bulk part*/
	double I3 = fStretch.Det();
	double I1 = fStretch.Trace();
	energy += fgamma/fbeta*(pow(I3, -fbeta) - 1.0) + fgamma*(I1 -3.0);
	
	return(energy);
}

/* modulus */
const dMatrixT& WLC::C_IJKL(void)
{
    /* stretch tensor */
    Compute_C(fC);
    if (NumSD() == 2)
    {
      fStretch[0] = fC[0];
      fStretch[1] = fC[1];
      fStretch[2] = 1.0;
    
      fStretch[3] = 0.0;
      fStretch[4] = 0.0;
      fStretch[5] = fC[2];
    }
	else fStretch = fC; 
	
//  	cout <<setprecision(16)<< "\nfStretch: "<<fStretch<<endl;
	/*initializing*/
	fModulus3D = 0.0; 
	
	/*entropic part*/
	Compute_R(fR, fStretch);
	double p0 = 0.0625*fN*fk*fT/fA*fR0*fR0*fR0*fR0;
	
	for (int i = 0; i < fn; i++) {
		const double& r = fR[i]; 
		double p1 = (fL-r);
		double p4 = p0/r*(3*fL-r)/(p1*p1*p1); 
		const dSymMatrixT& M = fM[i];
	
		fModulus3D(0,0) += p4*M[0]*M[0];
		fModulus3D(0,1) += p4*M[0]*M[1];
		fModulus3D(0,2) += p4*M[0]*M[2];
		fModulus3D(0,3) += p4*M[0]*M[3];
		fModulus3D(0,4) += p4*M[0]*M[4];
		fModulus3D(0,5) += p4*M[0]*M[5];

		fModulus3D(1,1) += p4*M[1]*M[1];
		fModulus3D(1,2) += p4*M[1]*M[2];
		fModulus3D(1,3) += p4*M[1]*M[3];
		fModulus3D(1,4) += p4*M[1]*M[4];
		fModulus3D(1,5) += p4*M[1]*M[5];

		fModulus3D(2,2) += p4*M[2]*M[2];
		fModulus3D(2,3) += p4*M[2]*M[3];
		fModulus3D(2,4) += p4*M[2]*M[4];
		fModulus3D(2,5) += p4*M[2]*M[5];

		fModulus3D(3,3) += p4*M[3]*M[3];
		fModulus3D(3,4) += p4*M[3]*M[4];
		fModulus3D(3,5) += p4*M[3]*M[5];

		fModulus3D(4,4) += p4*M[4]*M[4];
		fModulus3D(4,5) += p4*M[4]*M[5];

		fModulus3D(5,5) += p4*M[5]*M[5];
	}

	/*repulsion part*/
	Compute_eigs(fEigs, fStretch);
	double p1 = sqrt(2.0*fA/fL);
	double p2 = 0.5*fN*fk*fT/(fA*fA*sqrt(2.0*fL/fA)) * (p1 - 0.25 + 0.25/((1-p1)*(1-p1)));

	dArrayT coeff(3);
	double p5 = p2 * fl1*fl1/(fEigs[0]*fEigs[0]);
	double p6 = p2 * fl2*fl2/(fEigs[1]*fEigs[1]);
	double p7 = p2 * fl3*fl3/(fEigs[2]*fEigs[2]);

	const dSymMatrixT& M5 = fM[4];
	const dSymMatrixT& M6 = fM[5];
	const dSymMatrixT& M7 = fM[6];
	
	fModulus3D(0,0) += p5*M5[0]*M5[0] + p6*M6[0]*M6[0] + p7*M7[0]*M7[0];
	fModulus3D(0,1) += p5*M5[0]*M5[1] + p6*M6[0]*M6[1] + p7*M7[0]*M7[1];
	fModulus3D(0,2) += p5*M5[0]*M5[2] + p6*M6[0]*M6[2] + p7*M7[0]*M7[2];
	fModulus3D(0,3) += p5*M5[0]*M5[3] + p6*M6[0]*M6[3] + p7*M7[0]*M7[3];
	fModulus3D(0,4) += p5*M5[0]*M5[4] + p6*M6[0]*M6[4] + p7*M7[0]*M7[4];
	fModulus3D(0,5) += p5*M5[0]*M5[5] + p6*M6[0]*M6[5] + p7*M7[0]*M7[5];

	fModulus3D(1,1) += p5*M5[1]*M5[1] + p6*M6[1]*M6[1] + p7*M7[1]*M7[1];
	fModulus3D(1,2) += p5*M5[1]*M5[2] + p6*M6[1]*M6[2] + p7*M7[1]*M7[2];
	fModulus3D(1,3) += p5*M5[1]*M5[3] + p6*M6[1]*M6[3] + p7*M7[1]*M7[3];
	fModulus3D(1,4) += p5*M5[1]*M5[4] + p6*M6[1]*M6[4] + p7*M7[1]*M7[4];
	fModulus3D(1,5) += p5*M5[1]*M5[5] + p6*M6[1]*M6[5] + p7*M7[1]*M7[5];
	
	fModulus3D(2,2) += p5*M5[2]*M5[2] + p6*M6[2]*M6[2] + p7*M7[2]*M7[2];
	fModulus3D(2,3) += p5*M5[2]*M5[3] + p6*M6[2]*M6[3] + p7*M7[2]*M7[3];
	fModulus3D(2,4) += p5*M5[2]*M5[4] + p6*M6[2]*M6[4] + p7*M7[2]*M7[4];
	fModulus3D(2,5) += p5*M5[2]*M5[5] + p6*M6[2]*M6[5] + p7*M7[2]*M7[5];

	fModulus3D(3,3) += p5*M5[3]*M5[3] + p6*M6[3]*M6[3] + p7*M7[3]*M7[3];
	fModulus3D(3,4) += p5*M5[3]*M5[4] + p6*M6[3]*M6[4] + p7*M7[3]*M7[4];
	fModulus3D(3,5) += p5*M5[3]*M5[5] + p6*M6[3]*M6[5] + p7*M7[3]*M7[5];

	fModulus3D(4,4) += p5*M5[4]*M5[4] + p6*M6[4]*M6[4] + p7*M7[4]*M7[4];
	fModulus3D(4,5) += p5*M5[4]*M5[5] + p6*M6[4]*M6[5] + p7*M7[4]*M7[5];

	fModulus3D(5,5) += p5*M5[5]*M5[5] + p6*M6[5]*M6[5] + p7*M7[5]*M7[5]; 

	/*bulk part*/
    double I3 = fStretch.Det();
	dSymMatrixT& iC = fStretch.Inverse();
	double p3 = pow(I3, -fbeta)*4.0*fgamma;

	fModulus3D(0,0) += p3*(1.0+fbeta)*iC[0]*iC[0];
	fModulus3D(0,1) += p3*(iC[5]*iC[5] + fbeta*iC[0]*iC[1]);
	fModulus3D(0,2) += p3*(iC[4]*iC[4] + fbeta*iC[0]*iC[2]);
	fModulus3D(0,3) += p3*(iC[5]*iC[4] + fbeta*iC[0]*iC[3]);
	fModulus3D(0,4) += p3*(1.0+fbeta)*iC[0]*iC[4];
	fModulus3D(0,5) += p3*(1.0+fbeta)*iC[0]*iC[5];

	fModulus3D(1,1) += p3*(1.0+fbeta)*iC[1]*iC[1];
	fModulus3D(1,2) += p3*(iC[3]*iC[3] + fbeta*iC[1]*iC[2]);
	fModulus3D(1,3) += p3*(1.0+fbeta)*iC[1]*iC[3];
	fModulus3D(1,4) += p3*(iC[5]*iC[3] + fbeta*iC[1]*iC[4]);
	fModulus3D(1,5) += p3*(1.0+fbeta)*iC[1]*iC[5];

	fModulus3D(2,2) += p3*(1.0+fbeta)*iC[2]*iC[2];
	fModulus3D(2,3) += p3*(1.0+fbeta)*iC[2]*iC[3];
	fModulus3D(2,4) += p3*(1.0+fbeta)*iC[2]*iC[4];
	fModulus3D(2,5) += p3*(iC[3]*iC[4] + fbeta*iC[2]*iC[5]);
	
	fModulus3D(3,3) += p3*(0.5*iC[1]*iC[2] + (0.5+fbeta)*iC[3]*iC[3]);
	fModulus3D(3,4) += p3*(0.5*iC[5]*iC[2] + (0.5*fbeta)*iC[3]*iC[4]);
	fModulus3D(3,5) += p3*(0.5*iC[1]*iC[4] + (0.5*fbeta)*iC[3]*iC[5]);

	fModulus3D(4,4) += p3*(0.5*iC[0]*iC[2] + (0.5*fbeta)*iC[4]*iC[4]);
	fModulus3D(4,5) += p3*(0.5*iC[0]*iC[3] + (0.5*fbeta)*iC[5]*iC[4]);
	
	fModulus3D(5,5) += p3*(0.5*iC[0]*iC[1] + (0.5*fbeta)*iC[5]*iC[5]);
	
	/*copy lower half from symmetry*/
	
	fModulus3D(1,0) = fModulus3D(0,1); fModulus3D(2,0) = fModulus3D(0,2); fModulus3D(3,0) = fModulus3D(0,3); 
	fModulus3D(4,0) = fModulus3D(4,0); fModulus3D(5,0) = fModulus3D(0,5);
	fModulus3D(2,1) = fModulus3D(1,2); fModulus3D(3,1) = fModulus3D(1,3); fModulus3D(4,1) = fModulus3D(1,4);
	fModulus3D(5,1) = fModulus3D(1,5);
	fModulus3D(3,2) = fModulus3D(2,3); fModulus3D(4,2) = fModulus3D(2,4); fModulus3D(5,2) = fModulus3D(2,5);
	fModulus3D(4,3) = fModulus3D(3,4); fModulus3D(5,3) = fModulus3D(3,5);
	fModulus3D(5,4) = fModulus3D(4,5); 
	
    if (NumSD() == 2)
    {
      fModulus[0] = fModulus3D[0];
      fModulus[1] = fModulus3D[1];
      fModulus[2] = fModulus3D[5];

      fModulus[3] = fModulus3D[6];
      fModulus[4] = fModulus3D[7];
      fModulus[5] = fModulus3D[11];

      fModulus[6] = fModulus3D[30];
      fModulus[7] = fModulus3D[31];
      fModulus[8] = fModulus3D[35];
    }
    else fModulus = fModulus3D; 
//	cout<< "\nfModulus3D: "<<fModulus3D<<endl;
    return fModulus;
}

/* stresses */
const dSymMatrixT& WLC::S_IJ(void)
{

    /* stretch tensor */
    Compute_C(fC);
    if (NumSD() == 2)
    {
      fStretch[0] = fC[0];
      fStretch[1] = fC[1];
      fStretch[2] = 1.0;
    
      fStretch[3] = 0.0;
      fStretch[4] = 0.0;
      fStretch[5] = fC[2];
    }
	else fStretch = fC;
	fStress3D = 0.0;
	
	/*entropic part*/
	Compute_R(fR, fStretch);
	
//	cout << "\nfR: "<<fR<<endl;
	
	double p0 = 0.25*fN*fk*fT/fA;
	for (int i = 0; i < fn; i++) {
		const double r = fR[i];
		double p2 = (r/fL - 0.25 + 0.25/((1-r/fL)*(1-r/fL))); 
		const dSymMatrixT& M = fM[i];
		fStress3D.AddScaled(p0*p2*fR0*fR0/r, M);
	}
	/*repulsion part*/
	Compute_eigs(fEigs, fStretch);
//	cout << "\nfEigs: "<<fEigs<<endl;
	double p1 = sqrt(2.0*fA/fL);
	double p2 = -0.25*fN*fk*fT/(fA*fA*sqrt(2.0*fL/fA)) * (p1 - 0.25 + 0.25/((1-p1)*(1-p1)));
	double p5 = p2*(fl1*fl1/fEigs[0]);
	double p6 = p2*(fl2*fl2/fEigs[1]);
	double p7 = p2*(fl3*fl3/fEigs[2]);
	
	fStress3D.AddScaled(p5, fM[4]);
	fStress3D.AddScaled(p6, fM[5]);
	fStress3D.AddScaled(p7, fM[6]);

	/*bulk part*/
	double I3 = fStretch.Det();
	dSymMatrixT& bulk = fStretch.Inverse();
	bulk *= -pow(I3, -fbeta);
	bulk.PlusIdentity(1.0);
	fStress3D.AddScaled(2.0*fgamma, bulk);
//	cout << "\nfStress3D: "<<fStress3D<<endl;
	
    if (NumSD() == 2)
    {
        fStress[0] = fStress3D[0];
        fStress[1] = fStress3D[1];
        fStress[2] = fStress3D[5];
    }
    else fStress = fStress3D;

	return fStress;
}

/* material description */
const dMatrixT& WLC::c_ijkl(void)
{
    /* deformation gradient */
    const dMatrixT& Fmat = F();
  
    /* transform */
    fModulus.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, C_IJKL()));
    return fModulus;	
}

const dSymMatrixT& WLC::s_ij(void)
{
    /* deformation gradient */
    const dMatrixT& Fmat = F();
  
    /* transform */
    fStress.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, S_IJ()));
    return fStress;
}


void WLC::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FSSolidMatT::DefineParameters(list);

	/* common limit */
	LimitT positive(0.0, LimitT::Lower);

	/* general parameters */
	ParameterT persistence_length(ParameterT::Double, "persistence_length");
	persistence_length.AddLimit(positive);
	list.AddParameter(persistence_length);	
	ParameterT N_links(ParameterT::Double, "N_links");
	ParameterT temperature_K(ParameterT::Double, "temperature_K");
	N_links.AddLimit(positive);
	temperature_K.AddLimit(positive);
	list.AddParameter(N_links);
	list.AddParameter(temperature_K);

	/* unit cell */
	ParameterT l1(ParameterT::Double, "unit_cell_dimension_l1");
	ParameterT a_x(ParameterT::Double,"unit_cell_orientation_a_x");
	ParameterT a_y(ParameterT::Double,"unit_cell_orientation_a_y");
	ParameterT a_z(ParameterT::Double,"unit_cell_orientation_a_z");
	
	ParameterT l2(ParameterT::Double, "unit_cell_dimension_l2");
	ParameterT b_x(ParameterT::Double,"unit_cell_orientation_b_x");
	ParameterT b_y(ParameterT::Double,"unit_cell_orientation_b_y");
	ParameterT b_z(ParameterT::Double,"unit_cell_orientation_b_z");

	ParameterT l3(ParameterT::Double, "unit_cell_dimension_l3");
	ParameterT c_x(ParameterT::Double,"unit_cell_orientation_c_x");
	ParameterT c_y(ParameterT::Double,"unit_cell_orientation_c_y");
	ParameterT c_z(ParameterT::Double,"unit_cell_orientation_c_z");

	l1.AddLimit(positive);
	l2.AddLimit(positive);
	l3.AddLimit(positive);

	list.AddParameter(l1);
	list.AddParameter(a_x);
	list.AddParameter(a_y);
	list.AddParameter(a_z);

	list.AddParameter(l2);
	list.AddParameter(b_x);
	list.AddParameter(b_y);
	list.AddParameter(b_z);

	list.AddParameter(l3);
	list.AddParameter(c_x);
	list.AddParameter(c_y);
	list.AddParameter(c_z);
	
	/*bulk parameters*/
	ParameterT gamma(ParameterT::Double, "bulk_response_gamma");
	ParameterT beta(ParameterT::Double,"bulk_response_beta");
	gamma.AddLimit(positive);
	beta.AddLimit(positive);
	list.AddParameter(gamma);
	list.AddParameter(beta);
}

void WLC::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSSolidMatT::TakeParameterList(list);

	/* dimension work space */
	fC.Dimension(NumSD());
	fStretch.Dimension(3);
	fEigs.Dimension(3);
	fStress.Dimension(NumSD());
	fModulus.Dimension(dSymMatrixT::NumValues(NumSD()));
	fStress3D.Dimension(3);
	fModulus3D.Dimension(6);
	f_e1.Dimension(3);
	f_e2.Dimension(3);
	f_e3.Dimension(3);
	fP.Dimension(3);

	/* general parameters */
	fN = list.GetParameter("N_links");
	fA = list.GetParameter("persistence_length");
	fT = list.GetParameter("temperature_K");

	/* elastic properties */
	fl1 = list.GetParameter("unit_cell_dimension_l1");
	f_e1[0] = list.GetParameter("unit_cell_orientation_a_x");
	f_e1[1] = list.GetParameter("unit_cell_orientation_a_y");
	f_e1[2] = list.GetParameter("unit_cell_orientation_a_z");

	fl2 = list.GetParameter("unit_cell_dimension_l2");
	f_e2[0] = list.GetParameter("unit_cell_orientation_b_x");
	f_e2[1] = list.GetParameter("unit_cell_orientation_b_y");
	f_e2[2] = list.GetParameter("unit_cell_orientation_b_z");

	fl3 = list.GetParameter("unit_cell_dimension_l3");
	f_e3[0] = list.GetParameter("unit_cell_orientation_b_x");
	f_e3[1] = list.GetParameter("unit_cell_orientation_b_y");
	f_e3[2] = list.GetParameter("unit_cell_orientation_b_z");

	fR0 = sqrt(fl1*fl1+fl2*fl2+fl3*fl3)*0.5;
	fL = fR0*fR0/(2*fA);

	/* bulk parameters */
	fgamma = list.GetParameter("bulk_response_gamma");
	fbeta = list.GetParameter("bulk_response_beta");
	
	/*Set Representative Bond Vectors and Structure tensor*/
	fR.Dimension(fn);
	fM.Dimension(fn+3);
	for (int i = 0; i < fn+3; i++) 
		fM[i].Dimension(3);
		
	double iR0 =1.0/fR0;
	
	/*P1, M1*/
	fP.SetToCombination(0.5*fl1*iR0, f_e1, 0.5*fl2*iR0, f_e2, 0.5*fl3*iR0, f_e3);
	fM[0].Outer(fP, 1.0);
	/*P2, M2*/
	fP.SetToCombination(0.5*fl1*iR0, f_e1, 0.5*fl2*iR0, f_e2, -0.5*fl3*iR0, f_e3);
	fM[1].Outer(fP, 1.0);
	/*P3, M3*/
	fP.SetToCombination(0.5*fl1*iR0, f_e1, -0.5*fl2*iR0, f_e2, 0.5*fl3*iR0, f_e3);
	fM[2].Outer(fP, 1.0);
	/*P4, M4*/
	fP.SetToCombination(0.5*fl1*iR0, f_e1, -0.5*fl2*iR0, f_e2, -0.5*fl3*iR0, f_e3);
	fM[3].Outer(fP, 1.0);
		
	/*M5-M8*/
	fM[4].Outer(f_e1, 1.0);
	fM[5].Outer(f_e2, 1.0);
	fM[6].Outer(f_e3, 1.0);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

inline void WLC::Compute_R(dArrayT& r, const dSymMatrixT& C) 
{		
	for (int i = 0; i < fn; i++) 
		r[i] = fR0*sqrt(fM[i].ScalarProduct(C));
}

inline void WLC::Compute_eigs(dArrayT& eigs, const dSymMatrixT& C) 
{
	for (int i = 0; i < 3; i++)
		eigs[i] = fM[i+fn].ScalarProduct(C);
}
