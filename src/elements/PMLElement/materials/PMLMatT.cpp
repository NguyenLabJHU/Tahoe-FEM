/* $Id: PMLMatT.cpp,v 1.9 2004-07-15 08:28:02 paklein Exp $ */
/* created: TDN (5/31/2001) */
#include "PMLMatT.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "ifstreamT.h"
#include "LinearDecreaseT.h"
#include "PMLT.h"

using namespace Tahoe;

/* constructor */
PMLMatT::PMLMatT(ifstreamT& in, const MaterialSupportT& support, const PMLT& element):
	ParameterInterfaceT("PML_material"),
//	SolidMaterialT(support),
//	Material2DT(in),
	IsotropicT(in),
	fPMLElement(element),
	fNumSD(NumSD()),
	fNumDOF(NumDOF()),
	fStress(NumSD()),
	fModulus(dSymMatrixT::NumValues(NumSD())),
	fdt(support.TimeStep())
{
 	int code;
	in >> code;
	switch (code)
	{
		case (klinear):
		{
			double A, L;
			in >> A >> L;
			fFuna = new LinearDecreaseT(A,L);
			in >> A >> L;
			fFunb = new LinearDecreaseT(A,L);
 			break;
		}
	    case (kquadratic):
	    {
			double A, L;
			in >> A >> L;
			fFuna = new LinearDecreaseT(A,L);
			in >> A >> L;
			fFunb = new LinearDecreaseT(A,L);
 			break;
	    }
        default:	
	    {
	      throw ExceptionT::kBadInputValue;
	      break;
	    }
	}
	in >> fRefCoorda;
	in >> fRefCoordb;

	int numstress = fNumSD*(fNumSD+1)/2;
		
	/*allocates storage for state variable array*/
	fnstatev = 0;
	fnstatev += numstress;   /*current sigma_a*/
	fnstatev += numstress;   /*current sigma_b*/
	fnstatev += numstress;   /*current sigma0_a*/
	fnstatev += numstress;   /*current sigma0_b*/
	
	fnstatev += numstress;   /*previous sigma_a*/
	fnstatev += numstress;   /*previous sigma_b*/
	fnstatev += numstress; 	 /*previous sigma0_a*/
	fnstatev += numstress;   /*previous sigma0_b*/
	
	fstatev.Dimension(fnstatev);
	double* pstatev = fstatev.Pointer();
	/* assign pointers to current and preceding blocks of state variable array */
	
	fStressa_n1.Set(fNumSD, pstatev);        
	pstatev += numstress;
	fStressb_n1.Set(fNumSD, pstatev);	
	pstatev += numstress;
	fStress0a_n1.Set(fNumSD, pstatev);
	pstatev += numstress;
	fStress0b_n1.Set(fNumSD, pstatev);
	pstatev += numstress;
	
	fStressa_n.Set(fNumSD, pstatev);        
	pstatev += numstress;
	fStressb_n.Set(fNumSD, pstatev);	
	pstatev += numstress;
	fStress0a_n.Set(fNumSD, pstatev);
	pstatev += numstress;
	fStress0b_n.Set(fNumSD, pstatev);
	pstatev += numstress;
	
}	

void PMLMatT::PointInitialize(void)
{
	/* allocate element storage */
	if (CurrIP() == 0)
	{
		ElementCardT& element = CurrentElement();
		element.Dimension(0, fnstatev*NumIP());
	
	/* initialize internal variables to 0.0*/
		element.DoubleData() = 0.0;
	}

	/* store results as last converged */
	if (CurrIP() == NumIP() - 1) 
		UpdateHistory();
}
 
void PMLMatT::UpdateHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables */
		Load(element, ip);
	
		/* assign "current" to "last" */	
		fStressa_n = fStressa_n1;
		fStressb_n = fStressb_n1;
		fStress0a_n = fStress0a_n1;
		fStress0b_n = fStress0b_n1;
	
		/* write to storage */
		Store(element, ip);
	}
}

void PMLMatT::ResetHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables*/
		Load(element, ip);
	
		/* assign "last" to "current" */
		fStressa_n1 = fStressa_n;
		fStressb_n1 = fStressb_n;
		fStress0a_n1 = fStress0a_n;
		fStress0b_n1 = fStress0b_n;
		
		/* write to storage */
		Store(element, ip);
	}
}


void PMLMatT::Load(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pd = d_array.Pointer(fnstatev*ip);
	double* pdr = fstatev.Pointer();
	for (int i = 0; i < fnstatev; i++)
		*pdr++ = *pd++;
}

void PMLMatT::Store(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pdr = fstatev.Pointer();
	double* pd = d_array.Pointer(fnstatev*ip);
	for (int i = 0; i < fnstatev; i++)
		*pd++ = *pdr++;
}
	
/* spatial description */
const dMatrixT& PMLMatT::c_ijkl(void) 
{ 
ExceptionT::GeneralFail("PMLMatT::c_ijkl", "out of date");
#if 0	
	ComputeModuli2D(fModulus, fConstraintOption);
#endif
	return(fModulus);
}

const dSymMatrixT& PMLMatT::s_ij(void)
{
	dArrayT ip_coords(fNumSD);
	/*calculate normalized time step i.e. dt/tau*/
	const double& dampa = DampFacta(ip_coords);
	const double& dampb = DampFactb(ip_coords);
	
	double cscale_a = 1.0/(1.0+dampa*fdt*0.5);
	double cscale_b = 1.0/(1.0+dampb*fdt*0.5);
	double sscale_a = (1.0-dampa*fdt*0.5)*cscale_a;
	double sscale_b = (1.0-dampb*fdt*0.5)*cscale_b;
	
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	
	/* obtain displacement gradients*/
	const dMatrixT& gradU = fPMLElement.GradU();
		
 //TEMP - revised Update/ComputeOutput sequence may make this
 //       check unnecessary	
	if (1 /* || fRunState == GlobalT::kFormRHS */)
	{
		const dMatrixT& modulus = c_ijkl(); 
//Q:What is the ordering of GradU?  GradU(i,j) = Ui,j?
		fStress0a_n1(0,0) = cscale_a*modulus(0,0)*gradU(0,0);
		fStress0a_n1(1,1) = cscale_a*modulus(0,1)*gradU(0,0);
		fStress0a_n1(0,1) = cscale_a*modulus(2,2)*gradU(1,0);
		
		fStress0b_n1(0,0) = cscale_b*modulus(1,0)*gradU(1,1);
		fStress0b_n1(1,1) = cscale_b*modulus(1,1)*gradU(1,1);
		fStress0b_n1(0,1) = cscale_b*modulus(2,2)*gradU(0,1);
		
		fStressa_n1.SetToCombination(sscale_a,fStressa_n,1.0,fStress0a_n1,1.0,fStress0a_n);
		fStressb_n1.SetToCombination(sscale_b,fStressb_n,1.0,fStress0b_n1,1.0,fStress0b_n);
		
		
		/*Store State Variables*/
		Store(element, CurrIP());
	}
	else 
	{
		/*calculate overstress*/	
		fStressa_n1 = fStressa_n;
		fStressb_n1 = fStressb_n;
	}

	fStress+= fStressa_n1;
	fStress+= fStressb_n1;
	return (fStress);
}

const dMatrixT& PMLMatT::C_IJKL(void) 
{ 
	return(c_ijkl());
}

const dSymMatrixT& PMLMatT::S_IJ(void)
{
	return(s_ij());
}

const double& PMLMatT::DampFacta(dArrayT& ip_coords)
{
	double del=fabs(ip_coords[0]-fRefCoorda);

	fDampa = fFuna->Function(del);
	if (fDampa < 0)
	{
	  	cout<<"Dimension of viscous layer is too small";
	    throw ExceptionT::kBadInputValue;
	}
     return(fDampa);
}

const double& PMLMatT::DampFactb(dArrayT& ip_coords)
{
	double del=fabs(ip_coords[1]-fRefCoordb);

	fDampb = fFunb->Function(del);
	if (fDampb < 0)
	{
	  	cout<<"Dimension of viscous layer is too small";
	    throw ExceptionT::kBadInputValue;
	}
     return(fDampb);
}

	
double PMLMatT::StrainEnergyDensity(void)
{
	/*obtain displacement gradients*/
	const dMatrixT& gradU = fPMLElement.GradU();
	dSymMatrixT strain;
	strain.Symmetrize(gradU);
	const dSymMatrixT& stress = s_ij();
	
	double energy = 0;
	
	energy = 0.5* (strain(0,0)*stress(0,0)+strain(1,1)*stress(1,1));
	energy += strain(0,1)*stress(0,1);
	
	return(0.5*energy);
}

