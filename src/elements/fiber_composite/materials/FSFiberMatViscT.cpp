/* $Id: FSFiberMatViscT.cpp,v 1.1 2006-08-10 01:35:44 thao Exp $ */
/* created: paklein (06/09/1997) */
#include "FSFiberMatViscT.h"
#include "FSFiberMatSupportT.h"
#include "iArray2DT.h"

using namespace Tahoe;

/* constructor */
FSFiberMatViscT::FSFiberMatViscT(void):
	ParameterInterfaceT("vicoelastic_fiber_composite_material")
{
	/*Reset default*/
	fNumFibProcess = 1;
	fNumMatProcess = 1;
}


/* modulus */
const dMatrixT& FSFiberMatViscT::C_IJKL(void)
{
	/* stretch */
	Compute_C(fC);
	
	/*equilibrium contribution*/
	/*calculate eq. matrix contribution*/
	ComputeMatrixMod(fC, fStress, fModulus);
	
	/* eq. fiber contribution*/
	ComputeFiberStretch(fC, fFiberStretch);
	ComputeFiberMod(fFiberStretch, fFiberStress, fFiberMod);
				
	/* rotate and assemble eq. stress to lab coordinates */
	AssembleFiberStress(fFiberStress, fStress);
	
	/* rotate and assemble eq. modulus to lab coordinates */
	AssembleFiberModuli(fFiberMod, fModulus);

	/*calculate nonequilibrium contribution*/
	/*Load state variables (Cv and Cvn)*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());

	for (int i = 0; i < fNumMatProcess; i++)
	{
		/*calculate neq. matrix contribution*/
		ComputeMatrixMod(fC, fC_v[i], fStress, fModulus, i, dSymMatrixT::kAccumulate);
	}

	int j = fNumMatProcess;
	for (int i = 0; i < fNumFibProcess && fNumFibProcess > 0; i++)
	{
		/* neq. fiber contribution*/
		ComputeFiberStretch(fC_v[j], fFiberStretch_v);

		/*calculate SNEQ and dSNEQ/dC*/
		ComputeFiberMod(fFiberStretch, fFiberStretch_v, fFiberStress, fFiberMod, i);
				
		/* rotate neq. stress to lab coordinates and assemble in fStress */
		AssembleFiberStress(fFiberStress, fStress);

		/* rotate and assemble eq. modulus to lab coordinates */
		AssembleFiberModuli(fFiberMod, fModulus);
		
		/*compute dSNEQ/dCv delta Cv/delta C*/
		ComputeCalg(fFiberStretch, fFiberStretch_v, fCalg, i);
		
/*		dMatrixT test(6);
		AssembleFiberNonSymModuli(fCalg, test,dSymMatrixT::kOverwrite);
		cout <<"\nCalg: "<<test;
*/
		AssembleFiberNonSymModuli(fCalg, fModulus);
		j++;
	}
/*	for (int i = 0; i<6;i++)
		if( fModulus(i,i)<.1) cout << "\nMod: "<<fModulus;
*/
	return fModulus;
}
	
/* stress */
const dSymMatrixT& FSFiberMatViscT::S_IJ(void)
{
	
	/* stretch */
	Compute_C(fC);
	
	/*matrix contribution*/
	/*calculate matrix contribution*/
	ComputeMatrixStress(fC, fStress);

	/*fiber contribution*/
	ComputeFiberStretch(fC, fFiberStretch);
	ComputeFiberStress(fFiberStretch, fFiberStress);
	
	/* rotate and assemble stress to lab coordinates */
	AssembleFiberStress(fFiberStress, fStress);

	/*calculate nonequilibrium contribution*/
	/*Load state variables (Cv and Cvn)*/
    ElementCardT& element = CurrentElement();
    Load(element, CurrIP());

    if (fFSMatSupport->RunState() == GlobalT::kFormRHS)
    {	
		for (int i = 0; i < fNumMatProcess && fNumMatProcess > 0; i++)
		{	
			/*calculate neq. matrix contribution*/
			ComputeMatrixMod(fC, fC_v[i], fStress, fModulus, i, dSymMatrixT::kAccumulate);
		}
		int j = fNumMatProcess;
		
		for (int i = 0; i < fNumFibProcess && fNumFibProcess > 0; i++)
		{
			/* neq. fiber contribution*/
			if(Need_F_last())
			{
				const dMatrixT& F_last = F_mechanical_last();
				fC_n.MultATA(F_last);
			}
			Compute_Cv(fC_n, fC, fC_vn[j], fC_v[j], i);		
			AssembleFiberStress(fFiberStretch_v, fC_v[j], dSymMatrixT::kOverwrite);

			ComputeFiberStress(fFiberStretch, fFiberStretch_v, fFiberStress, i);
				
			/* rotate neq. stress to lab coordinates and assemble in fStress */
			AssembleFiberStress(fFiberStress, fStress);
			
			j++;
		}
		Store(element, CurrIP());
	}
	else 
	{
		for (int i = 0; i < fNumMatProcess && fNumMatProcess > 0; i++)
		{
			/*calculate neq. matrix contribution*/
			ComputeMatrixStress(fC, fC_v[i], fStress,  i, dSymMatrixT::kAccumulate);
		}
		int j = fNumMatProcess;
		for (int i = 0; i < fNumFibProcess && fNumFibProcess > 0; i++)
		{			
			/* neq. fiber contribution*/
			ComputeFiberStretch(fC_v[i], fFiberStretch_v);
			ComputeFiberStress(fFiberStretch, fFiberStretch_v, fFiberStress, i);
				
			/* rotate neq. stress to lab coordinates and assemble in fStress */
			AssembleFiberStress(fFiberStress, fStress);

			j++;
		}
	}

	return(fStress);
}

/* material description */
const dMatrixT& FSFiberMatViscT::c_ijkl(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fModulus.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, C_IJKL()));
	return fModulus;
}
/**< \todo construct directly in material description */

const dSymMatrixT& FSFiberMatViscT::s_ij(void)
{
	/* deformation gradient */
	const dMatrixT& Fmat = F();
	
	/* transform */
	fStress.SetToScaled(1.0/Fmat.Det(), PushForward(Fmat, S_IJ()));
	return fStress;
}

/*initializes history variable */
void  FSFiberMatViscT::PointInitialize(void)
{
	/* allocate element storage */
	ElementCardT& element = CurrentElement();	
	if (CurrIP() == 0)
	{
		ElementCardT& element = CurrentElement();
		element.Dimension(0, fnstatev*NumIP());
	
		/* initialize internal variables to identity*/
		for (int ip = 0; ip < NumIP(); ip++)
		{
		      /* load state variables */
		      Load(element, ip);
		      
			  int numprocess = fNumFibProcess+fNumMatProcess;
			  for (int i = 0; i < numprocess; i++)
			  {
				fC_vn[i].Identity();
				fC_v[i].Identity();
			  }

		      /* write to storage */
		      Store(element, ip);
		}
	}
}
 
void FSFiberMatViscT::UpdateHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables */
		Load(element, ip);
	
		int numprocess = fNumFibProcess+fNumMatProcess;
		/* assign "current" to "last" */	
		for (int i = 0; i < numprocess; i++)
			fC_vn[i] = fC_v[i];

		/* write to storage */
		Store(element, ip);
	}
}

void FSFiberMatViscT::ResetHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables*/
		Load(element, ip);
	
		/* assign "last" to "current" */
		
		/* write to storage */
		Store(element, ip);
	}
}

void FSFiberMatViscT::Load(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pd = d_array.Pointer(fnstatev*ip);
	double* pdr = fstatev.Pointer();
	for (int i = 0; i < fnstatev; i++)
		*pdr++ = *pd++;
}
void FSFiberMatViscT::Store(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pdr = fstatev.Pointer();
	double* pd = d_array.Pointer(fnstatev*ip);
	for (int i = 0; i < fnstatev; i++)
		*pd++ = *pdr++;

}

/* accept parameter list */
void FSFiberMatViscT::SetStateVariables(const int numprocess)
{

	fC_v.Dimension(numprocess);
	fC_vn.Dimension(numprocess);

	int ndof = 3;
	int numstress = dSymMatrixT::NumValues(ndof);

	fnstatev = 0;
	fnstatev += numstress;   /*current C_v*/
	fnstatev += numstress;   /*last C_vn*/

	fnstatev *= numprocess;
	
	fstatev.Dimension(fnstatev);
	double* pstatev = fstatev.Pointer();
		
	/* assign pointers to current and last blocks of state variable array */
	for (int i = 0; i < numprocess; i++)
	{
		fC_v[i].Set(ndof, pstatev);
		pstatev += numstress;
		fC_vn[i].Set(ndof, pstatev);
		pstatev += numstress;
	}
}


/* accept parameter list */
void FSFiberMatViscT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSFiberMatT::TakeParameterList(list);
	/* allocate memory */
	/*2D fiber stress and modulus*/
	fC_n.Dimension(fNumSD);
	
	fFiberStretch_n.Dimension(fNumSD-1);
	fFiberStretch_v.Dimension(fNumSD-1);
	fFiberStretch_vn.Dimension(fNumSD-1);

	/*Dimension work spaces*/
	fCalg.Dimension(fNumFibStress);

	int numprocess = fNumFibProcess+fNumMatProcess;
	if (numprocess > 0) SetStateVariables(numprocess);
}


/*********************************************************************************************
 *protected                                                                                  *
 *********************************************************************************************/
void FSFiberMatViscT::AssembleFiberNonSymModuli(const dMatrixT& cf, dMatrixT& mod, const int fillmode)
{
	const double& c11 = cf(0,0);
	const double& c22 = cf(1,1);
	const double& c33 = cf(2,2);
	
	const double& c23 = cf(1,2);
	const double& c32 = cf(2,1);
	const double& c13 = cf(0,2);
	const double& c31 = cf(2,0);
	const double& c12 = cf(0,1);
	const double& c21 = cf(1,0);
	
	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();

	const double x1 = Fibers(0,0);
	const double x2 = Fibers(0,1);
	const double x3 = Fibers(0,2);
	
	const double y1 = Fibers(1,0);
	const double y2 = Fibers(1,1);
	const double y3 = Fibers(1,2);
	
	/*Rotate moduli from local frame (defined by fibrils) to global cartesian frame*/
	/*C_IJKL = QIa QJb QKc QLd Cf_abcd*/

	if (fillmode == dSymMatrixT::kOverwrite)
	mod = 0.0;
	
	mod(0,0) += c11*x1*x1*x1*x1 + 2*c13*x1*x1*x1*y1 + 2*c31*x1*x1*x1*y1 + c12*x1*x1*y1*y1 + c21*x1*x1*y1*y1 
		+ 4*c33*x1*x1*y1*y1 + 2*c23*x1*y1*y1*y1 + 2*c32*x1*y1*y1*y1 + c22*y1*y1*y1*y1;

	mod(1,1) += c11*x2*x2*x2*x2 + 2*c13*x2*x2*x2*y2 + 2*c31*x2*x2*x2*y2 + c12*x2*x2*y2*y2 + c21*x2*x2*y2*y2 
		+ 4*c33*x2*x2*y2*y2 + 2*c23*x2*y2*y2*y2 + 2*c32*x2*y2*y2*y2 + c22*y2*y2*y2*y2;

	mod(2,2) += c11*x3*x3*x3*x3 + 2*c13*x3*x3*x3*y3 + 2*c31*x3*x3*x3*y3 + c12*x3*x3*y3*y3 + c21*x3*x3*y3*y3 
		+ 4*c33*x3*x3*y3*y3 + 2*c23*x3*y3*y3*y3 + 2*c32*x3*y3*y3*y3 + c22*y3*y3*y3*y3;

	mod(3,3) += c11*x2*x2*x3*x3 + c13*x2*x3*x3*y2 + c31*x2*x3*x3*y2 + c33*x3*x3*y2*y2 + c13*x2*x2*x3*y3 
		+ c31*x2*x2*x3*y3 + c12*x2*x3*y2*y3 + c21*x2*x3*y2*y3 + 2*c33*x2*x3*y2*y3 + c23*x3*y2*y2*y3 + c32*x3*y2*y2*y3 
		+ c33*x2*x2*y3*y3 + c23*x2*y2*y3*y3 + c32*x2*y2*y3*y3 + c22*y2*y2*y3*y3;
		
	mod(4,4) += c11*x1*x1*x3*x3 + c13*x1*x3*x3*y1 + c31*x1*x3*x3*y1 + c33*x3*x3*y1*y1 + c13*x1*x1*x3*y3 
		+ c31*x1*x1*x3*y3 + c12*x1*x3*y1*y3 + c21*x1*x3*y1*y3 + 2*c33*x1*x3*y1*y3 + c23*x3*y1*y1*y3 + c32*x3*y1*y1*y3 
		+ c33*x1*x1*y3*y3 + c23*x1*y1*y3*y3 + c32*x1*y1*y3*y3 + c22*y1*y1*y3*y3;

	mod(5,5) += c11*x1*x1*x2*x2 + c13*x1*x2*x2*y1 + c31*x1*x2*x2*y1 + c33*x2*x2*y1*y1 + c13*x1*x1*x2*y2 
		+ c31*x1*x1*x2*y2 + c12*x1*x2*y1*y2 + c21*x1*x2*y1*y2 + 2*c33*x1*x2*y1*y2 + c23*x2*y1*y1*y2 + c32*x2*y1*y1*y2 
		+ c33*x1*x1*y2*y2 + c23*x1*y1*y2*y2 + c32*x1*y1*y2*y2 + c22*y1*y1*y2*y2;
		
		
	mod(0,1) += c11*x1*x1*x2*x2 + 2*c31*x1*x2*x2*y1 + c21*x2*x2*y1*y1 + 2*c13*x1*x1*x2*y2 + 4*c33*x1*x2*y1*y2 
		+ 2*c23*x2*y1*y1*y2 + c12*x1*x1*y2*y2 + 2*c32*x1*y1*y2*y2 + c22*y1*y1*y2*y2;

	mod(0,2) += c11*x1*x1*x3*x3 + 2*c31*x1*x3*x3*y1 + c21*x3*x3*y1*y1 + 2*c13*x1*x1*x3*y3 + 4*c33*x1*x3*y1*y3 
		+ 2*c23*x3*y1*y1*y3 + c12*x1*x1*y3*y3 + 2*c32*x1*y1*y3*y3 + c22*y1*y1*y3*y3;

	mod(0,3) += c11*x1*x1*x2*x3 + 2*c31*x1*x2*x3*y1 + c21*x2*x3*y1*y1 + c13*x1*x1*x3*y2 + 2*c33*x1*x3*y1*y2 + c23*x3*y1*y1*y2 
		+ c13*x1*x1*x2*y3 + 2*c33*x1*x2*y1*y3 + c23*x2*y1*y1*y3 + c12*x1*x1*y2*y3 + 2*c32*x1*y1*y2*y3 + c22*y1*y1*y2*y3;

	mod(0,4) += c11*x1*x1*x1*x3 + c13*x1*x1*x3*y1 + 2*c31*x1*x1*x3*y1 + c21*x1*x3*y1*y1 + 2*c33*x1*x3*y1*y1 + c23*x3*y1*y1*y1 
		+ c13*x1*x1*x1*y3 + c12*x1*x1*y1*y3 + 2*c33*x1*x1*y1*y3 + c23*x1*y1*y1*y3 + 2*c32*x1*y1*y1*y3 + c22*y1*y1*y1*y3;

	mod(0,5) += c11*x1*x1*x1*x2 + c13*x1*x1*x2*y1 + 2*c31*x1*x1*x2*y1 + c21*x1*x2*y1*y1 + 2*c33*x1*x2*y1*y1 + c23*x2*y1*y1*y1 
		+ c13*x1*x1*x1*y2 + c12*x1*x1*y1*y2 + 2*c33*x1*x1*y1*y2 + c23*x1*y1*y1*y2 + 2*c32*x1*y1*y1*y2 + c22*y1*y1*y1*y2;


	mod(1,0) += c11*x1*x1*x2*x2 + 2*c13*x1*x2*x2*y1 + c12*x2*x2*y1*y1 + 2*c31*x1*x1*x2*y2 + 4*c33*x1*x2*y1*y2 
		+ 2*c32*x2*y1*y1*y2 + c21*x1*x1*y2*y2 + 2*c23*x1*y1*y2*y2 + c22*y1*y1*y2*y2;

	mod(1,2) += c11*x2*x2*x3*x3 + 2*c31*x2*x3*x3*y2 + c21*x3*x3*y2*y2 + 2*c13*x2*x2*x3*y3 + 4*c33*x2*x3*y2*y3 
		+ 2*c23*x3*y2*y2*y3 + c12*x2*x2*y3*y3 + 2*c32*x2*y2*y3*y3 + c22*y2*y2*y3*y3;

	mod(1,3) += c11*x2*x2*x2*x3 + c13*x2*x2*x3*y2 + 2*c31*x2*x2*x3*y2 + c21*x2*x3*y2*y2 + 2*c33*x2*x3*y2*y2 + c23*x3*y2*y2*y2 
		+ c13*x2*x2*x2*y3 + c12*x2*x2*y2*y3 + 2*c33*x2*x2*y2*y3 + c23*x2*y2*y2*y3 + 2*c32*x2*y2*y2*y3 + c22*y2*y2*y2*y3;

	mod(1,4) += c11*x1*x2*x2*x3 + c13*x2*x2*x3*y1 + 2*c31*x1*x2*x3*y2 + 2*c33*x2*x3*y1*y2 + c21*x1*x3*y2*y2 + c23*x3*y1*y2*y2 
		+ c13*x1*x2*x2*y3 + c12*x2*x2*y1*y3 + 2*c33*x1*x2*y2*y3 + 2*c32*x2*y1*y2*y3 + c23*x1*y2*y2*y3 + c22*y1*y2*y2*y3;
	
	mod(1,5) += c11*x1*x2*x2*x2 + c13*x2*x2*x2*y1 + c13*x1*x2*x2*y2 + 2*c31*x1*x2*x2*y2 + c12*x2*x2*y1*y2 + 2*c33*x2*x2*y1*y2 
		+ c21*x1*x2*y2*y2 + 2*c33*x1*x2*y2*y2 + c23*x2*y1*y2*y2 + 2*c32*x2*y1*y2*y2 + c23*x1*y2*y2*y2 + c22*y1*y2*y2*y2;


	mod(2,0) += c11*x1*x1*x3*x3 + 2*c13*x1*x3*x3*y1 + c12*x3*x3*y1*y1 + 2*c31*x1*x1*x3*y3 + 4*c33*x1*x3*y1*y3 
		+ 2*c32*x3*y1*y1*y3 + c21*x1*x1*y3*y3 + 2*c23*x1*y1*y3*y3 + c22*y1*y1*y3*y3;

	mod(2,1) += c11*x2*x2*x3*x3 + 2*c13*x2*x3*x3*y2 + c12*x3*x3*y2*y2 + 2*c31*x2*x2*x3*y3 + 4*c33*x2*x3*y2*y3 
		+ 2*c32*x3*y2*y2*y3 + c21*x2*x2*y3*y3 + 2*c23*x2*y2*y3*y3 + c22*y2*y2*y3*y3;

	mod(2,3) += c11*x2*x3*x3*x3 + c13*x3*x3*x3*y2 + c13*x2*x3*x3*y3 + 2*c31*x2*x3*x3*y3 + c12*x3*x3*y2*y3 + 2*c33*x3*x3*y2*y3 
		+ c21*x2*x3*y3*y3 + 2*c33*x2*x3*y3*y3 + c23*x3*y2*y3*y3 + 2*c32*x3*y2*y3*y3 + c23*x2*y3*y3*y3 + c22*y2*y3*y3*y3;

	mod(2,4) += c11*x1*x3*x3*x3 + c13*x3*x3*x3*y1 + c13*x1*x3*x3*y3 + 2*c31*x1*x3*x3*y3 + c12*x3*x3*y1*y3 + 2*c33*x3*x3*y1*y3 
		+ c21*x1*x3*y3*y3 + 2*c33*x1*x3*y3*y3 + c23*x3*y1*y3*y3 + 2*c32*x3*y1*y3*y3 + c23*x1*y3*y3*y3 + c22*y1*y3*y3*y3;

	mod(2,5) += c11*x1*x2*x3*x3 + c13*x2*x3*x3*y1 + c13*x1*x3*x3*y2 + c12*x3*x3*y1*y2 + 2*c31*x1*x2*x3*y3 + 2*c33*x2*x3*y1*y3 
		+ 2*c33*x1*x3*y2*y3 + 2*c32*x3*y1*y2*y3 + c21*x1*x2*y3*y3 + c23*x2*y1*y3*y3 + c23*x1*y2*y3*y3 + c22*y1*y2*y3*y3;


	mod(3,0) += c11*x1*x1*x2*x3 + 2*c13*x1*x2*x3*y1 + c12*x2*x3*y1*y1 + c31*x1*x1*x3*y2 + 2*c33*x1*x3*y1*y2 + c32*x3*y1*y1*y2 
		+ c31*x1*x1*x2*y3 + 2*c33*x1*x2*y1*y3 + c32*x2*y1*y1*y3 + c21*x1*x1*y2*y3 + 2*c23*x1*y1*y2*y3 + c22*y1*y1*y2*y3;

	mod(3,1) += c11*x2*x2*x2*x3 + 2*c13*x2*x2*x3*y2 + c31*x2*x2*x3*y2 + c12*x2*x3*y2*y2 + 2*c33*x2*x3*y2*y2 + c32*x3*y2*y2*y2 
		+ c31*x2*x2*x2*y3 + c21*x2*x2*y2*y3 + 2*c33*x2*x2*y2*y3 + 2*c23*x2*y2*y2*y3 + c32*x2*y2*y2*y3 + c22*y2*y2*y2*y3;

	mod(3,2) += c11*x2*x3*x3*x3 + c31*x3*x3*x3*y2 + 2*c13*x2*x3*x3*y3 + c31*x2*x3*x3*y3 + c21*x3*x3*y2*y3 + 2*c33*x3*x3*y2*y3 
		+ c12*x2*x3*y3*y3 + 2*c33*x2*x3*y3*y3 + 2*c23*x3*y2*y3*y3 + c32*x3*y2*y3*y3 + c32*x2*y3*y3*y3 + c22*y2*y3*y3*y3;

	mod(3,4) += c11*x1*x2*x3*x3 + c13*x2*x3*x3*y1 + c31*x1*x3*x3*y2 + c33*x3*x3*y1*y2 + c13*x1*x2*x3*y3 + c31*x1*x2*x3*y3 
		+ c12*x2*x3*y1*y3 + c33*x2*x3*y1*y3 + c21*x1*x3*y2*y3 + c33*x1*x3*y2*y3 + c23*x3*y1*y2*y3 + c32*x3*y1*y2*y3 + c33*x1*x2*y3*y3 
		+ c32*x2*y1*y3*y3 + c23*x1*y2*y3*y3 + c22*y1*y2*y3*y3;

	mod(3,5) += c11*x1*x2*x2*x3 + c13*x2*x2*x3*y1 + c13*x1*x2*x3*y2 + c31*x1*x2*x3*y2 + c12*x2*x3*y1*y2 + c33*x2*x3*y1*y2 
		+ c33*x1*x3*y2*y2 + c32*x3*y1*y2*y2 + c31*x1*x2*x2*y3 + c33*x2*x2*y1*y3 + c21*x1*x2*y2*y3 + c33*x1*x2*y2*y3 
		+ c23*x2*y1*y2*y3 + c32*x2*y1*y2*y3 + c23*x1*y2*y2*y3 + c22*y1*y2*y2*y3;


	mod(4,0) += c11*x1*x1*x1*x3 + 2*c13*x1*x1*x3*y1 + c31*x1*x1*x3*y1 + c12*x1*x3*y1*y1 + 2*c33*x1*x3*y1*y1 
		+ c32*x3*y1*y1*y1 + c31*x1*x1*x1*y3 + c21*x1*x1*y1*y3 + 2*c33*x1*x1*y1*y3 + 2*c23*x1*y1*y1*y3 + c32*x1*y1*y1*y3 
		+ c22*y1*y1*y1*y3;

	mod(4,1) += c11*x1*x2*x2*x3 + c31*x2*x2*x3*y1 + 2*c13*x1*x2*x3*y2 + 2*c33*x2*x3*y1*y2 + c12*x1*x3*y2*y2 + c32*x3*y1*y2*y2 
		+ c31*x1*x2*x2*y3 + c21*x2*x2*y1*y3 + 2*c33*x1*x2*y2*y3 + 2*c23*x2*y1*y2*y3 + c32*x1*y2*y2*y3 + c22*y1*y2*y2*y3;

	mod(4,2) += c11*x1*x3*x3*x3 + c31*x3*x3*x3*y1 + 2*c13*x1*x3*x3*y3 + c31*x1*x3*x3*y3 + c21*x3*x3*y1*y3 + 2*c33*x3*x3*y1*y3 
		+ c12*x1*x3*y3*y3 + 2*c33*x1*x3*y3*y3 + 2*c23*x3*y1*y3*y3 + c32*x3*y1*y3*y3 + c32*x1*y3*y3*y3 + c22*y1*y3*y3*y3;

	mod(4,3) += c11*x1*x2*x3*x3 + c31*x2*x3*x3*y1 + c13*x1*x3*x3*y2 + c33*x3*x3*y1*y2 + c13*x1*x2*x3*y3 + c31*x1*x2*x3*y3 
		+ c21*x2*x3*y1*y3 + c33*x2*x3*y1*y3 + c12*x1*x3*y2*y3 + c33*x1*x3*y2*y3 + c23*x3*y1*y2*y3 + c32*x3*y1*y2*y3 + c33*x1*x2*y3*y3 
		+ c23*x2*y1*y3*y3 + c32*x1*y2*y3*y3 + c22*y1*y2*y3*y3;

	mod(4,5) += c11*x1*x1*x2*x3 + c13*x1*x2*x3*y1 + c31*x1*x2*x3*y1 + c33*x2*x3*y1*y1 + c13*x1*x1*x3*y2 + c12*x1*x3*y1*y2 
		+ c33*x1*x3*y1*y2 + c32*x3*y1*y1*y2 + c31*x1*x1*x2*y3 + c21*x1*x2*y1*y3 + c33*x1*x2*y1*y3 + c23*x2*y1*y1*y3 
		+ c33*x1*x1*y2*y3 + c23*x1*y1*y2*y3 + c32*x1*y1*y2*y3 + c22*y1*y1*y2*y3;
	
	mod(5,0) += c11*x1*x1*x1*x2 + 2*c13*x1*x1*x2*y1 + c31*x1*x1*x2*y1 + c12*x1*x2*y1*y1 + 2*c33*x1*x2*y1*y1 + c32*x2*y1*y1*y1 
		+ c31*x1*x1*x1*y2 + c21*x1*x1*y1*y2 + 2*c33*x1*x1*y1*y2 + 2*c23*x1*y1*y1*y2 + c32*x1*y1*y1*y2 + c22*y1*y1*y1*y2;

	mod(5,1) += c11*x1*x2*x2*x2 + c31*x2*x2*x2*y1 + 2*c13*x1*x2*x2*y2 + c31*x1*x2*x2*y2 + c21*x2*x2*y1*y2 + 2*c33*x2*x2*y1*y2 
		+ c12*x1*x2*y2*y2 + 2*c33*x1*x2*y2*y2 + 2*c23*x2*y1*y2*y2 + c32*x2*y1*y2*y2 + c32*x1*y2*y2*y2 + c22*y1*y2*y2*y2;

	mod(5,2) += c11*x1*x2*x3*x3 + c31*x2*x3*x3*y1 + c31*x1*x3*x3*y2 + c21*x3*x3*y1*y2 + 2*c13*x1*x2*x3*y3 + 2*c33*x2*x3*y1*y3 
		+ 2*c33*x1*x3*y2*y3 + 2*c23*x3*y1*y2*y3 + c12*x1*x2*y3*y3 + c32*x2*y1*y3*y3 + c32*x1*y2*y3*y3 + c22*y1*y2*y3*y3;

	mod(5,3) += c11*x1*x2*x2*x3 + c31*x2*x2*x3*y1 + c13*x1*x2*x3*y2 + c31*x1*x2*x3*y2 + c21*x2*x3*y1*y2 + c33*x2*x3*y1*y2 + c33*x1*x3*y2*y2 
		+ c23*x3*y1*y2*y2 + c13*x1*x2*x2*y3 + c33*x2*x2*y1*y3 + c12*x1*x2*y2*y3 + c33*x1*x2*y2*y3 + c23*x2*y1*y2*y3 + c32*x2*y1*y2*y3 
		+ c32*x1*y2*y2*y3 + c22*y1*y2*y2*y3;

	mod(5,4) += c11*x1*x1*x2*x3 + c13*x1*x2*x3*y1 + c31*x1*x2*x3*y1 + c33*x2*x3*y1*y1 + c31*x1*x1*x3*y2 + c21*x1*x3*y1*y2 + c33*x1*x3*y1*y2 
		+ c23*x3*y1*y1*y2 + c13*x1*x1*x2*y3 + c12*x1*x2*y1*y3 + c33*x1*x2*y1*y3 + c32*x2*y1*y1*y3 + c33*x1*x1*y2*y3 
		+ c23*x1*y1*y2*y3 + c32*x1*y1*y2*y3 + c22*y1*y1*y2*y3;
}
