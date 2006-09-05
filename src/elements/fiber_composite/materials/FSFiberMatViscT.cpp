/* $Id: FSFiberMatViscT.cpp,v 1.2 2006-09-05 23:10:23 thao Exp $ */
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
		
		AssembleFiberModuli(fCalg, fModulus);
		j++;
	}
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
	
	fFiberStretch_n.Dimension(fNumSD);
	fFiberStretch_v.Dimension(fNumSD);
	fFiberStretch_vn.Dimension(fNumSD);

	/*Dimension work spaces*/
	fCalg.Dimension(fNumFibStress);

	int numprocess = fNumFibProcess+fNumMatProcess;
	if (numprocess > 0) SetStateVariables(numprocess);
}

