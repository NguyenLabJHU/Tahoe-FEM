/* $Id: FSFiberMatViscT.cpp,v 1.3 2006-10-20 20:02:38 thao Exp $ */
/* created: paklein (06/09/1997) */
#include "FSFiberMatViscT.h"
#include "FSFiberMatSupportT.h"
#include "ParameterContainerT.h"
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

		/* rotate neq. stress to lab coordinates and assemble in fStress */
		AssembleFiberStress(fFiberStress, fStress);

		/*calculate SNEQ and dSNEQ/dC*/
		ComputeFiberMod(fFiberStretch, fFiberStretch_v, fFiberStress, fFiberMod, i);				

		/* rotate and assemble neq. modulus to lab coordinates */
		AssembleFiberModuli(fFiberMod, fModulus);

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
			ComputeMatrixCv(fC, fC_vn[i], fC_v[i], i);
			ComputeMatrixStress(fC, fC_v[i], fStress, i, dSymMatrixT::kAccumulate);
		}
		int j = fNumMatProcess;
		
		for (int i = 0; i < fNumFibProcess && fNumFibProcess > 0; i++)
		{
			/* neq. fiber contribution*/
/*			if(Need_F_last())
			{
				const dMatrixT& F_last = F_mechanical_last();
				fC_n.MultATA(F_last);
			}
			Compute_Cv(fC_n, fC, fC_vn[j], fC_v[j], i);		
			*/
			/*rotate stretches in fiber plane*/
			ComputeFiberStretch(fC_vn[j], fFiberStretch_vn);
			ComputeFiberStretch(fC_v[j], fFiberStretch_v);

			Compute_Cv(fFiberStretch, fFiberStretch_vn, fFiberStretch_v, i);		
			/*Rotate fFiberStretch_v from local fiber coords to lab coordinates to get C_v*/		
			AssembleFiberStress(fFiberStretch_v, fC_v[j], dSymMatrixT::kOverwrite);

			/*compute fiber stress*/
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
			ComputeFiberStretch(fC_v[j], fFiberStretch_v);
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
void FSFiberMatViscT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSFiberMatT::TakeParameterList(list);
	/*Dimension work spaces*/
	/*C at time step n
	fC_n.Dimension(fNumSD); 
	fFiberStretch_n.Dimension(fNumSD); */

//	fQ.Dimension(fNumSD);
	/*2D fiber stress and modulus*/
//	fFiberStretch.Dimension(fNumSD-1);
//	fFiberStress.Dimension(fNumSD-1);
//	fFiberMod.Dimension(fNumFibStress);

	/*viscous fiber stretch at time step n, viscous stretch at time step n and vn*/
	fFiberStretch_v.Dimension(fNumSD);
	fFiberStretch_vn.Dimension(fNumSD);

	/*Dimension work spaces*/
	fCalg.Dimension(fNumFibStress);	
}

/* accept parameter list */
void FSFiberMatViscT::SetStateVariables(const int numprocess)
{
	/*dimension state variables*/
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

/************************************ private **********************************************/
/* construct symmetric rank-4 mixed-direction tensor (6.1.44) */
void FSFiberMatViscT::MixedRank4_2D(const dArrayT& a, const dArrayT& b, 
	dMatrixT& rank4_ab) const
{
#if __option(extended_errorcheck)
	if (a.Length() != 2 ||
	    b.Length() != 2 ||
	    rank4_ab.Rows() != 3 ||
	    rank4_ab.Cols() != 3) throw ExceptionT::kSizeMismatch;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11;

//	z1 = A(1.);
//	z2 = A(2.);
//	z3 = B(1.);
//	z4 = B(2.);

	z1 = a[0];
	z2 = a[1];
	z3 = b[0];
	z4 = b[1];

	z5 = z1*z1;
	z6 = z2*z2;
	z7 = z3*z3;
	z8 = 2.*z1*z2*z3*z4;
	z9 = z4*z4;
	z3 = 2.*z3*z4;
	z4 = z3*z5;
	z3 = z3*z6;
	z10 = 2.*z1*z2*z7;
	z11 = 2.*z5*z7;
	z7 = z6*z7;
	z1 = 2.*z1*z2*z9;
	z2 = z5*z9;
	z5 = 2.*z6*z9;
	z4 = z10 + z4;
	z1 = z1 + z3;
	z2 = z2 + z7 + z8;
	z3 = 0.5*z4;
	z1 = 0.5*z1;
	z2 = 0.5*z2;

	//{{z11, z8, z3}, 
	// {z8, z5, z1}, 
	// {z3, z1, z2}}

	double* p = rank4_ab.Pointer();
	*p++ = z11;
    *p++ = z8;
    *p++ = z3;
    *p++ = z8;
    *p++ = z5;
    *p++ = z1;
    *p++ = z3;
    *p++ = z1;
    *p   = z2;
}

void FSFiberMatViscT::MixedRank4_3D(const dArrayT& a, const dArrayT& b, 
	dMatrixT& rank4_ab) const
{
#if __option(extended_errorcheck)
	if (a.Length() != 3 ||
	    b.Length() != 3 ||
	    rank4_ab.Rows() != 6 ||
	    rank4_ab.Cols() != 6) throw ExceptionT::kSizeMismatch;
#endif

	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24;
	double z25, z26, z27, z28, z29, z30, z31, z32, z33, z34, z35, z36;
	double z37, z38, z39, z40, z41;

//	z1 = A(1.);
//	z2 = A(2.);
//	z3 = A(3.);
//	z4 = B(1.);
//	z5 = B(2.);
//	z6 = B(3.);

	z1 = a[0];	
	z2 = a[1];
	z3 = a[2];
	z4 = b[0];
	z5 = b[1];
	z6 = b[2];
	z7 = z1*z1;
	z8 = z2*z2;
	z9 = z3*z3;
	z10 = 2.*z1*z4;
	z11 = z10*z2;
	z12 = z2*z4;
	z13 = z1*z12;
	z14 = z4*z4;
	z15 = z10*z5;
	z15 = z15*z3;
	z16 = z11*z5;
	z17 = z12*z5;
	z18 = 2.*z17;
	z17 = z17*z3;
	z18 = z18*z3;
	z19 = z1*z4*z5;
	z19 = z19*z3;
	z20 = z5*z5;
	z11 = z11*z6;
	z13 = z13*z6;
	z10 = z10*z3*z6;
	z12 = z12*z3*z6;
	z21 = 2.*z12;
	z22 = z1*z2*z5*z6;
	z23 = 2.*z22;
	z24 = z1*z3*z5*z6;
	z25 = 2.*z24;
	z26 = 2.*z2*z3*z5*z6;
	z27 = z6*z6;
	z28 = 2.*z1*z14;
	z29 = 2.*z1*z2;
	z30 = z14*z2;
	z11 = z11 + z15;
	z15 = z1*z20*z3;
	z31 = 2.*z2*z20*z3;
	z18 = z18 + z23;
	z21 = z21 + z25;
	z23 = z1*z2*z27;
	z1 = 2.*z1*z27*z3;
	z25 = 2.*z2*z27*z3;
	z2 = z2*z28;
	z28 = z28*z3;
	z29 = z20*z29;
	z3 = z3*z30;
	z30 = 2.*z14*z7;
	z32 = z20*z7;
	z33 = z27*z7;
	z34 = 2.*z4*z5*z7;
	z35 = 2.*z4*z6*z7;
	z7 = z5*z6*z7;
	z36 = z14*z8;
	z37 = 2.*z20*z8;
	z38 = z27*z8;
	z39 = 2.*z4*z5*z8;
	z40 = z4*z6*z8;
	z8 = 2.*z5*z6*z8;
	z14 = z14*z9;
	z20 = z20*z9;
	z27 = 2.*z27*z9;
	z41 = z4*z5*z9;
	z4 = 2.*z4*z6*z9;
	z5 = 2.*z5*z6*z9;
	z6 = 0.5*z11;
	z9 = 0.5*z18;
	z11 = 0.5*z21;
	z2 = z2 + z34;
	z18 = z28 + z35;
	z3 = z13 + z19 + z3 + z7;
	z7 = z16 + z32 + z36;
	z13 = z29 + z39;
	z15 = z15 + z17 + z22 + z40;
	z8 = z31 + z8;
	z14 = z10 + z14 + z33;
	z17 = z20 + z26 + z38;
	z12 = z12 + z23 + z24 + z41;
	z1 = z1 + z4;
	z4 = z25 + z5;
	z2 = 0.5*z2;
	z5 = 0.5*z18;
	z3 = 0.5*z3;
	z7 = 0.5*z7;
	z13 = 0.5*z13;
	z15 = 0.5*z15;
	z8 = 0.5*z8;
	z14 = 0.5*z14;
	z17 = 0.5*z17;
	z12 = 0.5*z12;
	z1 = 0.5*z1;
	z4 = 0.5*z4;
	
	//{{z30, z16, z10,  z6,  z5,  z2}, 
	// {z16, z37, z26,  z8,  z9, z13}, 
	// {z10, z26, z27,  z4,  z1, z11}, 
	// { z6,  z8,  z4, z17, z12, z15}, 
	// { z5,  z9,  z1, z12, z14,  z3},
	// { z2, z13, z11, z15,  z3,  z7}}
	
	double* p = rank4_ab.Pointer();
    *p++ = z30;
    *p++ = z16;
    *p++ = z10;
    *p++ = z6;
    *p++ = z5;
    *p++ = z2;
    *p++ = z16;
    *p++ = z37;
    *p++ = z26;
    *p++ = z8;
    *p++ = z9;
    *p++ = z13;
    *p++ = z10;
    *p++ = z26;
    *p++ = z27;
    *p++ = z4;
    *p++ = z1;
    *p++ = z11;
    *p++ = z6;
    *p++ = z8;
    *p++ = z4;
    *p++ = z17;
    *p++ = z12;
    *p++ = z15;
    *p++ = z5;
    *p++ = z9;
    *p++ = z1;
    *p++ = z12;
    *p++ = z14;
    *p++ = z3;
    *p++ = z2;
    *p++ = z13;
    *p++ = z11;
    *p++ = z15;
    *p++ = z3;
    *p  = z7;
}
