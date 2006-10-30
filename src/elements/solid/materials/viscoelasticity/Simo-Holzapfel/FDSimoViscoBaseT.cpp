/* $Id: FDSimoViscoBaseT.cpp,v 1.1 2006-10-30 23:32:06 thao Exp $ */
/* created:   TDN (5/31/2001) */
#include "FDSimoViscoBaseT.h"

#include "fstreamT.h"
#include "ExceptionT.h"

using namespace Tahoe;

const int kNumOutputVar = 2;
static const char* Labels[kNumOutputVar] = {"r_dil","r_dev"};

FDSimoViscoBaseT::FDSimoViscoBaseT(ifstreamT& in,  
				   const FSMatSupportT& support):
	ParameterInterfaceT("FDSimoViscoBaseT")
//	FSSolidMatT(in, support)
{
ExceptionT::GeneralFail("FDSimoViscoBaseT::FDSimoViscoBaseT", "out of date");
#if 0
	int nsd = NumSD();
        int numstress = (nsd*(nsd+1))/2;
		
	/*allocates storage for state variable array*/
	fnstatev = 0;
	fnstatev += numstress;   /*current deviatoric overstress*/
	fnstatev += numstress;   /*current deviatoric inelastic stress*/
	fnstatev ++;			 /*current mean over stress*/
	fnstatev ++; 			 /*current mean inelastic stress*/
	
	fnstatev += numstress;   /*preceding deviatoric overstress*/
	fnstatev += numstress;   /*preceding deviatoric inelastic stress*/
	fnstatev ++;			 /*preceding mean overstress*/
	fnstatev ++; 			 /*preceding mean inelastic stress*/
	
	fstatev.Dimension(fnstatev);
	double* pstatev = fstatev.Pointer();
	/* assign pointers to current and preceding blocks of state variable array */
	
	fdevQ.Set(nsd, pstatev);        
	pstatev += numstress;
	fdevSin.Set(nsd, pstatev);	
	pstatev += numstress;
	fmeanQ.Set(1, pstatev);
	pstatev ++;
	fmeanSin.Set(1, pstatev);
	pstatev ++;
	
	fdevQ_n.Set(nsd, pstatev);        
	pstatev += numstress;
	fdevSin_n.Set(nsd, pstatev);	
	pstatev += numstress;
	fmeanQ_n.Set(1, pstatev);
	pstatev ++;
	fmeanSin_n.Set(1, pstatev);
#endif
}	

void FDSimoViscoBaseT::PointInitialize(void)
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
 
void FDSimoViscoBaseT::UpdateHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables */
		Load(element, ip);
	
		/* assign "current" to "last" */	
		fdevQ_n = fdevQ;
		fdevSin_n = fdevSin;
		fmeanQ_n = fmeanQ;
		fmeanSin_n = fmeanSin;

		/* write to storage */
		Store(element, ip);
	}
}

void FDSimoViscoBaseT::ResetHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load state variables*/
		Load(element, ip);
	
		/* assign "last" to "current" */
		fdevQ = fdevQ_n;
		fdevSin = fdevSin_n;
		fmeanQ = fmeanQ_n;
		fmeanSin = fmeanSin_n;
		
		/* write to storage */
		Store(element, ip);
	}
}

void FDSimoViscoBaseT::Load(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pd = d_array.Pointer(fnstatev*ip);
	double* pdr = fstatev.Pointer();
	for (int i = 0; i < fnstatev; i++)
		*pdr++ = *pd++;
}

void FDSimoViscoBaseT::Store(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pdr = fstatev.Pointer();
	double* pd = d_array.Pointer(fnstatev*ip);
	for (int i = 0; i < fnstatev; i++)
		*pd++ = *pdr++;
}

int FDSimoViscoBaseT::NumOutputVariables() const {return kNumOutputVar;}

void FDSimoViscoBaseT::OutputLabels(ArrayT<StringT>& labels) const
{
	//allocates space for labels
	labels.Dimension(kNumOutputVar);
	
	//copy labels
	for (int i = 0; i< kNumOutputVar; i++)
		labels[i] = Labels[i];
}
	
void FDSimoViscoBaseT::ComputeOutput(dArrayT& output)
{
  // obtains ratio of elastic strain to total strain
	
	ElementCardT& element = CurrentElement();
	Load(element, CurrIP());
	
	output[0] = fmeanQ[0];

	/*equivalent deviatoric stress*/
	double third = 1.0/3.0;
	output[1] = sqrt(2.0*third*fdevQ.ScalarProduct());
}	

double FDSimoViscoBaseT::StrainEnergyDensity(void)
{
        cout << "\nStrain Energy is undefined.";
	throw ExceptionT::kGeneralFail;
}
