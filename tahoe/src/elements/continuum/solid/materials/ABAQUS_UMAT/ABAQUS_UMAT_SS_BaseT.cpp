/* $Id: ABAQUS_UMAT_SS_BaseT.cpp,v 1.5 2004-06-17 07:40:11 paklein Exp $ */
#include "ABAQUS_UMAT_SS_BaseT.h"

#ifdef __F2C__

#include <ctype.h>
#include <float.h>

#include "ifstreamT.h"
#include "ContinuumElementT.h" //needed for ip coordinates
#include "SSMatSupportT.h"

#include "ThermalDilatationT.h"

using namespace Tahoe;

/* debugging */
//#define DEBUG 1
#undef DEBUG

/* constructor */
ABAQUS_UMAT_SS_BaseT::ABAQUS_UMAT_SS_BaseT(ifstreamT& in, const SSMatSupportT& support):
	SSSolidMatT(in, support),
	fApproxModulus(false),
	fNumElasticIterations(-1),
	fTangentType(GlobalT::kSymmetric),
	fStress(NumSD()),
	fIPCoordinates(NumSD()),
	fPressure(0.0),
	fmat_nsd(NumSD()),
	fsym_mat_nsd(dSymMatrixT::int2DimensionT(NumSD()))
{
	const char caller[] = "ABAQUS_UMAT_SS_BaseT::ABAQUS_UMAT_SS_BaseT";

	/* UMAT dimensions */
	ndi = 0;
	nshr = 0;
	ntens = 0;
	nstatv = 0;

	/* read non-ABAQUS parameters */
	in >> fNumElasticIterations;

	/* read ABAQUS-format input */
	bool nonsym = false;	
	Read_ABAQUS_Input(in, fUMAT_name, fProperties, fDensity, nstatv, nonsym);
	if (nonsym)
		fTangentType = GlobalT::kNonSymmetric;

#if DEBUG
StringT UMAT_file;
UMAT_file.Root(in.filename());
UMAT_file.Append(".UMAT.log");
flog.open(UMAT_file);
flog.precision(DBL_DIG);
flog.setf(ios::showpoint);
flog.setf(ios::right, ios::adjustfield);
flog.setf(ios::scientific, ios::floatfield);
#endif
}

/* print parameters */
void ABAQUS_UMAT_SS_BaseT::Print(ostream& out) const
{
	/* inherited */
	SSSolidMatT::Print(out);
	IsotropicT::Print(out);
	
	/* class parameters */
	out << " Number of elastic iterations. . . . . . . . . . = " << fNumElasticIterations << '\n';
	
	/* write properties array */
	out << " Number of ABAQUS UMAT internal variables. . . . = " << nstatv << '\n';
	out << " Number of ABAQUS UMAT properties. . . . . . . . = " << fProperties.Length() << '\n';
	PrintProperties(out);
}

/* disable multiplicative thermal strains */
void ABAQUS_UMAT_SS_BaseT::Initialize(void)
{
	const char caller[] = "ABAQUS_UMAT_SS_BaseT::Initialize";

	/* inherited */
	SSSolidMatT::Initialize();

	/* notify */
	if (fThermal->IsActive())
		cout << "\n ABAQUS_UMAT_SS_BaseT::Initialize: thermal strains must\n"
		     <<   "    be handled within the UMAT\n" << endl;
	
	/* UMAT dimensions */
	ndi = 3; // always 3 direct components
	int nsd = NumSD();
	if (nsd == 2)
		nshr = 1;
	else if (nsd == 3)
		nshr = 3;
	else
		ExceptionT::GeneralFail(caller);
	ntens = ndi + nshr;

	/* modulus storage */
	if (fTangentType == GlobalT::kDiagonal)
		fModulusDim = ntens;
	else if (fTangentType == GlobalT::kSymmetric)
	{
		if (nsd == 2) fModulusDim = 10;
		else if (nsd == 3) fModulusDim = 21;
		else ExceptionT::GeneralFail(caller);
	}
	else if (fTangentType == GlobalT::kNonSymmetric)
		fModulusDim = ntens*ntens;
	else
		ExceptionT::GeneralFail(caller);

	/* storage block size (per ip) */
	fBlockSize = 0;
	fBlockSize += ntens;       // fstress
	fBlockSize += ntens;       // fstrain
	fBlockSize += 3;           // fsse_pd_cd
	fBlockSize += nstatv;      // fstatv
	fBlockSize += fModulusDim; // fmodulus
	fBlockSize += ntens;       // fstress_last
	fBlockSize += ntens;       // fstrain_last
	fBlockSize += 3;           // fsse_pd_cd_last
	fBlockSize += nstatv;      // fstatv_last
	
	/* argument array */
	fArgsArray.Dimension(fBlockSize);

	/* assign pointers */
	doublereal* parg = fArgsArray.Pointer();
	fstress.Set(ntens, parg);        parg += ntens;
	fstrain.Set(ntens, parg);        parg += ntens;
	fsse_pd_cd.Set(3, parg);         parg += 3;
	fstatv.Set(nstatv, parg);        parg += nstatv;
	fmodulus.Set(fModulusDim, parg); parg += fModulusDim;
	fstress_last.Set(ntens, parg);   parg += ntens;
	fstrain_last.Set(ntens, parg);   parg += ntens;
	fsse_pd_cd_last.Set(3, parg);    parg += 3;
	fstatv_last.Set(nstatv, parg);

	/* UMAT array arguments */
	fddsdde.Dimension(ntens);
	fddsdde = 0.0;
	fdstran.Dimension(ntens);
	fdstran = 0.0;
	fdrot.Dimension(3);   // always 3
	fdrot.Identity();
	fdfgrd0.Dimension(3); // always 3
	fdfgrd0.Identity();
	fdfgrd1.Dimension(3); // always 3
	fdfgrd1.Identity();
	fcoords.Dimension(nsd);
}

/* materials initialization */
bool ABAQUS_UMAT_SS_BaseT::NeedsPointInitialization(void) const { return true; }
void ABAQUS_UMAT_SS_BaseT::PointInitialize(void)
{
	/* allocate element storage */
	if (CurrIP() == 0)
	{
		ElementCardT& element = CurrentElement();
		element.Dimension(0, fBlockSize*NumIP());
	
		/* initialize */
		element.DoubleData() = 0.0;
	}

	/* call UMAT - time signals initialization */
	Call_UMAT(0.0, 0.0, 0, 0);

	/* store results as last converged */
	if (CurrIP() == NumIP() - 1) UpdateHistory();
}

/* update/reset internal variables */
void ABAQUS_UMAT_SS_BaseT::UpdateHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load stored data */
		Load(element, ip);
	
		/* assign "current" to "last" */	
		fstress_last    = fstress;
		fstrain_last    = fstrain;
		fsse_pd_cd_last = fsse_pd_cd;
		fstatv_last     = fstatv;

		/* write to storage */
		Store(element, ip);
	}
}

void ABAQUS_UMAT_SS_BaseT::ResetHistory(void)
{
	/* current element */
	ElementCardT& element = CurrentElement();	
	for (int ip = 0; ip < NumIP(); ip++)
	{
		/* load stored data */
		Load(element, ip);
	
		/* assign "last" to "current" */
		fstress    = fstress_last;
		fstrain    = fstrain_last;
		fsse_pd_cd = fsse_pd_cd_last;
		fstatv     = fstatv_last;

		/* write to storage */
		Store(element, ip);
	}
}

/* spatial description */
const dMatrixT& ABAQUS_UMAT_SS_BaseT::c_ijkl(void)
{
	const char caller[] = "ABAQUS_UMAT_SS_BaseT::c_ijkl";

	/* load stored data */
	Load(CurrentElement(), CurrIP());

	int nsd = NumSD();

#if __option(extended_errorcheck)
	if (nsd == 2) {
		if (ntens != 4)	
			ExceptionT::SizeMismatch(caller);
	} else if (nsd == 3) {
		if (ntens != 6)	
			ExceptionT::SizeMismatch(caller);
	} else
		ExceptionT::GeneralFail(caller);
#endif

	if (fTangentType == GlobalT::kDiagonal)
	{
		if (nsd == 2)
		{
			fModulus(0,0) = double(fmodulus[0]); // 11
			fModulus(1,1) = double(fmodulus[1]); // 22
			fModulus(2,2) = double(fmodulus[3]); // 12
		}
		else
		{
			fModulus(0,0) = double(fmodulus[0]); // 11
			fModulus(1,1) = double(fmodulus[1]); // 22
			fModulus(2,2) = double(fmodulus[2]); // 33
			fModulus(3,3) = double(fmodulus[5]); // 23
			fModulus(4,4) = double(fmodulus[4]); // 13
			fModulus(5,5) = double(fmodulus[3]); // 12
		}
	}
	else if (fTangentType == GlobalT::kSymmetric)
	{
		if (nsd == 2)
		{
			fModulus(0,0) = double(fmodulus[0]);

			fModulus(1,0) = fModulus(0,1) = double(fmodulus[1]);
			fModulus(1,1) = double(fmodulus[2]);

			fModulus(2,0) = fModulus(0,2) = double(fmodulus[6]);
			fModulus(2,1) = fModulus(1,2) = double(fmodulus[7]);
			fModulus(2,2) = double(fmodulus[9]);
		}
		else
		{
			fModulus(0,0) = double(fmodulus[0]);

			fModulus(1,0) = fModulus(0,1) = double(fmodulus[1]);
			fModulus(1,1) = double(fmodulus[2]);

			fModulus(2,0) = fModulus(0,2) = double(fmodulus[3]);
			fModulus(2,1) = fModulus(1,2) = double(fmodulus[4]);
			fModulus(2,2) = double(fmodulus[5]);

			fModulus(3,0) = fModulus(0,3) = double(fmodulus[15]);
			fModulus(3,1) = fModulus(1,3) = double(fmodulus[16]);
			fModulus(3,2) = fModulus(2,3) = double(fmodulus[17]);
			fModulus(3,3) = double(fmodulus[20]);
		
			fModulus(4,0) = fModulus(0,4) = double(fmodulus[10]);
			fModulus(4,1) = fModulus(1,4) = double(fmodulus[11]);
			fModulus(4,2) = fModulus(2,4) = double(fmodulus[12]);
			fModulus(4,3) = fModulus(3,4) = double(fmodulus[19]);
			fModulus(4,4) = double(fmodulus[14]);
		
			fModulus(5,0) = fModulus(0,5) = double(fmodulus[6]);
			fModulus(5,1) = fModulus(1,5) = double(fmodulus[7]);
			fModulus(5,2) = fModulus(2,5) = double(fmodulus[8]);
			fModulus(5,3) = fModulus(3,5) = double(fmodulus[18]);
			fModulus(5,4) = fModulus(4,5) = double(fmodulus[13]);		
			fModulus(5,5) = double(fmodulus[9]);
		}
	}
	else if (fTangentType == GlobalT::kNonSymmetric)
	{
		if (nsd == 2) 
		{
			double* mod_tahoe = fModulus.Pointer();
			*mod_tahoe++ = double(fmodulus[0]);
			*mod_tahoe++ = double(fmodulus[1]);
			*mod_tahoe++ = double(fmodulus[3]);

			*mod_tahoe++ = double(fmodulus[4]);
			*mod_tahoe++ = double(fmodulus[5]);
			*mod_tahoe++ = double(fmodulus[7]);

			*mod_tahoe++ = double(fmodulus[12]);
			*mod_tahoe++ = double(fmodulus[13]);
			*mod_tahoe   = double(fmodulus[15]);
		}
		else
		{
			int tahoe2abaqus[6] = {0,1,2,5,4,3};
			double* mod_tahoe = fModulus.Pointer();
			for (int i = 0; i < 6; i++)
			{
				doublereal* mod_abaqus = fmodulus.Pointer(ntens*tahoe2abaqus[i]);
				for (int j = 0; j < 6; j++)
					*mod_tahoe++ = double(mod_abaqus[tahoe2abaqus[j]]);
			}
		}
	}
	else 
		ExceptionT::GeneralFail(caller);

#if DEBUG
flog << setw(10) << "element: " << MaterialSupport().CurrElementNumber()+1 << '\n';
flog << setw(10) << "     ip: " << CurrIP()+1 << '\n';
flog << setw(10) << "modulus:\n" << fModulus << endl;
#endif

	return fModulus;
}

const dSymMatrixT& ABAQUS_UMAT_SS_BaseT::s_ij(void)
{
	/* call UMAT */
	if (MaterialSupport().RunState() == GlobalT::kFormRHS ||
	    (fApproxModulus && MaterialSupport().RunState() == GlobalT::kFormLHS))
	{
		double  t = fSSMatSupport->Time();
		double dt = fSSMatSupport->TimeStep();
		int  step = fSSMatSupport->StepNumber();
		int  iter = fSSMatSupport->IterationNumber();
		if (iter < fNumElasticIterations-1) /* treat as 'first' iteration */
			iter = -1;
		
		Call_UMAT(t, dt, step, iter);
	}
	else
		/* load stored data */
		Load(CurrentElement(), CurrIP());

	/* copy/convert stress */
	ABAQUS_to_dSymMatrixT(fstress.Pointer(), fStress);
	fPressure = fStress.Trace()/3.0;
	return fStress;
}

/* returns the strain energy density for the specified strain */
double ABAQUS_UMAT_SS_BaseT::StrainEnergyDensity(void)
{
	/* load stored data */
	Load(CurrentElement(), CurrIP());

	/* pull from storage */
	return double(fsse_pd_cd[0]);
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int ABAQUS_UMAT_SS_BaseT::NumOutputVariables(void) const
{
	/* set material output variables/labels */
	if (fOutputIndex.Length() == 0)
	{
		//TEMP - better place for this?
		ABAQUS_UMAT_SS_BaseT* tmp = (ABAQUS_UMAT_SS_BaseT*) this;
		tmp->SetOutputVariables(tmp->fOutputIndex, tmp->fOutputLabels);
	}

	return fOutputIndex.Length();
}

void ABAQUS_UMAT_SS_BaseT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(fOutputLabels.Length());
	for (int i = 0; i < labels.Length(); i++)
		labels[i] = fOutputLabels[i];
}

void ABAQUS_UMAT_SS_BaseT::ComputeOutput(dArrayT& output)
{
	/* check */
	if (output.Length() != fOutputIndex.Length())
	{
		cout << "\n ABAQUS_UMAT_SS_BaseT::ComputeOutput: not enough space to return\n"
		     <<   "     output variables: given " << output.Length()
		     << ". expecting " << fOutputIndex.Length() << "." << endl;
		throw ExceptionT::kSizeMismatch;
	}

	/* load stored data */
	Load(CurrentElement(), CurrIP());

	/* collect variables */
	for (int i = 0; i < fOutputIndex.Length(); i++)
		output[i] = double(fstatv[fOutputIndex[i]]);
}

/***********************************************************************
* Protected
***********************************************************************/

/* I/O functions */
void ABAQUS_UMAT_SS_BaseT::PrintName(ostream& out) const
{
	/* inherited */
	SSSolidMatT::PrintName(out);
	out << "    ABAQUS user material: " << fUMAT_name << '\n';
}

void ABAQUS_UMAT_SS_BaseT::PrintProperties(ostream& out) const
{
	/* just write numbered list */
	int d_width = OutputWidth(out, fProperties.Pointer());
	for (int i = 0; i < fProperties.Length(); i++)	
		out << setw(kIntWidth) << i+1
		    << setw(  d_width) << fProperties[i] << '\n';
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* load element data for the specified integration point */
void ABAQUS_UMAT_SS_BaseT::Load(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pd = d_array.Pointer(fBlockSize*ip);
	doublereal* pdr = fArgsArray.Pointer();
	for (int i = 0; i < fBlockSize; i++)
		*pdr++ = doublereal(*pd++);
}

void ABAQUS_UMAT_SS_BaseT::Store(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	doublereal* pdr = fArgsArray.Pointer();
	double* pd = d_array.Pointer(fBlockSize*ip);
	for (int i = 0; i < fBlockSize; i++)
		*pd++ = double(*pdr++);
}

/* make call to the UMAT */
void ABAQUS_UMAT_SS_BaseT::Call_UMAT(double t, double dt, int step, int iter)
{
	/* load stored data */
	Load(CurrentElement(), CurrIP());

	/* set stored variables to values at beginning of increment */
	Reset_UMAT_Increment();

	/* compute strain/rotated stress */
	Set_UMAT_Arguments();

	/* map UMAT arguments */
	doublereal* stress = fstress.Pointer();             // i/o: Cauchy stress
	doublereal* statev = fstatv.Pointer();              // i/o: state variables
	doublereal* ddsdde = fddsdde.Pointer();             //   o: constitutive Jacobian
	doublereal  sse = fsse_pd_cd[0];                    // i/o: specific elastic strain energy
	doublereal  spd = fsse_pd_cd[1];                    // i/o: plastic dissipation
	doublereal  scd = fsse_pd_cd[2];                    // i/o: creep dissipation

	// for fully-coupled only
	doublereal  rpl;                                    // o: volumetric heat generation
	doublereal* ddsddt = NULL;                          // o: stress-temperature variation
	doublereal* drplde = NULL;                          // o: rpl-strain variation
	doublereal  drpldt;                                 // o: rpl-temperature variation

	doublereal* stran  = fstrain.Pointer();             // i: total integrated strain
	doublereal* dstran = fdstran.Pointer();             // i: strain increment
	doublereal  time[2];                                // i: {step time, total time} at the beginning of increment
	time[0] = time[1]  = doublereal(t);
	doublereal  dtime  = doublereal(dt);                // i: time step
	doublereal  temp   = 0.0;                           // i: temperature at start
	doublereal  dtemp  = 0.0;                           // i: temperature increment
	doublereal* predef = NULL;                          // i: pre-defined field variables
	doublereal* dpred  = NULL;                          // i: increment of pre-defined field variables
	char*       cmname = fUMAT_name.Pointer();          // i: UMAT name
	doublereal* props  = fProperties.Pointer();         // i: material properties array
	integer     nprops = integer(fProperties.Length()); // i: number of material properties
	doublereal* coords = fcoords.Pointer();             // i: coordinates of the integration point
	doublereal* drot   = fdrot.Pointer();               // i: rotation increment matrix
	doublereal  pnewdt = dtime;                         // o: suggested time step (automatic time integration)
	doublereal  celent;                                 // i: characteristic element length
	doublereal* dfgrd0 = fdfgrd0.Pointer();             // i: deformation gradient at the beginning of the increment
	doublereal* dfgrd1 = fdfgrd1.Pointer();             // i: deformation gradient at the end of the increment
	integer     noel   = integer(CurrElementNumber());  // i: element number
	integer     npt    = integer(CurrIP());             // i: integration point number
	integer     layer  = 1;                             // i: layer number (composites/layered solids)
	integer     kspt   = 1;                             // i: section point
	integer     kstep  = integer(step);                 // i: step number
	integer     kinc   = integer(iter);                 // i: increment number
	ftnlen      cmname_len = strlen(fUMAT_name);      // f2c: length of cmname string

#ifdef DEBUG
if (MaterialSupport().RunState() == GlobalT::kFormRHS) {
	int d_width = OutputWidth(flog, stress);
	flog << "\n THE INPUT\n";
	flog << setw(10) << "   time: " << setw(d_width) << time[0]  << '\n';
	flog << setw(10) << "   iter: " << MaterialSupport().IterationNumber() << '\n';
	flog << setw(10) << "element: " << MaterialSupport().CurrElementNumber()+1 << '\n';
	flog << setw(10) << "     ip: " << CurrIP()+1 << '\n';
	flog << setw(10) << " stress: " << fstress.no_wrap() << '\n';
	flog << setw(10) << " strain: " << fstrain.no_wrap() << '\n';
	flog << setw(10) << "dstrain: " << fdstran.no_wrap() << '\n';
	flog << setw(10) << "  state:\n";
	flog << fstatv.wrap(5) << '\n';
}
#endif

	/* call UMAT wrapper */
	UMAT(stress, statev, ddsdde, &sse, &spd, &scd, &rpl, ddsddt, drplde,
		&drpldt, stran, dstran, time, &dtime, &temp, &dtemp, predef, dpred, cmname,
		&ndi, &nshr, &ntens, &nstatv, props, &nprops, coords, drot, &pnewdt, &celent,
		dfgrd0, dfgrd1, &noel, &npt, &layer, &kspt, &kstep, &kinc, cmname_len);

	/* check for step cut */
	if (fabs(dtime) > kSmall && pnewdt/dtime < 0.55) {
	
		bool re_run_UMAT = false;
		if (re_run_UMAT)
			Call_UMAT(t, dt, step, iter);
	
		ExceptionT::BadJacobianDet("ABAQUS_UMAT_SS_BaseT::Call_UMAT", "material signaled step cut");
	}

#ifdef DEBUG
if (MaterialSupport().RunState() == GlobalT::kFormRHS) {
	flog << " THE OUTPUT\n";
	flog << setw(10) << " stress: " << fstress.no_wrap() << '\n';
	flog << setw(10) << " state:\n";
	flog << fstatv.wrap(5) << endl;
}
#endif

	/* update strain */
	fstrain += fdstran;

	/* store modulus */
	Store_UMAT_Modulus();
	
	/* write to storage */
	Store(CurrentElement(), CurrIP());
}

/* set variables to last converged */
void ABAQUS_UMAT_SS_BaseT::Reset_UMAT_Increment(void)
{
	/* assign "last" to "current" */
	fstress    = fstress_last;
	fstrain    = fstrain_last;
	fsse_pd_cd = fsse_pd_cd_last;
	fstatv     = fstatv_last;
}

/* set stress/strain arguments */
void ABAQUS_UMAT_SS_BaseT::Set_UMAT_Arguments(void)
{
	/* integration point coordinates */
	ContinuumElement().IP_Coords(fIPCoordinates);	
	fcoords[0] = doublereal(fIPCoordinates[0]);
	fcoords[1] = doublereal(fIPCoordinates[1]);
	if (NumSD() == 3)
		fcoords[2] = doublereal(fIPCoordinates[2]);

	/* deformation gradients */
	fmat_nsd.Identity();
	dMatrixT_to_ABAQUS(fmat_nsd, fdfgrd0);
	dMatrixT_to_ABAQUS(fmat_nsd, fdfgrd1);

	/* incremental rotation */
	dMatrixT_to_ABAQUS(fmat_nsd, fdrot);

	/* total integrated strain */
	const dSymMatrixT& e_int = e_last();
	dSymMatrixT_to_ABAQUS(e_int, fstrain.Pointer(), true);

	/* incremental strain */
	fsym_mat_nsd.DiffOf(e(), e_int);
	dSymMatrixT_to_ABAQUS(fsym_mat_nsd, fdstran.Pointer(), true);

	/* initialize Jacobian */
	fddsdde = 0.0;
}

/* store the modulus */
void ABAQUS_UMAT_SS_BaseT::Store_UMAT_Modulus(void)
{
	if (fTangentType == GlobalT::kDiagonal)
	{
		/* take diagonal values */
		for (int i = 0; i < fmodulus.Length(); i++)
			fmodulus[i] = fddsdde(i,i);
	}
	else if (fTangentType == GlobalT::kSymmetric)
	{
		/* columns */
		int dex = 0;
		for (int j = 0; j < ntens; j++)
			for (int i = 0; i <= j; i++)
				if (i == j)
					fmodulus[dex++] = fddsdde(i,j);
				else
					fmodulus[dex++] = 0.5*(fddsdde(i,j) + fddsdde(j,i));
	}
	else if (fTangentType == GlobalT::kNonSymmetric)
	{
		/* store everything */
		fmodulus = fddsdde;
	}
	else 
		ExceptionT::GeneralFail("ABAQUS_UMAT_SS_BaseT::Store_UMAT_Modulus");	
}

#endif /* __F2C__ */
