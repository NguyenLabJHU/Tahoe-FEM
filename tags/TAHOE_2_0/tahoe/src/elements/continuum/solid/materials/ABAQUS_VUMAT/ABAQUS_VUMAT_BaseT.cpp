/* $Id: ABAQUS_VUMAT_BaseT.cpp,v 1.25 2004-07-15 08:26:37 paklein Exp $ */
#include "ABAQUS_VUMAT_BaseT.h"

#ifdef __F2C__

#include <ctype.h>
#include <float.h>


#include "ContinuumElementT.h"

#include "SpectralDecompT.h"
#include "ThermalDilatationT.h"

#define VUMAT_DEBUG 0

using namespace Tahoe;

/* constructor */
ABAQUS_VUMAT_BaseT::ABAQUS_VUMAT_BaseT(ifstreamT& in, const FSMatSupportT& support):
	FSSolidMatT(in, support),
	fTangentType(GlobalT::kSymmetric),
	fModulus(dSymMatrixT::NumValues(NumSD())),
	fStress(NumSD()),
	fIPCoordinates(NumSD()),
	fPressure(0.0),
	fDecomp(NULL),
	fF_rel(NumSD()),
	fROld(NumSD()),
	fRNew(NumSD()),
	fA_nsd(NumSD()),
	fU1(NumSD()), fU2(NumSD()), fU1U2(NumSD()), fUOld(NumSD()), fUNew(NumSD())
{
	/* read ABAQUS-format input */
	nstatv = 0;
	bool nonsym = false;
	Read_ABAQUS_Input(in, fVUMAT_name, fProperties, fDensity, nstatv, nonsym);

	/* set (material tangent) modulus tensor (fixed) */
	double   Young = double(fProperties[0]);
	double Poisson = double(fProperties[1]);
	Set_E_nu(Young, Poisson);
		
	/* VUMAT dimensions */
	ndi = 3; // always 3 direct components
	int nsd = NumSD();
	if (nsd == 2)
		nshr = 1;
	else if (nsd == 3)
		nshr = 3;
	else
		throw ExceptionT::kGeneralFail;
	ntens = ndi + nshr;

	/* modulus storage */
	if (fTangentType == GlobalT::kDiagonal)
		fModulusDim = ntens;
	else if (fTangentType == GlobalT::kSymmetric)
	{
		if (nsd == 2) fModulusDim = 10;
		else if (nsd == 3) fModulusDim = 21;
		else throw ExceptionT::kGeneralFail;
	}
	else if (fTangentType == GlobalT::kNonSymmetric)
		fModulusDim = ntens*ntens;
	else
		throw ExceptionT::kGeneralFail;

	/* storage block size (per ip) */
	fBlockSize = 0;
	fBlockSize += ntens;       // fstress
	fBlockSize += ntens;       // fstrain
	//fBlockSize += 3;           // fsse_pd_cd
	/* may not need that (above) */
	fBlockSize += nstatv;      // fstatv
	//fBlockSize += fModulusDim; // fmodulus
	fBlockSize += ntens;       // fstress_last
	fBlockSize += ntens;       // fstrain_last
	//fBlockSize += 3;           // fsse_pd_cd_last
	/* may not need this (above) */
	fBlockSize += nstatv;      // fstatv_last
	
	/* argument array */
	fArgsArray.Dimension(fBlockSize);

	/* assign pointers */
	doublereal* parg = fArgsArray.Pointer();
	fstress.Set(ntens, parg);        parg += ntens;
	fstrain.Set(ntens, parg);        parg += ntens;
	//fsse_pd_cd.Set(3, parg);         parg += 3;
	fstatv.Set(nstatv, parg);        parg += nstatv;
	//fmodulus.Set(fModulusDim, parg); parg += fModulusDim;
	fstress_last.Set(ntens, parg);   parg += ntens;
	fstrain_last.Set(ntens, parg);   parg += ntens;
	//fsse_pd_cd_last.Set(3, parg);    parg += 3;
	fstatv_last.Set(nstatv, parg);
	
	/* VUMAT array arguments */
	//fddsdde.Dimension(ntens);
	fdstran.Dimension(ntens);
	fdstran = 0.0;
	fdrot.Dimension(3);   // always 3
	fdrot.Identity();
	fdfgrd0.Dimension(3); // always 3
	fdfgrd0.Identity();
	fdfgrd1.Dimension(3); // always 3
	fdfgrd1.Identity();
	fcoords.Dimension(nsd);

	/* initialize other VUMAT array arguments */
	fROld = 0.0;
	fRNew = 0.0;
	fRelSpin = 0.0;
	fUOld = 0.0;
	fUNew = 0.0;

	/* spectral decomp */
	fDecomp = new SpectralDecompT(NumSD());
	if (!fDecomp) throw ExceptionT::kOutOfMemory;

//DEBUG
#if VUMAT_DEBUG
StringT VUMAT_file;
VUMAT_file.Root(in.filename());
VUMAT_file.Append(".VUMAT.log");
flog.open(VUMAT_file);
flog.precision(DBL_DIG);
flog.setf(ios::showpoint);
flog.setf(ios::right, ios::adjustfield);
flog.setf(ios::scientific, ios::floatfield);
#endif
}

/* destructor */
ABAQUS_VUMAT_BaseT::~ABAQUS_VUMAT_BaseT(void)
{
	delete fDecomp;
	fDecomp = NULL;
}

/* print parameters */
void ABAQUS_VUMAT_BaseT::Print(ostream& out) const
{
	/* inherited */
	FSSolidMatT::Print(out);
	
	/* write properties array */
	out << " Number of ABAQUS VUMAT internal variables. . . . = " << nstatv << '\n';
	out << " Number of ABAQUS VUMAT properties. . . . . . . . = " << fProperties.Length() << '\n';
	PrintProperties(out);
}

/* disable multiplicative thermal strains */
void ABAQUS_VUMAT_BaseT::Initialize(void)
{
	/* inherited */
	FSSolidMatT::Initialize();

	/* notify */
	if (fThermal->IsActive())
		cout << "\n ABAQUS_VUMAT_BaseT::Initialize: thermal strains must\n"
		     <<   "    be handled within the VUMAT\n" << endl;
	
	/* disable thermal transform */
	//SetFmodMult(NULL);	
}


/* materials initialization */
bool ABAQUS_VUMAT_BaseT::NeedsPointInitialization(void) const { return true; }
void ABAQUS_VUMAT_BaseT::PointInitialize(void)
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
	double dt = fFSMatSupport->TimeStep();
	Call_VUMAT(0.0, dt, 0, 0);

	/* store results as last converged */
	if (CurrIP() == NumIP() - 1) UpdateHistory();
}

/* update/reset internal variables */
void ABAQUS_VUMAT_BaseT::UpdateHistory(void)
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
		//fsse_pd_cd_last = fsse_pd_cd;
		fstatv_last     = fstatv;

		/* write to storage */
		Store(element, ip);
	}
}

void ABAQUS_VUMAT_BaseT::ResetHistory(void)
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
		//fsse_pd_cd = fsse_pd_cd_last;
		fstatv     = fstatv_last;

		/* write to storage */
		Store(element, ip);
	}
}

/* spatial description */
const dMatrixT& ABAQUS_VUMAT_BaseT::c_ijkl(void)
{
	/* very approximate */
	return ABAQUS_VUMAT_BaseT::C_IJKL();
}

const dSymMatrixT& ABAQUS_VUMAT_BaseT::s_ij(void)
{
	/* call VUMAT */
	if (MaterialSupport().RunState() == GlobalT::kFormRHS)
	{
		double  t = fFSMatSupport->Time();
		double dt = fFSMatSupport->TimeStep();
		int  step = fFSMatSupport->StepNumber();
		int  iter = fFSMatSupport->IterationNumber();
		Call_VUMAT(t, dt, step, iter);
	}
	else
		/* load stored data */
		Load(CurrentElement(), CurrIP());

	/* copy/convert stress */
	ABAQUS_to_dSymMatrixT(fstress.Pointer(), fStress);
	return fStress;
}

/* material description */
const dMatrixT& ABAQUS_VUMAT_BaseT::C_IJKL(void)
{
	/* assuming constant and isotropic */
	if (NumSD() == 2)
		ComputeModuli2D(fModulus, Material2DT::kPlaneStrain);
	else
		ComputeModuli(fModulus);

	return fModulus;
}

const dSymMatrixT& ABAQUS_VUMAT_BaseT::S_IJ(void)
{
	/* Cauchy stress */
	const dSymMatrixT& s = ABAQUS_VUMAT_BaseT::s_ij();

	/* spatial -> material */
	fStress.SetToScaled(F().Det(), PullBack(F(), s));	
	return fStress;
}

/* returns the strain energy density for the specified strain */
double ABAQUS_VUMAT_BaseT::StrainEnergyDensity(void)
{
	/* load stored data */
	Load(CurrentElement(), CurrIP());

	/* pull from storage */
	//return double(fsse_pd_cd[0]);
	return 0.0;
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int ABAQUS_VUMAT_BaseT::NumOutputVariables(void) const
{
	/* set material output variables/labels */
	if (fOutputIndex.Length() == 0)
	{
		//TEMP - better place for this?
		ABAQUS_VUMAT_BaseT* tmp = (ABAQUS_VUMAT_BaseT*) this;
		tmp->SetOutputVariables(tmp->fOutputIndex, tmp->fOutputLabels);
	}

	return fOutputIndex.Length();
}

void ABAQUS_VUMAT_BaseT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Dimension(fOutputLabels.Length());
	for (int i = 0; i < labels.Length(); i++)
		labels[i] = fOutputLabels[i];
}

void ABAQUS_VUMAT_BaseT::ComputeOutput(dArrayT& output)
{
	/* check */
	if (output.Length() != fOutputIndex.Length())
	{
		cout << "\n ABAQUS_VUMAT_BaseT::ComputeOutput: not enough space to return\n"
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
void ABAQUS_VUMAT_BaseT::PrintName(ostream& out) const
{
	/* inherited */
	FSSolidMatT::PrintName(out);
	out << "    ABAQUS user material: " << fVUMAT_name << '\n';
}

void ABAQUS_VUMAT_BaseT::PrintProperties(ostream& out) const
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
void ABAQUS_VUMAT_BaseT::Load(const ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	const dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	const double* pd = d_array.Pointer(fBlockSize*ip);
	doublereal* pdr = fArgsArray.Pointer();
	for (int i = 0; i < fBlockSize; i++)
		*pdr++ = doublereal(*pd++);
}

void ABAQUS_VUMAT_BaseT::Store(ElementCardT& element, int ip)
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
void ABAQUS_VUMAT_BaseT::Call_VUMAT(double t, double dt, int step, int iter)
{	
	/* load stored data */
	Load(CurrentElement(), CurrIP());

	/* set stored variables to values at beginning of increment */
	Reset_VUMAT_Increment();

	/* compute strain/rotated stress */
	Set_VUMAT_Arguments();

	/* map VUMAT arguments */
	doublereal* stressold = fstress_last.Pointer();     // i: Cauchy stress - rotated
	doublereal* statevold = fstatv_last.Pointer();      // i: state variables
	//doublereal* ddsdde = fddsdde.Pointer();             //   o: constitutive Jacobian
	//doublereal  sse = fsse_pd_cd[0];                    // i/o: specific elastic strain energy
	//doublereal  spd = fsse_pd_cd[1];                    // i/o: plastic dissipation
	//doublereal  scd = fsse_pd_cd[2];                    // i/o: creep dissipation

	// for fully-coupled only
//	doublereal  rpl;                                    // o: volumetric heat generation
	doublereal* ddsddt = NULL;                          // o: stress-temperature variation
	doublereal* drplde = NULL;                          // o: rpl-strain variation
//	doublereal  drpldt;                                 // o: rpl-temperature variation

	doublereal* stran  = fstrain.Pointer();             // i: total integrated strain
	doublereal* dstran = fdstran.Pointer();             // i: strain increment
	doublereal  time[2];                                // i: {step time, total time} at the beginning of increment
	doublereal  stime = doublereal(t);
	doublereal  totime = doublereal(t);
	time[0] = time[1]  = doublereal(t);
	doublereal  dtime  = doublereal(dt);                // i: time step
	doublereal  temp   = 0.0;                           // i: temperature at start
	doublereal  dtemp  = 0.0;                           // i: temperature increment
	doublereal* predef = NULL;                          // i: pre-defined field variables
	doublereal* dpred  = NULL;                          // i: increment of pre-defined field variables
	char*       cmname = fVUMAT_name.Pointer();          // i: UMAT name
	doublereal* props  = fProperties.Pointer();         // i: material properties array
	integer     nprops = integer(fProperties.Length()); // i: number of material properties
	doublereal* coords = fcoords.Pointer();             // i: coordinates of the integration point
	doublereal* drot   = fdrot.Pointer();               // i: rotation increment matrix
//	doublereal  pnewdt;                                 // o: suggested time step (automatic time integration)
	doublereal  celent;                                 // i: characteristic element length
	doublereal* dfgrd0 = fdfgrd0.Pointer();             // i: deformation gradient at the beginning of the increment
	doublereal* dfgrd1 = fdfgrd1.Pointer();             // i: deformation gradient at the end of the increment
	integer     noel   = integer(CurrElementNumber());  // i: element number
	integer     npt    = integer(CurrIP());             // i: integration point number
	integer     layer  = 1;                             // i: layer number (composites/layered solids)
	integer     kspt   = 1;                             // i: section point
	integer     kstep  = integer(step);                 // i: step number
	integer     kinc   = integer(iter);                 // i: increment number
	ftnlen      cmname_len = strlen(fVUMAT_name);        // f2c: length of cmname string
	// below were added by Harold for VUMAT
	integer     lanneal = 0;                            // i: whether this analysis describes an annealing process
	integer     nfieldv = 0;                            // i: number of user defined field varibles - default to 0
	integer     nblock = 1;                             // i: number of material points to be processed - usually 1 IP
	doublereal  enerInelasOld = 0.0;                    
	doublereal  enerInelasNew = 0.0;                    // i: dissipated internal energy per unit mass - set to 0
	doublereal  enerInternOld = 0.0;
	doublereal  enerInternNew = 0.0;                    // i: these are set to 0 because BCJ does not define these
	doublereal  tempOld = 0.0;
	doublereal  tempNew = 0.0;                          // i: these are set to 0 because BCJ VUMAT uses the state
	                                                    //    variables arrays (SV) to track the temperature evolution
	doublereal* stretchold = fUOld2.Pointer();           // i: Stretch tensor at beginning of increment
	doublereal* stretchnew = fUNew2.Pointer();           // i: Stretch tensor at end of increment
	doublereal* relspininc = fRelSpin.Pointer();        // i: Relative spin increment
	doublereal  density = fDensity;                   // i: Density of material
	doublereal* stressnew = fstress.Pointer();          // o: This is the stress to be updated
	doublereal* statevnew = fstatv.Pointer();           // o: This is the state variable array to be updated

//DEBUG
#if VUMAT_DEBUG
int d_width = OutputWidth(flog, fstress.Pointer());
flog << " THE INPUT\n";
flog << setw(10) << "time:" << setw(d_width) << time[0]  << '\n';
flog << setw(10) << " stress: " << fstress.no_wrap() << '\n';
flog << setw(10) << " strain: " << fstrain.no_wrap() << '\n';
flog << setw(10) << "dstrain: " << fdstran.no_wrap() << '\n';
flog << setw(10) << "  state:\n";
flog << fstatv.wrap(5) << '\n';
#endif
//DEBUG

	/* call VUMAT wrapper */
       VUMAT(&nblock, &ndi, &nshr, &nstatv, &nfieldv, &nprops, &lanneal, &stime, &totime, &dtime, cmname, coords,
       &celent, props, &density, dstran, relspininc, &tempOld, stretchold, dfgrd0, predef, stressold, statevold,
       &enerInternOld, &enerInelasOld, &tempNew, stretchnew, dfgrd1, dpred, stressnew, statevnew, 
       &enerInternNew, &enerInelasNew);
 
//DEBUG
#if VUMAT_DEBUG
flog << " THE OUTPUT\n";
flog << setw(10) << " stress: " << fstress.no_wrap() << '\n';
flog << setw(10) << " state:\n" << '\n';
flog << fstatv.wrap(5) << endl;
#endif
//DEBUG

	/* update strain */
	fstrain += fdstran;
	
	/* write to storage */
	Store(CurrentElement(), CurrIP());
}

/* set variables to last converged */
void ABAQUS_VUMAT_BaseT::Reset_VUMAT_Increment(void)
{
  ///* assign "last" to "current" */
  //fstress    = fstress_last;
  //fstrain    = fstrain_last;
	//fsse_pd_cd = fsse_pd_cd_last;
  //fstatv     = fstatv_last;
  /* instead, assign "current" to "last" */
  fstress_last = fstress;
  fstrain_last = fstrain;
  fstatv_last = fstatv;
}

/* set stress/strain arguments */
void ABAQUS_VUMAT_BaseT::Set_VUMAT_Arguments(void)
{
	/* integration point coordinates */
	ContinuumElement().IP_Coords(fIPCoordinates);	
	fcoords[0] = doublereal(fIPCoordinates[0]);
	fcoords[1] = doublereal(fIPCoordinates[1]);
	if (NumSD() == 3)
		fcoords[2] = doublereal(fIPCoordinates[2]);

	/* deformation gradient at beginning of increment */
	fA_nsd = F_total_last();
	dMatrixT_to_ABAQUS(fA_nsd, fdfgrd0);

	/* stretch at beginning of increment */
	bool perturb_repeated_roots = false;
	fDecomp->PolarDecomp(fA_nsd, fROld, fUOld, perturb_repeated_roots);
	fUOld2 = fUOld;

	/* deformation gradient at end of increment */
	const dMatrixT& F_n = F();
	dMatrixT_to_ABAQUS(F_n, fdfgrd1);

	/* stretch at end of increment */
	fDecomp->PolarDecomp(F_n, fRNew, fUNew, perturb_repeated_roots);
	fUNew2 = fUNew;

	/* relative deformation gradient */
	fA_nsd.Inverse();
	fF_rel.MultAB(F_n, fA_nsd);

	/* polar decomposition - intermediate configuration */
	fDecomp->PolarDecomp(fF_rel, fA_nsd, fU1, perturb_repeated_roots);

	/* incremental rotation */
	dMatrixT_to_ABAQUS(fA_nsd, fdrot);

	/* Compute the relative spin here */
	// BLAH....

	/* incremental strain */
	fU2 = fU1;
	fU1.PlusIdentity(-1.0);
	fU2.PlusIdentity( 1.0);
	fU2.Inverse();
	fU1U2.MultAB(fU1, fU2);
	if (NumSD() == 2)
	{
		fdstran[0] = 2.0*doublereal(fU1U2[0]); // 11
		fdstran[1] = 2.0*doublereal(fU1U2[1]); // 22
		fdstran[3] = 2.0*doublereal(fU1U2[2]); // 12
	}
	else
	{
		fdstran[0] = 2.0*doublereal(fU1U2[0]); // 11
		fdstran[1] = 2.0*doublereal(fU1U2[1]); // 22
		fdstran[2] = 2.0*doublereal(fU1U2[2]); // 33
		fdstran[5] = 2.0*doublereal(fU1U2[3]); // 23
		fdstran[4] = 2.0*doublereal(fU1U2[4]); // 13
		fdstran[3] = 2.0*doublereal(fU1U2[5]); // 12
	}

	/* total integrated strain */
	ABAQUS_to_dSymMatrixT(fstrain.Pointer(), fU1);
	fU2.MultQBQT(fA_nsd, fU1);
	dSymMatrixT_to_ABAQUS(fU2, fstrain.Pointer(), true);

	/* rotate LAST stress to current configuration, instead of current stress */
	ABAQUS_to_dSymMatrixT(fstress_last.Pointer(), fU1);
	fU2.MultQBQT(fA_nsd, fU1);
	dSymMatrixT_to_ABAQUS(fU2, fstress_last.Pointer(), true);
}
#endif /* __F2C__ */
