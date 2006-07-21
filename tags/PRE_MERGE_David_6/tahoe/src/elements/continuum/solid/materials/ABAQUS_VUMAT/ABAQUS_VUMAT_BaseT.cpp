/* $Id: ABAQUS_VUMAT_BaseT.cpp,v 1.27 2005-08-04 07:06:02 paklein Exp $ */
#include "ABAQUS_VUMAT_BaseT.h"

#ifdef __F2C__

#include <ctype.h>
#include <float.h>

#include "ContinuumElementT.h"
#include "SpectralDecompT.h"
#include "ThermalDilatationT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"

using namespace Tahoe;

/* constructor */
ABAQUS_VUMAT_BaseT::ABAQUS_VUMAT_BaseT(void):
	ParameterInterfaceT("ABAQUS_VUMAT_material"),
	fDebug(false),
	fTangentType(GlobalT::kSymmetric),
	fPressure(0.0),
	fDecomp(NULL)
{

}

/* destructor */
ABAQUS_VUMAT_BaseT::~ABAQUS_VUMAT_BaseT(void) {
	delete fDecomp;
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
		fstress_last = fstress;
		fstrain_last = fstrain;
		fstatv_last  = fstatv;

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
		fstress = fstress_last;
		fstrain = fstrain_last;
		fstatv  = fstatv_last;

		/* write to storage */
		Store(element, ip);
	}
}

const dSymMatrixT& ABAQUS_VUMAT_BaseT::s_ij(void)
{
	/* call VUMAT */
	if (MaterialSupport().RunState() == GlobalT::kFormRHS ||
	    MaterialSupport().RunState() == GlobalT::kFormLHS)
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

/* returns the strain energy density for the specified strain */
double ABAQUS_VUMAT_BaseT::StrainEnergyDensity(void) {
	return 0.0;
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int ABAQUS_VUMAT_BaseT::NumOutputVariables(void) const {
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
		ExceptionT::SizeMismatch("ABAQUS_VUMAT_BaseT::ComputeOutput",
			"output array should be length %d not %d", fOutputIndex.Length(), output.Length());

	/* load stored data */
	Load(CurrentElement(), CurrIP());

	/* collect variables */
	for (int i = 0; i < fOutputIndex.Length(); i++)
		output[i] = double(fstatv[fOutputIndex[i]]);
}

/* describe the parameters needed by the interface */
void ABAQUS_VUMAT_BaseT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FSIsotropicMatT::DefineParameters(list);

	ParameterT debug(fDebug, "debug");
	debug.SetDefault(fDebug);
	list.AddParameter(debug);

	/* file with UMAT materials parameters */
	list.AddParameter(ParameterT::Word, "VUMAT_parameter_file");
}

/* accept parameter list */
void ABAQUS_VUMAT_BaseT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "ABAQUS_VUMAT_BaseT::TakeParameterList";	

	/* inherited */
	FSIsotropicMatT::TakeParameterList(list);

	fDebug = list.GetParameter("debug");

	/* dimension work space */
	int nsd = NumSD();
	fIPCoordinates.Dimension(nsd);
	fF_rel.Dimension(nsd);
	fROld.Dimension(nsd);
	fRNew.Dimension(nsd);
	fA_nsd.Dimension(nsd);
	fU1.Dimension(nsd);
	fU2.Dimension(nsd);
	fU1U2.Dimension(nsd);
	fUOld.Dimension(nsd);
	fUNew.Dimension(nsd);

	/* open VUMAT parameters file */
	StringT path;
	path.FilePath(MaterialSupport().InputFile());
	StringT params = list.GetParameter("VUMAT_parameter_file");
	params.ToNativePathName();
	params.Prepend(path);
	ifstreamT in('#', params);
	if (!in.is_open())
		ExceptionT::GeneralFail(caller, "could not open file \"%s\"",
			params.Pointer());

	/* read ABAQUS-format input */
	nstatv = 0;
	bool nonsym = false;
	Read_ABAQUS_Input(in, fVUMAT_name, fProperties, fDensity, nstatv, nonsym);
		
	/* VUMAT dimensions */
	ndi = 3; // always 3 direct components
	if (nsd == 2)
		nshr = 1;
	else if (nsd == 3)
		nshr = 3;
	else
		ExceptionT::GeneralFail(caller, "unexpected dimension %d", nsd);
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

	/* may not need that (above) */
	fBlockSize += nstatv;      // fstatv
	fBlockSize += ntens;       // fstress_last
	fBlockSize += ntens;       // fstrain_last

	/* may not need this (above) */
	fBlockSize += nstatv;      // fstatv_last
	
	/* argument array */
	fArgsArray.Dimension(fBlockSize);

	/* assign pointers */
	doublereal* parg = fArgsArray.Pointer();
	fstress.Set(ntens, parg);        parg += ntens;
	fstrain.Set(ntens, parg);        parg += ntens;
	fstatv.Set(nstatv, parg);        parg += nstatv;
	fstress_last.Set(ntens, parg);   parg += ntens;
	fstrain_last.Set(ntens, parg);   parg += ntens;
	fstatv_last.Set(nstatv, parg);
	
	/* VUMAT array arguments */
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
	fDecomp = new SpectralDecompT(nsd);
	if (!fDecomp) ExceptionT::OutOfMemory(caller);

	/* write properties array */
	ofstreamT& out = MaterialSupport().Output();
	out << " Number of ABAQUS VUMAT internal variables. . . . = " << nstatv << '\n';
	out << " Number of ABAQUS VUMAT properties. . . . . . . . = " << fProperties.Length() << '\n';
	PrintProperties(out);

	/* notify */
	if (fThermal->IsActive())
		cout << "\n ABAQUS_VUMAT_BaseT::Initialize: thermal strains must\n"
		     <<   "    be handled within the VUMAT\n" << endl;

	/* set material output variables/labels */
	SetOutputVariables(fOutputIndex, fOutputLabels);
}

/***********************************************************************
* Protected
***********************************************************************/

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
	/* assign "last" to "current" */
	fstress = fstress_last;
	fstrain = fstrain_last;
	fstatv  = fstatv_last;
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
