/* $Id: ABAQUS_UMAT_BaseT.cpp,v 1.1.1.1 2001-01-29 08:20:26 paklein Exp $ */
/* created: paklein (05/14/2000)                                          */

#include "ABAQUS_UMAT_BaseT.h"

#ifdef __F2C__

#include <ctype.h>
#include <float.h>

#include "fstreamT.h"
#include "FEManagerT.h"
#include "ElasticT.h"
#include "ShapeFunctionT.h"
#include "SpectralDecompT.h"
#include "ThermalDilatationT.h"

/* constructor */
ABAQUS_UMAT_BaseT::	ABAQUS_UMAT_BaseT(ifstreamT& in, const ElasticT& element):
	FDStructMatT(in, element),
	fRunState(ContinuumElement().RunState()),
	fTangentType(GlobalT::kSymmetric),
	fLocLastDisp(element.LastDisplacements()),
	fModulus(dSymMatrixT::NumValues(NumSD())),
	fStress(NumSD()),
	fIPCoordinates(NumSD()),
	fDecomp(NULL),
	fF_rel(NumSD()),
	fA_nsd(NumSD()),
	fU1(NumSD()), fU2(NumSD()), fU1U2(NumSD())
{
	/* read ABAQUS-format input */
	nstatv = 0;
	Read_ABAQUS_Input(in);
	
	/* UMAT dimensions */
	ndi = 3; // always 3 direct components
	int nsd = NumSD();
	if (nsd == 2)
		nshr = 1;
	else if (nsd == 3)
		nshr = 3;
	else
		throw eGeneralFail;
	ntens = ndi + nshr;

	/* modulus storage */
	if (fTangentType == GlobalT::kDiagonal)
		fModulusDim = ntens;
	else if (fTangentType == GlobalT::kSymmetric)
	{
		if (nsd == 2) fModulusDim = 10;
		else if (nsd == 3) fModulusDim = 21;
		else throw eGeneralFail;
	}
	else if (fTangentType == GlobalT::kNonSymmetric)
		fModulusDim = ntens*ntens;
	else
		throw eGeneralFail;

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
	fArgsArray.Allocate(fBlockSize);

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
	fddsdde.Allocate(ntens);
	fdstran.Allocate(ntens);
	fdstran = 0.0;
	fdrot.Allocate(3);   // always 3
	fdrot.Identity();
	fdfgrd0.Allocate(3); // always 3
	fdfgrd0.Identity();
	fdfgrd1.Allocate(3); // always 3
	fdfgrd1.Identity();
	fcoords.Allocate(nsd);
	
	/* spectral decomp */
	fDecomp = new SpectralDecompT(NumSD());
	if (!fDecomp) throw eOutOfMemory;

//DEBUG
#if 1
flog.open("UMAT.log");
flog.precision(DBL_DIG);
flog.setf(ios::showpoint);
flog.setf(ios::right, ios::adjustfield);
flog.setf(ios::scientific, ios::floatfield);
#endif
}

/* destructor */
ABAQUS_UMAT_BaseT::~ABAQUS_UMAT_BaseT(void)
{
	delete fDecomp;
	fDecomp = NULL;
}

/* print parameters */
void ABAQUS_UMAT_BaseT::Print(ostream& out) const
{
	/* inherited */
	FDStructMatT::Print(out);
	
	/* write properties array */
	out << " Number of ABAQUS UMAT internal variables. . . . = " << nstatv << '\n';
	out << " Number of ABAQUS UMAT properties. . . . . . . . = " << fProperties.Length() << '\n';
	PrintProperties(out);
}

/* disable multiplicative thermal strains */
void ABAQUS_UMAT_BaseT::Initialize(void)
{
	/* inherited */
	FDStructMatT::Initialize();

	/* notify */
	if (fThermal->IsActive())
		cout << "\n ABAQUS_UMAT_BaseT::Initialize: thermal strains must\n"
		     <<   "    be handled within the UMAT\n" << endl;
	
	/* disable thermal transform */
	SetFmodMult(NULL);	
}


/* materials initialization */
bool ABAQUS_UMAT_BaseT::NeedsPointInitialization(void) const { return true; }
void ABAQUS_UMAT_BaseT::PointInitialize(void)
{
	/* allocate element storage */
	if (CurrIP() == 0)
	{
		ElementCardT& element = CurrentElement();
		element.Allocate(0, fBlockSize*NumIP());
	
		/* initialize */
		element.DoubleData() = 0.0;
	}

	/* call UMAT - time signals initialization */
	Call_UMAT(0.0, 0.0, 0, 0);

	/* store results as last converged */
	if (CurrIP() == NumIP() - 1) UpdateHistory();
}

/* update/reset internal variables */
void ABAQUS_UMAT_BaseT::UpdateHistory(void)
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

void ABAQUS_UMAT_BaseT::ResetHistory(void)
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

/* required parameter flags */
bool ABAQUS_UMAT_BaseT::NeedDisp(void) const     { return true; }
bool ABAQUS_UMAT_BaseT::NeedLastDisp(void) const { return true; }
bool ABAQUS_UMAT_BaseT::NeedVel(void) const      { return true; }

/* spatial description */
const dMatrixT& ABAQUS_UMAT_BaseT::c_ijkl(void)
{
	/* load stored data */
	Load(CurrentElement(), CurrIP());

	int nsd = NumSD();
	if (nsd != 2 && nsd != 3) throw eGeneralFail;
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
		cout << "\n ABAQUS_UMAT_BaseT::c_ijkl: index mapping for nonsymmetric\n"
		     <<   "     tangent is not implemented" << endl;
		throw eGeneralFail;
	}
	else throw eGeneralFail;	

	return fModulus;
}

const dSymMatrixT& ABAQUS_UMAT_BaseT::s_ij(void)
{
	/* call UMAT */
	if (fRunState == GlobalT::kFormRHS)
	{
		const FEManagerT& fe_man = ContinuumElement().FEManager();
		double  t = fe_man.Time();
		double dt = fe_man.TimeStep();
		int  step = fe_man.StepNumber();
		int  iter = fe_man.IterationNumber();
		Call_UMAT(t, dt, step, iter);
	}
	else
		/* load stored data */
		Load(CurrentElement(), CurrIP());

	/* copy/convert stress */
	ABAQUS_to_dSymMatrixT(fstress.Pointer(), fStress);
	return fStress;
}

/* material description */
const dMatrixT& ABAQUS_UMAT_BaseT::C_IJKL(void)
{
	/* set deformation gradient */
	F();

	/* transform */
	return c_to_C(ABAQUS_UMAT_BaseT::c_ijkl());
}

const dSymMatrixT& ABAQUS_UMAT_BaseT::S_IJ(void)
{
//NOTE: F(), the current deformation gradient must be the last
//      one computed in s_ij()

	/* transform */
	return  s_to_S(ABAQUS_UMAT_BaseT::s_ij());
}

/* returns the strain energy density for the specified strain */
double ABAQUS_UMAT_BaseT::StrainEnergyDensity(void)
{
	/* load stored data */
	Load(CurrentElement(), CurrIP());

	/* pull from storage */
	return double(fsse_pd_cd[0]);
}

/* returns the number of variables computed for nodal extrapolation
* during for element output, ie. internal variables. Returns 0
* by default */
int ABAQUS_UMAT_BaseT::NumOutputVariables(void) const
{
	/* set material output variables/labels */
	if (fOutputIndex.Length() == 0)
	{
		//TEMP - better place for this?
		ABAQUS_UMAT_BaseT* tmp = (ABAQUS_UMAT_BaseT*) this;
		tmp->SetOutputVariables(tmp->fOutputIndex, tmp->fOutputLabels);
	}

	return fOutputIndex.Length();
}

void ABAQUS_UMAT_BaseT::OutputLabels(ArrayT<StringT>& labels) const
{
	labels.Allocate(fOutputLabels.Length());
	for (int i = 0; i < labels.Length(); i++)
		labels[i] = fOutputLabels[i];
}

void ABAQUS_UMAT_BaseT::ComputeOutput(dArrayT& output)
{
	/* check */
	if (output.Length() != fOutputIndex.Length())
	{
		cout << "\n ABAQUS_UMAT_BaseT::ComputeOutput: not enough space to return\n"
		     <<   "     output variables: given " << output.Length()
		     << ". expecting " << fOutputIndex.Length() << "." << endl;
		throw eSizeMismatch;
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
void ABAQUS_UMAT_BaseT::PrintName(ostream& out) const
{
	/* inherited */
	FDStructMatT::PrintName(out);
	out << "    ABAQUS user material: " << fUMAT_name << '\n';
}

void ABAQUS_UMAT_BaseT::PrintProperties(ostream& out) const
{
	/* just write numbered list */
	int d_width = OutputWidth(out, fProperties.Pointer());
	for (int i = 0; i < fProperties.Length(); i++)	
		out << setw(kIntWidth) << i+1
		    << setw(  d_width) << fProperties[i] << '\n';
}

/* conversion functions */
void ABAQUS_UMAT_BaseT::dMatrixT_to_ABAQUS(const dMatrixT& A,
	nMatrixT<doublereal>& B) const
{
#if __option(extended_errorcheck)
	/* expecting ABAQUS matrix to be 3D always */
	if (B.Rows() != 3 ||
	    B.Cols() != 3) throw eGeneralFail;
#endif

	if (NumSD() == 2)
	{
		B(0,0) = doublereal(A(0,0));
		B(1,0) = doublereal(A(1,0));
		B(0,1) = doublereal(A(0,1));
		B(1,1) = doublereal(A(1,1));
	}
	else
	{
		doublereal* pB = B.Pointer();
		double* pA = A.Pointer();
		*pB++ = doublereal(*pA++);
		*pB++ = doublereal(*pA++);
		*pB++ = doublereal(*pA++);
		*pB++ = doublereal(*pA++);
		*pB++ = doublereal(*pA++);
		*pB++ = doublereal(*pA++);
		*pB++ = doublereal(*pA++);
		*pB++ = doublereal(*pA++);
		*pB   = doublereal(*pA);
	}	
}

void ABAQUS_UMAT_BaseT::ABAQUS_to_dSymMatrixT(const doublereal* pA,
	dSymMatrixT& B) const
{
	double* pB = B.Pointer();
	if (NumSD() == 2)
	{
		*pB++ = double(pA[0]); // 11
		*pB++ = double(pA[1]); // 22
		*pB   = double(pA[3]); // 12
	}
	else
	{
		*pB++ = double(pA[0]); // 11
		*pB++ = double(pA[1]); // 22
		*pB++ = double(pA[2]); // 33
		*pB++ = double(pA[5]); // 23
		*pB++ = double(pA[4]); // 13
		*pB   = double(pA[3]); // 12
	}
}

void ABAQUS_UMAT_BaseT::dSymMatrixT_to_ABAQUS(const dSymMatrixT& A,
	doublereal* pB) const
{
	double* pA = A.Pointer();
	if (NumSD() == 2)
	{
		*pB++ = doublereal(pA[0]); // 11
		*pB++ = doublereal(pA[1]); // 22
		*pB++;                     // 33 - leave unchanged
		*pB   = doublereal(pA[2]); // 12
	}
	else
	{
		*pB++ = doublereal(pA[0]); // 11
		*pB++ = doublereal(pA[1]); // 22
		*pB++ = doublereal(pA[2]); // 33
		*pB++ = doublereal(pA[5]); // 12
		*pB++ = doublereal(pA[4]); // 13
		*pB   = doublereal(pA[3]); // 23
	}
}

/***********************************************************************
* Private
***********************************************************************/

/* read ABAQUS format input - reads until first line
* 	not containing a comment of keyword
*
* looks for kewords:
*    *MATERIAL
*    *USER MATERIAL
*    *DEPVAR
*/
void ABAQUS_UMAT_BaseT::Read_ABAQUS_Input(ifstreamT& in)
{
	/* disable comment skipping */
	int skip = in.skip_comments();
	char marker;
	if (skip)
	{
		marker = in.comment_marker();
		in.clear_marker();
	}

	/* required items */
	bool found_MATERIAL = false;
	bool found_USER     = false;

// when to stop
//    need at least material name
//    allow for no depvar and no constants
//    allow *DENSITY, *EXPANSION, *DAMPING <- mass damping

	/* advance to first keyword */
	bool found_keyword = Next_ABAQUS_Keyword(in);
	while (found_keyword)
	{
		StringT next_word;
		Read_ABAQUS_Word(in, next_word);
	
		/* first keyword */
		if (!found_MATERIAL)
		{
			if (next_word == "MATERIAL")
			{
				found_MATERIAL = true;
				/* scan for parameters */
				while (Skip_ABAQUS_Symbol(in, ','))
				{
					StringT next_parameter;
					Read_ABAQUS_Word(in, next_parameter);
					if (next_parameter == "NAME")
					{
						/* skip '=' */
						if (!Skip_ABAQUS_Symbol(in, '=')) throw eBadInputValue;
						
						/* read name */
						Read_ABAQUS_Word(in, fUMAT_name, false);
					}
					else
						cout << "\n ABAQUS_UMAT_BaseT::Read_ABAQUS_Input: skipping parameter:"
						     << next_parameter << endl;
				}
			}
			else
			{
				cout << "\n ABAQUS_UMAT_BaseT::Read_ABAQUS_Input: first keyword must be\n"
				     <<   "     *MATERIAL not *" << next_word << endl;
				throw eBadInputValue;
			}
		}
		/* other keywords */
		else if (next_word == "USER")
		{
			Read_ABAQUS_Word(in, next_word);
			if (next_word == "MATERIAL")
			{
				found_USER = true;
				/* scan for parameters */
				while (Skip_ABAQUS_Symbol(in, ','))
				{
					StringT next_parameter;
					Read_ABAQUS_Word(in, next_parameter);
					if (next_parameter == "CONSTANTS")
					{
						/* skip '=' */
						if (!Skip_ABAQUS_Symbol(in, '=')) throw eBadInputValue;

						int nprops = -1;
						in >> nprops;
						if (nprops < 0)
						{
							cout << "\n ABAQUS_UMAT_BaseT::Read_ABAQUS_Input: error reading "
							     << next_parameter << ": " << nprops << endl;
							throw eBadInputValue;
						}
						
						/* read properties */
						fProperties.Allocate(nprops);
						in.clear_line();
						Skip_ABAQUS_Comments(in);
						for (int i = 0; i < nprops && in.good(); i++)
						{	
							in >> fProperties[i];
							Skip_ABAQUS_Symbol(in, ',');
						}
						
						/* stream error */
						if (!in.good())
						{
							cout << "\n ABAQUS_UMAT_BaseT::Read_ABAQUS_Input: error reading "
							     << next_parameter << endl;
							throw eBadInputValue;
						}
					}
					else
						cout << "\n ABAQUS_UMAT_BaseT::Read_ABAQUS_Input: skipping parameter:"
						     << next_parameter << endl;
				}
			
			}			
			else
			{
				cout << "\n ABAQUS_UMAT_BaseT::Read_ABAQUS_Input: expecting MATERIAL after\n"
				     <<   "     keyword *USER not " << next_word << endl;
				throw eBadInputValue;
			}
		
		
		}
		else if (next_word == "DEPVAR")
		{
			nstatv = -1;
			in >> nstatv;
			if (nstatv < 0)
			{
				cout << "\n ABAQUS_UMAT_BaseT::Read_ABAQUS_Input: error reading "
				     << next_word << endl;
				throw eBadInputValue;
			}

			/* clear trailing comma */
			Skip_ABAQUS_Symbol(in, ',');
		}
		/* skipping all others */
		else
			cout << "\n ABAQUS_UMAT_BaseT::Read_ABAQUS_Input: skipping keyword *"
			     << next_word << endl;

		/* next keyword */
		found_keyword = Next_ABAQUS_Keyword(in);
	}

	/* check all data found */
	if (!found_USER || !found_MATERIAL)
	{
		cout << '\n';
		if (!found_MATERIAL)
			cout << " ABAQUS_UMAT_BaseT::Read_ABAQUS_Input: missing keyword: *MATERIAL\n";
		if (!found_USER)
			cout << " ABAQUS_UMAT_BaseT::Read_ABAQUS_Input: missing keyword: *USER MATERIAL\n";
		cout.flush();
		throw eBadInputValue;
	}

	/* restore comment skipping */
	if (skip) in.set_marker(marker);
}

/* advance to next non-comment line */
bool ABAQUS_UMAT_BaseT::Next_ABAQUS_Keyword(ifstreamT& in) const
{
//rules: (1) keyword lines have a '*' in column 1
//       (2) comment lines have a '*' in column 1 and 2

	while (in.good())
	{
		/* advance to next '*' or alpanumeric character */
		char c = in.next_char();
		while (c != '*' && !isgraph(c) && in.good())
		{
			in.get(c);	
			c = in.next_char();
		}

		/* stream error */
		if (!in.good()) break;

		/* comment or keyword? */
		if (c == '*')
		{
			/* next */
			in.get(c);
			c = in.next_char();
			
			/* comment */
			if (c == '*')
				in.clear_line();
			/* keyword */
			else
				return true;
		}
		/* other alphanumeric character */
		else
		{
			in.putback(c);
			return false;
		}
	}
	
	/* stream error */
	in.clear();
	return false;
}

void ABAQUS_UMAT_BaseT::Skip_ABAQUS_Comments(ifstreamT& in)
{
	while (in.good())
	{
		/* advance to next '*' or alpanumeric character */
		char c = in.next_char();
		while (c != '*' && !isgraph(c) && in.good())
		{
			in.get(c);	
			c = in.next_char();
		}

		/* stream error */
		if (!in.good()) break;
	
		/* comment or keyword? */
		if (c == '*')
		{
			in.get(c);
			c = in.next_char();
			
			/* comment */
			if (c == '*')
				in.clear_line();
			/* keyword */
			else
				return;
		}
		/* other alphanumeric character */
		else
			return;
	}

	/* stream error */
	in.clear();
}

bool ABAQUS_UMAT_BaseT::Skip_ABAQUS_Symbol(ifstreamT& in, char c) const
{
	/* advance stream */
	char a = in.next_char();

	if (!in.good())
	{
		in.clear();
		return false;
	}
	else if (a == c)
	{
		in.get(a);
		return true;
	}
	else
		return false;
}

void ABAQUS_UMAT_BaseT::Read_ABAQUS_Word(ifstreamT& in, StringT& word, bool to_upper) const
{
	const int max_word_length = 255;
	char tmp[max_word_length + 1];
	
	int count = 0;
	char c;
	/* skip whitespace to next word */
	c = in.next_char();
	
	in.get(c);
	while (count < max_word_length && isalnum(c))
	{
		tmp[count++] = c;
		in.get(c);
	}
	in.putback(c);
	tmp[count] = '\0';

	/* copy in */
	word = tmp;
	if (to_upper) word.ToUpper();
}

/* load element data for the specified integration point */
void ABAQUS_UMAT_BaseT::Load(ElementCardT& element, int ip)
{
	/* fetch internal variable array */
	dArrayT& d_array = element.DoubleData();

	/* copy/convert */
	double* pd = d_array.Pointer(fBlockSize*ip);
	doublereal* pdr = fArgsArray.Pointer();
	for (int i = 0; i < fBlockSize; i++)
		*pdr++ = doublereal(*pd++);
}

void ABAQUS_UMAT_BaseT::Store(ElementCardT& element, int ip)
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
void ABAQUS_UMAT_BaseT::Call_UMAT(double t, double dt, int step, int iter)
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
	doublereal  pnewdt;                                 // o: suggested time step (automatic time integration)
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

//DEBUG
int d_width = OutputWidth(flog, fstress.Pointer());
if (CurrIP() == 0 && (step == 1 || step == 26))
{
flog << " THE INPUT\n";
flog << setw(10) << "time:" << setw(d_width) << time[0]  << '\n';
flog << setw(10) << " stress: " << fstress.no_wrap() << '\n';
flog << setw(10) << " strain: " << fstrain.no_wrap() << '\n';
flog << setw(10) << "dstrain: " << fdstran.no_wrap() << '\n';
flog << setw(10) << "  state:\n";
flog << fstatv.wrap(5) << '\n';
}
//DEBUG

	/* call UMAT wrapper */
	UMAT(stress, statev, ddsdde, &sse, &spd, &scd, &rpl, ddsddt, drplde,
		&drpldt, stran, dstran, time, &dtime, &temp, &dtemp, predef, dpred, cmname,
		&ndi, &nshr, &ntens, &nstatv, props, &nprops, coords, drot, &pnewdt, &celent,
		dfgrd0, dfgrd1, &noel, &npt, &layer, &kspt, &kstep, &kinc, cmname_len);

//DEBUG
if (false && CurrIP() == 0 && (step == 1 || step == 26))
{
flog << " THE OUTPUT\n";
flog << setw(10) << " stress: " << fstress.no_wrap() << '\n';
flog << setw(10) << " state:\n" << '\n';
flog << fstatv.wrap(5) << endl;
}
//DEBUG

	/* update strain */
	fstrain += fdstran;

	/* store modulus */
	Store_UMAT_Modulus();
	
	/* write to storage */
	Store(CurrentElement(), CurrIP());
}

/* set variables to last converged */
void ABAQUS_UMAT_BaseT::Reset_UMAT_Increment(void)
{
	/* assign "last" to "current" */
	fstress    = fstress_last;
	fstrain    = fstrain_last;
	fsse_pd_cd = fsse_pd_cd_last;
	fstatv     = fstatv_last;
}

/* set stress/strain arguments */
void ABAQUS_UMAT_BaseT::Set_UMAT_Arguments(void)
{
	/* integration point coordinates */
	const ShapeFunctionT& shape_function = ShapeFunction();
	shape_function.IPCoords(fIPCoordinates);	
	fcoords[0] = doublereal(fIPCoordinates[0]);
	fcoords[1] = doublereal(fIPCoordinates[1]);
	if (NumSD() == 3)
		fcoords[2] = doublereal(fIPCoordinates[2]);

	/* deformation gradient at beginning of increment */
	fA_nsd = F(fLocLastDisp);
	dMatrixT_to_ABAQUS(fA_nsd, fdfgrd0);
	
	/* deformation gradient at end of increment */
	const dMatrixT& F_n = F();
	dMatrixT_to_ABAQUS(F_n, fdfgrd1);

	/* relative deformation gradient */
	fA_nsd.Inverse();
	fF_rel.MultAB(F_n, fA_nsd);

	/* polar decomposition */
	fDecomp->PolarDecomp(fF_rel, fA_nsd, fU1, true);

	/* incremental rotation */
	dMatrixT_to_ABAQUS(fA_nsd, fdrot);
	
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
	dSymMatrixT_to_ABAQUS(fU2, fstrain.Pointer());

	/* rotate stress to current configuration */
	ABAQUS_to_dSymMatrixT(fstress.Pointer(), fU1);
	fU2.MultQBQT(fA_nsd, fU1);
	dSymMatrixT_to_ABAQUS(fU2, fstress.Pointer());
}

/* store the modulus */
void ABAQUS_UMAT_BaseT::Store_UMAT_Modulus(void)
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
	else throw eGeneralFail;	
}
#endif /* __F2C__ */
