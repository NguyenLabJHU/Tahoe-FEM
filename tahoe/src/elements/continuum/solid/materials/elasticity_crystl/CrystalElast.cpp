/* $Id: CrystalElast.cpp,v 1.10 2004-07-15 08:28:16 paklein Exp $ */
#include "CrystalElast.h"
#include "CrystalElastMat.h"
#include "CrystalElastLat.h"
#include "Utils.h"

#include "FSSolidMatT.h"
#include "StringT.h"

using namespace Tahoe;

/* number of elastic material properties used in computations */
const int kNumMatProp = 3;

/* initialization flag value */
const int kIsInit = 1;

/* spatial dimensions of the problem */
const int kNSD = 3;

CrystalElast::CrystalElast(ifstreamT& in, const FSMatSupportT& support) :
	ParameterInterfaceT("crystial_elasticity"),
//	FDHookeanMatT(in, support),
  // pointers to supporting classes
  	fCrystalElastLat (NULL),
  	fCrystalElastMat (NULL),
  // material properties
        fMatProp      (kNumMatProp)
{
ExceptionT::GeneralFail("CrystalElast::CrystalElast", "out of date");
#if 0
  // input file
  StringT filename;
  in >> filename;

  // generate relative path in native format
  filename.ToNativePathName();
  StringT path;
  path.FilePath(in.filename());
  filename.Prepend(path);

  OpenExternal(fInput, filename, "CrystalElast data");
  if (in.skip_comments())
    fInput.set_marker(in.comment_marker());

  // read number of crystals per integration point 
  fInput >> fNumGrain;

  // read crystal orientations
  SetLatticeOrientation();

  // set crystal elasticity type
  SetCrystalElasticity();
#endif
}

CrystalElast::~CrystalElast()
{
  delete fCrystalElastLat;
  delete fCrystalElastMat;
}

void CrystalElast::Initialize()
{
  // allocate space for all elements
  AllocateElements();

  // initialize state variables in all elements
  InitializeCrystalVariables();
}

  // print input values
#if 0
  out << " FDcrystal data:\n";
  out << "    Number of grains   . . . . . . . . . . . . . = " << fNumGrain    << "\n";
  out << "    Lattice Orientation code . . . . . . . . . . = " << fODFCode << "\n";
  out << "    Initial Temperature   .  . . . . . . . . . . = " << fInit_Temp_DegC << "\n";
#endif

/* set (material) tangent modulus */
void CrystalElast::SetModulus(dMatrixT& modulus)
{
#pragma unused(modulus)
}

void CrystalElast::AllocateElements()
{
  // determine storage size
  int i_size = NumIP();
  int d_size = NumVariablesPerElement();

  // allocate space for all elements
  for (int elem = 0; elem < NumElements(); elem++)
    {
      // get pointer to element elem
      ElementCardT& element = ElementCard(elem);
      
      // construct element
      element.Dimension(i_size, d_size);
      
      // initialize values
      element.IntegerData() = kIsInit;
      element.DoubleData()  = 0.0;
    }
}

void CrystalElast::SetLatticeOrientation()
{
  // input code to distribute ODF at each IP/ELEM
  fInput >> fODFCode;

  // read orientation data
  fCrystalElastLat = new CrystalElastLat(*this);
  if (!fCrystalElastLat) throwMemoryError("CrystalElast::SetLatticeOrientation");

  // number of elements and integration points
  int numelem = NumElements();
  int numint = NumIP();

  // allocate array for euler angles at integration point
  fangles.Dimension(fNumGrain);
  for (int i = 0; i < fNumGrain; i++)
    fangles[i].Dimension(3);

  // allocate array to hold crystal orientations
  fEuler.Dimension(numelem);
  for (int i = 0; i< numelem; i++)
    {
      fEuler[i].Dimension(numint, fNumGrain);
      for (int j = 0; j < numint; j++)
	for (int k = 0; k < fNumGrain; k++)
	  fEuler[i](j,k).Dimension(3);
    }

  // assign orientation angles to each IP/ELEM
  fCrystalElastLat->AssignEulerAngles(fODFCode, numelem, numint,
				    fNumGrain, fEuler);
}

void CrystalElast::SetCrystalElasticity()
{
  // input initial temperature
  fInput >> fInit_Temp_DegC;

  // select crystal elasticity type
  fCrystalElastMat = new CrystalElastMat(*this);
  if (!fCrystalElastMat) throwMemoryError("CrystalElast::SetCrystalElasticity");
}
