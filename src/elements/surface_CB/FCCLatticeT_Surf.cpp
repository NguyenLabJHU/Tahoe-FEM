/* $Id: FCCLatticeT_Surf.cpp,v 1.1 2005-06-29 22:38:16 hspark Exp $ */
#include "FCCLatticeT_Surf.h"
#include "ParameterContainerT.h"

using namespace Tahoe;

/* number of atoms per shell */
/* Construct shells corresponding to atoms exactly on free surface */
/* Are actually counting ALL bonds; need to split energy in half for surface clusters */
const int atoms_per_shell[] = {8, 5, 12};
const int atoms_in_shells[] = {8, 13, 25};
static int AtomsInShells(int nshells) {
	if (nshells < 0 || nshells > 3) ExceptionT::OutOfRange();
	return atoms_in_shells[nshells-1];
};
const double sqrt2 = sqrt(2.0);
const double sqrt3 = sqrt(3.0);

/* constructor */
FCCLatticeT_Surf::FCCLatticeT_Surf(int nshells):
	ParameterInterfaceT("CB_lattice_FCC"),
	fNumShells(nshells)
{

}

/* information about subordinate parameter lists */
void FCCLatticeT_Surf::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	sub_list.AddSub("FCC_lattice_orientation", ParameterListT::Once, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FCCLatticeT_Surf::NewSub(const StringT& name) const
{
	if (name == "FCC_lattice_orientation")
	{
		ParameterContainerT* orientation = new ParameterContainerT(name);
		orientation->SetListOrder(ParameterListT::Choice);
	
		ParameterContainerT natural("FCC_natural");
		orientation->AddSub(natural);
		
		ParameterContainerT FCC110("FCC_110");
		ParameterT FCC110_type(ParameterT::Enumeration, "sense");
		FCC110_type.AddEnumeration("[1 0 0][0 1 1][0 -1 1]", 0);
		FCC110_type.AddEnumeration("[1 1 0][-1 1 0][0 0 1]", 1);
		FCC110_type.SetDefault(0);
		FCC110.AddParameter(FCC110_type);
		orientation->AddSub(FCC110);

		ParameterContainerT FCC111("FCC_111");
		ParameterT FCC111_type(ParameterT::Enumeration, "sense");
		FCC111_type.AddEnumeration("[-1 1 0][-1-1 2][ 1 1 1]", 0);
		FCC111_type.AddEnumeration("[ 1-1 0][ 1 1-2][ 1 1 1]", 1);
		FCC111_type.SetDefault(0);
		FCC111.AddParameter(FCC111_type);
		orientation->AddSub(FCC111);

		ParameterContainerT Euler_angles("FCC_Euler_angles");
		Euler_angles.AddParameter(ParameterT::Double, "theta");
		Euler_angles.AddParameter(ParameterT::Double, "phi");
		Euler_angles.AddParameter(ParameterT::Double, "psi");
		orientation->AddSub(Euler_angles);
	
		return orientation;
	}
	else /* inherited */
		return ParameterInterfaceT::NewSub(name);
}

/* accept parameter list */
void FCCLatticeT_Surf::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "FCCLatticeT_Surf::TakeParameterList";

	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);
	
	/* set Q */
	const ParameterListT& orientation = list.GetListChoice(*this, "FCC_lattice_orientation");
	dMatrixT Q;
	SetQ(orientation, Q);
	
	/* initialize bond table */
	Initialize(&Q);
}

/* set the transformation matrix for the given orientation */
void FCCLatticeT_Surf::SetQ(const ParameterListT& list, dMatrixT& Q)
{
	const char caller[] = "FCCLatticeT_Surf::SetQ";

	/* dimension */
	Q.Dimension(3);
	Q = 0.0;

	/* extract orientation */
	if (list.Name() == "FCC_natural")
		Q.Identity();
	else if (list.Name() == "FCC_110")
	{
		double cos45 = 0.5*sqrt2;

		int sense = list.GetParameter("sense");
		if (sense == 0) /* [1 0 0][0 1 1][0 -1 1] */ {
			Q(0,0) = 1.0;
			Q(1,1) = Q(2,2) = cos45;
			Q(1,2) =-cos45;
			Q(2,1) = cos45;	
		}
		else if (sense == 1) /* [1 1 0][-1 1 0][0 0 1] */ {
			Q(0,0) = Q(1,1) = cos45;
			Q(0,1) =-cos45;
			Q(1,0) = cos45;
			Q(2,2) = 1.0;
		}
		else
			ExceptionT::GeneralFail(caller, "unrecognized 110 sense %d", sense);
	}
	else if (list.Name() == "FCC_111")
	{
		int sense = list.GetParameter("sense");
		OrientationCodeT code = (sense == 1) ? kFCC3D111_b : kFCC3D111_a;

		double rt2b2 = sqrt2/2.0;
		double rt3b3 = sqrt3/3.0;
		double rt6b6 = (sqrt2*sqrt3)/6.0;
		double rt23  = sqrt2/sqrt3;
		if (code == kFCC3D111_a)
		{
			Q(0,0) =-rt2b2;
			Q(0,1) =-rt6b6;
			Q(0,2) = rt3b3;
			
			Q(1,0) = rt2b2;
			Q(1,1) =-rt6b6;
			Q(1,2) = rt3b3;
			
			Q(2,0) = 0.0;
			Q(2,1) = rt23;
			Q(2,2) = rt3b3;
		}
		else /* kFCC3D111_b */
		{
			Q(0,0) = rt2b2;
			Q(0,1) = rt6b6;
			Q(0,2) = rt3b3;
			
			Q(1,0) =-rt2b2;
			Q(1,1) = rt6b6;
			Q(1,2) = rt3b3;
			
			Q(2,0) = 0.0;
			Q(2,1) =-rt23;
			Q(2,2) = rt3b3;
		}
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized orientation \"%s\"", list.Name().Pointer());
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* initialize bond table values */
/* For surface cluster corresponding to atoms exactly on a free surface */
/* Free surface corresponds to "left" surface */
void FCCLatticeT_Surf::LoadBondTable(void)
{
	/* dimension work space */
	int num_bonds = AtomsInShells(fNumShells);
	fBondCounts.Dimension(num_bonds);
	fDefLength.Dimension(num_bonds);
	fBonds.Dimension(num_bonds, 3);

	/* initialize */
  	fBondCounts = 1;
  	fDefLength = 0.0; 

  	double bonddata1[8*3] =
	{ 1.0/sqrt2, 1.0/sqrt2,       0.0,
	  1.0/sqrt2,-1.0/sqrt2,       0.0,
	  1.0/sqrt2,       0.0, 1.0/sqrt2,
	  1.0/sqrt2,       0.0,-1.0/sqrt2,
	        0.0, 1.0/sqrt2, 1.0/sqrt2,
	        0.0,-1.0/sqrt2, 1.0/sqrt2,
	        0.0, 1.0/sqrt2,-1.0/sqrt2,
			0.0,-1.0/sqrt2,-1.0/sqrt2};

  	double bonddata2[5*3] =
	{sqrt2,   0.0,   0.0,
	   0.0, sqrt2,   0.0,
	   0.0,   0.0, sqrt2,
	   0.0,-sqrt2,   0.0,
	   0.0,   0.0,-sqrt2};

  	double bonddata3[12*3] =
	{     sqrt2, 1.0/sqrt2, 1.0/sqrt2,
      1.0/sqrt2,     sqrt2, 1.0/sqrt2,
      1.0/sqrt2, 1.0/sqrt2,     sqrt2,
          sqrt2, 1.0/sqrt2,-1.0/sqrt2,
      1.0/sqrt2,     sqrt2,-1.0/sqrt2,
      1.0/sqrt2, 1.0/sqrt2,    -sqrt2,
          sqrt2,-1.0/sqrt2, 1.0/sqrt2,
      1.0/sqrt2,    -sqrt2, 1.0/sqrt2,
      1.0/sqrt2,-1.0/sqrt2,     sqrt2,
          sqrt2,-1.0/sqrt2,-1.0/sqrt2,
      1.0/sqrt2,    -sqrt2,-1.0/sqrt2,
      1.0/sqrt2,-1.0/sqrt2,    -sqrt2};

	double* shells[3];
	shells[0] = bonddata1;
	shells[1] = bonddata2;
	shells[2] = bonddata3;

	int bond = 0;
	for (int i = 0; i < fNumShells; i++)
	{
		dArray2DT bonds(atoms_per_shell[i], 3, shells[i]);
		for (int j = 0; j < bonds.MajorDim(); j++)
		{
			fBonds(bond,0) = bonds(j,0);
			fBonds(bond,1) = bonds(j,1);
			fBonds(bond,2) = bonds(j,2);
			bond++;
		}
	}
}
