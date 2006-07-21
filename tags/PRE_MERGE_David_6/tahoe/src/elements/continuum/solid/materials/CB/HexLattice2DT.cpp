/* $Id: HexLattice2DT.cpp,v 1.6 2005-08-30 07:53:40 jzimmer Exp $ */
#include "HexLattice2DT.h"
#include "ParameterContainerT.h"

using namespace Tahoe;

/* number of atoms per shell */
const int atoms_per_shell[] = {3, 3, 3, 6, 3};
const int atoms_in_shells[] = {3, 6, 9, 15, 18};
static int AtomsInShells(int nshells) {
	if (nshells < 0 || nshells > 5) ExceptionT::OutOfRange();
	return atoms_in_shells[nshells-1];
};
const double sqrt2 = sqrt(2.0);
const double sqrt3 = sqrt(3.0);

/* constructor */
HexLattice2DT::HexLattice2DT(int nshells):
        ParameterInterfaceT("CB_lattice_Hex"),
	fNumShells(nshells)
{

}

/* information about subordinate parameter lists */
void HexLattice2DT::DefineSubs(SubListT& sub_list) const
{
        /* inherited */
        ParameterInterfaceT::DefineSubs(sub_list);

        sub_list.AddSub("Hex_lattice_orientation", ParameterListT::Once, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* HexLattice2DT::NewSub(const StringT& name) const
{
        if (name == "Hex_lattice_orientation")
        {
                ParameterContainerT* orientation = new ParameterContainerT(name);
                orientation->SetListOrder(ParameterListT::Choice);

                ParameterContainerT natural("HEX2D_natural");
                orientation->AddSub(natural);

                ParameterContainerT HEX2D90("HEX2D_90");
                orientation->AddSub(HEX2D90);

                ParameterContainerT Rotation_angle("HEX2D_Rotation_angle");
                Rotation_angle.AddParameter(ParameterT::Double, "phi");
                orientation->AddSub(Rotation_angle);

                return orientation;
        }
        else /* inherited */
                return ParameterInterfaceT::NewSub(name);
}

/* accept parameter list */
void HexLattice2DT::TakeParameterList(const ParameterListT& list)
{
        const char caller[] = "HexLattice2DT::TakeParameterList";

        /* inherited */
        ParameterInterfaceT::TakeParameterList(list);

        /* set Q */
        const ParameterListT& orientation = list.GetListChoice(*this, "Hex_lattice_orientation");
        dMatrixT Q;
        SetQ(orientation, Q);

        /* initialize bond table */
        Initialize(&Q);
}

/* set the transformation matrix for the given orientation */
void HexLattice2DT::SetQ(const ParameterListT& list, dMatrixT& Q)
{
        const char caller[] = "HexLattice2DT::SetQ";

        /* dimension */
        Q.Dimension(2);
        Q = 0.0;

        /* extract orientation */
        if (list.Name() == "HEX2D_natural")
                Q.Identity();
        else if (list.Name() == "HEX2D_90")
        {
                Q(0,1) = 1.0;
		Q(1,0) = -1.0;
        }
        else if (list.Name() == "HEX2D_Rotation_angle")
	{
        	double phi = list.GetParameter("phi");
		double cosphi = cos(phi);
		double sinphi = sin(phi);
                Q(0,0) = cosphi; Q(0,1) = -sinphi;
		Q(1,0) = sinphi; Q(1,1) = cosphi;
	}
	else
                ExceptionT::GeneralFail(caller, "unrecognized orientation \"%s\"", list.Name().Pointer());
}

/* initialize bond table values */
void HexLattice2DT::LoadBondTable(void)
{
	/* dimension work space */
	int num_bonds = AtomsInShells(fNumShells);
	fBondCounts.Dimension(num_bonds);
	fDefLength.Dimension(num_bonds);
	fBonds.Dimension(num_bonds, 2);

	/* initialize */
  	fBondCounts = 1;
  	fDefLength = 0.0; 
  
  	double bonddata1[3*2] =
  		{ 1.0, 0.0,
  		  0.5, sqrt3/2.0,
  		 -0.5, sqrt3/2.0};

  	double bonddata2[3*2] =
  		{ 0.0, sqrt3,
  		  1.5, sqrt3/2.0,
  		 -1.5, sqrt3/2.0};

  	double bonddata3[3*2] =
  		{ 2.0, 0.0,
  		  1.0, sqrt3,
  		 -1.0, sqrt3};

  	double bonddata4[6*2] =
  		{ 2.5, sqrt3/2.0,
  		  2.0, sqrt3,
  		  0.5, 3.0*sqrt3/2.0,
  		 -0.5, 3.0*sqrt3/2.0,
  		 -2.0, sqrt3,
  		 -2.5, sqrt3/2.0};

  	double bonddata5[3*2] =
  		{ 3.0, 0.0,
  		  1.5, 3.0*sqrt3/2.0,
  		 -1.5, 3.0*sqrt3/2.0};

	double* shells[5];
	shells[0] = bonddata1;
	shells[1] = bonddata2;
	shells[2] = bonddata3;
	shells[3] = bonddata4;
	shells[4] = bonddata5;

	int bond = 0;
	for (int i = 0; i < fNumShells; i++)
	{
		dArray2DT bonds(atoms_per_shell[i], 2, shells[i]);
		for (int j = 0; j < bonds.MajorDim(); j++)
		{
			fBonds(bond,0) = bonds(j,0);
			fBonds(bond,1) = bonds(j,1);
			bond++;
		}
	}
}
