/* $Id */
#include "ParadynEAMT.h"
#include "FBD_EAMGlue.h"
#include "FCCLatticeT.h"
#include "ofstreamT.h"
#include "ifstreamT.h"
#include "StringT.h"
#include "C1FunctionT.h"

using namespace Tahoe;

int main(void)
{
	/* output file */
	ofstreamT out;
	StringT parameter_file = "auu3";
	
	/* potentials */
	ParadynEAMT ParadynEAM(parameter_file);
	dMatrixT Q(3);
	Q.Identity();
	FCCLatticeT fcc_lattice(Q, 5);
	ifstreamT in(parameter_file);
	FBD_EAMGlue eam_glue(fcc_lattice, in);

	/* plot range */
	double r_min = 2.0;
	double r_max = 10.0;
	double dr = 0.05;

	/* pair force */
	out.open("pairforce.out");
	EAMPropertyT::PairEnergyFunction pair_z = ParadynEAM.getPairEnergy();
	EAMPropertyT::PairForceFunction Dpair_z = ParadynEAM.getPairForce();
	const C1FunctionT* pair_potential = eam_glue.PairPotential();
	for (double r = r_min; r < r_max; r += dr) {
		out << r << " " 
		    << 2.0*pair_z(r, NULL, NULL)*Dpair_z(r, NULL, NULL)/r << " " 
		    << pair_potential->DFunction(r) << '\n';
	}
	out.close();

	/* electron density */
	out.open("electrondensity.out");
	EAMPropertyT::EDEnergyFunction e_density_1 = ParadynEAM.getElecDensEnergy();
	const C1FunctionT* e_density_2 = eam_glue.ElectronDensity();
	for (double r = r_min; r < r_max; r += dr) {
		out << r << " " << e_density_1(r, NULL, NULL) 
		         << " " << e_density_2->Function(r) << '\n';
	}
	out.close();

	return 0;
}

