/* $Id: ErcolessiAdamsAl.cpp,v 1.2.4.1 2002-10-17 04:37:50 paklein Exp $ */
/* created: paklein (12/02/1996)                                          */
/* ErcolessiAdamsAl.cpp                                                   */

#include "ErcolessiAdamsAl.h"
#include <iostream.h> //TEMP
#include "CubicSplineT.h"

/* lattice parameters - angstrom */
//const double kLatticeParameterAl = 4.032; //given by E&A

using namespace Tahoe;

const double kLatticeParameterAl = 4.03515412; //zero stress

/* Constructor */
ErcolessiAdamsAl::ErcolessiAdamsAl(CBLatticeT& lattice):
	EAM(lattice)
{

}

/*
* Unstressed lattice parameter.
*/
double ErcolessiAdamsAl::LatticeParameter(void) const
{
	return kLatticeParameterAl;
}

/**********************************************************************
* Private
**********************************************************************/

/*
* Set the spline data - called by the constructor
*/
void ErcolessiAdamsAl::SetPairPotential(void)
{
	double rdata[17] = {
	      .202111069753385E+01,
.227374953472558E+01,
.252638837191732E+01,
.277902720910905E+01,
.303166604630078E+01,
.328430488349251E+01,
.353694372068424E+01,
.378958255787597E+01,
.404222139506771E+01,
.429486023225944E+01,
.454749906945117E+01,
.480013790664290E+01,
.505277674383463E+01,
.530541558102636E+01,
.555805441821810E+01,
.555807968210182E+01,
.555810494598553E+01};

	double coeffdata[18*4] = {
1.9601647220e+00, -5.7473949602e+00,  0.0000000000e+00,  0.0000000000e+00,
1.9601647220e+00, -5.7473949602e+00,  0.0000000000e+00,  1.0826322480e+01,
6.8272424075e-01, -3.6743799755e+00,  8.2054485672e+00, -8.1107746287e+00,
1.4737082454e-01, -1.0813942069e+00,  2.0581585545e+00, -1.5102172234e+00,
-1.8818823586e-02, -3.3062774345e-01,  9.1353998483e-01, -8.4098547808e-01,
-5.7601190269e-02, -3.0067565749e-02,  2.7614320499e-01, -2.7363970686e-01,
-5.1984649964e-02,  5.7065133239e-02,  6.8747152942e-02, -2.7630034947e-01,
-3.7635248485e-02,  3.8895779795e-02, -1.4066544407e-01, -3.6400424951e-02,
-3.7373787969e-02, -3.9149250943e-02, -1.6825392717e-01,  3.0191247769e-01,
-5.3135103012e-02, -6.6354260524e-02,  6.0570524720e-02,  1.7031110191e-01,
-6.3286498356e-02, -3.1383022794e-03,  1.8965212096e-01, -1.7586429056e-01,
-5.4810362384e-02,  5.9014336145e-02,  5.6361671353e-02, -6.1099800233e-02,
-3.7288923234e-02,  7.5793292790e-02,  1.0053123902e-02, -8.6121577486e-02,
-1.8887651763e-02,  6.4382425731e-02, -5.5219841678e-02,  1.8247395890e-02,
-5.8523936253e-03,  3.9975068334e-02, -4.1389839039e-02, -9.9540787058e-02,
0.0000000000e+00,  1.7218055004e-06, -1.1683344513e-01,  1.9268851408e+03,
0.0000000000e+00, -4.9194442861e-07,  2.9208361287e-02, -3.8537702835e+02,
0.0000000000e+00,  2.4597221430e-07,  0.0000000000e+00,  0.0000000000e+00};

	dArrayT		knots(17, rdata);
	dArray2DT	coeff(18, 4, coeffdata);

	fPairPotential = new CubicSplineT(knots, coeff);
	if (!fPairPotential) throw ExceptionT::kOutOfMemory;
}

void ErcolessiAdamsAl::SetEmbeddingEnergy(void)
{
	double rdata[13] = {
	     .000000000000000e+00,
.100000000000000e+00,
.200000000000000e+00,
.300000000000000e+00,
.400000000000000e+00,
.500000000000000e+00,
.600000000000000e+00,
.700000000000000e+00,
.800000000000000e+00,
.900000000000000e+00,
.100000000000000e+01,
.110000000000000e+01,
.120000000000000e+01};
	
	double coeffdata[14*4] = {
0.0000000000e+00, -1.3579223135e+01,  0.0000000000e+00,  0.0000000000e+00,
0.0000000000e+00, -1.3579223135e+01,  0.0000000000e+00,  2.1838907206e+02,
-1.1395332414e+00, -7.0275509732e+00,  6.5516721617e+01, -2.6997747547e+02,
-1.4570985981e+00, -2.0235309139e+00, -1.5476521023e+01,  6.5083819652e+01,
-1.7491330800e+00, -3.1663205290e+00,  4.0486248721e+00, -4.3243371611e+00,
-2.0296032214e+00, -2.4863256694e+00,  2.7513237238e+00, -1.3006986074e+00,
-2.2520232497e+00, -1.9750818829e+00,  2.3611141416e+00, -1.3102432489e+00,
-2.4272305398e+00, -1.5421663520e+00,  1.9680411669e+00,  1.0046998653e+01,
-2.5517197647e+00, -8.4714815905e-01,  4.9821407627e+00, -1.8603215371e+01,
-2.6052163883e+00, -4.0881646766e-01, -5.9882384874e-01,  8.1073297593e+00,
-2.6439789438e+00, -2.8536134462e-01,  1.8333750791e+00, -2.8975209314e+00,
-2.6570788484e+00, -5.6119567576e-03,  9.6411879962e-01,  2.3573620997e+00,
-2.6456414940e+00,  2.5793266616e-01,  1.6713274295e+00, -5.5710914317e+00,
-2.6087060445e+00,  4.2506540911e-01,  0.0000000000e+00,  0.0000000000e+00};

	dArrayT		knots(13, rdata);
	dArray2DT	coeff(14, 4, coeffdata);

	fEmbeddingEnergy = new CubicSplineT(knots, coeff);
	if (!fEmbeddingEnergy) throw ExceptionT::kOutOfMemory;
}

void ErcolessiAdamsAl::SetElectronDensity(void)
{
	double rdata[17] = {
	      .202111069753385e+01,
.227374953472558e+01,
.252638837191732e+01,
.277902720910905e+01,
.303166604630078e+01,
.328430488349251e+01,
.353694372068424e+01,
.378958255787597e+01,
.404222139506771e+01,
.429486023225944e+01,
.454749906945117e+01,
.480013790664290e+01,
.505277674383463e+01,
.530541558102636e+01,
.555805441821810e+01,
.555807968210182e+01,
.555810494598553e+01};

	double coeffdata[18*4] = {
8.6567462371e-02,  3.6025048447e-02,  0.0000000000e+00,  0.0000000000e+00,
8.6567462371e-02,  3.6025048447e-02,  0.0000000000e+00, -1.9518226033e-01,
9.2521470294e-02, -1.3482841090e-03, -1.4793185787e-01,  2.1466179050e-01,
8.6200312383e-02, -3.4991692861e-02,  1.4763857554e-02, -1.2581353417e-01,
7.6273629275e-02, -5.1622513351e-02, -8.0592297374e-02,  1.5877939335e-01,
6.0648184127e-02, -6.1941059273e-02,  3.9749226544e-02, -5.7886447813e-02,
4.6603095959e-02, -5.2940710128e-02, -4.1238680501e-03,  5.6372567884e-02,
3.3874013885e-02, -4.4230237420e-02,  3.8601831949e-02, -1.1821919537e-01,
2.3257266171e-02, -4.7362104137e-02, -5.0998448209e-02,  1.7785573044e-01,
1.0904640549e-02, -3.9074818039e-02,  8.3801346568e-02, -7.0229646571e-02,
5.2491060568e-03, -1.0179381296e-02,  3.0573137830e-02, -4.4139322443e-02,
3.9170241914e-03, -3.1832170890e-03, -2.8807834598e-03,  9.5396962632e-03,
3.0827777629e-03, -2.8121597607e-03,  4.3495098535e-03, -9.1648251086e-03,
2.5021474535e-03, -2.3693224239e-03, -2.5966624220e-03, -1.6472775667e-02,
1.4722051380e-03, -6.8355509934e-03, -1.5081651091e-02,  7.5493085827e-02,
0.0000000000e+00, -6.2096654791e-07,  4.2135805167e-02, -6.9492820985e+02,
0.0000000000e+00,  1.7741901366e-07, -1.0533951294e-02,  1.3898564204e+02,
0.0000000000e+00, -8.8709506830e-08,  0.0000000000e+00,  0.0000000000e+00};

	dArrayT		knots(17, rdata);
	dArray2DT	coeff(18, 4, coeffdata);

	fElectronDensity = new CubicSplineT(knots, coeff);
	if (!fElectronDensity) throw ExceptionT::kOutOfMemory;
}
