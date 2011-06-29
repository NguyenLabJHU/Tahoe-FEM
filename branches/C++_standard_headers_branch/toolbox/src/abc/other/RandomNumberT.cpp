/* $Id: RandomNumberT.cpp,v 1.7 2005-06-01 20:11:35 cjkimme Exp $ */
#include "RandomNumberT.h"
#include "ifstreamT.h"
#include <math.h>
#include "dArrayT.h"
#include "dArray2DT.h"

using namespace Tahoe;

/* line length */
const int kLineLength = 254;

const int im = 714025;
const int ia = 1366;
const int ic = 150889;

const double rminv = 1.0/double(im);

RandomNumberT::RandomNumberT(ifstreamT& in)
{
#pragma unused(in)
	ExceptionT::GeneralFail("RandomNumberT::RandomNumberT","not implemented yet");
}

RandomNumberT::RandomNumberT(DistributionT type)
{
	randomType = type;
	switch(randomType)
    {
    	case kUniform:
    	{
			randFunc = &RandomNumberT::UniformRandom;
			uniformType = randomType;
			break;
		}
		case kGaussian: //It's paradyn's till I put in another one
		{
			randFunc = &RandomNumberT::GaussianRandom;
			uniformType = kUniform;
			uniformFunc = &RandomNumberT::UniformRandom;
			break;
		}
    	case kParadynUniform:
		{
			randFunc = &RandomNumberT::ParadynUniformRandom;
			uniformType = randomType;
			break;
		}
		case kParadynGaussian:
		{
			randFunc = &RandomNumberT::GaussianRandom;
			uniformType = kParadynUniform;
			uniformFunc = &RandomNumberT::ParadynUniformRandom;
			break;
		}
	}
}

/* random numbers on [0..1) */
double RandomNumberT::UniformRandom(void)
{
	fseed = (fseed*ia + ic) % im;
	return double(fseed)*rminv;
}

/* Return a Gaussian random number with zero mean and unit variance */
double RandomNumberT::GaussianRandom(void)
{  
	static bool is_saved = false;
	static double g_save;
	if (is_saved) {
		is_saved = false;
		return g_save;
	} else {
		double rsq = 0.0, a1, a2, fac;
		do { 
			a1 =  2.0*(this->*uniformFunc)() - 1.;
			a2 =  2.0*(this->*uniformFunc)() - 1.;
			rsq = a1*a1 + a2*a2;
		} while (rsq > 1.0);
		fac = sqrt(-2.0*log(rsq)/rsq);
		g_save = a1*fac;;
		is_saved = true;
		return a2*fac;
	}
}

/* random number generator from Paradyn */
double RandomNumberT::ParadynUniformRandom(void)
{
	dseed = fmod(da*dseed,drm);

	return dseed/drm;

}

/* set the parameters */
void RandomNumberT::sRand(long seed, long a, long rm)
{
	if (seed == 0) 
		ExceptionT::GeneralFail("RandomNumberT::sRand","0 provided as seed");
    fseed = seed;
	fa = a;
	frm = rm;
	dseed = double(seed);
	da = double(a);
	drm = double(rm);

	switch(uniformType)
	{
		case kUniform:
		{
			break ;
		}
		case kParadynUniform:
		{
			break ;
		}
	}

	return ;

}

/* fill an array with random numbers */
dArrayT& RandomNumberT::RandomArray(dArrayT& fillArray)
{

	double* f = fillArray.Pointer();

	for (int i = 0; i < fillArray.Length(); i++)
		*f++ = Rand();

	return fillArray;

}

dArray2DT& RandomNumberT::RandomArray(dArray2DT& fillArray)
{

	double* f = fillArray.Pointer();

	for (int i = 0; i < fillArray.Length(); i++)
		*f++ = Rand();

	return fillArray;

}


long RandomNumberT::RandSeed(void)
{
	if (uniformType == kUniform)
		return fseed;
	else
		return long(dseed);
}
