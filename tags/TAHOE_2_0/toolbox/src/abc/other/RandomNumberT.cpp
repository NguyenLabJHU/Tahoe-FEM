/* $Id: RandomNumberT.cpp,v 1.6 2003-04-24 20:49:06 cjkimme Exp $ */
#include "RandomNumberT.h"
#include "ifstreamT.h"
#include <math.h>
#include "dArrayT.h"
#include "dArray2DT.h"

using namespace Tahoe;

/* line length */
const int kLineLength = 254;
const double TWOPI = 6.283185307179586476925286;

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

/* Implemeneted in .h file for now */
//double RandomNumberT::Rand(void)
//{
//  return (this->*randFunc)();
//}

/* a slightly more high-falutin' uniform random number generator */
double RandomNumberT::UniformRandom(void)
{
	long k = fseed/12773;
	fseed = fa*(fseed-k*12773)-2836*k;
	if (fseed < 0)
		fseed += frm;

	return 1./frm*fseed;
}

/* Return a Gaussian random number with zero mean and unit variance */
double RandomNumberT::GaussianRandom(void)
{  
	double a1 = (this->*uniformFunc)();
	double a2 = (this->*uniformFunc)();
	return sqrt(-2.0*log(a1))*cos(TWOPI*a2);
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
