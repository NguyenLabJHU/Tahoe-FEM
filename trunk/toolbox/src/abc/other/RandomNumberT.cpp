/* $Id: RandomNumberT.cpp,v 1.2 2003-03-14 16:31:06 paklein Exp $ */
#include "RandomNumberT.h"
#include "ifstreamT.h"
#include <math.h>
#include "dArrayT.h"

using namespace Tahoe;

/* line length */
const int kLineLength = 254;

RandomNumberT::RandomNumberT(ifstreamT& in)
{
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
    case kGaussian:
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

//double RandomNumberT::Rand(void)
//{
//  return (this->*randFunc)();
//}

/* a more high-falutin' uniform random number generator */
/* obviously, it's not implemented yet */
double RandomNumberT::UniformRandom(void)
{
  return 0.;
}

/* Return a Gaussian random number with zero mean and unit variance */
double RandomNumberT::GaussianRandom(void)
{  
  return sqrt(-2.0*log((this->*uniformFunc)()))*cos(6.283185308*(this->*uniformFunc)());
}

/* "minimal standard" random generator from Numerical Recipes. It's
 * the one used by Paradyn, too.
 */
double RandomNumberT::ParadynUniformRandom(void)
{

  long k = fseed/12773;
  fseed = fa*(fseed-k*12773)-2836*k;
  if (fseed < 0)
    fseed += frm;

  return 1./frm*fseed;

}

/* set the parameters */
void RandomNumberT::sRand(int seed, int a, int rm)
{

  if (seed == 0) 
    ExceptionT::GeneralFail("RandomNumberT::sRand","0 provided as seed");
  fseed = seed;
  fa = a;
  frm = rm;
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

  double *f = fillArray.Pointer();

  for (int i = 0; i < fillArray.Length(); i++)
    *f++ = Rand();

  return fillArray;

}
