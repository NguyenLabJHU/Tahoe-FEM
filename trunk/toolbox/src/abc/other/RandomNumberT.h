/* $Id: RandomNumberT.h,v 1.2 2003-03-14 16:31:06 paklein Exp $ */
#ifndef _RANDOM_NUMBER_T_H_
#define _RANDOM_NUMBER_T_H_

/* Environmental */
//#include "Environment.h"

//#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class dArrayT;
class ifstreamT;

class RandomNumberT
{
public:
  enum DistributionT 
  {
    kUniform = 0,
    kGaussian = 1,
    kParadynUniform = 2,
    kParadynGaussian = 3
  };

	/** \name constructors */
	/*@{*/
	RandomNumberT(DistributionT type = kUniform);
	/*@}*/

	/*@{*/
	RandomNumberT(ifstreamT& in);
	/*@}*/
	
 public:
	void sRand(int seed, int a = 16807, int rm = 2147483647);

	double Rand(void);

	dArrayT& RandomArray(dArrayT& fillArray);

        // Need some way to write seed to restart. 

 private:

	double UniformRandom(void);

	double GaussianRandom(void);

	double ParadynUniformRandom(void);

	long fseed, fa, frm; 
	DistributionT randomType, uniformType;

	double (RandomNumberT::*randFunc)(void);
	double (RandomNumberT::*uniformFunc)(void);

};

/* Is this fast or slow? */
inline double RandomNumberT::Rand(void)
{
  return (this->*randFunc)();
}

}//namespace Tahoe
#endif /* _RANDOM_NUMBER_T_H_ */
