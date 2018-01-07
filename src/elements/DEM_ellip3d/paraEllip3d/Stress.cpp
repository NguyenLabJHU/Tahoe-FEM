#ifdef STRESS_STRAIN

#include "Stress.h"
#include "const.h"
#include <iostream>

namespace dem {

  Stress::Stress()
  //:coord{0,0,0}, stress{0,0,0,0,0,0}, stressRate{0,0,0,0,0,0}, ...
  {
    setZero();
  } 

  void Stress::setZero() {
    for (int i = 0; i < 3; ++i) {
      coord[i] = 0;
      spin[i] = 0;
      stressEigenValue[i] = 0;
      stressRateEigenValue[i] = 0;
      rateOfDeformEigenValue[i] = 0;
    }

    for (int i = 0; i < 6; ++i) {
      stress[i] = 0;
      stressRate[i] = 0;
      stretch[i] = 0;
      rateOfDeform[i] = 0;
      greenStrain[i] = 0;
      eulerStrain[i] = 0;
    }

    for (int i = 0; i < 9; ++i) {
      OldroStressRate[i] = 0;
      TruesStressRate[i] = 0;
      deformGradient[i] = 0;
      rotation[i] = 0;
      velocityGradient[i] = 0;  
      norm[i] = 0;
      stressEigenVector[i] = 0;
      stressRateEigenVector[i] = 0;
      rateOfDeformEigenVector[i] = 0;
    }

  } 

  void Stress::print(std::ostream &ofs) const {
    ofs << std::setw(OWID) << coord[0] << std::setw(OWID) << coord[1] << std::setw(OWID) << coord[2]

	<< std::setw(OWID) << stress[0]<< std::setw(OWID) << stress[1]<< std::setw(OWID) << stress[2]
	<< std::setw(OWID) << stress[3]<< std::setw(OWID) << stress[4]<< std::setw(OWID) << stress[5]

	<< std::setw(OWID) << stressRate[0]<< std::setw(OWID) << stressRate[1]<< std::setw(OWID) << stressRate[2]
	<< std::setw(OWID) << stressRate[3]<< std::setw(OWID) << stressRate[4]<< std::setw(OWID) << stressRate[5]

	<< std::setw(OWID) << OldroStressRate[0]<< std::setw(OWID) << OldroStressRate[1]<< std::setw(OWID) << OldroStressRate[2]
	<< std::setw(OWID) << OldroStressRate[3]<< std::setw(OWID) << OldroStressRate[4]<< std::setw(OWID) << OldroStressRate[5]
	<< std::setw(OWID) << OldroStressRate[6]<< std::setw(OWID) << OldroStressRate[7]<< std::setw(OWID) << OldroStressRate[8]

	<< std::setw(OWID) << TruesStressRate[0]<< std::setw(OWID) << TruesStressRate[1]<< std::setw(OWID) << TruesStressRate[2]
	<< std::setw(OWID) << TruesStressRate[3]<< std::setw(OWID) << TruesStressRate[4]<< std::setw(OWID) << TruesStressRate[5]
	<< std::setw(OWID) << TruesStressRate[6]<< std::setw(OWID) << TruesStressRate[7]<< std::setw(OWID) << TruesStressRate[8]

	<< std::setw(OWID) << deformGradient[0]<< std::setw(OWID) << deformGradient[1]<< std::setw(OWID) << deformGradient[2]
	<< std::setw(OWID) << deformGradient[3]<< std::setw(OWID) << deformGradient[4]<< std::setw(OWID) << deformGradient[5]
	<< std::setw(OWID) << deformGradient[6]<< std::setw(OWID) << deformGradient[7]<< std::setw(OWID) << deformGradient[8]

	<< std::setw(OWID) << rotation[0]<< std::setw(OWID) << rotation[1]<< std::setw(OWID) << rotation[2]
	<< std::setw(OWID) << rotation[3]<< std::setw(OWID) << rotation[4]<< std::setw(OWID) << rotation[5]
	<< std::setw(OWID) << rotation[6]<< std::setw(OWID) << rotation[7]<< std::setw(OWID) << rotation[8]

	<< std::setw(OWID) << stretch[0]<< std::setw(OWID) << stretch[1]<< std::setw(OWID) << stretch[2]
	<< std::setw(OWID) << stretch[3]<< std::setw(OWID) << stretch[4]<< std::setw(OWID) << stretch[5]

	<< std::setw(OWID) << greenStrain[0]<< std::setw(OWID) << greenStrain[1]<< std::setw(OWID) << greenStrain[2]
	<< std::setw(OWID) << greenStrain[3]<< std::setw(OWID) << greenStrain[4]<< std::setw(OWID) << greenStrain[5]

	<< std::setw(OWID) << eulerStrain[0]<< std::setw(OWID) << eulerStrain[1]<< std::setw(OWID) << eulerStrain[2]
	<< std::setw(OWID) << eulerStrain[3]<< std::setw(OWID) << eulerStrain[4]<< std::setw(OWID) << eulerStrain[5]

	<< std::setw(OWID) << velocityGradient[0]<< std::setw(OWID) << velocityGradient[1]<< std::setw(OWID) << velocityGradient[2]
	<< std::setw(OWID) << velocityGradient[3]<< std::setw(OWID) << velocityGradient[4]<< std::setw(OWID) << velocityGradient[5]
	<< std::setw(OWID) << velocityGradient[6]<< std::setw(OWID) << velocityGradient[7]<< std::setw(OWID) << velocityGradient[8]

	<< std::setw(OWID) << rateOfDeform[0]<< std::setw(OWID) << rateOfDeform[1]<< std::setw(OWID) << rateOfDeform[2]
	<< std::setw(OWID) << rateOfDeform[3]<< std::setw(OWID) << rateOfDeform[4]<< std::setw(OWID) << rateOfDeform[5]

	<< std::setw(OWID) << spin[0]<< std::setw(OWID) << spin[1]<< std::setw(OWID) << spin[2]

	<< std::setw(OWID) << norm[0]<< std::setw(OWID) << norm[1]<< std::setw(OWID) << norm[2]
	<< std::setw(OWID) << norm[3]<< std::setw(OWID) << norm[4]<< std::setw(OWID) << norm[5]
	<< std::setw(OWID) << norm[6]<< std::setw(OWID) << norm[7]<< std::setw(OWID) << norm[8]

	<< std::setw(OWID) << (stressEigenValue[0] + stressEigenValue[1] + stressEigenValue[2]) / 3
	<< std::setw(OWID) << sqrt((pow(stressEigenValue[0]-stressEigenValue[1],2) + pow(stressEigenValue[1]-stressEigenValue[2],2) + pow(stressEigenValue[2]-stressEigenValue[0],2)) / 2)

	<< std::setw(OWID) << stressEigenValue[0]<< std::setw(OWID) << stressEigenValue[1]<< std::setw(OWID) << stressEigenValue[2]
	<< std::setw(OWID) << stressEigenVector[0]<< std::setw(OWID) << stressEigenVector[1]<< std::setw(OWID) << stressEigenVector[2]
	<< std::setw(OWID) << stressEigenVector[3]<< std::setw(OWID) << stressEigenVector[4]<< std::setw(OWID) << stressEigenVector[5]
	<< std::setw(OWID) << stressEigenVector[6]<< std::setw(OWID) << stressEigenVector[7]<< std::setw(OWID) << stressEigenVector[8]

	<< std::setw(OWID) << stressRateEigenValue[0]<< std::setw(OWID) << stressRateEigenValue[1]<< std::setw(OWID) << stressRateEigenValue[2]
	<< std::setw(OWID) << stressRateEigenVector[0]<< std::setw(OWID) << stressRateEigenVector[1]<< std::setw(OWID) << stressRateEigenVector[2]
	<< std::setw(OWID) << stressRateEigenVector[3]<< std::setw(OWID) << stressRateEigenVector[4]<< std::setw(OWID) << stressRateEigenVector[5]
	<< std::setw(OWID) << stressRateEigenVector[6]<< std::setw(OWID) << stressRateEigenVector[7]<< std::setw(OWID) << stressRateEigenVector[8]

	<< std::setw(OWID) << rateOfDeformEigenValue[0]<< std::setw(OWID) << rateOfDeformEigenValue[1]<< std::setw(OWID) << rateOfDeformEigenValue[2]
	<< std::setw(OWID) << rateOfDeformEigenVector[0]<< std::setw(OWID) << rateOfDeformEigenVector[1]<< std::setw(OWID) << rateOfDeformEigenVector[2]
	<< std::setw(OWID) << rateOfDeformEigenVector[3]<< std::setw(OWID) << rateOfDeformEigenVector[4]<< std::setw(OWID) << rateOfDeformEigenVector[5]
	<< std::setw(OWID) << rateOfDeformEigenVector[6]<< std::setw(OWID) << rateOfDeformEigenVector[7]<< std::setw(OWID) << rateOfDeformEigenVector[8]

	<< std::endl;
  }

}

#endif
