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
    density = 0;
    voidRatio = 0;
    angle = 0;
    Jacobian = 0;
    JDot = 0;

    for (int i = 0; i < 3; ++i) {
      coord[i] = 0;
      spin[i] = 0;
      stressEigenValue[i] = 0;
      stressRateEigenValue[i] = 0;
      rateOfDeformEigenValue[i] = 0;
      unitVec[i] = 0;
      OldroEigenValue[i] = 0;
      TruesEigenValue[i] = 0;
      EDotEigenValue[i] = 0;
      eDotEigenValue[i] = 0;
    }

    for (int i = 0; i < 6; ++i) {
      fabric[i] = 0;
      stress[i] = 0;
      stressRate[i] = 0;
      stretch[i] = 0;
      rateOfDeform[i] = 0;
      greenStrain[i] = 0;
      eulerStrain[i] = 0;
      EDot[i] = 0;
      eDot[i] = 0;
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
      OldroEigenVector[i] = 0;
      TruesEigenVector[i] = 0;
      EDotEigenVector[i] = 0;
      eDotEigenVector[i] = 0;
    }

  } 

  void Stress::print(std::ostream &ofs) const {
    double xx, yy, zz, xy, xz, yz;
    xx = greenStrain[0];
    yy = greenStrain[1];
    zz = greenStrain[2];
    xy = greenStrain[3];
    xz = greenStrain[4];
    yz = greenStrain[5];
    double greenVolumeStrain = (xx + yy + zz) / 3;
    double greenShearStrain = sqrt(2.0)/3 * sqrt(pow(xx-yy,2) + pow(yy-zz,2) + pow(zz-xx,2) + 6.0*(xy*xy + yz*yz + xz*xz));

    xx = eulerStrain[0];
    yy = eulerStrain[1];
    zz = eulerStrain[2];
    xy = eulerStrain[3];
    xz = eulerStrain[4];
    yz = eulerStrain[5];
    double eulerVolumeStrain = (xx + yy + zz) / 3;
    double eulerShearStrain = sqrt(2.0)/3 * sqrt(pow(xx-yy,2) + pow(yy-zz,2) + pow(zz-xx,2) + 6.0*(xy*xy + yz*yz + xz*xz));

    ofs << std::setw(OWID) << coord[0] << std::setw(OWID) << coord[1] << std::setw(OWID) << coord[2]

	<< std::setw(OWID) << density << std::setw(OWID) << voidRatio

	<< std::setw(OWID) << fabric[0] << std::setw(OWID) << fabric[1] << std::setw(OWID) << fabric[2]
	<< std::setw(OWID) << fabric[3] << std::setw(OWID) << fabric[4] << std::setw(OWID) << fabric[5]

	<< std::setw(OWID) << stress[0] << std::setw(OWID) << stress[1] << std::setw(OWID) << stress[2]
	<< std::setw(OWID) << stress[3] << std::setw(OWID) << stress[4] << std::setw(OWID) << stress[5]

	<< std::setw(OWID) << stressRate[0] << std::setw(OWID) << stressRate[1] << std::setw(OWID) << stressRate[2]
	<< std::setw(OWID) << stressRate[3] << std::setw(OWID) << stressRate[4] << std::setw(OWID) << stressRate[5]

	<< std::setw(OWID) << OldroStressRate[0] << std::setw(OWID) << OldroStressRate[1] << std::setw(OWID) << OldroStressRate[2]
	<< std::setw(OWID) << OldroStressRate[3] << std::setw(OWID) << OldroStressRate[4] << std::setw(OWID) << OldroStressRate[5]
	<< std::setw(OWID) << OldroStressRate[6] << std::setw(OWID) << OldroStressRate[7] << std::setw(OWID) << OldroStressRate[8]

	<< std::setw(OWID) << TruesStressRate[0] << std::setw(OWID) << TruesStressRate[1] << std::setw(OWID) << TruesStressRate[2]
	<< std::setw(OWID) << TruesStressRate[3] << std::setw(OWID) << TruesStressRate[4] << std::setw(OWID) << TruesStressRate[5]
	<< std::setw(OWID) << TruesStressRate[6] << std::setw(OWID) << TruesStressRate[7] << std::setw(OWID) << TruesStressRate[8]

	<< std::setw(OWID) << deformGradient[0] << std::setw(OWID) << deformGradient[1] << std::setw(OWID) << deformGradient[2]
	<< std::setw(OWID) << deformGradient[3] << std::setw(OWID) << deformGradient[4] << std::setw(OWID) << deformGradient[5]
	<< std::setw(OWID) << deformGradient[6] << std::setw(OWID) << deformGradient[7] << std::setw(OWID) << deformGradient[8]

	<< std::setw(OWID) << rotation[0] << std::setw(OWID) << rotation[1] << std::setw(OWID) << rotation[2]
	<< std::setw(OWID) << rotation[3] << std::setw(OWID) << rotation[4] << std::setw(OWID) << rotation[5]
	<< std::setw(OWID) << rotation[6] << std::setw(OWID) << rotation[7] << std::setw(OWID) << rotation[8]

	<< std::setw(OWID) << stretch[0] << std::setw(OWID) << stretch[1] << std::setw(OWID) << stretch[2]
	<< std::setw(OWID) << stretch[3] << std::setw(OWID) << stretch[4] << std::setw(OWID) << stretch[5]

	<< std::setw(OWID) << greenStrain[0] << std::setw(OWID) << greenStrain[1] << std::setw(OWID) << greenStrain[2]
	<< std::setw(OWID) << greenStrain[3] << std::setw(OWID) << greenStrain[4] << std::setw(OWID) << greenStrain[5]
	<< std::setw(OWID) << greenVolumeStrain << std::setw(OWID) << greenShearStrain

	<< std::setw(OWID) << eulerStrain[0] << std::setw(OWID) << eulerStrain[1] << std::setw(OWID) << eulerStrain[2]
	<< std::setw(OWID) << eulerStrain[3] << std::setw(OWID) << eulerStrain[4] << std::setw(OWID) << eulerStrain[5]
	<< std::setw(OWID) << eulerVolumeStrain << std::setw(OWID) << eulerShearStrain

	<< std::setw(OWID) << Jacobian

	<< std::setw(OWID) << velocityGradient[0] << std::setw(OWID) << velocityGradient[1] << std::setw(OWID) << velocityGradient[2]
	<< std::setw(OWID) << velocityGradient[3] << std::setw(OWID) << velocityGradient[4] << std::setw(OWID) << velocityGradient[5]
	<< std::setw(OWID) << velocityGradient[6] << std::setw(OWID) << velocityGradient[7] << std::setw(OWID) << velocityGradient[8]

	<< std::setw(OWID) << rateOfDeform[0] << std::setw(OWID) << rateOfDeform[1] << std::setw(OWID) << rateOfDeform[2]
	<< std::setw(OWID) << rateOfDeform[3] << std::setw(OWID) << rateOfDeform[4] << std::setw(OWID) << rateOfDeform[5]

	<< std::setw(OWID) << spin[0] << std::setw(OWID) << spin[1] << std::setw(OWID) << spin[2]

	<< std::setw(OWID) << norm[0] << std::setw(OWID) << norm[1] << std::setw(OWID) << norm[2]
	<< std::setw(OWID) << norm[3] << std::setw(OWID) << norm[4] << std::setw(OWID) << norm[5]
	<< std::setw(OWID) << norm[6] << std::setw(OWID) << norm[7] << std::setw(OWID) << norm[8]

	<< std::setw(OWID) << (stressEigenValue[0] + stressEigenValue[1] + stressEigenValue[2]) / 3
	<< std::setw(OWID) << sqrt((pow(stressEigenValue[0]-stressEigenValue[1],2) + pow(stressEigenValue[1]-stressEigenValue[2],2) + pow(stressEigenValue[2]-stressEigenValue[0],2)) / 2)

	<< std::setw(OWID) << stressEigenValue[0] << std::setw(OWID) << stressEigenValue[1] << std::setw(OWID) << stressEigenValue[2]
	<< std::setw(OWID) << stressEigenVector[0] << std::setw(OWID) << stressEigenVector[1] << std::setw(OWID) << stressEigenVector[2]
	<< std::setw(OWID) << stressEigenVector[3] << std::setw(OWID) << stressEigenVector[4] << std::setw(OWID) << stressEigenVector[5]
	<< std::setw(OWID) << stressEigenVector[6] << std::setw(OWID) << stressEigenVector[7] << std::setw(OWID) << stressEigenVector[8]

	<< std::setw(OWID) << stressRateEigenValue[0] << std::setw(OWID) << stressRateEigenValue[1] << std::setw(OWID) << stressRateEigenValue[2]
	<< std::setw(OWID) << stressRateEigenVector[0] << std::setw(OWID) << stressRateEigenVector[1] << std::setw(OWID) << stressRateEigenVector[2]
	<< std::setw(OWID) << stressRateEigenVector[3] << std::setw(OWID) << stressRateEigenVector[4] << std::setw(OWID) << stressRateEigenVector[5]
	<< std::setw(OWID) << stressRateEigenVector[6] << std::setw(OWID) << stressRateEigenVector[7] << std::setw(OWID) << stressRateEigenVector[8]

	<< std::setw(OWID) << rateOfDeformEigenValue[0] << std::setw(OWID) << rateOfDeformEigenValue[1] << std::setw(OWID) << rateOfDeformEigenValue[2]
	<< std::setw(OWID) << rateOfDeformEigenVector[0] << std::setw(OWID) << rateOfDeformEigenVector[1] << std::setw(OWID) << rateOfDeformEigenVector[2]
	<< std::setw(OWID) << rateOfDeformEigenVector[3] << std::setw(OWID) << rateOfDeformEigenVector[4] << std::setw(OWID) << rateOfDeformEigenVector[5]
	<< std::setw(OWID) << rateOfDeformEigenVector[6] << std::setw(OWID) << rateOfDeformEigenVector[7] << std::setw(OWID) << rateOfDeformEigenVector[8]

	<< std::setw(OWID) << unitVec[0]*angle << std::setw(OWID) << unitVec[1]*angle << std::setw(OWID) << unitVec[2]*angle
	<< std::setw(OWID) << angle

	<< std::setw(OWID) << - spin[2] << std::setw(OWID) << spin[1] << std::setw(OWID) << -spin[0]
	<< std::setw(OWID) << sqrt(spin[0]*spin[0] + spin[1]*spin[1] + spin[2]*spin[2])

	<< std::setw(OWID) << EDot[0] << std::setw(OWID) << EDot[1] << std::setw(OWID) << EDot[2]
	<< std::setw(OWID) << EDot[3] << std::setw(OWID) << EDot[4] << std::setw(OWID) << EDot[5]

	<< std::setw(OWID) << eDot[0] << std::setw(OWID) << eDot[1] << std::setw(OWID) << eDot[2]
	<< std::setw(OWID) << eDot[3] << std::setw(OWID) << eDot[4] << std::setw(OWID) << eDot[5]
	
	<< std::setw(OWID) << JDot

	<< std::setw(OWID) << OldroEigenValue[0] << std::setw(OWID) << OldroEigenValue[1] << std::setw(OWID) << OldroEigenValue[2]
	<< std::setw(OWID) << OldroEigenVector[0] << std::setw(OWID) << OldroEigenVector[1] << std::setw(OWID) << OldroEigenVector[2]
	<< std::setw(OWID) << OldroEigenVector[3] << std::setw(OWID) << OldroEigenVector[4] << std::setw(OWID) << OldroEigenVector[5]
	<< std::setw(OWID) << OldroEigenVector[6] << std::setw(OWID) << OldroEigenVector[7] << std::setw(OWID) << OldroEigenVector[8]

	<< std::setw(OWID) << TruesEigenValue[0] << std::setw(OWID) << TruesEigenValue[1] << std::setw(OWID) << TruesEigenValue[2]
	<< std::setw(OWID) << TruesEigenVector[0] << std::setw(OWID) << TruesEigenVector[1] << std::setw(OWID) << TruesEigenVector[2]
	<< std::setw(OWID) << TruesEigenVector[3] << std::setw(OWID) << TruesEigenVector[4] << std::setw(OWID) << TruesEigenVector[5]
	<< std::setw(OWID) << TruesEigenVector[6] << std::setw(OWID) << TruesEigenVector[7] << std::setw(OWID) << TruesEigenVector[8]

	<< std::setw(OWID) << EDotEigenValue[0] << std::setw(OWID) << EDotEigenValue[1] << std::setw(OWID) << EDotEigenValue[2]
	<< std::setw(OWID) << EDotEigenVector[0] << std::setw(OWID) << EDotEigenVector[1] << std::setw(OWID) << EDotEigenVector[2]
	<< std::setw(OWID) << EDotEigenVector[3] << std::setw(OWID) << EDotEigenVector[4] << std::setw(OWID) << EDotEigenVector[5]
	<< std::setw(OWID) << EDotEigenVector[6] << std::setw(OWID) << EDotEigenVector[7] << std::setw(OWID) << EDotEigenVector[8]

	<< std::setw(OWID) << eDotEigenValue[0] << std::setw(OWID) << eDotEigenValue[1] << std::setw(OWID) << eDotEigenValue[2]
	<< std::setw(OWID) << eDotEigenVector[0] << std::setw(OWID) << eDotEigenVector[1] << std::setw(OWID) << eDotEigenVector[2]
	<< std::setw(OWID) << eDotEigenVector[3] << std::setw(OWID) << eDotEigenVector[4] << std::setw(OWID) << eDotEigenVector[5]
	<< std::setw(OWID) << eDotEigenVector[6] << std::setw(OWID) << eDotEigenVector[7] << std::setw(OWID) << eDotEigenVector[8]

	<< std::endl;
  }

}

#endif
