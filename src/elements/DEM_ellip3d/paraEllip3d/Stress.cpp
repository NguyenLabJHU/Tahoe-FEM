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
    }

    for (int i = 0; i < 6; ++i) {
      stress[i] = 0;
      stressRate[i] = 0;
      stretch[i] = 0;
      velocityGradient[i] = 0;  
      rateOfDeform[i] = 0;
    }

    for (int i = 0; i < 9; ++i) {
      OldroStressRate[i] = 0;
      TruesStressRate[i] = 0;
      deformGradient[i] = 0;
      rotation[i] = 0;
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

	<< std::setw(OWID) << velocityGradient[0]<< std::setw(OWID) << velocityGradient[1]<< std::setw(OWID) << velocityGradient[2]
	<< std::setw(OWID) << velocityGradient[3]<< std::setw(OWID) << velocityGradient[4]<< std::setw(OWID) << velocityGradient[5]
	<< std::setw(OWID) << velocityGradient[6]<< std::setw(OWID) << velocityGradient[7]<< std::setw(OWID) << velocityGradient[8]

	<< std::setw(OWID) << rateOfDeform[0]<< std::setw(OWID) << rateOfDeform[1]<< std::setw(OWID) << rateOfDeform[2]
	<< std::setw(OWID) << rateOfDeform[3]<< std::setw(OWID) << rateOfDeform[4]<< std::setw(OWID) << rateOfDeform[5]

	<< std::setw(OWID) << spin[0]<< std::setw(OWID) << spin[1]<< std::setw(OWID) << spin[2]
	<< std::endl;
  }

}

#endif
