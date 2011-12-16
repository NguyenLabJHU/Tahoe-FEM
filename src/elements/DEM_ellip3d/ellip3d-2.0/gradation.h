#ifndef GRADATION_H
#define GRADATION_H

#include "realtypes.h"
#include <vector>

namespace dem {

  class gradation{
    
  public:
    gradation()
      :sievenum(),percent(),size(),ptcl_ratio_ba(),ptcl_ratio_ca()
      {}

    gradation(int _sn, std::vector<REAL>& _v1, std::vector<REAL>& _v2, REAL _ba, REAL _ca)
      :sievenum(_sn), percent(_v1), size(_v2), ptcl_ratio_ba(_ba), ptcl_ratio_ca(_ca)
      {}
    
    int getSieveNum() { return sievenum; }
    std::vector<REAL>& getPercent() {return percent;}
    std::vector<REAL>& getSize() {return size;}
    REAL getPtclRatioBA() {return ptcl_ratio_ba;}
    REAL getPtclRatioCA() {return ptcl_ratio_ca;}
    void setPtclRatioBA(REAL ba) {ptcl_ratio_ba = ba;}
    void setPtclRatioCA(REAL ca) {ptcl_ratio_ca = ca;}

    REAL getMaxPtclRadius() {return size[0];}
    REAL getMinPtclRadius() {return size[sievenum-1] * ptcl_ratio_ca;}

  private:
    int  sievenum; // sievenum == percent.size() == size.size()
    std::vector<REAL> percent;
    std::vector<REAL> size;
    REAL ptcl_ratio_ba; // ratio of radius b to radius a
    REAL ptcl_ratio_ca; // ratio of radius c to radius a
  };    

} // namespace dem

#endif
