#ifndef GRADATION_H
#define GRADATION_H

#include "realtypes.h"
#include <vector>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>

namespace dem {

  class Gradation {
    
  public:
    Gradation()
      :sievenum(),percent(),size(),ptcl_ratio_ba(),ptcl_ratio_ca()
      {}
    
    Gradation(int sn, std::vector<REAL> v1, std::vector<REAL> v2, REAL ba, REAL ca)
      :sievenum(sn), percent(v1), size(v2), ptcl_ratio_ba(ba), ptcl_ratio_ca(ca)
      {}
    
    int getSieveNum() const { return sievenum; }
    std::vector<REAL>& getPercent() {return percent;}
    const std::vector<REAL>& getPercent() const {return percent;}
    std::vector<REAL>& getSize() {return size;}
    const std::vector<REAL>& getSize() const {return size;}
    REAL getPtclRatioBA() const {return ptcl_ratio_ba;}
    REAL getPtclRatioCA() const {return ptcl_ratio_ca;}
    void setPtclRatioBA(REAL ba) {ptcl_ratio_ba = ba;}
    void setPtclRatioCA(REAL ca) {ptcl_ratio_ca = ca;}

    REAL getPtclMaxRadius() const {return size[0];}
    REAL getPtclMinRadius() const {return size[sievenum-1] * ptcl_ratio_ca;}

  private:
    int  sievenum; // sievenum == percent.size() == size.size()
    std::vector<REAL> percent;
    std::vector<REAL> size;
    REAL ptcl_ratio_ba; // ratio of radius b to radius a
    REAL ptcl_ratio_ca; // ratio of radius c to radius a
    
    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & sievenum;
      ar & percent;
      ar & size;
      ar & ptcl_ratio_ba;
      ar & ptcl_ratio_ca;
    }
    
  };    
  
} // namespace dem

#endif
