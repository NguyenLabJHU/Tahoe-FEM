#include "ExceptionT.h"
#include "FSThermoMechMatT.h"
#include "RGSplitT2.h"

namespace Tahoe {

    FSThermoMechMatT::FSThermoMechMatT(void):
    ParameterInterfaceT("FSThermoMechMatT"),
    fSpectralDecompSpat(3),
    fFSThermoMechSupport(NULL)
    {
        /*initialize*/
        fNumProcess = 1;
        fNumR = 1;
    }


    void FSThermoMechMatT::SetFSThermoMechSupport(const FSThermoMechSupportT* support)
    {
      FSSolidMatT::SetFSMatSupport(support);
      fFSThermoMechSupport  =support;
    }

    const dSymMatrixT& FSThermoMechMatT:: d_ij(void)
    {
        fCoupledModulus=0.0;
        return fCoupledModulus;
    }
   
    const dSymMatrixT& FSThermoMechMatT::h_ij(void)
    {
        fhij=0.0;
        return fhij;
    }
    
    const dArrayT& FSThermoMechMatT:: q_i(void)
    {
       const dMatrixT& kij=Conductivity();
        //check this function
        dArrayT fT_i=TemperatureGradient();
        fHeatFlux[0]=kij(0,0)*fT_i[0]+kij(0,1)*fT_i[1]+kij(0,2)*fT_i[2];
        fHeatFlux[1]=kij(1,0)*fT_i[0]+kij(1,1)*fT_i[1]+kij(1,2)*fT_i[2];
        fHeatFlux[2]=kij(2,0)*fT_i[0]+kij(2,1)*fT_i[1]+kij(2,2)*fT_i[2];
        return fHeatFlux;
    }
    
   double FSThermoMechMatT:: b1(void)
    {
        return 0.0;
    }
    
    double FSThermoMechMatT:: heatres(void)
    {
        return 0.0;
    }

    double FSThermoMechMatT:: dtime(void)
    {
        return 0.0;
    }

    const dMatrixT& FSThermoMechMatT::Conductivity(void)
    
    {   fkij=0.0;
        fkij(0,0)=fk;
        fkij(1,1)=fk;
        fkij(2,2)=fk;
        return fkij;
    }
    
    
    const dArrayT FSThermoMechMatT::TemperatureGradient(void)
    {
        fTemperatureGradient=fFSThermoMechSupport->TemperatureGradient();
        return fTemperatureGradient;
    }
    
    const dArrayT FSThermoMechMatT::TemperatureGradient(int ip) 
    {
        fTemperatureGradient=fFSThermoMechSupport->TemperatureGradient(ip);
        return fTemperatureGradient;
    }

	void FSThermoMechMatT::DefineSubs(SubListT& sub_list) const
  {
	//FSSolidMatT::DefineSubs(sub_list);
     RGSplitT2:: DefineSubs(sub_list);
    return;
  }

    void FSThermoMechMatT::Initialize(void)
    {
        int nsd=NumSD();
        fkij.Dimension(nsd);
        fhij.Dimension(nsd);
        fCoupledModulus.Dimension(nsd);
        fHeatFlux.Dimension(nsd);
    }


} //namespace Tahoe
