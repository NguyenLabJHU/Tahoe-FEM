#ifndef _FE_DE_MANAGER_H_
#define _FE_DE_MANAGER_H_

/* base class  */
#include "FEManagerT.h"

/* direct member */
#include "FBC_CardT.h"
#include "DEManagerT.h"

namespace Tahoe {

class FEDEManagerT: public FEManagerT
{
public:

    /** constructor */
    FEDEManagerT(const StringT& input_file, ofstreamT& output, CommunicatorT& comm,
		 const ArrayT<StringT>& argv, TaskT task);
    
    /** destructor */
    virtual ~FEDEManagerT(void);

    /** accept parameter list */
    virtual void TakeParameterList(const ParameterListT& list);

    /** solve all the time sequences */
    virtual void Solve(void);

    /* perform DEM simulation and calculate ghost forces */
    void DemComputeStep(void);

    /* update ghost particles coordinates based on FE mesh deformation  */
    void GhostDisplaceStep(void);

    /** execute the solution procedure */
    virtual ExceptionT::CodeT SolveStep(void);

    virtual void FormRHS(int group, ArrayT<FBC_CardT>& fGhostFBC);

    ArrayT<FBC_CardT>& GetGhostFBC() {return fGhostFBC;};

    DEManagerT& GetDEManager() {return fDEManager;};

protected:

    DEManagerT fDEManager;
    iArray2DT  fGhostElemSet;
    ArrayT<FBC_CardT> fGhostFBC;
};

} /* namespace Tahoe */

#endif /*_FE_DE_MANAGER_H_ */
