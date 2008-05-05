#ifndef _DE_MANAGER_H_
#define _DE_MANAGER_H_

#include "assembly.h"
#include "particle.h"

#include "ArrayT.h"
#include "FBC_CardT.h"
#include "ModelManagerT.h"
#include "NodeManagerT.h"

namespace Tahoe {

class DEManagerT: public dem::assembly
{
public:

    /** constructor */
    DEManagerT();
    
    /** destructor */
    virtual ~DEManagerT(void);

    virtual void TakeParameter(char* database, iArray2DT& fGhostElemSet);

    double GetTimeStep() { return TimeStep; };

    int GetNumStep() { return NumStep; };

    int GetNumPrint() { return NumPrint; };

    void MapToParentDomain(ModelManagerT* fModelManager, iArray2DT& fGhostElemSet);

    virtual void Run();

    void GhostForce(ArrayT<FBC_CardT>& fGhostFBC, ModelManagerT* fModelManager, iArray2DT& fGhostElemSet);

    void PrintFE(NodeManagerT* fNodeManager);

    void GhostDisplace(NodeManagerT* fNodeManager, iArray2DT& fGhostElemSet);
    

protected:
    double TimeStep;
    int NumStep;
    int NumPrint;
    std::list<dem::particle*> fGhostParticleList;

    int NumSD;
    int NumElemNodes;
    iArray2DT fElementSet;
    GeometryT::CodeT fGeometryCode;

    int RunTimes;
};

} /* namespace Tahoe */

#endif /*_DE_MANAGER_H_ */
