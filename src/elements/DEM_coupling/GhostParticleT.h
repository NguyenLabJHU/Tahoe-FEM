#ifndef _GhostParticle_H_
#define _GhostParticle_H_

#include "vec.h"
#include "particle.h"

#include "dArrayT.h"

namespace Tahoe {

class GhostParticleT: public dem::particle
{
public:

    /** constructor */
    GhostParticleT(int id, int tp, dem::vec dim, dem::vec position, dem::vec dirca, dem::vec dircb, dem::vec dircc, int elem_group, int elem_num);
    
    /** destructor */
    ~GhostParticleT(void);

    int GetElementGroup() { return ElemGroup; };

    int GetElementNum() { return ElemNum; };

    void SetParentCoord(dArrayT& Source);
    
    dArrayT& GetParentCoord();

    const dArrayT& GetParentCoord() const;

protected:

    int ElemGroup;
    int ElemNum;
    dArrayT ParentCoord; // coordinates in parent domain
    
};

} /* namespace Tahoe */

#endif /*_GhostParticle_H_ */
