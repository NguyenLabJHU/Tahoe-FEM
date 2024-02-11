#include "GhostParticleT.h"

using namespace dem;

namespace Tahoe {

GhostParticleT::GhostParticleT(int id, int tp, vec dim, vec position, vec dirca, vec dircb, vec dircc, int elem_group, int elem_num)
    :particle(id, tp, dim, position, dirca, dircb, dircc), ElemGroup(elem_group), ElemNum(elem_num)
{
}

void GhostParticleT::SetParentCoord(dArrayT& Source) 
{
    ParentCoord = Source;
}

dArrayT& GhostParticleT::GetParentCoord() 
{
    return ParentCoord;
}

const dArrayT& GhostParticleT::GetParentCoord() const
{
    return ParentCoord;
}

} /* namespace Tahoe */

