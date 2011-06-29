//
// $Id: FloatingPointT.cpp,v 1.2 2008-07-14 23:33:53 lxmota Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.1  2008/07/14 17:50:46  lxmota
// Initial sources.
//
//

//
// 2001/10/02 22:45:15 by Jaroslaw Knap
// Imported sources.
//

#include <fenv.h>

#if defined(sgi)
#include <sigfpe.h>
#endif // sgi

#include "FloatingPointT.h"

bool Tahoe::FloatingPointT::active_ = false;
unsigned  Tahoe::FloatingPointT::oldMask_ = Tahoe::emptyMask_;

//
//
//

namespace Tahoe {

  FloatingPointT::FloatingPointT()
  {

    if (active_) return;

    oldMask_ = getCurrentMask();
    active_  = true;
    return;

  }

  FloatingPointT::~FloatingPointT()
  {

    if (!active_) return;

    setMask(oldMask_);
    active_ = false;
    return;

  }


}


//
// all architectures with <fenv.h> available
//

#if !defined(sgi)

//
// fully ISO/IEC C99 compliant
//

#if defined(HAVE_FESETTRAPENABLE)

// set traps

namespace Tahoe {

  void FloatingPointT::trapInexact()
  {

    fesettrapenable(FE_INEXACT);
    return;

  }

  void FloatingPointT::trapDivbyzero()
  {

    fesettrapenable(FE_DIVBYZERO);
    return;

  }

  void FloatingPointT::trapUnderflow()
  {

    fesettrapenable(FE_UNDERFLOW);
    return;

  }

  void FloatingPointT::trapOverflow()
  {

    fesettrapenable(FE_OVERFLOW);
    return;

  }

  void FloatingPointT::trapInvalid()
  {

    fesettrapenable(FE_INVALID);
    return;

  }

  // get current trap mask

  unsigned FloatingPointT::getCurrentMask()
  {

    unsigned currentMask = emptyMask_;
    int currentTraps = fegettrapenable();

    if (currentTraps & FE_INEXACT)   currentMask |= inexactMask_;
    if (currentTraps & FE_DIVBYZERO) currentMask |= divbyzeroMask_;
    if (currentTraps & FE_UNDERFLOW) currentMask |= underflowMask_;
    if (currentTraps & FE_OVERFLOW)  currentMask |= overflowMask_;
    if (currentTraps & FE_INVALID)   currentMask |= invalidMask_;

    return currentMask;

  }

  // set mask

  void FloatingPointT::setMask(unsigned mask)
  {

    int currentTraps = 0;

    if (mask & inexactMask_)   currentTraps |= FE_INEXACT;
    if (mask & divbyzeroMask_) currentTraps |= FE_DIVBYZERO;
    if (mask & underflowMask_) currentTraps |= FE_UNDERFLOW;
    if (mask & overflowMask_)  currentTraps |= FE_OVERFLOW;
    if (mask & invalidMask_)   currentTraps |= FE_INVALID;

    fesettrapenable(currentTraps);

    return;

  }

}

#else

#if defined(linux)

//
// subset of ISO/IEC C99; (linux)
//

void Tahoe::FloatingPointT::trapInexact()
{
  feenableexcept( fegetexcept() | FE_INEXACT );
  return;
}

void Tahoe::FloatingPointT::trapDivbyzero()
{
  feenableexcept( fegetexcept() | FE_DIVBYZERO );
  return;
}

void Tahoe::FloatingPointT::trapUnderflow()
{
  feenableexcept( fegetexcept() | FE_UNDERFLOW );
  return;
}

void Tahoe::FloatingPointT::trapOverflow()
{
  feenableexcept( fegetexcept() | FE_OVERFLOW );
  return;
}

void Tahoe::FloatingPointT::trapInvalid()
{
  feenableexcept( fegetexcept() | FE_INVALID );
  return;
}

unsigned Tahoe::FloatingPointT::getCurrentMask()
{
  return fegetexcept();
}

// set mask

void Tahoe::FloatingPointT::setMask(unsigned mask)
{
  feenableexcept( mask );
  return;
}

#else

//
// dummy interfaces
//

void Tahoe::FloatingPointT::trapInexact()
{

  return;

}

void Tahoe::FloatingPointT::trapDivbyzero()
{

  return;

}

void Tahoe::FloatingPointT::trapUnderflow()
{

  return;

}

void Tahoe::FloatingPointT::trapOverflow()
{

  return;

}

void Tahoe::FloatingPointT::trapInvalid()
{

  return;

}

unsigned Tahoe::FloatingPointT::getCurrentMask()
{

  return emptyMask_;

}

// set mask

void Tahoe::FloatingPointT::setMask(unsigned mask)
{

  return;

}

#endif // linux

#endif // HAVE_FESETTRAPENABLE

#endif // !sgi

//
// IRIX
//

#if defined(sgi)


void Tahoe::FloatingPointT::trapInexact()
{

  // not available on IRIX
  return;

}

void Tahoe::FloatingPointT::trapDivbyzero()
{

  sigfpe_[_DIVZERO].abort=1;
  handle_sigfpes(_ON, _EN_DIVZERO, 0, _ABORT_ON_ERROR, 0);
  return;

}

void Tahoe::FloatingPointT::trapUnderflow()
{

  sigfpe_[_UNDERFL].abort=1;
  handle_sigfpes(_ON, _EN_UNDERFL, 0, _ABORT_ON_ERROR, 0);
  return;

}

void Tahoe::FloatingPointT::trapOverflow()
{

  sigfpe_[_OVERFL].abort=1;
  handle_sigfpes(_ON, _EN_OVERFL, 0, _ABORT_ON_ERROR, 0);
  return;

}

void Tahoe::FloatingPointT::trapInvalid()
{

  sigfpe_[_INVALID].abort=1;
  handle_sigfpes(_ON, _EN_INVALID, 0, _ABORT_ON_ERROR, 0);
  return;

}

unsigned Tahoe::FloatingPointT::getCurrentMask()
{

  return emptyMask_;

}

// set mask

void Tahoe::FloatingPointT::setMask(unsigned mask)
{

  return;

}


#endif // sgi

