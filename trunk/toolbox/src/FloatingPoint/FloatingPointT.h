//
// $Id: FloatingPointT.h,v 1.1 2008-07-14 17:50:46 lxmota Exp $
//
// $Log: not supported by cvs2svn $
//

//
// 2001/10/02 22:45:15 by Jaroslaw Knap
// Imported sources.
//

#if !defined(_FloatingPointT_h_)
#define _FloatingPointT_h_

namespace Tahoe {

  const unsigned emptyMask_     = 0x00;
  const unsigned inexactMask_   = 0x01;
  const unsigned divbyzeroMask_ = 0x02;
  const unsigned underflowMask_ = 0x04;
  const unsigned overflowMask_  = 0x08;
  const unsigned invalidMask_   = 0x10;

  class FloatingPointT {

  public:

    FloatingPointT();
    virtual ~FloatingPointT();

    void trapInexact();
    void trapDivbyzero();
    void trapUnderflow();
    void trapOverflow();
    void trapInvalid();


  private:

    static bool active_;
    static unsigned mask_;
    static unsigned oldMask_;

    // obtain and store current masks status;
    unsigned getCurrentMask();

    // set mask
    void setMask(unsigned mask);

  };

}

#endif // _FloatingPointT_h_
