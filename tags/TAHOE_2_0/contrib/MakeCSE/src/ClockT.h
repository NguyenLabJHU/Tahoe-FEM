// file: ClockT.h

#ifndef _SAWCLOCK_H_
#define _SAWCLOCK_H_

#include <time.h>
#include "ios_fwd_decl.h"
#include <iomanip.h>
#include "ArrayT.h"
#include "StringT.h"

using namespace Tahoe;

class ClockT
{
 public:
  ClockT (void) : fInitial (0) { };
  
  void Set (ArrayT<StringT>& names)
    { 
      fTime.Allocate (names.Length());
      fName.Allocate (names.Length());
      for (int i=0; i < names.Length(); i++)
	fName[i] = names[i];
      fTime = 0; 
    };

  void Initial (void) { fInitial = clock(); };

  void Sum (int index)
    { 
      clock_t t1 = clock();
      double secs = double (t1 - fInitial)/CLOCKS_PER_SEC;
      fTime[index] += secs/60;
      Initial();
    };

  void Print (ostream& out)
    {
      // find longest name
      int length = fName[0].Length();
      for (int i=1; i < fName.Length(); i++)
	if (fName[i].Length() > length) length = fName[i].Length();

      out << '\n';
      for (int j=0; j < fName.Length(); j++)
	out << setw(length) << fName[j] << ": " << fTime[j] << '\n';
    }

 private:
  double fInitial;
  ArrayT<double> fTime;
  ArrayT<StringT> fName;
};

#endif
