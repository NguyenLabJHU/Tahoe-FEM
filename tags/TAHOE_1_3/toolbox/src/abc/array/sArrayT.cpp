#include "sArrayT.h"
#include <iomanip.h>

using namespace Tahoe;

sArrayT::sArrayT (void) { }
sArrayT::sArrayT (int length) : ArrayT<StringT> (length) { }
sArrayT::sArrayT (const sArrayT& source) : ArrayT<StringT>(source) { }
sArrayT::sArrayT (const ArrayT<StringT>& source) : ArrayT<StringT>(source) { }

sArrayT& sArrayT::operator=(const sArrayT& RHS) 
{ 
  ArrayT<StringT>::operator=(RHS); 
  return *this;
}
sArrayT& sArrayT::operator=(const StringT& value) 
{ 
  ArrayT<StringT>::operator=(value); 
  return *this;
}

int sArrayT::HasValue (const StringT& value) const
{
  const StringT* p = Pointer();
  for (int i=0; i < Length(); i++, p++)
    {
      int length = p->StringLength();
      int vlength = value.StringLength();
      int test = (length > vlength) ? length : vlength;
      if (strncmp (*p, value.Pointer(), test) == 0)
	return 1;
    }
  return 0;
}

int sArrayT::HasValue (const StringT& value, int& index) const
{
  index  = -1;
  const StringT* p = Pointer();
  for (int i = 0; i < Length() && index < 0; i++, p++)
    {
      int length = p->StringLength();
      int vlength = value.StringLength();
      int test = (length > vlength) ? length : vlength;
      if (strncmp (*p, value.Pointer(), test) == 0)
	{
	  index = i;
	  return 1;
	}
    }
  return 0;			
}

void sArrayT::WriteWrapped (ostream& out, int wrap) const
{
  int width = 0;
  const StringT *p = Pointer();
  for (int i=0; i < Length(); i++, p++)
    {
      int ewidth = p->Length();
      if (ewidth > width) 
	width = ewidth;
    }

  p = Pointer();
  for (int j=0, count = 0; j < Length(); j++)
    {
      if (count == wrap)
	{
	  out << '\n';
	  count = 0;
	}
      count ++;
      out << setw (width) << *p++;
    }
}
