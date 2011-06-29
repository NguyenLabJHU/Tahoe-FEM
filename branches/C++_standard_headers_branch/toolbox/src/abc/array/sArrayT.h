#ifndef _SARRAYT_H_
#define _SARRAYT_H_

#include "StringT.h"
#include "ArrayT.h"

namespace Tahoe {

class sArrayT : public ArrayT<StringT>
{
 public:
  /* constructor*/
  sArrayT(void);
  explicit sArrayT(int length);
  sArrayT(const sArrayT& source);
  sArrayT(const ArrayT<StringT>& source);

  /* assignment operators */
  sArrayT& operator=(const sArrayT& RHS);
  sArrayT& operator=(const StringT& value);

  /** return 1 if the value is present, 0 otherwise */
  int HasValue (const StringT& value) const;

  /** find the first occurence of a value in the array
   * \param value value to match
   * \param index returns with the index of the first occurence, if found, -1 otherwise
   * \return 1 if value found, 0 otherwise */
  int HasValue (const StringT& value, int& index) const;

  /** write wrapped, using width = largest component's width */
  void WriteWrapped (ostream& out, int wrap) const;
};

} // namespace Tahoe

#endif
