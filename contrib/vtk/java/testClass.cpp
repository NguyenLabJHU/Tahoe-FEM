/* $Id: testClass.cpp,v 1.3 2002-07-24 19:01:28 recampb Exp $ */
#include "testClass.h"

testClass::testClass(int a): a_(a)
{ }

void testClass::Print(ostream& out)
{
	out << "\n testClass::Print: field value is " << a_ << endl;
}

void testClass::SetA(int x)
{
a_ = x;
}
  
