/* $Id: testClass.cpp,v 1.2 2002-07-20 01:59:11 paklein Exp $ */
#include "testClass.h"

testClass::testClass(int a): a_(a)
{ }

void testClass::Print(ostream& out)
{
	out << "\n testClass::Print: field value is " << a_ << endl;
}
