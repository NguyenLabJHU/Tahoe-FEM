/* $Id: testClass.cpp,v 1.1 2002-07-17 23:28:43 paklein Exp $ */
#include "testClass.h"

testClass::testClass(int a): a_(a)
{ }

void testClass::Print(ostream& out)
{
	out << a_ << endl;
}
