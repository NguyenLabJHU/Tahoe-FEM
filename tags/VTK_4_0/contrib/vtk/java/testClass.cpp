/* $Id: testClass.cpp,v 1.4 2002-07-29 21:11:50 recampb Exp $ */
#include "testClass.h"
#include "iConsoleT.h"

testClass::testClass(int a): a_(a)
{}

void testClass::Print(ostream& out)
{
	out << "\n testClass::Print: field value is " << a_ << endl;
}

void testClass::SetA(int x)
{
a_ = x;
}
  
int testClass::GetA(){
  return a_;
}


