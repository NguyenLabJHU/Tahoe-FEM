// $Id: testImpl.cpp,v 1.3 2002-07-19 17:28:25 recampb Exp $
//#include <StubPreamble.h>
// not needed for 1.3.1

#include "test.h"
#include "testClass.h"


JNIEXPORT void JNICALL Java_test_Print(JNIEnv *, jobject javaObj)
{
  //testClass* tc = (testClass*) unhand(javaObj)->cpp_ptr;
  //testClass* tc = (testClass*) javaObj->cpp_ptr;
  //	tc->Print(cout);
  cout << "Hello world" << endl;
  return;

}
