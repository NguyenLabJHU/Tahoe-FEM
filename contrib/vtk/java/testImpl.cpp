// $Id: testImpl.cpp,v 1.1 2002-07-17 23:28:43 paklein Exp $
#include <StubPreamble.h>
#include "test.h"
#include "testClass.h"

JNIEXPORT void JNICALL Java_test_initCppSide(JNIEnv *, jobject javaObj)
{
	testClass* tc = new testClass(5);
	unhand(javaObj)->cpp_ptr = (long) tc;
}

JNIEXPORT void JNICALL Java_test_Print(JNIEnv *, jobject javaObj)
{
	testClass* tc = (testClass*) unhand(javaObj)->cpp_ptr;
	tc->Print(cout);
}
