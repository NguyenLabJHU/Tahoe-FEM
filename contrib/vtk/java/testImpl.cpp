// $Id: testImpl.cpp,v 1.4 2002-07-20 01:59:11 paklein Exp $
//#include <StubPreamble.h>
// not needed for 1.3.1

#include "test.h"
#include "testClass.h"

JNIEXPORT void JNICALL Java_test_InitCpp(JNIEnv * env, jobject obj)
{
	jclass cls = env->GetObjectClass(obj);
	jfieldID fid = env->GetFieldID(cls, "cpp_obj", "J");
	if (fid == 0) {
    	return;
  	}

	int val = -99;
	cout << "\n Java_test_InitCpp: storing a " << -99 << endl;
  	testClass* p = new testClass(val);
	env->SetLongField(obj, fid, jlong(p));
}

JNIEXPORT void JNICALL Java_test_Print(JNIEnv * env, jobject obj)
{
	cout << "\n Java_test_Print: Hello!" << endl;

	jclass cls = env->GetObjectClass(obj);
	jfieldID fid = env->GetFieldID(cls, "cpp_obj", "J");
	if (fid == 0) {
    	return;
  	}
	
	jlong p_long = env->GetLongField(obj, fid);
	testClass* p = (testClass*) p_long;
	p->Print(cout);
	return;
}
