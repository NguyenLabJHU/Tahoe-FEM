// $Id: testImpl.cpp,v 1.6 2002-07-29 21:11:50 recampb Exp $
//#include <StubPreamble.h>
// not needed for 1.3.1

#include "test.h"
#include "testClass.h"
// #include "vtkRenderWindow.h"

#include "iArrayT.h"
#include "iConsoleT.h"
#include "iConsoleObjectT.h"

using namespace Tahoe;

JNIEXPORT void JNICALL Java_test_InitCpp(JNIEnv * env, jobject obj)
{
	jclass cls = env->GetObjectClass(obj);
	jfieldID fid = env->GetFieldID(cls, "cpp_obj", "J");
	jfieldID fid2 = env->GetFieldID(cls, "console", "J");
	if (fid == 0) {
    	return;
  	}

	int val = -99;
	//int val2 = 100;
	cout << "\n Java_test_InitCpp: storing a " << -99 << endl;


	iArrayT test(10);
	test.SetValueToPosition();
	test++;
	cout << test.wrap(3) << endl;

	iConsoleObjectT *a = new iConsoleObjectT;
	iConsoleObjectT	*b = new iConsoleObjectT;
	iConsoleObjectT	*c = new iConsoleObjectT;
	iConsoleObjectT	*d = new iConsoleObjectT;
	a->iSetName("a");
	b->iSetName("b");
	c->iSetName("c");
	d->iSetName("d");

	c->iAddSub(*d);
	a->iAddSub(*b);
	a->iAddSub(*c);

	StringT log_file = "testClass.log";
	//iConsoleT console(log_file, *b);
	iConsoleT* console = new iConsoleT(log_file, *b);
	//console_(log_file, *b);

  	testClass* p = new testClass(val);
	env->SetLongField(obj, fid, jlong(p));
	env->SetLongField(obj, fid2, jlong(console)); 
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

JNIEXPORT void JNICALL Java_test_SetMinSc(JNIEnv * env, jobject obj, jint x)
  {
	jclass cls = env->GetObjectClass(obj);
	jfieldID fid = env->GetFieldID(cls, "cpp_obj", "J");
	if (fid == 0) {
    	return;
  	}

	
	jlong p_long = env->GetLongField(obj, fid);
	testClass* p = (testClass*) p_long;
	p->SetA(x);
	env->SetLongField(obj, fid, jlong(p));
	p->Print(cout);
	return;
	

  }

JNIEXPORT jint JNICALL Java_test_GetMinSc(JNIEnv * env, jobject obj)
{
  jclass cls = env->GetObjectClass(obj);
  jfieldID fid = env->GetFieldID(cls, "cpp_obj", "J");
  if (fid == 0) {
    return 0;
  }
  
  jlong p_long = env->GetLongField(obj, fid);
  testClass* p = (testClass*) p_long;
  jint temp = p->GetA();
  return temp;

}

JNIEXPORT void JNICALL Java_test_AddScope(JNIEnv * env, jobject obj, jstring s)
{
  jclass cls = env->GetObjectClass(obj);
  jfieldID fid = env->GetFieldID(cls, "console", "J");
  if (fid == 0) {
    return;
  }
  
  jlong p_long = env->GetLongField(obj, fid);
  iConsoleT* p = (iConsoleT*) p_long;
  StringT test = "test";
  StringT& test2 = (StringT&)test;
  
  cout << "Add scope" << endl;
  return;
  

}
