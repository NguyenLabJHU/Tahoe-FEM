// $Id: testImpl.TESTING.cpp,v 1.1 2002-09-22 18:13:29 paklein Exp $
#include "test.h"
#include "VTKConsoleT.h"
#include "iConsoleT.h"

#include "iArrayT.h"

JNIEXPORT void JNICALL Java_test_InitCpp(JNIEnv * env, jobject obj)
{
//TEMP
cout << " sizeof(jlong): " << sizeof(jlong) << endl;
cout << " sizeof(void*): " << sizeof(void*) << endl;

	jclass cls = env->GetObjectClass(obj);

  	iArrayT* p = new iArrayT(10);
  	p->SetValueToPosition();
	jfieldID fid = env->GetFieldID(cls, "cpp_obj", "J");
	env->SetLongField(obj, fid, jlong(p));
	
#if 1
	/* construct VTK console object */
	StringT log_file = "testClass.log";	
	ArrayT<StringT> arguments;
	VTKConsoleT* vtk_console = new VTKConsoleT(arguments);	
	iConsoleT* console = new iConsoleT(log_file, *vtk_console, NULL, false);
	jfieldID fid2 = env->GetFieldID(cls, "console", "J");
	env->SetLongField(obj, fid2, jlong(console));
#endif
	
//TEMP
cout << " Java_test_InitCpp: DONE" << endl;
}

JNIEXPORT void JNICALL Java_test_Print(JNIEnv * env, jobject obj)
{
	jclass cls = env->GetObjectClass(obj);
	jfieldID fid = env->GetFieldID(cls, "cpp_obj", "J");
	if (fid == 0) {
    	return;
  	}
  	
	jlong p_long = env->GetLongField(obj, fid);
	iArrayT* p = (iArrayT*) p_long;
	cout << " array:\n" << (*p) << endl;
}

JNIEXPORT void JNICALL Java_test_SetMinSc(JNIEnv * env, jobject obj, jint x)
  {
return;
  }

JNIEXPORT jint JNICALL Java_test_GetMinSc(JNIEnv * env, jobject obj)
{
    return 0;
}

JNIEXPORT void JNICALL Java_test_SetScope(JNIEnv * env, jobject obj, jstring s)
{
    return;
}

JNIEXPORT void JNICALL Java_test_Interact(JNIEnv * env, jobject obj)
{
     return;
}

JNIEXPORT void JNICALL Java_test_CommandLine(JNIEnv * env, jobject obj)
{
     return;
}

JNIEXPORT void JNICALL Java_test_DoCommand(JNIEnv *env, 
	jobject obj, jstring command, jstring arguments)
{
	return;
}
