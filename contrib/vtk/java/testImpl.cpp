// $Id: testImpl.cpp,v 1.7 2002-08-08 16:50:32 paklein Exp $
//#include <StubPreamble.h>
// not needed for 1.3.1

#include "test.h"
#include "testClass.h"
// #include "vtkRenderWindow.h"

#include "iArrayT.h"
#include "iConsoleT.h"
#include "iConsoleObjectT.h"
#include "CommandSpecT.h"

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
	iConsoleT* console = new iConsoleT(log_file, *b, NULL, false);
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

	/* test communication with the console */
	fid = env->GetFieldID(cls, "console", "J");
	if (fid == 0) {
    	return;
  	}
	p_long = env->GetLongField(obj, fid);
	iConsoleT* p_console = (iConsoleT*) p_long;

	cout << " current scope name: " << p_console->Scope() << endl;

	/* try "help" command */	
	StringT empty_line;
	const CommandSpecT* help_command = p_console->iResolveCommand("help", empty_line);
	if (!help_command) throw;
	if (p_console->iDoCommand(*help_command, empty_line))
		cout << " help OK\n" << endl;
	else
		cout << " help NOT OK\n" << endl;

	/* try "read" command */
	StringT read_arguments = "file.test";
	const CommandSpecT* read_command = p_console->iResolveCommand("read", read_arguments);
	if (!read_command) throw;
	cout << "read command: ";
	read_command->WriteCommand(cout);
	cout << endl;
	if (p_console->iDoCommand(*read_command, empty_line))
		cout << " read OK\n" << endl;
	else
		cout << " read NOT OK\n" << endl;

	/* try "echo" command */
	StringT echo_arguments;
	echo_arguments.Append("X");
	echo_arguments.Append(" ", 5);
	echo_arguments.Append("Y");
	echo_arguments.Append(" ", 10);
	const CommandSpecT* echo_command = p_console->iResolveCommand("echo", echo_arguments);
	if (p_console->iDoCommand(*echo_command, empty_line))
		cout << " echo OK\n" << endl;
	else
		cout << " echo NOT OK\n" << endl;
		
	/* call command of current console scope */
	iConsoleObjectT& current_scope = p_console->Current();
	cout << "current console object name: " << current_scope.iName() << endl;
	echo_command = current_scope.iResolveCommand("echo", echo_arguments);
	if (echo_command && current_scope.iDoCommand(*echo_command, empty_line))
		cout << " echo OK\n" << endl;
	else
		cout << " echo NOT OK\n" << endl;

	const CommandSpecT* list_command = current_scope.iResolveCommand("list", empty_line);
	if (list_command && current_scope.iDoCommand(*list_command, empty_line))
		cout << " list OK\n" << endl;
	else
		cout << " list NOT OK\n" << endl;

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
