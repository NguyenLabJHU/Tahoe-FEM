// $Id: testImpl.cpp,v 1.14 2002-08-14 00:20:47 paklein Exp $
#include "test.h"
#include "testClass.h"
// #include "vtkRenderWindow.h"

#include "iArrayT.h"
#include "iConsoleT.h"
#include "iConsoleObjectT.h"
#include "CommandSpecT.h"
#include "AutoArrayT.h"
#include "VTKConsoleT.h"

using namespace Tahoe;

static iConsoleObjectT& GetCurrentConsoleObject(JNIEnv * env, jobject obj)
{
	jclass cls = env->GetObjectClass(obj);
	jfieldID fid = env->GetFieldID(cls, "console", "J");
  	if (fid == 0) {
		cout << " GetCurrentConsoleObject: error resolving field ID" << endl;
		throw;
	}
  
	jlong p_long = env->GetLongField(obj, fid);
	iConsoleT* p_console = (iConsoleT*) p_long;
	return p_console->Current();
}

JNIEXPORT void JNICALL Java_test_InitCpp(JNIEnv * env, jobject obj)
{
	jclass cls = env->GetObjectClass(obj);
	jfieldID fid = env->GetFieldID(cls, "cpp_obj", "J");
	jfieldID fid2 = env->GetFieldID(cls, "console", "J");
	jfieldID fid3 = env->GetFieldID(cls, "consoleObjects", "J");
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

#if 0
	iConsoleObjectT *a = new iConsoleObjectT;
	iConsoleObjectT	*b = new iConsoleObjectT;
	iConsoleObjectT	*c = new iConsoleObjectT;
	iConsoleObjectT	*d = new iConsoleObjectT;
	a->iSetName("Root");
	b->iSetName("0.0.frame");
	c->iSetName("0.body");
 	d->iSetName("0.body");
	
	a->iAddSub(*b);
	a->iAddSub(*c);
	b->iAddSub(*d);
#endif

	StringT log_file = "testClass.log";
	//iConsoleT* console = new  iConsoleT(log_file, *b, NULL, false);

#if 1
	/* construct VTK console object */
	ArrayT<StringT> arguments;
	VTKConsoleT* vtk_console = new VTKConsoleT(arguments);
	
	iConsoleT* console = new iConsoleT(log_file, *vtk_console, NULL, false);
#endif

  	testClass* p = new testClass(val);
	env->SetLongField(obj, fid, jlong(p));
	env->SetLongField(obj, fid2, jlong(console)); 

#if 1
	/* call command of current console scope */
	iConsoleObjectT& current_scope = console->Current();
	cout << "current console object name: " << current_scope.iName() << endl;

	/* try loading a body */
	StringT empty_line;
	StringT read_arguments = "billet.io0.geo";
	const CommandSpecT* read_command = current_scope.iResolveCommand("AddBody", read_arguments);
	if (read_command && current_scope.iDoCommand(*read_command, empty_line))
		cout << " AddBody OK\n" << endl;
	else
		cout << " AddBody NOT OK\n" << endl;

	read_arguments = "x 20";
//	const CommandSpecT* rotate_command = current_scope.iResolveCommand("Interactive", empty_line);
//	current_scope.iDoCommand(*rotate_command, empty_line);
	const CommandSpecT* update_command = current_scope.iResolveCommand("Update", empty_line);

	//const CommandSpecT* interact_command = current_scope.iResolveCommand("Interactive", empty_line);
	if (update_command && current_scope.iDoCommand(*update_command, empty_line))
		cout << " Update OK\n" << endl;
	else
		cout << " Update NOT OK\n" << endl;
#endif
}

JNIEXPORT void JNICALL Java_test_Print(JNIEnv * env, jobject obj)
{
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
	echo_arguments.Append("\"X");
	echo_arguments.Append(" ", 5);
	echo_arguments.Append(" Y");
	echo_arguments.Append(" ", 10);
	echo_arguments.Append("\"");	
	const CommandSpecT* echo_command = p_console->iResolveCommand("echo", echo_arguments);
	if (p_console->iDoCommand(*echo_command, empty_line))
		cout << " echo OK\n" << endl;
	else
		cout << " echo NOT OK\n" << endl;

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

JNIEXPORT void JNICALL Java_test_SetScope(JNIEnv * env, jobject obj, jstring s)
{
  jclass cls = env->GetObjectClass(obj);
  jfieldID fid = env->GetFieldID(cls, "console", "J");
  if (fid == 0) {
    return;
  }
  
  jlong p_long = env->GetLongField(obj, fid);
  iConsoleT* p_console = (iConsoleT*) p_long;

	
  StringT empty_line;
  StringT scope = "::Root";
  //StringT temp = s;
//   if (!(s.equals("Console Root")))
//      scope.Append(": ");
  
  p_console->iResolveCommand("::Root:0.0.frame", empty_line);


  cout << " current scope name: " << p_console->Scope() << endl;

  
}

JNIEXPORT void JNICALL Java_test_Interact(JNIEnv * env, jobject obj)
{
   jclass cls = env->GetObjectClass(obj);
   jfieldID fid = env->GetFieldID(cls, "console", "J");
   if (fid == 0) {
     return;
   }
  
   jlong p_long = env->GetLongField(obj, fid);
   iConsoleT* p_console = (iConsoleT*) p_long;
   iConsoleObjectT& current_scope = p_console->Current();
  
  	/* start synchronized block */
  	env->MonitorEnter(obj);
  
   StringT empty_line;
   StringT read_args = "x 20";
   const CommandSpecT* interact_command = current_scope.iResolveCommand("Interactive", empty_line);
	if (interact_command && current_scope.iDoCommand(*interact_command, empty_line))
		cout << " Interactive OK\n" << endl;
 	else
 		cout << " Interactive NOT OK\n" << endl;

  	/* end synchronized block */
  	env->MonitorExit(obj);

// 	jclass cls = env->GetObjectClass(obj);
// 	jfieldID fid = env->GetFieldID(cls, "console", "J");
// 	jlong p_long = env->GetLongField(obj, fid);	

// 	/* test communication with the console */

// 	if (fid == 0) {
//     	return;
//   	}
// 	p_long = env->GetLongField(obj, fid);
// 	iConsoleT* p_console = (iConsoleT*) p_long;

// 	cout << " current scope name: " << p_console->Scope() << endl;

// 	/* try "help" command */	
// 	StringT empty_line;
// 	const CommandSpecT* help_command = p_console->iResolveCommand("help", empty_line);
// 	if (!help_command) throw;
// 	if (p_console->iDoCommand(*help_command, empty_line))
// 		cout << " help OK\n" << endl;
// 	else
// 		cout << " help NOT OK\n" << endl;

// 	/* try "read" command */
// 	StringT read_arguments = "file.test";
// 	const CommandSpecT* read_command = p_console->iResolveCommand("read", read_arguments);
// 	if (!read_command) throw;
// 	cout << "read command: ";
// 	read_command->WriteCommand(cout);
// 	cout << endl;
// 	if (p_console->iDoCommand(*read_command, empty_line))
// 		cout << " read OK\n" << endl;
// 	else
// 		cout << " read NOT OK\n" << endl;

// 	/* try "echo" command */
// 	StringT echo_arguments;
// 	echo_arguments.Append("\"X");
// 	echo_arguments.Append(" ", 5);
// 	echo_arguments.Append(" Y");
// 	echo_arguments.Append(" ", 10);
// 	echo_arguments.Append("\"");	
// 	const CommandSpecT* echo_command = p_console->iResolveCommand("echo", echo_arguments);
// 	if (p_console->iDoCommand(*echo_command, empty_line))
// 		cout << " echo OK\n" << endl;
// 	else
// 		cout << " echo NOT OK\n" << endl;
		
	
// #if 1
// 	/* call command of current console scope */
// 	iConsoleObjectT& current_scope = p_console->Current();
// 	cout << "current console object name: " << current_scope.iName() << endl;

// 	/* try loading a body */
// 	read_arguments = "billet.io0.geo";
// 	read_command = current_scope.iResolveCommand("AddBody", read_arguments);
// 	if (read_command && current_scope.iDoCommand(*read_command, empty_line))
// 		cout << " AddBody OK\n" << endl;
// 	else
// 		cout << " AddBody NOT OK\n" << endl;

//  	//StringT 
// 	read_arguments = "x 20";
// 	//StringT empty_line;
//  	const CommandSpecT* rotate_command = current_scope.iResolveCommand("Rotate", read_arguments);
// 	current_scope.iDoCommand(*rotate_command, empty_line);
// 	const CommandSpecT* update_command = current_scope.iResolveCommand("Update", empty_line);
	
	
// 	//const CommandSpecT* interact_command = current_scope.iResolveCommand("Interactive", empty_line);
// 	if (update_command && current_scope.iDoCommand(*update_command, empty_line))
// 		cout << " Update OK\n" << endl;
// 	else
// 		cout << " Update NOT OK\n" << endl;
// #endif

// 	return;
}

JNIEXPORT void JNICALL Java_test_CommandLine(JNIEnv * env, jobject obj)
{
   jclass cls = env->GetObjectClass(obj);
   jfieldID fid = env->GetFieldID(cls, "console", "J");
   if (fid == 0) {
     return;
   }
  
   jlong p_long = env->GetLongField(obj, fid);
   iConsoleT* p_console = (iConsoleT*) p_long;

	/* let console run */
	p_console->DoInteractive();

	cout << "\n Control returned to GUI" << endl;
}

JNIEXPORT void JNICALL Java_test_DoCommand(JNIEnv *env, 
	jobject obj, jstring command, jstring arguments)
{
    const char *command_str = env->GetStringUTFChars(command, 0);
    const char *arguments_str = env->GetStringUTFChars(arguments, 0);

	cout << "   command " << "(" << strlen(command_str) <<  "):" << command_str << '\n'
	     << " arguments " << "(" << strlen(arguments_str) <<  "):" << arguments_str << endl;

	/* get current console object */
	iConsoleObjectT& current = GetCurrentConsoleObject(env, obj);

	/* resolve command */
	StringT com = command_str;
	StringT arg = arguments_str;	
	const CommandSpecT* command_spec = current.iResolveCommand(com, arg);
	StringT empty_line;
	if (command_spec && current.iDoCommand(*command_spec, empty_line))
		cout << " " << command_str << " OK\n" << endl;
 	else
		cout << " " << command_str << " FAILED\n" << endl;

    env->ReleaseStringUTFChars(command, command_str);
    env->ReleaseStringUTFChars(arguments, arguments_str);
}
