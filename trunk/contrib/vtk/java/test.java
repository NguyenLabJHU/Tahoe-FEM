// $Id: test.java,v 1.1 2002-07-17 23:28:42 paklein Exp $
import java.io.*;

public class test {

	// looks like an unnamed static function
	static {
		System.loadLibrary("testClass");
	}
		
	// "native" function to construct C++ object
	test() {
		initCppSide();
	}
	
	public native void Print();
	public native void initCppSide();
	
	// pointer to the C++ object
	private long cpp_ptr;	
}
