// $Id: test.java,v 1.3 2002-07-20 01:59:11 paklein Exp $
import java.io.*;

public class test {

  // looks like an unnamed static function
  static {
	  System.loadLibrary("testClass");
	}
  
  public native void InitCpp();  
  public native void Print();
  
	long cpp_obj;
}
