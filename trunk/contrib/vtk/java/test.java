// $Id: test.java,v 1.2 2002-07-19 17:28:25 recampb Exp $
import java.io.*;

public class test {

  // looks like an unnamed static function
  static {
	  System.loadLibrary("testClass");
	}
  
  public native void Print();
  
  public static void main(String[] args){
    new test().Print();
  }
}
