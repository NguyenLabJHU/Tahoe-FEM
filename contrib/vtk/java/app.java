import java.io.*;

public class app {

	public static void main(String args[]) throws IOException {

	    System.out.print("\n Hello world\n");
	    
	    // make test object
		test a = new test();

		// construct the internal C++ object
	    a.InitCpp();
	    
	    // see if it prints
	    a.Print();
	}
}
