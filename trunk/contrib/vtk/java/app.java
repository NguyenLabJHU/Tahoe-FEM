import java.io.*;
import javax.swing.*;
import java.awt.event.*;
import java.awt.*;

public class app {

  protected test a;

  public app(){
    a = new test();
    // construct the internal C++ object
    a.InitCpp();
    
    
    // see if it prints
    a.Print();

  }


public static void main(String args[]) throws IOException {

	    System.out.print("\n Hello world\n");

                JFrame frame = new JFrame("VTK for Tahoe");

                frame.addWindowListener(new WindowAdapter() {
		  public void windowClosing(WindowEvent e) {System.exit(0);}
                });
		
		JLabel emptyLabel = new JLabel("TEST");
		emptyLabel.setPreferredSize(new Dimension(175, 100));
		emptyLabel.setHorizontalAlignment(JLabel.CENTER);
		//frame.getContentPane().add(emptyLabel, BorderLayout.CENTER);		

		frame.getContentPane().add(new test());
                frame.pack();
                frame.setVisible(true);

	}



  
}
