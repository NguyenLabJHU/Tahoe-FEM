// $Id: test.java,v 1.5 2002-07-24 20:55:54 recampb Exp $
import java.io.*;
import javax.swing.*;
import java.awt.event.*;
import java.awt.*;
import javax.swing.tree.*;

public class test extends JPanel implements ActionListener {

  // looks like an unnamed static function
  static {
	  System.loadLibrary("testClass");
	  //System.loadLibrary("toolbox");
	}
  
  public native void InitCpp();  
  public native void Print();
  public native void SetMinSc(int x);

 
  long cpp_obj;
  long console;
  int newNodeSuffix = 1;
  protected JButton testButton, b2;
  protected JTextField minScalarTF;

  public test(){

    InitCpp();
    Print();
    JTabbedPane tabbedPane = new JTabbedPane();

    JPanel rootPanel = new JPanel(false);
    JPanel framePanel = new JPanel(false);
    JPanel bodyPanel = new JPanel(false);
    JPanel bodyVarPanel = new JPanel(false);
    tabbedPane.addTab("Root Commands", rootPanel);
    tabbedPane.addTab("Frame Commands", framePanel);
    tabbedPane.addTab("Body Commands", bodyPanel);
    tabbedPane.addTab("Body Variables", bodyVarPanel);
    bodyVarPanel.setLayout(new GridLayout(2,2,20,5));

    minScalarTF = new JTextField("-99", 8);
    JLabel minScalarL = new JLabel("min_Scalar_Range");
    JButton bodyVarOKButton = new JButton("OK");
    JButton bodyVarCancelButton = new JButton("Cancel");
    bodyVarOKButton.addActionListener(this);
    bodyVarOKButton.setActionCommand("BodyVarOK");

    bodyVarPanel.add(minScalarL);
    bodyVarPanel.add(minScalarTF);
    bodyVarPanel.add(bodyVarOKButton);
    bodyVarPanel.add(bodyVarCancelButton);


    setLayout(new GridLayout(1, 2)); 
    add(tabbedPane);
 
    testButton = new JButton ("Exit");
    testButton.addActionListener(this);
    testButton.setActionCommand("Exit");
    rootPanel.add(testButton);

    b2 = new JButton ("Print value");
    b2.addActionListener(this);
    b2.setActionCommand("Print value");
    rootPanel.add(b2);
    

    final DynamicTree treePanel = new DynamicTree();
    populateTree(treePanel);

   

 JButton addButton = new JButton("Add");
        addButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                treePanel.addObject("New Node " + newNodeSuffix++);
            }
        });

        JButton removeButton = new JButton("Remove");
        removeButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                treePanel.removeCurrentNode();
            }
        });

        JButton clearButton = new JButton("Clear");
        clearButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                treePanel.clear();
            }
        });

        //Lay everything out.
        //setLayout(new BorderLayout());
        treePanel.setPreferredSize(new Dimension(300, 150));
        //add(treePanel, BorderLayout.CENTER);

        JPanel panel = new JPanel();
        panel.setLayout(new GridLayout(0,1));
        panel.add(addButton);
        panel.add(removeButton);
        panel.add(clearButton);
        add(panel, BorderLayout.EAST);


    //setLayout(new BorderLayout());
    //treePanel.setPreferredSize(new Dimension(300,150));
    add(treePanel);
  }

  public void actionPerformed(ActionEvent e){
    if (e.getActionCommand().equals("Exit")){
      System.exit(0);  
    }

    else if (e.getActionCommand().equals("Print value")){
      Print();

    }
    
    else if (e.getActionCommand().equals("BodyVarOK")){
      SetMinSc(Integer.parseInt(minScalarTF.getText()));
      
    }
    
  }


    public void populateTree(DynamicTree treePanel) {
        String p1Name = new String("Parent 1");
        String p2Name = new String("Parent 2");
        String c1Name = new String("Child 1");
        String c2Name = new String("Child 2");

        DefaultMutableTreeNode p1, p2;

        p1 = treePanel.addObject(null, p1Name);
        p2 = treePanel.addObject(null, p2Name);

        treePanel.addObject(p1, c1Name);
        treePanel.addObject(p1, c2Name);

        treePanel.addObject(p2, c1Name);
        treePanel.addObject(p2, c2Name);
    }



} // class test
