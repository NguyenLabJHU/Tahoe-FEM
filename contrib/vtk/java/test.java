// $Id: test.java,v 1.6 2002-07-29 21:11:50 recampb Exp $
import java.io.*;
import javax.swing.*;
import java.awt.event.*;
import java.awt.*;
import javax.swing.tree.*;
import javax.swing.event.*;

public class test extends JPanel implements ActionListener {

  // looks like an unnamed static function
  static {
	  System.loadLibrary("testClass");
	  //System.loadLibrary("toolbox");
	}
  
  public native void InitCpp();  
  public native void Print();
  public native void SetMinSc(int x);
  public native int GetMinSc();
  public native void AddScope(String s);

 
  long cpp_obj;
  long console;

  int newNodeSuffix = 1;
  protected JButton testButton, b2;
  protected JTextField minScalarTF;

  public test(){

    InitCpp();
    Print();
    final JTabbedPane tabbedPane = new JTabbedPane();

    final JPanel rootPanel = new JPanel();
    final JPanel framePanel = new JPanel();
    JPanel bodyPanel = new JPanel();
    JPanel bodyVarPanel = new JPanel();

    bodyVarPanel.setLayout(new GridLayout(2,2,20,5));

    minScalarTF = new JTextField(Integer.toString(GetMinSc()), 8);
    JLabel minScalarL = new JLabel("min_Scalar_Range");
    JButton bodyVarOKButton = new JButton("OK");
    JButton bodyVarCancelButton = new JButton("Cancel");
    bodyVarOKButton.addActionListener(this);
    bodyVarOKButton.setActionCommand("BodyVarOK");

    bodyVarPanel.add(minScalarL);
    bodyVarPanel.add(minScalarTF);
    bodyVarPanel.add(bodyVarOKButton);
    bodyVarPanel.add(bodyVarCancelButton);

    final JToolBar bodyToolBar = new JToolBar();
    JToolBar bodyDataToolBar = new JToolBar();
    JToolBar frameToolBar = new JToolBar();
    JToolBar rootToolBar = new JToolBar();

    JButton ShowContoursButton = new JButton ("ShowContours");
    JButton HideContoursButton = new JButton ("HideContours");
    JButton ShowCuttingButton = new JButton ("ShowCuttingPlane");
    JButton HideCuttingButton = new JButton ("HideCuttingPlane");

    ShowContoursButton.addActionListener(this);
    HideContoursButton.addActionListener(this);    
    ShowCuttingButton.addActionListener(this);
    HideCuttingButton.addActionListener(this);

    ShowContoursButton.setActionCommand("ShowContours");
    HideContoursButton.setActionCommand("HideContours");
    ShowCuttingButton.setActionCommand("ShowCutting");
    HideCuttingButton.setActionCommand("HideCutting");


    bodyToolBar.add(ShowContoursButton);
    bodyToolBar.add(HideContoursButton);
    bodyToolBar.add(ShowCuttingButton);
    bodyToolBar.add(HideCuttingButton);
    
    
    tabbedPane.addTab("Root Commands", rootPanel);
    tabbedPane.addTab("Frame Commands", framePanel);
    tabbedPane.addTab("Body Commands", bodyToolBar);
    tabbedPane.addTab("Body Variables", bodyVarPanel);
    

    setLayout(new GridLayout(1, 2)); 
     add(tabbedPane);
     //add(rootPanel);
     //add(bodyToolBar);
    
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
    //treePanel.setup();
    treePanel.getTree().addTreeSelectionListener(new TreeSelectionListener() {
                public void valueChanged(TreeSelectionEvent e) {
                    DefaultMutableTreeNode node = (DefaultMutableTreeNode)
                                       treePanel.getTree().getLastSelectedPathComponent();
                    
                    if (node == null) return;

                    Object nodeInfo = node.getUserObject();
		    if (((String)nodeInfo).equals( "0.body"))
		      tabbedPane.setSelectedComponent(bodyToolBar);
		    else if (((String)nodeInfo).equals("Console Root"))
		      tabbedPane.setSelectedComponent(rootPanel);
		    else if (((String)nodeInfo).equals("0.0.frame"))
		      tabbedPane.setSelectedComponent(framePanel);

		    System.out.println((String)nodeInfo);
		    
		}
                
    });
    
    
    JButton addButton = new JButton("Add");
        addButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                treePanel.addObject("New Node " + newNodeSuffix++);
		AddScope("New Node " + newNodeSuffix++);
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
    
    else if (e.getActionCommand().equals("ShowContours")){
      System.out.println("ShowContours");
    }
    
    else if (e.getActionCommand().equals("HideContours")){
      System.out.println("HideContours");
    }

    else if (e.getActionCommand().equals("ShowCutting")){
      System.out.println("ShowCuttingPlane");
    }
    else if (e.getActionCommand().equals("HideCutting")){
      System.out.println("HideCuttingPlane");
    }

  }


    public void populateTree(DynamicTree treePanel) {
        String f1Name = new String("0.0.frame");
	String b1Name = new String("0.body");
        String bd1Name = new String("0.body");

        DefaultMutableTreeNode p1, p2;

        p1 = treePanel.addObject(null, f1Name);
	p2 = treePanel.addObject(null, bd1Name);
        treePanel.addObject(p1, b1Name);

    }



} // class test
