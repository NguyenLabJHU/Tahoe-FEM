// $Id: test.java,v 1.12 2002-08-13 10:11:20 paklein Exp $
import java.io.*;
import javax.swing.*;
import java.awt.event.*;
import java.awt.*;
import javax.swing.tree.*;
import javax.swing.event.*;

public class test extends JPanel implements ActionListener, ChangeListener {

  // looks like an unnamed static function
  static {
	  System.loadLibrary("testClass");
	  //System.loadLibrary("toolbox");
	}
  
  public native void InitCpp();  
  public native void Print();
  public native void SetMinSc(int x);
  public native int GetMinSc();
  public native void SetScope(String s);
  public native void Interact();
  public native void CommandLine();
  //public native void Rotate(int x, int y , int z);

 
  long cpp_obj;
  long console;
  long consoleObjects;

  int newNodeSuffix = 1;
  protected JButton testButton, b2, flipBookButton, nextTimeButton, prevTimeButton, selectTimeButton, addButton, removeButton, clearButton;
  protected JTextField minScalarTF, selectTimeTF, panXTF, panYTF, rotateXTF, rotateYTF, rotateZTF, zoomTF, selectTimeTFF, panXTFF, panYTFF, rotateXTFF, rotateYTFF,  rotateZTFF, zoomTFF;
  protected JPanel leftPanel, leftRootPanel, leftFramePanel, leftBodyPanel;
  protected JScrollPane leftRootScrollPane, leftFrameScrollPane, leftBodyScrollPane;
  protected DynamicTree treePanel;
  protected JTabbedPane tabbedPane;
  protected JSlider selectTimeSlider, selectTimeSliderF;
  protected JSplitPane /* splitPane, */ splitPane2;
  private boolean tabBool;
  protected JFrame parentFrame;

  public test(JFrame parentFrame){
    this.parentFrame = parentFrame;
    InitCpp();
    Print();
    tabbedPane = new JTabbedPane();
    tabbedPane.addChangeListener(this);
    tabBool = false;

    final JPanel rootPanel = new JPanel();
    final JPanel framePanel = new JPanel();
    JPanel bodyPanel = new JPanel();
    JPanel bodyVarPanel = new JPanel();

    JMenuBar menuBar = new JMenuBar();
    JMenu file = new JMenu("File");
    JMenu help = new JMenu("Help");
    JMenuItem exitItem = new JMenuItem("Exit");
    exitItem.addActionListener(this);
    exitItem.setActionCommand("Exit");
    JMenuItem addBodyItem = new JMenuItem("Add Body");
    JCheckBoxMenuItem frameNumsItem = new JCheckBoxMenuItem("Frame Numbers", false);
    JMenuItem removeBodyItem = new JMenuItem("Remove Body");
    JMenuItem saveItem = new JMenuItem("Save");
    JMenuItem saveFlipItem = new JMenuItem("Save Flip Book");
    JMenuItem windowSizeItem= new JMenuItem("Window Size");

    GridBagConstraints gbc = new GridBagConstraints();    
    gbc.anchor=gbc.NORTHWEST;
    
    file.add(addBodyItem);
    file.add(removeBodyItem);
    file.add(frameNumsItem);
    file.add(saveItem);
    file.add(saveFlipItem);
    file.add(windowSizeItem);
    file.add(exitItem);

    menuBar.add(file);
    menuBar.add(help);


    setLayout(new GridBagLayout());
   
    add(menuBar, new GridBagConstraints(0,0,4,1,0.0,0.05, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(0,0,0,0),0,0));
    
    add(tabbedPane,new GridBagConstraints(0,1,4,1,0.0,0.0,GridBagConstraints.NORTHWEST,GridBagConstraints.HORIZONTAL,new Insets(0,0,0,0),0,0));
    
    
    //bodyVarPanel.setLayout(new GridLayout(2,2,20,5));


    minScalarTF = new JTextField(Integer.toString(GetMinSc()), 8);
    JLabel minScalarL = new JLabel("min_Scalar_Range");
    JButton bodyVarOKButton = new JButton("OK");
    JButton bodyVarCancelButton = new JButton("Cancel");
    bodyVarOKButton.addActionListener(this);
    bodyVarOKButton.setActionCommand("BodyVarOK");

    bodyVarPanel.setLayout(new GridBagLayout()); 
    gbc.gridx=0; gbc.gridy=0; gbc.gridwidth=1;
    bodyVarPanel.add(minScalarL, gbc);
    gbc.gridx=1; gbc.gridy=0; gbc.gridwidth=1;
    bodyVarPanel.add(minScalarTF, gbc);

    JPanel buttonPanel=new JPanel();
    buttonPanel.add(bodyVarOKButton);
    buttonPanel.add(bodyVarCancelButton);
    gbc.gridx=0; gbc.gridy=1; gbc.gridwidth=2; gbc.anchor=gbc.CENTER;
    bodyVarPanel.add(buttonPanel,gbc);


    final JToolBar bodyToolBar = new JToolBar();
    JToolBar bodyDataToolBar = new JToolBar();
    JToolBar frameToolBar = new JToolBar();
    JToolBar rootToolBar = new JToolBar();

    //JButton ShowContoursButton = new JButton ("ShowContours");
    JButton ShowContoursButton = new JButton(new ImageIcon("../contButtonIconSmall.jpg"));
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
    

    //setLayout(new GridLayout(1, 2)); 
    //add(tabbedPane);
     //add(rootPanel);
     //add(bodyToolBar);

     rootPanel.setLayout(new GridBagLayout());
 
     //gbc.insets = new Insets(4,4,20,20);

    testButton = new JButton ("Exit");
    testButton.addActionListener(this);
    testButton.setActionCommand("Exit");
    gbc.gridx=0; gbc.gridy=1; gbc.gridwidth=1;
    rootPanel.add(testButton, gbc);


    b2 = new JButton ("Print value");
    b2.addActionListener(this);
    b2.setActionCommand("Print value");
    gbc.gridx=1; gbc.gridy=1; gbc.gridwidth=1;
    rootPanel.add(b2, gbc);
    
    JButton commandButton = new JButton("Command Line");
    commandButton.addActionListener(this);
    commandButton.setActionCommand("CommandLine");
    gbc.gridx=2;
    rootPanel.add(commandButton, gbc);

    flipBookButton = new JButton("Flip Book");
    flipBookButton.addActionListener(this);
    flipBookButton.setActionCommand("Flip Book");
    gbc.gridx=0; gbc.gridy=0; gbc.gridwidth=1;
    rootPanel.add(flipBookButton, gbc);

    nextTimeButton = new JButton("Next Time Step");
    nextTimeButton.addActionListener(this);
    nextTimeButton.setActionCommand("Next Time Step");
    gbc.gridx=1; gbc.gridy=0; gbc.gridwidth=1;
    rootPanel.add(nextTimeButton, gbc);

    prevTimeButton = new JButton("Previous Time Step");
    prevTimeButton.addActionListener(this);
    prevTimeButton.setActionCommand("Previous Time Step");
    gbc.gridx=2; gbc.gridy=0; gbc.gridwidth=1;
    rootPanel.add(prevTimeButton, gbc);   

    selectTimeSlider = new JSlider(0,10,0);
    selectTimeSlider.addChangeListener(this);
    JLabel selectTimeLabel = new JLabel("Time Step");
    selectTimeButton = new JButton("Select Time Step");
    selectTimeButton.addActionListener(this);
    selectTimeButton.setActionCommand("Select Time Step");
    gbc.gridx=3; gbc.gridy=0; gbc.gridwidth=1;
    rootPanel.add(selectTimeLabel, gbc);

    selectTimeTF = new JTextField("0", 5);
    selectTimeTF.setEnabled(true);
    gbc.gridx=4; gbc.gridy=0; gbc.gridwidth=1;
    rootPanel.add(selectTimeSlider, gbc);
    gbc.gridx=5; gbc.gridy=0; gbc.gridwidth=1;
    rootPanel.add(selectTimeTF, gbc);


    /* frame top panel stuff */
    framePanel.setLayout(new GridBagLayout());
    JButton flipBookButtonF = new JButton("Flip Book");
    flipBookButtonF.addActionListener(this);
    flipBookButtonF.setActionCommand("Flip Book");
    gbc.gridx=0; gbc.gridy=0; gbc.gridwidth=1;
    framePanel.add(flipBookButtonF, gbc);

    JButton nextTimeButtonF = new JButton("Next Time Step");
    nextTimeButtonF.addActionListener(this);
    nextTimeButtonF.setActionCommand("Next Time Step");
    gbc.gridx=1; gbc.gridy=0; gbc.gridwidth=1;
    framePanel.add(nextTimeButtonF, gbc);

    JButton prevTimeButtonF = new JButton("Previous Time Step");
    prevTimeButtonF.addActionListener(this);
    prevTimeButtonF.setActionCommand("Previous Time Step");
    gbc.gridx=2; gbc.gridy=0; gbc.gridwidth=1;
    framePanel.add(prevTimeButtonF, gbc);   

    selectTimeSliderF = new JSlider(0,10,0);
    selectTimeSliderF.addChangeListener(this);
    JLabel selectTimeLabelF = new JLabel("Time Step");
    JButton selectTimeButtonF = new JButton("Select Time Step");
    selectTimeButtonF.addActionListener(this);
    selectTimeButtonF.setActionCommand("Select Time Step");
    gbc.gridx=3; gbc.gridy=0; gbc.gridwidth=1;
    framePanel.add(selectTimeLabelF, gbc);

    selectTimeTFF = new JTextField("0", 5);
    selectTimeTFF.setEnabled(true);
    selectTimeTFF.addActionListener(this);
    gbc.gridx=4; gbc.gridy=0; gbc.gridwidth=1;
    framePanel.add(selectTimeSliderF, gbc);
    gbc.gridx=5; gbc.gridy=0; gbc.gridwidth=1;
    framePanel.add(selectTimeTFF, gbc);
 


    treePanel = new DynamicTree();
    populateTree(treePanel);

    JPanel rightPanel = new JPanel();
    rightPanel.setLayout(new GridLayout(0,1));

    addButton = new JButton("Add");
    addButton.addActionListener(this);
    addButton.setActionCommand("Add");

    removeButton = new JButton("Remove");
    removeButton.addActionListener(this);
    removeButton.setActionCommand("Remove");
    
    JButton clearButton = new JButton("Clear");
    clearButton.addActionListener(this);
    clearButton.setActionCommand("Clear");

    rightPanel.add(addButton);
    rightPanel.add(removeButton);
    rightPanel.add(clearButton);

    treePanel.add(rightPanel);





    //treePanel.setup();
    treePanel.getTree().addTreeSelectionListener(new TreeSelectionListener() {
                public void valueChanged(TreeSelectionEvent e) {
                    DefaultMutableTreeNode node = (DefaultMutableTreeNode)
                                       treePanel.getTree().getLastSelectedPathComponent();
                    
                    if (node == null) return;

                    Object nodeInfo = node.getUserObject();
		    if (((String)nodeInfo).equals( "0.body")){
		      tabbedPane.setSelectedComponent(bodyToolBar);
		      SetScope((String)nodeInfo);

		    }
		    else if (((String)nodeInfo).equals("Console Root")){
		      tabbedPane.setSelectedComponent(rootPanel);
		      //remove(leftPanel);
		      //add(leftRootPanel,new GridBagConstraints(0,2,1,1,0.05,1.0,GridBagConstraints.NORTHWEST,GridBagConstraints.BOTH,new Insets(0,0,0,0),0,0));
		      //splitPane.setLeftComponent(leftRootPanel);
		      SetScope((String)nodeInfo);
		      updateUI();

		    }
		    else if (((String)nodeInfo).equals("0.0.frame")){
		      tabbedPane.setSelectedComponent(framePanel);
		      //remove(leftPanel);
		      //add(leftRootPanel,new GridBagConstraints(0,2,1,1,0.05,1.0,GridBagConstraints.NORTHWEST,GridBagConstraints.BOTH,new Insets(0,0,0,0),0,0));
		      //splitPane.setLeftComponent(leftFrameScrollPane);
		      SetScope((String)nodeInfo);
		      updateUI();

		    }

		    System.out.println((String)nodeInfo);
		    
		}
                
    });
    
   
    /* left root panel stuff */


    JButton panButton = new JButton("Pan");
    panButton.addActionListener(this);
    panButton.setActionCommand("Pan");
    
    JButton rotateButton = new JButton("Rotate");
    rotateButton.addActionListener(this);
    rotateButton.setActionCommand("Rotate");
    
    JButton zoomButton = new JButton("Zoom");
    zoomButton.addActionListener(this);
    zoomButton.setActionCommand("Zoom");

    JLabel rotXLabel = new JLabel("X");
    JLabel rotYLabel = new JLabel("Y");
    JLabel rotZLabel = new JLabel("Z");

    JLabel panXLabel = new JLabel("X");
    JLabel panYLabel = new JLabel("Y");
    JLabel zoomFactorLabel = new JLabel("Factor");

    panXTF = new JTextField("0", 5);
    panYTF = new JTextField("0", 5);
    rotateXTF = new JTextField("0", 5);
    rotateYTF = new JTextField("0", 5);
    rotateZTF = new JTextField("0", 5);
    zoomTF = new JTextField("0", 5);

    JButton resetViewButton = new JButton("Reset View");
    JButton interactiveButton = new JButton("Interactive");
    JButton updateButton = new JButton("Update");
    interactiveButton.addActionListener(this);
    interactiveButton.setActionCommand("Interact");
    resetViewButton.addActionListener(this);
    resetViewButton.setActionCommand("Reset View");
    
     
        //treePanel.setPreferredSize(new Dimension(300, 150));
        //add(treePanel, BorderLayout.CENTER);

//     JPanel rotPanel = new JPanel();
//     rotPanel.setLayout(new GridLayout(3,2));
//     rotPanel.add(rotXLabel);
//     rotPanel.add(rotateXTF);
//     rotPanel.add(rotYLabel);
//     rotPanel.add(rotateYTF);
//     rotPanel.add(rotZLabel);
//     rotPanel.add(rotateZTF);
//     JPanel rotPanel2 = new JPanel();
//     rotPanel2.setLayout(new GridLayout(2,1));
//     rotPanel2.add(rotateButton);
//     rotPanel2.add(rotPanel);
    

	leftPanel = new JPanel();
	gbc.anchor=gbc.NORTHWEST;
        //leftPanel.setLayout(new GridLayout(0,1));
	leftPanel.setLayout(new GridBagLayout());
// 	gbc.gridx=0; gbc.gridy=0; gbc.gridwidth=1;
// 	leftPanel.add(addButton, gbc);
// 	gbc.gridx=1; gbc.gridy=0; gbc.gridwidth=1;
// 	leftPanel.add(removeButton, gbc);
// 	gbc.gridx=0; gbc.gridy=1; gbc.gridwidth=1;
// 	leftPanel.add(clearButton, gbc);

	leftRootPanel = new JPanel();
	leftRootScrollPane = new JScrollPane(leftRootPanel);
	gbc.anchor = gbc.NORTHEAST;
	gbc.insets=new Insets(4,4,4,4);
	leftRootPanel.setLayout(new GridBagLayout());
	gbc.gridx=1; gbc.gridy=0; gbc.gridwidth=1;
	leftRootPanel.add(rotateButton, gbc);
	//leftRootPanel.add(rotPanel2, gbc);
	gbc.gridx=0; gbc.gridy=1; gbc.gridwidth=1;
	leftRootPanel.add(rotXLabel, gbc);
	gbc.gridx=1; gbc.gridy=1; gbc.gridwidth=1;
	leftRootPanel.add(rotateXTF, gbc);
	gbc.gridx=0; gbc.gridy=2; gbc.gridwidth=1;
	leftRootPanel.add(rotYLabel, gbc);
	gbc.gridx=1; gbc.gridy=2; gbc.gridwidth=1;
	leftRootPanel.add(rotateYTF, gbc);
	gbc.gridx=0; gbc.gridy=3; gbc.gridwidth=1;
	leftRootPanel.add(rotZLabel, gbc);
	gbc.gridx=1; gbc.gridy=3; gbc.gridwidth=1;
	leftRootPanel.add(rotateZTF, gbc);

	leftRootPanel.add(new JSeparator(JSeparator.HORIZONTAL),new GridBagConstraints(0,4,2,1,0.0,0.0,GridBagConstraints.CENTER,GridBagConstraints.HORIZONTAL,new Insets(0,0,0,0),0,0));
	gbc.gridx=1; gbc.gridy=5; gbc.gridwidth=1;
	leftRootPanel.add(panButton, gbc);
	gbc.gridx=0; gbc.gridy=6; gbc.gridwidth=1;
	leftRootPanel.add(panXLabel, gbc);
	gbc.gridx=1; gbc.gridy=6; gbc.gridwidth=1;
	leftRootPanel.add(panXTF, gbc);
	gbc.gridx=0; gbc.gridy=7; gbc.gridwidth=1;
	leftRootPanel.add(panYLabel, gbc);
	gbc.gridx=1; gbc.gridy=7; gbc.gridwidth=1;
	leftRootPanel.add(panYTF, gbc);

	leftRootPanel.add(new JSeparator(JSeparator.HORIZONTAL),new GridBagConstraints(0,8,2,1,0.0,0.0,GridBagConstraints.CENTER,GridBagConstraints.HORIZONTAL,new Insets(0,0,0,0),0,0));
	
	gbc.gridx=1; gbc.gridy = 9;
	leftRootPanel.add(zoomButton, gbc);
	gbc.gridx=0; gbc.gridy = 10;
	leftRootPanel.add(zoomFactorLabel, gbc);
	gbc.gridx=1; gbc.gridy = 10;
	leftRootPanel.add(zoomTF, gbc);

	leftRootPanel.add(new JSeparator(JSeparator.HORIZONTAL),new GridBagConstraints(0,11,2,1,0.0,0.0,GridBagConstraints.CENTER,GridBagConstraints.HORIZONTAL,new Insets(0,0,0,0),0,0));
	
	gbc.gridx=0; gbc.gridy = 12; gbc.gridwidth=2;
	leftRootPanel.add(resetViewButton, gbc);
// 	gbc.gridx=0; gbc.gridy = 13;
// 	leftRootPanel.add(addButton, gbc);
// 	gbc.gridx=0; gbc.gridy = 14;
// 	leftRootPanel.add(removeButton, gbc);
// 	gbc.gridx=0; gbc.gridy = 15;
// 	leftRootPanel.add(clearButton, gbc);
	gbc.gridx=0; gbc.gridy = 13;
	leftRootPanel.add(interactiveButton, gbc);
	gbc.gridx=0; gbc.gridy = 14;
	leftRootPanel.add(updateButton, gbc);

	/* left frame panel stuff */
	JButton panButtonF = new JButton("Pan");
	panButtonF.addActionListener(this);
	panButtonF.setActionCommand("Pan");
	
	JButton rotateButtonF = new JButton("Rotate");
	rotateButtonF.addActionListener(this);
	rotateButtonF.setActionCommand("Rotate");
	
	JButton zoomButtonF = new JButton("Zoom");
	zoomButtonF.addActionListener(this);
	zoomButtonF.setActionCommand("Zoom");
	
	JLabel rotXLabelF = new JLabel("X");
	JLabel rotYLabelF = new JLabel("Y");
	JLabel rotZLabelF = new JLabel("Z");
	
	JLabel panXLabelF = new JLabel("X");
	JLabel panYLabelF = new JLabel("Y");
	JLabel zoomFactorLabelF = new JLabel("Factor");
	
	panXTFF = new JTextField("0", 5);
	panYTFF = new JTextField("0", 5);
	rotateXTFF = new JTextField("0", 5);
	rotateYTFF = new JTextField("0", 5);
	rotateZTFF = new JTextField("0", 5);
	zoomTFF = new JTextField("0", 5);
	
	JButton resetViewButtonF = new JButton("Reset View");
	JButton interactiveButtonF = new JButton("Interactive");
	JButton updateButtonF = new JButton("Update");
	
	JLabel bgLabel = new JLabel("Background Color");
	JRadioButton bgBlackButton = new JRadioButton("Black");
	bgBlackButton.addActionListener(this);
	bgBlackButton.setActionCommand("Black BG");
	bgBlackButton.setSelected(true);
	JRadioButton bgWhiteButton = new JRadioButton("White");
	bgWhiteButton.addActionListener(this);
	bgWhiteButton.setActionCommand("White BG");
	ButtonGroup bgButtonGroup = new ButtonGroup();
	bgButtonGroup.add(bgBlackButton);
	bgButtonGroup.add(bgWhiteButton);

	JLabel repLabel = new JLabel("Representation");
	JRadioButton repSurfButton = new JRadioButton("Surface");
	repSurfButton.setSelected(true);
	repSurfButton.addActionListener(this);
	repSurfButton.setActionCommand("Surf");
	JRadioButton repWireButton = new JRadioButton("Wire");
	repWireButton.addActionListener(this);
	repWireButton.setActionCommand("Wire");
	JRadioButton repPointButton = new JRadioButton("Points");
	repPointButton.addActionListener(this);
	repPointButton.setActionCommand("Points");
	ButtonGroup repButtons = new ButtonGroup();
	repButtons.add(repSurfButton);
	repButtons.add(repWireButton);
	repButtons.add(repPointButton);

	JLabel axesLabel = new JLabel("Axes");
	JLabel colorBarLabel = new JLabel ("Color Bar");
	JLabel nodeNumsLabel = new JLabel ("Node Numbers");
	JLabel elemNumsLabel = new JLabel ("Element Numbers");
	JRadioButton axesShowButton = new JRadioButton("Show");
	JRadioButton axesHideButton = new JRadioButton("Hide");
	ButtonGroup axesButtons = new ButtonGroup();
	JRadioButton colorBarShowButton = new JRadioButton("Show");
	JRadioButton colorBarHideButton = new JRadioButton("Hide");
	ButtonGroup colorBarButtons = new ButtonGroup();
	JRadioButton nodeNumsShowButton = new JRadioButton("Show");
	JRadioButton nodeNumsHideButton = new JRadioButton("Hide");
	ButtonGroup nodeNumsButtons = new ButtonGroup();
	JRadioButton elemNumsShowButton = new JRadioButton("Show");
	JRadioButton elemNumsHideButton = new JRadioButton("Hide");
	ButtonGroup elemNumsButtons = new ButtonGroup();
	axesButtons.add(axesShowButton);
	axesButtons.add(axesHideButton);
	colorBarButtons.add(colorBarShowButton);
	colorBarButtons.add(colorBarHideButton);
	nodeNumsButtons.add(nodeNumsShowButton);
	nodeNumsButtons.add(nodeNumsHideButton);
	elemNumsButtons.add(elemNumsShowButton);
	elemNumsButtons.add(elemNumsHideButton);
	axesHideButton.setSelected(true);
	colorBarHideButton.setSelected(true);
	nodeNumsHideButton.setSelected(true);
	elemNumsHideButton.setSelected(true);

	axesShowButton.addActionListener(this);
	axesHideButton.addActionListener(this);
	axesShowButton.setActionCommand("Show Axes");
	axesHideButton.setActionCommand("Hide Axes");
	colorBarShowButton.addActionListener(this);
	colorBarHideButton.addActionListener(this);
	colorBarShowButton.setActionCommand("Show Color Bar");
	colorBarHideButton.setActionCommand("Hide Color Bar");
	nodeNumsShowButton.addActionListener(this);
	nodeNumsHideButton.addActionListener(this);
	nodeNumsShowButton.setActionCommand("Show Node Nums");
	nodeNumsHideButton.setActionCommand("Hide Node Nums");
	elemNumsShowButton.addActionListener(this);
	elemNumsHideButton.addActionListener(this);
	elemNumsShowButton.setActionCommand("Show Elem Nums");
	elemNumsHideButton.setActionCommand("Hide Elem Nums");


	
	leftFramePanel = new JPanel();
	leftFrameScrollPane = new JScrollPane(leftFramePanel);
	gbc.anchor = gbc.NORTHEAST;
	gbc.insets=new Insets(4,4,4,4);
	leftFramePanel.setLayout(new GridBagLayout());
	gbc.gridx=1; gbc.gridy=0; gbc.gridwidth=1;
	leftFramePanel.add(rotateButtonF, gbc);
	gbc.gridx=0; gbc.gridy=1; gbc.gridwidth=1;
	leftFramePanel.add(rotXLabelF, gbc);
	gbc.gridx=1; gbc.gridy=1; gbc.gridwidth=1;
	leftFramePanel.add(rotateXTFF, gbc);
	gbc.gridx=0; gbc.gridy=2; gbc.gridwidth=1;
	leftFramePanel.add(rotYLabelF, gbc);
	gbc.gridx=1; gbc.gridy=2; gbc.gridwidth=1;
	leftFramePanel.add(rotateYTFF, gbc);
	gbc.gridx=0; gbc.gridy=3; gbc.gridwidth=1;
	leftFramePanel.add(rotZLabelF, gbc);
	gbc.gridx=1; gbc.gridy=3; gbc.gridwidth=1;
	leftFramePanel.add(rotateZTFF, gbc);

	leftFramePanel.add(new JSeparator(JSeparator.HORIZONTAL),new GridBagConstraints(0,4,2,1,0.0,0.0,GridBagConstraints.CENTER,GridBagConstraints.HORIZONTAL,new Insets(0,0,0,0),0,0));
	gbc.gridx=1; gbc.gridy=5; gbc.gridwidth=1;
	leftFramePanel.add(panButtonF, gbc);
	gbc.gridx=0; gbc.gridy=6; gbc.gridwidth=1;
	leftFramePanel.add(panXLabelF, gbc);
	gbc.gridx=1; gbc.gridy=6; gbc.gridwidth=1;
	leftFramePanel.add(panXTFF, gbc);
	gbc.gridx=0; gbc.gridy=7; gbc.gridwidth=1;
	leftFramePanel.add(panYLabelF, gbc);
	gbc.gridx=1; gbc.gridy=7; gbc.gridwidth=1;
	leftFramePanel.add(panYTFF, gbc);

	leftFramePanel.add(new JSeparator(JSeparator.HORIZONTAL),new GridBagConstraints(0,8,2,1,0.0,0.0,GridBagConstraints.CENTER,GridBagConstraints.HORIZONTAL,new Insets(0,0,0,0),0,0));
	
	gbc.gridx=1; gbc.gridy = 9;
	leftFramePanel.add(zoomButtonF, gbc);
	gbc.gridx=0; gbc.gridy = 10;
	leftFramePanel.add(zoomFactorLabelF, gbc);
	gbc.gridx=1; gbc.gridy = 10;
	leftFramePanel.add(zoomTFF, gbc);

	leftFramePanel.add(new JSeparator(JSeparator.HORIZONTAL),new GridBagConstraints(0,11,2,1,0.0,0.0,GridBagConstraints.CENTER,GridBagConstraints.HORIZONTAL,new Insets(0,0,0,0),0,0));
	
	gbc.gridx=0; gbc.gridy = 12; gbc.gridwidth=2;
	leftFramePanel.add(resetViewButtonF, gbc);
// 	gbc.gridx=0; gbc.gridy = 13;
// 	leftFramePanel.add(addButtonF, gbc);
// 	gbc.gridx=0; gbc.gridy = 14;
// 	leftFramePanel.add(removeButtonF, gbc);
// 	gbc.gridx=0; gbc.gridy = 15;
// 	leftFramePanel.add(clearButtonF, gbc);
	gbc.gridx=0; gbc.gridy = 13;
	leftFramePanel.add(interactiveButtonF, gbc);
	gbc.gridx=0; gbc.gridy = 14;
	leftFramePanel.add(updateButtonF, gbc);
	gbc.gridx=0; gbc.gridy = 15; gbc.gridwidth=1; gbc.anchor = gbc.NORTHWEST;
	leftFramePanel.add(bgLabel, gbc);
	gbc.gridx=0; gbc.gridy = 16;
	leftFramePanel.add(bgBlackButton, gbc);
	gbc.gridx=1; gbc.gridy = 16;
	leftFramePanel.add(bgWhiteButton, gbc);
	gbc.gridx=0; gbc.gridy = 17;
	leftFramePanel.add(repLabel, gbc);
	gbc.gridx=0; gbc.gridy = 18;
	leftFramePanel.add(repSurfButton, gbc);
	gbc.gridx=1; gbc.gridy = 18;
	leftFramePanel.add(repWireButton, gbc);
	gbc.gridx=0; gbc.gridy = 19;
	leftFramePanel.add(repPointButton, gbc);
	gbc.gridx=0; gbc.gridy = 20;
	leftFramePanel.add(axesLabel, gbc);
	gbc.gridx=0; gbc.gridy = 21;
	leftFramePanel.add(axesShowButton, gbc);
	gbc.gridx=1; gbc.gridy = 21;
	leftFramePanel.add(axesHideButton, gbc);
	gbc.gridx=0; gbc.gridy = 22;
	leftFramePanel.add(colorBarLabel, gbc);
	gbc.gridx=0; gbc.gridy = 23;
	leftFramePanel.add(colorBarShowButton, gbc);
	gbc.gridx=1; gbc.gridy = 23;
	leftFramePanel.add(colorBarHideButton, gbc);
	gbc.gridx=0; gbc.gridy = 24;
	leftFramePanel.add(nodeNumsLabel, gbc);
	gbc.gridx=0; gbc.gridy = 25;
	leftFramePanel.add(nodeNumsShowButton, gbc);
	gbc.gridx=1; gbc.gridy = 25;
	leftFramePanel.add(nodeNumsHideButton, gbc);
	gbc.gridx=0; gbc.gridy = 26;
	leftFramePanel.add(elemNumsLabel, gbc);
	gbc.gridx=0; gbc.gridy = 27;
	leftFramePanel.add(elemNumsShowButton, gbc);
	gbc.gridx=1; gbc.gridy = 27;
	leftFramePanel.add(elemNumsHideButton, gbc); 



	/* left body panel stuff */
      
	JButton interactiveButtonB = new JButton("Interactive");
	JButton updateButtonB = new JButton("Update");
	JButton cuttingButtonB = new JButton("Show Cutting Plane");
	interactiveButtonB.addActionListener(this);
	interactiveButtonB.setActionCommand("Interactive");
	updateButtonB.addActionListener(this);
	updateButtonB.setActionCommand("Update");
	cuttingButtonB.addActionListener(this);
	cuttingButtonB.setActionCommand("ShowCutting");

	JLabel repLabelB = new JLabel("Representation");
	JRadioButton repSurfButtonB = new JRadioButton("Surface");
	repSurfButtonB.setSelected(true);
	repSurfButtonB.addActionListener(this);
	repSurfButtonB.setActionCommand("Surf");
	JRadioButton repWireButtonB = new JRadioButton("Wire");
	repWireButtonB.addActionListener(this);
	repWireButtonB.setActionCommand("Wire");
	JRadioButton repPointButtonB = new JRadioButton("Points");
	repPointButtonB.addActionListener(this);
	repPointButtonB.setActionCommand("Points");
	ButtonGroup repButtonsB = new ButtonGroup();
	repButtonsB.add(repSurfButton);
	repButtonsB.add(repWireButton);
	repButtonsB.add(repPointButton);

	JLabel nodeNumsLabelB = new JLabel ("Node Numbers");
	JLabel elemNumsLabelB = new JLabel ("Element Numbers");

	JRadioButton nodeNumsShowButtonB = new JRadioButton("Show");
	JRadioButton nodeNumsHideButtonB = new JRadioButton("Hide");
	ButtonGroup nodeNumsButtonsB = new ButtonGroup();
	JRadioButton elemNumsShowButtonB = new JRadioButton("Show");
	JRadioButton elemNumsHideButtonB = new JRadioButton("Hide");
	ButtonGroup elemNumsButtonsB = new ButtonGroup();

	nodeNumsButtonsB.add(nodeNumsShowButtonB);
	nodeNumsButtonsB.add(nodeNumsHideButtonB);
	elemNumsButtonsB.add(elemNumsShowButtonB);
	elemNumsButtonsB.add(elemNumsHideButtonB);

	nodeNumsHideButtonB.setSelected(true);
	elemNumsHideButtonB.setSelected(true);

	nodeNumsShowButtonB.addActionListener(this);
	nodeNumsHideButtonB.addActionListener(this);
	nodeNumsShowButtonB.setActionCommand("Show Node Nums");
	nodeNumsHideButtonB.setActionCommand("Hide Node Nums");
	elemNumsShowButtonB.addActionListener(this);
	elemNumsHideButtonB.addActionListener(this);
	elemNumsShowButtonB.setActionCommand("Show Elem Nums");
	elemNumsHideButtonB.setActionCommand("Hide Elem Nums");


	JLabel contoursLabelB = new JLabel ("Contours");
	JLabel glyphsLabelB = new JLabel ("Glyphs");

	JRadioButton contoursShowButtonB = new JRadioButton("Show");
	JRadioButton contoursHideButtonB = new JRadioButton("Hide");
	ButtonGroup contoursButtonsB = new ButtonGroup();
	JRadioButton glyphsShowButtonB = new JRadioButton("Show");
	JRadioButton glyphsHideButtonB = new JRadioButton("Hide");
	ButtonGroup glyphsButtonsB = new ButtonGroup();

	contoursButtonsB.add(contoursShowButtonB);
	contoursButtonsB.add(contoursHideButtonB);
	glyphsButtonsB.add(glyphsShowButtonB);
	glyphsButtonsB.add(glyphsHideButtonB);

	contoursHideButtonB.setSelected(true);
	glyphsHideButtonB.setSelected(true);

	contoursShowButtonB.addActionListener(this);
	contoursHideButtonB.addActionListener(this);
	contoursShowButtonB.setActionCommand("Show Contours");
	contoursHideButtonB.setActionCommand("Hide Contours");
	glyphsShowButtonB.addActionListener(this);
	glyphsHideButtonB.addActionListener(this);
	glyphsShowButtonB.setActionCommand("Show Glyphs");
	glyphsHideButtonB.setActionCommand("Hide Glyphs");


	leftBodyPanel = new JPanel();
	leftBodyScrollPane = new JScrollPane(leftBodyPanel);
	gbc.anchor = gbc.NORTHEAST;
	gbc.insets=new Insets(4,4,4,4);
	leftBodyPanel.setLayout(new GridBagLayout());
 	gbc.gridx=1; gbc.gridy=0; gbc.gridwidth=1;
 	leftBodyPanel.add(interactiveButtonB, gbc);
 	gbc.gridx=1; gbc.gridy=1; gbc.gridwidth=1;
 	leftBodyPanel.add(updateButtonB, gbc);
	gbc.gridx=0; gbc.gridy=2; gbc.gridwidth=1;gbc.anchor = gbc.NORTHWEST;
	leftBodyPanel.add(repLabelB, gbc);
 	gbc.gridx=0; gbc.gridy=3; gbc.gridwidth=1;
 	leftBodyPanel.add(repSurfButtonB, gbc);
 	gbc.gridx=1; gbc.gridy=3; gbc.gridwidth=1;
 	leftBodyPanel.add(repWireButtonB, gbc);
 	gbc.gridx=0; gbc.gridy=4; gbc.gridwidth=1;
 	leftBodyPanel.add(repPointButtonB, gbc);
 	gbc.gridx=0; gbc.gridy=5; gbc.gridwidth=1;
 	leftBodyPanel.add(nodeNumsLabelB, gbc);

// 	leftBodyPanel.add(new JSeparator(JSeparator.HORIZONTAL),new GridBagConstraints(0,4,2,1,0.0,0.0,GridBagConstraints.CENTER,GridBagConstraints.HORIZONTAL,new Insets(0,0,0,0),0,0));
	gbc.gridx=0; gbc.gridy=6; gbc.gridwidth=1;
 	leftBodyPanel.add(nodeNumsShowButtonB, gbc);
 	gbc.gridx=1; gbc.gridy=6; gbc.gridwidth=1;
 	leftBodyPanel.add(nodeNumsHideButtonB, gbc);
 	gbc.gridx=0; gbc.gridy=7; gbc.gridwidth=1;
 	leftBodyPanel.add(elemNumsLabelB, gbc);
 	gbc.gridx=0; gbc.gridy=8; gbc.gridwidth=1;
 	leftBodyPanel.add(elemNumsShowButtonB, gbc);
 	gbc.gridx=1; gbc.gridy=8; gbc.gridwidth=1;
 	leftBodyPanel.add(elemNumsHideButtonB, gbc);

// 	leftBodyPanel.add(new JSeparator(JSeparator.HORIZONTAL),new GridBagConstraints(0,8,2,1,0.0,0.0,GridBagConstraints.CENTER,GridBagConstraints.HORIZONTAL,new Insets(0,0,0,0),0,0));
	
 	gbc.gridx=0; gbc.gridy = 9;
 	leftBodyPanel.add(contoursLabelB, gbc);
 	gbc.gridx=0; gbc.gridy = 10;
 	leftBodyPanel.add(contoursShowButtonB, gbc);
 	gbc.gridx=1; gbc.gridy = 10;
 	leftBodyPanel.add(contoursHideButtonB, gbc);

// 	leftFramePanel.add(new JSeparator(JSeparator.HORIZONTAL),new GridBagConstraints(0,11,2,1,0.0,0.0,GridBagConstraints.CENTER,GridBagConstraints.HORIZONTAL,new Insets(0,0,0,0),0,0));
	
// 	gbc.gridx=0; gbc.gridy = 12; gbc.gridwidth=2;
// 	leftFramePanel.add(resetViewButtonF, gbc);
// // 	gbc.gridx=0; gbc.gridy = 13;
// // 	leftFramePanel.add(addButtonF, gbc);
// // 	gbc.gridx=0; gbc.gridy = 14;
// // 	leftFramePanel.add(removeButtonF, gbc);
// // 	gbc.gridx=0; gbc.gridy = 15;
// // 	leftFramePanel.add(clearButtonF, gbc);
// 	gbc.gridx=0; gbc.gridy = 13;
// 	leftFramePanel.add(interactiveButtonF, gbc);
// 	gbc.gridx=0; gbc.gridy = 14;
// 	leftFramePanel.add(updateButtonF, gbc);
// 	gbc.gridx=0; gbc.gridy = 15; gbc.gridwidth=1; gbc.anchor = gbc.NORTHWEST;
// 	leftFramePanel.add(bgLabel, gbc);
// 	gbc.gridx=0; gbc.gridy = 16;
// 	leftFramePanel.add(bgBlackButton, gbc);
// 	gbc.gridx=1; gbc.gridy = 16;
// 	leftFramePanel.add(bgWhiteButton, gbc);
// 	gbc.gridx=0; gbc.gridy = 17;
// 	leftFramePanel.add(repLabel, gbc);
// 	gbc.gridx=0; gbc.gridy = 18;
// 	leftFramePanel.add(repSurfButton, gbc);
// 	gbc.gridx=1; gbc.gridy = 18;
// 	leftFramePanel.add(repWireButton, gbc);
// 	gbc.gridx=0; gbc.gridy = 19;
// 	leftFramePanel.add(repPointButton, gbc);
// 	gbc.gridx=0; gbc.gridy = 20;
// 	leftFramePanel.add(axesLabel, gbc);
// 	gbc.gridx=0; gbc.gridy = 21;
// 	leftFramePanel.add(axesShowButton, gbc);
// 	gbc.gridx=1; gbc.gridy = 21;
// 	leftFramePanel.add(axesHideButton, gbc);
// 	gbc.gridx=0; gbc.gridy = 22;
// 	leftFramePanel.add(colorBarLabel, gbc);
// 	gbc.gridx=0; gbc.gridy = 23;
// 	leftFramePanel.add(colorBarShowButton, gbc);
// 	gbc.gridx=1; gbc.gridy = 23;
// 	leftFramePanel.add(colorBarHideButton, gbc);
// 	gbc.gridx=0; gbc.gridy = 24;
// 	leftFramePanel.add(nodeNumsLabel, gbc);
// 	gbc.gridx=0; gbc.gridy = 25;
// 	leftFramePanel.add(nodeNumsShowButton, gbc);
// 	gbc.gridx=1; gbc.gridy = 25;
// 	leftFramePanel.add(nodeNumsHideButton, gbc);
// 	gbc.gridx=0; gbc.gridy = 26;
// 	leftFramePanel.add(elemNumsLabel, gbc);
// 	gbc.gridx=0; gbc.gridy = 27;
// 	leftFramePanel.add(elemNumsShowButton, gbc);
// 	gbc.gridx=1; gbc.gridy = 27;
// 	leftFramePanel.add(elemNumsHideButton, gbc); 
       

// 	leftPanel.add(addButton);
//         leftPanel.add(removeButton);
//         leftPanel.add(clearButton);




	JPanel renderPanel = new JPanel();
	JLabel renderPic = new JLabel(new ImageIcon("../presCutCont3.jpg"));
	renderPanel.add(renderPic);
	renderPic.setPreferredSize(new Dimension(700,850));

//Create a split pane with the two scroll panes in it.
			//splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
            //                           leftRootPanel, renderPanel);
            //splitPane.setOneTouchExpandable(true);
            //splitPane.setDividerLocation(150);

            //Provide minimum sizes for the two components in the split pane
            Dimension minimumSize = new Dimension(150, 250);
            leftRootPanel.setMinimumSize(minimumSize);
            renderPanel.setMinimumSize(minimumSize);

//             splitPane2 = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
//                                        splitPane, treePanel);

//PAKLEIN
//	    splitPane2 = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, leftRootPanel, treePanel);
	    splitPane2 = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, leftRootScrollPane, treePanel);
            splitPane2.setOneTouchExpandable(true);
            splitPane2.setDividerLocation(400);

            //Provide minimum sizes for the two components in the split pane
            Dimension minimumSize2 = new Dimension(300, 250);
            treePanel.setMinimumSize(minimumSize);
	    //splitPane.setMinimumSize(minimumSize2);
            //renderPanel.setMinimumSize(minimumSize);
    //setLayout(new BorderLayout());
	    //treePanel.setPreferredSize(new Dimension(150,150));
	//add(treePanel);
// 	add(treePanel,new GridBagConstraints(1,2,1,1,0.05,1.0,GridBagConstraints.NORTHWEST,GridBagConstraints.BOTH,new Insets(0,0,0,0),0,0));
// 	add(renderPanel,new GridBagConstraints(2,2,1,1,0.8,1.0,GridBagConstraints.NORTHWEST,GridBagConstraints.BOTH,new Insets(0,0,0,0),0,0));
// 	add(leftRootPanel,new GridBagConstraints(0,2,1,1,0.05,1.0,GridBagConstraints.NORTHWEST,GridBagConstraints.BOTH,new Insets(0,0,0,0),0,0));
// 	add(new JSeparator(JSeparator.VERTICAL),new GridBagConstraints(1,2,1,3,0.0,0.0,GridBagConstraints.CENTER,GridBagConstraints.VERTICAL,new Insets(0,0,0,0),0,0));
	    //	add(splitPane,new GridBagConstraints(0,2,1,1,0.95,1.0,GridBagConstraints.NORTHWEST,GridBagConstraints.BOTH,new Insets(0,0,0,0),0,0));
	add(splitPane2,new GridBagConstraints(0,2,1,1,1.0,.95,GridBagConstraints.NORTHWEST,GridBagConstraints.BOTH,new Insets(0,0,0,0),0,0));
  }
  
  public void actionPerformed(ActionEvent e){
    if (e.getActionCommand().equals("Exit")){
      System.exit(0);  
    }

    else if (e.getActionCommand().equals("Print value")){
      Print();
      //leftPanel.remove(addButton);
      //updateUI();
      
    }
    
    else if (e.getActionCommand().equals("CommandLine")){
    CommandLine();
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
       JDialog temp = new JDialog(parentFrame, true);
       temp.setTitle("Show Cutting Plane");
       temp.setSize(400, 300);
       JPanel cutPanel = new JPanel();
       JButton cutOKButton = new JButton("OK");
       JButton cutCancelButton = new JButton("Cancel");
       cutPanel.add(cutOKButton);
       cutPanel.add(cutCancelButton);
       temp.getContentPane().add(cutPanel);
   
       temp.pack();
       temp.show();
       
    }
    else if (e.getActionCommand().equals("HideCutting")){
      System.out.println("HideCuttingPlane");
    }

    else if (e.getActionCommand().equals("Flip Book")){
      System.out.println("Flip Book");
    }
    else if (e.getActionCommand().equals("Next Time Step")){
      System.out.println("Next Time Step");
    }
    else if (e.getActionCommand().equals("Previous Time Step")){
      System.out.println("Previous Time Step");
    }
    else if (e.getActionCommand().equals("Select Time Step")){
      System.out.println("Select Time Step");
    }
    else if (e.getActionCommand().equals("Add")){
	treePanel.addObject("New Node " + newNodeSuffix++);
	//AddScope("New Node " + newNodeSuffix++);
    }
    else if (e.getActionCommand().equals("Remove")){
      treePanel.removeCurrentNode();
    }
    else if (e.getActionCommand().equals("Clear")){
       treePanel.clear();
    }
    else if (e.getActionCommand().equals("Interact")){
      Interact();
    }
    else if (e.getActionCommand().equals("Black BG")){
      System.out.println("Black BG");
    }
    else if (e.getActionCommand().equals("White BG")){
      System.out.println("White BG");
    }    
    else if (e.getActionCommand().equals("Surf")){
      System.out.println("Surf");
    }
    else if (e.getActionCommand().equals("Wire")){
      System.out.println("Wire");
    }
    else if (e.getActionCommand().equals("Points")){
      System.out.println("Points");
    }
    else if (e.getActionCommand().equals("Show Axes")){
      System.out.println("Show Axes");
    }
    else if (e.getActionCommand().equals("Hide Axes")){
      System.out.println(e.getActionCommand());
    }
    else if (e.getActionCommand().equals("Show Color Bar")){
      System.out.println(e.getActionCommand());
    }
    else if (e.getActionCommand().equals("Hide Color Bar")){
      System.out.println(e.getActionCommand());
    }
    else if (e.getActionCommand().equals("Show Node Nums")){
      System.out.println(e.getActionCommand());
    }
    else if (e.getActionCommand().equals("Hide Node Nums")){
      System.out.println(e.getActionCommand());
    }
    else if (e.getActionCommand().equals("Show Elem Nums")){
      System.out.println(e.getActionCommand());
    }
    else if (e.getActionCommand().equals("Hide Elem Nums")){
      System.out.println(e.getActionCommand());
    }

  }

public void stateChanged(ChangeEvent e)
{
	if(e.getSource()==selectTimeSlider){
		selectTimeTF.setText(Integer.toString((int)selectTimeSlider.getValue()));
	}
	else if(e.getSource()==selectTimeSliderF){
		selectTimeTFF.setText(Integer.toString((int)selectTimeSliderF.getValue()));
	}

	else if (e.getSource()==tabbedPane){
		if (tabbedPane.getTitleAt(tabbedPane.getSelectedIndex()).equals("Root Commands")){

			/* for some reason this event occurs when the app is started */
			if (tabBool)
				splitPane2.setLeftComponent(leftRootScrollPane);
			else { 
				tabBool=true;
				}
	}

	else if (tabbedPane.getTitleAt(tabbedPane.getSelectedIndex()).equals("Root Commands"))
		splitPane2.setLeftComponent(leftRootScrollPane);
	else if (tabbedPane.getTitleAt(tabbedPane.getSelectedIndex()).equals("Frame Commands"))
		splitPane2.setLeftComponent(leftFrameScrollPane);
	else if (tabbedPane.getTitleAt(tabbedPane.getSelectedIndex()).equals("Body Commands"))
		splitPane2.setLeftComponent(leftBodyScrollPane);
		updateUI();
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
