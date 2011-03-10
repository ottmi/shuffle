/*
	Output Selector is the GUI panel where users decide whether they want to output a file
	If they choose to, they are required to input the compatibilty siglevel and/or MNIC by 
	which sites will be included and the name the file should be called
*/
import java.awt.*;
import javax.swing.*;
import java.awt.event.*;
import java.lang.*;
import java.util.*;
import java.io.*;

public class OutputSelector extends JPanel implements ActionListener {
	
	private LinkedList<Double> siglevels; // compatibility significance level - only sites with a POC less than or equal to this will be included in the output alignment
	private double compthresh; // compatibility threshold - only sites with an overall compatibility greater than or equal to this will be included in the output alignment
	private int maxchanges; //  maximum MNIC - only sites with an MNIC value less than this will be included in the output alignment
	private int entropy; // maximum entropy - only sites with an entropy value less than this will be included in the output alignment
	private int gens; // number of generations the randomization procedure should be conducted for
	private boolean output; // whether or not an alignment file is to be outputted
	private boolean summary; // whether or not a file is to be outputted
	private boolean draw; // whether or not an compatibility plot is to be drawn
	private boolean nss; // whether or not the neighbour similarity score is to be calculated
	//private boolean PC;
	private JPanel outPanel;
	private JRadioButton yes, no, yes2, no2, yes3, no3, yes4, no4;
	private JCheckBox siglevel_button, compthresh_button, changes_button, entropy_button;
	private JTextField siglevel_box, compthresh_box, gens_box, changes_box, filename_box, entropy_box;
	
	public OutputSelector() {
		super(new BorderLayout());
		output=false;
		siglevels = new LinkedList();
		compthresh = 0;
		maxchanges = 21; 
		entropy=0;
		//PC=false;
		gens=100;
		summary=false;
		draw=true;
		nss=false;
				
		// creates radio buttons for choosing if the file should be outputted
		no = new JRadioButton("No");
		no.setActionCommand("0");
		no.setSelected(true);
		yes = new JRadioButton("Yes");
		yes.setActionCommand("1");
		no2 = new JRadioButton("No");
		no2.setActionCommand("2");
		no2.setSelected(true);
		yes2 = new JRadioButton("Yes");
		yes2.setActionCommand("3");
		no3 = new JRadioButton("No");
		no3.setActionCommand("4");
		no3.setSelected(true);
		yes3 = new JRadioButton("Yes");
		yes3.setActionCommand("5");
		no4 = new JRadioButton("No");
		no4.setActionCommand("6");
		yes4 = new JRadioButton("Yes");
		yes4.setActionCommand("7");
		yes4.setSelected(true);
		
		
		ButtonGroup g1 = new ButtonGroup();
		g1.add(no);
		g1.add(yes);
		
		ButtonGroup g2 = new ButtonGroup();
		g2.add(no2);
		g2.add(yes2);
		
		ButtonGroup g3 = new ButtonGroup();
		g3.add(no3);
		g3.add(yes3);
		
		ButtonGroup g4 = new ButtonGroup();
		g4.add(no4);
		g4.add(yes4);
		
		// creates radio buttons for choosing if the file outputted should be created based on a 
		// compatibility siglevel and/ or a MNIC max
		// each radio button has a corresponding textfield where users can input the siglevel/MNIC
		siglevel_button = new JCheckBox("POC");
		siglevel_button.setActionCommand("1");
		siglevel_button.setEnabled(false);
		siglevel_box = new JTextField("0.05");
		siglevel_box.setEnabled(false);
		compthresh_button =  new JCheckBox("Min Compatibility (%)");
		compthresh_button.setActionCommand("0");
		compthresh_button.setEnabled(false);
		compthresh_box = new JTextField("0");
		compthresh_box.setEnabled(false);
		gens_box = new JTextField("100");
		gens_box.setEnabled(false);
		changes_button = new JCheckBox("Max MNIC");
		changes_button.setActionCommand("2");
		changes_button.setEnabled(false);
		changes_box = new JTextField();
		changes_box.setEnabled(false);
		filename_box = new JTextField("outfile");
		filename_box.setEnabled(false);
		entropy_button = new JCheckBox("Max Entropy");
		entropy_button.setActionCommand("3");
		entropy_button.setEnabled(false);
		entropy_box = new JTextField();
		entropy_box.setEnabled(false);
		
		/*
			Specific output Listener
		*/
		class OutputListener implements ActionListener
		{
			// called when the siglevel/MNIC button is changed
			public void actionPerformed(ActionEvent e) 
			{
				updateState();

			}
		}
		
		// add action listener to buttons
		no.addActionListener(this);
		yes.addActionListener(this);
		no2.addActionListener(this);
		yes2.addActionListener(this);
		no3.addActionListener(this);
		yes3.addActionListener(this);
		no4.addActionListener(this);
		yes4.addActionListener(this);
		siglevel_button.addActionListener(new OutputListener());
		compthresh_button.addActionListener(new OutputListener());
		changes_button.addActionListener(new OutputListener());
		entropy_button.addActionListener(new OutputListener());
		//probcomp_button.addActionListener(new OutputListener());
		
		// add buttons and textfields to panel
		outPanel = new JPanel(new GridLayout(16,2));
		outPanel.add(new JLabel("Draw Plot?"));
		outPanel.add(new JLabel(""));
		outPanel.add(yes4);
		outPanel.add(no4);
		outPanel.add(new JLabel("Output Alignment?"));
		outPanel.add(new JLabel(""));
		outPanel.add(yes);
		outPanel.add(no);
		outPanel.add(new JLabel("Include sites:"));
		outPanel.add(new JLabel(""));
		outPanel.add(compthresh_button);
		outPanel.add(compthresh_box);
		outPanel.add(siglevel_button);
		outPanel.add(new JLabel(""));
		outPanel.add(new JLabel("#	POC Threshold(s):"));
		outPanel.add(siglevel_box);
		outPanel.add(new JLabel("#	Randomizations:"));
		outPanel.add(gens_box);
		outPanel.add(changes_button);
		outPanel.add(changes_box);
		outPanel.add(entropy_button);
		outPanel.add(entropy_box);
		outPanel.add(new JLabel("File name:"));
		outPanel.add(filename_box);
		outPanel.add(new JLabel("Output Site Summary?"));
		outPanel.add(new JLabel(""));
		outPanel.add(yes2);
		outPanel.add(no2);
		outPanel.add(new JLabel("Calc. Neighbour Similarity?"));
		outPanel.add(new JLabel(""));
		outPanel.add(yes3);
		outPanel.add(no3);

		add(outPanel, BorderLayout.CENTER);
		
	}
	
	/*
		Called by the action listeners.
		If a radio button is changed, updates the state coinicide with the new status
		Eg. if the 'yes' radio button is selected, the siglevel and MNIC radio buttons become activated
	*/
	public void updateState()
	{
		if(output){
			siglevel_button.setEnabled(true);
			compthresh_button.setEnabled(true);
			changes_button.setEnabled(true);
			entropy_button.setEnabled(true);
			filename_box.setEnabled(true);
			if(compthresh_button.isSelected()){
				compthresh_box.setEnabled(true);
			}
			else{
				compthresh_box.setEnabled(false);
			}
			if(siglevel_button.isSelected()){
				siglevel_box.setEnabled(true);
				gens_box.setEnabled(true);
				//probcomp_button.setEnabled(false);
			}
			else {
				siglevel_box.setEnabled(false);
				gens_box.setEnabled(false);
			}
			if(changes_button.isSelected()){
				changes_box.setEnabled(true);
			}
			else{
				changes_box.setEnabled(false);
			}
			if(entropy_button.isSelected()){
				entropy_box.setEnabled(true);
			}
			else{
				entropy_box.setEnabled(false);
			}
		}
		else{
			siglevel_button.setEnabled(false);
			compthresh_button.setEnabled(false);
			changes_button.setEnabled(false);
			entropy_button.setEnabled(false);
			siglevel_box.setEnabled(false);
			gens_box.setEnabled(false);
			changes_box.setEnabled(false);
			entropy_box.setEnabled(false);
			filename_box.setEnabled(false);
		}		
	}
	
	/*
		Main Action Listener
	*/
	public void actionPerformed(ActionEvent e)
	{
		if(Integer.parseInt(e.getActionCommand()) < 2){
			if(Integer.parseInt(e.getActionCommand()) == 0)
			{	
				output = false;
			}
			else output = true;
			updateState();
		}else{
			if(Integer.parseInt(e.getActionCommand()) == 2){summary = false;}
			else if(Integer.parseInt(e.getActionCommand()) == 3){summary = true;}
			else if(Integer.parseInt(e.getActionCommand()) == 4){nss = false;}
			else if(Integer.parseInt(e.getActionCommand()) == 5){nss = true;}
			else if(Integer.parseInt(e.getActionCommand()) == 6){draw = false;}
			else if(Integer.parseInt(e.getActionCommand()) == 7){draw = true;}
		}
	}

	
	// Helper method: Returns whether the a new alignment is to be outputted
	public boolean getOutput() { return output;}
	public boolean getSummary() { return summary;}
	public boolean getNSS() { return nss;}
	public boolean getDraw() { return draw;}	
	/*
		Helper method: Returns the current siglevel of the new alignment to be outputted
	*/
	public LinkedList<Double> getPOCThreshs() 
	{ 
		if(siglevel_button.isSelected()){
			String str = siglevel_box.getText();
			String[] strlist = str.split(",");
			for(int i=0;i<strlist.length;i++){
				try{
					siglevels.add(new Double(strlist[i].trim()));
				}
				catch(NumberFormatException nfe){
				}
			}
		}
		return siglevels;
	}
	
		/*
		Helper method: Returns the current raw compatibility threshold of the new alignment to be outputted
	*/
	public double getCompThresh() 
	{ 
		String str = compthresh_box.getText();
		try{
			compthresh = new Double(str);
		}
		catch(NumberFormatException nfe){
		}
		return compthresh;
	}
	
		/*
		Helper method: Returns the number of randomizations to be used to determine site significance
	*/
	public int getGenerations() 
	{ 
		String str = gens_box.getText();
		try{
			gens = new Integer(str);
		}
		catch(NumberFormatException nfe){
		}
		return gens;
	}
	
	/*
		Helper method: Returns the current max MNIC of the new alignment to be outputted
	*/
	public int getMaxChanges() 
	{
		if(siglevel_button.isEnabled()){
			String str = changes_box.getText();
			try{
				maxchanges = new Integer(str);
			}
			catch(NumberFormatException nfe){
			}
		}
		return maxchanges;
	}
	
	/*
		Helper method: Returns the current max Entropy of the new alignment to be outputted
	*/
	public int getEntropy() 
	{
		if(entropy_button.isEnabled()){
			String str = entropy_box.getText();
			try{
				entropy = new Integer(str);
			}
			catch(NumberFormatException nfe){
			}
		}
		return entropy;
	}
	
	
	/*
		returns the name of the file to be outputted
	*/
	public String getFilename() 
	{ 
		return filename_box.getText();
	}
	
	
	public void resetSiglevels()
	{
		siglevels = new LinkedList();
	}
}
